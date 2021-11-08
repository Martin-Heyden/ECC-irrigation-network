function [sys] = simulate_system(N,sys_choice,contr_set,params,Dist,T,x0)
%Simulate an irrigation network
%   N number of nodes, 
%   sys_choice is the model to be used for each pool (1 - pool 1, otherwise
%      pool 2)
%   contr_set is the choise of controller (1 struct, 2 P, 3 third order LQ
%   params is the design parameters for the controller
%   Dist is the planned disturbances (NxT matrix)
%   T is the simulation horizon
%   x0 is the initial condition

if nargin == 6
    x0 = zeros(N,1);
end

z = tf('z',60);
%%% Data for the two pools %%%
%%pool 9%%
%first order data - used for control design
tau_additional = 10; %Same for both pools
tau_first_9 = 2;
b9 = 0.069;
c9 = 0.063;

%third order data
ci = [0.137 -0.155 0.053];
cip1 = [-0.190 0.333 -0.175];%c_{i+1}
alfa1 = 0.978;
alfa2 = 0.468;
tau_i = 3;
g9 = [ci(1)*z^(-tau_i) + ci(2)*z^(-tau_i-1)+ci(3)*z^(-tau_i-2), ...
	cip1(1) + cip1(2)*z^(-1)+cip1(3)*z^(-2)]...
	/(z-1-alfa1*(1-2*z^(-1)+z^(-2))-alfa2*(1-z^(-1)));
%%pool 10%%

%First order%
tau_first_10 = 15;
b10 = 0.0142 * 1.5;
c10 = 0.0156; 
%Third oder%
ci = [0.134 -0.244 0.114]; %todo
cip1 = [-0.101 0.185 -0.087];%c_{i+1}
alfa1 = 0.314;
alfa2 = 0.814;
tau_i = 16;
g10 = [ci(1)*z^(-tau_i) + ci(2)*z^(-tau_i-1)+ci(3)*z^(-tau_i-2), ...
	cip1(1) + cip1(2)*z^(-1)+cip1(3)*z^(-2)]...
	/(z-1-alfa1*(1-2*z^(-1)+z^(-2))-alfa2*(1-z^(-1)));
%%% System  generation %%%
g = {g9,g10};
b_vec = zeros(1,N);
c_vec = zeros(1,N);
tau = zeros(1,N);
for i = 1:N
	if sys_choice(i) == 1
		b_vec(i) = b9;
		c_vec(i) = c9;
		tau(i) = tau_first_9;
	else
		b_vec(i) = b10;
		c_vec(i) = c10;
		tau(i) = tau_first_10;
	end
end




%%% Low pass Filter Design %%%
fc = 0.003/2/pi;
fs = 1/60;
[A,B,C,D] = butter(4,fc/(fs/2));
lp = ss(A,B,C,D,60);


%%% Controller Synthesis %%%
if  contr_set == 1
	r = params.r;
	qvec = ones(N,1);
	R1 = params.R1;
	R2 = params.R2;
	kf = Kalman_filter(N,x0,tau,tau_additional,R1,R2);
	contr = StructuredController(N, qvec, r, tau,tau_additional,params.H,kf,b_vec,c_vec);
	D_c = contr.scale_dist(Dist);%needed due to the change of variables.
	sys = System_Discrete(N,g,sys_choice, lp,x0, 0);
elseif contr_set ==2
	D_c = Dist;
	contr = Pcontroller(N,tau,tau_additional,b_vec,c_vec,params.k);
	sys = System_Discrete(N, g,sys_choice, lp, x0, 0);
else
	sys = System_Discrete(N, g, sys_choice, ss(1), x0, 1); %use minimal realization, no filter
	q = ones(N,1);
	r = params.r;
	r_du = params.r_du;
	H = params.H;
	contr =  LQGcontroller(N, sys, q, r, r_du, H);
	
    %Filter dist as there is no filter in the system for the LQ controller
    D_f_state = zeros(length(lp.A),N);
    for t = 2:T
        for i = 1:N
            D_f_state(:,i) = lp.A*D_f_state(:,i)+lp.B*Dist(i,t);
            Dist(i,t) = lp.C*D_f_state(:,i)+lp.D*Dist(i,t);
        end
    end
	D_c = Dist;
end


for t = 1:T
    D_contr = D_c(:,t:end);

	if contr_set == 1 %Stuct or LQ -> Output feedback
		measurement = sys.get_levels();
		[contr, u]  = contr.controller_sample(measurement,D_contr);
    end
    if contr_set == 2
        measurement = sys.get_levels();
		[contr, u]  = contr.calculate_u(measurement,D_contr);
    end
 	if contr_set == 3 %LQ -> State feedback
 		state = sys.get_state();
 		u = contr.calculate_u(state,D_contr);
    end
    sys = sys.apply_input(u,Dist(:,t));
end
end



