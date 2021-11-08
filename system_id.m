%System id and Fig 4

set = 2; %which input to plot for

%%% Low pass filter %%%
fc = 0.003/2/pi;
fs = 1/60
[A,B,C,D] = butter(3,fc/(fs/2));
lp = ss(A,B,C,D,60)
z = tf('z',60);


%%% Pool 9 %%%
%Third order system
ci = [0.137 -0.155 0.053];
cip1 = [-0.190 0.333 -0.175];%c_{i+1}
alfa1 = 0.978;
alfa2 = 0.468;
tau_i = 3;
g9 = [ci(1)*z^(-tau_i) + ci(2)*z^(-tau_i-1)+ci(3)*z^(-tau_i-2), ...
	cip1(1) + cip1(2)*z^(-1)+cip1(3)*z^(-2)]...
	/(z-1-alfa1*(1-2*z^(-1)+z^(-2))-alfa2*(1-z^(-1)));
sys_9 = g9*lp;
%First order coefficients
b9 = 0.069;
c9 = -0.063; 

%%% Pool 10 %%%
%Third order system
ci = [0.134 -0.244 0.114];
cip1 = [-0.101 0.185 -0.087];%c_{i+1}
alfa1 = 0.314;
alfa2 = 0.814;
tau_i = 16;
g10 = [ci(1)*z^(-tau_i) + ci(2)*z^(-tau_i-1)+ci(3)*z^(-tau_i-2), ...
	cip1(1) + cip1(2)*z^(-1)+cip1(3)*z^(-2)]...
	/(z-1-alfa1*(1-2*z^(-1)+z^(-2))-alfa2*(1-z^(-1)));
sys_10 = g10*lp;
tau_10 = 14;
%First order coefficients
b10 = 0.0142*1.5;
c10 = -0.0156; 

%%% Simulation %%%

u = [ones(25,1);zeros(50,1);-ones(25,1);zeros(40,1)];
z_vec = zeros(size(u));
T = (0:(length(u)-1))*60;
tau_9 = 3;%0:3;
tau_add_cand = 1:20;
err_mat = zeros(2,size(tau_add_cand,2));
z = tf('z',60);
%loop over candidate delays, 
% Look at dynamics for outflow and find joint tau_bar
for j = 1:length(tau_add_cand)
    tau_additional = tau_add_cand(j);
    sys_first_9 = [b9*z^(-tau_9), c9]/(z-1)*z^(-tau_additional);
    %simulate response
    y_third_9_1 = lsim(sys_9,[z_vec, u],T);
    y_first_9_1 = lsim(sys_first_9,[z_vec, u],T);
    %calculate least square error
    err_mat(1,j) = sum((y_third_9_1-y_first_9_1).^2);
    sys_first_10 = [b10*z^(-tau_10), c10]/(z-1)*z^(-tau_additional);
    %simulate response
    y_third_10_1 = lsim(sys_10,[z_vec, u],T);
    y_first_10_1 = lsim(sys_first_10,[z_vec, u],T);
    %calculate least square error
    err_mat(2,j) = sum((y_third_10_1-y_first_10_1).^2);
end

%%% Finding best tau_add %%%
%normalize each system
for i = 1:2
    err_mat(i,:) = err_mat(i,:)/max(err_mat(i,:));
end
total_error = sum(err_mat);
[~,i] = min(total_error);
best_tau_add = tau_add_cand(i)

%%% finding tau_i %%%
tau_cand_9 = 0:6;
err_9 = zeros(1,length(tau_cand_9));
for i = 1:length(tau_cand_9)
    tau_9 = tau_cand_9(i); 
    sys_first_9 = [b9*z^(-tau_9), c9]/(z-1)*z^(-best_tau_add);
    y_third_9_2= lsim(sys_9,[u, z_vec],T);
    y_first_9_2 = lsim(sys_first_9,[u, z_vec],T);
    err_9(i) = sum((y_third_9_2-y_first_9_2).^2);    
end
[~,i] = min(err_9);
best_tau_9 = tau_cand_9(i)

tau_cand_10 = 5:18;
err_10 = zeros(1,length(tau_cand_10));
for i = 1:length(tau_cand_10)
    tau_10 = tau_cand_10(i);
    sys_first_10 = [b10*z^(-tau_10), c10]/(z-1)*z^(-best_tau_add);
    y_third_10_2= lsim(sys_10,[u, z_vec],T);
    y_first_10_2 = lsim(sys_first_10,[u, z_vec],T);
    err_10(i) = sum((y_third_10_2-y_first_10_2).^2);
end
[~,i] = min(err_10);
best_tau_10 = tau_cand_10(i)


%%% Plotting %%%
if set == 1
	outflow  = u';
	inflow = z_vec';
else
	outflow  = z_vec';
	inflow = u';
end
yg1_9 = lsim(g9,[inflow;outflow]);
yg1lp_9  = lsim(g9*lp,[inflow;outflow]);
sys_first = [b9*z^(-best_tau_9), c9]/(z-1)*z^(-best_tau_add);
yfirst_9 = lsim(sys_first,[inflow;outflow]);
subplot(1,2,1)
hold off
plot(yg1_9,'b','LineWidth',0.1)
hold on
plot(yg1lp_9,'r','LineWidth',2)
plot([yfirst_9],'k','LineWidth',2)
legend( 'Third order', 'Third Order low pass', 'First order')
title('Pool 9');


yg1_10 = lsim(g10,[inflow;outflow]);
yg1lp_10  = lsim(g10*lp,[inflow;outflow]);
sys_first = [b10*z^(-best_tau_10), c10]/(z-1)*z^(-best_tau_add);
yfirst_10 = lsim(sys_first,[inflow;outflow]);
subplot(1,2,2)
hold off
plot(yg1_10,'b','LineWidth',0.1)
hold on
plot(yg1lp_10,'r','LineWidth',2)
plot([yfirst_10],'k','LineWidth',2)
legend( 'Third order', 'Third Order low pass', 'First order')
title('Pool 10');

T = 0:(length(yg1_10)-1);