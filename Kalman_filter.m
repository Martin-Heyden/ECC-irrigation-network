classdef Kalman_filter
	%Estimates all first order states (including transportation states)
	% (but only the levels are updated based on measurements)
    properties
		L
		N
		priori
        sigma
		tau_bar
		A
		B
		state_vec
	end
	
    methods
		function obj = Kalman_filter(N,x0,tau,tau_bar,R1,R2)
			obj.N = N;
			L = zeros(N,1);
			obj.sigma = cumsum([0,tau]);
			obj.tau_bar = tau_bar;
			[A,B] = generate_state_space(tau, tau_bar);
			obj.A = A;
			obj.B = B;
			for i = 1:N
				%R1 statenoise, R2 measurement noise 
				[~,~,L_tr] = dare(1,1,R1,R2);
				L(i) = L_tr; %Scalar
			end
			obj.L = L;
			obj.priori = [x0;zeros(length(A)-N,1)];
		end
			
		%dist is assumed to be d[t-tau_bar]
		function [obj, new_est] = update_estimate(obj, y,input,dist)
			posteriori = obj.priori;
			posteriori(1:obj.N) = obj.priori(1:obj.N) + obj.L.*(y-obj.priori(1:obj.N)); % Local
			obj.priori = obj.A*posteriori+obj.B*input + [dist;zeros(obj.sigma(end)+obj.N*obj.tau_bar,1)];
			new_est = obj.priori;
			obj.state_vec = [obj.state_vec new_est(1:obj.N)];


		end


    end
end