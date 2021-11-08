function [A,B,C,D] = generate_state_space(tau_vec, tau_bar)
%   First N states are the node levels, Then follows tau_vec(1)+tau_bar
%   states for the memory of u_1 in the order u_1(t-tau_i-tau_bar),
%   u_1(t-tau_i-tau_bar + 1)... u_1(t-1) etc. 

%used in kalman filter for keeping track of the necessary old inputs.

N = length(tau_vec);
nbr_states = N+sum(tau_vec)+N*tau_bar;
A = zeros(nbr_states,nbr_states);
B = zeros(nbr_states,N);
C = zeros(N,nbr_states);
C(1:N,1:N) = eye(N);
D = zeros(N,1)

%%% A matrix %%%
A(1:N,1:N) = eye(N);
current_state = N+1;
for i = 1:N %loop over inputs
	A(i,current_state) = 1; %inflow to node i
	for j = 1:tau_vec(i)+tau_bar-1
		A(current_state,current_state+1) = 1;
		current_state=current_state+1;
		if j == tau_vec(i) 
			if i <N
				A(i+1,current_state) = -1; %Outflow from node i
			end
		end
	end
	B(current_state,i) = 1; %Save u_i[t-1]
	current_state = current_state+1; %u_i[t-1]
	
end

end

