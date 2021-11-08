classdef System_Discrete
    properties
		N
		State_vec %Cell for all the states
		A
		B
		C
		level_trajectory;
		input_trajectory;
		filter_state_u
		filter_state_D
		LP
		sys_choice
		nbr_of_states;
	end
	
    methods
		%N: number of systems
		%sys: The dynamics of each system (discrete time transfer function)
		%LP: Low pass filter (discrete time state space)
		function obj = System_Discrete(N, sys_vec, sys_choice, LP, x0, min_real)
			obj.N = N;
			obj.A = cell(length(sys_vec),1);
			obj.B = cell(length(sys_vec),1);
			obj.C = cell(length(sys_vec),1);
			obj.sys_choice = sys_choice;
			for i = 1:length(sys_vec)
				if min_real
					sys_ss = ss(sys_vec{i},'min');
				else
					sys_ss = ss(sys_vec{i});
				end
				obj.A{i} = sys_ss.A;
				obj.B{i} = sys_ss.B;
				obj.C{i} = sys_ss.C;
			end
			obj.State_vec = cell(N,1);
			obj.nbr_of_states = 0;
			for i = 1:N
				n = sys_choice(i);
				%Solving for x such that Cx = x0 and Ax = x
				x = [obj.C{n};eye(size(obj.A{n}))-obj.A{n}]\[x0(i);zeros(length(obj.A{n}),1)];
				obj.State_vec{i} = x;
				obj.nbr_of_states = obj.nbr_of_states+length(x);
			end
			obj.level_trajectory = x0;
			obj.input_trajectory = [];
			obj.LP = LP;
			%%% Filter states %%%
			obj.filter_state_u = cell(length(LP.A),1);
			obj.filter_state_D = cell(length(LP.A),1);
			for i = 1:N
				obj.filter_state_u{i} = zeros(length(LP.A),1);
				obj.filter_state_D{i} = zeros(length(LP.A),1);
			end
			
		end
		
		function obj = apply_input(obj,u,d)
			%%% filter inputs %%%
			u_filt = zeros(size(u));
			d_filt = zeros(size(d));
			for i = 1:obj.N
				obj.filter_state_u{i} = obj.LP.A*obj.filter_state_u{i} + obj.LP.B*(u(i));
				u_filt(i) = obj.LP.C*obj.filter_state_u{i} + obj.LP.D*(u(i));
				obj.filter_state_D{i} = obj.LP.A*obj.filter_state_D{i} + obj.LP.B*(d(i));
				d_filt(i) = obj.LP.C*obj.filter_state_D{i} + obj.LP.D*(d(i));
			end
			for i = 1:obj.N
				n = obj.sys_choice(i);
				inflow = u_filt(i);
				if i ~= 1 %Not most downstream node
					outflow = -d_filt(i)+ u_filt(i-1); 
				else %Most downstream node
					outflow = -d_filt(i);
				end
				obj.State_vec{i} = obj.A{n} * obj.State_vec{i} + obj.B{n}*[inflow;outflow];
			end
			obj.level_trajectory = [obj.level_trajectory obj.get_levels()];
			obj.input_trajectory = [obj.input_trajectory u_filt];
		end
		
		function x = get_levels(obj)
			x = zeros(obj.N,1);
			for i = 1:obj.N
				x(i) = obj.C{obj.sys_choice(i)}*obj.State_vec{i};
			end
		end
		
		function x = get_state(obj)
			x = zeros(obj.nbr_of_states,1);
			j = 1;
			for i = 1:obj.N
				state = obj.State_vec{i};
				l = length(state);
				x(j:j+l-1) = state;
				j = j+l;
			end
		end

    end
end