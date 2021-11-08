classdef LQGcontroller
    properties
		K
		Kv
		H
		A
		B
		S
		N
		B_d
		last_u;
		sys_choice;
    end
    methods
        function u = calculate_u(obj,state,d)
            u = obj.K*[state;obj.last_u;zeros(obj.N,1)]; %Adding delta u states
            if ~isempty(d)
                if size(d,2)< obj.H+1 %Pad with zero if dist is too small
                    d = [d, zeros(size(d,1),obj.H+1-size(d,2))];
                end
                D = zeros(length(obj.K),1);
                for i = obj.H:-1:1
					dyn_d = [];
					for j = 1:obj.N
						n = obj.sys_choice(j);
						dyn_d = [dyn_d;-d(j,i)*obj.B_d{n}(:,2)];
					end
					dyn_d = [dyn_d;zeros(2*obj.N,1)]; %Adding delta u states.
					D = (obj.A+obj.B*obj.K)'*D - obj.S*dyn_d;
                end
                u = u+obj.Kv*D;              
			end
			obj.last_u = u;
        end
        
        function obj = LQGcontroller(N, sys, q, r , r_du, H)
			
			%%% Create state space for entire system %%%
			A = zeros(sys.nbr_of_states);
			B = zeros(sys.nbr_of_states,N);
            s_i = 1
			for i =1:N
				%Indices for the states corresponding to subsystem i
				n = sys.sys_choice(i);
				sub_sys_size = length(sys.A{n});
				e_i = s_i+sub_sys_size-1;
				A(s_i:e_i,s_i:e_i) = sys.A{n};
				if i ==1 %Only inflow
					B(s_i:e_i,i) = sys.B{n}(:,1);
				else %Inflow and outflow
					B(s_i:e_i,i-1) = sys.B{n}(:,2);
					B(s_i:e_i,i) = sys.B{n}(:,1);
                end
                s_i = s_i + sub_sys_size;
			end
			%%% Add states to penalize u_i[t]-u_i[t-1] %%%
			%introduce N states corresponding to z_i[t+1] = u_i[t]
			%and N states corresponding to dx_i[t+1] = u_i[t] - z_i[t]
			A_d = zeros(2*N,2*N);
			B_d = zeros(2*N,N);
			for i = 1:N
				A_d(N+i,i) = -1;
				B_d(i,i) = 1;
				B_d(N+i,i) = 1;
			end
			A = [A, zeros(sys.nbr_of_states,2*N);
				zeros(2*N,sys.nbr_of_states), A_d];
			B = [B;B_d];
				
			
			
			%%% Penalty matrices %%%
			Qext = zeros(sys.nbr_of_states);
            
            s_i = 1;
			for i = 1:N
				n = sys.sys_choice(i);
				sub_sys_size = length(sys.A{sys.sys_choice(i)});
				e_i = s_i+sub_sys_size-1;
				Qext(s_i:e_i,s_i:e_i) = sys.C{n}'*q(i)*sys.C{n}; %The output is the cost
                s_i = s_i + sub_sys_size;
			end
			Q_d = diag([zeros(1,N),r_du]);
			Qext = [Qext,zeros(size(Qext,1),size(Q_d,2));
				zeros(size(Qext,1),size(Q_d,2))',Q_d];
			R = diag(r);
			[S,~,K] = dare(A,B,Qext,R);
			obj.K = -K;
            %%% Feedforward %%%
			obj.Kv = (B'*S*B+R)\B';
            obj.H = H;
            obj.A = A;
            obj.B = B;
            obj.S = S;    
			obj.N = N;
			obj.B_d = sys.B %the effect of disturbances on the states
			
			obj.last_u = zeros(N,1);
			obj.sys_choice = sys.sys_choice;
        end
    end
end