classdef Pcontroller
    properties
        N
        K
        k_ff
        tau
		b_vec
		c_vec
		uprev
    end
    methods
        function  [obj, u] = calculate_u(obj,x,d) %x is given as volume
            u = zeros(obj.N,1);

            
           if ~isempty(d)
               if size(d,2)< 1+max(obj.tau) %Pad with zero if dist is too small
                        d = [d, zeros(obj.N,max(obj.tau)+1-size(d,2))];
               end
               for i = obj.N:-1:1
                   %feedforward due to offtake disturbance
                   u(i) = u(i) - obj.k_ff*obj.c_vec(i)/obj.b_vec(i)*d(i,1+obj.tau(i));
               end
		   end
		   for i = 1:obj.N
			   u(i) = u(i)-obj.K(i)*x(i);
			   if i ~= 1
				   u(i) = u(i) + obj.k_ff*obj.uprev(i-1)*obj.c_vec(i)/obj.b_vec(i); 
			   end
		   end
		   obj.uprev = u;
        end
        function obj = Pcontroller(N, tau,tau_additional, b_vec,c_vec,k,k_ff) 
			if nargin == 5
				k = 4;
                k_ff = 1
            elseif nargin == 6
                 k_ff = 1
            end
            obj.N = N;
            obj.tau = tau;
			obj.b_vec = b_vec;
			obj.c_vec = c_vec;
            for i = 1:N
               obj.K(i) = pi/2/(tau(i)+tau_additional)/b_vec(i)/k;
            end
            obj.k_ff = k_ff;
			obj.uprev = zeros(N,1);
		end
    end
end