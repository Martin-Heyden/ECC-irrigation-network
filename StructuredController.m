classdef StructuredController
    properties
        K
        gamma_N
        H
        tau
		tau_bar
        N
        r
        g
		d_mem
		kf
		state
		b_hat
		c
    end
    methods
		function [dist] = scale_dist(obj,dist)
			for j = 1:size(dist,2) %First node
				dist(1,j) = dist(1,j)*obj.c(1);
			end
			for i = 2:size(dist,1)
				for j = 1:size(dist,2)
					dist(i,j) = dist(i,j)*obj.b_hat(i-1);
				end
			end
		end
		
		%dist is assumed to be scaled (so that it only need to be done
		%once)
		function [obj, u] = controller_sample(obj,measurement,dist)
			%scaling of measurement
			for i = 2:obj.N
				measurement(i) = measurement(i)*obj.b_hat(i-1)/obj.c(i);
			end

			%calculate input
			u = obj.calculate_u(obj.state,dist); %using the a priori estimate that was calculated last iteration
			%update kalman filter.
			[obj.kf, obj.state] = obj.kf.update_estimate(measurement,u,obj.d_mem(:,1));
			%update d_mem
			if isempty(dist)
				obj.d_mem = [obj.d_mem(:,2:end) zeros(obj.N,1)];
			else
				obj.d_mem = [obj.d_mem(:,2:end) dist(:,1)];
			end
			%scaling of outputs
			for i = 1:obj.N
				u(i) = u(i)/obj.b_hat(i);
			end
		end
		
        function u = calculate_u(obj,state,dist)
            %%% Effect from states %%%
            u = obj.K*state;
            %%% effect from disurbance %%%
            if ~isempty(dist) %  dist
                %Pad with zeros if needed
                q = size(dist,2);
				if q<(1+sum(obj.tau)+obj.H) %Pad with zero if dist is too small
					dist = [dist, zeros(obj.N,1+sum(obj.tau)+obj.H-q)];
				end
				dist(:,obj.H+1:end) = 0; %only use first H dist (for fair comparison with LQ.
                %Calculate shifted disturbances
                dist_horizon = [obj.tau(1:end-1),obj.H];
                sigma = [0 cumsum(dist_horizon)];
                Dsum = zeros(obj.N,1);%sum of D_i[t] for sigma_i <= t < sigma_{i+1}
                Dfirst = zeros(obj.N,1); % D_i[sigma_i]
                
				%Node 1 to N-1
                for i = 1:obj.N-1 
                    for j = sigma(i):(sigma(i+1)-1)
                        for k = 1:i
                            Dsum(i) =  Dsum(i)+dist(k,j+1-sigma(k));
                            if j == sigma(i)%disturbance about to arrive, D_i[sigma_i]
                                Dfirst(i) =  Dfirst(i)+dist(k,j+1-sigma(k));
                            end
                        end
                    end
				end     
				%%% Node N (which needs to be handled differently due to the
				%%% production.
                DN_weighted = 0;
                weight = 1;
                for j = 1:obj.H
                    if j > obj.tau(obj.N)+1 %j = 1 corresponds to t, so j = tau_N+2 corresponds to t = tau_N+1
                       weight = weight*obj.g;
                    end
                    for k = 1:obj.N
                        DN_weighted = DN_weighted+weight*dist(k,j+sigma(end-1)-sigma(k));
						if j == 1 %Disturbance about to arrive
							Dfirst(obj.N) = Dfirst(obj.N) + dist(k,j+sigma(end-1)-sigma(k));
						end
					end
				end
			else
				Dsum = zeros(obj.N,1);
                Dfirst = zeros(obj.N,1);
				DN_weighted = 0;
				dist = zeros(obj.N,1);
			end
			%The effect of disurbance is given by
			%-gamma_i/q_i(D_i + sum_{j<i}sum_{t = ...} D_i[t])
			%gamma_i/q_i are the elements of K
			for i = 1:obj.N-1
				%sum_{d=0->tau_i-1}D_i[t+sigma_i+d] +
				%    sum_{d=1->tau_bar}d_i[t-s]
				Dtemp = Dsum + sum(obj.d_mem,2);
				Dtemp(i) = Dtemp(i) + Dfirst(i+1); %Dfirst(i+1) corresponds to D_{i+1}[t+sigma_i]
				Dtemp(i+1) = sum(obj.d_mem(i+1,:));%no sum_{d=0->tau_i-1}D_i[t+sigma_i+d]
				u(i) = u(i) + obj.K(i,1:obj.N)*Dtemp + dist(i+1,1);
			end
			%%% Production %%%
			u(obj.N) = u(obj.N) +obj.K(obj.N,1)*(sum(Dsum(1:end-1))+ ...
				DN_weighted + sum(sum(obj.d_mem,2)));
 			for i = 1:obj.N-1 
 				u(i) = u(i)*1;
 			end
        end
        
        function obj = StructuredController(N, q_vec, r_prod, tau, tau_bar, H, kf, b_vec, c_vec)
           obj.N = N;
           obj.H = H;
           obj.tau = tau;
		   obj.tau_bar = tau_bar;
           
           sigma = cumsum([0, tau]);
           nbr_states = N+sum(tau)+N*tau_bar;
		   obj.d_mem = zeros(N,tau_bar);
		   obj.kf = kf;
		   obj.state = kf.priori;
		   
		   %%% Change of Variables %%%
		   b_hat_vec = b_vec;
		   for i = 2:N
			   b_hat_vec(i) = b_vec(i)/c_vec(i)*b_hat_vec(i-1);
		   end
		   obj.b_hat =b_hat_vec;
		   obj.c = c_vec;
		   
		   qvec = ones(N,1);
		   for i = 2:N
			   qvec(i) = (c_vec(i)/b_hat_vec(i-1))^2*q_vec(i); 
		   end
		   q_vec = qvec;
		   r_prod = r_prod/b_hat_vec(N)^2;
		   obj.r=r_prod;  
		   %%% Calculation of Control Parameters
           gamma = zeros(N,1);
           gamma(1) = q_vec(1);
           for i = 2:N
               gamma(i) = q_vec(i)*gamma(i-1)/(q_vec(i) + gamma(i-1));
           end
           obj.gamma_N = gamma(N);
          %%% Internal Transportation %%%
          obj.K = zeros(N,nbr_states);
          for k = 2:N
              obj.K(k-1,1:k-1) = -gamma(k)/q_vec(k); %z_{k-1}[t]
			  for i =1:k-1
				  %sum_{s=t_bar+1->tau_i+tau_bar} u_i[t-s]
				  first_state = N+1+sigma(i)+(i-1)*tau_bar;
				  last_state = first_state+(tau(i)-1);
				  obj.K(k-1,first_state:last_state) = -gamma(k)/q_vec(k); 
			  end
			  %sum{1->t_bar} u_{k-1}[t-s]
			  first_state = N+1+sigma((k-1)+1)+((k-1)-1)*tau_bar;
			  last_state = first_state+tau_bar-1;
			  obj.K(k-1,first_state:last_state)  = -1;
              obj.K(k-1,k) = (1-gamma(k)/q_vec(k)); %z_k[t]
			  %sum_{s=tau_k}^{tau_k+tau_bar} u_k[t-s]
			  first_state = N+1+sigma(k)+(k-1)*tau_bar;
			  last_state = first_state+tau_bar;
              obj.K(k-1,first_state:last_state) = (1-gamma(k)/q_vec(k));
		  end
		  %%% Production %%%
          [~,~,G] = dare(1,1,gamma(N),r_prod);
		  for i = 1:N-1
			  obj.K(N,i) = -G; %z_i[t]
			  %sum_{s=tau_bar+1 ->tau_bar + tau_i} u_i[t-s]
			  first_state =  N+1+sigma(i)+(i-1)*tau_bar;
			  last_state = first_state+(tau(i)-1);
			  obj.K(N,first_state:last_state) = -G;
		  end
          obj.K(N,N) = -G;  %z_N[t]
		  %sum_{s=1->tau_bar+tau_N} u_N[t-s]
		  obj.K(N,N+1+sigma(N)+(N-1)*tau_bar:end) = -G;
          
            Xr = obj.gamma_N/2+sqrt(obj.gamma_N*obj.r+obj.gamma_N^2/4);%same as dare(1,1,gamma(N),rho(N));
            X =  Xr-obj.gamma_N; 
            obj.g = X/(X+obj.gamma_N);
        end
    end
end