%used for Fig 6 Left

N_range = 2:1:10;
nbr_sim = 1;


LQ_cost = zeros(size(N_range));
P_cost = zeros(size(N_range));
struc_cost = zeros(size(N_range));

k_range = [1 2 4 6 8];

for i = 1: length(N_range)
	N = N_range(i);
	for dist_node = N-1
		for j = 1:nbr_sim
			sys_choice = ones(1,N);
			T = 500;
			
			D = zeros(N,T);
			D(dist_node,100:300) = -0.5;
			contr_set = 1;
			params.r = 0.3;
			params.R1 = 1;
			params.R2 = 100;
			params.H = 300;
			sys_S = simulate_system(N,sys_choice,contr_set,params,D,T);
			struc_cost(i) = struc_cost(i)+sum_of_squares(sys_S.level_trajectory) + ...
                sum_of_squares(sys_S.input_trajectory(N,:))*params.r;
			
			contr_set =2;
            best_cost = inf;
            best_k = 0;
            
            for k = k_range
                params.k =k;
                sys_P = simulate_system(N,sys_choice,contr_set,params,D,T);
                cost = P_cost(i)+sum_of_squares(sys_P.level_trajectory) + ...
                    sum_of_squares(sys_P.input_trajectory(N,:))*params.r;
                if cost <best_cost
                    best_cost = cost;
                    best_k = k;
                end
            end
            P_cost(i) = best_cost;
            best_k

			contr_set =3;
			params_LQ.r = 0.00003*ones(N,1);
			params_LQ.r_du = 0.0*ones(1,N);
			params_LQ.r(end) = 0.3;
			params_LQ.H = 300;		
			sys_LQ = simulate_system(N,sys_choice,contr_set,params_LQ,D,T);
			LQ_cost(i) = LQ_cost(i)+sum_of_squares(sys_LQ.level_trajectory) + ...
                sum_of_squares(sys_LQ.input_trajectory(N,:))*params.r;
		end
	end
end

figure(1)
hold off
plot(N_range,struc_cost,'-d');
hold on
plot(N_range,P_cost,'-d');
plot(N_range,LQ_cost,'-d');
legend('struct','P','LQ');



