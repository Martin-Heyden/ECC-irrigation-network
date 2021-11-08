%Fig 8 left
r_range_struct = [ 1 4 8 16 32 64 124, 248]
r_range_LQ = [0.1 0.2 0.5 1 2 4 6 8 16]
N = 10;

LQ_state_cost = zeros(1,length(r_range_LQ));
LQ_in_cost = zeros(1,length(r_range_LQ));

struc_state_cost = zeros(1,length(r_range_struct));
struc_in_cost = zeros(1,length(r_range_struct));

    dist_node = 5;
    sys_choice = ones(1,N);
    T = 500;

    D = zeros(N,T);
    D(dist_node,100:300) = -0.05;
for i = 1:length(r_range_struct)
    contr_set = 1;
    params.r = r_range_struct(i);
    params.R1 = 1;
    params.R2 = 100;
    params.H = 200;
    sys_S = simulate_system(N,sys_choice,contr_set,params,D,T);
    struc_state_cost(i) = sum_of_squares(sys_S.level_trajectory);
    struc_in_cost(i) =   sum_of_squares(sys_S.input_trajectory);
end
for i = 1:length(r_range_LQ)
    contr_set =3;
    params_LQ.r = r_range_LQ(i)*ones(N,1);
    params_LQ.r_du = 0.0*ones(1,N);
    params_LQ.H = 200;		
    sys_LQ = simulate_system(N,sys_choice,contr_set,params_LQ,D,T);
    LQ_state_cost(i) = sum_of_squares(sys_LQ.level_trajectory);
    LQ_in_cost(i) =   sum_of_squares(sys_LQ.input_trajectory);
end

figure(1)
hold off
plot(struc_state_cost,struc_in_cost,'d');
hold on
plot(LQ_state_cost,LQ_in_cost,'d');
legend('struct','LQ');
xlabel('State Deviation');
ylabel('Input Deviation');


