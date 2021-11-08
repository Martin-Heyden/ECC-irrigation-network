%Time response plot

rng(123)

contr_set = 1;
N = 5;
x0 = 5*[-1 0 0 0 1]';
if  contr_set == 1
	params.r = 0.3;
	params.R1 = 1;
	params.R2 = 100;
	params.H = 300
elseif contr_set ==2
	params.k = 4;
    params.k_ff = 1;
else
	params.r = 0.0003*ones(N,1);
	params.r_du = 0.00*ones(1,N);
	params.r(end) = 0.3;
	params.H = 300;
end
sys_choice = [1 2 1 2 1]
T = 600;

D = zeros(N,T);
D(1,250:450) = -1/0.069;

sys = simulate_system(N,sys_choice,contr_set,params,D,T,x0)

%%% Visualization %%%
figure(contr_set)


clf
subplot(2,1,1)
plot(0:1:T,sys.level_trajectory,'LineWidth', 2)
hold on
title('Node Levels')
legend('Node 1','Node 2', 'Node 3', 'Node 4', 'Node 5')
grid on
subplot(2,1,2)
plot(0:(T-1),sys.input_trajectory, 'LineWidth', 2)
hold on
title('Node Inputs')
grid on
sum_of_squares(sys.level_trajectory)*1+sum_of_squares(sys.input_trajectory(N,:))*0.3





