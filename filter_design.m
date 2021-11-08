%Filter design (and fig 3)


z = tf('z',60)
w1 = {1e-4,1e-1}
w1f = {1e-4,1.7e-2}
%%% pool 9 %%%
ci = [0.137 -0.155 0.053];
cip1 = [-0.190 0.333 -0.175];%c_{i+1}
alfa1 = 0.978;
alfa2 = 0.468;
tau_i = 3;
g9 = [ci(1)*z^(-tau_i) + ci(2)*z^(-tau_i-1)+ci(3)*z^(-tau_i-2), ...
	cip1(1) + cip1(2)*z^(-1)+cip1(3)*z^(-2)]...
	/(z-1-alfa1*(1-2*z^(-1)+z^(-2))-alfa2*(1-z^(-1)))

%%% pool 10 %%%
w2 = {1e-4,1e-1}
w2f = {1e-4,1.2e-2}
ci = [0.134 -0.244 0.114];
cip1 = [-0.101 0.185 -0.087];%c_{i+1}
alfa1 = 0.314;
alfa2 = 0.814;
tau_i = 16;
g10 = [ci(1)*z^(-tau_i) + ci(2)*z^(-tau_i-1)+ci(3)*z^(-tau_i-2), ...
	cip1(1) + cip1(2)*z^(-1)+cip1(3)*z^(-2)]...
	/(z-1-alfa1*(1-2*z^(-1)+z^(-2))-alfa2*(1-z^(-1)))



%%% Low pass Filter Design %%%
%angular frequency 0.003 rad/s, sample rate 1 (/minute)
fc = 0.003/2/pi;
fs = 1/60
[A,B,C,D] = butter(3,fc/(fs/2));
lp = ss(A,B,C,D,60)
figure(2)
bode(lp)


%%% Plotting %%%
subplot(1,2,1)
hold off
[mag1p, phase, w1p] = bode(g9,w1)
loglog(w1p,squeeze(mag1p(1,1,:)),'linewidth',2)
hold on
[mag1fp, phase, w1fp] = bode(g9*lp,w1f)
loglog(w1fp,squeeze(mag1fp(1,1,:)),'linewidth',2)
legend('Third Order','Third Order + lowpass')

subplot(1,2,2)
hold off
[mag2p, phase, w2p] = bode(g10,w2)
loglog(w2p,squeeze(mag2p(1,1,:)),'linewidth',2)
hold on
[mag2fp, phase, w2fp] = bode(g10*lp,w2f)
loglog(w2fp,squeeze(mag2fp(1,1,:)),'linewidth',2)
legend('Third Order','Third Order + lowpass')
