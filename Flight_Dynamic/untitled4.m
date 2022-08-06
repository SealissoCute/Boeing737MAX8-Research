%
twopi = 4. *asin(1.)
CM elevator = -0.9823 3
Mdelta = CM elevator*Q*S wing*c mac/Iy
% redefine the values of
Moment_w =-0.5046; Moment_q = -0.2971; Moment_wot = -0.0210;
%
Malpha = Moment_w*v
Mq = Moment_q*v
Ma_dot = Moment_wdot
A_short = [0 1; Malpha (Mq+Ma dot)]
B = [0; -Mdelta]
C = [1 0;0 1 ]
D = [0;0]
lambdal= poly(A_short)
[u vy] = eig(A_short)
%
t_half = 0.69/abs(real(vy(1,1)))
period = twopi/imag(vy(1,1))
ss(A_short,B,C,D)
sys_ss = ss(A_short, B, C, D)
%
% the transfer function alpha/delta
%
sys_tf0l=tf(sys_ss)
%
% figure(1)
% title('step response'):
% step(sys tf01,35)
% SS=stepinfo(systf01,35)
% grid on
%
