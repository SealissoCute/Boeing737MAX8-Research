function [] = call_pend()
tspan=0:0.01:100; ; % set time interval
y0=[15/57.3,0]; % set initial conditions
[t,y]=ode45(@pend,tspan,y0);
plot(t,y(:,:))
title('Pendulum')
xlabel('Time')
legend('theta ( t )','omega ( t )')
function dydt = pend(t,y)
G=9.8; L=2; % set constants
y1=y(1); % get y1
y2=y(2); % get y2
dydt = [y2 ; -G/L*sin(y1);];
end
end