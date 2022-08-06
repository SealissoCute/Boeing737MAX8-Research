global g L
g = 9.81; % acceleration of gravity, m/s^2
L = 1; % length of the pendulum string, m
t = 0:0.1:100;

y_ic = [15/360*2*pi 0];
[t,y] = ode45(@fun,t,y_ic);
figure(1)
plot(t,y(:,1))
xlabel('time( sec )') 
ylabel('delta theta')
grid on

function dxdt=fun(t,y)
global g L
dxdt = zeros(2,1);
dxdt(1) = y(2);
dxdt(2) = -(g/L)*sin(y(1));
end

