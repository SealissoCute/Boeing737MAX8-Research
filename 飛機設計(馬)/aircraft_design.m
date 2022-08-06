clc
clear
%%
% Parameters
ryo = 1.225;
Vs = 7.0;
CL_max = 1.2;
%
%%
% Wing loading
WoverS = 0.5*ryo*(Vs^2)*CL_max
% Maxmium Speed

%%
% Plot
x1 = [0:70];
y1 = zeros(1,71)+1/7;
x2 = zeros(1,2)+WoverS;
y2 = [0:1];
plot(x1,y1,x2,y2)
%
