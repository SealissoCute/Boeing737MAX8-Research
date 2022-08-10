clc;clear
format long
CruiseAltitude = 36000;
FieldAltitude = 5000;
Range = 2300;
[a,rho]=Standard_Atmosphere(FieldAltitude);
ft_s_to_kt = 0.592484;
WoverS = 0:10:200;
W_TO = 182200;
c = 0.0199;
d = 0.7531;
S_wet = 10^(c+d*log10(W_TO));
% From Table 3.4 Correlation Coefficients For Parasite Area Versus Wetted Area
cf_1 = 0.002; a_1 = -2.6990; b_1 = 1;
cf_2 = 0.003; a_2 = -2.5229; b_2 = 1;
cf_3 = 0.004; a_3 = -2.3979; b_3 = 1;
f_1 = 10^(a_1+b_1*log10(S_wet));
f_2 = 10^(a_2+b_2*log10(S_wet));
f_3 = 10^(a_3+b_3*log10(S_wet));
AR = 10;
e =0.85;
%% FAR25 TAKEOFF DISTANCE SIZING
figure()
hold on
for CL_max_TO = 1.6:0.2:2.2
    ToverW = (0.009640.*WoverS)/CL_max_TO;
    plot(WoverS,ToverW)
end
legend('1.6','1.8','2.0','2.2');
xlabel('W/S');
ylabel('T/W');
title('FAR25 TAKEOFF DISTANCE SIZING');
hold off

%% FAR25 LANDING DISTANCE SIZING
figure()
hold on
for CL_max_L = 1.8:0.2:2.8
    xline(5000/0.3/1.69/2*rho*CL_max_L/ft_s_to_kt^2)
end
hold off

%% FAR25 CLIMB RATE SIZING
% For Take-off Climb
CD0_1 = f_1/S_wet;
CD0_2 = f_2/S_wet;
CD0_3 = f_3/S_wet;
syms CL
CD = vpa(CD0_2 + CL^2/(pi*AR*e));
%% DIRECT CLIMB SIZING

%% CRUISE SPEED SIZING
   [a,rho]=Standard_Atmosphere(CruiseAltitude);
   CruiseSpeed_Mach = 0.79;
   CruiseSpeed = CruiseSpeed_Mach*a;
   WoverS = 0:10:200;
   C_D0 = 0.0184; %p.145&182 low speed,clean drag polar
   delta_C_D0 = 0.0001*2.5; % p.166 figure 3.32
   C_D0_modification = C_D0 + delta_C_D0;
   q_overline = 0.5*rho*CruiseSpeed^2;
   ToverW_cruise_reqd = C_D0_modification*q_overline./WoverS + WoverS./(q_overline*pi*AR*e);
   ToverW_TO = ToverW_cruise_reqd./0.23;
   plot(WoverS,ToverW_TO)
%%
function [a,rho]=Standard_Atmosphere(h)
% Standard Atmosphere (SI Units)
% [C,a,P,rho,g,mu]=Standard_Atmosphere(h)
%
% bibliography :
% [1] Yunus A. Cengel & John M.Cimbala°ßFLUID MECHANICS ...
% Fundamentals and  Applications°®, McGraw-Hill., p897.
% [2] John J. Bertin & Russell M. Cummings°ßAERODYNAMICS FOR ENGINEERS°®,
% 5th Edition, Pearson Education International., p21-p43.
% [3] WARREN F.PHILLIPS,°ßMECHANICS of FLIGHT°®, 2nd Edition, John Wiley.
% p10-p14.
% [4] John D.Anderson,"Modern Compressibe Flow", third Edition, McGraw-Hill
% ., p585-p613.
%
% input arguments:
%    h = Geometric altitude. (default : sea level)
%
% output arguments:
%    Rankine = The temperature in Rankine scale.
%    a = Speed of sound.
%    P = The standard atmosphere at h.
%    rho = Density.
%    g = Is the gravitational acceleration at height h above sea level.
%    mu = Coefficient of viscosity.
%
% example :
%    % Plot C v.s geometrix altitude and P v.s geometrix altitude.
%    >> [C,a,P,rho,g,mu]=Standard_Atmosphere(100:100:90000);
%    >> figure, subplot(1,2,1);
%    >> plot(C,100:100:90000);
%    >> set(gca,'XTick',-100:20:20,'YTick',0:10000:100000,...
%            'YTickLabel',0:10:100,'DataAspectRatio',[1 650 1]);
%    >> title('Standard Atmosphere'); grid on;
%    >> xlabel('Temperature (Celsius)'); ylabel('Geometrix Altitude (Km)');
%    >> subplot(1,2,2);
%    >> plot(P,100:100:90000);
%    >> set(gca,'XTick',0:50000:150000,'YTick',0:10000:100000,...
%            'XTickLabel',0:50:150,'YTickLabel',0:10:100,...
%            'DataAspectRatio',[1 .65 1]);
%    >> title('Standard Atmosphere'); grid on;
%    >> xlabel('Pressure (kPa)'); ylabel('Geometrix Altitude (Km)');

% Last change: 2022/07/25 14:00 pm

% Unit exchange
ft_to_m = 0.3048;       % ft to m
m_s_to_mph = 2.236936;  % m/s to mph
m_s_to_kt = 1.943844;   % m/s to kt
m_s_to_ft_s = 3.280840; % m/s to ft/s
ft_s_to_kt = 0.592484;  % ft/s to kt
kg_to_slug = 0.068522;
%
h = h*ft_to_m; % unit:m

% default values
if ~exist('h','var'), h = 0; end;  % sea level

% Gravitational acceleration
% go = Is the standard gravitational acceleration.
% re = Is the Earth's mean radius.
go = 9.806645;           % m/s^2
re = 6356766;            % m
g = go*(re./(re+h)).^2;  % m/s^2

% Geopotential Altitude
% Z = Geopotential altitude.
% Zi = Is the minimum geopotential altitude in the range.
Z = (re*h)./(re+h);                                        % m
Zi = [0 11000 20000 32000 47000 52000 61000 79000 90000];  % m
[Zig Zg] = meshgrid(Zi,Z);

% Temperature in Celsius
% B = The lapse rate (Temperature Gradient);in SI units [K/m].
% Bi = The lapse rate (Temperature Gradient) for the range;
% Ti = Inital Temperature (absolute) inthe range ; K
% To = Is the sea_level temperature (absolute).
% T = Is temperature in K.
Bi = [-0.0065,0,0.001,0.0028,0,-0.002,-0.004,0];
Ti = [288.150,216.650,216.650,228.650,270.650,270.650,252.650,180.650];
Big = meshgrid(Bi,Z);
Tig = meshgrid(Ti,Z);
B = Big(Zig(:,1:8) <= Zg(:,1:8) &  Zg(:,2:9) < Zig(:,2:9));     % K/m
To = Tig(Zig(:,1:8) <= Zg(:,1:8) &  Zg(:,2:9) < Zig(:,2:9));
T = To'+B'.*...
    (Z-(Zig(Zig(:,1:8) <= Zg(:,1:8) &  Zg(:,2:9) < Zig(:,2:9)))');
Rankine = T*1.8; % unit:Rankine
% standard atmosphere pressure
% R = The gas constant for air.
% The gas constant for air in SI units [(N*m)/(kg*K)].
R = 287.0528;
[~,n] = max(...
    double(Zig(:,1:7) <= Zg(:,1:7) &  Zg(:,2:8) < Zig(:,2:8)),[],2);
N = 0;
P = zeros(size(Z));
for n = n',
    N = N+1;
    if B(N)==0,
        P(N) = Pressure(n+1,go,R,Zi,Ti,Bi)*...
            exp((-go*(Z(N)-Zi(n)))/(R*Ti(n)));
    else
        P(N) = Pressure(n+1,go,R,Zi,Ti,Bi)*...
            ((T(N))/Ti(n))^(-go/(R*Bi(n)));
    end
end

% Recursion function.
    function Pi = Pressure(n,go,R,Zi,Ti,Bi)
        n = n-1;
        if (n > 1)
            if Bi(n-1)==0,
                Pi = Pressure(n,go,R,Zi,Ti,Bi)*...
                    exp((-go*(Zi(n-1+1)-Zi(n-1)))/(R*Ti(n-1)));
            else
                Pi = Pressure(n,go,R,Zi,Ti,Bi)*...
                    ((Ti(n-1)+Bi(n-1)*(Zi(n-1+1)-Zi(n-1)))/...
                    Ti(n-1))^(-go/(R*Bi(n-1)));
            end
        else
            % Standard atmosphere pressure at sea level.
            Pi = 1.01325*10^5;
        end
    end

% Density
rho = P./(R.*T);
rho = rho*kg_to_slug*ft_to_m^3;
% Coefficient of viscosity
if nargout > 5
    % Viscosity in SI units [kg/(s*m)].
    mu = 1.458*(10^-6)*((T.^1.5)./(T+110.4));
end

% molecular energy.
% Cv = Constant volime ; Cv = e(internal energy)/T.
% Cp = Constant pressure ; Cp = h(enthalpy)/T.
% e = Internal energy.
% e_tr = Translational energy.
% e_rot = Rotational energy.
% e_vib = Vibrational energy.
e_tr = (3/2).*R.*T;
e_rot = R.*T;
e_vib = (1/2).*R.*T;
if ( T >= 600 ) % when the air temperature reaches 600K or higher .
    e = e_tr+e_rot+e_vib;
    Cv = e./T;
    Cp = (e+R.*T)./T;
else
    e = e_tr+e_rot;
    Cv = e./T;
    Cp = (e+R.*T)./T;
end

% gamma = Define Cp(constant pressure)/Cv(constant volime).
gamma = Cp./Cv;

% Speed of sound .
a = sqrt(gamma.*R.*T);
a = a*m_s_to_ft_s; % unit:ft/s
end