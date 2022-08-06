% M-file for Midterm Exam of Advanced Dynamics Of Flight, fall 2017
%
%% Calculate the mean aerodynamic chord
%
% Wing
% Chord root 
  Croot_w = 0.3055;
% Chord tip 
  Ctip_w = 0.3055;
% Taper ratio
  lambda_w = Ctip_w /Croot_w;
% Mean Aerodynamic Chord
  mac_w =(2/3)*((1+lambda_w+lambda_w^2)/(1+lambda_w))*Croot_w
% Tail
%Chord root 
  Croot_t = 0.1833;
% Chord tip 
  Ctip_t = 0.1833;
% Taper ratio
  lambda_t = Ctip_t /Croot_t;
% Mean Aerodynamic Chord
  mac_t =(2/3)*((1+lambda_t+lambda_t^2)/(1+lambda_t))*Croot_t
%
%% (a) Determine the mounting angles of wing and tail at design point
%
rad2deg = 180/pi;
% initial setting
mass = 8.0; % mass of UAV (kg) 
S_w = 1.6804; % reference wing area (m2)
b_w = 5.5; % wing span (m)
AR_w = b_w^2/S_w; % aspect ratio of wing
h = 100; % m, cruising height
[C,a,P,rho,g,muc]=Standard_Atmosphere(h); % subroutine for Standard Atmosphere
v_0 = 9.04; % level flight velocity (m/s)
CL_alpha_w = 0.09823*rad2deg; % lift coefficient slope of wing
alpha_CL0_w =  -6.45739; % Zero lift angle of attack (degree)
CL_trim = (mass*g)/(0.5*rho*v_0^2*S_w); % lift coefficient
% long ppt pg22
i_w = (CL_trim/CL_alpha_w)*rad2deg + alpha_CL0_w % mounted angle ( degree ) Wing
% long ppt pg23
epsilon_d0 = (0.6/AR_w)*CL_trim*rad2deg ; % downwash angle (degree)
i_t = epsilon_d0 % mounted angle ( degree ) Tail
%
%% (b)The location of the aerodynamic center of the wing relative to the center of gravity
%
Cm_ac_w = -0.0114698; % moment coefficient of aerodynamic center of Wing
l_w = Cm_ac_w/CL_trim*mac_w % c.g. to wing location, from long ppt pg 38 
%
%% (c)The pitch stability derivative for the wing-tail combination
%
CL_alpha_t = 0.06311*rad2deg; % lift coefficient slope of Tail
eta_t = 1; % tail efficiency
S_t = 0.16804;% reference tail area (m2)
l_t = 1.1778+l_w; % c.g. to tail location 
epsilon_d_alpha = (0.6/AR_w)*CL_alpha_w % Differentiating '' epsilon_d0 '' with respect to angle of attack
Cm_alpha = -(l_w/mac_w)*CL_alpha_w-(l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha) % pitch stability derivative 
CL_alpha = CL_alpha_w + (S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha);
SM = -Cm_alpha/CL_alpha*100
%
%% (d)  the trimmed angle of attack and elevator deflection for airspeed of 8~10 m/sec.
%
epsilon_e = 0.4; % elevator efficiency
Cm_ac_delta = -0.01174*rad2deg;
CL_0 = CL_trim;
Cm_0 = 0;
v_1 = 8:0.01:10;
CL_trim_1 = (mass.*g)./(0.5.*rho.*v_1.^2.*S_w);
CL_delta = (S_t./S_w).*eta_t.*CL_alpha_t.*epsilon_e;
Cm_delta = (mac_t./mac_w).*(S_t./S_w).*eta_t.*Cm_ac_delta-(l_t./mac_w).*(S_t./S_w).*eta_t.*CL_alpha_t.*epsilon_e;
d0 = CL_alpha.*Cm_delta - CL_delta.*Cm_alpha;
alpha_trim =((CL_trim_1 - CL_0).*Cm_delta + CL_delta.*Cm_0)./d0.*rad2deg;
delta_trim = (-(CL_trim_1 - CL_0).*Cm_alpha - CL_alpha.*Cm_0)./d0.*rad2deg;
% plot solution
figure(1)
plot(v_1,alpha_trim)
hold on
plot(v_1,delta_trim,'color','r')
hold off,grid on
axis([8 10 -20 20])
title('\alpha_{trim} versus \delta_{trim} ')
xlabel('Velocity ( m/s )'),ylabel('Degree')
legend('\alpha_{trim}','\delta_{trim}')
%
% (e) If the static margin is required to be 10%, what is the length of lt.
% use SM find lnp/mac
% use lnp/mac*CL_alpha find neigative CM_alpha
% sub into long ppt pg27 CM_alpha and find new lt
l_np_over_mac_w = 0.1;
n1 = l_np_over_mac_w*CL_alpha-(l_w/mac_w)*CL_alpha_w
d1 = (S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha)
l_t_over_mac_w = n1/d1
l_t_new = l_t_over_mac_w*mac_w
%
%
%% (f)  the trimmed angle of attack and elevator deflection with given SM of 10% for airspeed of 8~10 m/sec.
%
epsilon_e = 0.4; % elevator efficiency
Cm_ac_delta = -0.01174*rad2deg;
CL_0 = CL_trim;
Cm_0 = 0;
l_t = l_t_new; replace lt with lt_new
v_1 = 8:0.01:10;
CL_trim_1 = (mass.*g)./(0.5.*rho.*v_1.^2.*S_w);
CL_delta = (S_t./S_w).*eta_t.*CL_alpha_t.*epsilon_e;
Cm_delta = (mac_t./mac_w).*(S_t./S_w).*eta_t.*Cm_ac_delta-(l_t./mac_w).*(S_t./S_w).*eta_t.*CL_alpha_t.*epsilon_e;
Cm_alpha = -(l_w/mac_w)*CL_alpha_w-(l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha) % pitch stability derivative 
d0 = CL_alpha.*Cm_delta - CL_delta.*Cm_alpha;
alpha_trim_1 =((CL_trim_1 - CL_0).*Cm_delta + CL_delta.*Cm_0)./d0.*rad2deg;
delta_trim_1 = (-(CL_trim_1 - CL_0).*Cm_alpha - CL_alpha.*Cm_0)./d0.*rad2deg;
% plot solution
figure(2)
plot(v_1,alpha_trim_1,'b')
hold on
plot(v_1,delta_trim_1,'color','r')
hold off,grid on
axis([8 10 -20 20])
title('\alpha_{trim} versus \delta_{trim} at SM = 10%')
xlabel('Velocity ( m/s )'),ylabel('Degree')
legend('\alpha_{trim}','\delta_{trim}')
%
% compare the results
figure(3)
plot(v_1,alpha_trim)
hold on
plot(v_1,delta_trim,'color','r')
hold on
plot(v_1,alpha_trim_1,'--')
hold on
plot(v_1,delta_trim_1,'--','color','r')
hold off,grid on
%
function [C,a,P,rho,g,mu]=Standard_Atmosphere(h)
% =========================== Separator ==================================
% Graduate Student : Boming, Lin
% Department of Aerospace Engineering
% Tamkang University
% E-mail : 601430092@s01.tku.edu.tw
% Last change: 2013/07/24 12:08 am
% =========================== Separator ==================================
% Standard Atmosphere, SI Units, Revised for Version 1.2.
%     The first input arguments are required;
%     the another have default values.
% =========================== Separator ==================================
% Symbol:
% a = Speed of sound.
% B = The lapse rate (Temperature Gradient);in SI units [K/m].
% Bi = The lapse rate (Temperature Gradient) for the range;
% C = The temperature in Celsius.
% Cv = Constant volime ; Cv = e(internal energy)/T.
% Cp = Constant pressure ; Cp = h(enthalpy)/T.
% e = Internal energy.
% e_tr = Translational energy.
% e_rot = Rotational energy.
% e_vib = Vibrational energy.
% g = Is the gravitational acceleration at height h above sea level.
% go = Is the standard gravitational acceleration.
% gamma = Define Cp(constant pressure)/Cv(constant volime).
% h = Geometric altitude.
% mu = Coefficient of viscosity.
% P = The standard atmosphere at h.
% re = Is the Earth's mean radius.
% R = The gas constant for air.
% T = Is temperature in K.
% Ti = Inital Temperature (absolute) inthe range ; K
% To = Is the sea_level temperature (absolute).
% rho = Density.
% Z = Geopotential altitude.
% Zi = Is the minimum geopotential altitude in the range.
% =========================== Separator ==================================
% ------------------- Gravitational acceleration -------------------------
go=9.806645;     % m/s^2
re=6356766;      % m
g=go*(re/(re+h))^2;
% ---------------------- Geopotential Altitude ---------------------------
% Geopotential Altitude :
Z=(re*h)/(re+h); % m
Zi=[0 11000 20000 32000 47000 52000 61000 79000 90000];  % m
% ---------------------- Temperature in Celsius --------------------------
% in SI units [K/m].
Bi=[-0.0065,0,0.001,0.0028,0,-0.002,-0.004,0];
Ti=[288.150,216.650,216.650,228.650,270.650,270.650,252.650,180.650];
N=0;
for n=1:8
    N=N+1;
    if Zi(n)<=Z && Z<Zi(n+1)
        B=Bi(n);   % K/m 
        To=Ti(n);  % K
        break
    end
end
T=To+B*(Z-Zi(N));
C=T-273.15;
% ----------------------- Standard Atmosphere ----------------------------
% The gas constant for air in SI units [(N*m)/(kg*K)].
R=287.0528;
% Density in SI units [kg/m^3].
Pi=zeros(1,N);
Pi(1)=1.01325*10^5;
for i=1:N-1
    if Bi(i)==0
        Pi(i+1)=Pi(i)*exp((-go*(Zi(i+1)-Zi(i)))/(R*Ti(i)));
    else
        Pi(i+1)=...
            Pi(i)*((Ti(i)+Bi(i)*(Zi(i+1)-Zi(i)))/Ti(i))^(-go/(R*Bi(i)));
    end
end
P=Pi(N)*((To+B*(Z-Zi(N)))/To)^(-go/(R*B));
% ----------------------------- Density ----------------------------------
rho=P/(R*T);
% -------------------- Coefficient of viscosity---------------------------
if nargout > 5
    % Viscosity in SI units [kg/(s*m)].
    mu=1.458*(10^-6)*((T^1.5)/(T+110.4));
end
% -------------------------- Speed of sound ------------------------------
% molecular energy.
e_tr=(3/2)*R*T;
e_rot=R*T;
e_vib=(1/2)*R*T;
if ( T >= 600 ) % when the air temperature reaches 600K or higher .
    e=e_tr+e_rot+e_vib;
    Cv=e/T;
    Cp=(e+R*T)/T;
else
    e=e_tr+e_rot;
    Cv=e/T;
    Cp=(e+R*T)/T;
end
gamma=Cp/Cv;
% Speed of sound .
a=sqrt(gamma*R*T);
% =========================== Separator ==================================
% bibliography :
% [1] Yunus A. Cengel & John M.Cimbala¢X?FLUID MECHANICS ...
% Fundamentals and  Applications¢X?, McGraw-Hill., p897.
% [2] John J. Bertin & Russell M. Cummings¢X?AERODYNAMICS FOR ENGINEERS¢X?,
% 5th Edition, Pearson Education International., p21-p43.
% [3] WARREN F.PHILLIPS,¢X?MECHANICS of FLIGHT¢X?, 2nd Edition, John Wiley.
% p10-p14.
% [4] John D.Anderson,"Modern Compressibe Flow", third Edition, McGraw-Hill
% ., p585-p613.
% =========================== Separator ==================================
end