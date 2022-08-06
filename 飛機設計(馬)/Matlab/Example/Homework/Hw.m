% m-file of the first exam, 2017/08/07
%
clear
clc
rad2deg = 90./asin(1);
altitude = 0;
[C,a,P,rho,g,mu]=Standard_Atmosphere(altitude);
%
% Problem 1
%
% NACA 747A415
% from the airfoil data, 
  alpha_CL0_w = -2.0
%  
% at alpha = -8 deg Cl = -0.7, Cm_c/4 = -0.025
%    alpha = +8 deg Cl = 1.0, Cm_c/4 = -0.05
  a_0 = (1.0 - (-0.7))/ (8 - (-8))*rad2deg;
  m_0 = (-0.05 - (-0.025))/(8 - (-8))*rad2deg;
  d_ac_over_c = -m_0/a_0 
  Cm_ac = 1.0 * d_ac_over_c + (-0.05)  % the values of Cl and Cm_c/4 at alpha = 8 deg.
%
% Problem 2
%
% The airfoil of tail is NACA 0009
%from the airfoil data, 
%  
% at alpha = -8 deg Cl = -0.9, Cm_c/4 = 0
%    alpha = +8 deg Cl = 0.8, Cm_c/4 =  0
  a_0_009 = (0.8 - (-0.9))/ (8 - (-8))*rad2deg;
  m_0_009 = (0 - 0)/(8 - (-8))*rad2deg;
  d_ac_over_c_009 = -m_0_009/a_0_009 ;
  Cm_ac_009 = 0.8 * d_ac_over_c_009 + 0  ;% the values of Cl and Cm_c/4 at alpha = 8 deg.
%
% Calculate the epsilon_e from Cl_delta_t and a_0_009 for Cl at alpha = -8. deg, delta = 60 deg
   Cl_delta_t = (.4 - (-.9))/60 * rad2deg
%
% (a) Calculate the arodynamic characteristic of wing and tail
%
S_w = 20;
AR_w = 9;
lamda_w = 0.8;
b_w = sqrt(S_w * AR_w);
c_r_w = S_w/(1 + lamda_w)/(b_w/2)
mac_w = (2/3)*(1 + lamda_w + 2. * lamda_w * lamda_w)/(1 + lamda_w)
CL_alpha_w = a_0 /(1+ a_0/(pi*0.95*AR_w))
Cm_ac_w = Cm_ac
%
S_t = 2;
AR_t = 6;
b_t = sqrt(S_t * AR_t);
mac_t = S_t / b_t;
CL_alpha_t = a_0_009 /(1+ a_0_009/(pi*0.95*AR_t))
%
%  Calculate the epsilon_e from Cl_delta_t and CL_alpha_t for Cl at alpha = -8. deg, delta = 60 deg
   Cl_delta_t = (.4 - (-.9))/60 * rad2deg
   epsilon_e = Cl_delta_t /CL_alpha_t
%  Calculate Cm_ac_delta from Cm at alpha = -8. deg, delta = 60 deg
   Cm_ac_delta = -0.225/60 * rad2deg
%   
l_t_2_l_w = 6;
%
% (b) Determine the mounting angles of wing and tail at design point
%
mass = 1800;
%rho = 1.225
v_0 = 115 * 0.514; % from knots to m/sec
CL_trim = (mass*9.8)/(0.5*rho*v_0^2*S_w);
i_w = (CL_trim/CL_alpha_w)*rad2deg + alpha_CL0_w
epsilon_d0 = (0.6/AR_w)*CL_trim*rad2deg;
i_t = epsilon_d0
%
% (c)The location of the aerodynamic center of the wing relative to the center of gravity
%
l_w = Cm_ac_w/CL_trim*mac_w
%
% (d)The pitch stability derivative for the wing-tail combination
%
eta_t = 1;
l_t = l_t_2_l_w +l_w;
epsilon_d_alpha = (0.6/AR_w)*CL_alpha_w;
Cm_alpha = -(l_w/mac_w)*CL_alpha_w-(l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha)
CL_alpha = CL_alpha_w + (S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha);
%
% (e) determine the static margin
%
SM = -Cm_alpha/CL_alpha*100
%
% (f)  the trimmed angle of attack and elevator deflection for airspeed from 80 to 140 knots.
%
CL_0 = CL_trim;
Cm_0=Cm_ac_w-l_w/mac_w*CL_alpha_w*(i_w-alpha_CL0_w)/rad2deg ...
    -S_t*l_t/S_w/mac_w*CL_alpha_t*(i_t-epsilon_d0)/rad2deg;
CL_delta = (S_t/S_w)*eta_t*CL_alpha_t*epsilon_e;
Cm_delta = (mac_t/mac_w)*(S_t/S_w)*eta_t*Cm_ac_delta-(l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*epsilon_e;
%
figure(1)
i = 1;
for v_1_1 = 80 : 5 : 140
v_1 = v_1_1 * 0.514;
 CL_trim_1 = (mass*9.8)/(0.5*rho*v_1^2*S_w);
 d0 = CL_alpha*Cm_delta - CL_delta*Cm_alpha;
 alpha_trim(i) =((CL_trim_1 - CL_0)*Cm_delta + CL_delta*Cm_0)/d0*rad2deg;
 delta_trim(i) = (-(CL_trim_1 - CL_0)*Cm_alpha - CL_alpha*Cm_0)/d0*rad2deg;
 i = i+1;
end
temp_V=80:5:140;
    plot(temp_V,alpha_trim,'r',temp_V,delta_trim,'b')
    xlabel('Velocity (knots)')
    ylabel('Angle (deg)')
    legend('AOA','Delta e')
    fprintf('At 80 knots AOA =%f deg  Elevator angle =%f deg\n',alpha_trim(1),delta_trim(1))
    fprintf('At 140 knots AOA =%f deg  Elevator angle =%f deg\n',alpha_trim(end),delta_trim(end))
 %
% (g) repeating (f) when the CG is moved to the location of the a.c. of the wing.
%
l_w = 0;
l_t = l_t_2_l_w +l_w;
epsilon_d_alpha = (0.6/AR_w)*CL_alpha_w;
Cm_alpha =-(l_w/mac_w)*CL_alpha_w-(l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha);
CL_alpha = CL_alpha_w + (S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha);
Cm_0=Cm_ac_w-l_w/mac_w*CL_alpha_w*(i_w-alpha_CL0_w)/rad2deg ...
   -S_t*l_t/S_w/mac_w*CL_alpha_t*(i_t-epsilon_d0)/rad2deg;
 CL_delta = (S_t/S_w)*eta_t*CL_alpha_t*epsilon_e; 
 Cm_delta = (mac_t/mac_w)*(S_t/S_w)*eta_t*Cm_ac_delta-(l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*epsilon_e;
 %
 figure(2)
 %
 i = 1; 
for v_1_1 = 80 : 5 : 140 
    v_1 = v_1_1 * 0.514;
 CL_trim_1 = (mass*9.8)/(0.5*rho*v_1^2*S_w); 
 d0 = CL_alpha*Cm_delta -CL_delta*Cm_alpha;
 alpha_trim(i) =((CL_trim_1 - CL_0)*Cm_delta +CL_delta*Cm_0)/d0*rad2deg; 
 delta_trim(i) = (-(CL_trim_1 - CL_0)*Cm_alpha - CL_alpha*Cm_0)/d0*rad2deg; 
 i = i+1;
end
temp_V=80:5:140;
    plot(temp_V,alpha_trim,'r',temp_V,delta_trim,'b')
    xlabel('Velocity (knots)')
    ylabel('Angle (deg)')
    legend('AOA','Delta e')
    fprintf('At 80 knots AOA =%f deg  Elevator angle =%f deg\n',alpha_trim(1),delta_trim(1))
    fprintf('At 140 knots AOA =%f deg  Elevator angle =%f deg\n',alpha_trim(end),delta_trim(end))
%  
% (h) plot the elevator deflection required to trim the airplane at 80 knots as function
% of the location of the center of gravity measured relative to the aerodynamic
% center of the main wing
%
figure(3)
%
x_cg_0 = -Cm_ac_w/CL_trim*mac_w;
x_t = l_t_2_l_w;
epsilon_d_alpha = (0.6/AR_w)*CL_alpha_w;
Cm_alpha =-(l_w/mac_w)*CL_alpha_w-(l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha);
CL_alpha = CL_alpha_w + (S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha);
Cm_0=Cm_ac_w-l_w/mac_w*CL_alpha_w*(i_w-alpha_CL0_w)/rad2deg ...
   -S_t*l_t/S_w/mac_w*CL_alpha_t*(i_t-epsilon_d0)/rad2deg;
 CL_delta = (S_t/S_w)*eta_t*CL_alpha_t*epsilon_e; 
 Cm_delta = (mac_t/mac_w)*(S_t/S_w)*eta_t*Cm_ac_delta-(l_t/mac_w)*(S_t/S_w)*epsilon_e;
i = 1; 
v_1_1 = 80;
CL_trim_1 = (mass*9.8)/(0.5*rho*v_1^2*S_w); 
for x_cg = -1. : 0.01 : 0.5
%
 Cm_0=Cm_ac_w+x_cg/mac_w*CL_alpha_w*(i_w-alpha_CL0_w)/rad2deg ...
   -S_t*(x_t - x_cg)/S_w/mac_w*CL_alpha_t*(i_t-epsilon_d0)/rad2deg;
 Cm_alpha =(x_cg/mac_w)*CL_alpha_w-((x_t - x_cg)/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha);
 Cm_delta = (mac_t/mac_w)*(S_t/S_w)*eta_t*Cm_ac_delta-((x_t - x_cg)/mac_w); 
 d0 = CL_alpha*Cm_delta -CL_delta*Cm_alpha;
 alpha_trim(i) =((CL_trim_1 - CL_0)*Cm_delta +CL_delta*Cm_0)/d0*rad2deg; 
 delta_trim(i) = (-(CL_trim_1 - CL_0)*Cm_alpha - CL_alpha*Cm_0)/d0*rad2deg; 
 i = i+1;
end
temp_xcg=-1. : 0.01 : 0.5;
    plot(temp_xcg,alpha_trim,'r',temp_xcg,delta_trim,'b')
    xlabel('x_cg (m)')
    ylabel('Angle (deg)')
    legend('AOA','Delta e')
    %fprintf('At 80 knots AOA =%f deg  Elevator angle =%f deg\n',alpha_trim(1),delta_trim(1))
    %fprintf('At 140 knots AOA =%f deg  Elevator angle =%f deg\n',alpha_trim(end),delta_trim(end))
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
% [1] Yunus A. Cengel & John M.Cimbala¡§FLUID MECHANICS ...
% Fundamentals and  Applications¡¨, McGraw-Hill., p897.
% [2] John J. Bertin & Russell M. Cummings¡§AERODYNAMICS FOR ENGINEERS¡¨,
% 5th Edition, Pearson Education International., p21-p43.
% [3] WARREN F.PHILLIPS,¡§MECHANICS of FLIGHT¡¨, 2nd Edition, John Wiley.
% p10-p14.
% [4] John D.Anderson,"Modern Compressibe Flow", third Edition, McGraw-Hill
% ., p585-p613.
% =========================== Separator ==================================
end