%% m-file of Example LO#1
%
clear
clc
%%
 rad2deg = 90./asin(1);
 mph2ftpersec = 5280./3600;
% 
 weight = 2700; % airplaen wieght lb
%
% Flight condition
 v_0 = 140   ;  % mile per hour
 v_0 = v_0*mph2ftpersec; % ft/sec
 h_0 = 0 ;      % sea level
 [C,a,P,rho,g,mu]=Standard_Atmosphere(h_0); % in MKS system
 rho = rho * 0.0624279/32 ;% slug/cubic feet
%
% Give the wing and tail characteristics:
%
 S_w = 180   ;        % wing surface area, sq_ft
 b_w = 33    ;        % wing span, ft
 AR_w = b_w^2 /S_w ;  % wing aspect ratio
 CL_alpha_w = 4.44 ;  % wing lfit slope per radian
 alpha_CL0_w = -2.20; % wing zero lift AOA (deg)
% 
% Keep the same as Example LO1
%
 Cm_ac_w = -0.053 ;   % wing A.C. moment coefficient 
%
 S_t = 36  ;          % tail surface area, sq_ft
 b_t = 12 ;           % tail span, ft  
 AR_t = b_t^2 / S_t ; % tail aspect ratio
 CL_alpha_t = 3.97 ;  % tail lfit slope per radian
 eta_t = 1   ;        % tail efficiency
 epsilon_e = 0.60 ;   % elevator efectiveness
 Cm_ac_delta = -0.55; % tail A.C. moment coefficient due to elevator deflection
%
 l_t2l_w = 15   ;     % distance from tail a.c. to wing a.c., ft
%
%% Calculate the mean aerodynamic chord
%
%  Given Trapezoid wing with taper rato lambda
%
   lambda = 0.8
   C_r = 2*S_w/(b_w*(1+lambda))
   mac_w = 2/3*((1+lambda+lambda^2)/(1+lambda))*C_r
%   
% mac_w = S_w / b_w ; mac_t = S_t/ b_t ;% assume rectangular planform
%
%% (a) Determine the mounting angles of wing and tail at design point
%
CL_trim = weight/(0.5*rho*v_0^2*S_w)
i_w = (CL_trim/CL_alpha_w)*rad2deg + alpha_CL0_w
epsilon_d0 = (0.6/AR_w)*CL_trim*rad2deg;
i_t = epsilon_d0
%
%% (b)The location of the aerodynamic center of the wing relative to the center of gravity
% Assuming the minimum drag axis of fuselage is aligned with the FRL
%
l_w = Cm_ac_w/CL_trim*mac_w
%
%% (c) The pitch stability derivative for the wing-tail combination and the static margin
%
eta_t = 1;
l_t = l_t2l_w +l_w;
epsilon_d_alpha = 0.6*CL_alpha_w/AR_w;
%
Cm_alpha = -(l_w/mac_w)*CL_alpha_w - ... 
           (l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha)
CL_alpha = CL_alpha_w + (S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha);
% Example LO#2
 SM = -Cm_alpha/CL_alpha*100;

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