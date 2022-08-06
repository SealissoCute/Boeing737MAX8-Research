% Homework #9 and #10
%  
clear all
clc
clf
%
clear;
clc;
%
% at alitude 100
altitude = 100
[C,speedofsound,P,rho,g,mu]=Standard_Atmosphere(altitude)
%
M=17.2
S_wing = 3.28
c_mac = 0.6
v0 = 8.6
Iy=18.492
t_ref = c_mac/(2.*v0)
%
CX_u = -0.1638 ; CX_alpha = 0.6856; % CX_q =0.; CX_alphadot = 0.;
CZ_u = -2.2930 ; CZ_alpha = -5.2629; CZ_q = -9.6668; CZ_alphadot = -3.5286;
CM_u = 0.; CM_alpha = -2.1741 ; CM_q = -30.3462; CM_alphadot = -11.0770;
%
Q=0.5*rho*v0^2;
X_u=CX_u*(1/v0)*Q*S_wing/M
X_w=CX_alpha*Q*S_wing/(M*v0)

Z_u=CZ_u*(1/v0)*Q*S_wing/M
Z_w=CZ_alpha*Q*S_wing/(M*v0)
Z_wdot=CZ_alphadot*(c_mac/2*v0)*Q*S_wing/(M*v0)
Z_q=CZ_q*(c_mac/(2*v0))*Q*S_wing/M
% Z_elevator=CZ_elevator*Q*S_wing/M

Moment_u=CM_u*(1/v0)*Q*S_wing*c_mac/Iy
Moment_w=CM_alpha*Q*S_wing*c_mac/(Iy*v0)
Moment_wdot= CM_alphadot*(c_mac/(2*v0))*Q*S_wing*c_mac/(Iy*v0)
Moment_q=CM_q*(c_mac/(2*v0))*Q*S_wing*c_mac/Iy
%Moment_elevator=CM_elevator*Q*S_wing*c_mac/Iy
%
% (a) the longitudinal system matrix
%
A = [X_u,X_w,0,-9.8;
    Z_u,Z_w,v0,0;
    Moment_u+Moment_wdot*Z_u,Moment_w+Moment_wdot*Z_w,Moment_q+Moment_wdot*v0,0;
    0,0,1,0]
%
% (b) the eigenvalues of the system
%
[u,d]=eig(A)
%
% u -- eigenvectors, d -- eigenvalues
%
% (c) the time-to-half and period of each mode
%
% Short period mode
%
lambda_sp = d( 1, 1 )
w_sp = imag( lambda_sp )
eta_sp = real( lambda_sp )  
P_sp = 2 * pi /  w_sp
T_half_sp = 0.69 / abs( eta_sp  )
N_half_sp = 0.11 * w_sp / abs( eta_sp )
%
% Phugoid mode
%
lambda_P = d( 3, 3 );
w_P = imag( lambda_P );
eta_P = real( lambda_P );  
P_P = 2 * pi /  w_P
T_half_P = 0.69 / abs( eta_P  )
N_half_P = 0.11 * w_P / abs( eta_P )
%
% (d) eigenvector with delta\theta = 1
% (e) phase angles relative with \theta
%
% short period eigenvalue and eigenvector
%
angle_u=phase(u(:,1))*57.3 
Amplitude_u=[abs(u(1,1))/v0;abs(u(2,1))/v0;abs(u(3,1))*t_ref;abs(u(4,1))] 
Amplitude_U_over_theta=Amplitude_u/abs(u(4,1)) 
angle_u_minus_theta=angle_u-angle_u(4)
%
% long-period eigenvalue and eigenvector
%
angle_u_l=phase(u(:,3))*57.3 
Amplitude_u_l=[abs(u(1,3))/v0;abs(u(2,3))/v0;abs(u(3,3))*t_ref;abs(u(4,3))] 
Amplitude_U_over_theta_l=Amplitude_u_l/abs(u(4,3)) 
angle_u_minus_theta=angle_u_l-angle_u_l(4)
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
% [1] Yunus A. Cengel & John M.Cimbala��FLUID MECHANICS ...
% Fundamentals and  Applications��, McGraw-Hill., p897.
% [2] John J. Bertin & Russell M. Cummings��AERODYNAMICS FOR ENGINEERS��,
% 5th Edition, Pearson Education International., p21-p43.
% [3] WARREN F.PHILLIPS,��MECHANICS of FLIGHT��, 2nd Edition, John Wiley.
% p10-p14.
% [4] John D.Anderson,"Modern Compressibe Flow", third Edition, McGraw-Hill
% ., p585-p613.
% =========================== Separator ==================================
end



