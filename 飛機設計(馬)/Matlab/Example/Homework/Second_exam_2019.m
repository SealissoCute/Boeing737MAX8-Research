%% m-file of the second exam Part II, 2018/08/21
%
clear
clc
%% (a) Calculate the mean aerodynamic chord
% mac_w = 0.6000 ; mac_t = 0.4860
%
c_r_w = 0.79;
c_t_w = 0.355;
lambda_w = c_t_w /c_r_w;
mac_w =(2/3)*((1+lambda_w+lambda_w^2)/(1+lambda_w))*c_r_w
c_r_t = 0.6;
c_t_t = 0.35;
lambda_t = c_t_t /c_r_t;
mac_t =(2/3)*((1+lambda_t+lambda_t^2)/(1+lambda_t))*c_r_t
%
%% (b) Determine the mounting angles of wing and tail at design point
% i_w = 4.5455 ; i_t = 3.8960
%
rad2deg = 90./asin(1);
mass = 17.2;
S_w = 3.28;
b_w = 5.73;
AR_w = b_w^2/S_w;
rho = 1.225;
v_0 = 8.6;
CL_alpha_w = 5.181;
alpha_CL0_w = -8.0;
CL_trim = (mass*9.8)/(0.5*rho*v_0^2*S_w)
i_w = (CL_trim/CL_alpha_w)*rad2deg + alpha_CL0_w
epsilon_d0 = (0.6/AR_w)*CL_trim*rad2deg;
i_t = epsilon_d0
%
%% (c)The location of the aerodynamic center of the wing relative to the center of gravity
% Assuming the minimum drag axis of fuselage is aligned with the FRL
% l_w = -0.1047
%
Cm_ac_w = -0.198;
l_w = Cm_ac_w/CL_trim*mac_w
%
%% (d) The pitch stability derivative for the wing-tail combination and the static margin
% Considering the effects of the fuselage
%
CL_alpha_t = 4.43;
eta_t = 1;
S_t = 1.14;
l_t = 1.785+l_w;
epsilon_d_alpha = (0.6/10.01)*5.181;
%
S_f = 0.05; % Maximum cross-sectional area of the fuselage
x_Sf = 1;   % The position of maximum cross-sectional area from nose
x_com = 1.8;    % The position of center of mass from nose
l_f = x_Sf/2-x_com; % Length from CG to center of presure of fuselage
d_f = 2*sqrt(S_f/pi);
c_f = 4;    % Length of fuselage
%
Cm_alpha = -(l_w/mac_w)*CL_alpha_w - (l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha) - 2*(S_f/S_w)*(l_f/mac_w)*(1-1.76*sqrt((d_f/c_f)^3))
CL_alpha = CL_alpha_w + (S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha) + 2*(S_f/S_w)*(1-1.76*sqrt((d_f/c_f)^3));
SM = -Cm_alpha/CL_alpha*100
%
%% (e) Plot the trimmed angle of attack and elevator deflection for airspeed from 7 to 10 m/sec
%
eta_e = 0.15;
Cm_ac_delta = -0.607;
CL_0 = CL_trim;
Cm_0 = 0;
%
v_1 = 7:0.1:10;
CL_trim_1 = (mass*9.8)./(0.5*rho*v_1.^2*S_w);
CL_delta = (S_t/S_w)*eta_t*CL_alpha_t*eta_e;
Cm_delta = (mac_t/mac_w)*(S_t/S_w)*eta_t*Cm_ac_delta-(l_t/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*eta_e;
d0 = CL_alpha*Cm_delta - CL_delta*Cm_alpha;
alpha_trim = ((CL_trim_1 - CL_0)*Cm_delta + CL_delta*Cm_0)/d0*rad2deg;
delta_trim = (-(CL_trim_1 - CL_0)*Cm_alpha - CL_alpha*Cm_0)/d0*rad2deg;
%
figure(1)
[Ax,Line_alpha,Line_delta]= plotyy(v_1,alpha_trim,v_1,delta_trim);
set(Line_delta,'LineStyle','--')
title('Trimmed Angles of Attack & Elevator Deflections')
xlabel('Airspeed (m/s)')
ylabel(Ax(1),'Degree (\alpha_{trim})')
ylabel(Ax(2),'Degree (\delta_{trim})')
legend('\alpha_{trim}','\delta_{trim}','Location','east')
grid on
%
%% (f) Plot the trimmed angle of attack and elevator deflection for airspeed
% from 7 to 10 m/sec when the CG is moved to the location of the a.c. of
% the wing
%
l_w_new = 0;
l_t_new = 1.785+l_w_new;
l_f_new = x_Sf/2-x_com-l_w; % Length from CG to center of presure of fuselage
%
Cm_alpha_new = -(l_w_new/mac_w)*CL_alpha_w - (l_t_new/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(1-epsilon_d_alpha) - 2*(S_f/S_w)*(l_f_new/mac_w)*(1-1.76*sqrt((d_f/c_f)^3));
Cm_delta_new = (mac_t/mac_w)*(S_t/S_w)*eta_t*Cm_ac_delta-(l_t_new/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*eta_e;
d0_new = CL_alpha*Cm_delta_new - CL_delta*Cm_alpha_new;
Cm_0_new = Cm_ac_w - (l_w_new/mac_w)*CL_alpha_w*(i_w-alpha_CL0_w) - (l_t_new/mac_w)*(S_t/S_w)*eta_t*CL_alpha_t*(i_t-epsilon_d0);
alpha_trim_new = ((CL_trim_1 - CL_0)*Cm_delta_new + CL_delta*Cm_0_new)/d0_new*rad2deg;
delta_trim_new = (-(CL_trim_1 - CL_0)*Cm_alpha_new - CL_alpha*Cm_0_new)/d0_new*rad2deg;
%
figure(2)
[Ax_new,Line_alpha_new,Line_delta_new]= plotyy(v_1,alpha_trim_new,v_1,delta_trim_new);
set(Line_delta_new,'LineStyle','--')
title('Trimmed Angles of Attack & Elevator Deflections (CG at a.c. of the wing)')
xlabel('Airspeed (m/s)')
ylabel(Ax_new(1),'Degree (\alpha_{trim\_new})')
ylabel(Ax_new(2),'Degree (\delta_{trim\_new})')
legend('\alpha_{trim\_new}','\delta_{trim\_new}','Location','east')
grid on
%
%% The yaw stability derivative for the plane
% Considering the effects of the fuselage
% 
S_v = 0.5;
b_v = 1;
c_r_v = 0.56;
c_t_v = 0.28;
CL_alpha_v = 3.4;
l_v = l_t;
h_v = 0.06;
eta_v = 1;
epsilon_s_beta = -0.1;
epsilon_s0 = 0;
epsilon_r = 0.56;
Cm_delta_r = -0.51;
%
Cn_beta_v = (S_v/S_w)*(l_v/b_w)*eta_v*CL_alpha_v*(1-epsilon_s_beta)
Cn_beta_f = 2*(S_f/S_w)*(l_f/b_w)*(1-1.76*(sqrt(d_f/c_f))^3)
% Cn_beta = Cn_beta_v+2*Cn_beta_f
Cn_beta = Cn_beta_v+ Cn_beta_f