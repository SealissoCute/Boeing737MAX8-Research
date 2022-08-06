%
% Wing Design Example
%
% Problem Statement:
% Design a wing for a normal category GA aircraft with the following features:
%    S = 18.1m2, m = 1800 kg, Vc = 130 knot (at sea level), Vs = 60 knot
% Assume the aircraft has a monoplane high wing and employs a split flap.
%
 clc   
     S = 18.1 
     Mass = 1800 % Kg, assuming average mass and maximum take-off mass
     Mass_ave = Mass
     Mass_TO = Mass
     %g = 9.81 % m/sec2,  gravity acceleration
     V_c = 130 % cruising speed, Knots, at altitude 3000 m
     %rho = 1.225 % kg/m3, density at sea level
     V_s = 60 % stall speed, Knots, sea level
% crusing altitude     
     h = 3000
     [Temp,speedofSound,Pressure,rho_3000,g,mu_3000]=Standard_Atmosphere(h);
%     
% The number of wings and wing vertical position are stated in the problem statement, 
% so we do not need to investigate these two parameters.
%
% 1. Dihedral angle:
%   Since the aircraft is a high-wing, low-subsonic, and mono-wing aircraft, 
%  based on table given, a ?5 deg of anhedral is selected. This value will 
% be revised and optimized when other aircraft components are designed during 
% lateral stability analysis.
%
% 2. Sweep angle: 
%   The aircraft is a low-subsonic prop-driven normal category aircraft. 
% To keep the cost low in the manufacturing process, we select no sweep angle 
% at 50% of wing chord. However, we may need to taper the wing; hence the leading edge 
% and trailing edge may have sweep angles.
%
% 3. Airfoil Selection:
%   To be fast in the wing design, we select an airfoil from NACA selections.
% The design of an airfoil is beyond the scope of this book. The selection 
% process of an airfoil for the wing requires some calculation as follows.
%
% Section¡¦s ideal lift coefficient:
%
  CLC_w = 2.*Mass_ave*g/(rho_3000*(V_c*0.514)^2*S)
  Cl_i = CLC_w /0.9
%
% Section¡¦s maximum lift coefficient (see level):
%
 h = 0 % take_off at sea level
 [Temp,speedofSound,Pressure,rho_0,g,mu]=Standard_Atmosphere(h);
  Clmax_w = 2.*Mass_TO*g/(rho_0*(V_s*0.514)^2*S)
  Clmax_gross =  Clmax_w/0.9
%
% The aircraft has a split flap, and the split flap generates a 
% DCL of 0.36 when deflected 30 deg. 
%
  D_Clmax_HLD = 0.36
  Clmax =  Clmax_gross - D_Clmax_HLD 
%
% NACA airfoil section 64(2)-415 is chosen.
%
%  4. The lift distribution of the wing (AR = 7, £f = 0.8, £\t  = ?1.5, iw = 3. deg)
%
N = 9; % (number of segments - 1)
%S = 25; % m^2
AR = 7; % Aspect ratio
lambda = 0.8; % Taper ratio
alpha_twist = -1.5; % Twist angle (deg)
i_w = 1.86; % wing setting angle (deg)
a_2d = 6.3; % lift curve slope (1/rad)
alpha_0 = -3.; % zero-lift angle of attack (deg)
b = sqrt(AR*S); % wing span (m)
bf_b=0.6; %flap-to-wing span ratio
%Croot = (1.5*(1+lambda)*MAC)/(1+lambda+lambda^2); % root chord (m)
Croot = 2*S/(1+lambda)/b
theta = pi/(2*N):pi/(2*N):pi/2;
alpha = i_w+alpha_twist:-alpha_twist/(N-1):i_w;
% segment¡¦s angle of attack
z = (b/2)*cos(theta);
c = Croot * (1 - (1-lambda)*cos(theta)); % Mean Aerodynamics
%Chord at each segment (m)
mu = c * a_2d / (4 * b);
LHS = mu .* (alpha-alpha_0)/57.3; % Left Hand Side
% Solving N equations to find coefficients A(i):
for i=1:N
for j=1:N
B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1)) /sin(theta(i)));
end
end
A=B\transpose(LHS);
for i = 1:N
sum1(i) = 0;
sum2(i) = 0;
for j = 1 : N
sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
end
end
CL = 4*b*sum2 ./ c;
CL1=[0 CL(1) CL(2) CL(3) CL(4) CL(5) CL(6) CL(7) CL(8) CL(9)];
y_s=[b/2 z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9)];
plot(y_s,CL1,'-o')
grid
title('Lift distribution')
xlabel('Semi-span location (m)')
ylabel ('Lift coefficient')
A
CL_wing = pi * AR * A(1)
delta = 3*(A(3)/A(1))^2+5*(A(5)/A(1))^2+7*(A(7)/A(1))^2+9*(A(9)/A(1))^2
CD_wing_i = CL_wing^2/(pi*AR)*(1+delta)
%
%  5. Flap parameters: A flap is usually employed during take-off and landing operations.
% The design of the flap based on the take-off requirements and adjust it for the
% landing requirements. The take-off speed for a GA aircraft is about 20% faster than
% the stall speed 
%
   V_TO = 1.2 * V_s * 0.514
   CL_max_W = 2*Mass_TO*g/(rho_0*V_TO^2*S)
%
cf_c = 0.2
D_alpha_0 = -1.15*cf_c*13
%
N = 9; % (number of segments - 1)
%S = 25; % m^2
AR = 7; % Aspect ratio
lambda = 0.8; % Taper ratio
alpha_twist = -1.5; % Twist angle (deg)
i_w = 8.9; % wing setting angle (deg)
a_2d = 6.3; % lift curve slope (1/rad)
a_0 = -3; % flap up zero-lift angle of attack (deg)
a_0_fd = -6; % flap down zero-lift angle of attack (deg)
b = sqrt(AR*S); % wing span (m)
% Croot = (1.5*(1+lambda)*MAC)/(1+lambda+lambda^2); % root chord (m)
Croot = 2*S/(1+lambda)/b
theta = pi/(2*N):pi/(2*N):pi/2;
alpha = i_w+alpha_twist:-alpha_twist/(N-1):i_w;
% segment¡¦s angle of attack
for i=1:N
if (i/N)>(1-bf_b)
alpha_0(i)=a_0_fd; %flap down zero lift AOA
else
alpha_0(i)=a_0; %flap up zero lift AOA
end
end
%
z = (b/2)*cos(theta);
c = Croot * (1 - (1-lambda)*cos(theta)); % Mean Aerodynamics
%Chord at each segment (m)
mu = c * a_2d / (4 * b);
LHS = mu .* (alpha-alpha_0)/57.3; % Left Hand Side
% Solving N equations to find coefficients A(i):
for i=1:N
for j=1:N
B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1)) /sin(theta(i)));
end
end
A=B\transpose(LHS);
for i = 1:N
sum1(i) = 0;
sum2(i) = 0;
for j = 1 : N
sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
end
end
CL = 4*b*sum2 ./ c;
CL1=[0 CL(1) CL(2) CL(3) CL(4) CL(5) CL(6) CL(7) CL(8) CL(9)];
y_s=[b/2 z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9)];
% plot(y_s,CL1,'-o')
% grid
% title('Lift distribution')
% xlabel('Semi-span location (m)')
% ylabel ('Lift coefficient')
% A
CL_TO_wing = pi * AR * A(1)
delta = 3*(A(3)/A(1))^2+5*(A(5)/A(1))^2+7*(A(7)/A(1))^2+9*(A(9)/A(1))^2;
CD_TO_wing_i = CL_TO_wing^2/(pi*AR)*(1+delta);
%
% 6. Other wing parameters. To determine the other wing parameters (i.e., wing span
% (b), root chord (Cr), tip chord (Ct), and MAC),
%
  b = sqrt(AR*S)
  Croot = 2*S/(1+lambda)/b
  Ctip = Croot * lambda
  MAC = Croot *(2/3)*(1+lambda+lambda^2)/(1+lambda)
  b_f = bf_b*b
  c_f = cf_c * MAC
  Reynold = rho_3000*V_c*0.514*MAC/mu_3000
  
  %
  
  
  
