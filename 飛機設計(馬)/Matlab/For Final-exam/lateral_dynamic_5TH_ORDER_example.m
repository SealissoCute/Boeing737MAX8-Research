% Lateral dynamics (5th-order) example - Dutch modes
%
rho = 0.00088;
u0 = 733;
S = 1667;
b = 108;
mass = 100000/32.2;
mu = 2.*mass/(rho*S*b);
two_mu = 2.* mu;
%
ix = 3.69;
iz = 9.22;
ixz = -0.39;
Cap_I = ix*iz - ixz*ixz;
%
CL0 = 0.25;
Cyb = -0.168;
Cyp = 0;
Cyr = 0.192;
Clb = -0.022- 0.10*CL0;
Clp = -0.43;
Clr = 0.0077 + 0.25*CL0;
Cnb = 0.036 + 0.04*CL0*CL0;
Cnp = 0.008 - 0.10*CL0;
Cnr = -0.116-0.020*CL0*CL0;
%
% 4th order system matrix
a11=Cyb/two_mu;
a12=Cyp/two_mu;
a13= Cyr/two_mu-1;
a14=CL0/two_mu;
a15=0;
%
a21= (iz*Clb+ixz*Cnb)/Cap_I;
a22=(iz*Clp+ixz*Cnp)/Cap_I;
a23=(iz*Clr+ixz*Cnr)/Cap_I;
a24=0;
a25=0;
%
a31=(ixz*Clb+ix*Cnb)/Cap_I;
a32=(ixz*Clp+ix*Cnp)/Cap_I;
a33=(ixz*Clr+ix*Cnr)/Cap_I;
a34=0;
a35=0;
%
a41=0;
a42=1;
a43=0;
a44=0;
a45=0;
%
a51=0;
a52=0;
a53=1;
a54=0;
a55=0;
%
A=[a11 a12 a13 a14 a15; 
   a21 a22 a23 a24 a25;
   a31 a32 a33 a34 a35; 
   a41 a42 a43 a44 a45;
   a51 a52 a53 a54 a55] 

B=[0 0;0 0;0 0;0 0;0 0]

C=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1]

D=[0 0;0 0;0 0;0 0;0 0]

[u v] = eig(A)
%
fprintf('dutch roll mode\n')
angle_u = phase(u(:,4))*57.3
Amplitude_u = abs(u(:,4))
Amplitude_U_over_beta = Amplitude_u/abs(u(1,4))
angle_u_minus_beta = angle_u-angle_u(1)
%
G = ss(A,B,C,D)
x0 = [1; 0; 0; 0; 0]
initial(G, x0, 400.)
grid