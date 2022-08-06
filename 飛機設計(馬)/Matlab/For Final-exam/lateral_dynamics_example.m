% Lateral dynamics example
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
%
A = [Cyb/two_mu Cyp/two_mu Cyr/two_mu-1 CL0/two_mu; 
    (iz*Clb+ixz*Cnb)/Cap_I (iz*Clp+ixz*Cnp)/Cap_I (iz*Clr+ixz*Cnr)/Cap_I 0; 
    (ixz*Clb+ix*Cnb)/Cap_I (ixz*Clp+ix*Cnp)/Cap_I (ixz*Clr+ix*Cnr)/Cap_I 0; 
     0 1 0 0]
    
[u v] = eig(A)
%
%%
fprintf('dutch roll mode\n') 
angle_u = phase(u(:,2))*57.3 
Amplitude_u = abs(u(:,2)) 
Amplitude_U_over_phi = Amplitude_u/abs(u(4,2)) 
angle_u_minus_phi = angle_u-angle_u(4)
%