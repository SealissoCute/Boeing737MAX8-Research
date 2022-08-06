clc;
clear;
% Input matrix
C = [-0.6645 0.6645 -0.3420; -0.7472 -0.5817 0.3214; 0.0146 0.4691 0.8830];
%
% Calculate Euler angles
theta = asin(-C(1,3));
phi = atan2(C(2,3),C(3,3));
psi = atan2(C(1,2),C(1,1));
rad2deg = 360/(2*pi);
theta_deg = theta*rad2deg
phi_deg = phi*rad2deg
psi_deg = psi*rad2deg
%
% Calculate quarternion
e0 = 0.5*sqrt(1+cos(theta)*cos(psi)+cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(theta))
e1 = (sin(phi)*cos(theta)+sin(phi)*cos(psi)-cos(phi)*sin(theta)*sin(psi))/(4*e0)
e2 = (sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)+sin(theta))/(4*e0)
e3 = (cos(theta)*sin(psi)+cos(phi)*sin(psi)-sin(phi)*sin(theta)*cos(psi)) /(4*e0)
%