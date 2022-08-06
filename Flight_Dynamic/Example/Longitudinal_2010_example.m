clear; 
A=[-0.045 0.036 0 -32.2; -0.369 -2.02 176 0; 0.0019 -0.0396 -2.948 0;0 0 1 0] 
[u v]=eig(A) 
% 
% short period eigenvalue and eigenvector 
% 
fprintf('short period\n')
angle_u=phase(u(:,1))*57.3 
Amplitude_u=[abs(u(1,1))/176;abs(u(2,1))/176;abs(u(3,1))*0.016;abs(u(4,1))] 
Amplitude_U_over_theta=Amplitude_u/abs(u(4,1)) 
angle_u_minus_theta=angle_u-angle_u(4) 
% 
% long-period eigenvalue and eigenvector 
% 
fprintf('Long period\n') 
angle_u_l=phase(u(:,3))*57.3 
Amplitude_u_l=[abs(u(1,3))/176;abs(u(2,3))/176;abs(u(3,3))*0.016;abs(u(4,3))] 
Amplitude_U_over_theta_l=Amplitude_u_l/abs(u(4,3)) 
angle_u_minus_theta=angle_u_l-angle_u_l(4)
