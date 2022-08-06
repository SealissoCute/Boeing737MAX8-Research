clc;
clear;
%% longitudinal system matrix
A=[-0.163 0.6821 0 -9.8; -2.2813 -5.236 8.6 0; 0.4894 -0.0838 -6.8998 0;0 0 1 0]

%% (a) Obtain the eigenvectors of the system,let delta_theta = 1.
[u v]=eig(A);
fprintf('eigenvectors\n')
u

%% (b) What are the phase angle of delta_u, delta_alpha, delta_q, behind delta_theta?

% short period 
fprintf('short period\n')
angle_u_1=phase(u(:,1))*57.3;
angle_u_2=phase(u(:,2))*57.3;
angle_u_1_minus_theta=angle_u_1-angle_u_1(4) 
angle_u_2_minus_theta=angle_u_2-angle_u_2(4)

% long-period 
fprintf('Long period\n') 
angle_u_3=phase(u(:,3))*57.3;
angle_u_4=phase(u(:,4))*57.3;
angle_u_3_minus_theta=angle_u_3-angle_u_3(4) 
angle_u_4_minus_theta=angle_u_4-angle_u_4(4) 
 
