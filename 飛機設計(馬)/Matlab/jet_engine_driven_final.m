% Aircraft Design for Jet Transport
%Length transfer
%1 nm= 1.15 sm
%velocity
%1 kts= 1.15 mph
clear;clc
A=-.4179;
B=1.1446;
M_ff=0.35108;
C=1-(1)*(1-M_ff)-0.005;
D=2000;
%% 
% M_ff= 0.4276;
W_To_guess=1;
n=0;
while 1
n=n+1;
W_To_guess=W_To_guess+1;
W_tfo=double(0.005*W_To_guess); 
W_F=(1-M_ff)*W_To_guess;  % Fuel used
W_E_tent=double(W_To_guess-W_F-W_tfo-D);
W_E_real=double(10^((log10(W_To_guess)-A)/B));
error=double((abs(W_E_real-W_E_tent)/W_E_real))*100;
F=double(-B*W_To_guess^2*((C*W_To_guess*(1-B)-D)^-1)*(1+0)*M_ff);
if abs(error)<1e-4
    break
elseif n>=100000000
    break
end
end
disp('----------------------------------------------------')
disp('Iteration Result')
string=['W_to = ',num2str(W_To_guess),' lbs'];
disp(string);
string1=['W_E = ',num2str(W_E_tent),' lbs'];
disp(string1);
string2=['W_F = ',num2str(W_F),' lbs'];
disp(string2);
string3=['Error = ',num2str(error),' %'];
disp(string3);
string3=['F = ',num2str(F),' lbs'];
disp(string3);
string=['Iteration = ',num2str(n)];
disp(string)

%% 

R=3500;
L_over_D=5.5;
c_j_cruise=1.3;
L_over_D_loiter=8;
c_j_loiter=0.7;
velocity=1781.89;
E=1;
Wto_DR=F*c_j_cruise*(velocity*L_over_D)^-1;
Wto_c_j=F*R*(velocity*L_over_D)^-1;
Wto_v=-F*R*c_j_cruise*(velocity^2*L_over_D)^-1;
Wto_L_over_D=-F*R*c_j_cruise*(velocity*L_over_D^2)^-1;
Wto_E=F*c_j_loiter*(L_over_D_loiter)^-1;
Wto_c_j_loiter=F*E*(L_over_D_loiter)^-1;
Wto_L_over_D_loiter=-F*E*c_j_loiter*(L_over_D_loiter)^-2;
disp('----------------------------------------------------')
disp('sensitivity')
string=['Wto_DR = ',num2str(Wto_DR),' lbs/nm'];
disp(string);
string1=['Wto_c_j = ',num2str(Wto_c_j),' lbs/lbs/hp/hr'];
disp(string1);
string2=['Wto_v = ',num2str(Wto_v),' lbs/kt'];
disp(string2);
string3=['Wto_L_over_D = ',num2str(Wto_L_over_D),' lbs'];
disp(string3);
string3=['Wto_E = ',num2str(Wto_E),' lbs'];
disp(string3);
string3=['Wto_c_j_loiter = ',num2str(Wto_c_j_loiter),' lbs'];
disp(string3);
string3=['Wto_L_over_D_loiter = ',num2str(Wto_L_over_D_loiter),' lbs'];
disp(string3);