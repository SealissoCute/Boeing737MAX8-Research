%% Reverse Engineering of Boeing 737 MAX 8
% Clear previous data

format long
clc
clear
%% Unit exchange
%
ft_to_m = 0.3048;       % ft to m
m_s_to_mph = 2.236936;  % m/s to mph
m_s_to_kt = 1.943844;   % m/s to kt
m_s_to_ft_s = 3.280840; % m/s to ft/s
ft_s_to_kt = 0.592484;  % ft/s to kt

%% Parameters
% Take-off weight (unit: lbs)
W_TO = 182043.622463998;

% Parameters at CruiseAltitude
CruiseAltitude = 40000; % unit: ft
[a,rho]=Standard_Atmosphere(CruiseAltitude);
a_CruiseAltitude = a;
rho_CruiseAltitude = rho;
CruiseSpeed_Mach = 0.79;
CruiseSpeed = CruiseSpeed_Mach*a_CruiseAltitude;
q_overline = 0.5*rho_CruiseAltitude*CruiseSpeed^2;

% Parameters at FieldAltitude
FieldAltitude = 5000; % unit: ft
[a,rho]=Standard_Atmosphere(FieldAltitude);
rho_FieldAltitude = rho;

%
WoverS = 0:10:200;
c = 0.0199;
d = 0.7531;
S_wet = 10^(c+d*log10(W_TO));

% From Table 3.4 Correlation Coefficients For Parasite Area Versus Wetted Area
cf_2 = 0.003; a_2 = -2.5229; b_2 = 1;
cf_3 = 0.004; a_3 = -2.3979; b_3 = 1;
f_2 = 10^(a_2+b_2*log10(S_wet));
f_3 = 10^(a_3+b_3*log10(S_wet));

%
delta_CD0_TOflaps = 0.015;  % From p.127, Table 3.6
delta_CDO_Lflaps = 0.065;   % From p.127, Table 3.6
delta_CD0_LG = 0.02;        % From p.127, Table 3.6CD_0_clean
e_TOflaps = 0.8;            % From p.127, Table 3.6
e_Lflaps = 0.75;            % From p.127, Table 3.6
W_L = W_TO*0.84;            % 0.84 is from p.107, Table 3.3
S = W_TO/100;               % W/S = 100
CD_0_clean = f_2/S;         % Take cf = 0.003

%
C_D0 = 0.0184; % p.145&182 low speed,clean drag polar
delta_C_D0 = 0.0001*2.5; % p.166 figure 3.32
C_D0_modification = C_D0 + delta_C_D0;

%
W_TO_wiki = 182200; % unit: lb
b = 117.833; % unit: ft^2
S_wiki = 1370; % unit: ft
AR = b^2/S_wiki; % unit: ft
S_wiki_TO = 1370*1.3; % unit: ft^2
S_wet_wiki = 10^(c+d*log10(W_TO_wiki));
StaticThrust_TO = 29317; % unit: lbs
WoverS_TO_wiki = W_TO_wiki/S_wiki_TO; % unit: lb/ft^2
ToverW_TO_wiki = StaticThrust_TO*2/W_TO_wiki; % unit: lb/lb

%%
figure()
hold on
for CL_max_TO = 1.6:0.2:3
    ToverW = (0.009640.*WoverS)/CL_max_TO;
    plot(WoverS,ToverW,'color',[0 0.4470 0.7410]); % blue
end

for CL_max_L = 1.8:0.2:3
    WoverS_landing = FieldAltitude/0.3/1.69/2*rho_FieldAltitude*CL_max_L/(ft_s_to_kt^2)/0.85;
    ToverW_landing = [0 1.6];
    WoverS_landing = [WoverS_landing WoverS_landing];
    plot(WoverS_landing,ToverW_landing,'color',[0.4660 0.6740 0.1880]); % green
end

for e_clean = 0.8:0.05:0.85              % From p.127, Table 3.6
    ToverW_cruise_reqd = C_D0_modification*q_overline./WoverS + WoverS./(q_overline*pi*AR*e_clean);
    ToverW_TO = ToverW_cruise_reqd./0.191;
    plot(WoverS,ToverW_TO,'color',[0.9290 0.6940 0.1250]); % orange
end

% FAR25.111 OEI (P.145)
CL_TO_max = 2;                     % From Table 3.1
CL = CL_TO_max/1.2^2;              % at 1.2 V_stall_TO
LoverD = CL/(CD_0_clean+delta_CD0_TOflaps+delta_CD0_LG+CL^2/(pi*AR*e_TOflaps)); % CL/CD_TO_GearDown
ToverW_TO = 2*(1/LoverD+0.012);    % CGR>0.012
ToverW_TO1 = ToverW_TO/0.8; % 50°F效應(除以0.8)

% FAR25.121 OEI
CL = CL_TO_max/1.1^2; % V_LOF = 1.1 V_stall_TO
LoverD = CL/(CD_0_clean+delta_CD0_TOflaps+delta_CD0_LG+CL^2/(pi*AR*e_TOflaps));
ToverW_TO = 2*(1/LoverD); % CGR>0
ToverW_TO2 = ToverW_TO/0.8; % 50°F效應(除以0.8)

% FAR25.121 OEI
CL = CL_TO_max/1.2^2; % at 1.2 V_stall_TO
LoverD = CL/(CD_0_clean+delta_CD0_TOflaps+CL^2/(pi*AR*e_TOflaps));
ToverW_TO = 2*(1/LoverD+0.024); % CGR>0.024
ToverW_TO3 = ToverW_TO/0.8;

% FAR25.121 OEI
CL_max = 1.4; % From Table 3.1
CL = CL_max/1.25^2; % at 1.25 V_stall
LoverD = CL/(CD_0_clean + CL^2/(pi*AR*0.85));
ToverW_TO = 2*(1/LoverD+0.012); % CGR>0.012
ToverW_TO4 = ToverW_TO/0.94/0.8; % 最大推力校正(除以0.94), 50°F效應(除以0.8)

% FAR25.119 AEO
CL_max_L = 2.8; % From Table 3.1
CL = CL_max_L/1.3^2; % at 1.3 V_stall_L
LoverD = CL/(CD_0_clean+delta_CDO_Lflaps+delta_CD0_LG+CL^2/(pi*AR*e_Lflaps));
ToverW_L = 1/LoverD+0.032; % CGR>0.032
ToverW_TO5 = ToverW_L*(W_L/W_TO)/0.8;

% FAR25.121 OEI
CL_max_A = 2.4; % From Table 3.1
CL = CL_max_A/1.5^2; % at 1.5 V_stall_A
LoverD = CL/((CD_0_clean+delta_CD0_TOflaps+CD_0_clean+delta_CDO_Lflaps)/2+delta_CD0_LG+CL^2/(pi*AR*e_Lflaps));
ToverW_L = 2*(1/LoverD+0.021); % CGR>0.021
ToverW_TO6 = ToverW_L*(W_L/W_TO)/0.8; 

yline(ToverW_TO1,'color',[0.4940 0.1840 0.5560])
yline(ToverW_TO2,'color',[0.4940 0.1840 0.5560])
yline(ToverW_TO3,'color',[0.4940 0.1840 0.5560])
yline(ToverW_TO4,'color',[0.4940 0.1840 0.5560])
yline(ToverW_TO5,'color',[0.4940 0.1840 0.5560])
yline(ToverW_TO6,'color',[0.4940 0.1840 0.5560])

plot(WoverS_TO_wiki,ToverW_TO_wiki,'rx')

title('Matching result for sizing of Boeing 737MAX8')
xlabel('(W/S)_{TO}');
ylabel('(T/W)_{TO}');

hold off

%% FAR25 TAKEOFF DISTANCE SIZING
figure()
hold on
for CL_max_TO = 1.6:0.2:2.2
    ToverW = (0.009640.*WoverS)/CL_max_TO;
    plot(WoverS,ToverW)
end
legend('1.6','1.8','2.0','2.2');
xlabel('W/S');
ylabel('T/W');
title('FAR25 TAKEOFF DISTANCE SIZING');
hold off

%% FAR25 CLIMB RATE SIZING
delta_CD0_TOflaps = 0.015;  % From p.127, Table 3.6
delta_CDO_Lflaps = 0.065;   % From p.127, Table 3.6
delta_CD0_LG = 0.02;        % From p.127, Table 3.6
e_TOflaps = 0.8;            % From p.127, Table 3.6
e_Lflaps = 0.75;            % From p.127, Table 3.6
W_L = W_TO*0.84;            % 0.84 is from p.107, Table 3.3
WoverS = 100;               % from example in the book ( 0:200 )
S = W_TO/WoverS;
CD_0_clean = f_2/S;         % Take cf = 0.003

% CD_clean = CD_0_clean + CL^2/(pi*AR*e)
% CD_TO_GearUp = CD_0_clean+delta_CD0_TOflaps+CL^2/(pi*AR*e_TOflaps)
% CD_TO_GearDown = CD_0_clean+delta_CD0_TOflaps+delta_CD0_LG+CL^2/(pi*AR*e_TOflaps)
% CD_L_GearUp = CD_0_clean+delta_CDO_Lflaps+CL^2/(pi*AR*e_Lflaps)
% CD_L_GearDown = CD_0_clean+delta_CDO_Lflaps+delta_CD0_LG+CL^2/(pi*AR*e_Lflaps)

% FAR25.111 OEI (P.145)
CL_TO_max = 2;                     % From Table 3.1
CL = CL_TO_max/1.2^2;              % at 1.2 V_stall_TO
LoverD = CL/(CD_0_clean+delta_CD0_TOflaps+delta_CD0_LG+CL^2/(pi*AR*e_TOflaps)); % CL/CD_TO_GearDown
ToverW_TO = 2*(1/LoverD+0.012);    % CGR>0.012
ToverW_TO1 = ToverW_TO/0.8; % 50°F效應(除以0.8)
% FAR25.121 OEI
CL = CL_TO_max/1.1^2; % V_LOF = 1.1 V_stall_TO
LoverD = CL/(CD_0_clean+delta_CD0_TOflaps+delta_CD0_LG+CL^2/(pi*AR*e_TOflaps));
ToverW_TO = 2*(1/LoverD); % CGR>0
ToverW_TO2 = ToverW_TO/0.8; % 50°F效應(除以0.8)
% FAR25.121 OEI
CL = CL_TO_max/1.2^2; % at 1.2 V_stall_TO
LoverD = CL/(CD_0_clean+delta_CD0_TOflaps+CL^2/(pi*AR*e_TOflaps));
ToverW_TO = 2*(1/LoverD+0.024); % CGR>0.024
ToverW_TO3 = ToverW_TO/0.8;
% FAR25.121 OEI
CL_max = 1.4; % From Table 3.1
CL = CL_max/1.25^2; % at 1.25 V_stall
LoverD = CL/(CD_0_clean + CL^2/(pi*AR*e));
ToverW_TO = 2*(1/LoverD+0.012); % CGR>0.012
ToverW_TO4 = ToverW_TO/0.94/0.8; % 最大推力校正(除以0.94), 50°F效應(除以0.8)
% FAR25.119 AEO
CL_max_L = 2.8; % From Table 3.1
CL = CL_max_L/1.3^2; % at 1.3 V_stall_L
LoverD = CL/(CD_0_clean+delta_CDO_Lflaps+delta_CD0_LG+CL^2/(pi*AR*e_Lflaps));
ToverW_L = 1/LoverD+0.032; % CGR>0.032
ToverW_TO5 = ToverW_L*(W_L/W_TO)/0.8;
% FAR25.121 OEI
CL_max_A = 2.4; % From Table 3.1
CL = CL_max_A/1.5^2; % at 1.5 V_stall_A
LoverD = CL/((CD_0_clean+delta_CD0_TOflaps+CD_0_clean+delta_CDO_Lflaps)/2+delta_CD0_LG+CL^2/(pi*AR*e_Lflaps));
ToverW_L = 2*(1/LoverD+0.021); % CGR>0.021
ToverW_TO6 = ToverW_L*(W_L/W_TO)/0.8
WoverS = 0:10:200;
hold on
yline(ToverW_TO1)
yline(ToverW_TO2)
yline(ToverW_TO3)
yline(ToverW_TO4)
yline(ToverW_TO5)
yline(ToverW_TO6)
hold off

% Table 3.6 First Estimates for ΔCD_0 and e With Flaps and Gear Down
% Configuration     ΔCD_0           e
% Clean             0               0.80 - 0.85
% Take-off flaps    0.010 - 0.020   0.75 - 0.80
% Landing flaps     0.055 - 0.075   0.70 - 0.75
% Landing Gear      0.015 - 0.025   no effect

%%
function [a,rho]=Standard_Atmosphere(h)