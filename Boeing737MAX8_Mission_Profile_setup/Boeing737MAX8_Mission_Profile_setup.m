%% Reverse Engineering of Boeing 737 MAX 8
% Clear previous data

format long
clc
clear

%% Setup output file directory (Notice: Need to be customized at different computer)
prompt = ' Which OS are you using? (Windows/MacOS) Ans:';
OS = input(prompt,'s');
if strcmp(OS,'MacOS')
    Boeing737MAX8_Mission_Profile_setup_WorkspaceSavedDirectory = 'Boeing737MAX8_Mission_Profile_setup/Boeing737MAX8_Mission_Profile_setup.mat';

elseif strcmp(OS,'Windows') 
    Boeing737MAX8_Mission_Profile_setup_WorkspaceSavedDirectory = 'G:\飛設\Boeing737MAX8-Research\Boeing737MAX8_Mission_Profile_setup\Boeing737MAX8_Mission_Profile_setup.mat';
    
else
    error;
end
%% Start time record
%
time = now;
date = datetime(time,'ConvertFrom','datenum');
disp('----------------------------------------------------')
string_StartTime=[' StartTime: ' ,datestr(date)];
disp(string_StartTime)
disp('----------------------------------------------------')

%% Unit exchange
%
ft_to_m = 0.3048;       % ft to m
m_s_to_mph = 2.236936;  % m/s to mph
m_s_to_kt = 1.943844;   % m/s to kt
m_s_to_ft_s = 3.280840; % m/s to ft/s
ft_s_to_kt = 0.592484;  % ft/s to kt


%% Mission profile parameters
% Payload weight (unit: lbs)
W_PL = (175+30)*180;        % 180 Passengers at 175 lbs each and 30 lbs of baggage each

% Crew weight (unit: lbs)
W_crew = (175+30)*8;        % 2 Pilots and 6 flight attendents at 175 lbs each and 30 lbs of baggage each

%
D = W_crew + W_PL;

% CruiseAltitude (unit: ft)
CruiseAltitudeMin = 35000;
CruiseAltitudeMax = 41000;
CruiseAltitudeInterval = 1000;
CruiseAltitudeMatrix = [CruiseAltitudeMin:CruiseAltitudeInterval:CruiseAltitudeMax]; %

% Range (unit: nm)
RangeMin = 2000;
RangeMax = 3500;
RangeInterval = 100;
RangeMatrix = [RangeMin:RangeInterval:RangeMax]; %

% LoverD_Cruise
LoverD_CruiseMin = 13;
LoverD_CruiseMax = 14;
LoverD_CruiseInterval = 0.1;
LoverD_CruiseMatrix = [LoverD_CruiseMin:LoverD_CruiseInterval:LoverD_CruiseMax]; %

% LoverD_Loiter
LoverD_LoiterMin = 17;
LoverD_LoiterMax = 18;
LoverD_LoiterInterval = 0.1;
LoverD_LoiterMatrix = [LoverD_LoiterMin:LoverD_LoiterInterval:LoverD_LoiterMax]; %

% c_j_cruise
c_j_cruiseMin = 0.5;
c_j_cruiseMax = 0.55;
c_j_cruiseInterval = 0.01;
c_j_cruiseMatrix = [c_j_cruiseMin:c_j_cruiseInterval:c_j_cruiseMax]; %

% c_j_loiter
c_j_loiterMin = 0.55;
c_j_loiterMax = 0.6;
c_j_loiterInterval =0.01;
c_j_loiterMatrix = [c_j_loiterMin:c_j_loiterInterval:c_j_loiterMax]; %

%
Endurance = 0.5;            % Loiter, unit: hr
AverageClimbRate = 2500;    % unit: fpm
CruiseSpeed_Mach = 0.79;    % unit: Mach
AlternateCruiseSpeed = 250; % unit: kts
AlternateRange = 100;       % unit: nm


% Regression Line Constants A and B of Equation
% Airplane type: Transport jet
A = 0.0833;
B = 1.0383;

% W_TO & W_OE from wikipedia (unit: lbs)
W_OE_wiki = 99360;
W_TO_wiki = 182200;
W_E_wiki = W_OE_wiki - W_TO_wiki*0.005 - W_crew;
W_E_wiki2W_TO_guess = 10^(A+B*log10(W_E_wiki));
W_E_real_wiki = 10^((log10(W_TO_wiki)-A)/B);
error_wiki = abs(W_E_real_wiki-W_E_wiki)/W_E_real_wiki;

%% Fuel Fraction parameters
%
W1_W_TO_guess_ratio = 0.990;                                           % Engine start and warm up, from Table 2.1
W2_W1_ratio = 0.990;                                                   % Taxi, from Table 2.1
W3_W2_ratio = 0.995;                                                   % Take-off, from Table 2.1
W4_W3_ratio = 0.980;                                                   % Climb, from Table 2.1
W7_W6_ratio = 0.990;                                                   % Decent,from Table 2.1
W8_W7_ratio = 1/(exp(AlternateRange/(AlternateCruiseSpeed/0.9*10)));   % Fly to alternate and descend, from Brequet's range equation, L/D = 10, c_j = 0.9
W9_W8_ratio = 0.992;                                                   % Landing, Taxi and Shutdown, From Table 2.1

%% InputParametersMatrix setup section
% InputParametersMatrix sizing
InputParametersMatrix_row = width(CruiseAltitudeMatrix)*width(RangeMatrix)*width(LoverD_CruiseMatrix)...
    *width(LoverD_LoiterMatrix)*width(c_j_cruiseMatrix)*width(c_j_loiterMatrix);
InputParametersMatrix_column = 14;
InputParametersMatrixTemp = zeros(InputParametersMatrix_row,InputParametersMatrix_column);
InputParametersMatrix = zeros(InputParametersMatrix_row,InputParametersMatrix_column);


% Create InputParametersMatrix
n=1;
for CruiseAltitude = CruiseAltitudeMin:CruiseAltitudeInterval:CruiseAltitudeMax
    for Range = RangeMin:RangeInterval:RangeMax
        for LoverD_Cruise = LoverD_CruiseMin:LoverD_CruiseInterval:LoverD_CruiseMax
            for LoverD_Loiter = LoverD_LoiterMin:LoverD_LoiterInterval:LoverD_LoiterMax
                for c_j_cruise = c_j_cruiseMin:c_j_cruiseInterval:c_j_cruiseMax
                    for c_j_loiter = c_j_loiterMin:c_j_loiterInterval:c_j_loiterMax
                        InputParametersMatrixTemp(n,1) = CruiseAltitude;
                        InputParametersMatrixTemp(n,2) = Range;
                        InputParametersMatrixTemp(n,3) = LoverD_Cruise;
                        InputParametersMatrixTemp(n,4) = LoverD_Loiter;
                        InputParametersMatrixTemp(n,5) = c_j_cruise;
                        InputParametersMatrixTemp(n,6) = c_j_loiter;
                        n=n+1;
                    end
                end
            end
        end
    end
    string_CruiseAltitude=[' Parameters at CruiseAltitude = ',num2str(CruiseAltitude),' ft are created'];
    disp(string_CruiseAltitude)
end


%% Parallel computing CruiseSpeed/AverageClimbSpeed/ClimbTime/CruiseRange/W5_W4_ratio/W6_W5_ratio/M_ff
parfor row = 1:InputParametersMatrix_row
    % Temporary matrix for parallel computing
    temp = zeros(1,InputParametersMatrix_column);

    % Read data form InputParametersMatrixtemp
    CruiseAltitude = InputParametersMatrixTemp(row,1);
    Range = InputParametersMatrixTemp(row,2);
    LoverD_Cruise = InputParametersMatrixTemp(row,3);
    LoverD_Loiter = InputParametersMatrixTemp(row,4);
    c_j_cruise  = InputParametersMatrixTemp(row,5);
    c_j_loiter = InputParametersMatrixTemp(row,6);

    % Calculate parameters
    [a]=Standard_Atmosphere(CruiseAltitude);                                     % unit:Imperial system
    CruiseSpeed = CruiseSpeed_Mach*a*ft_s_to_kt;                                 % unit: kts
    AverageClimbSpeed = CruiseSpeed*0.6;                                         % unit: kts
    ClimbTime = CruiseAltitude/AverageClimbRate;                                 % Climb time, unit: minute
    ClimbRange = AverageClimbSpeed*(ClimbTime/60);                               % Climb range, unit: nm
    CruiseRange = Range - ClimbRange;                                            % Cruise range, unit: nm
    W5_W4_ratio = 1/(exp(CruiseRange/(CruiseSpeed/c_j_cruise*LoverD_Cruise)));   % Cruise, from Breguet's range equation
    W6_W5_ratio = 1/(exp(0.5/(1/c_j_loiter*LoverD_Loiter)));                     % Loiter, from Breguet's endurance equation
    M_ff = W1_W_TO_guess_ratio*W2_W1_ratio*W3_W2_ratio*W4_W3_ratio*...
        W5_W4_ratio*W6_W5_ratio*W7_W6_ratio*W8_W7_ratio*W9_W8_ratio
    C = 1-(1-M_ff)-0.005;
 
    % Output data
    temp(1) = CruiseAltitude;
    temp(2) = Range;
    temp(3) = LoverD_Cruise;
    temp(4) = LoverD_Loiter;
    temp(5) = c_j_cruise;
    temp(6) = c_j_loiter;
    temp(7) = CruiseSpeed;
    temp(8) = AverageClimbSpeed;
    temp(9) = ClimbTime;
    temp(10) = ClimbRange;
    temp(11) = CruiseRange;
    temp(12) = W5_W4_ratio;
    temp(12) = W6_W5_ratio;
    temp(13) = M_ff;
    temp(14) = C;

    % Output calculate result into InputParametersMatrix
    InputParametersMatrix(row, :) = temp;

end

%% Plot W_E_real/W_E_tent_min/W_E_tent_max
%
C_min = min(InputParametersMatrix(:,14));
C_max = max(InputParametersMatrix(:,14));

%
C_min = min(InputParametersMatrix(:,14));
C_max = max(InputParametersMatrix(:,14));

%
x_W_TO_guess = 0:100:500000;
y_W_E_real = 10.^((log10(x_W_TO_guess)-A)/B);
y_W_E_tent_min = C_min.*x_W_TO_guess - D;
y_W_E_tent_max = C_max.*x_W_TO_guess - D;

% Find the lower and upper bound of W_To_guess
syms x
W_TO_guess_min = vpasolve( A + B*log10(C_max*x - D) - log10(x) == 0 );
W_TO_guess_max = vpasolve( A + B*log10(C_min*x - D) - log10(x) == 0 );
W_TO_guess_LowerBound = floor(W_TO_guess_min);
W_to_guess_UpperBound = ceil(W_TO_guess_max);
%
x1 = [182044 182044];
y1 = [-0.5*10^5 96809];
x2 = [182200 182200];
y2 = [-0.5*10^5 96889];
x3 = [0 182044];
y3 = [96809 96809];
x4 = [0 182200];
y4 = [96889 96889];

hold on
plot(x_W_TO_guess,y_W_E_real)
plot(x_W_TO_guess,y_W_E_tent_min)
plot(x_W_TO_guess,y_W_E_tent_max)
plot(x1,y1,'--b')
plot(x2,y2,'--m')
plot(x3,y3,'--b')
plot(x4,y4,'--m')
line1 = xline(double(W_TO_guess_min),'--');
line2 = xline(double(W_TO_guess_max),'--');
line1.LabelVerticalAlignment = 'bottom';
line2.LabelVerticalAlignment = 'bottom';
xlabel('W_T_Oguess');
ylabel('W_E');
legend('W_Ereal','W_Etent_m_i_n','W_Etent_m_a_x');
hold off

% section saved
save(Boeing737MAX8_Mission_Profile_setup_WorkspaceSavedDirectory);
time = now;
date = datetime(time,'ConvertFrom','datenum');
string_RecordTime=[' RecordTime: ',datestr(date)];
string_Workspace_saved=[' Workspace is saved'];
disp('----------------------------------------------------')
disp(string_Workspace_saved)
disp(string_RecordTime)

%%
function [a]=Standard_Atmosphere(h)
% Standard Atmosphere (SI Units)
% [C,a,P,rho,g,mu]=Standard_Atmosphere(h)
%
% bibliography :
% [1] Yunus A. Cengel & John M.Cimbala°ßFLUID MECHANICS ...
% Fundamentals and  Applications°®, McGraw-Hill., p897.
% [2] John J. Bertin & Russell M. Cummings°ßAERODYNAMICS FOR ENGINEERS°®,
% 5th Edition, Pearson Education International., p21-p43.
% [3] WARREN F.PHILLIPS,°ßMECHANICS of FLIGHT°®, 2nd Edition, John Wiley.
% p10-p14.
% [4] John D.Anderson,"Modern Compressibe Flow", third Edition, McGraw-Hill
% ., p585-p613.
%
% input arguments:
%    h = Geometric altitude. (default : sea level)
%
% output arguments:
%    Rankine = The temperature in Rankine scale.
%    a = Speed of sound.
%    P = The standard atmosphere at h.
%    rho = Density.
%    g = Is the gravitational acceleration at height h above sea level.
%    mu = Coefficient of viscosity.
%
% example :
%    % Plot C v.s geometrix altitude and P v.s geometrix altitude.
%    >> [C,a,P,rho,g,mu]=Standard_Atmosphere(100:100:90000);
%    >> figure, subplot(1,2,1);
%    >> plot(C,100:100:90000);
%    >> set(gca,'XTick',-100:20:20,'YTick',0:10000:100000,...
%            'YTickLabel',0:10:100,'DataAspectRatio',[1 650 1]);
%    >> title('Standard Atmosphere'); grid on;
%    >> xlabel('Temperature (Celsius)'); ylabel('Geometrix Altitude (Km)');
%    >> subplot(1,2,2);
%    >> plot(P,100:100:90000);
%    >> set(gca,'XTick',0:50000:150000,'YTick',0:10000:100000,...
%            'XTickLabel',0:50:150,'YTickLabel',0:10:100,...
%            'DataAspectRatio',[1 .65 1]);
%    >> title('Standard Atmosphere'); grid on;
%    >> xlabel('Pressure (kPa)'); ylabel('Geometrix Altitude (Km)');

% Last change: 2022/07/25 14:00 pm

% Unit exchange
ft_to_m = 0.3048;       % ft to m
m_s_to_mph = 2.236936;  % m/s to mph
m_s_to_kt = 1.943844;   % m/s to kt
m_s_to_ft_s = 3.280840; % m/s to ft/s
ft_s_to_kt = 0.592484;  % ft/s to kt

%
h = h*ft_to_m; % unit:m

% default values
if ~exist('h','var'), h = 0; end;  % sea level

% Gravitational acceleration
% go = Is the standard gravitational acceleration.
% re = Is the Earth's mean radius.
go = 9.806645;           % m/s^2
re = 6356766;            % m
g = go*(re./(re+h)).^2;  % m/s^2

% Geopotential Altitude
% Z = Geopotential altitude.
% Zi = Is the minimum geopotential altitude in the range.
Z = (re*h)./(re+h);                                        % m
Zi = [0 11000 20000 32000 47000 52000 61000 79000 90000];  % m
[Zig Zg] = meshgrid(Zi,Z);

% Temperature in Celsius
% B = The lapse rate (Temperature Gradient);in SI units [K/m].
% Bi = The lapse rate (Temperature Gradient) for the range;
% Ti = Inital Temperature (absolute) inthe range ; K
% To = Is the sea_level temperature (absolute).
% T = Is temperature in K.
Bi = [-0.0065,0,0.001,0.0028,0,-0.002,-0.004,0];
Ti = [288.150,216.650,216.650,228.650,270.650,270.650,252.650,180.650];
Big = meshgrid(Bi,Z);
Tig = meshgrid(Ti,Z);
B = Big(Zig(:,1:8) <= Zg(:,1:8) &  Zg(:,2:9) < Zig(:,2:9));     % K/m
To = Tig(Zig(:,1:8) <= Zg(:,1:8) &  Zg(:,2:9) < Zig(:,2:9));
T = To'+B'.*...
    (Z-(Zig(Zig(:,1:8) <= Zg(:,1:8) &  Zg(:,2:9) < Zig(:,2:9)))');
Rankine = T*1.8; % unit:Rankine
% standard atmosphere pressure
% R = The gas constant for air.
% The gas constant for air in SI units [(N*m)/(kg*K)].
R = 287.0528;
[~,n] = max(...
    double(Zig(:,1:7) <= Zg(:,1:7) &  Zg(:,2:8) < Zig(:,2:8)),[],2);
N = 0;
P = zeros(size(Z));
for n = n',
    N = N+1;
    if B(N)==0,
        P(N) = Pressure(n+1,go,R,Zi,Ti,Bi)*...
            exp((-go*(Z(N)-Zi(n)))/(R*Ti(n)));
    else
        P(N) = Pressure(n+1,go,R,Zi,Ti,Bi)*...
            ((T(N))/Ti(n))^(-go/(R*Bi(n)));
    end
end

% Recursion function.
    function Pi = Pressure(n,go,R,Zi,Ti,Bi)
        n = n-1;
        if (n > 1)
            if Bi(n-1)==0,
                Pi = Pressure(n,go,R,Zi,Ti,Bi)*...
                    exp((-go*(Zi(n-1+1)-Zi(n-1)))/(R*Ti(n-1)));
            else
                Pi = Pressure(n,go,R,Zi,Ti,Bi)*...
                    ((Ti(n-1)+Bi(n-1)*(Zi(n-1+1)-Zi(n-1)))/...
                    Ti(n-1))^(-go/(R*Bi(n-1)));
            end
        else
            % Standard atmosphere pressure at sea level.
            Pi = 1.01325*10^5;
        end
    end

% Density
rho = P./(R.*T);

% Coefficient of viscosity
if nargout > 5
    % Viscosity in SI units [kg/(s*m)].
    mu = 1.458*(10^-6)*((T.^1.5)./(T+110.4));
end

% molecular energy.
% Cv = Constant volime ; Cv = e(internal energy)/T.
% Cp = Constant pressure ; Cp = h(enthalpy)/T.
% e = Internal energy.
% e_tr = Translational energy.
% e_rot = Rotational energy.
% e_vib = Vibrational energy.
e_tr = (3/2).*R.*T;
e_rot = R.*T;
e_vib = (1/2).*R.*T;
if ( T >= 600 ) % when the air temperature reaches 600K or higher .
    e = e_tr+e_rot+e_vib;
    Cv = e./T;
    Cp = (e+R.*T)./T;
else
    e = e_tr+e_rot;
    Cv = e./T;
    Cp = (e+R.*T)./T;
end

% gamma = Define Cp(constant pressure)/Cv(constant volime).
gamma = Cp./Cv;

% Speed of sound .
a = sqrt(gamma.*R.*T);
a = a*m_s_to_ft_s; % unit:ft/s
end

%% Appendix

% Table 2.1 Suggested Fuel-Fractions For Several Mission Phses
% 1.Engine Start,Warm-up;   2.Taxi;   3.Take-off;   4.Climb;       5.Descent;    6.Landing Taxi,Shutdown
%   0.988,                    0.988,    0.988,        0.995,         0.995,        0.995   Homebuilt
%   0.995,                    0.997,    0.998,        0.992,         0.993,        0.993   Single Engine
%   0.992,                    0.996,    0.966,        0.990,         0.992,        0.992   Twin Engine
%   0.996,                    0.995,    0.996,        0.998,         0.999,        0.998   Agricultural
%   0.990,                    0.995,    0.995,        0.980,         0.990,        0.992   Business Jets
%   0.990,                    0.995,    0.995,        0.985,         0.985,        0.995   Regional TBP's
%   0.990,                    0.990,    0.995,        0.980,         0.990,        0.992   Transport Jets

% Table 2.2 Suggested Values For L/D, c_j, η_p And For c_p For Several Mission Phases
% 1.Cruise - L/D,         c_j,       c_p,       η_p;    2.Loiter - L/D,         c_j,       c_p,       η_p
%            8.0-10.0,               0.6-0.8,   0.7,               10.0-12.0,              0.5-0.7,   0.6    Homebuilt
%            8.0-10.0,               0.5-0.7,   0.8,               10.0-12.0,              0.5-0.7,   0.7    Single Engine
%            8.0-10.0,               0.5-0.7,   0.82,              9.0-11.0,               0.5-0.7,   0.72   Twin Engin
%            5.0-7.0,                0.5-0.7,   0.82,              8.0-10.0,               0.5-0.7,   0.72   Agricultural
%            10.0-12.0,   0.5-0.9,                                 12.0-14.0,   0.4-0.6,                     Business Jets
%            11.0-13.0,              0.4-0.6,   0.85,              14.0-16.0,              0.5-0.7,   0.77   Regional TBP's
%            13.0-15.0,   0.5-0.9,                                 14.0-18.0,   0.4-0.6,                     Transport Jets

% Table 2.15 Regression Line Constants A and B of Equation

