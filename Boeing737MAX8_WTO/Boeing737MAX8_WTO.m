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
    Boeing737MAX8_WTO_WorkspaceSavedDirectory = 'Boeing737MAX8_WTO/Boeing737MAX8_WTO.mat';

elseif strcmp(OS,'Windows') 
    Boeing737MAX8_Mission_Profile_setup_WorkspaceSavedDirectory = 'G:\飛設\Boeing737MAX8-Research\Boeing737MAX8_Mission_Profile_setup\Boeing737MAX8_Mission_Profile_setup.mat';
    Boeing737MAX8_WTO_WorkspaceSavedDirectory = 'G:\飛設\Boeing737MAX8-Research\Boeing737MAX8_WTO\Boeing737MAX8_WTO.mat';

else
    error;
end

load(Boeing737MAX8_Mission_Profile_setup_WorkspaceSavedDirectory);
%% Start time record
%
time = now;
date = datetime(time,'ConvertFrom','datenum');
disp('----------------------------------------------------')
string_StartTime=[' StartTime: ' ,datestr(date)];
disp(string_StartTime)
disp('----------------------------------------------------')

%% Numerical approximation of W_TO_guess
% ResultMatrixApporx sizing
W_TO_Apporx_row = InputParametersMatrix_row;
W_TO_Approx_column = 12;
W_TO_Approx = zeros(W_TO_Apporx_row,W_TO_Approx_column);

parfor row = 1:W_TO_Apporx_row
    % Temporary matrix for parallel computing
    temp = zeros(1,Result_column);

    % Read data
    CruiseAltitude = InputParametersMatrix(row,1);
    Range = InputParametersMatrix(row,2);
    LoverD_Cruise = InputParametersMatrix(row,3);
    LoverD_Loiter = InputParametersMatrix(row,4);
    c_j_cruise  = InputParametersMatrix(row,5);
    c_j_loiter = InputParametersMatrix(row,6);
    M_ff = InputParametersMatrix(row,13);
    C = InputParametersMatrix(row,14);
  
    for W_TO_guess = 165991:182200 % W_TO_guess_LowerBound:W_TO_wiki
        W_E_real = 10^((log10(W_TO_guess)-A)/B);
        W_E_tent = C*W_TO_guess - D;
        error = abs(W_E_tent - W_E_real)/ W_E_real;
        if error < 0.005
            W_E_error =  abs(W_E_real - W_E_wiki)/W_E_real;

            % Output data
            temp(1) = CruiseAltitude;
            temp(2) = Range;
            temp(3) = LoverD_Cruise;
            temp(4) = LoverD_Loiter;
            temp(5) = c_j_cruise;
            temp(6) = c_j_loiter;
            temp(7) = M_ff;
            temp(8) = C;
            temp(9) = W_TO_guess;
            temp(10) = W_E_tent;
            temp(11) = W_E_real;
            temp(12) = W_E_error;

            % Output calculate result into Result matrix
            W_TO_Approx(row, :) = temp;
            break
        end
    end
end

save(Boeing737MAX8_WTO_WorkspaceSavedDirectory);
time = now;
date = datetime(time,'ConvertFrom','datenum');
string_RecordTime=[' RecordTime: ',datestr(date)];
string_Workspace_saved=[' Workspace is saved'];
disp('----------------------------------------------------')
disp(string_Workspace_saved)
disp(string_RecordTime)

%% Check the amount of numerical approximation solutions
n = 0;
% ResultMatrixApporxSolutions sizing
W_TO_ApproxSolutions = zeros(W_TO_Approx_column);

for row = 1:W_TO_Apporx_row
    if W_TO_Approx(row,8) > 0 && W_TO_Approx(row,8) < W_TO_wiki
        if W_TO_Approx(row,12) < 0.005 | W_TO_Approx(row,11) < W_E_wiki
            n = n+1;
            W_TO_ApproxSolutions(n,:) = W_TO_Approx(row,:);
        end
    end
end

disp('----------------------------------------------------')
string_solutions=[' There are ',num2str(n),' numerical approximation of W_TO_guess less than W_TO_wiki.'];
disp(string_solutions)

%% Numerical solution of W_TO_guess
% ResultMatrix sizing
W_TO_row = height(W_TO_ApproxSolutions);
W_TO_column = width(W_TO_ApproxSolutions) ;
W_TO = zeros(W_TO_row,W_TO_column);

x1 = sym('x1', [1,W_TO_row]);
parfor row = 1:height(W_TO_ApproxSolutions)
    % Temporary matrix for parallel computing
    temp = zeros(1,Result_column);
   
            % Read data
            CruiseAltitude = W_TO_ApproxSolutions(row,1);
            Range = W_TO_ApproxSolutions(row,2);
            LoverD_Cruise = W_TO_ApproxSolutions(row,3);
            LoverD_Loiter = W_TO_ApproxSolutions(row,4);
            c_j_cruise  = W_TO_ApproxSolutions(row,5);
            c_j_loiter = W_TO_ApproxSolutions(row,6);
            M_ff = W_TO_ApproxSolutions(row,7);
            C = W_TO_ApproxSolutions(row,8);
    
            % vpasolve
            W_TO_guess = vpasolve( A + B*log10(C*x1(row) - D) - log10(x1(row)) == 0 );
            
            % Computing
            W_E_real = 10^((log10(W_TO_guess)-A)/B);
            W_E_tent = C*W_TO_guess-D;
            W_E_error =  abs(W_E_real - W_E_wiki)/W_E_real;
    
            % Output data
            temp(1) = CruiseAltitude;
            temp(2) = Range;
            temp(3) = LoverD_Cruise;
            temp(4) = LoverD_Loiter;
            temp(5) = c_j_cruise;
            temp(6) = c_j_loiter;
            temp(7) = M_ff;
            temp(8) = C;
            temp(9) = W_TO_guess;
            temp(10) = W_E_tent;
            temp(11) = W_E_real;
            temp(12) = W_E_error;

            % Output calculate result into Result matrix
            W_TO(row, :) = temp;
end

% section saved
save(WorkspaceSavedDirectory);
time = now;
date = datetime(time,'ConvertFrom','datenum');
string_RecordTime4=[' RecordTime: ',datestr(date)];
string_Workspace_saved4=[' Workspace is saved 4'];
disp('----------------------------------------------------')
disp(string_Workspace_saved4)
disp(string_RecordTime4)

%% Sensitivity section
%
m = 0;

% Open TransportJet_WTO_CheatingVersion_result.txt
fid = fopen(SelectedResultOutputDirectory,'wt');
%
for row = 1:1:height(W_TO)
    if W_TO(row,1) > 0 && W_TO(row,11) < error_wiki/250 % need to be correct e.g. error_wiki = abs(W_E_real_wiki-W_E_wiki)/W_E_real_wiki
        % Read data
        CruiseAltitude = InputParametersMatrix(row,1);
        Range = InputParametersMatrix(row,2);
        LoverD_Cruise = InputParametersMatrix(row,3);
        LoverD_Loiter = InputParametersMatrix(row,4);
        c_j_cruise  = InputParametersMatrix(row,5);
        c_j_loiter = InputParametersMatrix(row,6);
        CruiseSpeed = InputParametersMatrix(row,7);
        M_ff = W_TO(row,7);
        C = InputParametersMatrix(row,14);
        W_TO_guess = W_TO(row,8);

        % Sensitivity calculate
        W_TO = W_TO_guess;
        F=-B*(W_TO^2)*((C*W_TO*(1-B)-D)^-1)*(1+0)*M_ff;
        % W_TO over W_PL
        W_TO_over_W_PL = B*W_TO*(D-C*(1-B)*W_TO)^-1;
        % W_TO over W_E
        W_TO_over_W_E = B*W_TO*(10^((log10(W_TO)-A)/B))^-1;
        % W_TO over Range
        W_TO_over_Range = F*c_j_cruise*(CruiseSpeed*LoverD_Cruise)^-1;
        % W_TO over Endurance
        W_TO_over_Endurance = F*c_j_loiter*LoverD_Loiter^-1;
        % W_TO over Cruise speed
        W_TO_over_CriuseSpeed = -F*Range*c_j_cruise*(CruiseSpeed^2*LoverD_Cruise)^-1;
        % W_TO over c_j_Range
        W_TO_over_c_j_Range = F*Range*(CruiseSpeed*LoverD_Cruise)^-1;
        % W_TO over L/D_Range
        W_TO_over_LoverD_Range = -F*Range*c_j_cruise*(CruiseSpeed*LoverD_Cruise^2)^-1;
        % W_TO over c_j_Loiter
        W_TO_over_c_j_Loiter = F*Endurance*LoverD_Loiter^-1;
        % W_TO over L/D_Loiter
        W_TO_over_LoverD_Loiter = -F*Endurance*c_j_loiter*LoverD_Loiter^-2;

        % Output sensitivity calculate result
        W_TO(row,12) = W_TO_over_W_PL;
        W_TO(row,13) = W_TO_over_W_E ;
        W_TO(row,14) = W_TO_over_Range;
        W_TO(row,15) = W_TO_over_Endurance;
        W_TO(row,16) = W_TO_over_CriuseSpeed;
        W_TO(row,17) = W_TO_over_c_j_Range;
        W_TO(row,18) = W_TO_over_LoverD_Range;
        W_TO(row,19) = W_TO_over_c_j_Loiter;
        W_TO(row,20) = W_TO_over_LoverD_Loiter;

        % W_TO_guess Iteration Result
        string_CruiseAltitude=[' 1.CruiseAltitude = ',num2str(W_TO(row,1)),' ft'];
        string_Range=[' 2.Range = ',num2str(W_TO(row,2)),' nm'];
        string_LoverD_Cruise=[' 3.L/D Cruise = ',num2str(W_TO(row,3))];
        string_LoverD_Loiter=[' 4.L/D Loiter = ',num2str(W_TO(row,4))];
        string_c_j_cruise=[' 5.c_j Cruise = ',num2str(W_TO(row,5))];
        string_c_j_loiter=[' 6.c_j Loiter = ',num2str(W_TO(row,6))];
        string_M_ff=[' 7.M_ff = ',num2str(W_TO(row,7))];
        string_W_TO_guess=[' 8.W_TO_guess = ',num2str(W_TO(row,8)),' lbs'];
        string_W_E_tent=[' 9.W_E_tent = ',num2str(W_TO(row,9)),' lbs'];
        string_W_E_real=[' 10.W_E_real = ',num2str(W_TO(row,10)),' lbs'];
        string_W_E_error=[' 11.W_E_error = ',num2str(W_TO(row,11)*100),' %% (compare with W_E_wiki = W_OE_wiki - W_TO_wiki*0.005 - W_crew)/250'];

        % Sensitivity Result
        string_W_TO_over_W_PL=[' 12.W_TO_over_W_PL = ',num2str(W_TO(row,12))];
        string_W_TO_over_W_E=[' 13.W_TO_over_W_E = ',num2str(W_TO(row,13))];
        string_W_TO_over_Range=[' 14.W_TO_over_Range = ',num2str(W_TO(row,14)),' lbs/nm'];
        string_W_TO_over_Endurance=[' 15.W_TO_over_Endurance = ',num2str(W_TO(row,15)),' lbs/hr'];
        string_W_TO_over_CriuseSpeed=[' 16.W_TO_over_CriuseSpeed = ',num2str(W_TO(row,16)),' lbs/kt'];
        string_W_TO_over_c_j_Range=[' 17.W_TO_over_c_j_Range = ',num2str(W_TO(row,17)),' lbs/lbs/lbs/hr'];
        string_W_TO_over_LoverD_Range=[' 18.W_TO_over_LoverD_Range = ',num2str(W_TO(row,18)),' lbs'];
        string_W_TO_over_LoverD_Range=[' 18.W_TO_over_LoverD_Range = ',num2str(W_TO(row,18)),' lbs'];
        string_W_TO_over_c_j_Loiter=[' 19.W_TO_over_c_j_Loiter = ',num2str(W_TO(row,19)),' lbs/lbs/lbs/hr'];
        string_W_TO_over_LoverD_Loiter=[' 20.W_TO_over_LoverD_Loiter = ',num2str(W_TO(row,20)),' lbs'];

        % Print result in txt
        fprintf(fid,' ----------------------------------------------------' );
        fprintf(fid,'\n');
        fprintf(fid,' W_TO_guess Iteration Result' );
        fprintf(fid,'\n');
        fprintf(fid,string_CruiseAltitude );
        fprintf(fid,'\n');
        fprintf(fid,string_Range );
        fprintf(fid,'\n');
        fprintf(fid,string_LoverD_Cruise );
        fprintf(fid,'\n');
        fprintf(fid,string_LoverD_Loiter );
        fprintf(fid,'\n');
        fprintf(fid,string_c_j_cruise );
        fprintf(fid,'\n');
        fprintf(fid,string_c_j_loiter );
        fprintf(fid,'\n');
        fprintf(fid,string_M_ff );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_guess );
        fprintf(fid,'\n');
        fprintf(fid,string_W_E_tent );
        fprintf(fid,'\n');
        fprintf(fid,string_W_E_real );
        fprintf(fid,'\n');
        fprintf(fid,string_W_E_error );
        fprintf(fid,'\n');
        fprintf(fid,' Sensitivity Result' );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_over_W_PL );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_over_W_E );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_over_Range );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_over_Endurance );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_over_CriuseSpeed );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_over_c_j_Range );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_over_LoverD_Range );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_over_c_j_Loiter );
        fprintf(fid,'\n');
        fprintf(fid,string_W_TO_over_LoverD_Loiter );
        fprintf(fid,'\n');

        m=m+1;
    end
end

disp('----------------------------------------------------')
string_solutions=[' There are ',num2str(m),' solutions printed in txt file.'];
disp(string_solutions)
fprintf(fid,'----------------------------------------------------');
fprintf(fid,'\n');
fprintf(fid,string_solutions );

% Close TransportJet_WTO_CheatingVersion_result.txt
fclose(fid);

% section saved
save(WorkspaceSavedDirectory);
time = now;
date = datetime(time,'ConvertFrom','datenum');
string_RecordTime5 = [' RecordTime: ',datestr(date)];
string_Workspace_saved5 = [' Workspace is saved 5'];
disp('----------------------------------------------------')
disp(string_Workspace_saved5)
disp(string_RecordTime5)

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

