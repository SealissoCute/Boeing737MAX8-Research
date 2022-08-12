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
%% Numerical approximation of W_TO_guess
% ResultMatrixApporx sizing
W_TO_Apporx_row = InputParametersMatrix_row;
W_TO_Approx_column = 12;
W_TO_Approx = zeros(W_TO_Apporx_row,W_TO_Approx_column);

parfor row = 1:1324224 % W_TO_Apporx_row
    % Temporary matrix for parallel computing
    temp = zeros(1,W_TO_Approx_column);

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

for row = 1:1324224 % W_TO_Apporx_row
    if W_TO_Approx(row,9) > 0 && W_TO_Approx(row,9) < W_TO_wiki
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
parfor row = 1: W_TO_row
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
save(Boeing737MAX8_WTO_WorkspaceSavedDirectory);
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
