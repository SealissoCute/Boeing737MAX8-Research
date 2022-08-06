%% m-file of Example LO#4
%
clear
%%
 rad2deg = 90./asin(1);
 mph2ftpersec = 5280./3600;
% 
 weight = 2700; % airplaen wieght lb
%
% Flight condition
 v_0 = 120;     % mile per hour
 v_0 = v_0*mph2ftpersec; % ft/sec
 h_0 = 0;       % sea level
 [C,a,P,rho,g,mu]=Standard_Atmosphere(h_0); % in MKS system
 rho = rho * 0.0624279/32; % slug/cubic feet
%
% Give the wing and tail characteristics:
%
 S_w = 180;           % wing surface area, sq_ft
 b_w = 33;            % wing span, ft
 AR_w = b_w^2 /S_w;   % wing aspect ratio
 CL_alpha_w = 4.44;   % wing lfit slope per radian
 alpha_CL0_w = -2.20; % wing zero lift AOA (deg)
 lambda_w = 0.40;     % wing taper ratio
 Lambda_w = 10;       % wing quater-chord swept (deg)
 Cm_ac_w = -0.053;    % wing A.C. moment coefficient 
%
 l_t2l_w = 15;        % distance from tail a.c. to wing a.c., ft
 X = l_t2l_w;
 Y =3.4/(b_w/2);      % the horizontal tail were mounted 3.4 ft above 
                     % the aerodynamic center of wing
 Z = 0;
%
%
AR = AR_w;
RT = lambda_w;
CL = weight/(0.5*rho*v_0^2*S_w);
bw = b_w;
Planform = [];
Lambda_c4 = Lambda_w;
[Kv,Kb,Kp,Ks]=PLL_Downwash(AR,RT,CL,bw,X,Y,Z,Planform,Lambda_c4);
disp([Kv,Kb,Kp,Ks])
%
function [C,a,P,rho,g,mu]=Standard_Atmosphere(h)
% =========================== Separator ==================================
% Graduate Student : Boming, Lin
% Department of Aerospace Engineering
% Tamkang University
% E-mail : 601430092@s01.tku.edu.tw
% Last change: 2013/07/24 12:08 am
% =========================== Separator ==================================
% Standard Atmosphere, SI Units, Revised for Version 1.2.
%     The first input arguments are required;
%     the another have default values.
% =========================== Separator ==================================
% Symbol:
% a = Speed of sound.
% B = The lapse rate (Temperature Gradient);in SI units [K/m].
% Bi = The lapse rate (Temperature Gradient) for the range;
% C = The temperature in Celsius.
% Cv = Constant volime ; Cv = e(internal energy)/T.
% Cp = Constant pressure ; Cp = h(enthalpy)/T.
% e = Internal energy.
% e_tr = Translational energy.
% e_rot = Rotational energy.
% e_vib = Vibrational energy.
% g = Is the gravitational acceleration at height h above sea level.
% go = Is the standard gravitational acceleration.
% gamma = Define Cp(constant pressure)/Cv(constant volime).
% h = Geometric altitude.
% mu = Coefficient of viscosity.
% P = The standard atmosphere at h.
% re = Is the Earth's mean radius.
% R = The gas constant for air.
% T = Is temperature in K.
% Ti = Inital Temperature (absolute) inthe range ; K
% To = Is the sea_level temperature (absolute).
% rho = Density.
% Z = Geopotential altitude.
% Zi = Is the minimum geopotential altitude in the range.
% =========================== Separator ==================================
% ------------------- Gravitational acceleration -------------------------
go=9.806645;     % m/s^2
re=6356766;      % m
g=go*(re/(re+h))^2;
% ---------------------- Geopotential Altitude ---------------------------
% Geopotential Altitude :
Z=(re*h)/(re+h); % m
Zi=[0 11000 20000 32000 47000 52000 61000 79000 90000];  % m
% ---------------------- Temperature in Celsius --------------------------
% in SI units [K/m].
Bi=[-0.0065,0,0.001,0.0028,0,-0.002,-0.004,0];
Ti=[288.150,216.650,216.650,228.650,270.650,270.650,252.650,180.650];
N=0;
for n=1:8
    N=N+1;
    if Zi(n)<=Z && Z<Zi(n+1)
        B=Bi(n);   % K/m 
        To=Ti(n);  % K
        break
    end
end
T=To+B*(Z-Zi(N));
C=T-273.15;
% ----------------------- Standard Atmosphere ----------------------------
% The gas constant for air in SI units [(N*m)/(kg*K)].
R=287.0528;
% Density in SI units [kg/m^3].
Pi=zeros(1,N);
Pi(1)=1.01325*10^5;
for i=1:N-1
    if Bi(i)==0
        Pi(i+1)=Pi(i)*exp((-go*(Zi(i+1)-Zi(i)))/(R*Ti(i)));
    else
        Pi(i+1)=...
            Pi(i)*((Ti(i)+Bi(i)*(Zi(i+1)-Zi(i)))/Ti(i))^(-go/(R*Bi(i)));
    end
end
P=Pi(N)*((To+B*(Z-Zi(N)))/To)^(-go/(R*B));
% ----------------------------- Density ----------------------------------
rho=P/(R*T);
% -------------------- Coefficient of viscosity---------------------------
if nargout > 5
    % Viscosity in SI units [kg/(s*m)].
    mu=1.458*(10^-6)*((T^1.5)/(T+110.4));
end
% -------------------------- Speed of sound ------------------------------
% molecular energy.
e_tr=(3/2)*R*T;
e_rot=R*T;
e_vib=(1/2)*R*T;
if ( T >= 600 ) % when the air temperature reaches 600K or higher .
    e=e_tr+e_rot+e_vib;
    Cv=e/T;
    Cp=(e+R*T)/T;
else
    e=e_tr+e_rot;
    Cv=e/T;
    Cp=(e+R*T)/T;
end
gamma=Cp/Cv;
% Speed of sound .
a=sqrt(gamma*R*T);
% =========================== Separator ==================================
% bibliography :
% [1] Yunus A. Cengel & John M.Cimbala¡§FLUID MECHANICS ...
% Fundamentals and  Applications¡¨, McGraw-Hill., p897.
% [2] John J. Bertin & Russell M. Cummings¡§AERODYNAMICS FOR ENGINEERS¡¨,
% 5th Edition, Pearson Education International., p21-p43.
% [3] WARREN F.PHILLIPS,¡§MECHANICS of FLIGHT¡¨, 2nd Edition, John Wiley.
% p10-p14.
% [4] John D.Anderson,"Modern Compressibe Flow", third Edition, McGraw-Hill
% ., p585-p613.
% =========================== Separator ==================================
end
%
function [Kv,Kb,Kp,Ks,Kd]=PLL_Downwash... 
    (AR,RT,CL,bw,X,Y,Z,Planform,Lambda_c4,BIG_omeaga,Cl_aoa_bar,N)
% Estimating the Low-Speed Downwash Angle on an Aft Tail.
% [Kv,Kb,Kp,Ks,Kd]=PLL_Downwash...
%    (AR,RT,CL,bw,X,Y,Z,Planform,Lambda_c4,BIG_omeaga,Cl_aoa_bar,N)
%
% bibliography:
% [1] WARREN F.PHILLIPS,¡§MECHANICS of FLIGHT¡¨, 2nd Edition, John Wiley.
% 4.5 Estimating the Downwash Angle on an Aft Tail. p411-p421.
% [2] Phillips, W.F., Anderson, E.A., Jenkins, J.C., and Sunouchi, S.,
% (2002), ¡§Estimating the Low-Speed Downwash Angle on an Aft Tail¡¨
% Journa; of Aircraft, 40, 6.
% [3] McCormick, B.W., (1995), Aerodynamics, Aeronautics, and Flight
% Mechanics, 2nd ed., Wiley, New York.
%
%    The first six input arguments are required;
%    the other six have default values.
% 
% input arguments:
%    AR = Aspect ratio.
%    RT = Taper ratio.
%    CL = Lift coefficient.
%    bw = Span of the main wing (m or ft). 
%    X = The axial coordinate in the direction of the freestream (m or ft). 
%    Y = The upward normal coordinate (m or ft). 
%    Z = The spanwise coordinate(m or ft). (default Z = 0)
%    Planform = Planform of the main wing. 
%               (default Planform = 'wings with linear taper')
%    % Planform can be change :
%    % 'wings with linear taper' or 'rectangular wing' or 'elliptic wing'
%    Lambda_c4 = Quarter-chord sweep angle (degrees). 
%                (default Lambda_c4 = 0 deg)
%    BIG_omeaga = Maximum total washout (degrees). 
%                 (default BIG_omeaga = 0 deg)
%    Cl_aoa_bar = Airfoil lift slope. 
%                 (default an airfoil lift slope of 2*pi).    
%    N = Number of nonzero terms in the truncated Fourier series.
%        (default N = 99).
% 
% output arguments:
%    Kv = vortex strength factor in downwash computations.
%    Kb = vortex span factor in downwash coputation.
%    Kp = position factor in downwash computation.
%    Ks = wing sweep factor in downwash computations.
%    Kd = proportionality constant.
%
% example : 
%    % The wing has a AR = 6.05 ,RT = 0.4 and Quarter-chord sweep of 
%    % 10 (degrees). If the horizontal tail were mounted 3.4 feet above the
%    % aerodynamic center of the wing. Assume an airfoil lift slope of 2*pi.
%    % For this wing-tail combination, we have
%    % bw = 33 ft, CL_alpha = 4.44, l_h-l_w = 15 ft, CL = 0.4075.
%    % (The lift on the wing for 120 mph at sea level is 2,700 lbf.
%    % CL = 2700/(0.5*0.0023769*(120*5280/3600)^2*180) = 0.4075).
%    % (where CL_alpha is lift slope for the main wing, l_h is distance aft 
%    % of the center of gravity to aerodynamic center of either an aft
%    % horizontal tail and l_w is distance aft of the center of gravity to
%    % aerodynamic center of the main wing.)
%    >> AR = 6.05; RT = 0.4; bw = 33; CL_alpha = 4.44; X = 15; Y = 3.4;
%    >> Z = 0; Planform = []; Lambda_c4 = 10; CL = 0.4075;
%    >> [Kv,Kb,Kp,Ks]=PLL_Downwash...
%    (AR,RT,CL,bw,X,Y,Z,Planform,Lambda_c4);
%    >> disp([Kv,Kb,Kp,Ks]);
%
%    % The downwash at the position of the horizontal tail is
%    >> downwash = Kv*Kp*Ks/Kb*CL/AR;
%    >> disp(downwash);  % rad
%
%    % The downwash gradient at the tail is
%    >> downwash_gradient = Kv*Kp*Ks/Kb*CL_alpha/AR;
%    >> disp(downwash_gradient);
%
%    See also PLL_Sidewash

% Author : Po-Ming Lin
% Compiler : MATLAB Ver 7.12.0.635 (R2011a).
% edition : Revised for Version 1.3.2
% E-mail : 601430092@s01.tku.edu.tw
%
% Last change: 2014/04/15 09:55 am

% Check input arguments.
if isempty(CL) == 1, CL=1; end
if nargin<6, error('At least 6 input arguments required'), end
if nargin<7 || isempty(Z), Z=0; end
if nargin<8 || isempty(Planform), Planform='wings with linear taper'; end
if nargin<9 || isempty(Lambda_c4), Lambda_c4=0; end
if nargin<10 || isempty(BIG_omeaga),BIG_omeaga=0; end
if nargin<11 || isempty(Cl_aoa_bar), Cl_aoa_bar=2*pi; end
if nargin<12 || isempty(N), N=99; end    % N can be change
if AR<=0, error('Aspect ratio cannot be below zero or equal zero.'), end
if RT<0, error('Taper ratio cannot be below zero.'), end
if N<=0, error('Number of nonzero terms in the truncated.'), end
if ~(RT<=1 || RT==0), error('Taper ratio not applicable'), end
if Cl_aoa_bar<=0, error('Airfoil lift slope not applicable'), end

% theta => theta ; theta= [(i-1)*pi/(N-1)] ; i = 1,N
% c(theta)/b
[i,j] = meshgrid(2:(N-1),1:N);
theta=((i'-1)*pi)/(N-1);
c_of_theta_over_b=(2/(AR*(1+RT)))*(1-(1-RT)*abs(cos(theta)));

% C(i,j) = [4b/(Cl,a_bar*c(theta))+j/sin(theta)]*sin(j*theta)
% However,applying l'Hospital's rule gives.
C(1,:) = j(:,1).^2;
C(2:(N-1),1:N) = (((1./c_of_theta_over_b).*(4/Cl_aoa_bar))+...
               (j'./sin(theta))).*sin(j'.*theta);
C(N,:)=((-1).^(j(:,1)+1)).*(j(:,1).^2);

% where the Fourier coefficients , a(n) are obtained from :
% [C]{a}={1}
[a]=C\ones(1,N)';

n=1:N;
if BIG_omeaga == 0,
    % Fourier coefficients,An
    [A]=((a(n)./(pi*AR*a(1))).*CL)';
else
    % optimum washout distribution function.
    i=1:N;
    theta=((i-1)*pi)/(N-1);
    omeaga=abs(cos(theta));
    % where the Fourier coefficients , b(n) are obtained from :
    % {b}=[C]^-1{£sopt}
    [b]=C\omeaga';

    % Fourier coefficients,An
    [A]=((((a(n)*b(1))./a(1))-b(n)).*deg2rad(BIG_omeaga)+...
        (a(n)./(pi*AR*a(1))).*CL)';
end

% dimensionless aftward axial coordinate, X_bar
X_bar=X/(bw/2);
% dimensionless upward normal coordinate, Y_bar
Y_bar=Y/(bw/2);
% dimensionless leftward aerodynamic spanwise coordinate, Z_bar
Z_bar=Z/(bw/2);

switch lower(Planform)
    case {'wings with linear taper','rectangular wing'}
        n=2:N;
        % vortex strength factor in downwash computations, Kv
        Kv=1+sum((A(n).*sin(n.*pi*0.5))./A(1));
        % vortex span factor in downwash coputation, Kb
        Kb_N=sum((n.*A(n).*cos(n.*pi*0.5))./(((n.^2)-1)*A(1)));
        Kb_D=sum((A(n).*sin(n.*pi*0.5))./A(1));
        Kb=((pi/4)+Kb_N)/(1+Kb_D);
    case 'elliptic wing'
        Kv=1.0;
        Kb=pi/4;
end

% position factor in downwash computation, Kp
if X == inf
    Kp=4/(pi^2);
else
    Kp=((2*(Kb^2))/((pi^2)*((Y_bar^2)+(Kb^2))))*...
    (1+(X_bar*((X_bar^2)+2*(Y_bar^2)+(Kb^2)))/(((X_bar^2)+(Y_bar^2))*...
    sqrt((X_bar^2)+(Y_bar^2)+(Kb^2))));
end

if nargout>3 
    % dimensionless variable, r_bar.
    r_bar=norm([X_bar,Y_bar]); 
    % dimensionless variable, s_bar.
    s_bar=Kb*tan(deg2rad(Lambda_c4));
    % dimensionless variable, t_bar.
    t_bar=norm([(X_bar-s_bar),Y_bar,Kb]);
    % dimensionless variable, to_bar.
    to_bar=norm([X_bar,Y_bar,Kb]);
    % wing sweep factor in downwash computations, Ks
    if X == inf
    Ks=1.0;
    else
    Ks=(1+(X_bar-s_bar)/(t_bar)+(X_bar*(r_bar+t_bar)*(to_bar^2-X_bar^2))...
        /(r_bar*t_bar*(r_bar*t_bar+r_bar^2-X_bar*s_bar)))/...
        (1+(X_bar*(r_bar^2+to_bar^2-X_bar^2 ))/(r_bar^2*to_bar));
    end
end

if nargout>4
    if Lambda_c4 == 0
        if X ~= 0 && X ~= inf
            Kd=(Kv/pi^2)*(((Kb-Z_bar)/(Y_bar^2+(Kb-Z_bar)^2))*...
                (1+X_bar/(norm([X_bar,Y_bar,(Kb-Z_bar)],2)))+...
                (X_bar/(X_bar^2+Y_bar^2))*...
                (((Kb-Z_bar)/(norm([X_bar,Y_bar,(Kb-Z_bar)],2)))+...
                (Kb+Z_bar/(norm([X_bar,Y_bar,(Kb+Z_bar)],2))))+...
                ((Kb+Z_bar)/(Y_bar^2+(Kb+Z_bar)^2))*...
                (1+X_bar/(norm([X_bar,Y_bar,(Kb+Z_bar)],2))));
        elseif X == inf
            Kd=4*Kv/(pi^2);
        else
            Kd=(Kv/pi^2)*(((Kb-Z_bar)/(Y_bar^2+(Kb-Z_bar)^2))*...
                (1)+((Kb+Z_bar)/(Y_bar^2+(Kb+Z_bar)^2))*(1));
        end
    else
        b_prime=Kb*bw;
        s_prime=(b_prime/2)*tan(deg2rad(Lambda_c4));
        if X ~= 0 && X ~= inf
            Kd=-((Kv*2*bw)/(2*pi^2))*...
                (((-b_prime/2)/((X-s_prime)^2+Y^2+(b_prime/2)^2-...
                (X-s_prime)*norm([(X-s_prime),Y,(b_prime/2)],2)))+...
                (-(b_prime/2)*X*(norm([X,Y],2)+...
                norm([(X-s_prime),Y,(b_prime/2)],2)))/...
                ((X^2+Y^2)*((X-s_prime)^2+Y^2+(b_prime/2)^2)+...
                (X^2-X*s_prime+Y^2)*norm([X,Y],2)*...
                norm([(X-s_prime),Y,(b_prime/2)],2)));
        elseif X == inf
            Kd=4*Kv/(pi^2);
        else
            Kd=-((Kv*2*bw)/(2*pi^2))*...
                (-b_prime)/((-s_prime)^2+Y^2+(b_prime/2)^2-...
                (-s_prime)*norm([(-s_prime),Y,(b_prime/2)],2));
        end
    end
end

end