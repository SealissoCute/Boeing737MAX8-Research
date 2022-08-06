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