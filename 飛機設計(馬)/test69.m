% The following data are obtained from the design requirements:
Vs = 8.5; % The stall speed according to the requirements in
 % certification specification EASA CS VLA, m/s
Vc = 26.0; % The cruising speed according to the design requirements,m/s
Vmax = 33.8; % Calculated maximum speed, m/s
Vto = 2.915; % Calculated take-off speed, m/s
Vto = 11.05;
Vr = Vto; % Take-off rotation speed, m/s
hC = 350; % Normal service altitude/ceiling above sea level, m
hac = 5000; % Absolute ceiling altitude, m
Clmax = 1.6; % Maximum lift coefficienfor the preliminary design phase
e = 0.8; % Oswald efficiency factor
AR = 12; % Wing aspect ratio for the preliminary design phase
K = 0.0331741; % Calculated induced drag coefficient
g = 9.81; % Gravitational acceleration, m/s^2
Cd0 = 0.0245; % Zero lift-drag coefficient
Cd0to = 0.0835; % Zero lift-drag coefficient at take-off
Clto = 0.85; % Aircraft lift coefficient at take-off
Cdto = 0.10747; % Aircraft drag coefficient at take-off
Cdg = 0.03947; % Coefficient
Clr = Clto; % Lift coefficient at take-off rotation
nu = 0.08; % Drag coefficient for the launch unit
Sto = 2; % Launch unit length
rhosl = 1.225; % Air density at sea level
rhoc = 1.184; % Air density at a cruising altitude of 350 m above sea level
rhoac = 0.736; % Air density at absolute ceiling altitude
mupto = 0.55; % Propeller efficiency coefficient at take-off
mupac = 0.8; % Propeller efficiency coefficient at cruising altitude
LDmax = 11.5; % Lift drag value for the preliminary design faze
ROCAC = 0; % Rate of climb at absolute ceiling, m/s
ROCSC = 0.5; % Rate of climb at service ceiling, m/s
ROCCrC = 1.5; % Rate of climb at cruise ceiling, m/s
ROCCoC = 5; % Rate of climb at combat ceiling, m/s
% Stall speed.
WS = 1/2*rhosl*Vs^2*Clmax;
x1 = WS;
x2 = WS;
y1 = 0;
y2 = 1.5;
plot([x1,x2],[y1,y2],'-g')
text(55,1.2,'Stall speed')
axis([0 80 -0.5 1.5])
xlabel('W/S, N/m^2')
ylabel('W/P, N/W')
grid on
hold on
% Maximum speed.
WSms = 0:2:80;
WPvmax = mupac./((0.5*rhosl*Vmax^3*Cd0./WSms)+(((2*K)./(rhoc*(rhoc/rhosl)*Vmax)).*WSms));
plot(WSms,WPvmax,'--r')
text(10,-0.05,'Maximum speed')
% Take-off run.
WPsto = (((1-exp(0.6*rhosl*g*Cdg*Sto)./WSms))./(nu-(nu+Cdg/Clr).*(exp(0.6*rhosl*g*Cdg*Sto)./WSms))).*(mupto/Vto);
disp(WPsto)
plot(WSms,WPsto,'b--o')
text(5,1.2,'Take-off run')
% Rate of Climb.
WProc = 1./(3.6363+(sqrt(1.0969.*WSms)*0.1826));
plot(WSms,WProc,'*-c')
text(5,0.3,'Rate of clime')
% Cruise ceiling.
WPslc =(rhoc/rhosl)./((ROCCrC/mupac)+sqrt((2/(rhoc*sqrt(3*Cd0/K)))*WSms)*(1.115/(LDmax*mupac)));
plot(WSms,WPslc,'*-y')
text(10,0.5,'Cruise ceiling')