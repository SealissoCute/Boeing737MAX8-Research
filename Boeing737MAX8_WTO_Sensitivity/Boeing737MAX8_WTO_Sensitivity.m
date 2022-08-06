%% Reverse Engineering of Boeing 737 MAX 8
% Clear previous data
format long
clc
clear

%% Output file directory (Notice: Need to be customized at different computer)
WorkspaceSavedDirectory = 'Boeing737MAX8_WTO_Sensitivity/Boeing737MAX8_WTO_Sensitivity.mat';
SelectedResultOutputDirectory = 'Boeing737MAX8_WTO_Sensitivity/Boeing737MAX8_WTO_Sensitivity_Result/Boeing737MAX8_WTO_Sensitivity_Result.txt';
RecordTimeDirectory = 'Boeing737MAX8_WTO_Sensitivity/Boeing737MAX8_WTO_Sensitivity_RunTimeRecord/Boeing737MAX8_WTO_Sensitivity_RunTimeRecord.txt';

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

% CruiseAltitude (unit: ft)
CruiseAltitudeMin = 23000;
CruiseAltitudeMax = 41000;
CruiseAltitudeInterval = 1000;
CruiseAltitudeMatrix = [CruiseAltitudeMin:CruiseAltitudeInterval:CruiseAltitudeMax]; % 

% Range (unit: nm)
RangeMin = 2000;
RangeMax = 3500;
RangeInterval = 100;
RangeMatrix = [RangeMin]; % :RangeInterval:RangeMax

% LoverD_Cruise
LoverD_CruiseMin = 13;
LoverD_CruiseMax = 14;
LoverD_CruiseInterval = 0.1;
LoverD_CruiseMatrix = [LoverD_CruiseMin]; % :LoverD_CruiseInterval:LoverD_CruiseMax

% LoverD_Loiter
LoverD_LoiterMin = 17;
LoverD_LoiterMax = 18;
LoverD_LoiterInterval = 0.1;
LoverD_LoiterMatrix = [LoverD_LoiterMin]; % :LoverD_LoiterInterval:LoverD_LoiterMax

% c_j_cruise
c_j_cruiseMin = 0.5;
c_j_cruiseMax = 0.6;
c_j_cruiseInterval = 0.01;
c_j_cruiseMatrix = [c_j_cruiseMin]; % :c_j_cruiseInterval:c_j_cruiseMax

% c_j_loiter
c_j_loiterMin = 0.5;
c_j_loiterMax = 0.6;
c_j_loiterInterval =0.01;
c_j_loiterMatrix = [c_j_loiterMin]; % :c_j_loiterInterval:c_j_loiterMax

%
Endurance = 0.5;            % Loiter, unit: hr 
AverageClimbRate = 2500;    % unit: fpm
CruiseSpeed_Mach = 0.79;    % unit: Mach
AlternateCruiseSpeed = 250; % unit: kts
AlternateRange = 100;       % unit: nm

% W_TO & W_OE from wikipedia (unit: lbs)
W_OE_wiki = 99360;  
W_TO_wiki = 182200;
W_E_wiki = W_OE_wiki - W_TO_wiki*0.005 - W_crew;

% Regression Line Constants A and B of Equation 
% Airplane type: Transport jet
A = 0.0833;
B = 1.0383;

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
% ResultMatrix sizing
Result_row = width(CruiseAltitudeMatrix)*width(RangeMatrix)*width(LoverD_CruiseMatrix)...
             *width(LoverD_LoiterMatrix)*width(c_j_cruiseMatrix)*width(c_j_loiterMatrix);
Result_column = 20;
InputParametersMatrix_column = 14;
InputParametersMatrix = zeros(Result_row,InputParametersMatrix_column);
InputParametersMatrixTemp = zeros(Result_row,InputParametersMatrix_column);
ResultMatrixApprox = zeros(Result_row,Result_column);
ResultMatrix = zeros(Result_row,Result_column);

% Create InputParametersMatrix1
n=1;
for CruiseAltitude = CruiseAltitudeMin:CruiseAltitudeInterval:CruiseAltitudeMax
    for Range = RangeMin%:RangeInterval:RangeMax
        for LoverD_Cruise = LoverD_CruiseMin%:LoverD_CruiseInterval:LoverD_CruiseMax
            for LoverD_Loiter = LoverD_LoiterMin%:LoverD_LoiterInterval:LoverD_LoiterMax
                for c_j_cruise = c_j_cruiseMin%:c_j_cruiseInterval:c_j_cruiseMax
                    for c_j_loiter = c_j_loiterMin%:c_j_loiterInterval:c_j_loiterMax
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

% section saved
save(WorkspaceSavedDirectory);
time = now;
date = datetime(time,'ConvertFrom','datenum');
string_RecordTime1=[' RecordTime: ',datestr(date)];
string_Workspace_saved1=[' Workspace is saved 1'];
disp('----------------------------------------------------')
disp(string_Workspace_saved1)
disp(string_RecordTime1)

%% Parallel computing CruiseSpeed/AverageClimbSpeed/ClimbTime/CruiseRange/W5_W4_ratio/W6_W5_ratio/M_ff
parfor row1 = 1:height(InputParametersMatrix)                          
    % Temporary matrix for parallel computing
    temp1 = zeros(1,InputParametersMatrix_column);
    
    % Read data form InputParametersMatrix
    CruiseAltitude = InputParametersMatrixTemp(row1,1);
    Range = InputParametersMatrixTemp(row1,2);
    LoverD_Cruise = InputParametersMatrixTemp(row1,3);
    LoverD_Loiter = InputParametersMatrixTemp(row1,4);
    c_j_cruise  = InputParametersMatrixTemp(row1,5);
    c_j_loiter = InputParametersMatrixTemp(row1,6);
    
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
    temp1(1) = CruiseAltitude;
    temp1(2) = Range;
    temp1(3) = LoverD_Cruise;
    temp1(4) = LoverD_Loiter;
    temp1(5) = c_j_cruise;
    temp1(6) = c_j_loiter;
    temp1(7) = CruiseSpeed;
    temp1(8) = AverageClimbSpeed;
    temp1(9) = ClimbTime;
    temp1(10) = ClimbRange;
    temp1(11) = CruiseRange;
    temp1(12) = W5_W4_ratio;
    temp1(12) = W6_W5_ratio;
    temp1(13) = M_ff;
    temp1(14) = C;

    % Output calculate result into InputParametersMatrix
    InputParametersMatrix(row1, :) = temp1;  
    
end

% section saved
save(WorkspaceSavedDirectory);
time = now;
date = datetime(time,'ConvertFrom','datenum');
string_RecordTime2=[' RecordTime: ',datestr(date),];
string_Workspace_saved2=[' Workspace is saved 2'];
disp('----------------------------------------------------')
disp(string_Workspace_saved2)
disp(string_RecordTime2)

%% Plot W_E_real/W_E_tent_min/W_E_tent_max
%
D = W_crew + W_PL;
%
C_min = min(InputParametersMatrix(:,14));
C_max = max(InputParametersMatrix(:,14));
%
x_W_TO_guess = 0:100:500000;
y_W_E_real = 10.^((log10(x_W_TO_guess)-A)/B);
y_W_E_tent_min = C_min.*x_W_TO_guess - D;
y_W_E_tent_max = C_max.*x_W_TO_guess - D;
hold on
plot(x_W_TO_guess,y_W_E_real)
plot(x_W_TO_guess,y_W_E_tent_min)
plot(x_W_TO_guess,y_W_E_tent_max)
xlabel('W_T_Oguess')
ylabel('W_E')
legend('W_Ereal','W_Etent_m_i_n','W_Etent_m_a_x');
hold off

%% 
% Find the lower and upper bound of W_To_guess
syms x
W_TO_guess_min = vpasolve( A + B*log10(C_max*x - D) - log10(x) == 0 );
W_TO_guess_max = vpasolve( A + B*log10(C_min*x - D) - log10(x) == 0 );
W_TO_guess_LowerBound = floor(W_TO_guess_min);
W_to_guess_UpperBound = ceil(W_TO_guess_max);

%% Numerical approximation of W_TO_guess
%
parfor row2 = 1:height(ResultMatrixApprox)
    % Temporary matrix for parallel computing
    temp2 = zeros(1,Result_column);

    % Read data
    CruiseAltitude = InputParametersMatrix(row2,1);
    Range = InputParametersMatrix(row2,2);
    LoverD_Cruise = InputParametersMatrix(row2,3);
    LoverD_Loiter = InputParametersMatrix(row2,4);
    c_j_cruise  = InputParametersMatrix(row2,5);
    c_j_loiter = InputParametersMatrix(row2,6);
    M_ff = InputParametersMatrix(row2,13);
    C = InputParametersMatrix(row2,14);


    %
    for W_TO_guess = W_TO_guess_LowerBound:W_TO_wiki
        W_E_real = 10^((log10(W_TO_guess)-A)/B);
        W_E_tent = C*W_TO_guess - D;
        error = abs(W_E_tent - W_E_real)/ W_E_real;
        if error < 0.005
            % Output data
            temp2(1) = CruiseAltitude;
            temp2(2) = Range;
            temp2(3) = LoverD_Cruise;
            temp2(4) = LoverD_Loiter;
            temp2(5) = c_j_cruise;
            temp2(6) = c_j_loiter;
            temp2(7) = M_ff;
            temp2(8) = W_TO_guess;
            temp2(9) = W_E_tent;
            temp2(10) = W_E_real;

            % Output calculate result into Result matrix
            ResultMatrixApprox(row2, :) = temp2;
            break
        end
    end
end

save(WorkspaceSavedDirectory);
time = now;
date = datetime(time,'ConvertFrom','datenum');
string_RecordTime3=[' RecordTime: ',datestr(date)];
string_Workspace_saved3=[' Workspace is saved 3'];
disp('----------------------------------------------------')
disp(string_Workspace_saved3)
disp(string_RecordTime3)

%% Numerical solution of W_TO_guess
x1 = sym('x1', [1,height(ResultMatrix)]);
parfor row3 = 1:height(ResultMatrix)
    % Temporary matrix for parallel computing
    temp3 = zeros(1,Result_column);

    if ResultMatrixApprox(row3,8) > 0 && ResultMatrixApprox(row3,8) < W_TO_wiki  
        % Read data
        CruiseAltitude = InputParametersMatrix(row3,1);
        Range = InputParametersMatrix(row3,2);
        LoverD_Cruise = InputParametersMatrix(row3,3);
        LoverD_Loiter = InputParametersMatrix(row3,4);
        c_j_cruise  = InputParametersMatrix(row3,5);
        c_j_loiter = InputParametersMatrix(row3,6);
        M_ff = InputParametersMatrix(row3,13);
        C = InputParametersMatrix(row3,14);

        % Using Matlab vpasolve function
        W_TO_guess = vpasolve( A + B*log10(C*x1(row3) - D) - log10(x1(row3)) == 0 );

        % Computing 
        W_E_real = 10^((log10(W_TO_guess)-A)/B);
        W_E_tent = C*W_TO_guess-D;
        W_E_error =  abs(W_E_real - W_E_wiki)/W_E_real;

        % Output data
        temp3(1) = CruiseAltitude;
        temp3(2) = Range;
        temp3(3) = LoverD_Cruise;
        temp3(4) = LoverD_Loiter;
        temp3(5) = c_j_cruise;
        temp3(6) = c_j_loiter;
        temp3(7) = M_ff;
        temp3(8) = W_TO_guess;
        temp3(9) = W_E_tent;
        temp3(10) = W_E_real;
        temp3(11) = W_E_error;

        % Output calculate result into Result matrix
        ResultMatrix(row3, :) = temp3;
    end  
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
    for row = 1:1:height(ResultMatrix)
%         if ResultMatrix(row,11) < 0.000005
            % Read data
            CruiseAltitude = InputParametersMatrix(row,1);
            Range = InputParametersMatrix(row,2);
            LoverD_Cruise = InputParametersMatrix(row,3);
            LoverD_Loiter = InputParametersMatrix(row,4);
            c_j_cruise  = InputParametersMatrix(row,5);
            c_j_loiter = InputParametersMatrix(row,6);
            CruiseSpeed = InputParametersMatrix(row,7);
            M_ff = ResultMatrix(row,7);
            C = InputParametersMatrix(row,14);
            W_TO_guess = ResultMatrix(row,8);
             
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
            ResultMatrix(row,12) = W_TO_over_W_PL;
            ResultMatrix(row,13) = W_TO_over_W_E ;
            ResultMatrix(row,14) = W_TO_over_Range;
            ResultMatrix(row,15) = W_TO_over_Endurance;
            ResultMatrix(row,16) = W_TO_over_CriuseSpeed;
            ResultMatrix(row,17) = W_TO_over_c_j_Range;
            ResultMatrix(row,18) = W_TO_over_LoverD_Range;
            ResultMatrix(row,19) = W_TO_over_c_j_Loiter;
            ResultMatrix(row,20) = W_TO_over_LoverD_Loiter;

            % W_TO_guess Iteration Result
            string_CruiseAltitude=[' 1.CruiseAltitude = ',num2str(ResultMatrix(row,1)),' ft'];
            string_Range=[' 2.Range = ',num2str(ResultMatrix(row,2)),' nm'];
            string_LoverD_Cruise=[' 3.L/D Cruise = ',num2str(ResultMatrix(row,3))];
            string_LoverD_Loiter=[' 4.L/D Loiter = ',num2str(ResultMatrix(row,4))];
            string_c_j_cruise=[' 5.c_j Cruise = ',num2str(ResultMatrix(row,5))];
            string_c_j_loiter=[' 6.c_j Loiter = ',num2str(ResultMatrix(row,6))];
            string_M_ff=[' 7.M_ff = ',num2str(ResultMatrix(row,7))];
            string_W_TO_guess=[' 8.W_TO_guess = ',num2str(ResultMatrix(row,8)),' lbs'];
            string_W_E_tent=[' 9.W_E_tent = ',num2str(ResultMatrix(row,9)),' lbs'];
            string_W_E_real=[' 10.W_E_real = ',num2str(ResultMatrix(row,10)),' lbs'];
            string_W_E_error=[' 11.W_E_error = ',num2str(ResultMatrix(row,11)*100),' %% (compare with W_E_wiki = W_OE_wiki - W_TO_wiki*0.005  - W_crew)'];

            % Sensitivity Result
            string_W_TO_over_W_PL=[' 12.W_TO_over_W_PL = ',num2str(ResultMatrix(row,12))];
            string_W_TO_over_W_E=[' 13.W_TO_over_W_E = ',num2str(ResultMatrix(row,13))];
            string_W_TO_over_Range=[' 14.W_TO_over_Range = ',num2str(ResultMatrix(row,14)),' lbs/nm'];
            string_W_TO_over_Endurance=[' 15.W_TO_over_Endurance = ',num2str(ResultMatrix(row,15)),' lbs/hr'];
            string_W_TO_over_CriuseSpeed=[' 16.W_TO_over_CriuseSpeed = ',num2str(ResultMatrix(row,16)),' lbs/kt'];
            string_W_TO_over_c_j_Range=[' 17.W_TO_over_c_j_Range = ',num2str(ResultMatrix(row,17)),' lbs/lbs/lbs/hr'];
            string_W_TO_over_LoverD_Range=[' 18.W_TO_over_LoverD_Range = ',num2str(ResultMatrix(row,18)),' lbs'];
            string_W_TO_over_LoverD_Range=[' 18.W_TO_over_LoverD_Range = ',num2str(ResultMatrix(row,18)),' lbs'];
            string_W_TO_over_c_j_Loiter=[' 19.W_TO_over_c_j_Loiter = ',num2str(ResultMatrix(row,19)),' lbs/lbs/lbs/hr'];
            string_W_TO_over_LoverD_Loiter=[' 20.W_TO_over_LoverD_Loiter = ',num2str(ResultMatrix(row,20)),' lbs'];

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

            % Print result at command line
%             disp('----------------------------------------------------')           
%             disp(' W_TO_guess Iteration Result')
%             disp(string_CruiseAltitude)
%             disp(string_Range)
%             disp(string_LoverD_Cruise)
%             disp(string_LoverD_Loiter)
%             disp(string_c_j_cruise)
%             disp(string_c_j_loiter)
%             disp(string_M_ff)
%             disp(string_W_TO_guess)
%             disp(string_W_E_tent)
%             disp(string_W_E_real)
%             disp(string_W_E_error)
%             disp(' Sensitivity Result')
%             disp(string_W_TO_over_W_PL)
%             disp(string_W_TO_over_W_E)
%             disp(string_W_TO_over_Range) 
%             disp(string_W_TO_over_Endurance)            
%             disp(string_W_TO_over_CriuseSpeed)            
%             disp(string_W_TO_over_c_j_Range)            
%             disp(string_W_TO_over_LoverD_Range)           
%             disp(string_W_TO_over_c_j_Loiter)            
%             disp(string_W_TO_over_LoverD_Loiter)

            % Numbers of result are printed
            m=m+1;
%         end
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

%% Section record time
%
fid = fopen(RecordTimeDirectory,'wt');
    fprintf(fid,'----------------------------------------------------');
    fprintf(fid,'\n');
    fprintf(fid,string_StartTime);
    fprintf(fid,'\n');
    fprintf(fid,'----------------------------------------------------');
    fprintf(fid,'\n');
    fprintf(fid,string_Workspace_saved1);
    fprintf(fid,'\n');   
    fprintf(fid,string_RecordTime1);
    fprintf(fid,'\n');
    fprintf(fid,'----------------------------------------------------');
    fprintf(fid,'\n');
    fprintf(fid,string_Workspace_saved2);
    fprintf(fid,'\n');
    fprintf(fid,string_RecordTime2);
    fprintf(fid,'\n');
    fprintf(fid,'----------------------------------------------------');
    fprintf(fid,'\n');
    fprintf(fid,string_Workspace_saved3);
    fprintf(fid,'\n');
    fprintf(fid,string_RecordTime3);
    fprintf(fid,'\n');
    fprintf(fid,'----------------------------------------------------');
    fprintf(fid,'\n');
    fprintf(fid,string_Workspace_saved4);
    fprintf(fid,'\n');
    fprintf(fid,string_RecordTime4);
    fprintf(fid,'\n');

    %
    save(WorkspaceSavedDirectory);
    time = now;
    date = datetime(time,'ConvertFrom','datenum');
    string_EndTime=[' EndTime: ' ,datestr(date)];
    disp('----------------------------------------------------')
    disp(string_EndTime)

    fprintf(fid,'----------------------------------------------------');
    fprintf(fid,'\n');
    fprintf(fid,string_EndTime);
    
fclose(fid);

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
