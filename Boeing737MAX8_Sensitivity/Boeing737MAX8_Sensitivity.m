%% Reverse Engineering of Boeing 737 MAX 8
% Clear previous data

format long
clc
clear

%% Setup output file directory (Notice: Need to be customized at different computer)
prompt = ' Which OS are you using? (Windows/MacOS) Ans:';
OS = input(prompt,'s');
if strcmp(OS,'MacOS')
    Boeing737MAX8_WTO_WorkspaceSavedDirectory = 'Boeing737MAX8_WTO/Boeing737MAX8_WTO.mat';
    Boeing737MAX8_Sensitivity_WorkspaceSavedDirectory = 'Boeing737MAX8_Sensitivity/Boeing737MAX8_Sensitivity.mat';
    Boeing737MAX8_W_TO_Senitivity_WorkspaceSavedDirectory = 'Boeing737MAX8_Sensitivity/Boeing737MAX8_W_TO_Senitivity.mat';
    Boeing737MAX8_W_TO_Senitivity_OutputDirectory = 'Boeing737MAX8_Sensitivity/Boeing737MAX8_W_TO_Senitivity.txt;'
elseif strcmp(OS,'Windows') 
    Boeing737MAX8_WTO_WorkspaceSavedDirectory = 'G:\飛設\Boeing737MAX8-Research\Boeing737MAX8_WTO\Boeing737MAX8_WTO.mat';
    Boeing737MAX8_Sensitivity_WorkspaceSavedDirectory = 'G:\飛設\Boeing737MAX8-Research\Boeing737MAX8_Sensitivity\Boeing737MAX8_Sensitivity.mat';
    Boeing737MAX8_W_TO_Senitivity_WorkspaceSavedDirectory = 'G:\飛設\Boeing737MAX8-Research\Boeing737MAX8_Sensitivity\Boeing737MAX8_W_TO_Senitivity.mat';
    Boeing737MAX8_W_TO_Senitivity_OutputDirectory = 'G:\飛設\Boeing737MAX8-Research\Boeing737MAX8_Sensitivity\Boeing737MAX8_W_TO_Senitivity.txt';
else
    error;
end
load(Boeing737MAX8_WTO_WorkspaceSavedDirectory);
%% Start time record
%
time = now;
date = datetime(time,'ConvertFrom','datenum');
disp('----------------------------------------------------')
string_StartTime=[' StartTime: ' ,datestr(date)];
disp(string_StartTime)

%% Sensitivity section
% W_TO_Senitivity sizing
W_TO_Senitivity = zeros();

%
n = 0;

% Open TransportJet_WTO_CheatingVersion_result.txt
fid = fopen(Boeing737MAX8_W_TO_Senitivity_OutputDirectory,'wt');

for row = 1:W_TO_solutions_row
    if  W_TO_solutions(row,13) < error_wiki/250 
        n = n+1;
        % Read data
        CruiseAltitude = W_TO_solutions(row,1);
        Range = W_TO_solutions(row,2);
        LoverD_Cruise = W_TO_solutions(row,3);
        LoverD_Loiter = W_TO_solutions(row,4);
        c_j_cruise  = W_TO_solutions(row,5);
        c_j_loiter = W_TO_solutions(row,6);
        CruiseSpeed = W_TO_solutions(row,7);
        M_ff = W_TO_solutions(row,8);
        C = W_TO_solutions(row,9);
        W_TO = W_TO_solutions(row,10);
        W_E_tent = W_TO_solutions(row,11);
        W_E_real = W_TO_solutions(row,12);
        W_E_error = W_TO_solutions(row,13);

        % Sensitivity calculate
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

        % Output result
        W_TO_Senitivity(n,1) = CruiseAltitude;
        W_TO_Senitivity(n,2) = Range;
        W_TO_Senitivity(n,3) = LoverD_Cruise;
        W_TO_Senitivity(n,4) = LoverD_Loiter;
        W_TO_Senitivity(n,5) = c_j_cruise;
        W_TO_Senitivity(n,6) = c_j_loiter;
        W_TO_Senitivity(n,7) = CruiseSpeed;
        W_TO_Senitivity(n,8) = M_ff;
        W_TO_Senitivity(n,9) = W_TO;
        W_TO_Senitivity(n,10) = W_E_tent;
        W_TO_Senitivity(n,11) = W_E_real;
        W_TO_Senitivity(n,12) = W_E_error;
        W_TO_Senitivity(n,13) = W_TO_over_W_PL;
        W_TO_Senitivity(n,14) = W_TO_over_W_E ;
        W_TO_Senitivity(n,15) = W_TO_over_Range;
        W_TO_Senitivity(n,16) = W_TO_over_Endurance;
        W_TO_Senitivity(n,17) = W_TO_over_CriuseSpeed;
        W_TO_Senitivity(n,18) = W_TO_over_c_j_Range;
        W_TO_Senitivity(n,19) = W_TO_over_LoverD_Range;
        W_TO_Senitivity(n,20) = W_TO_over_c_j_Loiter;
        W_TO_Senitivity(n,21) = W_TO_over_LoverD_Loiter;

        % W_TO Iteration Result
        string_CruiseAltitude=[' 1.CruiseAltitude = ',num2str(W_TO_Senitivity(n,1)),' ft'];
        string_Range=[' 2.Range = ',num2str(W_TO_Senitivity(n,2)),' nm'];
        string_LoverD_Cruise=[' 3.L/D Cruise = ',num2str(W_TO_Senitivity(n,3))];
        string_LoverD_Loiter=[' 4.L/D Loiter = ',num2str(W_TO_Senitivity(n,4))];
        string_c_j_cruise=[' 5.c_j Cruise = ',num2str(W_TO_Senitivity(n,5))];
        string_c_j_loiter=[' 6.c_j Loiter = ',num2str(W_TO_Senitivity(n,6))];
        string_Cruisespeed=[' 7.Cruisespeed = ',num2str(W_TO_Senitivity(n,7)),' kt'];
        string_M_ff=[' 8.M_ff = ',num2str(W_TO_Senitivity(n,8))];
        string_W_TO_guess=[' 9.W_TO = ',num2str(W_TO_Senitivity(n,9)),' lbs'];
        string_W_E_tent=[' 10.W_E_tent = ',num2str(W_TO_Senitivity(n,10)),' lbs'];
        string_W_E_real=[' 11.W_E_real = ',num2str(W_TO_Senitivity(n,11)),' lbs'];
        string_W_E_error=[' 12.W_E_error = ',num2str(W_TO_Senitivity(n,12)*100),' %% (compare with W_E_wiki = W_OE_wiki - W_TO_wiki*0.005 - W_crew)'];

        % Sensitivity Result
        string_W_TO_over_W_PL=[' 13.W_TO_over_W_PL = ',num2str(W_TO_Senitivity(n,13))];
        string_W_TO_over_W_E=[' 14.W_TO_over_W_E = ',num2str(W_TO_Senitivity(n,14))];
        string_W_TO_over_Range=[' 15.W_TO_over_Range = ',num2str(W_TO_Senitivity(n,15)),' lbs/nm'];
        string_W_TO_over_Endurance=[' 16.W_TO_over_Endurance = ',num2str(W_TO_Senitivity(n,16)),' lbs/hr'];
        string_W_TO_over_CriuseSpeed=[' 17.W_TO_over_CriuseSpeed = ',num2str(W_TO_Senitivity(n,17)),' lbs/kt'];
        string_W_TO_over_c_j_Range=[' 18.W_TO_over_c_j_Range = ',num2str(W_TO_Senitivity(n,18)),' lbs/lbs/lbs/hr'];
        string_W_TO_over_LoverD_Range=[' 19.W_TO_over_LoverD_Range = ',num2str(W_TO_Senitivity(n,19)),' lbs'];
        string_W_TO_over_c_j_Loiter=[' 20.W_TO_over_c_j_Loiter = ',num2str(W_TO_Senitivity(n,20)),' lbs/lbs/lbs/hr'];
        string_W_TO_over_LoverD_Loiter=[' 21.W_TO_over_LoverD_Loiter = ',num2str(W_TO_Senitivity(n,21)),' lbs'];

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
        fprintf(fid,string_Cruisespeed );
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
    end
end

disp('----------------------------------------------------')
string_solutions=[' There are ',num2str(n),' solutions printed in txt file.'];
disp(string_solutions)
fprintf(fid,'----------------------------------------------------');
fprintf(fid,'\n');
fprintf(fid,string_solutions );

% Close TransportJet_WTO_CheatingVersion_result.txt
fclose(fid);

% section saved
save(Boeing737MAX8_Sensitivity_WorkspaceSavedDirectory);
save(Boeing737MAX8_W_TO_Senitivity_WorkspaceSavedDirectory,'W_TO_Senitivity');
time = now;
date = datetime(time,'ConvertFrom','datenum');
string_RecordTime=[' RecordTime: ',datestr(date)];
string_Workspace_saved=[' Workspace is saved'];
disp('----------------------------------------------------')
disp(string_Workspace_saved)
disp(string_RecordTime)
