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
W_TO_Approx_column = 13;
W_TO_Approx = zeros(W_TO_Apporx_row,W_TO_Approx_column);

parfor row = 1:W_TO_Apporx_row
    % Temporary matrix for parallel computing
    temp = zeros(1,W_TO_Approx_column);

    % Read data
    CruiseAltitude = InputParametersMatrix(row,1);
    Range = InputParametersMatrix(row,2);
    LoverD_Cruise = InputParametersMatrix(row,3);
    LoverD_Loiter = InputParametersMatrix(row,4);
    c_j_cruise  = InputParametersMatrix(row,5);
    c_j_loiter = InputParametersMatrix(row,6);
    CruiseSpeed = InputParametersMatrix(row,7);
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
            temp(7) = CruiseSpeed;
            temp(8) = M_ff;
            temp(9) = C;
            temp(10) = W_TO_guess;
            temp(11) = W_E_tent;
            temp(12) = W_E_real;
            temp(13) = W_E_error;

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
    if W_TO_Approx(row,10) > 0 && W_TO_Approx(row,10) < W_TO_wiki
        if W_TO_Approx(row,13) < 0.005 | W_TO_Approx(row,12) < W_E_wiki
            n = n+1;
            W_TO_ApproxSolutions(n,:) = W_TO_Approx(row,:);
        end
    end
end

disp('----------------------------------------------------')
string_solutions=[' There are ',num2str(n),' numerical approximation of W_TO_guess less than W_TO_wiki.'];
disp(string_solutions)

%% Numerical solution of W_TO_guess
% W_TO_solutions sizing
W_TO_solutions_row = height(W_TO_ApproxSolutions);
W_TO_sloutions_column = width(W_TO_ApproxSolutions) ;
W_TO_solutions = zeros(W_TO_solutions_row,W_TO_sloutions_column);

x = sym('x', [1,W_TO_solutions_row]);
parfor row = 1: W_TO_solutions_row
    % Temporary matrix for parallel computing
    temp = zeros(1,W_TO_sloutions_column);
   
            % Read data
            CruiseAltitude = W_TO_ApproxSolutions(row,1);
            Range = W_TO_ApproxSolutions(row,2);
            LoverD_Cruise = W_TO_ApproxSolutions(row,3);
            LoverD_Loiter = W_TO_ApproxSolutions(row,4);
            c_j_cruise  = W_TO_ApproxSolutions(row,5);
            c_j_loiter = W_TO_ApproxSolutions(row,6);
            CruiseSpeed = W_TO_ApproxSolutions(row,7);
            M_ff = W_TO_ApproxSolutions(row,8);
            C = W_TO_ApproxSolutions(row,9);
    
            % vpasolve
            W_TO_guess = vpasolve( A + B*log10(C*x(row) - D) - log10(x(row)) == 0 );
            
            % Computing
            W_E_real = 10^((log10(W_TO_guess)-A)/B);
            W_E_tent = C*W_TO_guess-D;
            W_E_error = abs(W_E_real - W_E_wiki)/W_E_real;
    
            % Output data
            temp(1) = CruiseAltitude;
            temp(2) = Range;
            temp(3) = LoverD_Cruise;
            temp(4) = LoverD_Loiter;
            temp(5) = c_j_cruise;
            temp(6) = c_j_loiter;
            temp(7) = CruiseSpeed;
            temp(8) = M_ff;
            temp(9) = C;
            temp(10) = W_TO_guess;
            temp(11) = W_E_tent;
            temp(12) = W_E_real;
            temp(13) = W_E_error;

            % Output calculate result into Result matrix
            W_TO_solutions(row, :) = temp;
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
