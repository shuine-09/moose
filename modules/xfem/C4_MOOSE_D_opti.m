%% C4_MOOSE_D_OPTI
% Find the best diffusion coefficients for oxygen in the Zr alpha and beta
% phases, for T=1000�C to 1500�C.
% The results of the model (oxide thicknes, alpha thickness, weight gain)
% are fit to the experimental data of Cathcart-Pawel, extended by UW
% experiment when available.

clear all;close all
%% Global parameters
temperature = '1400C';

%% Reference experimental data

%% Experiment UW-Madison

if strcmp(temperature, '1100C')
    time_exp_UW = [100 200 500]%1000];% 500 1300 1500];       % [s]
    d_exp_UW = [22 30 47] %64];
    alpha_exp_UW = [20 31 48] %63];
    wg_exp_UW = [3.99 5.83 9.33] %12.50];% 9.52 14.80 15.92];     % [mg/cm^2]
end
if strcmp(temperature, '1200C')
    time_exp_UW = [50 200]% 500 1000 1500];% 500 800 1300];               % [s]
    d_exp_UW = [24 45] %68 94 117];
    alpha_exp_UW = [28 59] %88 123 147];
    wg_exp_UW = [4.90 10.19] %15.45 21.34 25.91];% 14.85 15.81 22.52];       % [mg/cm^2]
end


%% Cathcart-Pawel experiment

if strcmp(temperature, '1000C')
    time_exp_CP = [1130.5 725.2 301.1 925.0 300.2 483.4 678.6 1104.1 915.5 475.7 493.4 299.4];      % [s]
    d_exp_CP = [40.5 32.6 22.5 37.4 20.9 26.5 29.9 38.6 34.5 25.9 27.1 21.3];       % [um]
    alpha_exp_CP = [23.7 23.1 15.2 24.6 14.4 18.0 20.8 26.1 24.6 17.6 19.5 13.8];   % [um]
    wg_exp_CP = [7.04 5.77 3.95 6.57 3.69 4.67 5.30 6.81 6.13 4.57 4.80 3.73];      % [mg/cm^2]
    time_exp_CP_MAXI = [367.4 866.5 99.4 640.5 234.4 512.4];    %[s]
    d_exp_CP_MAXI = [25.8 40.6 13.7 35.3 20.6 31.0];            % [um]
    alpha_exp_CP_MAXI = [17.6 26.0 10.2 22.6 14.1 19.3];        % [um]
    wg_exp_CP_MAXI = [4.547 7.104 2.436 6.173 3.632 5.410];     % [um]
    time_CP_max = 1000;
end
if strcmp(temperature, '1100C')
    time_exp_CP = [495.6 252.7 313.4 512.4 348.1 283.4 252.1 159.5 173.9 449.3 458.9];
    d_exp_CP = [45.7 33.5 36.6 47.1 37.2 40.8 33.3 26.0 29.4 45.7 44.6];
    alpha_exp_CP = [41.6 33.0 36.6 45.4 37.0 37.0 24.5 23.1 26.9 43.5 44.5];
    wg_exp_CP = [8.77 6.49 7.12 9.10 7.25 7.82 6.22 4.97 5.61 8.79 8.66];
    time_CP_max = 500;
end
if strcmp(temperature, '1200C')
    time_exp_CP = [236.4 160.2 126.6 280.2 234.6 56.1 171.8 66.1 111.2];
    d_exp_CP = [49.4 42.8 38.9 53.5 50.2 26.6 43.2 28.6 37.2];
    alpha_exp_CP = [53.5 44.7 41.0 54.9 51.6 28.2 44.0 30.9 39.7];
    wg_exp_CP = [10.03 8.59 7.81 10.78 10.09 5.33 8.66 5.75 7.46];
    time_CP_max = 300;
end
if strcmp(temperature, '1300C')
    time_exp_CP = [138.4 151.4 83.4 124.4 123.2 32.9 30.0 111.1 62.2 79.9 59.3 57.7 34.3];
    d_exp_CP = [56.2 59.7 44.9 51.7 53.1 28.9 27.9 50.1 38.2 41.2 36.4 36.5 28.0];
    alpha_exp_CP = [66.9 71.8 53.4 61.7 60.6 32.5 31.7 58.9 45.2 48.6 42.8 42.2 32.4];
    wg_exp_CP = [11.87 12.61 9.45 10.96 11.13 6.01 5.81 10.57 8.05 8.72 7.68 7.67 5.89];
    time_CP_max = 200;
end
if strcmp(temperature, '1400C')
    time_exp_CP = [72.7 62.2 29.7 28.6 45.7 42.6 14.3 13.5 11.2 10.4 56.4 54.2 44.8 43.4 27.5 29.8 74.1 78.9 29.3 26.5];
    %time_exp_CP = [72.7 62.2 29.7 28.6 45.7 42.6 56.4 54.2 44.8 43.4 27.5 29.8 74.1 78.9 29.3 26.5];
    d_exp_CP = [58.4 53.4 38.7 35.2 45.7 44.7 28.0 26.1 24.5 23.4 51.6 48.8 47.1 45.1 39.1 39.2 59.0 60.3 39.0 35.1];
    %d_exp_CP = [58.4 53.4 38.7 35.2 45.7 44.7 51.6 48.8 47.1 45.1 39.1 39.2 59.0 60.3 39.0 35.1];
    alpha_exp_CP = [78.3 70.9 51.0 49.2 62.9 60.4 36.8 35.7 31.8 30.7 67.5 65.7 62.2 60.4 51.1 51.4 77.3 77.6 50.4 47.2];
    %alpha_exp_CP = [78.3 70.9 51.0 49.2 62.9 60.4 67.5 65.7 62.2 60.4 51.1 51.4 77.3 77.6 50.4 47.2];
    wg_exp_CP = [12.84 11.74 8.45 7.84 10.12 9.85 6.08 5.74 5.31 5.09 11.29 10.79 10.29 9.92 8.47 8.53 12.92 13.18 8.46 7.73];
    %wg_exp_CP = [12.84 11.74 8.45 7.84 10.12 9.85 11.29 10.79 10.29 9.92 8.47 8.53 12.92 13.18 8.46 7.73];
    time_CP_max = 100;
end
if strcmp(temperature, '1500C')
    time_exp_CP = [31.8 28.9 33.0 31.2 47.0 42.6 22.8 21.9 7.6 8.2 53.2 9.2 8.9 49.0 50.0 42.4 37.8 13.9 13.2];
    %time_exp_CP = [31.8 28.9 33.0 31.2 47.0 42.6 22.8 21.9 53.2 49.0 50.0 42.4 37.8];
    d_exp_CP = [53.6 50.7 57.8 50.3 65.6 61.8 45.8 43.8 27.3 28.5 73.7 30.0 28.9 62.6 62.5 65.9 57.5 37.5 35.5];
    %d_exp_CP = [53.6 50.7 57.8 50.3 65.6 61.8 45.8 43.8 73.7 62.6 62.5 65.9 57.5];
    alpha_exp_CP = [72.1 68.1 75.7 74.2 88.2 83.8 63.6 62.8 38.7 38.5 97.6 40.7 40.8 93.9 91.4 85.6 80.3 48.3 48.6];
    %alpha_exp_CP = [72.1 68.1 75.7 74.2 88.2 83.8 63.6 62.8 97.6 93.9 91.4 85.6 80.3];
    wg_exp_CP = [11.98 11.34 12.76 11.53 14.65 13.84 10.29 9.93 6.13 6.34 16.29 6.68 6.50 14.42 14.34 14.51 12.97 8.25 7.93];
    %wg_exp_CP = [11.98 11.34 12.76 11.53 14.65 13.84 10.29 9.93 16.29 14.42 14.34 14.51 12.97];
    time_CP_max = 50;
end

%% Initialization

% Start with the values chosen by Leo 
Da_Matlab =[0.5 2.425 10.3 29.9063 75.25 170];        % Diffusion coefficients in the aplha phase for every T [�m�/s]
Db_Matlab =[30 81.8 250 687.1875 702 1200];         % Diffusion coefficients in the beta phase for every T [�m�/s]

if strcmp(temperature, '1000C')
    Da = Da_Matlab(1);
    Db = Db_Matlab(1);
end
if strcmp(temperature, '1100C')
    Da = Da_Matlab(2);
    Db = Db_Matlab(2);
end
if strcmp(temperature, '1200C')
    Da = Da_Matlab(3);
    Db = Db_Matlab(3);
end
if strcmp(temperature, '1300C')
    Da = Da_Matlab(4);
    Db = Db_Matlab(4);
end
if strcmp(temperature, '1400C')
    Da = Da_Matlab(5);
    Db = Db_Matlab(5);
end
if strcmp(temperature, '1500C')
    Da = Da_Matlab(6);
    Db = Db_Matlab(6);
end

% Define criteria for the optimization to stop
tolerance = 0.01;                             % Tolerance for each of the 3 NRMSE
Da_bound = [0.01*Da 100*Da];                  % Bounds of the Da search : 2 orders of magnitude below and above
Db_bound = [0.01*Db 100*Db];                  % Bounds of the Db search : 2 orders of magnitude below and above
max_iter = 100;                               % Maximum nb of iterations
residual_type = 5;                            % To choose which type of residual to consider

% Initial step 
step_Da = 7; %0.1*Da
step_Db = 70; %0.1*Db

% Create different matrices that will store everything needed

next_D = zeros(5,2);                          % Matrix with the 5 possible (Da,Db) values for next step
next_errors = zeros(5,5);                     % Matrix with the 5 types of residual for each of the next 5 values possible

D_history = zeros(max_iter+1,2);               % Matrix with all succesive values of (Da,Db)
D_history(1,:) = [Da Db];

error_history = zeros(max_iter+1,5);                % Matrix with all succesive values for various NRMSE

%% Initialization
% This file will be located and run in louis/projects/moose/modules/xfem

% Run the file with the specified input
command_in = strcat(['./xfem-opt -i test/tests/moving_interface/moving_oxide_C4_for_opti_2.i ' ...
       'UserObjects/velocity_a_b/diffusivity_beta='],num2str(Db),...
       ' UserObjects/velocity_a_b/diffusivity_alpha=',num2str(Da),...
       ' UserObjects/velocity_ox_a/diffusivity_alpha=',num2str(Da), ...
       ' Materials/diffusivity_beta/prop_values=',num2str(Db), ...
       ' Materials/diffusivity_alpha/prop_values=',num2str(Da));

system('make -j 4');  
system(command_in);

% Retrieve the output  csv file (postprocessors)
R1=1; %skip column names
C1=0;
outputFilename='test/tests/moving_interface/moving_oxide_C4_for_opti_2_out.csv';
data = csvread(outputFilename,R1,C1);
% Extract output of interest
time_MOOSE = data(1:length(data),1);    %% to change : put length(data,1) instead of 1481
alpha_MOOSE = data(1:length(data),2);
oxide_MOOSE = data(1:length(data),3);
%wg_J_MOOSE = data(1:length(data),9);
wg_C_MOOSE = data(1:length(data),7);

% Calls file that computes the errors (NRMSE)
C4_MOOSE_least_square;

% Put the initial outputs in the corresponding matrices
next_errors(1,1:3) = [wg_error oxide_error alpha_error];
sum_error = wg_error + oxide_error + alpha_error;
LC_error = 0.2*wg_error + 0.2*oxide_error + 0.6*alpha_error; 
next_errors(1,4:5) = [sum_error LC_error];

error_history(1,1:5) = next_errors(1,1:5);
residual = next_errors(1,residual_type);

progress = waitbar(0, sprintf('step %d/%d', 0, max_iter), 'Name', 'Optimizing ...');

%% Iteration
iter = 1;

while (iter <= max_iter) && (wg_error > tolerance || oxide_error > tolerance || alpha_error > tolerance) && ((Da > Da_bound(1) && Da < Da_bound(2)) || (Db > Db_bound(1) && Db < Db_bound(2))) && (step_Da > 1e-3*Da && step_Db > 1e-3*Db)  
    
    waitbar(iter / max_iter, progress, sprintf('step %d/%d', iter, max_iter));
    
    next_D(1,:) = [Da Db];
    next_D(2,:) = [min(Da+step_Da,Da_bound(2)) Db];
    next_D(3,:) = [max(Da-step_Da,Da_bound(1)) Db];
    next_D(4,:) = [Da min(Db+step_Db,Db_bound(2))];
    next_D(5,:) = [Da max(Db-step_Db,Db_bound(1))];
    
    for iter_D = 2:5
        Da = next_D(iter_D,1);
        Db = next_D(iter_D,2);
        
        command_in = strcat(['./xfem-opt -i test/tests/moving_interface/moving_oxide_C4_for_opti_2.i ' ...
            'UserObjects/velocity_a_b/diffusivity_beta='],num2str(Db),...
            ' UserObjects/velocity_a_b/diffusivity_alpha=',num2str(Da),...
            ' UserObjects/velocity_ox_a/diffusivity_alpha=',num2str(Da), ...
            ' Materials/diffusivity_beta/prop_values=',num2str(Db), ...
            ' Materials/diffusivity_alpha/prop_values=',num2str(Da));

        system(command_in);
        
        data = csvread(outputFilename,R1,C1);
        % Extract output of interest
        time_MOOSE = data(1:length(data),1);    %% to change : put length(data,1) instead of 1481
        alpha_MOOSE = data(1:length(data),2);
        oxide_MOOSE = data(1:length(data),3);
        %wg_J_MOOSE = data(1:length(data),9);
        wg_C_MOOSE = data(1:length(data),7);
        
        % Calls file that computes the errors (NRMSE)
        C4_MOOSE_least_square;
        
        % Put the initial outputs in the corresponding matrices
        next_errors(iter_D,1:3) = [wg_error oxide_error alpha_error]
        sum_error = wg_error + oxide_error + alpha_error
        LC_error = 0.2*wg_error + 0.2*oxide_error + 0.6*alpha_error
        next_errors(iter_D,4:5) = [sum_error LC_error];
    end
    
    [min_error, index_min_error] = min(next_errors(1:5,residual_type));
    
    if min_error < next_errors(1,residual_type)
        Da = next_D(index_min_error,1);
        Db = next_D(index_min_error,2);
        next_errors(1,:) = next_errors(index_min_error,:);
        D_history(iter+1,:) = [Da Db];
        error_history(iter+1,:) = next_errors(index_min_error,:);
    else
        step_Da = step_Da / 2
        step_Db = step_Db / 2
        Da = next_D(1,1);
        Db = next_D(1,2);
        D_history(iter+1,:) = [Da Db];
        error_history(iter+1,:) = next_errors(1,:);
    end
    
    residual = error_history(iter+1,residual_type);
    iter = iter + 1;
    
    %save(strcat('output_file/', file_name, '_HT_optimization_dicho'), 'D_mat', 'DaDb', 'Da_bound', 'Db_bound', 'step_Da', 'step_Db', '-regexp', 'diff_LS', 'diff_LS_mat', 'residual_idx', 'file_name');
    
end