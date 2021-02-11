%% C4_least_square
%%Determine a least square criteria to compare the CP (+UW) experimental
%%data
%%The criteria is the normalized root-mean-square deviation

%%Author : Louis Bailly-Salins
%%Email : baillysalins@wisc.edu

%%Last updated : 10/26/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Weight gain least square criteria

% C4 model vs CP exp
idx = zeros(1, length(wg_exp_CP));
wg_error = 0;

for k = 1:length(wg_exp_CP)
    t_curr = time_exp_CP(k);
    [val idx(k)] = min(abs(time_MOOSE - t_curr));
    wg_error = sqrt(wg_error^2 + ((wg_C_MOOSE(idx(k)) - wg_exp_CP(k)) / wg_exp_CP(k))^2 / length(wg_exp_CP));
end

% C4 model vs CP + UW exp
if exist('wg_exp_UW', 'var') == 1
    idx = zeros(1, length(wg_exp_CP)+length(wg_exp_UW));
    wg_error = 0;
    
    for k = 1:length(wg_exp_CP)
        t_curr = time_exp_CP(k);
        [val idx(k)] = min(abs(time_MOOSE - t_curr));
        wg_error = sqrt(wg_error^2 + ((wg_C_MOOSE(idx(k)) - wg_exp_CP(k)) / wg_exp_CP(k))^2 / (length(wg_exp_CP)+length(wg_exp_UW)));
    end
    for k = 1:length(wg_exp_UW)
        t_curr = time_exp_UW(k);
        [val idx(k)] = min(abs(time_MOOSE - t_curr));
        wg_error = sqrt(wg_error^2 + ((wg_C_MOOSE(idx(k)) - wg_exp_UW(k)) / wg_exp_UW(k))^2 / (length(wg_exp_CP)+length(wg_exp_UW)));
    end
end


%% Oxide thickness least square criteria

% C4 model vs CP exp
idx = zeros(1, length(d_exp_CP));
oxide_error = 0;

for k = 1:length(d_exp_CP)
    t_curr = time_exp_CP(k);
    [val idx(k)] = min(abs(time_MOOSE - t_curr));
    oxide_error = sqrt(oxide_error^2 + ((oxide_MOOSE(idx(k)) - d_exp_CP(k)) / d_exp_CP(k))^2 / length(d_exp_CP));
end

% C4 model vs CP + UW exp
if exist('d_exp_UW', 'var') == 1
    idx = zeros(1, length(d_exp_CP)+length(d_exp_UW));
    oxide_error = 0;
    
    for k = 1:length(d_exp_CP)
        t_curr = time_exp_CP(k);
        [val idx(k)] = min(abs(time_MOOSE - t_curr));
        oxide_error = sqrt(oxide_error^2 + ((oxide_MOOSE(idx(k)) - d_exp_CP(k)) / d_exp_CP(k))^2 / (length(d_exp_CP)+length(d_exp_UW)));
    end
    for k = 1:length(d_exp_UW)
        t_curr = time_exp_UW(k);
        [val idx(k)] = min(abs(time_MOOSE - t_curr));
        oxide_error = sqrt(oxide_error^2 + ((oxide_MOOSE(idx(k)) - d_exp_UW(k)) / d_exp_UW(k))^2 / (length(d_exp_CP)+length(d_exp_UW)));
    end
end


%% Alpha thickness least square criteria

% C4 model vs CP exp
idx = zeros(1, length(alpha_exp_CP));
alpha_error = 0;

for k = 1:length(alpha_exp_CP)
    t_curr = time_exp_CP(k);
    [val idx(k)] = min(abs(time_MOOSE - t_curr));
    alpha_error = sqrt(alpha_error^2 + ((alpha_MOOSE(idx(k)) - alpha_exp_CP(k)) / alpha_exp_CP(k))^2 / length(alpha_exp_CP));
end

% C4 model vs CP + UW exp
if exist('alpha_exp_UW', 'var') == 1
    idx = zeros(1, length(alpha_exp_CP)+length(alpha_exp_UW));
    alpha_error = 0;
    
    for k = 1:length(alpha_exp_CP)
        t_curr = time_exp_CP(k);
        [val idx(k)] = min(abs(time_MOOSE - t_curr));
        alpha_error = sqrt(alpha_error^2 + ((alpha_MOOSE(idx(k)) - alpha_exp_CP(k)) / alpha_exp_CP(k))^2 / (length(alpha_exp_CP)+length(alpha_exp_UW)));
    end
    for k = 1:length(alpha_exp_UW)
        t_curr = time_exp_UW(k);
        [val idx(k)] = min(abs(time_MOOSE - t_curr));
        alpha_error = sqrt(alpha_error^2 + ((alpha_MOOSE(idx(k)) - alpha_exp_UW(k)) / alpha_exp_UW(k))^2 / (length(alpha_exp_CP)+length(alpha_exp_UW)));
    end
end