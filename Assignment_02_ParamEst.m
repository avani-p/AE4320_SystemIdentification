clc
close all
clear
%%%%%% Part - 4 Aerodynamic Model Identification
% IDENTIFICATION DATASETS
% Longitudinal Cx, Cz, Cm  -- de2311_1.mat
% Lateral Cy, Cl, Cn -- da3211_2.mat dr2311_1.mat
id_datasets = {'da3211_2.mat','dr3211_1.mat','de3211_1.mat'};

N_iddataset = length(id_datasets);
N_repeat    = 20;

% VALIDATION DATASET -- last one used for validation
% Longitudinal -- dedoublet_1.mat
% Lateral -- dadoublet_1.mat dr3211_2.mat
val_dataset = {'dadoublet_1.mat','dr3211_2.mat','dedoublet_1.mat'};
N_valdataset = length(val_dataset);


% Output Noise Standard Deviation Values (sigma)
sigma_gps_pos   = 2.5;              % [m]   GPS position
sigma_gps_vel   = 0.02;             % [m/s] GPS velocity
sigma_gps_att   = deg2rad(0.05);    % [rad] GPS attitudes
sigma_tas       = 0.1;              % [m/s] True airspeed, V
sigma_alpha     = deg2rad(0.1);     % [rad] Angle of Attack, α
sigma_beta      = deg2rad(0.1);     % [rad] Sideslip Angle, β

% Input Noise Standard Deviation Values (sigma) 
% for Accelerometer and Rate Gyro
sigma_acc = 0.02;               % [m/s^2]
sigma_rgyro = deg2rad(0.003);   % [rad/s]

% from Assignment
Ixx  =  11187.8;     % [kg.m^2]
Iyy  =  22854.8;     % [kg.m^2]
Izz  =  31974.8;     % [kg.m^2]
Ixz  =   1930.1;     % [kg.m^2]
mass =     4500;     % [kg]
b    =  13.3250;     % [m]   wing span
S    =  24.9900;     % [m^2] wing area
c    =   1.9910;     % chord
h    =     7500;     % height [m]

% Air Density at height h = 7500 m
[~, ~, ~, rho, ~] = atmosisa(h);  % [kg/m^3]

% g at h = 7500 m
g  = 9.80665*(6.3781e6)^2 / (6.3781e6 + h)^2;  % [m/s^2]

%% 4.0 Get data from multiple Datasets
dataset_CX_OLS_param    = [];
dataset_CY_OLS_param    = [];
dataset_CZ_OLS_param    = [];
dataset_Cl_OLS_param    = [];
dataset_Cm_OLS_param    = [];
dataset_Cn_OLS_param    = [];

dataset_CX_OLS_n_param    = [];
dataset_CY_OLS_n_param    = [];
dataset_CZ_OLS_n_param    = [];
dataset_Cl_OLS_n_param    = [];
dataset_Cm_OLS_n_param    = [];
dataset_Cn_OLS_n_param    = [];

dataset_CX_RLS_param    = [];
dataset_CY_RLS_param    = [];
dataset_CZ_RLS_param    = [];
dataset_Cl_RLS_param    = [];
dataset_Cm_RLS_param    = [];
dataset_Cn_RLS_param    = [];

coeff_name = {'C_X', 'C_Y', 'C_Z', 'C_l', 'C_m', 'C_n'};

CX_param_names = {'C_{X_0}'; 'C_{X_\alpha}'; 'C_{X_{\alpha^2}}'; 'C_{X_q}';...
                    'C_{X_{\delta_e}}'; 'C_{X_{T_c}}'};
CY_param_names = {'C_{Y_0}'; 'C_{Y_\beta}'; 'C_{Y_p}'; 'C_{Y_r}'; ...
                    'C_{Y_{\delta_a}}'; 'C_{Y_{\delta_r}}'};
CZ_param_names = {'C_{Z_0}'; 'C_{Z_\alpha}'; 'C_{Z_q}';'C_{Z_{\delta_e}}';...
                    'C_{Z_{T_c}}'};
Cl_param_names = {'C_{l_0}'; 'C_{l_\beta}'; 'C_{l_p}'; 'C_{l_r}';...
                    'C_{l_{\delta_a}}'; 'C_{l_{\delta_r}}'};
Cm_param_names = {'C_{m_0}'; 'C_{m_\alpha}'; 'C_{m_q}';'C_{m_{\delta_e}}';...
                    'C_{m_{T_c}}'};
Cn_param_names = {'C_{n_0}'; 'C_{n_\beta}'; 'C_{n_p}'; 'C_{n_r}';...
                    'C_{n_{\delta_a}}'; 'C_{n_{\delta_r}}'};


for number = 1:N_iddataset
    dataset_CX_OLS_res      = [];
    dataset_CY_OLS_res      = [];
    dataset_CZ_OLS_res      = [];
    dataset_Cl_OLS_res      = [];
    dataset_Cm_OLS_res      = [];
    dataset_Cn_OLS_res      = [];
    dataset_CX_RLS_res      = [];
    dataset_CY_RLS_res      = [];
    dataset_CZ_RLS_res      = [];
    dataset_Cl_RLS_res      = [];
    dataset_Cm_RLS_res      = [];
    dataset_Cn_RLS_res      = [];

    for seed = 1:N_repeat
        dataset_name = append('dataWorkspace/workspace_',string(2*(number)-1),string(seed),'_',id_datasets(number));
        loaded_data = load(dataset_name);
        start_idx = loaded_data.start_idx;
        t = loaded_data.t;
        %% 4.1  Compute the aerodynamic forces and moments from the trajectory and sensor parameters
        %       found using the EKF/IEKF. Clearly indicate how you performed the numerical differentiation to
        %       calculate the angular acceleration components required for the moment calculations.
       
        % denominator = (1/2)ρ*V^2*S
        den = (0.5*rho*(loaded_data.V_estm.^2)*S);
    
        %%%% Dimensionless Aerodynamic Force Body Components (Equation 13) %%%%%
        C_X = (mass*loaded_data.Ax_estm)./den;
        C_Y = (mass*loaded_data.Ay_estm)./den;
        C_Z = (mass*loaded_data.Az_estm)./den;
        
        %%%% Dimensionless Aerodynamic Moment Components  (Equation 14) %%%%
        C_l = ((loaded_data.p_dot_estm.*Ixx) + (loaded_data.q_estm.*loaded_data.r_estm.*(Izz - Iyy))...
                - ((loaded_data.p_estm.*loaded_data.q_estm) + loaded_data.r_dot_estm)*Ixz)./(den*b);
        C_m = ((loaded_data.q_dot_estm.*Iyy) + (loaded_data.r_estm.*loaded_data.p_estm.*(Ixx - Izz))...
                + (loaded_data.p_estm.^2 - loaded_data.r_estm.^2)*Ixz)./(den*c);
        C_n = ((loaded_data.r_dot_estm.*Izz) + (loaded_data.p_estm.*loaded_data.q_estm.*(Iyy - Ixx))...
                + ((loaded_data.q_estm.*loaded_data.r_estm) - loaded_data.p_dot_estm)*Ixz)./(den*b);

        %% 4.2  Formulate an OLS parameter estimator for the parameters of the aerodynamic model given in
        %       the assignment. Clearly indicate the structure of your linear regression model and estimator. (5
        %       points)
        %  4.3  Estimate the parameters of your aerodynamic model using the OLS estimator. Clearly indicate
        %       how (and why) you used the different data files in this process.
        useNewStructure = false;
        [OLS_CX_param, OLS_CY_param, OLS_CZ_param, OLS_Cl_param, OLS_Cm_param,...
            OLS_Cn_param, OLS_A_CX, OLS_A_CY, OLS_A_CZ, OLS_A_Cl, OLS_A_Cm, OLS_A_Cn] = ...
            ParamEst_OLS(c, b, C_X, C_Y, C_Z, C_l, C_m, C_n, loaded_data.p_estm, loaded_data.q_estm,...
            loaded_data.r_estm, loaded_data.da, loaded_data.de, loaded_data.dr, loaded_data.Tc1,...
            loaded_data.V_estm, loaded_data.alpha_estm, loaded_data.beta_estm, useNewStructure);
        
        % Calculating Aerodynamic Force and Moments from Parameter Estimation
        CX_OLS = OLS_A_CX * OLS_CX_param;
        CY_OLS = OLS_A_CY * OLS_CY_param;
        CZ_OLS = OLS_A_CZ * OLS_CZ_param;
        Cl_OLS = OLS_A_Cl * OLS_Cl_param;
        Cm_OLS = OLS_A_Cm * OLS_Cm_param;
        Cn_OLS = OLS_A_Cn * OLS_Cn_param;
        
        % Calculing OLS Residual
        OLS_CX_res = C_X - CX_OLS;
        OLS_CY_res = C_Y - CY_OLS;
        OLS_CZ_res = C_Z - CZ_OLS;
        OLS_Cl_res = C_l - Cl_OLS;
        OLS_Cm_res = C_m - Cm_OLS;
        OLS_Cn_res = C_n - Cn_OLS;


        % Second Model Structure (Part 4.6)
        useNewStructure = true;
        [OLS_CX_param_new, OLS_CY_param_new, OLS_CZ_param_new, ...
            OLS_Cl_param_new, OLS_Cm_param_new,OLS_Cn_param_new, ...
            OLS_A_CX_new, OLS_A_CY_new, OLS_A_CZ_new, ...
            OLS_A_Cl_new, OLS_A_Cm_new, OLS_A_Cn_new] = ...
        ParamEst_OLS(c, b, C_X, C_Y, C_Z, C_l, C_m, C_n, loaded_data.p_estm, loaded_data.q_estm,...
            loaded_data.r_estm, loaded_data.da, loaded_data.de, loaded_data.dr, loaded_data.Tc1,...
            loaded_data.V_estm, loaded_data.alpha_estm, loaded_data.beta_estm, useNewStructure);
        
        % Calculating Aerodynamic Force and Moments (new)
        CX_OLS_new = OLS_A_CX_new * OLS_CX_param_new;
        CY_OLS_new = OLS_A_CY_new * OLS_CY_param_new;
        CZ_OLS_new = OLS_A_CZ_new * OLS_CZ_param_new;
        Cl_OLS_new = OLS_A_Cl_new * OLS_Cl_param_new;
        Cm_OLS_new = OLS_A_Cm_new * OLS_Cm_param_new;
        Cn_OLS_new = OLS_A_Cn_new * OLS_Cn_param_new;
        
        % Calculing OLS Residual (new)
        OLS_CX_new_res = C_X - CX_OLS_new;
        OLS_CY_new_res = C_Y - CY_OLS_new;
        OLS_CZ_new_res = C_Z - CZ_OLS_new;
        OLS_Cl_new_res = C_l - Cl_OLS_new;
        OLS_Cm_new_res = C_m - Cm_OLS_new;
        OLS_Cn_new_res = C_n - Cn_OLS_new;

%% 4.4  Formulate an alternative estimator (WLS, RLS, MLE) for the parameters of the aerodynamic
%       model given in the assignment, and show how the resulting parameter estimates differ (or not)
%       from those estimated using OLS.
        
        % Initial Estimate of parameter covariance matrix
        P0_CX = eye(length(OLS_CX_param));
        P0_CY = eye(length(OLS_CY_param));
        P0_CZ = eye(length(OLS_CZ_param));
        P0_Cl = eye(length(OLS_Cl_param));
        P0_Cm = eye(length(OLS_Cm_param));
        P0_Cn = eye(length(OLS_Cn_param));
        
        % Defining Forgetting Factor
        lambda = 1;
       
        RLS_CX_param = ParamEst_RLS(OLS_A_CX, C_X, lambda,OLS_CX_param, P0_CX);
        RLS_CY_param = ParamEst_RLS(OLS_A_CY, C_Y, lambda,OLS_CY_param, P0_CY);
        RLS_CZ_param = ParamEst_RLS(OLS_A_CZ, C_Z, lambda,OLS_CZ_param, P0_CZ);
        RLS_Cl_param = ParamEst_RLS(OLS_A_Cl, C_l, lambda,OLS_Cl_param, P0_Cl);
        RLS_Cm_param = ParamEst_RLS(OLS_A_Cm, C_m, lambda,OLS_Cm_param, P0_Cm);
        RLS_Cn_param = ParamEst_RLS(OLS_A_Cn, C_n, lambda,OLS_Cn_param, P0_Cn);
        

        % Calculating Aerodynamic Force and Moments from Parameter Estimation
        CX_RLS = OLS_A_CX * RLS_CX_param;
        CY_RLS = OLS_A_CY * RLS_CY_param;
        CZ_RLS = OLS_A_CZ * RLS_CZ_param;
        Cl_RLS = OLS_A_Cl * RLS_Cl_param;
        Cm_RLS = OLS_A_Cm * RLS_Cm_param;
        Cn_RLS = OLS_A_Cn * RLS_Cn_param;
        
        % Calculing RLS Residual
        RLS_CX_res = C_X - CX_RLS;
        RLS_CY_res = C_Y - CY_RLS;
        RLS_CZ_res = C_Z - CZ_RLS;
        RLS_Cl_res = C_l - Cl_RLS;
        RLS_Cm_res = C_m - Cm_RLS;
        RLS_Cn_res = C_n - Cn_RLS;

        %%% SAVE DATA TO BE USED IN ANALYSIS
        dataset_CX_OLS_res     = [dataset_CX_OLS_res     OLS_CX_res];
        dataset_CY_OLS_res     = [dataset_CY_OLS_res     OLS_CY_res];
        dataset_CZ_OLS_res     = [dataset_CZ_OLS_res     OLS_CZ_res];
        dataset_Cl_OLS_res     = [dataset_Cl_OLS_res     OLS_Cl_res];
        dataset_Cm_OLS_res     = [dataset_Cm_OLS_res     OLS_Cm_res];
        dataset_Cn_OLS_res     = [dataset_Cn_OLS_res     OLS_Cn_res];
        
        dataset_CX_OLS_param     = [dataset_CX_OLS_param     OLS_CX_param];
        dataset_CY_OLS_param     = [dataset_CY_OLS_param     OLS_CY_param];
        dataset_CZ_OLS_param     = [dataset_CZ_OLS_param     OLS_CZ_param];
        dataset_Cl_OLS_param     = [dataset_Cl_OLS_param     OLS_Cl_param];
        dataset_Cm_OLS_param     = [dataset_Cm_OLS_param     OLS_Cm_param];
        dataset_Cn_OLS_param     = [dataset_Cn_OLS_param     OLS_Cn_param];

        dataset_CX_RLS_res     = [dataset_CX_RLS_res     RLS_CX_res];
        dataset_CY_RLS_res     = [dataset_CY_RLS_res     RLS_CY_res];
        dataset_CZ_RLS_res     = [dataset_CZ_RLS_res     RLS_CZ_res];
        dataset_Cl_RLS_res     = [dataset_Cl_RLS_res     RLS_Cl_res];
        dataset_Cm_RLS_res     = [dataset_Cm_RLS_res     RLS_Cm_res];
        dataset_Cn_RLS_res     = [dataset_Cn_RLS_res     RLS_Cn_res];

        dataset_CX_RLS_param     = [dataset_CX_RLS_param     RLS_CX_param];
        dataset_CY_RLS_param     = [dataset_CY_RLS_param     RLS_CY_param];
        dataset_CZ_RLS_param     = [dataset_CZ_RLS_param     RLS_CZ_param];
        dataset_Cl_RLS_param     = [dataset_Cl_RLS_param     RLS_Cl_param];
        dataset_Cm_RLS_param     = [dataset_Cm_RLS_param     RLS_Cm_param];
        dataset_Cn_RLS_param     = [dataset_Cn_RLS_param     RLS_Cn_param];
        

        dataset_CX_OLS_n_param     = [dataset_CX_OLS_n_param     OLS_CX_param_new];
        dataset_CY_OLS_n_param     = [dataset_CY_OLS_n_param     OLS_CY_param_new];
        dataset_CZ_OLS_n_param     = [dataset_CZ_OLS_n_param     OLS_CZ_param_new];
        dataset_Cl_OLS_n_param     = [dataset_Cl_OLS_n_param     OLS_Cl_param_new];
        dataset_Cm_OLS_n_param     = [dataset_Cm_OLS_n_param     OLS_Cm_param_new];
        dataset_Cn_OLS_n_param     = [dataset_Cn_OLS_n_param     OLS_Cn_param_new];


        % Adding data to workspace struct
        loaded_data.C_X = C_X; loaded_data.C_Y = C_Y; loaded_data.C_Z = C_Z;
        loaded_data.C_l = C_l; loaded_data.C_m = C_m; loaded_data.C_n = C_n;

        loaded_data.RLS_CX_param = RLS_CX_param; loaded_data.RLS_CY_param = RLS_CY_param;
        loaded_data.RLS_CZ_param = RLS_CZ_param; loaded_data.RLS_Cl_param = RLS_Cl_param;
        loaded_data.RLS_Cm_param = RLS_Cm_param; loaded_data.RLS_Cn_param = RLS_Cn_param;

        loaded_data.OLS_CX_param = OLS_CX_param; loaded_data.OLS_CY_param = OLS_CY_param;
        loaded_data.OLS_CZ_param = OLS_CZ_param; loaded_data.OLS_Cl_param = OLS_Cl_param;
        loaded_data.OLS_Cm_param = OLS_Cm_param; loaded_data.OLS_Cn_param = OLS_Cn_param;

        newname_data = append(string(id_datasets(number)),'_0',string(seed));
        newname_data = erase(newname_data, '.mat');
        assignin('base', newname_data, loaded_data);
    end
    
    OLS_residuals = [ dataset_CX_OLS_res(:,end) dataset_CY_OLS_res(:,end)...
                      dataset_CZ_OLS_res(:,end) dataset_Cl_OLS_res(:,end)...
                      dataset_Cm_OLS_res(:,end) dataset_Cn_OLS_res(:,end)];

    RLS_residuals = [ dataset_CX_RLS_res(:,end) dataset_CY_RLS_res(:,end)...
                      dataset_CZ_RLS_res(:,end) dataset_Cl_RLS_res(:,end)...
                      dataset_Cm_RLS_res(:,end) dataset_Cn_RLS_res(:,end)];

    PlotResidual(t,OLS_residuals,'OLS Residual',RLS_residuals,'RLS Residual',coeff_name,id_datasets(number))
    
end

%% 4.5  Prove the statistical quality of your parameter
%       estimates. Which parameter(s) are you most/less sure about? Clearly indicate why!

% STATISTICAL ANALYSIS

% Calculaete Variance per Dataset to Conclude on one value with least
% variance
for number = 1:N_iddataset
    dataset_var_CX(:,number) =  var(dataset_CX_OLS_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    dataset_var_CY(:,number) =  var(dataset_CY_OLS_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    dataset_var_CZ(:,number) =  var(dataset_CZ_OLS_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    
    dataset_var_Cl(:,number) =  var(dataset_Cl_OLS_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    dataset_var_Cm(:,number) =  var(dataset_Cm_OLS_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    dataset_var_Cn(:,number) =  var(dataset_Cn_OLS_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);


    dataset_var_CX_n(:,number) =  var(dataset_CX_OLS_n_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    dataset_var_CY_n(:,number) =  var(dataset_CY_OLS_n_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    dataset_var_CZ_n(:,number) =  var(dataset_CZ_OLS_n_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    
    dataset_var_Cl_n(:,number) =  var(dataset_Cl_OLS_n_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    dataset_var_Cm_n(:,number) =  var(dataset_Cm_OLS_n_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
    dataset_var_Cn_n(:,number) =  var(dataset_Cn_OLS_n_param(:,N_repeat*(number-1)+1:N_repeat*(number)),0,2);
end

PlotVariances(dataset_var_CX, CX_param_names, dataset_var_CY, CY_param_names,dataset_var_CZ, CZ_param_names,...
    dataset_var_Cl, Cl_param_names, dataset_var_Cm, Cm_param_names, dataset_var_Cn, Cn_param_names);


% minimum variance and their dataset index
[min_var_dataset_CX, min_var_CX_idx] = min(dataset_var_CX,[],2);
[min_var_dataset_CY, min_var_CY_idx] = min(dataset_var_CY,[],2);
[min_var_dataset_CZ, min_var_CZ_idx] = min(dataset_var_CZ,[],2);
[min_var_dataset_Cl, min_var_Cl_idx] = min(dataset_var_Cl,[],2);
[min_var_dataset_Cm, min_var_Cm_idx] = min(dataset_var_Cm,[],2);
[min_var_dataset_Cn, min_var_Cn_idx] = min(dataset_var_Cn,[],2);

% minimum variances and index for @nd model structure (Part 4.6)
[min_var_dataset_CX_n, min_var_CX_idx_n] = min(dataset_var_CX_n,[],2);
[min_var_dataset_CY_n, min_var_CY_idx_n] = min(dataset_var_CY_n,[],2);
[min_var_dataset_CZ_n, min_var_CZ_idx_n] = min(dataset_var_CZ_n,[],2);
[min_var_dataset_Cl_n, min_var_Cl_idx_n] = min(dataset_var_Cl_n,[],2);
[min_var_dataset_Cm_n, min_var_Cm_idx_n] = min(dataset_var_Cm_n,[],2);
[min_var_dataset_Cn_n, min_var_Cn_idx_n] = min(dataset_var_Cn_n,[],2);


% Extracting Better Estimates of Median and lowest variance 
for k = 1:9
    % 2ND STRUCTURE
        
    if(k<6)
        Median_CZ_param(k,1)  = median(dataset_CZ_OLS_param(k,N_repeat*(min_var_CZ_idx(k)-1)+1:N_repeat*(min_var_CZ_idx(k))));
        Median_Cm_param(k,1)  = median(dataset_Cm_OLS_param(k,N_repeat*(min_var_Cm_idx(k)-1)+1:N_repeat*(min_var_Cm_idx(k))));
    end
    if(k<7)
        Median_CX_param(k,1)  = median(dataset_CX_OLS_param(k,N_repeat*(min_var_CX_idx(k)-1)+1:N_repeat*(min_var_CX_idx(k))));
        Median_CY_param(k,1)  = median(dataset_CY_OLS_param(k,N_repeat*(min_var_CY_idx(k)-1)+1:N_repeat*(min_var_CY_idx(k))));
        Median_Cl_param(k,1)  = median(dataset_Cl_OLS_param(k,N_repeat*(min_var_Cl_idx(k)-1)+1:N_repeat*(min_var_Cl_idx(k))));
        Median_Cn_param(k,1)  = median(dataset_Cn_OLS_param(k,N_repeat*(min_var_Cn_idx(k)-1)+1:N_repeat*(min_var_Cn_idx(k))));

    end
    if(k<8)
        Median_CX_param_n(k,1)  = median(dataset_CX_OLS_n_param(k,N_repeat*(min_var_CX_idx_n(k)-1)+1:N_repeat*(min_var_CX_idx_n(k))));
        Median_CZ_param_n(k,1)  = median(dataset_CZ_OLS_n_param(k,N_repeat*(min_var_CZ_idx_n(k)-1)+1:N_repeat*(min_var_CZ_idx_n(k))));
        Median_Cm_param_n(k,1)  = median(dataset_Cm_OLS_n_param(k,N_repeat*(min_var_Cm_idx_n(k)-1)+1:N_repeat*(min_var_Cm_idx_n(k))));
    end

    Median_CY_param_n(k,1)  = median(dataset_CY_OLS_n_param(k,N_repeat*(min_var_CY_idx_n(k)-1)+1:N_repeat*(min_var_CY_idx_n(k))));
    Median_Cl_param_n(k,1)  = median(dataset_Cl_OLS_n_param(k,N_repeat*(min_var_Cl_idx_n(k)-1)+1:N_repeat*(min_var_Cl_idx_n(k))));
    Median_Cn_param_n(k,1)  = median(dataset_Cn_OLS_n_param(k,N_repeat*(min_var_Cn_idx_n(k)-1)+1:N_repeat*(min_var_Cn_idx_n(k))));

end

PlotVariances(min_var_dataset_CX, CX_param_names, min_var_dataset_CY, CY_param_names,min_var_dataset_CZ, CZ_param_names,...
    min_var_dataset_Cl, Cl_param_names, min_var_dataset_Cm, Cm_param_names, min_var_dataset_Cn, Cn_param_names,false);

% PARAMETER COMPARISON

showPlots = true;
deleteOutliers = false;
% For Longitudinal Forces and Moments (only last 10 data set -- of de)

PlotParamComparison(dataset_CX_OLS_param(:,end-N_repeat+1:end),dataset_CX_RLS_param(:,end-N_repeat+1:end),...
                        CX_param_names,deleteOutliers,showPlots);
PlotParamComparison(dataset_CZ_OLS_param(:,end-N_repeat+1:end),dataset_CZ_RLS_param(:,end-N_repeat+1:end),...
                        CZ_param_names,deleteOutliers,showPlots);
PlotParamComparison(dataset_Cm_OLS_param(:,end-N_repeat+1:end),dataset_Cm_RLS_param(:,end-N_repeat+1:end),...
                        Cm_param_names,deleteOutliers,showPlots);


% For Lateral Forces and Moments (first 20 datasets -- of da and dr)
PlotParamComparison(dataset_CY_OLS_param(:,1:end-N_repeat),dataset_CY_RLS_param(:,1:end-N_repeat),...
                       CY_param_names,deleteOutliers,showPlots);
PlotParamComparison(dataset_Cl_OLS_param(:,1:end-N_repeat),dataset_Cl_RLS_param(:,1:end-N_repeat),...
                        Cl_param_names,deleteOutliers,showPlots);
PlotParamComparison(dataset_Cn_OLS_param(:,1:end-N_repeat),dataset_Cn_RLS_param(:,1:end-N_repeat),...
                        Cn_param_names,deleteOutliers,showPlots);


delete('workspace_ParameterEstimation.mat')
save('workspace_ParameterEstimation.mat');
cont = input('Continue to Validation? (y/n) [n]   ',"s");

%% 4.5  Validate your aerodynamic model by using the estimated model output and measured data.
%       Provide a model residual analysis. 
if(cont == 'y')
close all; clc;
% VALIDATION DATASET -- last one used for validation
% Longitudinal -- dedoublet_1.mat
% Lateral -- dadoublet_1.mat dr3211_2.mat

for iter=1:N_valdataset
    load(append('dataWorkspace/workspace_',string(2*(iter)),string(randi(N_repeat)),'_',string(val_dataset(iter))));
    
    % Claculation of True Values
    pdot = gradient(p, dt);
    qdot = gradient(q, dt);
    rdot = gradient(r, dt);
    den = (1/2 .* rho .* S .* (V_airdata.^2));
    
    CX_true = mass*Ax ./ den;
    CY_true = mass*Ay ./ den;
    CZ_true = mass*Az ./ den;
    
    Cl_true = (pdot.*Ixx + q.*r.*(Izz - Iyy)...
            - (p.*q + rdot).*Ixz) ./ (den.*b);
    Cm_true = (qdot.*Iyy + r.*p.*(Ixx - Izz)...
            + ((p.^2) - (r.^2)).*Ixz) ./ (den.*c);
    Cn_true = (rdot.*Izz + p.*q.*(Iyy - Ixx)...
            + (q.*r - pdot).*Ixz) ./ (den.*b);
    
    
    % Estimated Kalman Filter Values
    
    % denominator = (1/2)ρ*V^2*S
    den = (0.5*rho*(V_estm.^2)*S);
    
    %%%% Dimensionless Aerodynamic Force Body Components (Equation 13) %%%%%
    CX_est = (mass*Ax_estm)./den;
    CY_est = (mass*Ay_estm)./den;
    CZ_est = (mass*Az_estm)./den;
    
    %%%% Dimensionless Aerodynamic Moment Components  (Equation 14) %%%%
    Cl_est = ((p_dot_estm.*Ixx) + (q_estm.*r_estm.*(Izz - Iyy))...
            - ((p_estm.*q_estm) + r_dot_estm)*Ixz)./(den*b);
    Cm_est = ((q_dot_estm.*Iyy) + (r_estm.*p_estm.*(Ixx - Izz))...
            + (p_estm.^2 - r_estm.^2)*Ixz)./(den*c);
    Cn_est = ((r_dot_estm.*Izz) + (p_estm.*q_estm.*(Iyy - Ixx))...
            + ((q_estm.*r_estm) - p_dot_estm)*Ixz)./(den*b);
    
    % Mathematical OLS Values - Parameter Estimation Values
    % Getting A matrix
    val_A_CX = [ones(N, 1), alpha_estm, alpha_estm.^2, q_estm*c./V_estm, de, Tc1];
    val_A_CY = [ones(N, 1),  beta_estm, p_estm*b./(2*V_estm), r_estm*b./(2*V_estm), da, dr];
    val_A_CZ = [ones(N, 1), alpha_estm, q_estm*c./V_estm, de, Tc1];
    
    val_A_Cl = [ones(N, 1),  beta_estm, p_estm*b./(2*V_estm), r_estm*b./(2*V_estm), da, dr];
    val_A_Cm = [ones(N, 1), alpha_estm, q_estm*c./V_estm, de, Tc1];
    val_A_Cn = [ones(N, 1),  beta_estm, p_estm*b./(2*V_estm), r_estm*b./(2*V_estm), da, dr];
    
    % Obtaning the aerodynamic forces and moments from validation dataset
    CX_OLS = val_A_CX * Median_CX_param;
    CY_OLS = val_A_CY * Median_CY_param;
    CZ_OLS = val_A_CZ * Median_CZ_param;
    Cl_OLS = val_A_Cl * Median_Cl_param;
    Cm_OLS = val_A_Cm * Median_Cm_param;
    Cn_OLS = val_A_Cn * Median_Cn_param;
    
    data_C_true = [CX_true  CY_true CZ_true Cl_true Cm_true Cn_true];
    data_C_KF   = [CX_est   CY_est  CZ_est  Cl_est  Cm_est  Cn_est];
    data_C_OLS  = [CX_OLS   CY_OLS  CZ_OLS  Cl_OLS  Cm_OLS  Cn_OLS];

    PlotAerodynamicModel(t,data_C_true, data_C_KF, data_C_OLS, coeff_name, val_dataset(iter),false);

    data_C_KF_res = data_C_true - data_C_KF;
    data_C_OLS_res = data_C_true - data_C_OLS;
    PlotResidual(t,data_C_OLS_res,'OLS Residual',data_C_KF_res,'KF Residual',coeff_name,val_dataset(iter))

    %% 4.6  Provide an indication which model terms are most dominant, and which have the least
    %       influence. Introduce an alternative model structure and compare its validation with the
    %       validation given in 4.3.

    % A Matrix of new model structure
    val_A_CX_n = [ones(N,1)   alpha_estm  alpha_estm.^2   q_estm*c./V_estm        (q_estm*c./V_estm).^2       de                    Tc1];
    val_A_CY_n = [ones(N,1)   beta_estm   beta_estm.^2    p_estm*b./(2*V_estm)    (p_estm*b./(2*V_estm)).^2   r_estm*b./(2*V_estm)  (r_estm*b./(2*V_estm)).^2     da  dr];
    val_A_CZ_n = [ones(N,1)   alpha_estm  alpha_estm.^2   q_estm*c./V_estm        (q_estm*c./V_estm).^2       de                    Tc1];
    
    val_A_Cl_n = [ones(N,1)   beta_estm   beta_estm.^2    p_estm*b./(2*V_estm)    (p_estm*b./(2*V_estm)).^2   r_estm*b./(2*V_estm)  (r_estm*b./(2*V_estm)).^2     da  dr];
    val_A_Cm_n = [ones(N,1)   alpha_estm  alpha_estm.^2   q_estm*c./V_estm        (q_estm*c./V_estm).^2       de                    Tc1];
    val_A_Cn_n = [ones(N,1)   beta_estm   beta_estm.^2    p_estm*b./(2*V_estm)    (p_estm*b./(2*V_estm)).^2   r_estm*b./(2*V_estm)  (r_estm*b./(2*V_estm)).^2     da  dr];
    
    % Obtaning the aerodynamic forces and moments from validation dataset
    CX_OLS_n = val_A_CX_n * Median_CX_param_n;
    CY_OLS_n = val_A_CY_n * Median_CY_param_n;
    CZ_OLS_n = val_A_CZ_n * Median_CZ_param_n;
    Cl_OLS_n = val_A_Cl_n * Median_Cl_param_n;
    Cm_OLS_n = val_A_Cm_n * Median_Cm_param_n;
    Cn_OLS_n = val_A_Cn_n * Median_Cn_param_n;
    

    data_C_OLS_n = [CX_OLS_n  CY_OLS_n CZ_OLS_n Cl_OLS_n Cm_OLS_n Cn_OLS_n];
    PlotAerodynamicModel(t,data_C_true, data_C_OLS, data_C_OLS_n, coeff_name, val_dataset(iter),true);
    
    data_C_OLS_n_res = data_C_true - data_C_OLS_n;
    PlotResidual(t,data_C_OLS_n_res,'OLS Model 2 Residual',data_C_OLS_res,'OLS Residual',coeff_name,val_dataset(iter))

end
end