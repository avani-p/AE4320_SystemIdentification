clc
clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doParamEst = false; % false for KALMAN FILTER of one dataset only

if(~doParamEst)
    % For Testing the KALMAN FILTER ONLY
    % Part 2 and 3 of assignment
    datasets = {'da3211_2.mat','dadoublet_1.mat','dr3211_1.mat',...
        'dr3211_2.mat','de3211_1.mat','dedoublet_1.mat'};
    datasets = datasets(4); % Change the number as per your choice
    N_dataset = length(datasets);
    showPlots = true;
    saveWorkspace = false;
    N_repeat = 1;
else
    % For creating datasets of Parameter Estimation
    datasets = {'da3211_2.mat','dadoublet_1.mat','dr3211_1.mat',...
        'dr3211_2.mat','de3211_1.mat','dedoublet_1.mat'};
    % datasets = {'da3211_2.mat','dadoublet_1.mat','dr3211_1.mat',...
    %     'dr3211_2.mat','de3211_1.mat','dedoublet_1.mat'};
    N_dataset = length(datasets);
    showPlots = false;
    saveWorkspace = true;
    N_repeat = 20;  % Specify number of times that each dataset should repeat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONSTANTS
% Output Noise Standard Deviation Values (sigma)
sigma_gps_pos   = 2.5;              % [m]   GPS position
sigma_gps_vel   = 0.02;             % [m/s] GPS velocity
sigma_gps_att   = deg2rad(0.05);    % [rad] GPS attitudes
sigma_tas       = 0.1;              % [m/s] True airspeed, V
sigma_alpha     = deg2rad(0.1);     % [rad] Angle of Attack, α
sigma_beta      = deg2rad(0.1);     % [rad] Sideslip Angle, β

% Real Biases of Accelerometers and Rate Gyros
l_xr = 0.02;    % [m/s^2]
l_yr = 0.02;    % [m/s^2]
l_zr = 0.02;    % [m/s^2]
l_pr = deg2rad(0.003);   % [rad/s]
l_qr = deg2rad(0.003);   % [rad/s]
l_rr = deg2rad(0.003);   % [rad/s]

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

% Simulation Parameters
n               = 18;       % number of states
nm              = 12;       % number of measurements
m               = 6;        % number of inputs
dt              = 0.01;     % time step (s)
epsilon         = 1e-10;    % performance bound
doIEKF          = 1;        % set 1 for IEKF and 0 for EKF
maxIterations   = 100;


%%%%%%%%%%%%% RUN IEKF FOR ALL DATASETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for number = 1:N_dataset
    % Load data
    dataset = string(datasets(number));
    load(append('simdata2023/',dataset));    
    clear string

    N = length(t);
    
    %% Part 2 : Data Pre-processing
    
    %% 2.1 Generate aircraft position using the given simulated airspeed
    %      components in the navigation frame and assumed wind velocity. 
    
    % Given u_n, v_n and w_n are airspeed in NED (navigation) frame 
    % Converting u,v, w from NED (navigation) frame to body frame
    u_body = zeros(N,1);
    v_body = zeros(N,1);
    w_body = zeros(N,1);
    
    % Calculation of airspeed in BODY frame
    for i = 1:N
        [u_body(i), v_body(i), w_body(i)] = ned2body(phi(i), theta(i), ...
                                psi(i), u_n(i), v_n(i), w_n(i));
    end
    
    % Wind components in the earth reference frame (constant wind)
    W_xE =  2;
    W_yE = -8;
    W_zE =  1;
    
    % Correcting Velocity for windspeed --- Ground Speed Components
    % NED Frame
    u_gnd = u_n + W_xE;
    v_gnd = v_n + W_yE;
    w_gnd = w_n + W_zE;
    V_gnd = sqrt(u_gnd.^2 + v_gnd.^2 + w_gnd.^2);
    
    % Calculating position vectors by integration
    pos_X = cumtrapz(t,u_gnd);
    pos_Y = cumtrapz(t,v_gnd);
    pos_Z = cumtrapz(t,w_gnd);
    % position = cumtrapz(t,V_gnd);
    
    if (showPlots)
        % Position Plot
        figure()
        plot3(pos_X,-pos_Y,-pos_Z+h,LineWidth=2)
        xlabel('X (North)'); ylabel('-Y (West)'); zlabel('-Z (+7500) (Up)')
        title('Position of Aircraft wrt NED'); grid on
        
        figure()
        plot(t,V_gnd,t,abs(u_gnd),t,abs(v_gnd),t,abs(w_gnd),LineWidth=1.5)
        xlabel('T[s]'); ylabel(' speed [m/s]'); 
        legend('V','|u|','|v|','|w|'); ylim([-5 max(V_gnd)+10])
        title('Ground Velocity of the aircraft')
    end
    %% 2.2 Generate GPS, IMU and airdata measurements using Eqs. (5) and (6) 
    %      using the flight data provided.
    
    for seed = 1:N_repeat % REPEAT OF SAME DATASET FOR DIFFERENT NOISE FOR
                          % MORE DATA SAMPLES FOR PARAMETER ESTIMATION
    % Generating Sensor Noise
    noise_1 = randn(N,1);
    noise_2 = randn(N,1);
    noise_3 = randn(N,1);
    
    noise_gps_pos   = sigma_gps_pos * [noise_1 noise_2 noise_3];
    noise_gps_vel   = sigma_gps_vel * [noise_1 noise_2 noise_3];
    noise_gps_att   = sigma_gps_att * [noise_1 noise_2 noise_3];
    noise_tas       = sigma_tas     * noise_1;
    noise_alpha     = sigma_alpha   * noise_2;
    noise_beta      = sigma_beta    * noise_3;
    
    % Generating Acc and Gyro Noise
    noise_acc       = sigma_acc     * [noise_1 noise_2 noise_3];
    noise_rgyro     = sigma_rgyro   * [noise_1 noise_2 noise_3];
    
    % Equation 5 (GPS Measurements)
    x_gps   = pos_X;
    x_gps_m = x_gps + noise_gps_pos(:,1);
    y_gps   = pos_Y;
    y_gps_m = y_gps + noise_gps_pos(:,2);
    z_gps   = pos_Z;
    z_gps_m = z_gps + noise_gps_pos(:,3);
    
    c_th = cos(theta);
    c_ps = cos(psi);
    c_ph = cos(phi);
    s_th = sin(theta);
    s_ps = sin(psi);
    s_ph = sin(phi);
    
    u_gps    = (u_body.*c_th + (v_body.*s_ph + w_body.*c_ph).*s_th ).*c_ps ...
                  - (v_body.*c_ph - w_body.*s_ph).*s_ps + W_xE*ones(N,1);
    u_gps_m  = u_gps + noise_gps_vel(:,1);
    
    v_gps    = (u_body.*c_th + (v_body.*s_ph + w_body.*c_ph ).*s_th ).*s_ps ...
                  + (v_body.*c_ph - w_body.*s_ph).*c_ps + W_yE*ones(N,1);
    v_gps_m = v_gps + noise_gps_vel(:,2);
    
    w_gps    = -u_body.*s_th + (v_body.*s_ph + w_body.*c_ph).*c_th + W_zE*ones(N,1);
    w_gps_m  = w_gps + noise_gps_vel(:,3);
    
    phi_gps      = phi;
    phi_gps_m    = phi_gps + noise_gps_att(:,1);
    theta_gps    = theta;
    theta_gps_m  = theta_gps + noise_gps_att(:,2);
    psi_gps      = psi;
    psi_gps_m    = psi_gps + noise_gps_att(:,3);
    
    % Equation 6 (airdata measurements)
    V_airdata        = sqrt(u_body.^2 + v_body.^2 + w_body.^2);
    V_airdata_m      = V_airdata + noise_tas;
    
    alpha_airdata    = atan(w_body./u_body);
    alpha_airdata_m  = alpha_airdata + noise_alpha;
    
    beta_airdata     = atan(v_body./(sqrt(u_body.^2 + w_body.^2)));
    beta_airdata_m   = beta_airdata + noise_beta;
    
    % Equation 4 for IMU calculations
    
    A_xm = Ax + l_xr*ones(N,1) + noise_acc(:,1);
    A_ym = Ay + l_yr*ones(N,1) + noise_acc(:,2);
    A_zm = Az + l_zr*ones(N,1) + noise_acc(:,3);
    p_m  = p  + l_pr*ones(N,1) + noise_rgyro(:,1);
    q_m  = q  + l_qr*ones(N,1) + noise_rgyro(:,2);
    r_m  = r  + l_rr*ones(N,1) + noise_rgyro(:,3);
    
    
    % Determine GPS Position
    GPS_pos_X = cumtrapz(t,u_gps);
    GPS_pos_Y = cumtrapz(t,v_gps);
    GPS_pos_Z = cumtrapz(t,w_gps);
    
    % Determine measured GPS Position
    GPSm_pos_X = cumtrapz(t,u_gps_m);
    GPSm_pos_Y = cumtrapz(t,v_gps_m);
    GPSm_pos_Z = cumtrapz(t,w_gps_m);
    
    if (showPlots)
        % Plotting GPS Position
        figure()
        plot3(x_gps,-y_gps,-z_gps+h,LineWidth=1.50,Marker='s',MarkerIndices=1:500:N)
        hold on
        plot3(x_gps_m,-y_gps_m,-z_gps_m+h,LineWidth=1.50,Marker='o',MarkerIndices=200:500:N)
        xlabel('X (North)'); ylabel('-Y (West)'); zlabel('-Z (Up)')
        legend('GPS Pos','GPSm Pos'); grid on
        title('Position Comparison','Original: GPS positions and Measured GPS positions')
        
        figure()
        plot3(x_gps,-y_gps,-z_gps+h,LineWidth=1.50,Marker='s',MarkerIndices=1:500:N)
        hold on
        plot3(GPSm_pos_X,-GPSm_pos_Y,-GPSm_pos_Z+h,LineWidth=1.50,Marker='o',MarkerIndices=200:500:N)
        xlabel('X (North)'); ylabel('-Y (West)'); zlabel('-Z (Up)')
        legend('Og. GPS Pos','Int. GPSm Pos'); grid on
        title('Position Comparison','Original GPS positions vs Integrated meausred GPS positions')
    end
    
    %% 2.3 Design an integrated GPS/IMU/Airdata navigation system using 
    %      the accelerometers, rate gyros, GPS receiver and airdata sensors
    %      by constructing a navigation model with 18 states
    
    
    % constant
    syms sym_g
    
    % 6 inputs
    % inputs = [A_xm,A_ym,A_zm,p_m,q_m,r_m]
    syms sym_A_xm sym_A_ym sym_A_zm
    syms sym_p_m sym_q_m sym_r_m
    
    % 18 states
    % states = [x,y,z,u,v,w,phi,theta,psi,u_wind,v_wind,w_wind,l_xr,l_yr,l_zr,l_pr,l_qr,l_rr]
    syms sym_x sym_y sym_z 
    syms sym_u sym_v sym_w  
    syms sym_phi sym_theta sym_psi
    syms sym_WxE sym_WyE sym_WzE
    syms sym_bias_xr sym_bias_yr sym_bias_zr
    syms sym_bias_pr sym_bias_qr sym_bias_rr 
    
    % 12 measurements
    % measurements = [x,y,z,u,v,w,phi,theta,psi,V,alpha,beta]
    % noise vectors
    syms sym_wx sym_wy sym_wz sym_wp sym_wq sym_wr 
    
    [DFx_sym, DHx_sym, F_sym, G_sym, H_sym] = jacobian_matrices(sym_g, ...
            sym_A_xm, sym_A_ym, sym_A_zm,...                    %INPUTS
            sym_p_m, sym_q_m, sym_r_m, ...                      %INPUTS
            sym_x, sym_y, sym_z, sym_u, sym_v, sym_w, ...       %STATES
            sym_phi, sym_theta, sym_psi,...                     %STATES
            sym_WxE, sym_WyE, sym_WzE, ...                      %STATES
            sym_bias_xr, sym_bias_yr, sym_bias_zr, ...          %STATES
            sym_bias_pr, sym_bias_qr, sym_bias_rr); ...         %STATES
            
    
    %% Part - 3 Iterative Extended Kalman Filter
    
    %% 3.1 Create an Extended or Iterated Extended Kalman filters (EKF/IEKF) 
    %      for this aircraft flight path reconstruction problem. Clearly 
    %      present in the report the different steps in your EKF/IEKF.
    
    % True State
    X_k         = [ x_gps,          y_gps,          z_gps, ...
                    u_body,         v_body,         w_body, ...
                    phi_gps,        theta_gps,      psi_gps, ...
                    W_xE*ones(N,1), W_yE*ones(N,1), W_zE*ones(N,1),...
                    l_xr*ones(N,1), l_yr*ones(N,1), l_zr*ones(N,1),...
                    l_pr*ones(N,1), l_qr*ones(N,1), l_rr*ones(N,1)
                  ]';
    
    % Measurement State
    Z_k       = [   x_gps_m,        y_gps_m,            z_gps_m,...
                    u_gps_m,        v_gps_m,            w_gps_m, ...
                    phi_gps_m,      theta_gps_m,        psi_gps_m, ...
                    V_airdata_m,    alpha_airdata_m,    beta_airdata_m
                ]';
    
    % Initial Estimate
    x_0        = [  x_gps(1);       y_gps(1);       z_gps(1);...
                    u_body(1);      v_body(1);      w_body(1);...
                    phi_gps(1);     theta_gps(1);   psi_gps(1); ...
                    W_xE;           W_yE;           W_zE;...
                    l_xr;           l_yr;           l_zr;...
                    l_pr;           l_qr;           l_rr      
                ]; % initial true state
    
    % Airspeed measurement data in BODY frame
    [u_body_m_est, v_body_m_est, w_body_m_est] = ned2body(phi_gps_m(1), theta_gps_m(1), ...
        psi_gps_m(1), u_gps_m(1), v_gps_m(1), w_gps_m(1));
    
    E_x_0     = [  x_gps_m(1);       y_gps_m(1);       z_gps_m(1);...
                    u_body_m_est;     v_body_m_est;     w_body_m_est;...
                    phi_gps_m(1);     theta_gps_m(1);   psi_gps_m(1); ...
                    0;                0;                0;   ...
                    0;                0;                0;   ...
                    0;                0;                0      
                    % wind components and biases are unknown     
                ]; % initial estimate of optimal value of x_k1_k1
    
    % Initial Estimate of State Vector
    % x(0|0)  = E{x_0}
    x_k1_k1   = E_x_0;
    
    % Initial Estimate of Covariance Matrix for State Esimation Error
    % P(0|0)  = P(0)  
    P_k1_k1   = eye(n,n);  
    
    % Input Noise Covariance Matrix Q
    Q = diag([  sigma_acc^2,     sigma_acc^2,      sigma_acc^2,...
                sigma_rgyro^2,   sigma_rgyro^2,    sigma_rgyro^2]);
    
    % Output Noise Covariance Matrix R
    R = diag([ sigma_gps_pos^2,  sigma_gps_pos^2,  sigma_gps_pos^2,...
               sigma_gps_vel^2,  sigma_gps_vel^2,  sigma_gps_vel^2,...
               sigma_gps_att^2,  sigma_gps_att^2,  sigma_gps_att^2,...
               sigma_tas^2,      sigma_alpha^2,    sigma_beta^2   ]);
    inputs = [A_xm A_ym A_zm p_m q_m r_m];
    
    % Run Kalman Filter
    [XX_k1_k1, ~, ZZ_pred, iter_info] = IterativeExtendedKF(n, nm, t, g, inputs, ...
        Z_k, x_k1_k1, P_k1_k1, Q, R, doIEKF, epsilon, maxIterations);
    
    % calculate state estimation error (in real life this is unknown!)
    EstErr_x    = XX_k1_k1 - X_k; % Convergence
    EstErr_z    = ZZ_pred  - Z_k; % Innovation
    
    %% 3.2 Clearly show that the Kalman filter estimates for the accelerometer and gyro biases are correct
    %      or provide a sound scientific argumentation when they are not.
    if (showPlots)
        PlotBiases(t,XX_k1_k1(end-5:end,:),'IEKF Est. ',X_k(end-5:end,:), 'True ');
    end
    
    
    %% 3.3 Use your EKF/IEKF to estimate all aircraft states. Clearly show the difference between the 'raw'
    %      sensor measurements and the filtered measurements.
    
    % Unfiltered values with System (Kinematics) Dynamics
    X_unf       = zeros(n, N);
    X_unf(:, 1) = E_x_0;
    
    for k = 1:N-1
        Xk_unf    = X_unf(:, k);
        % time
        tk_unf    = t(k);
        tk1_unf   = t(k+1);
        % Acceleration Input
        A_xm_unf  = A_xm(k);
        A_ym_unf  = A_ym(k);
        A_zm_unf  = A_zm(k);
        % Angular Rates Input
        p_m_unf   = p_m(k);
        q_m_unf   = q_m(k);
        r_m_unf   = r_m(k);
    
        [~, Xk1_unf] = ode45(@(ode_t, states) kf_sys_dynamics(ode_t, states, g,...
            A_xm_unf, A_ym_unf, A_zm_unf, p_m_unf, q_m_unf, r_m_unf),...
            [tk_unf tk1_unf], Xk_unf');
        X_unf(:, k+1) = Xk1_unf(end, :)';
    end
        
    
    if (showPlots)
        PlotStates('States',t,XX_k1_k1,' IEKF Est.', X_unf,' Raw', X_k, ' True')
        PlotStates('Measurements',t,ZZ_pred, ' IEKF Meas.', Z_k, ' Raw Sensor Meas.')
    end
    
    if (~saveWorkspace) % IF NEED TO SAVE WORKSPACE -- this portion is not needed
        %% 3.4 Clearly demonstrate how the estimates of the aircrafts states and IMU sensor biases change if
        %      standard deviations for the angle of attack and side slip angle measurement noises are
        %      increased to 2.0 degrees and 3.5 degrees, respectively. 
        
        % New noice with  new standard deviations
        noise_alpha_n     = deg2rad(2.0) * noise_2;
        noise_beta_n      = deg2rad(3.5) * noise_3;
        
        % New measurement data with measurement noise
        
        alpha_airdata_m_n  = alpha_airdata + noise_alpha_n;
        beta_airdata_m_n   = beta_airdata + noise_beta_n;
        
        
        % Measurement State with new data
        Z_k_n       = [ x_gps_m,        y_gps_m,            z_gps_m,...
                        u_gps_m,       v_gps_m,             w_gps_m, ...
                        phi_gps_m,      theta_gps_m,        psi_gps_m, ...
                        V_airdata_m,    alpha_airdata_m_n,  beta_airdata_m_n
                      ]';
        
        % Output Noise Covariance Matrix R for new standard deviations
        R_n = diag([ sigma_gps_pos^2,  sigma_gps_pos^2,  sigma_gps_pos^2,...
                   sigma_gps_vel^2,  sigma_gps_vel^2,  sigma_gps_vel^2,...
                   sigma_gps_att^2,  sigma_gps_att^2,  sigma_gps_att^2,...
                   sigma_tas^2,      deg2rad(2.0)^2,   deg2rad(3.5)^2   ]);
        
        
        % Run Kalman Filter with new data
        [XX_k1_k1_n, ~, ~, ~] = IterativeExtendedKF(n, nm, t, g, inputs,...
            Z_k_n, x_k1_k1, P_k1_k1, Q, R_n, doIEKF, epsilon, maxIterations);
        
        if (showPlots)
            PlotStates('States',t,XX_k1_k1_n(1:end-6,:),' IEKF Est. (α,β Noise)', XX_k1_k1(1:end-6,:),' IEKF Est. (OG)', X_k(1:end-6,:), ' True')
            PlotBiases(t,XX_k1_k1_n(end-5:end,:),'IEKF Est. (α,β Noise) ', XX_k1_k1(end-5:end,:), 'IEKF Est. (OG) ')
        end

        %% 3.5 Show (graphically) that your Kalman filter functions properly by discussing 
        %     convergence and filter innovation.
        
       
        PlotError('Convergence',t,EstErr_x);
        PlotError('Innovation',t,EstErr_z);

    end


    % Discarding first 2 seconds of data because the bias estimates are
    % off during this time.
    start_idx = find(t==2);

    % Extract Required Estimated States from IEKF result
    u_estm      = XX_k1_k1(4,:)';
    v_estm      = XX_k1_k1(5,:)';
    w_estm      = XX_k1_k1(6,:)';
    l_xr_estm   = XX_k1_k1(13,:)';
    l_yr_estm   = XX_k1_k1(14,:)';
    l_zr_estm   = XX_k1_k1(15,:)';
    l_pr_estm   = XX_k1_k1(16,:)';
    l_qr_estm   = XX_k1_k1(17,:)';
    l_rr_estm   = XX_k1_k1(18,:)';
    
    % Estimated Value of Accelerations and Angular Rates
    Ax_estm = A_xm - l_xr_estm;
    Ay_estm = A_ym - l_yr_estm;
    Az_estm = A_zm - l_zr_estm;
    
    p_estm  = p_m  - l_pr_estm;
    q_estm  = q_m  - l_qr_estm;
    r_estm  = r_m  - l_rr_estm;
    
    % Estimated V, α and β from measurements
    V_estm      = ZZ_pred(10,:)';
    alpha_estm  = ZZ_pred(11,:)';
    beta_estm   = ZZ_pred(12,:)';
    
    % Estimate Angular Acceleration Components ̇p, ̇q, ̇r
    p_dot_estm  = gradient(smoothdata(p_estm),t);
    q_dot_estm  = gradient(smoothdata(q_estm),t);
    r_dot_estm  = gradient(smoothdata(r_estm),t);
    
    %%%%%%%%%%%%%%%%% Split ID and Validation Data %%%%%%%%%%%%%%
    % 35 percent identification data
    % 65 percent validation data
    id_idx = [(floor(N*0.05):floor(N*0.25)) (floor(0.5*N):floor(0.6*N)) (floor(0.85*N):floor(0.9*N))];
    val_idx = setdiff(1:N,id_idx);

    %%%%%%%%%%%%%% SAVE WORKSPACE FOR DATASETS %%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveWorkspace
        folder = 'dataWorkspace';
        workspacename = append(folder,'/workspace_',string(number),string(seed),'_',dataset);
        save(workspacename)
    end
    end %% END OF REPEAT OF SAME DATASET
end
