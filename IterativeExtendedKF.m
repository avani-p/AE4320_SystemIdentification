%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IterativeExtendedKF() is used to perform state estimation using
% (Iterative) Extended Kalman Filter on the input and measurement data
% provided.
%
% Gives output of the estimated state XX_k1_k1, Estimated covariance matrix
% diagonal elements PP_k1_k1, predicted state ZZ_pred, and info on the
% iterations performed by the IEKF iter_info.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XX_k1_k1, PP_k1_k1, ZZ_pred, iter_info] = IterativeExtendedKF(n, nm, ...
    t, g, inputs, Z_k, x_k1_k1, P_k1_k1, Q, R, doIEKF, epsilon, maxIterations)
    
    N = length(t);
    dt = t(2)-t(1);
    
    % Inputs
    A_xm = inputs(:,1);
    A_ym = inputs(:,2);
    A_zm = inputs(:,3);
    p_m  = inputs(:,4);
    q_m  = inputs(:,5);
    r_m  = inputs(:,6);

    % Initialize Iterated Extended Kalman filter
    XX_k1_k1                = zeros(n, N); % final states
    PP_k1_k1                = zeros(n, N);
    ZZ_pred                 = zeros(nm, N);
    
    iter_info.IEKFitcount   = zeros(N, 1);
    iter_info.IterError     = zeros(N, maxIterations);
    
    % Run the KALMAN FILTER
    tic;
    
    t_k             = 0; 
    t_k1            = dt;
    % for all N samples of the data
    for k = 1:N
        % Initialization of current iteration
        A_xm_k = A_xm(k);   A_ym_k = A_ym(k);   A_zm_k = A_zm(k);
        p_m_k  = p_m(k);    q_m_k  = q_m(k);    r_m_k  = r_m(k);

        % STEP 1: One Step Ahead Prediction (using ODE45)
        % x(k+1|k) (prediction)
        [~, x_k1_k] = ode45(@(ode_t, states) kf_sys_dynamics(ode_t, states,...
            g, A_xm_k, A_ym_k, A_zm_k, p_m_k, q_m_k, r_m_k),...
        [t_k t_k1], x_k1_k1);
        x_k1_k = x_k1_k(end,:)';
        
        % STEP 2: Calculate Jacobian
        [DFx, G, H, DHx] = kf_calc_jacobian(g, p_m_k, q_m_k, r_m_k,...
            x_k1_k(1),  x_k1_k(2),  x_k1_k(3),  x_k1_k(4),  x_k1_k(5), ...
            x_k1_k(6),  x_k1_k(7),  x_k1_k(8),  x_k1_k(9),  x_k1_k(10), ...
            x_k1_k(11), x_k1_k(12), x_k1_k(16), x_k1_k(17), x_k1_k(18), false);
        
        % z(k+1|k) (predicted observation)
        z_k1_k = H;
        ZZ_pred(:,k) = z_k1_k;% store this observation prediction, since later prediction observations
                              % are corrected using the actual observation
        % STEP 3: Discretize state transition and input matrix
        [Phi_k1k, Gamma_k1k] = c2d(DFx, G, dt);

        % STEP 4: Calculation of Covariance Matrix of State Prediction
        % Error
        P_k1_k    = Phi_k1k*P_k1_k1*(Phi_k1k.') + Gamma_k1k*Q*(Gamma_k1k.'); 
    
        % EKF or IEKF
        if(~doIEKF)
            
            % STEP 5: Calculation of KALMAN Gain Matrix
            % P_zz(k+1|k) (covariance matrix of innovation)
            P_zz            = (DHx*P_k1_k * (DHx.') + R);% covariance matrix of observation error
            
            % K(k+1) (gain)
            K_ekf           = mrdivide(P_k1_k *(DHx.'),P_zz);
                        
            % STEP 6: Measurement Update optimal state x(k+1|k+1) 
            x_k1_k1         = x_k1_k + K_ekf*(Z_k(:,k) - z_k1_k); 

            % STEP 7: Covariance Matrix of State Estimation Error
            P_k1_k1         = (eye(n) - K_ekf*DHx) * P_k1_k;  
        
        else
            % do the iterative part
            eta_i1    = x_k1_k;
            err     = 2*epsilon;
    
            iter    = 0;

            % Iterate till the error is less than performance bound
            while (err >= epsilon)
                % Break if iterations exceed the max value
                if (iter > maxIterations)
                    fprintf('Terminating IEKF: exceeded max iterations (%d)\n', maxIterations);
                    break
                end
                iter    = iter + 1;
                eta_i    = eta_i1;

                % STEP 5: Recalculate Jacobian of Measurement Equation
                % Construct the Jacobian H = d/dx(h(x))) with h(x) 
                % the observation model transition matrix 
                [~, ~, ~, DHx_new] = kf_calc_jacobian(g, p_m_k, q_m_k, r_m_k,...
                    eta_i(1),  eta_i(2),  eta_i(3),  eta_i(4),  eta_i(5), ...
                    eta_i(6),  eta_i(7),  eta_i(8),  eta_i(9),  eta_i(10), ...
                    eta_i(11), eta_i(12), eta_i(16), eta_i(17), eta_i(18), true);

                % Check observability of state
                % if (k == 1 && iter == 1)
                %     rankHF  = kf_calc_ObsRank(DHx_new, DFx);
                %     if (rankHF < n)
                %         warning('The current state is not observable; rank of Observability Matrix is %d, should be %d', rankHF, n);
                %     end
                % end
                % 
                % STEP 6: Calculation of KALMAN gain matrix
                % Observation and observation error predictions
                P_zz        = (DHx_new*P_k1_k*(DHx_new.') + R); % covariance matrix of observation error
                
                % calculate the Kalman gain matrix
                K_iekf      = mrdivide(P_k1_k*(DHx_new.'),P_zz);
             
                % STEP 7: Measurement and state estimate update
                % new observation state
                [~, ~, z_p, ~] = kf_calc_jacobian(g, p_m_k, q_m_k, r_m_k,...
                    eta_i(1),  eta_i(2),  eta_i(3),  eta_i(4),  eta_i(5), ...
                    eta_i(6),  eta_i(7),  eta_i(8),  eta_i(9),  eta_i(10), ...
                    eta_i(11), eta_i(12), eta_i(16), eta_i(17), eta_i(18), true);

                eta_i1        = x_k1_k + K_iekf*(Z_k(:,k) - z_p - DHx_new*(x_k1_k - eta_i));
                err         = norm((eta_i1 - eta_i), inf)/norm(eta_i, inf);
                iter_info.IterError(k,iter) = err;
            end
            
            iter_info.IEKFitcount(k)  = iter;
            % STEP 8:  Covariance Matrix of State Estimation Error
            P_k1_k1         = (eye(n) - K_iekf*DHx_new) * P_k1_k * ((eye(n) - K_iekf*DHx_new).')...
                                + K_iekf*R*(K_iekf.');  
    
            % Next Time Step
            x_k1_k1         = eta_i1;
           
        end

        t_k             = t_k1; 
        t_k1            = t_k + dt;

        % store results
        XX_k1_k1(:,k)     = x_k1_k1;
        PP_k1_k1(:,k)     = diag(P_k1_k1);
    end
    
    IEKFtime = toc;    
    fprintf('IEKF: Completed run with %d samples in %2.2f seconds.\n', N, IEKFtime);
end
