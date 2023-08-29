%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ParamEst_RLS_OneEq() is used to estimate the parameter of the aerodynamic
% moments and forces using sytem inputs, states, and measurements 
% and Recursive Least Square Method.
% 
% Gives output of estimated parameters of the polynomial equation based on
% A_matrix and the y values provided. (lamba = forgetting factor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RLS_param = ParamEst_RLS(A, y, lambda,THETA0, P0)
    N = length(y);
    % Initial Estimate of parameters and parameter covariance matrix
    THETA_k1 = THETA0;
    P_k1 = P0;

    for k = 0 : N-1
        THETA_k = THETA_k1;
        P_k    = P_k1;

        % STEP 2: Obtain Measurement data (x_k1 and y_k1)
        % x_k1 not needed as A matrix already contains the function values.
        y_k1 = y(k+1);

        % STEP 3: Formulate Regression Matrix/Vector for new data points
        a_k1 = A(k+1, :);

        % Step 4: Calculate KALMAN gain
        K_k1 = P_k * (a_k1') / (a_k1*P_k*(a_k1') + lambda);

        % Step 5: Calculate Parameter Update
        THETA_k1 = THETA_k + K_k1*(y_k1 - a_k1*THETA_k);

        % Step 6: Calculate Update Parameter Covariance Matrix
        P_k1 = P_k - K_k1*a_k1*P_k;
    end

    RLS_param = THETA_k1;

end