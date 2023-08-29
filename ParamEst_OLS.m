%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ParamEst_OLS() is used to estimate the parameter of the aerodynamic
% moments and forces using sytem inputs, states, and measurements and
% Ordinary Least Square Method.
% 
% Gives output of estimated parameters of - Cx, Cy, Cz, Cl, Cm, Cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OLS_CX_param, OLS_CY_param, OLS_CZ_param, OLS_Cl_param, OLS_Cm_param,...
    OLS_Cn_param, OLS_A_CX, OLS_A_CY, OLS_A_CZ, OLS_A_Cl, OLS_A_Cm, OLS_A_Cn] = ...
    ParamEst_OLS(c, b, C_X, C_Y, C_Z, C_l, C_m, C_n,... % CONSTANTS
    p_estm, q_estm, r_estm, da, de, dr, Tc, ... % SYS INPUTS
    V_estm, alpha_estm, beta_estm, useNewStructure) % SYS OUTPUTS
    
    N = length(p_estm);
    % STEP 2: Formulate Linear regression model structure 
    % p(x,THETA) =  A(x).THETA

    % The aerodynamic forces and moment are given (Equations 15 and 16)
    % CX = CX_0 + CX_alpha*alpha + CX_alpha2*alpha.^2 + CX_q*(q*c/V)     + CX_de*de + CX_Tc*Tc
    % CY = CY_0 + CY_beta*beta   + CY_p*(p*b/(2*V))   + CY_r*(r*b/(2*V)) + CY_da*da + CY_dr*dr
    % CZ = CZ_0 + CZ_alpha*alpha + CZ_q*(q*c/V)       + CZ_de*de         + CZ_Tc*Tc
    
    % Cl = Cl_0 + Cl_beta*beta   + Cl_p*(p*b/(2*V))   + Cl_r*(r*b/(2*V)) + Cl_da*da + Cl_dr*dr
    % Cm = Cm_0 + Cm_alpha*alpha + Cm_q*(q*c/V)       + Cm_de*de         + Cm_Tc*Tc
    % Cn = Cn_0 + Cn_beta*beta   + Cn_p*(p*b/(2*V))   + Cn_r*(r*b/(2*V)) + Cn_da*da + Cn_dr*dr
    
    % So, THETA = CX_0, CX_alpha, CX_alpha2, CX_q, CX_de, CX_Tc for Cx and
    % similarily for others
    % STEP 3: Formulate regression matrix according to equation

    if (~useNewStructure)
        OLS_A_CX = [ones(N,1)   alpha_estm  alpha_estm.^2           q_estm*c./V_estm     de     Tc];
        OLS_A_CY = [ones(N,1)   beta_estm   p_estm*b./(2*V_estm)    r_estm*b./(2*V_estm) da     dr];
        OLS_A_CZ = [ones(N,1)   alpha_estm  q_estm*c./V_estm        de                   Tc];
        
        OLS_A_Cl = [ones(N,1)   beta_estm   p_estm*b./(2*V_estm)    r_estm*b./(2*V_estm) da     dr];
        OLS_A_Cm = [ones(N,1)   alpha_estm  q_estm*c./V_estm        de                   Tc];
        OLS_A_Cn = [ones(N,1)   beta_estm   p_estm*b./(2*V_estm)    r_estm*b./(2*V_estm) da     dr];
    else
        OLS_A_CX = [ones(N,1)   alpha_estm  alpha_estm.^2   q_estm*c./V_estm        (q_estm*c./V_estm).^2       de                    Tc];
        OLS_A_CY = [ones(N,1)   beta_estm   beta_estm.^2    p_estm*b./(2*V_estm)    (p_estm*b./(2*V_estm)).^2   r_estm*b./(2*V_estm)  (r_estm*b./(2*V_estm)).^2     da  dr];
        OLS_A_CZ = [ones(N,1)   alpha_estm  alpha_estm.^2   q_estm*c./V_estm        (q_estm*c./V_estm).^2       de                    Tc];
        
        OLS_A_Cl = [ones(N,1)   beta_estm   beta_estm.^2    p_estm*b./(2*V_estm)    (p_estm*b./(2*V_estm)).^2   r_estm*b./(2*V_estm)  (r_estm*b./(2*V_estm)).^2     da  dr];
        OLS_A_Cm = [ones(N,1)   alpha_estm  alpha_estm.^2   q_estm*c./V_estm        (q_estm*c./V_estm).^2       de                    Tc];
        OLS_A_Cn = [ones(N,1)   beta_estm   beta_estm.^2    p_estm*b./(2*V_estm)    (p_estm*b./(2*V_estm)).^2   r_estm*b./(2*V_estm)  (r_estm*b./(2*V_estm)).^2     da  dr];
    end

    % STEP 4: Formulate Least Squares Estimator
    % THETA_cap = (A'.A)^-1 . A'.y
    OLS_CX_param = mldivide((OLS_A_CX.')*OLS_A_CX, OLS_A_CX.')*C_X;
    OLS_CY_param = mldivide((OLS_A_CY.')*OLS_A_CY, OLS_A_CY.')*C_Y;
    OLS_CZ_param = mldivide((OLS_A_CZ.')*OLS_A_CZ, OLS_A_CZ.')*C_Z;
    
    OLS_Cl_param = mldivide((OLS_A_Cl.')*OLS_A_Cl, OLS_A_Cl.')*C_l;
    OLS_Cm_param = mldivide((OLS_A_Cm.')*OLS_A_Cm, OLS_A_Cm.')*C_m;
    OLS_Cn_param = mldivide((OLS_A_Cn.')*OLS_A_Cn, OLS_A_Cn.')*C_n;

end