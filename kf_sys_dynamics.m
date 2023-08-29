%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xdot = kf_sys_dynamics(x) Calculates the system dynamics equation f(x,u,t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = kf_sys_dynamics(t, x, g, A_xm, A_ym, A_zm, p_m, q_m, r_m)

    u       = x(4);
    v       = x(5);
    w       = x(6);
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    u_wind  = x(10);
    v_wind  = x(11);
    w_wind  = x(12);
    bias_xr = x(13);
    bias_yr = x(14);
    bias_zr = x(15);
    bias_pr = x(16);
    bias_qr = x(17);
    bias_rr = x(18);

    cth = cos(theta);
    sth = sin(theta);
    cph = cos(phi);
    sph = sin(phi);
    cps = cos(psi);
    sps = sin(psi);
    tth = tan(theta);

    % system dynamics go here!
    xdot = [ ...
        ... x(dot)
        (u.*cth + (v.*sph+w.*cph).*sth).*cps - (v.*cph-w.*sph).*sps + u_wind;
        ... y(dot)
        (u.*cth + (v.*sph+w.*cph).*sth).*sps + (v.*cph-w.*sph).*cps + v_wind;
        ... z(dot)
        -u.*sth + (v.*sph+w.*cph).*cth + w_wind;

        ... u(dot)
        (A_xm - bias_xr) - g.*sth + (r_m - bias_rr).*v - (q_m - bias_qr).*w;
        ... v(dot)
        (A_ym - bias_yr) + g.*cth.*sph + (p_m - bias_pr).*w - (r_m - bias_rr).*u;
        ... w(dot)
        (A_zm - bias_zr) + g.*cth.*cph + (q_m - bias_qr).*u - (p_m - bias_pr).*v;

        ... phi(dot)
        (p_m - bias_pr) + (q_m - bias_qr).*sph.*tth + (r_m - bias_rr).*cph.*tth;
        ... theta(dot)
        (q_m - bias_qr).*cph - (r_m - bias_rr).*sph;
        ... psi(dot)
        (q_m - bias_qr).*(sph./cth) + (r_m - bias_rr).*(cph./cth);

        0; % u_wind (dot)
        0; % v_wind (dot)
        0; % w_wind (dot)

        0; % bias_xr (dot)
        0; % bias_yr (dot)
        0; % bias_zr (dot)
        0; % bias_pr (dot)
        0; % bias_qr (dot)
        0]; % bias_rr (dot)    

end
