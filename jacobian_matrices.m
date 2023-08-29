%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobian function to linearize and then obtain the non-linear system
% matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DFx, DHx, F, G, H] = jacobian_matrices(g, ...
    A_xm, A_ym, A_zm, p_m, q_m, r_m, ...
    x, y, z, u, v, w, phi, theta, psi, u_wind, v_wind, w_wind, ...
    bias_xr, bias_yr, bias_zr, bias_pr, bias_qr, bias_rr)

	% ORDER OF STATES
    states = [x y z u v w phi theta psi u_wind v_wind w_wind...
                bias_xr bias_yr bias_zr bias_pr bias_qr bias_rr];
    % inputs: A_xm, A_ym, A_zm, p_m, q_m, r_m
	% noise: wx, wy, wz, wp, wq, wr
 
%     syms x y z u v w p_m q_m r_m phi theta psi ...
%         A_xm A_ym A_zm bias_xr bias_yr bias_zr ...
%         bias_pr bias_qr bias_rr u_wind v_wind w_wind...
%         wx wy wz wp wq wr 
    
    cth = cos(theta);
    cps = cos(psi);
    cph = cos(phi);
    sth = sin(theta);
    sps = sin(psi);
    sph = sin(phi);
    tth = tan(theta);
    
    F = [ ...
        ... x(dot)
        (u*cth + (v*sph+w*cph)*sth)*cps - (v*cph-w*sph)*sps + u_wind;
        ... y(dot)
        (u*cth + (v*sph+w*cph)*sth)*sps + (v*cph-w*sph)*cps + v_wind;
        ... z(dot)
        -u*sth + (v*sph+w*cph)*cth + w_wind;

        ... u(dot)
        (A_xm - bias_xr) - g*sth + (r_m - bias_rr)*v - (q_m - bias_qr)*w;
        ... v(dot)
        (A_ym - bias_yr) + g*cth*sph + (p_m - bias_pr)*w - (r_m - bias_rr)*u;
        ... w(dot)
        (A_zm - bias_zr) + g*cth*cph + (q_m - bias_qr)*u - (p_m - bias_pr)*v;

        ... phi(dot)
        (p_m - bias_pr) + (q_m - bias_qr)*sph*tth + (r_m - bias_rr)*cph*tth;
        ... theta(dot)
        (q_m - bias_qr)*cph - (r_m - bias_rr)*sph;
        ... psi(dot)
        (q_m - bias_qr)*(sph/cth) + (r_m - bias_rr)*(cph/cth);

        0; % u_wind (dot)
        0; % v_wind (dot)
        0; % w_wind (dot)

        0; % bias_xr (dot)
        0; % bias_yr (dot)
        0; % bias_zr (dot)
        0; % bias_pr (dot)
        0; % bias_qr (dot)
        0]; % bias_rr (dot)
    
    G = [  0   0   0   0   0            0;
           0   0   0   0   0            0;
           0   0   0   0   0            0;
          -1   0   0   0   w           -v;
           0  -1   0  -w   0            u;
           0   0  -1   v  -u            0;
           0   0   0   -1  -sph*tth   -cph*tth;
           0   0   0   0   -cph         sph;
           0   0   0   0   -sph/cth   -cph/cth;
           0   0   0   0   0            0;
           0   0   0   0   0            0;
           0   0   0   0   0            0;
           0   0   0   0   0            0;
           0   0   0   0   0            0;
           0   0   0   0   0            0;
           0   0   0   0   0            0;
           0   0   0   0   0            0;
           0   0   0   0   0            0];
    
    
    H = [   x;
            y;
            z;
           (u*cth + (v*sph + w*cph)*sth)*cps - (v*cph - w*sph)*sps + u_wind;
           (u*cth + (v*sph + w*cph)*sth)*sps + (v*cph - w*sph)*cps + v_wind;
           -u*sth + (v*sph + w*cph)*cth + w_wind;
            phi;
            theta;
            psi;
            sqrt(u^2 + v^2 + w^2);
            atan(w/u);
            atan(v/(sqrt(u^2 + w^2)))];
    
    DFx = jacobian(F,states);
    DHx = jacobian(H,states);
end


