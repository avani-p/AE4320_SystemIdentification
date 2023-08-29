function [u_NED, v_NED, w_NED] = body2ned(phi, theta, psi, u, v, w)
    
    cth = cos(theta);
    cps = cos(psi);
    cph = cos(phi);
    sth = sin(theta);
    sps = sin(psi);
    sph = sin(phi);
    
    R_Body2NED = [cth*cps     sph*sth*cps - cph*sps      cph*sth*cps + sph*sps;
                  cth*sps     sph*sth*sps + cph*cps      cph*sth*sps - sph*cps;
                     -sth     sph*cth                    cph*cth             ];
    
    U = [u; v; w];
    X = R_Body2NED*U;
    
    u_NED = X(1);
    v_NED = X(2);
    w_NED = X(3);
end