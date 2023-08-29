function [u_BODY, v_BODY, w_BODY] = ned2body(phi, theta, psi, u, v, w)
    
    cth = cos(theta);
    cps = cos(psi);
    cph = cos(phi);
    sth = sin(theta);
    sps = sin(psi);
    sph = sin(phi);
    
    R_Body2NED = [cth*cps     sph*sth*cps - cph*sps      cph*sth*cps + sph*sps;
                  cth*sps     sph*sth*sps + cph*cps      cph*sth*sps - sph*cps;
                     -sth     sph*cth                    cph*cth             ];
    
    R_NED2Body = R_Body2NED';
    U = [u; v; w];
    X = R_NED2Body*U;
    
    u_BODY = X(1);
    v_BODY = X(2);
    w_BODY = X(3);
end