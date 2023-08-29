%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DFx, G, H, DHx] = kf_calc_jacobian(x) 
% Calculates the system dynamics equation h(x,u,t) g(x,u,t)
% and Jacobian DFx and DHx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DFx, G, H, DHx] = kf_calc_jacobian(sym_g, sym_p_m, sym_q_m, sym_r_m,...
    sym_x, sym_y, sym_z, sym_u, sym_v, sym_w, sym_phi, sym_theta, sym_psi,...
    sym_u_wind, sym_v_wind, sym_w_wind, sym_bias_pr, sym_bias_qr, sym_bias_rr, ...
    calc_only_H_DHx)
    
    G = 0;
    DFx = 0;

    if(~calc_only_H_DHx)
    
        G =[   
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0;
            -1,  0,  0,      0,                        sym_w,                       -sym_v;
             0, -1,  0, -sym_w,                            0,                        sym_u;
             0,  0, -1,  sym_v,                       -sym_u,                            0;
             0,  0,  0,     -1, -sin(sym_phi)*tan(sym_theta), -cos(sym_phi)*tan(sym_theta);
             0,  0,  0,      0,                -cos(sym_phi),                 sin(sym_phi);
             0,  0,  0,      0, -sin(sym_phi)/cos(sym_theta), -cos(sym_phi)/cos(sym_theta);
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0;
             0,  0,  0,      0,                            0,                            0
            ];

        DFx = [
            0, 0, 0, cos(sym_psi)*cos(sym_theta), cos(sym_psi)*sin(sym_phi)*sin(sym_theta) - cos(sym_phi)*sin(sym_psi), sin(sym_phi)*sin(sym_psi) + cos(sym_phi)*cos(sym_psi)*sin(sym_theta), sin(sym_psi)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) + cos(sym_psi)*sin(sym_theta)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)),                                                 -cos(sym_psi)*(sym_u*sin(sym_theta) - cos(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi))), - sin(sym_psi)*(sin(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) + sym_u*cos(sym_theta)) - cos(sym_psi)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)), 1, 0, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0, cos(sym_theta)*sin(sym_psi), cos(sym_phi)*cos(sym_psi) + sin(sym_phi)*sin(sym_psi)*sin(sym_theta), cos(sym_phi)*sin(sym_psi)*sin(sym_theta) - cos(sym_psi)*sin(sym_phi), sin(sym_psi)*sin(sym_theta)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)) - cos(sym_psi)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)),                                                 -sin(sym_psi)*(sym_u*sin(sym_theta) - cos(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi))),   cos(sym_psi)*(sin(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) + sym_u*cos(sym_theta)) - sin(sym_psi)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)), 0, 1, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,             -sin(sym_theta),                                          cos(sym_theta)*sin(sym_phi),                                          cos(sym_phi)*cos(sym_theta),                                                                       cos(sym_theta)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)),                                                               - sin(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) - sym_u*cos(sym_theta),                                                                                                                                                         0, 0, 0, 1,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,                           0,                                                sym_r_m - sym_bias_rr,                                                sym_bias_qr - sym_q_m,                                                                                                                              0,                                                                                                                           -sym_g*cos(sym_theta),                                                                                                                                                         0, 0, 0, 0, -1,  0,  0,      0,                        sym_w,                       -sym_v;
            0, 0, 0,       sym_bias_rr - sym_r_m,                                                                    0,                                                sym_p_m - sym_bias_pr,                                                                                              sym_g*cos(sym_phi)*cos(sym_theta),                                                                                                              -sym_g*sin(sym_phi)*sin(sym_theta),                                                                                                                                                         0, 0, 0, 0,  0, -1,  0, -sym_w,                            0,                        sym_u;
            0, 0, 0,       sym_q_m - sym_bias_qr,                                                sym_bias_pr - sym_p_m,                                                                    0,                                                                                             -sym_g*cos(sym_theta)*sin(sym_phi),                                                                                                              -sym_g*cos(sym_phi)*sin(sym_theta),                                                                                                                                                         0, 0, 0, 0,  0,  0, -1,  sym_v,                       -sym_u,                            0;
            0, 0, 0,                           0,                                                                    0,                                                                    0,                      cos(sym_phi)*tan(sym_theta)*(sym_q_m - sym_bias_qr) - sin(sym_phi)*tan(sym_theta)*(sym_r_m - sym_bias_rr),                       cos(sym_phi)*(sym_r_m - sym_bias_rr)*(tan(sym_theta)^2 + 1) + sin(sym_phi)*(sym_q_m - sym_bias_qr)*(tan(sym_theta)^2 + 1),                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,     -1, -sin(sym_phi)*tan(sym_theta), -cos(sym_phi)*tan(sym_theta);
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                  - cos(sym_phi)*(sym_r_m - sym_bias_rr) - sin(sym_phi)*(sym_q_m - sym_bias_qr),                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                -cos(sym_phi),                 sin(sym_phi);
            0, 0, 0,                           0,                                                                    0,                                                                    0,                  (cos(sym_phi)*(sym_q_m - sym_bias_qr))/cos(sym_theta) - (sin(sym_phi)*(sym_r_m - sym_bias_rr))/cos(sym_theta), (sin(sym_phi)*sin(sym_theta)*(sym_q_m - sym_bias_qr))/cos(sym_theta)^2 + (cos(sym_phi)*sin(sym_theta)*(sym_r_m - sym_bias_rr))/cos(sym_theta)^2,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0, -sin(sym_phi)/cos(sym_theta), -cos(sym_phi)/cos(sym_theta);
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                                                                                              0,                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                                                                                              0,                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                                                                                              0,                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                                                                                              0,                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                                                                                              0,                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                                                                                              0,                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                                                                                              0,                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                                                                                              0,                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                            0,                            0;
            0, 0, 0,                           0,                                                                    0,                                                                    0,                                                                                                                              0,                                                                                                                                               0,                                                                                                                                                         0, 0, 0, 0,  0,  0,  0,      0,                            0,                            0
            ];
    
    end

        H = [  
                                                                                                                                                               sym_x;
                                                                                                                                                               sym_y;
                                                                                                                                                               sym_z;
sym_u_wind - sin(sym_psi)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)) + cos(sym_psi)*(sin(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) + sym_u*cos(sym_theta));
sym_v_wind + sin(sym_psi)*(sin(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) + sym_u*cos(sym_theta)) + cos(sym_psi)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi));
                                                                        sym_w_wind - sym_u*sin(sym_theta) + cos(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi));
                                                                                                                                                             sym_phi;
                                                                                                                                                           sym_theta;
                                                                                                                                                             sym_psi;
                                                                                                                                 (sym_u^2 + sym_v^2 + sym_w^2)^(1/2);
                                                                                                                                                   atan(sym_w/sym_u);
                                                                                                                               atan(sym_v/(sym_u^2 + sym_w^2)^(1/2))
          ];


        DHx = [      
            1, 0, 0,                                                                            0,                                                                    0,                                                                            0,                                                                                                                              0,                                                                                               0,                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 1, 0,                                                                            0,                                                                    0,                                                                            0,                                                                                                                              0,                                                                                               0,                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 1,                                                                            0,                                                                    0,                                                                            0,                                                                                                                              0,                                                                                               0,                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0,                                                  cos(sym_psi)*cos(sym_theta), cos(sym_psi)*sin(sym_phi)*sin(sym_theta) - cos(sym_phi)*sin(sym_psi),         sin(sym_phi)*sin(sym_psi) + cos(sym_phi)*cos(sym_psi)*sin(sym_theta), sin(sym_psi)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) + cos(sym_psi)*sin(sym_theta)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)), -cos(sym_psi)*(sym_u*sin(sym_theta) - cos(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi))), - sin(sym_psi)*(sin(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) + sym_u*cos(sym_theta)) - cos(sym_psi)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)), 1, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0,                                                  cos(sym_theta)*sin(sym_psi), cos(sym_phi)*cos(sym_psi) + sin(sym_phi)*sin(sym_psi)*sin(sym_theta),         cos(sym_phi)*sin(sym_psi)*sin(sym_theta) - cos(sym_psi)*sin(sym_phi), sin(sym_psi)*sin(sym_theta)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)) - cos(sym_psi)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)), -sin(sym_psi)*(sym_u*sin(sym_theta) - cos(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi))),   cos(sym_psi)*(sin(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) + sym_u*cos(sym_theta)) - sin(sym_psi)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)), 0, 1, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0,                                                              -sin(sym_theta),                                          cos(sym_theta)*sin(sym_phi),                                                  cos(sym_phi)*cos(sym_theta),                                                                       cos(sym_theta)*(sym_v*cos(sym_phi) - sym_w*sin(sym_phi)),               - sin(sym_theta)*(sym_w*cos(sym_phi) + sym_v*sin(sym_phi)) - sym_u*cos(sym_theta),                                                                                                                                                         0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
            0, 0, 0,                                                                            0,                                                                    0,                                                                            0,                                                                                                                              1,                                                                                               0,                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0,                                                                            0,                                                                    0,                                                                            0,                                                                                                                              0,                                                                                               1,                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0,                                                                            0,                                                                    0,                                                                            0,                                                                                                                              0,                                                                                               0,                                                                                                                                                         1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0,                                    sym_u/(sym_u^2 + sym_v^2 + sym_w^2)^(1/2),                            sym_v/(sym_u^2 + sym_v^2 + sym_w^2)^(1/2),                                    sym_w/(sym_u^2 + sym_v^2 + sym_w^2)^(1/2),                                                                                                                              0,                                                                                               0,                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0,                                       -sym_w/(sym_u^2*(sym_w^2/sym_u^2 + 1)),                                                                    0,                                              1/(sym_u*(sym_w^2/sym_u^2 + 1)),                                                                                                                              0,                                                                                               0,                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, -(sym_u*sym_v)/((sym_u^2 + sym_w^2)^(3/2)*(sym_v^2/(sym_u^2 + sym_w^2) + 1)),      1/((sym_u^2 + sym_w^2)^(1/2)*(sym_v^2/(sym_u^2 + sym_w^2) + 1)), -(sym_v*sym_w)/((sym_u^2 + sym_w^2)^(3/2)*(sym_v^2/(sym_u^2 + sym_w^2) + 1)),                                                                                                                              0,                                                                                               0,                                                                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ];

end