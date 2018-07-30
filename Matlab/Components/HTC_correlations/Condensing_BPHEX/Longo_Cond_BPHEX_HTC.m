function [h_cond, Nu, flag ]  = Longo_Cond_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v, i_fg, G, Dh, DT_wall, phi, L, disp_flag) % verified
% Condensation heat transfer correlation published by Longo et al. in
% "A new computational procedure for refrigerant condensation inside herringbone-type Brazed Plate Heat Exchangers", International Journal of Heat and Mass Transfer, 2015

% RDickes - 20/07/2018
G_eq = G * ( (1 - x) + x * (rho_l/rho_v)^0.5);
Re_eq = G_eq*Dh/mu_l;
g = 9.81;
if Re_eq < 1600
    h_cond = phi*0.943*((k_l^3*rho_l^2*g*i_fg)/(mu_l*DT_wall*L))^0.25;
else
    h_cond = 1.875*phi*k_l/Dh*Re_eq^0.445*Pr_l^0.3333333;
end

Nu = h_cond*Dh/k_l;
flag = 1;
        
end