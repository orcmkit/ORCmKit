function [h_boiling, Nu, flag] = Han_Boiling_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v,  i_fg, G, DT_log, Qdot, honv_h, Dh, theta, pitch_co, disp_flag)
% Boiling heat transfer correlation published by Han et. al in
% "Experiments on the characteristics of evaporation of R410A in brazed plate heat exchangers with different geometric configurations"

flag = 1;
G_eq = G * ( (1 - x) + x * (rho_l/rho_v)^0.5);
Re_eq = G_eq*Dh/mu_l;
AU_tp = Qdot/DT_log;
Ge1 = 2.81*(pitch_co/Dh)^-0.041*(theta)^(-2.83);
Ge2 = 0.746*(pitch_co/Dh)^-0.082*(theta)^(0.61);
Bo_0 = 1;
f_Bo = @(xx) res_iter_Han_boiling(xx, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh);
opts = optimoptions('fsolve', 'display', 'none');
[Bo,err_Bo] = fsolve(f_Bo, Bo_0, opts);
[~, Nu, h, ~, ~, ~, ~] = iter_Han_boiling(Bo, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh);

if err_Bo > 1e-4
    if disp_flag
        display(['Han boiling: Wrong boiling number --> err_Bo = ' num2str(err_Bo) ' !!!'])
    end
    flag = -1;
end
h_boiling = h;


end

function res_Bo = res_iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
[res_Bo, ~, ~, ~, ~, ~, ~] = iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh);
end

function [res_Bo, Nu, h, U, A_tp, q, Bo] = iter_Han_boiling(Bo_g, Ge1, Ge2,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
Bo_g = max(Bo_g, 0);
Nu = Ge1*Re_eq^Ge2*Bo_g^0.3*Pr_l^0.4;
h = Nu*k_l/Dh;
U = (1/h +  1/honv_h)^-1;
A_tp = AU_tp/U;
q = Qdot/A_tp;
Bo = q/(G_eq*i_fg);
res_Bo = (Bo-Bo_g)/Bo_g;
end