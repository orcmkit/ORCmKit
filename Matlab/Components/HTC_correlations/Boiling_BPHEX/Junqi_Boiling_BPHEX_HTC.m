function [h_boiling, Nu, flag] = Junqi_Boiling_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v,  i_fg, G, DT_log, Qdot, honv_h, Dh, disp_flag)

% Boiling heat transfer correlation published by Junqi et. al for
% characterizing boiling process in BPHEX (validated with R245fa only!)
% sources : "Experimental investigation on heat transfer characteristics 
% of plat heat exchanger applied in organic Rankine cycle (ORC)"

Re_min = 250;
Re_max = 2500;
flag = 1;

G_eq = G * ( (1 - x) + x * (rho_l/rho_v)^0.5);
Re_eq = G_eq*Dh/mu_l;
AU_tp = Qdot/DT_log;
Bo_0 = 1;
f_Bo = @(xx) res_iter_Junqi_boiling(xx,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh);
opts = optimoptions('fsolve', 'display', 'none');
[Bo,err_Bo] = fsolve(f_Bo, Bo_0, opts);
[~, Nu, h, ~, ~, ~, ~] = iter_Junqi_boiling(Bo,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh);

if Re_eq < Re_min || Re_eq > Re_max
    if disp_flag
        display(['Junqi boiling: Out of validity range --> Re_eq = ' num2str(Re_eq) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = -1;
end

if err_Bo > 1e-4
    if disp_flag
        display(['Junqi boiling: Wrong boiling number --> err_Bo = ' num2str(err_Bo) ' !!!'])
    end
    flag = -5;
end

h_boiling = h;


end

function res_Bo = res_iter_Junqi_boiling(Bo_g,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
[res_Bo, ~, ~, ~, ~, ~, ~] = iter_Junqi_boiling(Bo_g,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh);
end

function [res_Bo, Nu, h, U, A_tp, q, Bo] = iter_Junqi_boiling(Bo_g,Re_eq, Pr_l, k_l, honv_h, AU_tp, Qdot, G_eq, i_fg, Dh)
Bo_g = max(Bo_g, 0);
Nu = 2.64*Re_eq^0.815*Bo_g^0.343*Pr_l^0.333;
h = Nu*k_l/Dh;
U = (1/h +  1/honv_h)^-1;
A_tp = AU_tp/U;
q = Qdot/A_tp;
Bo = q/(G_eq*i_fg);
res_Bo = (Bo-Bo_g)/Bo_g;
end