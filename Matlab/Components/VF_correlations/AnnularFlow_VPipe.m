function out = AnnularFlow_VPipe(G, P_su, T_su, C_oil, D, fluid_ref, fluid_lub, fit_ratio_rho, fit_DTP_zeta)
% RDickes - 11/04/2018

% Annular flow as described by Sithi.
% Oil-refrigerant solubility obtained by Grebner model fitted on Zhelty data
% Valid for vertical annular flow

%% NOMENCLATURE :

% C_ref =     Mdot_ref/Mdot_tot;  --> REF mass fraction in TOTAL flow
% C_oil =     Mdot_oil/Mdot_tot;  --> OIL mass fraction in TOTAL flow
% C_rv =     Mdot_rv/Mdot_tot;  --> VAP REF mass fraction in TOTAL flow
% C_rl =     Mdot_rl/Mdot_tot;  --> VAP REF mass fraction in TOTAL flow
% zeta_ref =  Mdot_rel/Mdot_liq;  --> REF mass fraction in LIQUID flow
% zeta_oil =  Mdot_oil/Mdot_liq;  --> OIL mass fraction in LIQUID flow

%% Transfer of inputs;
out.G = G;
out.T_su = T_su;
out.P_su = P_su;
out.D = D;
out.R = D/2;
out.fluid_ref = fluid_ref;
out.fluid_lub = fluid_lub;

%% Flow fraction factors and mixture properties
out.Tsat_pure = CoolProp.PropsSI('T', 'P', P_su, 'Q', 0, fluid_ref);
out.Tbubble_min = R245fa_POE_Tbubble(1-C_oil, P_su, out.Tsat_pure);
[~, ~, ~, ~, out.zeta_r, out.x, out.C_rl, out.C_rv] = TP_solubMixt(T_su, P_su, C_oil, fluid_ref, fluid_lub, out.Tbubble_min, out.Tsat_pure, fit_DTP_zeta);
out.C_oil = C_oil;
out.zeta_oil = 1-out.zeta_r;
out.res_zeta_r = 0;
[out.rho_vap, out.rho_rl, out.rho_oil, out.rho_liq, out.K] = R245fa_POE_density(T_su, P_su, out.zeta_r, fluid_ref, fluid_lub, fit_ratio_rho, out.Tbubble_min, out.Tsat_pure);
[out.mu_vap, out.mu_rl, out.mu_oil, out.mu_liq] = R245fa_POE_viscosity(T_su, P_su, C_oil, out.zeta_r, out.rho_liq, out.rho_oil, out.rho_rl, fluid_ref, fluid_lub, out.Tbubble_min, out.Tsat_pure);

%[out.Tbubble_min, out.Tsat_pure, ~, ~, ~, ~,  out.C_rl,  out.C_rv,  out.C_oil,  out.x,  out.zeta_r,  out.res_zeta_r, out.zeta_oil,   out.rho_vap,  out.rho_rl,  out.rho_oil,  out.rho_liq,  out.K,  out.mu_vap,  out.mu_rl,  out.mu_oil,  out.mu_liq] = R245fa_POE_mixtProp(T_su, P_su, C_oil, fluid_ref, fluid_lub, fit_ratio_rho);
out.C_r = 1- C_oil;
out.nu_liq = out.mu_liq/out.rho_liq;
out.nu_vap = out.mu_vap/out.rho_vap;

%% Flow distribution
out.G_liq = G*(1-out.C_rv);
out.G_vap = G*out.C_rv;
out.j_liq = out.G_liq/out.rho_liq;
out.j_vap = out.G_vap/out.rho_vap;
out.Mdot_vap = out.G_vap*(pi*out.R^2);
out.Mdot_liq = out.G_liq*(pi*out.R^2);

%% Deterministic annular flow model
f = @(x) res_LiquidFilm_detModel(x, out.D, out.G_vap, out.rho_vap, out.mu_vap, out.nu_vap, out.G_liq, out.mu_liq, out.rho_liq,out.Mdot_liq);
delta_by_R_min = 0; 
delta_by_R_max = 0;
stop = 0;
while not(stop)
    fmax = f(delta_by_R_max);
    if fmax > 0 && delta_by_R_max <0.9
        delta_by_R_min = delta_by_R_max;
        delta_by_R_max = delta_by_R_max +1e-3;
    else
        stop = 1;
    end
end
out.delta_by_R = fzero(f,[delta_by_R_min delta_by_R_max]);
[out.res_Mdot_liq, out.delta, out.alpha, out.u_vap, out.Re_vap, out.Re_liq, out.f_s, out.tau, out.res_tau, out.delta_plus, out.f_i, out.tau_i, out.dpdz, out.Mdot_liq_bis] = LiquidFilm_detModel(out.delta_by_R, D, out.G_vap, out.rho_vap, out.mu_vap, out.nu_vap, out.G_liq, out.mu_liq, out.rho_liq, out.Mdot_liq);

%% Other output Retentation
out.u_liq = out.G_liq/(out.rho_liq*(1-out.alpha));
out.slip = out.u_vap/out.u_liq;
out.oil_retention = pi*(out.R^2 - (out.R-out.delta)^2)*out.rho_liq*out.zeta_oil; %kg/m
out.alpha_bis = 1/(1+(1-out.C_rv)/out.C_rv*(out.rho_vap/out.rho_liq)*out.slip);
if max(abs([out.res_Mdot_liq out.res_tau  out.res_zeta_r])) < 1e-5
    out.flag = 1;
    out.Validity_delta_by_R = 1;
    out.Validity_Re_vap = 1;
    out.Validity_Re_liq = 1;
    if (0.5*out.delta_by_R > 0.07 || 0.5*out.delta_by_R < 0.05)
        out.flag = 0.5;
        out.Validity_delta_by_R = 0;
    end
    if (out.Re_vap > 210000 || out.Re_vap < 48000)
        out.flag = 0.5;
        out.Validity_Re_vap = 0;
    end
    if  (out.Re_liq > 10 || out.Re_liq < 0.3)
        out.flag = 0.5;
        out.Validity_Re_liq = 0;
    end
else
    out.flag = -1;
    disp('Prob annular flow')
end
out = orderfields(out);

end