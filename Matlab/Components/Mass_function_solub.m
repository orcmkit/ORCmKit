function out = Mass_function_solub(T_K, P_Pa , V, C_oil, fluid_r, fluid_lub, param)
fit_ratio_rho = param.fit_ratio_rho;
out.Tsat_pure = CoolProp.PropsSI('T', 'P', P_Pa, 'Q', 0, fluid_r);
out.Tbubble_min = R245fa_POE_Tbubble(1-C_oil, P_Pa, out.Tsat_pure);
[~, ~, ~, ~, out.zeta_r, out.x, out.C_rl, out.C_rv] = TP_solubMixt(T_K, P_Pa, C_oil, fluid_r, fluid_lub, out.Tbubble_min, out.Tsat_pure, param.fit_DTP_zeta);
out.C_oil = C_oil;
out.zeta_oil = 1-out.zeta_r;
out.res_zeta_r = 0;
[out.rho_vap, out.rho_rl, out.rho_oil, out.rho_liq, out.K] = R245fa_POE_density(T_K, P_Pa, out.zeta_r, fluid_r, fluid_lub, fit_ratio_rho, out.Tbubble_min, out.Tsat_pure);
[out.mu_vap, out.mu_rl, out.mu_oil, out.mu_liq] = R245fa_POE_viscosity(T_K, P_Pa, C_oil, out.zeta_r, out.rho_liq, out.rho_oil, out.rho_rl, fluid_r, fluid_lub, out.Tbubble_min, out.Tsat_pure);


if out.C_rv > 0 %there is a vapour phase
    out.alpha = VoidFraction_RefLub(out.C_rv, out.rho_vap, out.rho_liq, param);
else % there is no vapour phase
    out.alpha = 0;
end

out.M_rv = V*out.alpha*out.rho_vap;
out.M_rl = V*(1-out.alpha)*out.zeta_r*out.rho_liq;
out.M_lub = V*(1-out.alpha)*out.zeta_oil*out.rho_liq;
out.M_ref = out.M_rl + out.M_rv;
out.M_liq = out.M_rl + out.M_lub;
out.M_vap = out.M_rv;
out.M_mix = out.M_ref + out.M_lub;

end
