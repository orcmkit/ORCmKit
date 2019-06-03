function [M_rv, M_rl, M_lub, M_ref, M_liq, M_vap, M_mix, alpha_mean] = MassMeanInt_solub_Fast(T_K_1, P_Pa_1, zeta_r_1, C_rv_1, Tbubble_min_1, Tsat_pure_1, T_K_2, P_Pa_2, zeta_r_2, C_rv_2, Tbubble_min_2, Tsat_pure_2, V, fluid_r, fluid_lub, fit_ratio_rho, param)
%compute  flow densities
[rho_vap_1, rho_rl_1, rho_oil_1, rho_liq_1, ~] = R245fa_POE_density_Fast(T_K_1, P_Pa_1, zeta_r_1, fluid_r, fluid_lub, fit_ratio_rho, Tbubble_min_1, Tsat_pure_1, param);
[rho_vap_2, rho_rl_2, rho_oil_2, rho_liq_2, ~] = R245fa_POE_density_Fast(T_K_2, P_Pa_2, zeta_r_2, fluid_r, fluid_lub, fit_ratio_rho, Tbubble_min_2, Tsat_pure_2, param);
rho_vap_mean = 0.5*rho_vap_1 + 0.5*rho_vap_2;
rho_liq_mean = 0.5*rho_liq_1 + 0.5*rho_liq_2;
zeta_r_mean = 0.5*zeta_r_1+0.5*zeta_r_2;
if 0.5*C_rv_1 + 0.5*C_rv_2 > 0 %there is a vapour phase
    
    if strcmp(param.correlation.type_void_fraction, 'LM') || strcmp(param.correlation.type_void_fraction, 'Premoli') || strcmp(param.correlation.type_void_fraction, 'Hughmark')
        [mu_vap_1, ~, ~, mu_liq_1] = R245fa_POE_viscosity_Fast(T_K_1, P_Pa_1, param.C_oil, zeta_r_1, rho_liq_1, rho_oil_1, rho_rl_1, fluid_r, fluid_lub, Tbubble_min_1, Tsat_pure_1, param);
        [mu_vap_2, ~, ~, mu_liq_2] = R245fa_POE_viscosity_Fast(T_K_2, P_Pa_2, param.C_oil, zeta_r_2, rho_liq_2, rho_oil_2, rho_rl_2, fluid_r, fluid_lub, Tbubble_min_2, Tsat_pure_2, param);
        param.mu_v = 0.5*mu_vap_1 + 0.5*mu_vap_2;
        param.mu_l = 0.5*mu_liq_1 + 0.5*mu_liq_2;
        param.sig = CoolPropt_I_PQ(param.abs_lowLevel,param.CP_file, 0.5*P_Pa_1+0.5*P_Pa_2, 0); 
    end
    alpha_mean = VoidFraction_Integration(min(1,max(0,C_rv_1)), min(1,max(0,C_rv_2)), rho_vap_mean, rho_liq_mean, param);
else % there is no vapour phase
    alpha_mean = 0;
end

%alpha_mean = alpha_mean*param.C_fit_VFM;

M_rv = V*alpha_mean*rho_vap_mean;
M_rl = V*(1-alpha_mean)*zeta_r_mean*rho_liq_mean;
M_lub = V*(1-alpha_mean)*(1-zeta_r_mean)*rho_liq_mean;
M_ref = M_rl + M_rv;
M_liq = M_rl + M_lub;
M_vap = M_rv;
M_mix = M_ref + M_lub;

end
MassMeanInt_solub