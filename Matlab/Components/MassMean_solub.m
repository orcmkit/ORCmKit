function [M_rv, M_rl, M_lub, M_ref, M_liq, M_vap, M_mix, alpha_mean] = MassMean_solub(T_K_mean, P_Pa_mean , V, C_oil, fluid_r, fluid_lub, Tbubble_min, Tsat_pure, fit_ratio_rho, param)

%compute mean flow composition
[~, ~, ~, ~, zeta_r_mean, ~, ~, C_rv_mean] = TP_solubMixt(T_K_mean, P_Pa_mean, C_oil, fluid_r, fluid_lub, Tbubble_min, Tsat_pure); 
   
%compute mean flow densities
[rho_vap_mean, ~, ~, rho_liq_mean, ~] = R245fa_POE_density(T_K_mean, P_Pa_mean, zeta_r_mean, fluid_r, fluid_lub, fit_ratio_rho, Tbubble_min, Tsat_pure);

if C_rv_mean > 0 %there is a vapour phase
    switch param.correlation.type_void_fraction
        case 'Homogenous'
            alpha_mean = VoidFraction_homogenous(C_rv_mean, rho_vap_mean,  rho_liq_mean);
            
        case 'Zivi'
            alpha_mean = VoidFraction_Zivi(C_rv_mean, rho_vap_mean,  rho_liq_mean);
            
        case 'SlipRatio'
            alpha_mean = VoidFraction_SlipRatio(C_rv_mean, rho_vap_mean,  rho_liq_mean, param.S_ratio);
            
        case 'LockMart'
            alpha_mean = VoidFraction_Lockhart_Martinelli(C_rv_mean, rho_vap_mean,  rho_liq_mean, mu_vap_mean, mu_liq_mean);
    end
else % there is no vapour phase
    alpha_mean = 0;
end

M_rv = V*alpha_mean*rho_vap_mean;
M_rl = V*(1-alpha_mean)*zeta_r_mean*rho_liq_mean;
M_lub = V*(1-alpha_mean)*(1-zeta_r_mean)*rho_liq_mean;
M_ref = M_rl + M_rv;
M_liq = M_rl + M_lub;
M_vap = M_rv;
M_mix = M_ref + M_lub;

end
