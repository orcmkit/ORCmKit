function [rho_vap, rho_rl, rho_oil, rho_liq, K] = R245fa_POE_density(T_K, P_Pa, zeta_r, fluid_r, fluid_lub, fit_ratio_rho, Tbubble_min_K, Tsat_pure_K)

if T_K >= Tbubble_min_K % The temperature is sufficiently high to permit a vapour phase
    %[zeta_r, ~] = R245fa_POE_solubility(T_K , P_Pa, Tsat_pure_K);
    if abs(T_K-Tsat_pure_K)<1e-3
        rho_vap = CoolProp.PropsSI('D', 'T', T_K, 'Q', 1, fluid_r);
    else
        rho_vap = CoolProp.PropsSI('D', 'T', T_K, 'P', P_Pa, fluid_r);
    end
    rho_rl= CoolProp.PropsSI('D', 'T', T_K, 'Q', 0, fluid_r);
    rho_oil = PropsSI_ICP('D', 'T',T_K, 'P', P_Pa, fluid_lub);
    rho_ideal = rho_oil/(1+zeta_r*(rho_oil/rho_rl-1));
    K = fit_ratio_rho(max(285,min(T_K,422)), max(1e-4, min(0.9999,zeta_r)));
    rho_liq = rho_ideal/K;
    
else  % The temperature is too low and there is only a liquid phase
    if T_K >= Tsat_pure_K-1e-2
        rho_rl= CoolProp.PropsSI('D', 'T', T_K, 'Q', 0, fluid_r);
    else
        rho_rl = CoolProp.PropsSI('D', 'T', T_K, 'P', P_Pa, fluid_r);
    end
    rho_vap = CoolProp.PropsSI('D', 'T', T_K, 'Q', 1, fluid_r); %0;
    rho_oil = PropsSI_ICP('D', 'T',T_K, 'P', P_Pa, fluid_lub);
    rho_ideal = rho_oil/(1+zeta_r*(rho_oil/rho_rl-1));
    K = fit_ratio_rho(max(285,min(T_K,422)), max(1e-4, min(0.9999,zeta_r)));
    rho_liq = rho_ideal/K;
end

end