function [mu_vap, mu_rl, mu_oil, mu_liq] = R245fa_POE_viscosity(T_K, P_Pa, C_oil, zeta_r, rho_liq, rho_oil, rho_rl, fluid_r, fluid_lub, Tbubble_min_K, Tsat_pure_K)

if T_K >= Tbubble_min_K % The temperature is sufficiently high to permit a vapour phase
    zeta_oil = 1 - zeta_r;
    if abs(T_K-Tsat_pure_K)<1e-3
        mu_vap = CoolProp.PropsSI('V', 'T', T_K, 'Q', 1, fluid_r);
    else
        try
            mu_vap = CoolProp.PropsSI('V', 'T', T_K, 'P', P_Pa, fluid_r);
        catch
            mu_vap = CoolProp.PropsSI('V', 'T', Tsat_pure_K+5, 'P', P_Pa, fluid_r);
        end
    end
    mu_rl= CoolProp.PropsSI('V', 'T', T_K, 'Q', 0, fluid_r);
    mu_oil = PropsSI_ICP('V', 'T',T_K, 'P', P_Pa, fluid_lub);
    mu_liq = rho_liq*exp(zeta_oil*log(mu_oil/rho_oil) + zeta_r*log(mu_rl/rho_rl)); % Schroeder equation (from M Conde)
    
else  % The temperature is too low and there is only a liquid phase
    zeta_oil = C_oil;
    if T_K >= Tsat_pure_K-1e-2
        mu_rl= CoolProp.PropsSI('D', 'T', T_K, 'Q', 0, fluid_r);
    else
        mu_rl = CoolProp.PropsSI('D', 'T', T_K, 'P', P_Pa, fluid_r);
    end
    mu_vap = 0;
    mu_oil = PropsSI_ICP('V', 'T',T_K, 'P', P_Pa, fluid_lub);
    mu_liq = rho_liq*exp(zeta_oil*log(mu_oil/rho_oil) + zeta_r*log(mu_rl/rho_rl)); % Schroeder equaition (from M Conde)
    
end

end