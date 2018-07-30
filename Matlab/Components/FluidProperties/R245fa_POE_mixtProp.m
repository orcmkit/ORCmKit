function [Tbubble_min_K, Tsat_pure_K, h_mix, h_rl, h_rv, h_oil, C_rl, C_rv, C_oil, x, zeta_r, res_zeta_r, zeta_oil,  rho_vap, rho_rl, rho_oil, rho_liq, K, mu_vap, mu_rl, mu_oil, mu_liq] = R245fa_POE_mixtProp(T_K, P_Pa, C_oil, fluid_r, fluid_lub, fit_ratio_rho)

zeta_r_max = 1-C_oil;
Tsat_pure_K = CoolProp.PropsSI('T', 'P', P_Pa, 'Q', 0, fluid_r);
Tbubble_min_K = R245fa_POE_Tbubble(zeta_r_max, P_Pa, Tsat_pure_K);

if T_K >= Tbubble_min_K % The temperature is sufficiently high to permit a vapour phase
    
    [zeta_r, res_zeta_r] = R245fa_POE_solubility(T_K , P_Pa, Tsat_pure_K);
    zeta_oil = 1 - zeta_r;
    x = (1-zeta_r-C_oil)/(1-zeta_r-C_oil+zeta_r*C_oil);
    
    
    C_rl = (1-x)*(1-C_oil);
    C_rv = x*(1-C_oil);
    h_rl = CoolProp.PropsSI('H', 'T', T_K, 'Q', 0, fluid_r);
    if abs(T_K-Tsat_pure_K)<1e-3
        h_rv = CoolProp.PropsSI('H', 'T', T_K, 'Q', 1, fluid_r);
        rho_vap = CoolProp.PropsSI('D', 'T', T_K, 'Q', 1, fluid_r);
        mu_vap = CoolProp.PropsSI('V', 'T', T_K, 'Q', 1, fluid_r);
    else
        h_rv = CoolProp.PropsSI('H', 'T', T_K, 'P', P_Pa, fluid_r);
        rho_vap = CoolProp.PropsSI('D', 'T', T_K, 'P', P_Pa, fluid_r);
        try
            mu_vap = CoolProp.PropsSI('V', 'T', T_K, 'P', P_Pa, fluid_r);
        catch
            mu_vap = NaN;
        end
    end
    h_oil =  sf_PropsSI_bar('H', T_K, T_K, P_Pa, fluid_lub);
    h_mix = C_rl*h_rl + C_rv*h_rv + C_oil*h_oil;
    
    rho_rl= CoolProp.PropsSI('D', 'T', T_K, 'Q', 0, fluid_r);
    rho_oil = sf_PropsSI_bar('D', T_K, T_K, P_Pa, fluid_lub);
    rho_ideal = rho_oil/(1+zeta_r*(rho_oil/rho_rl-1));
    K = fit_ratio_rho(T_K, zeta_r);
    rho_liq = rho_ideal/K;
        
    mu_rl= CoolProp.PropsSI('V', 'T', T_K, 'Q', 0, fluid_r);
    mu_oil = sf_PropsSI_bar('V', T_K, T_K, P_Pa, fluid_lub);  
    mu_liq = rho_liq*exp(zeta_oil*log(mu_oil/rho_oil) + zeta_r*log(mu_rl/rho_rl)); % Schroeder equation (from M Conde)

    
else  % The temperature is too low and there is only a liquid phase
    x = 0;
    zeta_oil = C_oil;
    zeta_r = 1-C_oil;
    res_zeta_r = 0;
    C_rl = zeta_r;
    C_rv = 0;
    
    h_rv = 0;
    if T_K >= Tsat_pure_K-1e-2
        h_rl = CoolProp.PropsSI('H', 'T', T_K, 'Q', 0, fluid_r);
        rho_rl= CoolProp.PropsSI('D', 'T', T_K, 'Q', 0, fluid_r);
        mu_rl= CoolProp.PropsSI('V', 'T', T_K, 'Q', 0, fluid_r);
    else
        h_rl = CoolProp.PropsSI('H', 'T', T_K, 'P', P_Pa, fluid_r);
        rho_rl = CoolProp.PropsSI('D', 'T', T_K, 'P', P_Pa, fluid_r);
        mu_rl = CoolProp.PropsSI('V', 'T', T_K, 'P', P_Pa, fluid_r);

    end
    h_oil =  sf_PropsSI_bar('H', T_K, T_K, P_Pa, fluid_lub);
    h_mix = C_rl*h_rl + C_oil*h_oil;
    
    
    rho_vap = 0;
    rho_oil = sf_PropsSI_bar('D', T_K, T_K, P_Pa, fluid_lub);    
    rho_ideal = rho_oil/(1+zeta_r*(rho_oil/rho_rl-1));
    K = fit_ratio_rho(T_K, zeta_r);
    rho_liq = rho_ideal/K;
    
    mu_vap = 0;
    mu_oil = sf_PropsSI_bar('V', T_K, T_K, P_Pa, fluid_lub);  
    mu_liq = rho_liq*exp(zeta_oil*log(mu_oil/rho_oil) + zeta_r*log(mu_rl/rho_rl)); % Schroeder equaition (from M Conde)
end

end