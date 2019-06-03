function [M_ref, M_lub] = Mass_function_cycle_solub(in1, P_Pa , V, C_oil, fluid_r, fluid_lub, param)
if C_oil >0
    T_K = in1;
    fit_ratio_rho = param.fit_ratio_rho;
    Tsat_pure = CoolProp.PropsSI('T', 'P', P_Pa, 'Q', 0, fluid_r);
    Tbubble_min = R245fa_POE_Tbubble(1-C_oil, P_Pa, Tsat_pure);
    [~, ~, ~, ~, zeta_r, ~, ~, C_rv] = TP_solubMixt(T_K, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min, Tsat_pure, param.fit_DTP_zeta);
    zeta_oil = 1-zeta_r;
    [rho_vap, rho_rl, rho_oil , rho_liq, ~] = R245fa_POE_density(T_K, P_Pa, zeta_r, fluid_r, fluid_lub, fit_ratio_rho, Tbubble_min, Tsat_pure);
    
    if strcmp(param.type_void_fraction, 'LM') || strcmp(param.type_void_fraction, 'Premoli') || strcmp(param.type_void_fraction, 'Hughmark') || strcmp(param.type_void_fraction, 'AnnularFlow')
        [mu_vap, ~, ~, mu_liq] = R245fa_POE_viscosity(T_K, P_Pa, C_oil, zeta_r, rho_liq, rho_oil, rho_rl, fluid_r, fluid_lub, Tbubble_min, Tsat_pure);
        param.mu_v = mu_vap;
        param.mu_l = mu_liq;
        param.sig = CoolProp.PropsSI('I', 'P', P_Pa, 'Q', 0, fluid_r);
    end
    
    if C_rv > 0 %there is a vapour phase
        alpha = VoidFraction_RefLub(C_rv, rho_vap, rho_liq, param);
    else % there is no vapour phase
        alpha = 0;
    end
    
    M_lub = V*(1-alpha)*zeta_oil*rho_liq;
    M_ref = V*(1-alpha)*zeta_r*rho_liq + V*alpha*rho_vap;
    
else
    h_Jkg = in1;
    h_liq = CoolProp.PropsSI('H', 'P', P_Pa, 'Q', 0, fluid_r);
    h_vap = CoolProp.PropsSI('H', 'P', P_Pa, 'Q', 1, fluid_r);
    if h_Jkg < h_liq % liquid phase only
        alpha = 0;
        rho_vap = 0;
        rho_liq = CoolProp.PropsSI('D', 'P', P_Pa, 'H', h_Jkg, fluid_r);
    elseif h_Jkg > h_vap % vapour phase only
        alpha = 1;
        rho_vap = CoolProp.PropsSI('D', 'P', P_Pa, 'H', h_Jkg, fluid_r);
        rho_liq = 0;
    else
        q = min(1,max(0,CoolProp.PropsSI('Q', 'P', P_Pa, 'H', h_Jkg, fluid_r)));
        rho_liq = CoolProp.PropsSI('D', 'P', P_Pa, 'Q', 0, fluid_r);
        rho_vap = CoolProp.PropsSI('D', 'P', P_Pa, 'Q', 1, fluid_r);
        if strcmp(param.type_void_fraction, 'LM') || strcmp(param.type_void_fraction, 'Premoli') || strcmp(param.type_void_fraction, 'Hughmark')
            param.sig = CoolProp.PropsSI('I', 'Q', 0, 'P', P_Pa, fluid_r);
            param.mu_v = CoolProp.PropsSI('V','P',P_Pa,'Q',1,fluid_r);
            param.mu_l = CoolProp.PropsSI('V','P',P_Pa,'Q',0,fluid_r);
        end
        alpha = VoidFraction_RefLub(q, rho_vap, rho_liq, param);
    end
    M_lub = 0;
    M_ref = V*(1-alpha)*rho_liq + V*alpha*rho_vap;
end
end
