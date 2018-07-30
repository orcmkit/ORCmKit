function [h_mix, h_rl, h_oil, h_rv, zeta_r, x, C_rl, C_rv] = TP_solubMixt(T_K, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min_K, Tsat_pure_K, fit_DTP_zeta)
% Function that finds properties of a oil-refrigerant mixture based
% on the temperature and pressure
if T_K > Tbubble_min_K % liquid + two-phase
    
    %[zeta_r, res_zeta_r] = R245fa_POE_solubility(T_K , P_Pa, Tsat_pure_K);
    [zeta_r, res_zeta_r] = R245fa_POE_solubility_tab(T_K , P_Pa, Tsat_pure_K, fit_DTP_zeta);

    if abs(res_zeta_r) < 1e-4
        x = (1-zeta_r-C_oil)/(1-zeta_r-C_oil+zeta_r*C_oil);
        C_rl = (1-x)*(1-C_oil);
        C_rv = x*(1-C_oil);
        h_rl = CoolProp.PropsSI('H', 'T', T_K, 'Q', 0, fluid_r);
        if abs(T_K-Tsat_pure_K)<1e-3
            h_rv = CoolProp.PropsSI('H', 'T', T_K, 'Q', 1, fluid_r);
        else
            h_rv = CoolProp.PropsSI('H', 'T', T_K, 'P', P_Pa, fluid_r);
        end
        h_oil =  PropsSI_ICP('H', 'T', T_K, 'P', P_Pa, fluid_lub);
        h_mix = C_rl*h_rl + C_rv*h_rv + C_oil*h_oil;
    else
        h_mix = NaN;
        h_rl  = NaN;
        h_oil  = NaN;
        h_rv  = NaN; 
        zeta_r  = NaN; 
        x  = NaN; 
        C_rl  = NaN;
        C_rv  = NaN;
    end
    
else % liquid-only
    
    if T_K >= Tsat_pure_K-1e-2
        h_rl = CoolProp.PropsSI('H', 'T', T_K, 'Q', 0, fluid_r);
    else
        h_rl = CoolProp.PropsSI('H', 'T', T_K, 'P', P_Pa, fluid_r);
    end
    h_rv = 0;
    zeta_r = 1-C_oil;
    C_rl = 1-C_oil;
    C_rv = 0;
    x = 0;
    h_oil =  PropsSI_ICP('H', 'T', T_K, 'P', P_Pa, fluid_lub);
    h_mix = (1-C_oil)*h_rl + C_oil*h_oil;    
end

end