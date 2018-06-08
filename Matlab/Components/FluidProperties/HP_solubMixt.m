function [T, Tbubble_min, Tsat_pure, h_rl, h_oil, h_rv, zeta_r, x, C_rl, C_rv] = HP_solubMixt(C_oil, P_Pa, h_Jkg, fluid_r, fluid_lub)
% Function that finds properties of a oil-refrigerant mixture based
% on the enthalpy and pressure
zeta_r_max = 1-C_oil;
Tbubble_min = R245fa_POE_Tbubble(zeta_r_max, P_Pa);
Tsat_pure = CoolProp.PropsSI('T', 'P', P_Pa, 'Q', 0, fluid_r);
f = @(x) res_find_T_solubMixt(x, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min, Tsat_pure, h_Jkg);
Tlb = 0+273.15;
Tub = 153+273.15;
[ T, res_T] = fzero(f,[Tlb Tub]);
if abs(res_T) > 1e-2
    display(['Problem in finding T for C_oil = ' num2str(100*C_oil) ' h = ' num2str(h) ' J/kg - P = ' num2str(P_Pa) ' Pa'])
end
[~, h_rl, h_oil, h_rv, zeta_r, x, C_rl, C_rv] = TP_solubMixt(T, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min, Tsat_pure);
end

function res = res_find_T_solubMixt(T_K, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min, Tsat_pure, h_Jkg)
[h_mix, ~, ~, ~, ~, ~, ~, ~] = TP_solubMixt(T_K, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min, Tsat_pure);
res = (h_Jkg-h_mix)/h_Jkg;
end

