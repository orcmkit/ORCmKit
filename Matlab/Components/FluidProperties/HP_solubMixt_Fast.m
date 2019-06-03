function [T, h_rl, h_oil, h_rv, zeta_r, x, C_rl, C_rv] = HP_solubMixt_Fast(C_oil, P_Pa, h_Jkg, fluid_r, fluid_lub, Tbubble_min_K, Tsat_pure_K, fit_DTP_zeta, param)
% Function that finds properties of a oil-refrigerant mixture based
% on the enthalpy and pressure
f = @(x) res_find_T_solubMixt_Fast(x, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min_K, Tsat_pure_K, h_Jkg, fit_DTP_zeta,param);
Tlb = 0+273.15;
Tub = 153+273.15;
%options = optimset('TolX',1e-8);
[ T, res_T, exitflag] = fzero(f,[Tlb Tub]);%,options);

if exitflag<0 || abs(res_T)>1e-2
    disp(['flag = ' num2str(exitflag) ' - res_T = ' num2str(abs(res_T)) ' --> Problem in finding HP_solubMixt(' num2str(C_oil) ',' num2str(P_Pa) ',' num2str(h_Jkg) ',''' fluid_r ''',''' fluid_lub ''',' num2str(Tbubble_min_K) ',' num2str(Tsat_pure_K) ', fit_DT_Pbar_zeta2)'])
end
[~, h_rl, h_oil, h_rv, zeta_r, x, C_rl, C_rv] = TP_solubMixt_Fast(T, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min_K, Tsat_pure_K, fit_DTP_zeta, param);
end

function res = res_find_T_solubMixt_Fast(T_K, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min_K, Tsat_pure_K, h_Jkg, fit_DTP_zeta, param)
[h_mix, ~, ~, ~, ~, ~, ~, ~] = TP_solubMixt_Fast(T_K, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min_K, Tsat_pure_K, fit_DTP_zeta, param);
res = (h_Jkg-h_mix)/h_Jkg;
end

