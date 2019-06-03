function [T, h_mix, h_rl, h_oil, h_rv, zeta_r, C_rl, C_rv] = xP_solubMixt(C_oil, P_Pa, x_q, fluid_r, fluid_lub, Tbubble_min_K, Tsat_pure_K, fit_DTP_zeta)
% Function that finds properties of a oil-refrigerant mixture based
% on the refrigerant quality and pressure

zeta_r = (1- x_q - C_oil + x_q*C_oil)/(1- x_q + x_q*C_oil);
T = R245fa_POE_Tbubble(zeta_r, P_Pa, Tsat_pure_K);

[h_mix, h_rl, h_oil, h_rv, zeta_r, ~, C_rl, C_rv] = TP_solubMixt(T, P_Pa, C_oil, fluid_r, fluid_lub, Tbubble_min_K, Tsat_pure_K, fit_DTP_zeta);
           
end
