function [ zeta, res_zeta] = R245fa_POE_solubility2(T_K, P_Pa,fluid)
% RDickes - 11/04/2018
% Impletement correlations Thome's to retrieve solubility

f = @(x)res_R245fa_POE_Tbubble(x, T_K, P_Pa, fluid);
[ zeta_oil, res_zeta] = fzero(f,[1e-10 0.999999999999999999999999]);
zeta = 1-zeta_oil;
end

function res = res_R245fa_POE_Tbubble(x, T_K, P_Pa, fluid)
T_bubble = Thome_Tbubble(x, P_Pa, fluid);
res = T_K - T_bubble;
end