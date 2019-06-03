function [ zeta_r, res_zeta_r] = R245fa_POE_solubility_1(T_K, P_Pa, Tsat_pure_K)
% RDickes - 11/04/2018
% Impletement correlations from Grebner, tuned with Zhelny's data 

theta_1 = (T_K-Tsat_pure_K)/Tsat_pure_K;
if theta_1 >= 0
    f = @(x)res_R245fa_POE_solubility(x, theta_1, P_Pa);
    [ zeta_r, res_zeta_r] = fzero(f,[1e-10 0.999999999999999999999999]);
else
    zeta_r = 0;
    res_zeta_r = -2;
    disp('Problem zeta_r calculation');
end
end

function res = res_R245fa_POE_solubility(x, theta_1, P_Pa)
theta_2 = R245fa_POE_theta(x, P_Pa);
res = theta_1 - theta_2;
end