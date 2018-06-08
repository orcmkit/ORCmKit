function [ zeta_r, res_zeta_r] = R245fa_POE_solubility(T_K, P_Pa)
% RDickes - 11/04/2018
% Impletement correlations from Grebner, tuned with Zhelny's data 

T_sat = CoolProp.PropsSI('T', 'P', P_Pa, 'Q', 0, 'R245fa');
theta_1 = (T_K-T_sat)/T_sat;
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

f_theta = @(zeta_r, Pbar) (1-zeta_r)*(((-5.4713391e-02) + (5.6223340e-02)/zeta_r^0.5) + ((2.3650193e-02) + (-4.4703103e-02)/zeta_r^0.5 + (2.9199003e-02)/zeta_r + (-7.2069670e-03)/zeta_r^1.5 + (6.3295120e-04)/(zeta_r^2) )*Pbar);
theta_2 = f_theta(x, P_Pa/1e5);
res = theta_1 - theta_2;
end