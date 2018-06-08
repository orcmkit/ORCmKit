function Tbubble_K= R245fa_POE_Tbubble(zeta_r, P_Pa)
% RDickes - 11/04/2018
% Impletement correlations from Grebner, tuned with Zhelny's data 
Pbar = P_Pa/1e5;
theta = (1-zeta_r)*(((-5.4713391e-02) + (5.6223340e-02)/zeta_r^0.5) + ((2.3650193e-02) + (-4.4703103e-02)/zeta_r^0.5 + (2.9199003e-02)/zeta_r + (-7.2069670e-03)/zeta_r^1.5 + (6.3295120e-04)/(zeta_r^2) )*Pbar);
T_sat = CoolProp.PropsSI('T', 'P', P_Pa, 'Q', 0, 'R245fa');
Tbubble_K = (1+theta)*T_sat;
end
