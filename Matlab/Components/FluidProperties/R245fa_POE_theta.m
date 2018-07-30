function theta = R245fa_POE_theta(zeta_r, P_Pa)
theta = (1-zeta_r)*(((-5.4713391e-02) + (5.6223340e-02)/zeta_r^0.5) + ((2.3650193e-02) + (-4.4703103e-02)/zeta_r^0.5 + (2.9199003e-02)/zeta_r + (-7.2069670e-03)/zeta_r^1.5 + (6.3295120e-04)/(zeta_r^2) )*(P_Pa/1e5));
end
