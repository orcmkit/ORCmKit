function alpha = VoidFraction_Zivi(q, rho_v, rho_l)
S_zivi = (rho_v/rho_l)^(-1/3);
alpha = 1./(1+(((1-q)./q).*(rho_v/rho_l)*S_zivi));
end