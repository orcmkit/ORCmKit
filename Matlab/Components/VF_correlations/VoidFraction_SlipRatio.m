function alpha = VoidFraction_SlipRatio(q, rho_v, rho_l, S)
alpha = 1./(1+(((1-q)./q).*(rho_v/rho_l)*S));
end