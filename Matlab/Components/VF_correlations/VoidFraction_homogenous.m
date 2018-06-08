function alpha = VoidFraction_homogenous(q, rho_v, rho_l)
alpha = 1./(1+((1-q)./q).*(rho_v/rho_l));
end