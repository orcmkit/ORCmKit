function alpha = VoidFraction_Premoli(q, rho_v, rho_l, mu_l, sig, D, G)
Re_f = G*D/mu_l;
We_f = G^2*D/sig/rho_l;
alpha_h = VoidFraction_homogenous(q, rho_v, rho_l);
y = alpha_h./(1-alpha_h);
F1 = 1.578*(Re_f^-0.19)*(rho_l/rho_v)^0.22;
F2 = 0.0273*We_f*(Re_f^-0.51)*(rho_l/rho_v)^-0.08;
S_Premoli = 1+F1*sqrt(max(0,(y./(1+y.*F2)) -y.*F2));
alpha = 1./(1+(((1-q)./q).*(rho_v/rho_l).*S_Premoli));
end
