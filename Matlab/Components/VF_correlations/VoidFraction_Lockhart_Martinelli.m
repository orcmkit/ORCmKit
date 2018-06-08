function alpha = VoidFraction_Lockhart_Martinelli(q, rho_v, rho_l, mu_v, mu_l)
X_tt = (((1-q)/q)^0.9)*sqrt(((mu_l/mu_v)^0.1)*(rho_v/rho_l));
if X_tt <= 10
    alpha = (1+X_tt^0.8)^-0.378;
else
    alpha = 0.823-0.157*log(X_tt);
end
end