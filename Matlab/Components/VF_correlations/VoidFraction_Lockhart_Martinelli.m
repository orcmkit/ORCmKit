function alpha = VoidFraction_Lockhart_Martinelli(q, rho_v, rho_l, mu_v, mu_l)
X_tt = NaN*ones(size(q));
alpha =  NaN*ones(size(q));
for ii =1:length(q)
    X_tt(ii) = (((1-q(ii))/q(ii))^0.9)*((mu_l/mu_v)^0.1)*((rho_v/rho_l)^0.5);
    if X_tt(ii) <= 10
        alpha(ii) = (1+X_tt(ii)^0.8)^-0.378;
    else
        alpha(ii) = 0.823-0.157*log(X_tt(ii));
    end
end
end