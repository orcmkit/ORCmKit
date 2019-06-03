function alpha = VoidFraction_GrahamCondensation(q_vec, rho_v, G, D)
alpha = NaN*ones(size(q_vec));
ii= 0;
for q = q_vec
    ii = ii+1;
    Ft = ((q^3*G^2)/(9.81*rho_v^2*D*(1-q)))^0.5;
    if Ft>0.01032
        alpha(ii) = min(1,max(0,1-exp(-1-0.3*log(Ft)-0.0328*(log(Ft))^2)));
    else
        alpha(ii) = 0;
    end
end
end
