function x_di = Kim_DryOutIncipience(G, q, Dh, P_star, rho_l, rho_v, mu_l, sigma, i_fg)
Bo = q/G/i_fg;
We_lo = (Dh*G^2)/(rho_l*sigma);
Ca = (mu_l*G)/(rho_l*sigma);
x_di = 1.4*(We_lo^0.03)*(P_star^0.08) - 15*(Bo^0.15)*(Ca^0.35)*(rho_v/rho_l)^0.06;
end