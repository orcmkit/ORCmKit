function h_cond = Han_Cond_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v, G, Dh, pitch_co, theta)
% Condensation heat transfer correlation published by Han et. al in
% "The caracteristics of condensation in brazed plate heat exchangers with 
% different chevron angles", Journal of the Korean Physical Society, 2003

% RDickes - 26/04/2018

G_eq = G * ( (1-x) + x * (rho_l/rho_v)^0.5);
Re_eq = G_eq*Dh/mu_l;
Ge1 = 11.22*(pitch_co/Dh)^(-2.83)*(theta)^(-4.5);
Ge2 = 0.35*(pitch_co/Dh)^(0.23)*(theta)^(1.48);
Nu = Ge1*Re_eq^Ge2*Pr_l^0.33333333;
h_cond = Nu*k_l/Dh;
end