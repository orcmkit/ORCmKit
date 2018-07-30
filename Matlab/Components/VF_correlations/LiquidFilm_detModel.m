function [res, delta, alpha, u_v, Re_v, Re_lf, f_s, tau, res_tau, delta_plus, f_i, tau_i, dpdz, Mdot_l_bis] = LiquidFilm_detModel(delta_by_R_guess, D, G_v, rho_v, mu_v, nu_v, G_l, mu_l, rho_l, Mdot_l)
delta_by_R = delta_by_R_guess;  % guess on film relative thickness
R = D/2;
delta = delta_by_R*R;           % liquid film thickness
alpha = ((D-2*delta)/D)^2;      % void fraction
u_v = G_v/(rho_v*alpha);        % average vapour velocity
Re_v = rho_v*u_v*D/mu_v;        % Reynolds number of vapour phase
Re_lf = G_l*D/(4*mu_l);         % Reynolds number of liquid film flow
f_s = 0.046*Re_v^(-0.2) ;       % Smooth pipe friction factor

f = @(x) res_interfacial_shear_stress(x, delta, nu_v, rho_v, f_s, Re_v, Re_lf, u_v);
tau_0 = 2000;
opts = optimoptions(@fsolve,'Display','none');
tau = fsolve(f, tau_0,opts);    % interfacial shear stress
[res_tau, delta_plus, f_i, tau_i] = interfacial_shear_stress(tau, delta, nu_v, rho_v, f_s, Re_v, Re_lf, u_v);
dpdz = -(rho_v*9.81+(4*tau_i)/(D*sqrt(alpha))); % Pressure gradient
Mdot_l_bis = ((2*pi*rho_l/mu_l)*((tau_i*(R-delta)+((R-delta)^(2)/2)*(dpdz+rho_l*9.81))*((R^2-(R-delta)^2)/4 -((R-delta)^(2)/2)*log(1/(1-delta_by_R))))) -((pi*rho_l/(8*mu_l))*(dpdz+rho_l*9.81)*(R^2-(R-delta)^2)^2);

res = (Mdot_l - Mdot_l_bis)/Mdot_l;

end
