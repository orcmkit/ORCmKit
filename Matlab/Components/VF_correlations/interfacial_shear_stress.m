function [res, delta_plus, f_i, tau_i] = interfacial_shear_stress(tau_guess, delta, nu_v, rho_v, f_s, Re_v, Re_lf, u_v)
delta_plus = delta/(nu_v)*(tau_guess/rho_v)^(0.5) ;
f_i = f_s*(1+0.0784*(delta_plus)^(1.4)*Re_v^(-0.3)*Re_lf^(-0.3)); %interfacial shear factor
tau_i = 0.5*f_i*rho_v*u_v^2; %interfacial shear stress
res = (tau_i-tau_guess);
end
