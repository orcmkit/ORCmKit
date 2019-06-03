function [h_cond, Nu, flag ]  = Cavallini_Cond_pipe_HTC(x, mu_l, mu_v, rho_l, rho_v, k_l, Pr_l, i_fg, DT_wall, G, Dh, disp_flag) % verified
% Condensation heat transfer correlation published by Cavallini in
% "Condensation in Horizontal Smooth Tubes: A New Heat Transfer Model for
% Heat Exchanger Design", Heat Transfer Engineering, 2006

% RDickes - 25/07/2018
    x = min(0.999999,max(x, 1e-6));
    g = 9.81;
    C_T = 2.6; %1.6 for HC or 2.6 for other refrigerant
    X_tt = ((mu_l/mu_v)^0.1)*((rho_v/rho_l)^0.5)*((1-x)/x)^0.9; %LM factor
    J_v_T = (((7.5/(4.3*X_tt^1.111 + 1))^-3) + ((C_T)^-3) )^-0.333333333333333333;
    J_v = x*G/sqrt(g*Dh*rho_v*(rho_l-rho_v));
    Re_l = G*Dh/mu_l;
    
    h_lo = 0.023*Re_l^(0.8)*Pr_l^0.4*k_l/Dh;
    K = 1 + (1.128*x^0.817)*((rho_l/rho_v)^0.3685)*((mu_l/mu_v)^0.2363)*((1-mu_v/mu_l)^2.144)*(Pr_l^-0.1);
    h_A = h_lo*K;
    h_strat = 0.725*((1+0.741*((1-x)/x)^0.3321)^-1)*(((k_l^3*rho_l*(rho_l-rho_v)*g*i_fg)/(mu_l*Dh*DT_wall))^0.25)+((1-x^0.087)*h_lo);

    if J_v > J_v_T %delta_T-independent flow regime
        h_cond = h_A;
    elseif J_v <= J_v_T %delta_T-dependent flow regime
        h_cond = (h_A*(J_v_T/J_v)^0.8 - h_strat)*(J_v/J_v_T) + h_strat;
    end
    flag = 1;
    Nu = h_cond*Dh/k_l;
    if isnan(Nu)
        format long
        disp(['[h_cond, Nu, flag ]  = Cavallini_Cond_pipe_HTC(' num2str(x) ',' num2str(mu_l) ',' num2str(mu_v) ',' num2str(rho_l) ',' num2str(rho_v) ',' num2str(k_l) ',' num2str(Pr_l) ',' num2str(i_fg) ',' num2str(DT_wall) ',' num2str(G) ',' num2str(Dh) ',' num2str(disp_flag) ')'])
        scheiss
    end
end

