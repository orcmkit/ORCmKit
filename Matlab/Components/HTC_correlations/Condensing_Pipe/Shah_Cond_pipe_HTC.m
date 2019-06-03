function [h_cond, Nu, flag ]  = Shah_Cond_pipe_HTC(x, mu_l, k_l, Pr_l, p_star, G, Dh, disp_flag) % verified
% Condensation heat transfer correlation published by Shah in
% "A general correlation for heat transfer during film condensation inside pipes", International Journal of Heat and Mass Transfer, 1979

% RDickes - 20/07/2018
x = max(0.00001,min(x,0.99999));
Re_l = G*Dh/mu_l;
Nu_DB = 0.023*(Re_l^0.8)*(Pr_l^0.4); % Dittus and Boelter correlation
Nu = Nu_DB*( ((1-x)^0.8) + (((3.8)*(x^0.76)*((1-x)^0.04))/(p_star^0.38)) ); % Shah's correction

h_cond = Nu*k_l/Dh;

flag = 1;

end