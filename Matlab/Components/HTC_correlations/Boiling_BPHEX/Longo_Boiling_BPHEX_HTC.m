function h_boiling = Longo_Boiling_BPHEX_HTC(x, k_l, rho_l, Pr_l, mu_l, rho_v, G, Dh, phi, h_conv_h, DT_log, Qdot, P_su, fluid, info)
% Boiling heat transfer correlation for BPHEX published by Longo et. al in
% "A new model for refrigerant boiling inside Brazed Plate Heat Exchangers (BPHEs)"


G_eq = G * ( (1 - x) + x * (rho_l/rho_v)^0.5);
Re_eq = G_eq*Dh/mu_l;
AU_tp = Qdot/DT_log;

% Convective boiling HTC
h_conv_boiling = 0.122*phi*(k_l/Dh)*Re_eq^0.8*Pr_l^0.333333333333;

% Nucleat boiling HTC
P_crit = CoolProp.PropsSI('Pcrit', 'Q', 1, 'P', P_su, fluid);
p_star = P_su/P_crit;
sigma = CoolProp.PropsSI('I', 'Q', 0, 'P', 0.1*P_crit, fluid)*1000; %mN/m
P1 = (0.1*P_crit)+1e3; T1 = CoolProp.PropsSI('T', 'Q', 1, 'P', P1, fluid);
P2 = (0.1*P_crit)-1e3; T2 = CoolProp.PropsSI('T', 'Q', 1, 'P', P2, fluid);
dp_dTsat = (P1-P2)/(T1-T2)/1000;
P_f = dp_dTsat/sigma;
h0 = (3.58*P_f^0.6)*1000; % W/m².K --> cfr VDI page 765
C_ra = (0.4/0.4)^0.13333333;
F_p = (1.2*p_star^0.27) + (2.5 + 1/(1-p_star))*p_star;
q_0 = 20e3;
ub = 20e3;
f_q = @(xx) res_iter_Longo_nucleate_boiling(xx, h0, C_ra, F_p, h_conv_h, AU_tp, Qdot, phi, ub);
opts = optimoptions('fsolve', 'display', 'none');
[q,err_q] = fsolve(f_q, q_0./ub, opts);
q_nucleate_boiling = q*ub;
[res_q_nucleate_boiling, h_nucleate_boiling, U_nucleate_boiling, A_tp_nucleate_boiling, q_bis_nucleate_boiling]= iter_Longo_nucleate_boiling(q_nucleate_boiling, h0, C_ra, F_p, h_conv_h, AU_tp, Qdot, info);
if err_q > 1e-5
    display(['Longo boiling: Wrong heat flux --> err_q = ' num2str(err_q) ' !!!'])
end

% Final boiling HTC
h_boiling = max(h_conv_boiling, h_nucleate_boiling);


end

function res_q = res_iter_Longo_nucleate_boiling(x, h0, C_ra, F_p, h_conv_h, AU_tp, Qdot, phi, ub)
q = x*ub;
[res_q,~, ~, ~,~]= iter_Longo_nucleate_boiling(q, h0, C_ra, F_p, h_conv_h, AU_tp, Qdot, phi);
end

function [res_q, h, U, A_tp,q_bis]= iter_Longo_nucleate_boiling(q, h0, C_ra, F_p, h_conv_h, AU_tp, Qdot, phi)
h = 0.58*phi*h0*C_ra*F_p*(q/(20e3))^0.467;
U = (1/h +  1/h_conv_h)^-1;
A_tp = AU_tp/U;
q_bis = Qdot/A_tp;
res_q = (q-q_bis);
end
