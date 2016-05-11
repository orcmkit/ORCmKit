function res = FCT_Exp2_SemiEmp_res(x, ub, fluid, P_su, h_su, s_su, rho_su, h_ex_s, M_dot, V_s, r_v_in, P_ex, A_leak0, d_su, alpha, W_dot_loss_0, AU_su_n, M_dot_n, AU_ex_n, AU_amb, T_amb, C_loss, param)

% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems
%
% Remi Dickes - 11/05/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% FCT_Exp2_SemiEmp_res is a constitutive subfunction of ExpanderModel_SemiEmp.m
% It isolates the residuals of the semi-empirical model assuming a guess wall
% temperature.

%%

T_w = x.*ub;
out = FCT_Exp2_SemiEmp(T_w, fluid, P_su, h_su, s_su, rho_su, h_ex_s, M_dot, V_s, r_v_in, P_ex, A_leak0, d_su, alpha, W_dot_loss_0, AU_su_n, M_dot_n, AU_ex_n, AU_amb, T_amb, C_loss, param);
res = out.res;

end