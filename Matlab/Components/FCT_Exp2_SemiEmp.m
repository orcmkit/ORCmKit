function out = FCT_Exp2_SemiEmp(T_w, fluid, P_su, h_su, s_su, rho_su, h_ex_s, M_dot, V_s, r_v_in, P_ex, A_leak0, d_su, alpha, W_dot_loss_0, AU_su_n, M_dot_n, AU_ex_n, AU_amb, T_amb, C_loss, param)

% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems
%
% Remi Dickes - 11/05/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% FCT_Exp2_SemiEmp is a constitutive subfunction of ExpanderModel_SemiEmp.m
% It evaluates the semi-empirical model considering the supply
% conditions and a guess value of the wall temperature.


AU_su1 = AU_su_n*(M_dot/M_dot_n)^0.8;
AU_ex1 = AU_ex_n*(M_dot/M_dot_n)^0.8;
out.h_su1 = h_su;
out.P_su1 = max(CoolProp.PropsSI('P','S',s_su,'H',max(h_su - (M_dot/(pi*d_su^2/4*rho_su))^2/2, h_ex_s), fluid),P_ex+1);
out.DP_su = P_su-out.P_su1;
out.T_su1 = CoolProp.PropsSI('T','P',out.P_su1,'H',out.h_su1,fluid);
out.Q_su1 = CoolProp.PropsSI('Q','P',out.P_su1,'H',out.h_su1,fluid);
if out.Q_su1 < 0
    out.cp_su1 = CoolProp.PropsSI('C','P',out.P_su1,'H',out.h_su1,fluid);
else
    out.cp_su1 = CoolProp.PropsSI('C','P',out.P_su1,'Q',0,fluid);
end
out.epsilon_su1 = max(0,(1-exp(-AU_su1/(M_dot*out.cp_su1))));
out.Q_dot_su = max(0,out.epsilon_su1*M_dot*out.cp_su1*(out.T_su1 - T_w));
out.h_su2 = min(param.h_max,max(h_ex_s,out.h_su1 - out.Q_dot_su/M_dot));
out.s_su2 = CoolProp.PropsSI('S','H',out.h_su2,'P',out.P_su1,fluid);
out.rho_su2 = CoolProp.PropsSI('D','H',out.h_su2,'P',out.P_su1,fluid);
out.Q_su2 = CoolProp.PropsSI('Q','P',out.P_su1,'H',out.h_su2,fluid);
if  out.Q_su2 < 0
    out.gamma  = feval(param.gamma.gamma_PT_pol,[out.P_su1/1e5, CoolProp.PropsSI('T','P',out.P_su1,'H',out.h_su2,fluid)/1e2]);
else
    out.gamma = feval(param.gamma.gamma_PQ_pol,[out.P_su1/1e5, out.Q_su2]);
end
P_crit = out.P_su1*(2/(out.gamma+1))^(out.gamma/(out.gamma-1));
out.resgamma = 0;
out.P_thr = max(P_ex,P_crit);
out.rho_thr = CoolProp.PropsSI('D','P',out.P_thr,'S',out.s_su2,fluid);
out.C_thr = sqrt(2*(out.h_su2 - CoolProp.PropsSI('H','P',out.P_thr,'S',out.s_su2,fluid)));
out.M_dot_leak = A_leak0*out.C_thr*out.rho_thr;
out.M_dot_in = M_dot-out.M_dot_leak;
out.N_exp = 60*out.M_dot_in/(V_s*out.rho_su2);
out.rho_in = out.rho_su2/r_v_in;
try
    out.P_in = CoolProp.PropsSI('P','D',out.rho_in,'S',out.s_su2,fluid);
catch
    delta = 0.001;
    out.P_in = 0.5*CoolProp.PropsSI('P','D',out.rho_in*(1+delta),'S',out.s_su2,fluid)+0.5*CoolProp.PropsSI('P','D',out.rho_in*(1-delta),'S',out.s_su2,fluid);
end
out.h_in = CoolProp.PropsSI('H','D',out.rho_in,'P',out.P_in,fluid);
out.w_1 = out.h_su2-out.h_in;
out.w_2 = (out.P_in - P_ex)/out.rho_in;
out.h_ex2 = out.h_in - out.w_2;
out.W_dot_in = out.M_dot_in*(out.w_1+out.w_2);
out.W_dot_loss = alpha*out.W_dot_in + W_dot_loss_0 + C_loss*out.N_exp/60*2*pi;
out.W_dot = out.W_dot_in - out.W_dot_loss;
out.W_dot_s = M_dot*(h_su - h_ex_s);
out.epsilon_is = out.W_dot/(M_dot*(h_su - h_ex_s));
out.h_ex1 = max(min((out.M_dot_in*out.h_ex2 + out.M_dot_leak*out.h_su2)/M_dot, out.h_su2), out.h_ex2);
out.T_ex1 = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex1,fluid);
out.cp_ex1 = CoolProp.PropsSI('C','P',P_ex,'H',out.h_ex1,fluid);
out.epsilon_ex1 = max(0,(1-exp(-AU_ex1/(M_dot*out.cp_ex1))));
out.Q_dot_ex = max(0,out.epsilon_ex1*M_dot*out.cp_ex1*(T_w-out.T_ex1));
out.h_ex = min(out.h_ex1 + out.Q_dot_ex/M_dot,param.h_max);
out.Q_dot_amb = AU_amb*(T_w-T_amb);
out.resE = abs((out.Q_dot_su + out.W_dot_loss - out.Q_dot_ex - out.Q_dot_amb)/(out.Q_dot_su + out.W_dot_loss));
out.res = out.resE;

end