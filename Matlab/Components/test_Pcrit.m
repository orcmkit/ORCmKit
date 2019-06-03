function test_Pcrit
clc

% case study
fluid = 'R245fa';
T_su = 150+273.15;
P_su = 15e5;
P_ex = 2e5;
h_su = CoolProp.PropsSI('H','T', T_su, 'P',P_su,fluid);
s_su = CoolProp.PropsSI('S','T', T_su, 'P',P_su,fluid);
rho_su = CoolProp.PropsSI('D','T', T_su, 'P',P_su,fluid);
A_thr = 1e-5;

% % Method 1 : fictive gamma
f_gamma = @(x) gamma_fct(x, P_ex, P_su, s_su, rho_su, fluid);
lb = 0.1;
ub = 2;
gamma = zeroBrent ( lb, ub, 1e-6, 1e-6, f_gamma, 1e-6);
P_thr_1 = max(P_ex, P_su *(2/(gamma+1))^(gamma/(gamma-1)));
rho_thr_1 = CoolProp.PropsSI('D','P',P_thr_1,'S',s_su,fluid);
C_thr_1 = sqrt(2*(h_su - CoolProp.PropsSI('H','P',P_thr_1,'S',s_su,fluid)));
M_dot_thr_1= A_thr*C_thr_1*rho_thr_1;


% Method 2: speed of sound
f_leak_c_sound = @(x) sound_speed(x, s_su, h_su, fluid);
c_guess = CoolProp.PropsSI('A','H',h_su,'S',s_su,fluid);
lb = 0.5*c_guess;
ub = 2*c_guess;
c_sound_leak = zeroBrent ( lb, ub, 1e-6, 1e-6, f_leak_c_sound , 1e-6);
h_crit_2 = h_su - c_sound_leak^2/2;
P_crit_2 = CoolProp.PropsSI('P','H',h_crit_2,'S',s_su,fluid);
P_thr_2 = max(P_ex,P_crit_2);
rho_thr_2 = CoolProp.PropsSI('D','P',P_thr_2,'S',s_su,fluid);
C_thr_2 = sqrt(2*(h_su - CoolProp.PropsSI('H','P',P_thr_2,'S',s_su,fluid)));
M_dot_thr_2 = A_thr*C_thr_2*rho_thr_2;

disp(['[P_thr_1 = ' num2str(P_thr_1) '- P_thr_2 = ' num2str(P_thr_2) ']'])
disp(['[C_thr_1 = ' num2str(C_thr_1) '- C_thr_2 = ' num2str(C_thr_2) ']'])
disp(['[rho_thr_1 = ' num2str(rho_thr_1) '- rho_thr_2 = ' num2str(rho_thr_2) ']'])
disp(['[M_dot_thr_1 = ' num2str(M_dot_thr_1) '- M_dot_thr_2 = ' num2str(M_dot_thr_2) ']'])

end


function res = sound_speed(x, s_su, h_su, fluid)
c_guess = x;
h_crit = h_su - c_guess^2/2;
c_real = CoolProp.PropsSI('A', 'S', s_su, 'H', h_crit, fluid);
res = abs(c_guess-c_real)/c_real;
end


function res = gamma_fct(x, P_ex, P_su, s_su, rho_su, fluid)
gamma = max(0.1,min(x,3));
P_thr = max(P_ex, P_su *(2/(gamma+1))^(gamma/(gamma-1)));
rho_thr = CoolProp.PropsSI('D','P',P_thr,'S',s_su,fluid);
gamma_bis =  log10(P_su/P_thr)/log10(rho_su/rho_thr);
res = gamma_bis-gamma;
end