function out = CompressorNestedModel_3_vRD_2(fluid, P_su, T_su, M_dot, P_ex, T_amb, param)

% Compressor model based on Vincent Lemort's model 
%%
% The model inputs are:
%       - P_su: inlet pressure (P_r_in_cp)                          [bar]
%       - P_ex: outlet pressure (P_r_in_cp)                         [bar]
%       - T_su: inlet temperature (T_r_in_cp)                       [K]
%       - fluid: nature of the fluid (string)                       [-]
%       - N_rot: rotational speed                                   [rpm]
%       - T_amb : ambient temperature                               [K]
%       - param: structure variable containing the model parameters
%
% The model parameters provided in 'param' should contain the followings:
%           param.AU_amb, global heat transfer coefficient for ambient heat losses      [W/K];
%           param.AU_su_n, global heat transfer coefficient for supply heat transfer    [W/K];
%           param.AU_ex_n, global heat transfer coefficient for exhaust heat transfer	[W/K];
%           param.A_leak, nozzle cross section area for the leakage                     [m?];
%           param.W_dot_loss_0, constant losses term                                    [W];
%           param.alpha, proportional losses coefficient                                [-];
%           param.r_v_in_cp, built-in volumetric ratio                                  [-];
%           param.V_s, machine displacement volume                                      [m^3];
%           param.gamma, correlations to compute the isentropic component of the fluid

% The model outputs are:
%       - out: a structure variable which includes
%               - T_ex = exhaust temperature                    [K]
%               - M_dot = fluid mass flow rate                  [kg/s]
%               - W_dot = mechanical power                      [W]
%               - epsilon_is = isentropic efficiency            [-]
%




%% Input

if nargin == 0
    clear all
    close all
    clc
    
    fluid = 'R1234yf';                              % Nature of the fluid
    M_dot = 0.01;                                   % Mass flow rate
    P_su = 437500;                                  % Supply pressure        [Pa]
    T_su = 288.92;                                  % Supply temperature     [K]
    P_ex = 928600;                                  % Exhaust pressure       [Pa]
    T_amb = 298.1500;                               % Ambient temperature    [K]
    
    
    param.AU_amb = 5;
    param.AU_su_n = 10;
    param.AU_ex_n = 10;
    param.A_leak0 = 1e-8;
    param.W_dot_loss_0 = 10;
    param.alpha = 0.01;
    param.r_v_in = 2;
    param.M_dot_n = 0.09;
    param.V_s = 0.000034; 
    if not(isfield(param,'h_min'))
        param.h_min =  CoolProp.PropsSI('H','P',5e4,'T',253.15,fluid);
    end
    if not(isfield(param,'h_max'))
        param.h_max =  CoolProp.PropsSI('H','P',3e7,'T',400,fluid);
    end

end

%% MODELLING SECTION

% Modeling section of the code
s_su = CoolProp.PropsSI('S','P',P_su,'T',T_su,fluid);
rho_su = CoolProp.PropsSI('D','P',P_su,'T',T_su,fluid);
h_su = CoolProp.PropsSI('H','P',P_su,'T',T_su,fluid);
h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
T_ex_s = CoolProp.PropsSI('T','P',P_ex,'S',s_su,fluid);

x_T_guess = [0.9 0.8 0.7 0.2];
stop = 0;
k_iter = 0;
while not(stop) && k_iter<min(length(x_T_guess),20)
   k_iter = k_iter + 1;
   x0(1) = x_T_guess(k_iter)*min(T_su,T_amb) + (1-x_T_guess(k_iter))*max(T_amb,T_ex_s);  % Initial value for T_w 
   x0(2) = CoolProp.PropsSI('T','P',P_ex,'S',s_su,fluid); % Initial value for T_ex_2

    
   lb = [min(T_su,T_amb)-10     x0(2)-50];
   ub = [max(T_amb,T_ex_s)+50   x0(2)+50];
    
   options = optimoptions('fmincon','Display','none');
   [x, ~, flag] = fmincon(@(x)  FCT_Compr_SemiEmp_res(x, ub, fluid, M_dot, T_amb, P_su, T_su, h_su, P_ex, h_ex_s, param.V_s, param.r_v_in, param.A_leak0, param.alpha, param.W_dot_loss_0, param.AU_su_n, param.M_dot_n, param.AU_ex_n, param.AU_amb,param.h_min,param.h_max),x0./ub, [],[],[],[],lb./ub, ub./ub, [], options);
   if flag > 0 
      stop = 1;
   end
end
x = x.*ub;
int = FCT_Compr_SemiEmp(x(1), x(2), fluid, M_dot, T_amb, P_su, T_su, h_su, P_ex, h_ex_s, param.V_s, param.r_v_in, param.A_leak0, param.alpha, param.W_dot_loss_0, param.AU_su_n, param.M_dot_n, param.AU_ex_n, param.AU_amb,param.h_min,param.h_max);

out.M_dot = M_dot;
out.N_rot = int.N_rot;
out.W_dot = int.W_dot;
out.T_ex = CoolProp.PropsSI('T','H',int.h_ex,'P',P_ex,fluid);
out.VolEff = out.M_dot/(param.V_s*int.N_rot/60*rho_su);
out.Q_dot_amb = int.Q_dot_amb;
out.epsilon_is = int.epsilon_is;
if int.resE < 1e-4 && int.res_hex2 < 1e-4  
    out.flag = flag;
else
    out.flag = -1;
end
out.complete_results = int;
    
if out.flag <=0 
   warning('The model did not converge correctly - the solution may be wrong')
end
    
end

%% NESTED FUNCTIONS

function res = FCT_Compr_SemiEmp_res(x, ub, fluid, M_dot, T_amb, P_su, T_su, h_su, P_ex, h_ex_s, V_s, r_v_in, A_leak0, alpha, W_dot_loss_0, AU_su_n, M_dot_n, AU_ex_n, AU_amb,h_min,h_max)

x = x.*ub;
out = FCT_Compr_SemiEmp(x(1), x(2), fluid, M_dot, T_amb, P_su, T_su, h_su, P_ex, h_ex_s, V_s, r_v_in, A_leak0, alpha, W_dot_loss_0, AU_su_n, M_dot_n, AU_ex_n, AU_amb,h_min,h_max);
res = [out.resE; out.res_hex2];
res = norm(res);
end

function out = FCT_Compr_SemiEmp(T_w, T_ex2, fluid, M_dot, T_amb, P_su, T_su, h_su, P_ex, h_ex_s, V_s, r_v_in, A_leak0, alpha, W_dot_loss_0, AU_su_n, M_dot_n, AU_ex_n, AU_amb,~,~)

gamma = 1.03;
out.rho_su = CoolProp.PropsSI('D','P',P_su,'T',T_su,fluid);
out.c_p_su = CoolProp.PropsSI('C','P',P_su,'T',T_su,fluid);
out.C_dot_su = M_dot * out.c_p_su;
out.AU_su = AU_su_n * (M_dot/M_dot_n)^0.8;
out.NTU_su = out.AU_su / out.C_dot_su;
out.epsilon_su = 1 - exp(-out.NTU_su);
out.Q_dot_su = out.epsilon_su * out.C_dot_su *(T_w - T_su);
out.h_su1 = h_su + out.Q_dot_su / M_dot;
out.P_su1 = P_su; 
out.T_su1 = CoolProp.PropsSI('T','P',out.P_su1,'H',out.h_su1,fluid);

% Internal leakage: ex2 -> su1
out.h_ex2_bis = CoolProp.PropsSI('H','P',P_ex,'T',T_ex2,fluid);
out.s_ex2_bis = CoolProp.PropsSI('S','P',P_ex,'T',T_ex2,fluid);
out.s_thr = out.s_ex2_bis; 
out.P_thr_crit = P_ex*(2/(gamma+1))^(gamma/(gamma-1)); 
out.P_thr = max(out.P_su1,out.P_thr_crit);
out.rho_thr = CoolProp.PropsSI('D','P',out.P_thr,'S',out.s_thr,fluid);
h_thr = CoolProp.PropsSI('H','S',out.s_thr,'P',out.P_thr,fluid);
C_thr = sqrt(2*(out.h_ex2_bis - h_thr));
out.V_dot_leak = A_leak0 * C_thr;
out.M_dot_leak = A_leak0 * C_thr * out.rho_thr;

% Flow rate calculation: su2
out.M_dot_in = out.M_dot_leak + M_dot;
out.P_su2 = out.P_su1;
out.h_su2 = (M_dot * out.h_su1 + out.M_dot_leak*out.h_ex2_bis)/out.M_dot_in;
out.rho_su2 = CoolProp.PropsSI('D','H',out.h_su2,'P',out.P_su2,fluid);
out.s_su2 = CoolProp.PropsSI('S','H',out.h_su2,'P',out.P_su2,fluid);
out.N_rot = out.M_dot_in/V_s/out.rho_su2*60;

% Isentropic compression
out.s_in = out.s_su2; 
out.rho_in = out.rho_su2 * r_v_in; 
out.P_in = CoolProp.PropsSI('P','D',out.rho_in,'S',out.s_in,fluid);
out.h_in = CoolProp.PropsSI('H','D',out.rho_in,'P',out.P_in,fluid);
out.w_in1 = out.h_in - out.h_su2;
out.r_p_in = out.P_in/out.P_su2;

% Isochoric compression 
out.rho_ex2 = out.rho_in;
out.P_ex2 = P_ex;
out.T_ex2 = CoolProp.PropsSI('T','P',out.P_ex2,'D',out.rho_ex2,fluid);
out.h_ex2 = CoolProp.PropsSI('H','P',out.P_ex2,'D',out.rho_ex2,fluid);
out.s_ex2 = CoolProp.PropsSI('S','P',out.P_ex2,'D',out.rho_ex2,fluid);
out.w_in2 = (out.P_ex2 - out.P_in)/out.rho_in; 

% Total compression work and power
w_in = out.w_in1 + out.w_in2;
out.W_dot_in = out.M_dot_in * w_in;

% Electrical consuption
out.W_dot_loss = alpha*out.W_dot_in + W_dot_loss_0; %Constant electro-mechanical losses + losses proportional to internal power
out.W_dot = out.W_dot_in + out.W_dot_loss;


% Exhaust cooling-down: ex2 -> ex
c_p_ex = CoolProp.PropsSI('C','P',P_ex,'H',out.h_ex2,fluid);
AU_ex = AU_ex_n *(M_dot/M_dot_n)^0.8;
C_dot_ex = M_dot * c_p_ex;
NTU_ex = AU_ex / C_dot_ex;
epsilon_ex = 1 - exp(-NTU_ex);
out.Q_dot_ex = epsilon_ex * C_dot_ex * (T_w - out.T_ex2);
out.h_ex = out.h_ex2 + (out.Q_dot_ex / M_dot);


% Finctitious enveloppe heat balance
out.Q_dot_amb = AU_amb * (T_w - T_amb);
out.T_w = T_w;

%Global isentropic effectiveness
out.epsilon_is =(M_dot*(h_ex_s - h_su))/out.W_dot;

% Residuals
out.resE = abs((out.W_dot_loss - out.Q_dot_ex - out.Q_dot_su - out.Q_dot_amb)/(out.W_dot_loss));
out.res_hex2 = abs(out.h_ex2_bis - out.h_ex2)/out.h_ex2;

end

