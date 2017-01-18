function [out,TS, PV] = PistonExpander(fluid, P_su, h_su, N_exp, P_ex, T_amb, param)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 29/12/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be



%% DEMONSTRATION CASE
if nargin == 0
    
    % Define a demonstration case if ExpanderModel.mat is not executed externally
    fluid = 'R245fa';                           %Nature of the fluid
    N_exp = 2000;                               %Rotational speed       [rpm]
    P_su = 24e+05;                              %Supply pressure        [Pa]
    h_su = 5.0395e+05;                          %Supply enthalpy        [J/kg]
    P_ex = 3.028e+05;                           %Exhaust pressure       [Pa]
    T_amb = 298.1500;                           %Ambient temperature    [K]
    param.displayResults = 1;                   %Flag to control the resuslt display [0/1]
    param.displayFigure = 1;
    
    % Example of paramters for the piston expander model
    param.V = 1e-3;
    param.V_s = 4.551e-05;
    param.V_0 = 6.51e-06;
    param.f_a = 0.17;
    param.f_p = 0.87;
    param.rv_in_exp = 5.882;
    param.rv_in_comp = 6.082;
    param.A_leak = 2.376e-7;
    param.alpha_loss = 1e-2;
    param.W_dot_loss_0 = 50;
    param.C_loss = 0;
    param.M_dot_n = 0.013;
    param.AU_amb =  4;
    param.AU_su_n = 1;
    param.AU_ex_n = 2;
    param.d_su = 0.0016;
    param.d_ex = 0.01;
    %load('C:\Users\RDickes\Google Drive\PhD\MOR study\ORC\Experimental database\Sun2Power\OffDesign\gamma_R245fa.mat');
    %param.gamma.gamma_PQ_pol = gamma_PQ_R245fa; param.gamma.gamma_PT_pol = gamma_PT_R245fa;
end


tstart_exp = tic;

%% EXPANDER MODELING
% Modelling section of the code
if not(isfield(param,'displayResults'))
    param.displayResults = 0;
    %if nothing specified by the user, the results are not displayed by
    %default.
end

if not(isfield(param,'h_min'))
    param.h_min =  CoolProp.PropsSI('H','P',5e4,'T',253.15,fluid);
end
if not(isfield(param,'h_max'))
    param.h_max =  CoolProp.PropsSI('H','P',4e6,'T',500,fluid);
end

T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);



if P_su > P_ex && h_su > CoolProp.PropsSI('H','P',P_su,'Q',0.1,fluid);        
    %If the external conditions are viable (i.e. pressure ratio higher than
    %1 and in two phase conditions), we proceed to the modeling:
    
    
    % Default values assigned to the parameters if not specified
    if not(isfield(param,{'V'}))
        param.V = 0;
    end
    if not(isfield(param,{'alpha_loss'}))
        param.alpha_loss = 0;
    end
    if not(isfield(param,{'W_dot_loss'}))
        param.W_dot_loss = 0;
    end
    if not(isfield(param,{'C_loss'}))
        param.C_loss = 0;
    end
    if not(isfield(param,{'d_su'}))
        param.d_su = 1e2;
    end
    if not(isfield(param,{'d_ex'}))
        param.d_ex = 1e2;
    end
    if not(isfield(param,{'A_leak'}))
        param.A_leak = 0;
    end
    if not(isfield(param,{'AU_su_n'}))
        param.AU_su_n = 0;
    end
    if not(isfield(param,{'AU_ex_n'}))
        param.AU_ex_n = 0;
    end
    if not(isfield(param,{'AU_amb'}))
        param.AU_amb = 0;
    end
    
    ff_guess = [0.8 1.2 0.7 1.3 0.4 1.7 3 ]; %guesses on the filling factor to provide suitable initial point for the iteration
    x_T_guess = [0.7 0.9 0.95 0.8 0.99];%
    stop = 0;
    
    j = 1;
    while not(stop) && j <= length(x_T_guess)
        k =1;
        while not(stop) && k <= length(ff_guess)
            
            x0(1) = ff_guess(k)*param.V_s*param.f_a*N_exp/60*CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid); %initial value for M_dot
            x0(2) = x_T_guess(j)*T_su+(1-x_T_guess(j))*T_amb; %initial value for T_wall
            x0(3) = 0.01*P_su+0.99*P_ex; %initial value for P_ex_3
            x0(4) = CoolProp.PropsSI('H','P',x0(3),'S',CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid),fluid);% initial value for h_ex_3
            
            ub = 2*x0; % upper bound for fsolve
            options = optimset('Display','none');% 'TolX', 1e-8, 'TolFun', 1e-4);
            [x, res, flag] = fsolve(@(x)  Exp_PistonSemiEmp_res(x, ub, fluid, P_su, h_su, N_exp, param.V_s, param.V_0, param.f_a, param.f_p, param.rv_in_exp, param.rv_in_comp, P_ex, param.A_leak, param.d_su, param.alpha_loss, param.W_dot_loss, param.AU_su_n, param.M_dot_n, param.AU_ex_n, param.AU_amb, T_amb, param.C_loss, param.h_min, param.h_max, param), x0./ub, options);
            if norm(res) < 1e-4 %flag > 0
                stop = 1;
            end
            k = k + 1;
        end
        j = j + 1;
    end
    
    x = x.*ub;
    int = Exp_PistonSemiEmp(x(1), x(2), x(3), x(4), fluid, P_su, h_su, N_exp, param.V_s, param.V_0, param.f_a, param.f_p, param.rv_in_exp, param.rv_in_comp, P_ex, param.A_leak, param.d_su, param.alpha_loss, param.W_dot_loss, param.AU_su_n, param.M_dot_n, param.AU_ex_n, param.AU_amb, T_amb, param.C_loss, param.h_min, param.h_max, param);
    M_dot = x(1);
    W_dot = int.W_dot;
    h_ex = int.h_ex;
    epsilon_is = int.epsilon_is;
    FF = M_dot/(param.V_s*N_exp/60*CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid));
    Q_dot_amb = int.Q_dot_amb;
    if int.resE < 1e-4 && int.resMdot <1e-4 && int.resgamma <1e-4 && int.resP <1e-4 && int.resH3 <1e-4 && h_ex > param.h_min && h_ex < param.h_max
        out.flag = flag;
    else
        out.flag = -1;
    end
    out.complete_results = orderfields(int);
    orderfields(int)
else
    out.flag = -2;
end

if out.flag > 0;
    out.h_ex = h_ex;
    out.M_dot = M_dot;
    out.W_dot = W_dot;
    out.Q_dot_amb = Q_dot_amb;
    out.epsilon_is = epsilon_is;
    out.FF = FF;
    out.T_ex = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex,fluid);
    out.M = (CoolProp.PropsSI('D','H',h_su,'P',P_su,fluid)+CoolProp.PropsSI('D','H',out.h_ex,'P',P_ex,fluid))/2*param.V;    
else    
    out.M_dot = N_exp/60*param.V_s*rho_su;
    out.FF = 1;
    out.W_dot = out.M_dot*(h_su-h_ex_s);
    out.epsilon_is =1;
    out.Q_dot_amb = 0;
    out.h_ex = h_ex_s;
    out.T_ex = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex,fluid);
    out.M = (CoolProp.PropsSI('D','H',h_su,'P',P_su,fluid)+CoolProp.PropsSI('D','H',out.h_ex,'P',P_ex,fluid))/2*param.V;
end

if out.flag <=0 && param.displayResults == 1
    warning('The model did not converge correctly - the solution may be wrong')
end
out.time = toc(tstart_exp);

%% TS DIAGRAM and DISPLAY

% Generate the output variable PV 
i = 1;
PV.P(i) = P_su;          
PV.V(i) = param.V_0;      
i = i + 1;

PV.P(i) = out.complete_results.P_su3;
PV.V(i) = param.V_0;
i = i + 1;

PV.P(i) = out.complete_results.P_su3;
PV.V(i) = param.f_a*param.V_s;
i = i + 1;

for r_v = linspace(1,param.rv_in_exp, 20)
    PV.V(i) = param.f_a*param.V_s*r_v;     
    PV.P(i) = CoolProp.PropsSI('P','D', out.complete_results.rho_su3/r_v, 'S',out.complete_results.s_su3,fluid);
    i = i+1;
end

PV.P(i) = out.complete_results.P_ad1;
PV.V(i) = param.V_s;
i = i + 1;

PV.P(i) = out.complete_results.P_ex3;
PV.V(i) = param.V_s;
i = i + 1;

PV.P(i) = P_ex;
PV.V(i) = param.V_s;
i = i + 1;

PV.P(i) = out.complete_results.P_ex3;
PV.V(i) = param.V_s;
i = i + 1;

PV.P(i) = out.complete_results.P_ex3;
PV.V(i) = param.f_p*param.V_s;
i = i + 1;

for r_v = linspace(1,param.rv_in_comp, 20)
    PV.V(i) = param.f_p*param.V_s/r_v;     
    PV.P(i) = CoolProp.PropsSI('P','D', out.complete_results.rho_ex3*r_v, 'S',out.complete_results.s_ex3,fluid);
    i = i+1;
end

PV.P(i) = out.complete_results.P_ad2;
PV.V(i) = param.V_0;
i = i + 1;

PV.P(i) = out.complete_results.P_su3;
PV.V(i) = param.V_0;

% Generate the output variable TS 
s_su = CoolProp.PropsSI('S','H',h_su,'P',P_su, fluid);
s_ex = CoolProp.PropsSI('S','H', out.h_ex,'P', P_ex, fluid);
TS.T = [T_su out.T_ex];
TS.s = [s_su s_ex];

if param.displayFigure
    figure
    plot(PV.V, PV.P)
end

% If the param.displayResults flag is activated (=1), the results are 
% displayed on the command window
if param.displayResults ==1
    in.fluid = fluid;
    in.T_su = T_su;
    in.h_su = h_su;
    in.P_su = P_su;
    in.P_ex = P_ex;
    in.P_ex = P_ex;
    in.N_exp = N_exp;
    in.T_amb = T_amb;
    out = orderfields(out);
    if nargin ==0
        fprintf ( 1, '\n' );
        disp('-------------------------------------------------------')
        disp('--------------------   Demo Code   --------------------')
        disp('-------------------------------------------------------')
        fprintf ( 1, '\n' );
    end
    disp('Working conditions:')
    fprintf ( 1, '\n' );
    disp(in)
    disp('Results')
    fprintf ( 1, '\n' );
    disp(out)
    disp('Detailed Results')
    fprintf ( 1, '\n' );
    disp(out.complete_results)
end

end

%% NESTED FUNCTIONS

function res = Exp_PistonSemiEmp_res(x, ub, fluid, P_su, h_su, N_exp, V_s, V_0, f_a, f_p, rv_in_exp, rv_in_comp, P_ex, A_leak, d_su, alpha_loss, W_dot_loss, AU_su_n, M_dot_n, AU_ex_n, AU_amb, T_amb, C_loss, h_min, h_max, param)
x = x.*ub;
out = Exp_PistonSemiEmp(x(1), x(2), x(3), x(4), fluid, P_su, h_su, N_exp, V_s, V_0, f_a, f_p, rv_in_exp, rv_in_comp, P_ex, A_leak, d_su, alpha_loss, W_dot_loss, AU_su_n, M_dot_n, AU_ex_n, AU_amb, T_amb, C_loss, h_min, h_max, param);
res = [out.resE ; out.resMdot; out.resP; out.resH3];
end

function out = Exp_PistonSemiEmp(M_dot, T_w, P_ex3, h_ex3, fluid, P_su, h_su, N_exp, V_s, V_0, f_a, f_p, rv_in_exp, rv_in_comp, P_ex, A_leak, d_su, alpha_loss, W_dot_loss, AU_su_n, M_dot_n, AU_ex_n, AU_amb, T_amb, C_loss, h_min, h_max, param)
if isnan(T_w)
    T_w =  CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
end
M_dot = max(1e-5, M_dot);
out.T_w = T_w;
P_ex3 = min(inf, max(P_ex3, P_ex-1e5));

% Supply conditions
s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
out.h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
AU_su = AU_su_n*(M_dot/M_dot_n)^0.8;
AU_ex = AU_ex_n*(M_dot/M_dot_n)^0.8;

% su -> su1 : intake pressure drop
out.h_su1 = h_su;
out.P_su1 = max(CoolProp.PropsSI('P','S',s_su,'H',max(h_su - (M_dot/(pi*d_su^2/4*rho_su))^2/2, out.h_ex_s), fluid),P_ex+1);
out.T_su1 = CoolProp.PropsSI('T','P',out.P_su1,'H',out.h_su1,fluid);
out.DP_su = P_su-out.P_su1;

% su1 -> su2 : supply heat transfer
if CoolProp.PropsSI('Q','P',out.P_su1,'H',out.h_su1,fluid) == -1
    out.cp_su1 = CoolProp.PropsSI('C','P',out.P_su1,'H',out.h_su1,fluid);
else
    out.cp_su1 = CoolProp.PropsSI('C','P',out.P_su1,'Q',0,fluid);
end
out.epsilon_su1 = max(0,(1-exp(-AU_su/(M_dot*out.cp_su1))));
out.Q_dot_su = max(-inf,out.epsilon_su1*M_dot*out.cp_su1*(out.T_su1 - T_w));
out.h_su2 = min(h_max,max(max(out.h_ex_s, CoolProp.PropsSI('H','Q',0.1,'P',out.P_su1,fluid)),out.h_su1 - out.Q_dot_su/M_dot));
out.s_su2 = CoolProp.PropsSI('S','H',out.h_su2,'P',out.P_su1,fluid);
out.rho_su2 = CoolProp.PropsSI('D','H',out.h_su2,'P',out.P_su1,fluid);
out.P_su2 = out.P_su1;

% su2 -> leak : leakages
%out.Q_su2 = CoolProp.PropsSI('Q','P',out.P_su2,'H',out.h_su2,fluid);
%if  out.Q_su2 < 0
%    out.gamma  = feval(param.gamma.gamma_PT_pol,[out.P_su2/1e5, CoolProp.PropsSI('T','P',out.P_su2,'H',out.h_su2,fluid)/1e2]);
%else
%    out.gamma = feval(param.gamma.gamma_PQ_pol,[out.P_su2/1e5, out.Q_su2]);
%end
%out.gamma = max(1e-2,out.gamma);
%out.P_crit2 = out.P_su2*(2/(out.gamma+1))^(out.gamma/(out.gamma-1));

f_leak_c_sound = @(x) sound_speed(x, out.s_su2, out.h_su2, fluid);
c_guess = CoolProp.PropsSI('A','H',out.h_su2,'S',out.s_su2,fluid);
lb = 0.5*c_guess;
ub = 2*c_guess;
out.c_sound_leak = zeroBrent ( lb, ub, 1e-6, 1e-6, f_leak_c_sound );
out.h_crit = out.h_su2 - out.c_sound_leak^2/2;
out.P_crit = CoolProp.PropsSI('P','H',out.h_crit,'S',out.s_su2,fluid);
out.resgamma = 0;
out.P_thr = max(P_ex3,out.P_crit);
out.rho_thr = CoolProp.PropsSI('D','P',out.P_thr,'S',out.s_su2,fluid);
out.C_thr = sqrt(2*(out.h_su2 - CoolProp.PropsSI('H','P',out.P_thr,'S',out.s_su2,fluid)));
out.M_dot_leak = A_leak*out.C_thr*out.rho_thr;

% ex3 -> ad2 -> su2 : recompression
out.P_ex3 = P_ex3;
out.rho_ex3 = CoolProp.PropsSI('D','P',P_ex3,'H',h_ex3,fluid);
out.s_ex3 = CoolProp.PropsSI('S','P',P_ex3,'H',h_ex3,fluid);
out.M_dot_comp = N_exp/60*V_s*f_p*out.rho_ex3;
out.rho_ad2 = out.rho_ex3*rv_in_comp;
try
    out.P_ad2 = CoolProp.PropsSI('P','D',out.rho_ad2,'S',out.s_ex3,fluid);
catch
    delta = 0.0001;
    out.P_ad2 = 0.5*CoolProp.PropsSI('P','D',out.rho_ad2*(1+delta),'S',out.s_ex3,fluid)+0.5*CoolProp.PropsSI('P','D',out.rho_ad2*(1-delta),'S',out.s_ex3,fluid);
end
out.h_ad2 = CoolProp.PropsSI('H','D',out.rho_ad2,'P',out.P_ad2,fluid);
out.W_dot_3 = out.M_dot_comp*(out.h_ad2-h_ex3);
out.W_dot_4 = out.M_dot_comp/out.rho_ad2*(out.P_su2-out.P_ad2);

% su2 + ad2 -> su3 : swept mixing
out.M_dot_in = M_dot-out.M_dot_leak;
M_dot_exp_bis = out.M_dot_comp+out.M_dot_in;
out.h_su3 = (out.M_dot_comp*out.h_ad2 + out.M_dot_in*out.h_su2)/M_dot_exp_bis;
out.P_su3 = out.P_su2;
out.rho_su3 = CoolProp.PropsSI('D','P',out.P_su3,'H',out.h_su3,fluid);
out.s_su3 = CoolProp.PropsSI('S','P',out.P_su3,'H',out.h_su3,fluid);

% su3 -> ad1 -> ex3 : expansion
out.M_dot_exp = N_exp/60*V_s*f_a*out.rho_su3;
out.rho_ad1 = out.rho_su3/rv_in_exp;
try
    out.P_ad1 = CoolProp.PropsSI('P','D',out.rho_ad1,'S',out.s_su3,fluid);
catch
    delta = 0.0001;
    out.P_ad1 = 0.5*CoolProp.PropsSI('P','D',out.rho_ad1*(1+delta),'S',out.s_su3,fluid)+0.5*CoolProp.PropsSI('P','D',out.rho_ad1*(1-delta),'S',out.s_su3,fluid);
end
out.h_ad1 = CoolProp.PropsSI('H','D',out.rho_ad1,'P',out.P_ad1,fluid);
out.W_dot_1 = out.M_dot_exp*(out.h_su3-out.h_ad1);
out.W_dot_2 = out.M_dot_exp/out.rho_ad1*(out.P_ad1 - P_ex3);
out.W_dot_in = out.W_dot_1+out.W_dot_2-out.W_dot_3-out.W_dot_4;
out.h_ex3 = out.h_su2 - out.W_dot_in/out.M_dot_in;

% ex3 + leak -> ex2 : exhaust mixing
out.h_ex2 = max(min((out.M_dot_in*h_ex3 + out.M_dot_leak*out.h_su2)/M_dot, out.h_su2), out.h_ex3);

% ex2 -> ex1 : exhasut heat transfer
out.P_ex2 = P_ex3;
out.T_ex2 = CoolProp.PropsSI('T','P',out.P_ex2,'H',out.h_ex2,fluid);
out.cp_ex2 = CoolProp.PropsSI('C','P',out.P_ex2,'H',out.h_ex2,fluid);
out.epsilon_ex2 = max(-inf,(1-exp(-AU_ex/(M_dot*out.cp_ex2))));
out.Q_dot_ex = max(-inf,out.epsilon_ex2*M_dot*out.cp_ex2*(T_w-out.T_ex2));
out.h_ex1 = min(out.h_ex2 + out.Q_dot_ex/M_dot,h_max);
out.P_ex1 = out.P_ex2;

% ex1 -> ex : exhaust pressure drop
out.s_ex1 = CoolProp.PropsSI('S','P',out.P_ex1,'H',out.h_ex1,fluid);
out.rho_ex1 = CoolProp.PropsSI('D','P',out.P_ex1,'H',out.h_ex1,fluid);
P_ex_bis = max(CoolProp.PropsSI('P','S',out.s_ex1,'H',max(out.h_ex1 - (M_dot/(pi*param.d_ex^2/4*out.rho_ex1))^2/2, 0.5*out.h_ex_s), fluid),-inf);
out.h_ex = out.h_ex1;

% Energy balance
out.W_dot_loss = alpha_loss*out.W_dot_in + W_dot_loss + C_loss*N_exp/60*2*pi;
out.W_dot = out.W_dot_in - out.W_dot_loss;
out.W_dot_s = M_dot*(h_su - out.h_ex_s);
out.epsilon_is = out.W_dot/out.W_dot_s;
out.Q_dot_amb = AU_amb*(T_w-T_amb);

% Residuals
out.resP = 1-P_ex_bis/P_ex;
out.resH3 = 1- h_ex3/out.h_ex3;
out.resE = abs((out.Q_dot_su + out.W_dot_loss - out.Q_dot_ex - out.Q_dot_amb)/(out.Q_dot_su + out.W_dot_loss));
if M_dot_exp_bis == 0 && out.M_dot_exp == 0
    out.resMdot = 0;
elseif M_dot_exp_bis ~= 0 && - out.M_dot_exp == 0;
    out.resMdot = 1;
elseif M_dot_exp_bis == 0 && - out.M_dot_exp ~= 0;
    out.resMdot = 1;
else
    out.resMdot = abs((out.M_dot_exp - M_dot_exp_bis)/M_dot_exp_bis);
end

end

function res = sound_speed(x, s_su, h_su, fluid)
c_guess = x;
h_crit = h_su - c_guess^2/2;
c_real = CoolProp.PropsSI('A', 'S', s_su, 'H', h_crit, fluid);
res = abs(c_guess-c_real)/c_real;
end
