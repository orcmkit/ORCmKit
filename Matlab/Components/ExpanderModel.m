function [out,TS] = ExpanderModel(fluid, P_su, h_su, N_exp, P_ex, T_amb, Expander)
%fluid, P_su, h_su, N_exp, P_ex, T_amb
%% CODE INFORMATION
% Rémi Dickes - 26/01/2015 (University of Liège, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% ExpanderModel is a single matlab code implementing four different 
% modeling approaches to simulate a volumetric pump (see the Documentation)
%
% The inputs are the followings:
%       - P_su: inlet pressure of the WF                          	[Pa]
%       - h_su: inlet temperature of the WF                        	[K]
%       - P_ex: outlet pressure of the WF                          	[Pa]
%       - fluid: nature of the fluid                               	[string]
%       - N_exp: the pump speed                                    	[rpm]
%       - Expander: structure variable containing the model parameters
%
% The variable param depends of the type of model selected:
%
%       - if param.modelType = 'CstEff':
%           param.V_s , machine displacement volume                	[m³]
%           param.epsilon_is, isentropic efficiency                	[-]
%           param.epsilon_vol, volumetric efficiency               	[-]
% 
%       - if param.modelType = 'PolEff':
%           param.V_s , machine displacement volume                	[m³]
%           param.N_pp_nom, pump nominal shaft speed              	[rpm]
%           param.coeffPol_is, polynmial coef for epsilon_is        [-]
%           param.coeffPol_vol, polynmial coef for epsilon_vol      [-]
% 
%       - if param.modelType = 'SemiEmp':
%           param.V_s , machine displacement volume               	[m³]
%           param.A_leak, leakage surface area                     	[m²]
%           param.W_dot_loss, constant power losses                	[W]
%           param.K_0_loss, term for the proportional losses       	[-]
%
% See the documentation for further info or contact rdickes@ulg.ac.be


tstart_exp = tic;
%% DEMONSTRATION CASE

if nargin == 0
    clc
    clear all
    close all
    fluid = 'R245fa';
    N_exp = 2000;
    P_su = 1.523372634282432e+06;
    h_su = 3.654377176733231e+05; 
    T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid)
    T_sat = CoolProp.PropsSI('T','P',P_su,'Q',0,fluid)
    Q_su = CoolProp.PropsSI('Q','P',P_su,'H',h_su,fluid)
    P_ex = 2.585969394779355e+05;
    T_amb = 2.831500000000000e+02;
    path = 'C:\Users\RDickes\Google Drive\PhD\MOR study\ORC\Experimental database\Sun2Power';
    EXP_folder = [path '\Expander\'];
    load([EXP_folder, 'ParametersCalibration_EXP.mat'])
    Expander = EXP_SemiEmp;
%     Expander.V_s = 1.279908675799087e-05;
%     Expander.r_v_in = 2.19;
%     Expander.A_leak0 = 1.344675598144531e-06;
%     Expander.A_leak_n = 0;
%     Expander.d_su = 0.003227564086914;
%     Expander.alpha = 1.253662109374984e-05;
%     Expander.W_dot_loss_0 = 0;
%     Expander.C_loss = 7.952880859375330e-07;
%     Expander.AU_su_n = 50.033569335937500;
%     Expander.M_dot_n =  0.068378356905196;
%     Expander.P_su_n = 1.417898896236315e+06;
%     Expander.AU_ex_n = 94.017028808593750;
%     Expander.AU_amb = 0.674005126953125;
%     Expander.epsilon_s = 0.7;
%     Expander.FF = 1.2;
%     Expander.Pac = [0.894464014706633   0.843487146881312   2.496197210245431   0.698814176031069   1.750267746114657  -0.164103975268144   0.077850624779040  -0.173994287900423];
%     Expander.Dec = [0.986159987784339   0.026443495944529   0.123726882289625];
%     Expander.P_su_n = 10e5;
%     Expander.rp_n = 4;
%     Expander.aPole_s = [1.504191277585839   0.051680226738587  -0.000001018959447  -0.020713559079583   0.000000031728135   0.000000000000223];
%     Expander.aPolFF = [1.122789891166860  -0.050184374872146  -0.000000109428161   0.001845567200683   0.000000026842863   0.000000000000043];
%     Expander.modelType = 'SemEmp_rv_mc_lk_dp_ht';
%     Expander.displayResults =1;
end

%% INPUTS VERIFICATION
if not(isfield(Expander,'displayResults'))
    Expander.displayResults = 0;
end
tstart_exp = tic;

%% EXPANDER MODELING
if P_su > P_ex && h_su > CoolProp.PropsSI('H','P',P_su,'Q',0.01,fluid);
    switch Expander.modelType
        case 'CstEff'
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
            out.h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
            out.h_ex = h_su - Expander.epsilon_s*(h_su - out.h_ex_s);
            out.T_ex = CoolProp.PropsSI('T','H',out.h_ex,'P',P_ex,fluid);
            out.epsilon_s = Expander.epsilon_s;
            out.flag = 1;
            
        case 'CstEffVs'
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
            rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
            out.h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
            out.M_dot_in = N_exp/60*Expander.V_s*Expander.FF*rho_su;
            out.M_dot = out.M_dot_in;
            out.W_dot_in =  out.M_dot*(h_su-out.h_ex_s)*Expander.epsilon_s;
            out.W_dot = out.W_dot_in;
            out.Q_dot_amb = Expander.AU_amb*(T_su-T_amb);
            out.h_ex = h_su - (out.W_dot_in+out.Q_dot_amb)/out.M_dot;
            out.T_ex = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex,fluid);
            out.FF = Expander.FF;
            out.epsilon_s = Expander.epsilon_s;
            out.flag = 1;
            
        case 'Pol2Eff'
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
            rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
            out.h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
            aPole_s = Expander.aPole_s;
            aPolFF = Expander.aPolFF;
            if length(aPolFF) == 6
                FF = max(0.5, min(aPolFF(1) + aPolFF(2)*(P_su/P_ex) + aPolFF(3)*rho_su + aPolFF(4)*(P_su/P_ex)^2 + aPolFF(5)*(P_su/P_ex)*rho_su + aPolFF(6)*(rho_su)^2, 5));
            elseif length(aPolFF) == 10
                FF = max(0.3,min(aPolFF(1) + aPolFF(2)*(P_su/P_ex) + aPolFF(3)*rho_su + aPolFF(4)*(P_su/P_ex)^2 + aPolFF(5)*(P_su/P_ex)*rho_su + aPolFF(6)*(rho_su)^2 + aPolFF(7)*N_exp + aPolFF(8)*N_exp^2 + aPolFF(9)*N_exp*(P_su/P_ex) + aPolFF(10)*N_exp*rho_su,5));

            end
            if length(aPole_s) == 6
                epsilon_s = max(0.01,min(aPole_s(1) + aPole_s(2)*(P_su/P_ex) + aPole_s(3)*rho_su + aPole_s(4)*(P_su/P_ex)^2 + aPole_s(5)*(P_su/P_ex)*rho_su + aPole_s(6)*(rho_su)^2,1));
            elseif length(aPole_s) == 10
                epsilon_s = max(0.01,min(aPole_s(1) + aPole_s(2)*(P_su/P_ex) + aPole_s(3)*rho_su + aPole_s(4)*(P_su/P_ex)^2 + aPole_s(5)*(P_su/P_ex)*rho_su + aPole_s(6)*(rho_su)^2 + aPole_s(7)*N_exp + aPole_s(8)*N_exp^2 + aPole_s(9)*N_exp*(P_su/P_ex) + aPole_s(10)*N_exp*rho_su,1));
            end
            
            
            out.M_dot_in = N_exp/60*Expander.V_s*FF*rho_su;
            out.M_dot = out.M_dot_in;
            out.W_dot_in = out.M_dot*(h_su - out.h_ex_s)*epsilon_s;
            out.W_dot = out.W_dot_in;
            out.Q_dot_amb = max(0,Expander.AU_amb*(T_su - T_amb));
            out.h_ex = h_su - (out.W_dot_in + out.Q_dot_amb)/out.M_dot;
            out.T_ex = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex,fluid);
            out.FF = FF;
            out.epsilon_s = epsilon_s;
            out.flag = 1;
            
        case 'EmpEff'
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
            rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
            out.h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
            rp_0 = Expander.Pac(1);
            delta = Expander.Pac(2) + Expander.Pac(6)*(P_su-Expander.P_su_n)/Expander.P_su_n;
            rp_max = Expander.Pac(3) + Expander.Pac(7)*(P_su-Expander.P_su_n)/Expander.P_su_n;
            y_max = Expander.Pac(4) + Expander.Pac(8)*((P_su-Expander.P_su_n)/Expander.P_su_n);
            xi = Expander.Pac(5);
            B = delta/xi/y_max;
            E = (B*(rp_max-rp_0) - tan(pi/(2*xi)))/(B*(rp_max-rp_0) - atan(B*(rp_max-rp_0)));
            epsilon_s = max( 1e-5, min(y_max*sin(xi*atan(B*(P_su/P_ex-rp_0) - E*(B*(P_su/P_ex-rp_0) - atan(B*(P_su/P_ex-rp_0))))), 1));
            FF = max(0.8, min(5,Expander.Dec(1) + Expander.Dec(2)*(((P_su/P_ex)-Expander.rp_n)/Expander.rp_n) + Expander.Dec(3)*((P_su-Expander.P_su_n)/Expander.P_su_n)));
            out.M_dot_in = N_exp/60*Expander.V_s*FF*rho_su;
            out.M_dot = out.M_dot_in;
            out.W_dot_in = out.M_dot*(h_su - out.h_ex_s)*epsilon_s;
            out.W_dot = out.W_dot_in;
            out.Q_dot_amb = Expander.AU_amb*(T_su - T_amb);
            out.h_ex = h_su - (out.W_dot_in + out.Q_dot_amb)/out.M_dot;
            out.T_ex = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex,fluid);
            out.FF = FF;
            out.epsilon_s = epsilon_s;
            out.flag = 1;
            
        case 'SemEmp_rv_mc'
            
            if not(isfield(Expander,{'alpha'}))
                Expander.alpha = 0;
            end
            if not(isfield(Expander,{'W_dot_loss_0'}))
                Expander.W_dot_loss_0 = 0;
            end
            if not(isfield(Expander,{'C_loss'}))
                Expander.C_loss = 0;
            end
            
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
            rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
            out.M_dot_in = N_exp/60*Expander.V_s*rho_su;
            out.M_dot = out.M_dot_in;
            out.rho_in = rho_su/Expander.r_v_in;
            out.P_in = CoolProp.PropsSI('P','D',out.rho_in,'S',s_su,fluid);
            out.h_in = CoolProp.PropsSI('H','P',out.P_in,'S',s_su,fluid);
            out.w_1 = h_su- out.h_in;
            out.w_2 = (out.P_in - P_ex)/out.rho_in;
            out.W_dot_in = out.M_dot_in*(out.w_1+out.w_2);
            out.W_dot_loss = Expander.alpha*out.W_dot_in + Expander.W_dot_loss_0 + Expander.C_loss*N_exp/60*2*pi;
            out.W_dot = out.W_dot_in - out.W_dot_loss;
            out.h_ex = h_su - (out.w_1+out.w_2);
            out.T_ex = CoolProp.PropsSI('T','H',out.h_ex,'P',P_ex,fluid);
            out.h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
            out.W_dot_s = out.M_dot*(h_su - out.h_ex_s);
            out.epsilon_s = out.W_dot/out.W_dot_s;
            out.FF = 1;
            out.flag = 1;
            
        case 'SemEmp_rv_mc_lk'
            
            if not(isfield(Expander,{'alpha'}))
                Expander.alpha = 0;
            end
            if not(isfield(Expander,{'W_dot_loss_0'}))
                Expander.W_dot_loss_0 = 0;
            end
            if not(isfield(Expander,{'C_loss'}))
                Expander.C_loss = 0;
            end
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
            rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
            out.M_dot_in = N_exp/60*Expander.V_s*rho_su;
            out.rho_in = rho_su/Expander.r_v_in;
            out.P_in = CoolProp.PropsSI('P','D',out.rho_in,'S',s_su,fluid);
            out.h_in = CoolProp.PropsSI('H','P',out.P_in,'S',s_su,fluid);
            out.w_1 = h_su- out.h_in;
            out.w_2 = (out.P_in - P_ex)/out.rho_in;
            out.W_dot_in = out.M_dot_in*(out.w_1+out.w_2);
            out.W_dot_loss = Expander.alpha*out.W_dot_in + Expander.W_dot_loss_0 + Expander.C_loss*N_exp/60*2*pi;
            out.W_dot = out.W_dot_in - out.W_dot_loss;
            options = optimset('Display','off');
            [out.gamma, out.errgamma, out.flag_gamma] = fsolve (@(x) Nozzle(x, P_su, s_su, rho_su, P_ex, fluid),0.8,options);
            out.P_thr = max(P_ex,P_su*(2/(out.gamma+1))^(out.gamma/(out.gamma-1)));
            out.rho_thr = CoolProp.PropsSI('D','P',out.P_thr,'S',s_su,fluid);
            out.C_thr = sqrt(2*(h_su - CoolProp.PropsSI('H','P',out.P_thr,'S',s_su,fluid)));
            out.M_dot_leak = (Expander.A_leak0+Expander.A_leak_n*P_su/Expander.P_su_n)*out.C_thr*out.rho_thr;
            out.M_dot = out.M_dot_leak + out.M_dot_in;
            out.h_ex1 = h_su - (out.w_1+out.w_2);
            out.h_ex = (out.M_dot_leak*h_su + out.M_dot_in*out.h_ex1)/(out.M_dot);
            out.T_ex = CoolProp.PropsSI('T','H',out.h_ex,'P',P_ex,fluid);
            out.h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
            out.W_dot_s = out.M_dot*(h_su - out.h_ex_s);
            out.epsilon_s = out.W_dot/out.W_dot_s;
            out.FF = out.M_dot/(Expander.V_s*N_exp/60*rho_su);
            out.flag = 1;
            
        case 'SemEmp_rv_mc_lk_dp'
            
            if not(isfield(Expander,{'alpha'}))
                Expander.alpha = 0;
            end
            if not(isfield(Expander,{'W_dot_loss_0'}))
                Expander.W_dot_loss_0 = 0;
            end
            if not(isfield(Expander,{'C_loss'}))
                Expander.C_loss = 0;
            end
            x0 = 1.05*Expander.V_s*N_exp/60*CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
            options = optimset('Display','off', 'TolFun' , eps,  'TolFun', eps);
            [x, ~, flag] = fsolve(@(x) Exp_Rv_Mc_Lk_Dp_res(x, fluid, P_su, h_su, N_exp, Expander.V_s, Expander.r_v_in, P_ex, Expander.A_leak0, Expander.A_leak_n, Expander.d_su, Expander.alpha, Expander.W_dot_loss_0, Expander.P_su_n, Expander.C_loss), x0, options);
            out = Exp_Rv_Mc_Lk_Dp(x, fluid, P_su, h_su, N_exp, Expander.V_s, Expander.r_v_in, P_ex, Expander.A_leak0, Expander.A_leak_n, Expander.d_su, Expander.alpha, Expander.W_dot_loss_0, Expander.P_su_n, Expander.C_loss);
            out.M_dot = x;
            out.A_leak = Expander.A_leak0 + Expander.A_leak_n*(P_su/Expander.P_su_n);
            out.T_ex = CoolProp.PropsSI('T','H',out.h_ex,'P',P_ex,fluid);
            out.FF = out.M_dot/(Expander.V_s*N_exp/60*CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid));
            out.flag = flag;
            
        case 'SemEmp_rv_mc_lk_dp_ht'
            
            if not(isfield(Expander,{'alpha'}))
                Expander.alpha = 0;
            end
            if not(isfield(Expander,{'W_dot_loss_0'}))
                Expander.W_dot_loss_0 = 0;
            end
            if not(isfield(Expander,{'C_loss'}))
                Expander.C_loss = 0;
            end
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            ff_guess = [1 0.8 1.2 0.7 1.3 0.4 1.7];
            stop = 0;
            k =1;
            while not(stop) && k <= length(ff_guess)
                x0(1) = ff_guess(k)*Expander.V_s*N_exp/60*CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
                x0(2) = 0.9*T_su+0.1*T_amb;
                x0(3) = 1;
                ub = 2*x0;
                options = optimset('Display','off');
                [x, ~, flag] = fsolve(@(x)  Exp_Rv_Mc_Lk_Dp_Ht_res(x, ub, fluid, P_su, h_su, N_exp, Expander.V_s, Expander.r_v_in, P_ex, Expander.A_leak0, Expander.A_leak_n, Expander.d_su, Expander.alpha, Expander.W_dot_loss_0, Expander.AU_su_n, Expander.M_dot_n, Expander.P_su_n, Expander.AU_ex_n, Expander.AU_amb, T_amb, Expander.C_loss), x0./ub, options);
                if flag > 0
                    stop = 1;
                end
                k = k + 1;
            end
            x = x.*ub;
            out = Exp_Rv_Mc_Lk_Dp_Ht(x(1), x(2), x(3), fluid, P_su, h_su, N_exp, Expander.V_s, Expander.r_v_in, P_ex, Expander.A_leak0, Expander.A_leak_n, Expander.d_su, Expander.alpha, Expander.W_dot_loss_0, Expander.AU_su_n, Expander.M_dot_n, Expander.P_su_n, Expander.AU_ex_n, Expander.AU_amb, T_amb, Expander.C_loss);
            out.M_dot = x(1);
            out.T_w = x(2);
            out.gamma = x(3);
            out.A_leak = Expander.A_leak0 + Expander.A_leak_n*(P_su/Expander.P_su_n);
            if out.T_w <= T_su
                out.T_ex = CoolProp.PropsSI('T','H',out.h_ex,'P',P_ex,fluid);
            else
                out.T_ex = T_su;
            end
            out.FF = out.M_dot/(Expander.V_s*N_exp/60*CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid));
            if out.resE < 1e-4 && out.resMlk <1e-4
                out.flag = flag;
            else
                out.flag = -1;
            end
            
        otherwise
            warning('Error: model name is not valid');
            return
    end
else
    T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
    rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
    out.M_dot_in = N_exp/60*Expander.V_s*rho_su;
    out.M_dot = out.M_dot_in;
    out.W_dot_in =  0;
    out.W_dot = out.W_dot_in;
    out.Q_dot_amb = 0;
    out.h_ex = h_su;
    out.T_ex = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex,fluid);
    out.FF = 0;
    out.epsilon_s =0;
    out.flag = -2;
end

if out.flag <=0 && Expander.displayResults ==1
    warning('The model did not converge correctly - the solution may be wrong')
end
out.time = toc(tstart_exp);

%% TS DIAGRAM and DISPLAY

s_su = CoolProp.PropsSI('S','H',h_su,'P',P_su, fluid);
s_ex = CoolProp.PropsSI('S','H', out.h_ex,'P', P_ex, fluid);
TS.T = [T_su out.T_ex];
TS.s = [s_su s_ex];

if Expander.displayResults ==1
    in.fluid = fluid;
    in.T_su = T_su;
    in.h_su = h_su;
    in.P_su = P_su;
    in.P_ex = P_ex;
    in.P_ex = P_ex;
    in.N_exp = N_exp;
    in.T_amb = T_amb;
    in.modelType= Expander.modelType;
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
end

end

%% FUNCTIONS

function res = Nozzle(x, P_su, s_su, rho_su, P_ex, fluid)
P_thr = max(P_ex,P_su*(2/(x+1))^(x/(x-1)));
gamma_bis =  log10(P_su/P_thr)/log10(rho_su/CoolProp.PropsSI('D','P',P_thr,'S',s_su,fluid));
res = (gamma_bis-x)/x;
end

function res = Exp_Rv_Mc_Lk_Dp_Ht_res(x, ub, fluid, P_su, h_su, N_exp, V_s, r_v_in, P_ex, A_leak0, A_leak_n, d_su, alpha, W_dot_loss_0, AU_su_n, M_dot_n, P_su_n, AU_ex_n, AU_amb, T_amb, C_loss)
x = x.*ub;
out = Exp_Rv_Mc_Lk_Dp_Ht(x(1), x(2), x(3), fluid, P_su, h_su, N_exp, V_s, r_v_in, P_ex, A_leak0, A_leak_n, d_su, alpha, W_dot_loss_0, AU_su_n, M_dot_n, P_su_n, AU_ex_n, AU_amb, T_amb, C_loss);
res = [out.resE ; out.resMlk; out.resgamma];
end

function out = Exp_Rv_Mc_Lk_Dp_Ht(M_dot, T_w, gamma, fluid, P_su, h_su, N_exp, V_s, r_v_in, P_ex, A_leak0, A_leak_n, d_su, alpha, W_dot_loss_0, AU_su_n, M_dot_n, P_su_n, AU_ex_n, AU_amb, T_amb, C_loss)
if isnan(T_w)
    T_w =  CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
end
M_dot = max(1e-5, M_dot);
h_min = CoolProp.PropsSI('H','P',CoolProp.PropsSI('pmin','Q',0,'P',1e5,fluid)+1,'T',CoolProp.PropsSI('Tmin','Q',0,'P',1e5,fluid)+1,fluid);
h_max = CoolProp.PropsSI('H','P',CoolProp.PropsSI('pmax','Q',0,'P',1e5,fluid)-1,'T',CoolProp.PropsSI('Tmax','Q',0,'P',1e5,fluid)-1,fluid);

s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
out.h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
AU_su1 = AU_su_n*(M_dot/M_dot_n)^0.8;
AU_ex1 = AU_ex_n*(M_dot/M_dot_n)^0.8;
out.h_su1 = h_su;
out.P_su1 = max(CoolProp.PropsSI('P','S',s_su,'H',max(h_su - (M_dot/(pi*d_su^2/4*rho_su))^2/2, out.h_ex_s), fluid),P_ex+1);
out.T_su1 = CoolProp.PropsSI('T','P',out.P_su1,'H',out.h_su1,fluid);
if CoolProp.PropsSI('Q','P',out.P_su1,'H',out.h_su1,fluid) == -1
    out.cp_su1 = CoolProp.PropsSI('C','P',out.P_su1,'H',out.h_su1,fluid);
else
    out.cp_su1 = CoolProp.PropsSI('C','P',out.P_su1,'Q',0,fluid);
end
out.epsilon_su1 = max(0,(1-exp(-AU_su1/(M_dot*out.cp_su1))));
out.DP_su = P_su-out.P_su1;
out.Q_dot_su = max(0,out.epsilon_su1*M_dot*out.cp_su1*(out.T_su1 - T_w));
out.h_su2 = min(h_max,max(out.h_ex_s,out.h_su1 - out.Q_dot_su/M_dot));
out.s_su2 = CoolProp.PropsSI('S','H',out.h_su2,'P',out.P_su1,fluid);
out.rho_su2 = CoolProp.PropsSI('D','H',out.h_su2,'P',out.P_su1,fluid);
out.M_dot_in = N_exp/60*V_s*out.rho_su2;
out.rho_in = out.rho_su2/r_v_in;
out.h_in = CoolProp.PropsSI('H','D',out.rho_in,'S',out.s_su2,fluid);
out.P_in = CoolProp.PropsSI('P','D',out.rho_in,'H',out.h_in,fluid);
out.w_1 = out.h_su2-out.h_in;
out.w_2 = (out.P_in - P_ex)/out.rho_in;
out.h_ex2 = out.h_in - out.w_2;
out.W_dot_in = out.M_dot_in*(out.w_1+out.w_2);
out.W_dot_loss = alpha*out.W_dot_in + W_dot_loss_0 + C_loss*N_exp/60*2*pi;
out.W_dot = out.W_dot_in - out.W_dot_loss;
out.W_dot_s = M_dot*(h_su - out.h_ex_s);
out.epsilon_s = out.W_dot/out.W_dot_s;
out.M_dot_leak = M_dot - out.M_dot_in;
out.gamma = max(0.1,min(gamma,3));
out.P_thr = max(P_ex,out.P_su1*(2/(out.gamma+1))^(out.gamma/(out.gamma-1)));
out.rho_thr = CoolProp.PropsSI('D','P',out.P_thr,'S',out.s_su2,fluid);
out.gamma_bis =  log10(out.P_su1/out.P_thr)/log10(out.rho_su2/out.rho_thr);

out.resgamma = (out.gamma_bis-out.gamma)/out.gamma;


out.C_thr = sqrt(2*(out.h_su2 - CoolProp.PropsSI('H','P',out.P_thr,'S',out.s_su2,fluid)));
out.M_dot_leak_bis = (A_leak0+A_leak_n*P_su/P_su_n)*out.C_thr*out.rho_thr;
out.h_ex1 = max(min((out.M_dot_in*out.h_ex2 + out.M_dot_leak*out.h_su2)/M_dot, out.h_su2), out.h_ex2);
out.T_ex1 = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex1,fluid);
out.cp_ex1 = CoolProp.PropsSI('C','P',P_ex,'H',out.h_ex1,fluid);
out.epsilon_ex1 = max(0,(1-exp(-AU_ex1/(M_dot*out.cp_ex1))));
out.Q_dot_ex = max(0,out.epsilon_ex1*M_dot*out.cp_ex1*(T_w-out.T_ex1));
out.h_ex = min(out.h_ex1 + out.Q_dot_ex/M_dot,h_max);
out.Q_dot_amb = AU_amb*(T_w-T_amb);
out.resE = abs((out.Q_dot_su + out.W_dot_loss - out.Q_dot_ex - out.Q_dot_amb)/(out.Q_dot_su + out.W_dot_loss));
if out.M_dot_leak_bis == 0 && out.M_dot_leak == 0
    out.resMlk = 0;
elseif out.M_dot_leak ~= 0 && - out.M_dot_leak_bis == 0;
    out.resMlk = 1;
elseif out.M_dot_leak == 0 && - out.M_dot_leak_bis ~= 0;
    out.resMlk = 1;
else
    out.resMlk = abs((out.M_dot_leak - out.M_dot_leak_bis)/out.M_dot_leak_bis);
end

end

function res = Exp_Rv_Mc_Lk_Dp_res(x, fluid, P_su, h_su, N_exp, V_s, r_v_in, P_ex, A_leak0, A_leak_n, d_su, alpha, W_dot_loss_0, P_su_n, C_loss)

out = Exp_Rv_Mc_Lk_Dp(x, fluid, P_su, h_su, N_exp, V_s, r_v_in, P_ex, A_leak0, A_leak_n, d_su, alpha, W_dot_loss_0, P_su_n, C_loss);
res = out.resMlk;

end

function out = Exp_Rv_Mc_Lk_Dp(M_dot, fluid, P_su, h_su, N_exp, V_s, r_v_in, P_ex, A_leak0, A_leak_n, d_su, alpha, W_dot_loss_0, P_su_n, C_loss)

s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
out.h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
out.h_su1 = h_su;
out.P_su1 = CoolProp.PropsSI('P','S',s_su,'H',h_su - (M_dot/(pi*d_su^2/4*rho_su))^2/2,fluid);
out.DP_su = P_su-out.P_su1;
out.h_su2 = out.h_su1;
out.s_su2 = CoolProp.PropsSI('S','H',out.h_su2,'P',out.P_su1,fluid);
out.rho_su2 = CoolProp.PropsSI('D','H',out.h_su2,'P',out.P_su1,fluid);
out.M_dot_in = N_exp/60*V_s*out.rho_su2;
out.rho_in = out.rho_su2/r_v_in;
out.P_in = CoolProp.PropsSI('P','D',out.rho_in,'S',out.s_su2,fluid);
out.h_in = CoolProp.PropsSI('H','P',out.P_in,'S',out.s_su2,fluid);
out.w_1 = out.h_su2-out.h_in;
out.w_2 = (out.P_in - P_ex)/out.rho_in;
out.h_ex2 = out.h_in - out.w_2;
out.W_dot_in = out.M_dot_in*(out.w_1+out.w_2);
out.W_dot_loss = alpha*out.W_dot_in + W_dot_loss_0 + C_loss*N_exp/60*2*pi;
out.W_dot = out.W_dot_in - out.W_dot_loss;
out.W_dot_s = M_dot*(h_su - out.h_ex_s);
out.epsilon_s = out.W_dot/out.W_dot_s;
out.M_dot_leak = M_dot - out.M_dot_in;
fgamma=@(x) Nozzle(x, out.P_su1, out.s_su2, out.rho_su2, P_ex, fluid);
out.gamma = zeroBrent( 0.1, 2, 1e-6, 1e-6, fgamma );out.P_thr = max(P_ex,out.P_su1*(2/(out.gamma+1))^(out.gamma/(out.gamma-1)));
out.rho_thr = CoolProp.PropsSI('D','P',out.P_thr,'S',out.s_su2,fluid);
out.C_thr = sqrt(2*(out.h_su2 - CoolProp.PropsSI('H','P',out.P_thr,'S',out.s_su2,fluid)));
out.M_dot_leak_bis = (A_leak0+A_leak_n*P_su/P_su_n)*out.C_thr*out.rho_thr;
out.h_ex1 = (out.M_dot_in*out.h_ex2 + out.M_dot_leak*out.h_su2)/M_dot;
out.h_ex = out.h_ex1 ;
out.resMlk = abs((out.M_dot_leak - out.M_dot_leak_bis)/out.M_dot_leak_bis);

end

