function [out,TS] = ExpanderModel5(fluid, P_su, h_su, N_exp, P_ex, T_amb, param)
% fluid, P_su, h_su, N_exp, P_ex, T_amb, param
%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 27/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% ExpanderModel is a single matlab code implementing three different modelling
% approaches to simulate a volumetric expander (see the Documentation/ExpanderModel_MatlabDoc)
%
% The model inputs are:
%       - P_su: inlet pressure of the WF                          	[Pa]
%       - h_su: inlet temperature of the WF                        	[J/kg]
%       - P_ex: outlet pressure of the WF                          	[Pa]
%       - fluid: nature of the fluid (string)                       [-]
%       - N_exp: expander rotational speed                          [rpm]
%       - T_amb : ambien temperature                                [K]
%       - param: structure variable containing the model parameters
%
% The model paramters provided in 'param' depends of the type of model selected:
%       - if param.modelType = 'CstEff':
%           param.V_s , machine displacement volume                	[m3]
%           param.V, volume of the expander                         [m^3]
%           param.epsilon_is, isentropic efficiency                	[-]
%           param.FF, filling factor (volumetric efficiency)       	[-]
%           param.displayResults, flag to display the results or not[1/0]
%
%       - if param.modelType = 'PolEff':
%           param.V_s , machine displacement volume                	[m3]
%           param.V, volume of the expander                         [m^3]
%           param.coeffPol_is, polynmial coef for epsilon_is        [-]
%           param.coeffPol_ff, polynmial coef for FF      [-]
%           param.displayResults, flag to display the results or not[1/0]
%
%       - if param.modelType = 'SemiEmp':
%           param.V_s,  machine displacement volume                 [m^3];
%           param.V, volume of the expander                         [m^3]
%           param.alpha, proportional losses coefficient            [-];
%           param.W_dot_loss_0, constant losses term                [W];
%           param.C_loss, losses torque                             [Nm];
%           param.r_v_in, built-in volumetric ration                [-];
%           param.A_leak0, nozzle cross section area for the leakage[m²];
%           param.d_su, nozzle diameter for supply pressure drop    [m];
%           param.AU_su_n, global heat transfer coefficient for supply heat transfer    [W/K];
%           param.AU_ex_n, global heat transfer coefficient for exhaust heat transfer	[W/K];
%           param.AU_amb, global heat transfer coefficient for ambient heat losses      [W/K];
%           param.M_dot_n, nominal mass flow rate                   [kg/s];

%
% The model outputs are:
%       - out: a structure variable which includes
%               - T_ex =  exhaust temperature                    [K]
%               - h_ex =  exhaust enthalpy                       [J/kg]
%               - M_dot = fluid mass flow rate                   [kg/s]
%               - W_dot = mechanical power                       [W]
%               - Q_dot_amb = ambiant losses                     [W]
%               - epsilon_is = isentropic efficiency             [-]
%               - FF = filling factor (volumetric efficiency)    [-]
%               - time = the code computational time             [sec]
%               - flag = simulation flag                         [-1/1]
%               - M = mass of refrigerent inside the expaner     [kg]
%
%       - TS : a stucture variable which contains the vectors of temperature
%              and entropy of the fluid (useful to generate a Ts diagram 
%              when modelling the entire ORC system 
%
% See the documentation for further details or contact rdickes@ulg.ac.be


%% DEMONSTRATION CASE
if nargin == 0
    
    % Define a demonstration case if ExpanderModel.mat is not executed externally
    fluid = 'R245fa';                           %Nature of the fluid
    N_exp = 5000;                               %Rotational speed       [rpm]
    P_su = 16.753498330038136e+05;               %Supply pressure        [Pa]
    h_su = 4.052843743508205e+05;               %Supply enthalpy        [J/kg]
    P_ex = 2.471310061849047e+05;               %Exhaust pressure       [Pa]
    T_amb = 298.1500;              %Ambient temperature    [K]
    param.displayResults = 1;                   %Flag to control the resustl display [0/1]
    param.modelType = 'SemiEmp';                %Type of model          [CstEff, PolEff, SemiEmp]
    
    switch param.modelType
        case 'CstEff'
        % Example of paramters for modelType = CstEff
            param.epsilon_is = 0.7;
            param.FF = 1.2;
            param.AU_amb = 0.674005126953125;
            param.V_s = 1.279908675799087e-05;
            param.V = 1.492257e-3;
        case 'PolEff'
        % Example of paramters for modelType = PolEff
            param.coeffPol_is  = [1.504191277585839   0.051680226738587  -0.000001018959447  -0.020713559079583   0.000000031728135   0.000000000000223];
            param.coeffPol_ff = [1.122789891166860  -0.050184374872146  -0.000000109428161   0.001845567200683   0.000000026842863   0.000000000000043];
            param.AU_amb = 0.674005126953125;
            param.V_s = 1.279908675799087e-05;
            param.V = 1.492257e-3;
        case 'SemiEmp'
        % Example of paramters for modelType = SemiEmp
            param.V = 1.492257e-3;
            param.V_s = 1.279908675799087e-05;
            param.r_v_in = 2.19;
            param.A_leak0 = 1.344675598144531e-06;
            param.d_su = 0.003227564086914;
            param.alpha = 1.253662109374984e-05;
            param.W_dot_loss_0 = 0;
            param.C_loss = 7.952880859375330e-07;
            param.AU_su_n = 50.033569335937500;
            param.M_dot_n =  0.068378356905196;
            param.P_su_n = 1.417898896236315e+06;
            param.AU_ex_n = 94.017028808593750;
            param.AU_amb = 0.674005126953125;
            load('C:\Users\RDickes\Google Drive\PhD\MOR study\ORC\Experimental database\Sun2Power\OffDesign\gamma_R245fa.mat');
            param.gamma.gamma_PQ_pol = gamma_PQ_R245fa; param.gamma.gamma_PT_pol = gamma_PT_R245fa;

    end
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

if not(isfield(param,{'V'}))
    param.V = 0;
end

if P_su > P_ex && h_su > CoolProp.PropsSI('H','P',P_su,'Q',0.1,fluid);


    %If the external conditions are viable (i.e. pressure ratio higher than
    %1 and in two phase conditions), we proceed to the modeling:
             
    switch param.modelType
    %Select the proper model paradigm chosen by the user            
        case 'CstEff'
            M_dot = N_exp/60*param.V_s*param.FF*rho_su;
            W_dot =  M_dot*(h_su-h_ex_s)*param.epsilon_is;
            Q_dot_amb = max(0,param.AU_amb*(T_su - T_amb));
            FF = param.FF;
            epsilon_is = param.epsilon_is;           
            h_ex = h_su - (W_dot+Q_dot_amb)/M_dot;
            if h_ex > param.h_min && h_ex < param.h_max
                out.flag = 1;
            else
                out.flag = -1;
            end            
            
        case 'PolEff'
            if length(param.coeffPol_ff) == 6
                % if expander with constant speed -> quadratic function of Rp and rho_su
                FF = max(0.5, min(param.coeffPol_ff(1) + param.coeffPol_ff(2)*(P_su/P_ex) + param.coeffPol_ff(3)*rho_su + param.coeffPol_ff(4)*(P_su/P_ex)^2 + param.coeffPol_ff(5)*(P_su/P_ex)*rho_su + param.coeffPol_ff(6)*(rho_su)^2, 5));
            elseif length(param.coeffPol_ff) == 10
                % if expander with variable speed -> quadratic function of Rp, rho_su and N_exp
                FF = max(0.3,min(param.coeffPol_ff(1) + param.coeffPol_ff(2)*(P_su/P_ex) + param.coeffPol_ff(3)*rho_su + param.coeffPol_ff(4)*(P_su/P_ex)^2 + param.coeffPol_ff(5)*(P_su/P_ex)*rho_su + param.coeffPol_ff(6)*(rho_su)^2 + param.coeffPol_ff(7)*N_exp + param.coeffPol_ff(8)*N_exp^2 + param.coeffPol_ff(9)*N_exp*(P_su/P_ex) + param.coeffPol_ff(10)*N_exp*rho_su,5));
            end
            if length(param.coeffPol_is) == 6
                % if expander with constant speed -> quadratic function of Rp and rho_su
                epsilon_is = max(0.01,min(param.coeffPol_is(1) + param.coeffPol_is(2)*(P_su/P_ex) + param.coeffPol_is(3)*rho_su + param.coeffPol_is(4)*(P_su/P_ex)^2 + param.coeffPol_is(5)*(P_su/P_ex)*rho_su + param.coeffPol_is(6)*(rho_su)^2,1));
            elseif length(param.coeffPol_is) == 10
                % if expander with variable speed -> quadratic function of Rp, rho_su and N_exp
                epsilon_is = max(0.01,min(param.coeffPol_is(1) + param.coeffPol_is(2)*(P_su/P_ex) + param.coeffPol_is(3)*rho_su + param.coeffPol_is(4)*(P_su/P_ex)^2 + param.coeffPol_is(5)*(P_su/P_ex)*rho_su + param.coeffPol_is(6)*(rho_su)^2 + param.coeffPol_is(7)*N_exp + param.coeffPol_is(8)*N_exp^2 + param.coeffPol_is(9)*N_exp*(P_su/P_ex) + param.coeffPol_is(10)*N_exp*rho_su,1));
            end
            M_dot = N_exp/60*param.V_s*FF*rho_su;
            W_dot = M_dot*(h_su - h_ex_s)*epsilon_is;
            Q_dot_amb = max(0,param.AU_amb*(T_su - T_amb));
            h_ex = h_su - (W_dot + Q_dot_amb)/M_dot;
            if h_ex > param.h_min && h_ex < param.h_max
                out.flag = 1;
            else
                out.flag = -1;
            end   
                      
        case 'SemiEmp'
            % Default values assigned to the parameters if not specified
            if not(isfield(param,{'alpha'}))
                param.alpha = 0;
            end
            if not(isfield(param,{'W_dot_loss_0'}))
                param.W_dot_loss_0 = 0;
            end
            if not(isfield(param,{'C_loss'}))
                param.C_loss = 0;
            end
            if not(isfield(param,{'d_su'}))
                param.d_su = 1e2;
            end
            if not(isfield(param,{'A_leak0'}))
                param.A_leak0 = 0;
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
%             ff_guess = [1 0.8 1.2 0.7 1.3 0.4 1.7 3 ]; %guesses on the filling factor to provide suitable initial point for the iteration
%             stop = 0;
%             k =1;
%             while not(stop) && k <= length(ff_guess)
%             % Loop to permit multiple attempts to solve the implicit
%             % calculation of Exp_SemiEmp trough Exp_SemiEmp_res
%                 x0(1) = ff_guess(k)*param.V_s*N_exp/60*CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid); %initial value for M_dot
%                 x0(2) = 0.9*T_su+0.1*T_amb; %initial value for T_wall
%                 x0(3) = 1; %initial value for gamma
%                 ub = 2*x0; % upper bound for fsolve
%                 options = optimset('Display','off');
%                 [x, ~, flag] = fsolve(@(x)  Exp_SemiEmp_res(x, ub, fluid, P_su, h_su, N_exp, param.V_s, param.r_v_in, P_ex, param.A_leak0, param.d_su, param.alpha, param.W_dot_loss_0, param.AU_su_n, param.M_dot_n, param.AU_ex_n, param.AU_amb, T_amb, param.C_loss, param.h_min, param.h_max), x0./ub, options);
%                 if flag > 0
%                     stop = 1;
%                 end
%                 k = k + 1;
%             end
            
            ff_guess = [0.7 1.2 0.8 1.3 0.4 1.7 3 ]; %guesses on the filling factor to provide suitable initial point for the iteration
            x_T_guess = [0.7 0.95 0.8 0.99 0.9  ];
            stop = 0;
            
            j = 1;
            while not(stop) && j <= length(x_T_guess)
                k =1;
                while not(stop) && k <= length(ff_guess)
                        
                    % Loop to permit multiple attempts to solve the implicit
                    % calculation of Exp_SemiEmp trough Exp_SemiEmp_res
                    x0(1) = ff_guess(k)*param.V_s*N_exp/60*CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid); %initial value for M_dot
                    x0(2) = x_T_guess(j)*T_su+(1-x_T_guess(j))*T_amb; %initial value for T_wall
                    ub = 2*x0; % upper bound for fsolve
                    options = optimset('Display','none');% 'TolX', 1e-8, 'TolFun', 1e-4);
                    [x, res, flag] = fsolve(@(x)  Exp_SemiEmp_res(x, ub, fluid, P_su, h_su, N_exp, param.V_s, param.r_v_in, P_ex, param.A_leak0, param.d_su, param.alpha, param.W_dot_loss_0, param.AU_su_n, param.M_dot_n, param.AU_ex_n, param.AU_amb, T_amb, param.C_loss, param.h_min, param.h_max, param), x0./ub, options);
                    if norm(res) < 1e-4 %flag > 0 
                        stop = 1;
                    end
                    k = k + 1;
                end
                j = j + 1;
            end
            
            x = x.*ub;
            int = Exp_SemiEmp(x(1), x(2), fluid, P_su, h_su, N_exp, param.V_s, param.r_v_in, P_ex, param.A_leak0, param.d_su, param.alpha, param.W_dot_loss_0, param.AU_su_n, param.M_dot_n, param.AU_ex_n, param.AU_amb, T_amb, param.C_loss, param.h_min, param.h_max, param);
            M_dot = x(1);
            T_w = int.T_w;
            W_dot = int.W_dot;
            h_ex = int.h_ex;
            epsilon_is = int.epsilon_is;
            FF = M_dot/(param.V_s*N_exp/60*CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid));
            Q_dot_amb = int.Q_dot_amb;
            if int.resE < 1e-4 && int.resMlk <1e-4 && int.resgamma <1e-4 && h_ex > param.h_min && h_ex < param.h_max
                out.flag = flag;
            else
                out.flag = -1;
            end
            out.complete_results = int;
        otherwise
            warning('Error: model name is not valid');
            return
    end
else
    out.flag = -2;
end

if out.flag > 0;
    out.T_w = T_w;
    out.h_ex = h_ex;
    out.M_dot = M_dot;
    out.W_dot = W_dot;
    out.Q_dot_amb = Q_dot_amb;
    out.epsilon_is = epsilon_is;
    out.FF = FF;
    out.T_ex = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex,fluid);
    out.M = (CoolProp.PropsSI('D','H',h_su,'P',P_su,fluid)+CoolProp.PropsSI('D','H',out.h_ex,'P',P_ex,fluid))/2*param.V;    
else    
    out.T_w = NaN;
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

% Generate the output variable TS 
s_su = CoolProp.PropsSI('S','H',h_su,'P',P_su, fluid);
s_ex = CoolProp.PropsSI('S','H', out.h_ex,'P', P_ex, fluid);
TS.T = [T_su out.T_ex];
TS.s = [s_su s_ex];

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
    in.modelType= param.modelType;
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

%% NESTED FUNCTIONS
function res = Exp_SemiEmp_res(x, ub, fluid, P_su, h_su, N_exp, V_s, r_v_in, P_ex, A_leak0, d_su, alpha, W_dot_loss_0, AU_su_n, M_dot_n, AU_ex_n, AU_amb, T_amb, C_loss, h_min, h_max, param)
x = x.*ub;
out = Exp_SemiEmp(x(1), x(2), fluid, P_su, h_su, N_exp, V_s, r_v_in, P_ex, A_leak0, d_su, alpha, W_dot_loss_0, AU_su_n, M_dot_n, AU_ex_n, AU_amb, T_amb, C_loss, h_min, h_max, param);
res = [out.resE ; out.resMlk];
end

function out = Exp_SemiEmp(M_dot, T_w, fluid, P_su, h_su, N_exp, V_s, r_v_in, P_ex, A_leak0, d_su, alpha, W_dot_loss_0, AU_su_n, M_dot_n, AU_ex_n, AU_amb, T_amb, C_loss, h_min, h_max, param)
if isnan(T_w)
    T_w =  CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
end
M_dot = max(1e-5, M_dot);

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
out.h_su2 = min(h_max,max(max(out.h_ex_s, CoolProp.PropsSI('H','Q',0.1,'P',out.P_su1,fluid)),out.h_su1 - out.Q_dot_su/M_dot));
out.s_su2 = CoolProp.PropsSI('S','H',out.h_su2,'P',out.P_su1,fluid);
out.rho_su2 = CoolProp.PropsSI('D','H',out.h_su2,'P',out.P_su1,fluid);
out.M_dot_in = N_exp/60*V_s*out.rho_su2;
out.rho_in = out.rho_su2/r_v_in;
try
    out.P_in = CoolProp.PropsSI('P','D',out.rho_in,'S',out.s_su2,fluid);
catch
    delta = 0.0001;
    out.P_in = 0.5*CoolProp.PropsSI('P','D',out.rho_in*(1+delta),'S',out.s_su2,fluid)+0.5*CoolProp.PropsSI('P','D',out.rho_in*(1-delta),'S',out.s_su2,fluid);

end
out.h_in = CoolProp.PropsSI('H','D',out.rho_in,'P',out.P_in,fluid);
out.w_1 = out.h_su2-out.h_in;
out.w_2 = (out.P_in - P_ex)/out.rho_in;
out.h_ex2 = out.h_in - out.w_2;
out.W_dot_in = out.M_dot_in*(out.w_1+out.w_2);
out.W_dot_loss = alpha*out.W_dot_in + W_dot_loss_0 + C_loss*N_exp/60*2*pi;
out.W_dot = out.W_dot_in - out.W_dot_loss;
out.W_dot_s = M_dot*(h_su - out.h_ex_s);
out.epsilon_is = out.W_dot/out.W_dot_s;
out.M_dot_leak = M_dot - out.M_dot_in;

out.Q_su2 = CoolProp.PropsSI('Q','P',out.P_su1,'H',out.h_su2,fluid);
if  out.Q_su2 < 0
    out.gamma  = feval(param.gamma.gamma_PT_pol,[out.P_su1/1e5, CoolProp.PropsSI('T','P',out.P_su1,'H',out.h_su2,fluid)/1e2]);
else
    out.gamma = feval(param.gamma.gamma_PQ_pol,[out.P_su1/1e5, out.Q_su2]);
end
out.gamma = max(1e-2,out.gamma);
P_crit = out.P_su1*(2/(out.gamma+1))^(out.gamma/(out.gamma-1));
out.resgamma = 0;
out.P_thr = max(P_ex,P_crit);
out.rho_thr = CoolProp.PropsSI('D','P',out.P_thr,'S',out.s_su2,fluid);
out.C_thr = sqrt(2*(out.h_su2 - CoolProp.PropsSI('H','P',out.P_thr,'S',out.s_su2,fluid)));
out.M_dot_leak_bis = A_leak0*out.C_thr*out.rho_thr;

% out.gamma = max(0.1,min(gamma,3));
% out.P_thr = max(P_ex,out.P_su1*(2/(out.gamma+1))^(out.gamma/(out.gamma-1)));
% out.rho_thr = CoolProp.PropsSI('D','P',out.P_thr,'S',out.s_su2,fluid);
% out.gamma_bis =  log10(out.P_su1/out.P_thr)/log10(out.rho_su2/out.rho_thr);
% out.resgamma = (out.gamma_bis-out.gamma)/out.gamma;
% out.C_thr = sqrt(2*(out.h_su2 - CoolProp.PropsSI('H','P',out.P_thr,'S',out.s_su2,fluid)));
% out.M_dot_leak_bis = A_leak0*out.C_thr*out.rho_thr;

out.h_ex1 = max(min((out.M_dot_in*out.h_ex2 + out.M_dot_leak*out.h_su2)/M_dot, out.h_su2), out.h_ex2);
out.T_ex1 = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex1,fluid);
out.cp_ex1 = CoolProp.PropsSI('C','P',P_ex,'H',out.h_ex1,fluid);
out.epsilon_ex1 = max(0,(1-exp(-AU_ex1/(M_dot*out.cp_ex1))));
out.Q_dot_ex = max(0,out.epsilon_ex1*M_dot*out.cp_ex1*(T_w-out.T_ex1));
out.h_ex = min(out.h_ex1 + out.Q_dot_ex/M_dot,h_max);
out.Q_dot_amb = AU_amb*(T_w-T_amb);
out.resE =  abs((out.Q_dot_su + out.W_dot_loss - out.Q_dot_ex - out.Q_dot_amb)/(out.Q_dot_su + out.W_dot_loss));
% out.resE = abs((out.Q_dot_su  - out.Q_dot_ex - out.Q_dot_amb)/(out.Q_dot_ex+out.Q_dot_amb));

if out.M_dot_leak_bis == 0 && out.M_dot_leak == 0
    out.resMlk = 0;
elseif out.M_dot_leak ~= 0 && - out.M_dot_leak_bis == 0;
    out.resMlk = 1;
elseif out.M_dot_leak == 0 && - out.M_dot_leak_bis ~= 0;
    out.resMlk = 1;
else
    out.resMlk = abs((out.M_dot_leak - out.M_dot_leak_bis)/out.M_dot_leak_bis);
end
out.T_w = T_w;

%out.resMlk = 0;
out.resE = 0;
end

