function [out,TS] = PumpModel3_Solub(P_su, h_su, P_ex, fluid, fluid_lub, C_oil, N_pp, T_amb, param)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 27/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% PumpModel is a single matlab code implementing three different modelling
% approaches to simulate a volumetric pump (see the Documentation/PumpModel_MatlabDoc)
%
% The model inputs are:
%       - P_su: inlet pressure of the WF                          	[Pa]
%       - h_su: inlet temperature of the WF                        	[J/kg]
%       - P_ex: outlet pressure of the WF                          	[Pa]
%       - fluid: nature of the fluid (string)                       [-]
%       - N_pp: the pump speed                                    	[rpm]
%       - param: structure variable containing the model parameters
%
% The model paramters provided in 'param' depends of the type of model selected:
%       - if param.modelType = 'CstEff':
%           param.V_s , machine displacement volume                	[m3]
%           param.V, volume of the pump                             [m^3]
%           param.epsilon_is, isentropic efficiency                	[-]
%           param.epsilon_vol, volumetric efficiency               	[-]
%           param.displayResults, flag to display the results or not[1/0]
%
%       - if param.modelType = 'PolEff':
%           param.V_s , machine displacement volume                	[m3]
%           param.V, volume of the pump                             [m^3]
%           param.N_pp_nom, pump nominal shaft speed              	[rpm]
%           param.coeffPol_is, polynmial coef for epsilon_is        [-]
%           param.coeffPol_vol, polynmial coef for epsilon_vol      [-]
%           param.displayResults, flag to display the results or not[1/0]
%
%       - if param.modelType = 'SemiEmp':
%           param.V_s , machine displacement volume               	[m3]
%           param.V, volume of the pump                             [m^3]
%           param.A_leak, leakage surface area                     	[m2]
%           param.W_dot_loss, constant power losses                	[W]
%           param.K_0_loss, term for the proportional losses       	[-]
%           param.displayResults, flag to display the results or not [1/0]
%
% The model outputs are:
%       - out: a structure variable which includes
%               - T_ex =  exhaust temperature                    [K]
%               - h_ex =  exhaust enthalpy                       [J/kg]
%               - m_dot = fluid mass flow rate                   [kg/s]
%               - W_dot = mechanical power                       [W]
%               - epsilon_is = isentropic efficiency             [-]
%               - epsilon_vol = volumetric efficiency            [-]
%               - M = mass of fluid inside the pump              [kg]
%               - time = the code computational time             [sec]
%               - flag = simulation flag                         [-1/1]
%
%       - TS : a stucture variable which contains the vectors of temperature
%              and entropy of the fluid (useful to generate a Ts diagram 
%              when modelling the entire ORC system 
%
% See the documentation for further details or contact rdickes@ulg.ac.be


%% DEMONSTRATION CASE

if nargin == 0
    
    % Define a demonstration case if PumpModel.mat is not executed externally
    fluid = 'R134a';               %Nature of the fluid
    P_su = 4.0001e5;                %Supply pressure        [Pa]
    P_ex = 3.6510e+06*3;         %Exhaust pressure       [Pa]
    h_su = 2.6676e+05;              %Supply enthalpy        [J/kg]
    N_pp = 1500;                    %Rotational speed       [rpm]
    param.modelType = 'SemiEmp';    %Type of model          [CstEff, PolEff, SemiEmp]
    param.displayResults = 1;       %Flag to control the resustl display [0/1]
    T_amb = 298;
    param.V_s = 1e-6;               %Machine swepts volume  [m^3]
    param.V =1.4e-3;                %Volume inside the pump
    param.A_leak = 1e-7;
    param.W_dot_0_loss = 200;
    param.K_0_loss  = 0.2;
    param.AU = 1;
end

tstart_pp = tic;                    %Start to evaluate the simulation time

%% PUMP MODELING
% Modelling section of the code
if not(isfield(param, 'displayResults'))
    param.displayResults = 0;
    %if nothing specified by the user, the results are not displayed by
    %default.
end

%if not(isfield(param,'h_min'))
%    param.h_min =  CoolProp.PropsSI('H','P',4e6,'T',253.15,fluid);
%end
%if not(isfield(param,'h_max'))
%    param.h_max =  CoolProp.PropsSI('H','P',1e4,'T',500,fluid);
%end
%if not(isfield(param,'V'))
%    param.V =  0;
%end
%if not(isfield(param,'AU'))
%    param.AU =  0;
%end

if C_oil > 0
    Tsat_su = CoolProp.PropsSI('T', 'Q', 0, 'P', P_su, fluid);
    Tbubble_min_su = R245fa_POE_Tbubble(1-C_oil, P_su, Tsat_su);
    Tsat_ex = CoolProp.PropsSI('T', 'Q', 0, 'P', P_ex, fluid);
    Tbubble_min_ex = R245fa_POE_Tbubble(1-C_oil, P_ex, Tsat_ex);
    [T_su, ~, ~, ~, zeta_r_su, ~, ~, ~] = HP_solubMixt(C_oil, P_su, h_su, fluid, fluid_lub, Tbubble_min_su, Tsat_su, param.fit_DTP_zeta);
    s_su = CoolProp.PropsSI('S','P',P_su,'T',T_su,fluid);
    T_is_ex = CoolProp.PropsSI('T','P',P_ex,'S',s_su,fluid);
    [h_ex_s, ~, ~, ~, ~, ~, ~, ~] = TP_solubMixt(T_is_ex, P_ex, C_oil, fluid, fluid_lub, Tbubble_min_ex, Tsat_ex, param.fit_DTP_zeta);
    [~, ~, ~, rho_su, ~] = R245fa_POE_density(T_su, P_su, zeta_r_su, fluid, fluid_lub, param.fit_ratio_rho, Tbubble_min_su, Tsat_su);
    
else
    T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
    s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
    rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
    h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
end

if P_su < P_ex && N_pp > 0
    %If the external conditions are viable, we proceed to the modeling
    
    
    V_s = param.V_s;
    A_leak = param.A_leak;
    W_dot_loss = param.W_dot_0_loss;
    K_0_loss = param.K_0_loss;
    if isfield(param,'fit_NPSHr')
        NPSH_r = param.fit_NPSHr(N_pp);
        NPSH_a = (P_su-CoolProp.PropsSI('P', 'Q', 0, 'T', T_su, fluid))/rho_su/9.81;
        eta_cavitation =  max(0.1, min(1,1-exp(log(0.03)*(NPSH_a/NPSH_r))));
    else
        eta_cavitation = 1;
    end
    m_dot = max(1e-3, eta_cavitation*((N_pp/60*V_s*rho_su)-(A_leak*sqrt(2*rho_su*(P_ex-P_su)))));
    epsilon_vol = m_dot/(N_pp/60*V_s*rho_su);
    W_dot = W_dot_loss + K_0_loss*m_dot/rho_su*(P_ex-P_su);
    epsilon_is = (m_dot*(h_ex_s-h_su))/W_dot;
    AU = param.AU;
    Q_dot = AU*(T_su - T_amb);
    h_ex = 1/param.eta_hyd*(h_ex_s-h_su*(1-param.eta_hyd));%h_su+(W_dot-Q_dot)/m_dot;
    if h_ex > param.h_min && h_ex < param.h_max
        out.flag = 1;
    else
        out.flag = -1;
    end
    
else
    % If the external conditions are not viable, we fake a perfect machine
    % but we notice the user with a negative flag
    out.flag = -2;
end

if out.flag > 0
    out.h_ex = h_ex;
    if C_oil > 0
        [out.T_ex, ~, ~, ~, ~, ~, ~, ~] = HP_solubMixt(C_oil, P_ex, h_ex, fluid, fluid_lub, Tbubble_min_ex, Tsat_ex, param.fit_DTP_zeta);
    else
        out.T_ex = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex,fluid);
    end
    out.T_su = T_su;
    out.m_dot = m_dot;
    out.W_dot = W_dot;
    out.Q_dot = Q_dot;
    out.epsilon_is = epsilon_is;
    out.epsilon_vol = epsilon_vol;
    out.eta_cavitation = eta_cavitation;
else
    out.T_ex = T_su;
    out.h_ex = h_su;
    out.m_dot = N_pp/60* param.V_s*rho_su;
    out.W_dot = out.m_dot*(h_ex_s-h_su);
    out.Q_dot = 0;
    out.epsilon_is = 1;
    out.epsilon_vol = 1;
    out.T_su = T_su;
    out.eta_cavitation = 1;
end
out.time = toc(tstart_pp);

%% TS DIAGRAM and DISPLAY

% Generate the output variable TS 
%TS.T = [T_su out.T_ex];
%TS.s = [s_su CoolProp.PropsSI('S','H',out.h_ex,'P',P_ex,fluid)];
TS = NaN;
% If the param.displayResults flag is activated (=1), the results are displayed on the
% command window
if param.displayResults ==1
    in.fluid = fluid;
    in.N_pp = N_pp;
    in.T_su = T_su;
    in.h_su = h_su;
    in.P_su = P_su;
    in.P_su = P_ex;
    in.V_s = param.V_s;
    in.modelType= param.modelType;
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
    disp('Results:')
    disp(out)
end
end

