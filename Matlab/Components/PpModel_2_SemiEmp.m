function [out,TS] = PpModel_2_SemiEmp(P_su, h_su, P_ex, fluid, M_dot, param)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 01/06/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% "PpModel_2_SemiEmp.m" is a single matlab function implementing a constant-efficiency model
% of volumetric pump (see the Documentation/PpModel_2_SemiEmp_MatlabDoc). Unlike "PpModel_SemiEmp.m" 
% which imposes N_pp and deduces the mass flow rate, "PpModel_2_SemiEmp.m" imposes the mass flow rate
% and derives the corresponding pump rotational speed.
%
% The model inputs are:
%       - P_su: inlet pressure of the WF                          	[Pa]
%       - h_su: inlet temperature of the WF                        	[J/kg]
%       - P_ex: outlet pressure of the WF                          	[Pa]
%       - fluid: nature of the fluid (string)                       [-]
%       - M_dot: fluid mass flow rate                             	[kg/s]
%       - param: structure variable containing the model parameters
%
% The model paramters provided in 'param' should contain the followings:
%           param.V_s , machine displacement volume               	[m3]
%           param.V, volume of the pump                             [m^3]
%           param.A_leak, leakage surface area                     	[m2]
%           param.W_dot_loss, constant power losses                	[W]
%           param.K_0_loss, term for the proportional losses       	[-]
%           param.h_min, minimum enthalpy of the fluid              [J/kg]
%           param.h_max, maximum enthalpy of the fluid              [J/kg]
%
% The model outputs are:
%       - out: a structure variable which includes
%               - T_ex =  exhaust temperature                    [K]
%               - h_ex =  exhaust enthalpy                       [J/kg]
%               - N_pp = pump rotational speed                   [rpm]
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


%% DEMONSTRATION CASE -- COMMENT THIS SECTION IF EXTERNAL CALL FOR SPEED IMPROVEMENT
% if nargin == 0    
%     % Define a demonstration case if PumpModel.mat is not executed externally
%     fluid = 'R245fa';               %Nature of the fluid
%     P_su = 4.0001e5;                %Supply pressure        [Pa]
%     P_ex = 3.6510e+06*0.99;         %Exhaust pressure       [Pa]
%     h_su = 2.6676e+05;              %Supply enthalpy        [J/kg]
%     M_dot = 0.1;                    %Mass flow rate         [kg/s]
%     param.displayResults = 1;       %Flag to control the resustl display [0/1]   
%     param.V_s = 1e-6;               %Machine swepts volume  [m^3]
%     param.V =1.4e-3;                %Volume inside the pump [m^3]
%     param.A_leak = 0.5e-6;          %Leakage surface area   [m2]
%     param.W_dot_loss = 0;           %Constant power losses  [W]
%     param.K_0_loss = 0.1;           %Term for the proportional losses [-]
% end


%% PUMP MODELING
tstart_pp = tic;                    

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
            
if P_su < P_ex && M_dot > 0
    V_s = param.V_s;
    A_leak = param.A_leak;
    W_dot_loss = param.W_dot_0_loss;
    K_0_loss = param.K_0_loss;
    M_dot_leak = A_leak*sqrt(2*rho_su*(P_ex-P_su));
    M_dot_th = M_dot-M_dot_leak;
    N_pp = 60*M_dot_th/(V_s*rho_su);
    epsilon_vol = M_dot/(N_pp/60*V_s*rho_su);
    W_dot = W_dot_loss + K_0_loss*M_dot/rho_su*(P_ex-P_su);
    epsilon_is = (M_dot*(h_ex_s-h_su))/W_dot;
    h_ex = h_su+W_dot/M_dot;
    if h_ex > param.h_min && h_ex < param.h_max
        out.flag = 1;
    else
        out.flag = -1;
    end
else
    out.flag = -2;
end

if out.flag > 0
    out.h_ex = h_ex;
    out.T_ex = CoolProp.PropsSI('T','P',P_ex,'H',out.h_ex,fluid);
    out.N_pp = N_pp;
    out.W_dot = W_dot;
    out.epsilon_is = epsilon_is;
    out.epsilon_vol = epsilon_vol;
    out.M = (CoolProp.PropsSI('D','H',h_su,'P',P_su,fluid)+CoolProp.PropsSI('D','H',out.h_ex,'P',P_ex,fluid))/2*param.V;
else
    out.T_ex = T_su;
    out.h_ex = h_su;
    out.N_pp = 60*M_dot/(param.V_s*rho_su);
    out.W_dot = M_dot*(h_ex_s-h_su);
    out.epsilon_is = 1;
    out.epsilon_vol = 1;
    out.M =(CoolProp.PropsSI('D','H',h_su,'P',P_su,fluid)+CoolProp.PropsSI('D','H',out.h_ex,'P',P_ex,fluid))/2*param.V;
end

TS.T = [T_su out.T_ex];
TS.s = [s_su CoolProp.PropsSI('S','H',out.h_ex,'P',P_ex,fluid)];
out.time = toc(tstart_pp);

end

