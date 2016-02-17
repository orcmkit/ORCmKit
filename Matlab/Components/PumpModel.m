function [out,TS] = PumpModel(P_su, h_su, P_ex, fluid, N_pp, param)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Rémi Dickes - 26/01/2015 (University of Liège, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% PumpModel is a single matlab code implementing three different modeling
% approaches to simulate a volumetric pump (see the Documentation)
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
%           param.V_s , machine displacement volume                	[m³]
%           param.epsilon_is, isentropic efficiency                	[-]
%           param.epsilon_vol, volumetric efficiency               	[-]
%           param.displayResults, flag to display the results or not [1/0]
%
%       - if param.modelType = 'PolEff':
%           param.V_s , machine displacement volume                	[m³]
%           param.N_pp_nom, pump nominal shaft speed              	[rpm]
%           param.coeffPol_is, polynmial coef for epsilon_is        [-]
%           param.coeffPol_vol, polynmial coef for epsilon_vol      [-]
%           param.displayResults, flag to display the results or not [1/0]
%
%       - if param.modelType = 'SemiEmp':
%           param.V_s , machine displacement volume               	[m³]
%           param.A_leak, leakage surface area                     	[m²]
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
%               - time = the code computational time             [sec]
%               - flag = simulation flag                         [-1/1]
%
%
% See the documentation for further details or contact rdickes@ulg.ac.be


%% DEMONSTRATION CASE
% Define a demonstration case if PumpModel is executed without
% user-defined inputs

if nargin == 0
    fluid = 'R245fa';               %Nature of the fluid
    P_su = 4.0001e5;                %Supply pressure
    P_ex = 3.6510e+06*0.99;         %Exhaust pressure
    h_su = 2.6676e+05;              %Supply enthalpy
    param.modelType = 'CstEff';     %Type of model
    param.V_s = 1e-6;               %Machine swepts volume 
    N_pp = 1500;                    %Rotational speed 
    param.epsilon_is = 0.5;         %Cst isentropic efficiency
    param.epsilon_vol = 0.8;        %Cst volumetric efficiencyh
    param.displayResults = 1;       %Flag to control the resustl display
end

tstart_pp = tic;                    %Start to evaluate the simulation time

%% PUMP MODELING
% Modeling section of the code

if P_su < P_ex && N_pp > 0      
   %If the external conditions are viable, we proceed to the modeling
    
    switch param.modelType
        case 'CstEff'       
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
            rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
            h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
            V_s = param.V_s;
            epsilon_is = param.epsilon_is;
            epsilon_vol = param.epsilon_vol;
            m_dot = N_pp/60*V_s*epsilon_vol*rho_su;
            W_dot = m_dot*(h_ex_s-h_su)/epsilon_is;
            h_ex = h_su+W_dot/m_dot;
            T_ex = CoolProp.PropsSI('T','P',P_ex,'H',h_ex,fluid);
        case 'PolEff'
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
            rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
            h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
            V_s = param.V_s;
            N_pp_nom = param.N_pp_nom;
            a_is = param.coeffPol_is;
            a_vol = param.coeffPol_vol;
            epsilon_is = max(0.01,min(a_is(1) + a_is(2)*(P_ex/P_su) + a_is(3)*(N_pp/N_pp_nom) + a_is(4)*(P_ex/P_su)^2 + a_is(5)*(P_ex/P_su)*(N_pp/N_pp_nom) + a_is(6)*(N_pp/N_pp_nom)^2,1));
            epsilon_vol = max(0.01,min(a_vol(1) + a_vol(2)*(P_ex/P_su) + a_vol(3)*(N_pp/N_pp_nom) + a_vol(4)*(P_ex/P_su)^2 + a_vol(5)*(P_ex/P_su)*(N_pp/N_pp_nom) + a_vol(6)*(N_pp/N_pp_nom)^2,1));
            m_dot = N_pp/60*V_s*epsilon_vol*rho_su;
            W_dot = max(0,m_dot*(h_ex_s-h_su)/epsilon_is);
            h_ex = h_su + W_dot/m_dot;
            T_ex = CoolProp.PropsSI('T','P',P_ex,'H',h_ex,fluid);
        case 'SemiEmp'
            T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
            s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
            rho_su = CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
            h_ex_s = CoolProp.PropsSI('H','P',P_ex,'S',s_su,fluid);
            V_s = param.V_s;
            A_leak = param.A_leak;
            W_dot_loss = param.W_dot_0_loss;
            K_0_loss = param.K_0_loss;
            m_dot = max(1e-20, (N_pp/60*V_s*rho_su)-(A_leak*sqrt(2*rho_su*(P_ex-P_su))));
            epsilon_vol = m_dot/(N_pp/60*V_s*rho_su);
            W_dot = W_dot_loss + K_0_loss*m_dot/rho_su*(P_ex-P_su);
            epsilon_is = (m_dot*(h_ex_s-h_su))/W_dot;
            h_ex = h_su+W_dot/m_dot;
            T_ex = CoolProp.PropsSI('T','P',P_ex,'H',h_ex,fluid);
        otherwise
            disp('Error: type of pump model not valid');
    end  
    out.T_ex = T_ex;
    out.h_ex = h_ex;
    out.m_dot = m_dot;
    out.W_dot = W_dot;
    out.epsilon_is = epsilon_is;
    out.epsilon_vol = epsilon_vol;
    out.time = toc(tstart_pp);
    out.flag = 1;
else 
    % If the external conditions are not viable, we fake a perfect machine 
    % but we notice the user with a negative flag
    T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
    s_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
    rho_su = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
    out.T_ex = T_su;
    out.h_ex = h_su;
    out.m_dot =  N_pp/60*param.V_s*rho_su;
    out.W_dot = 0;
    out.epsilon_is = 0;
    out.epsilon_vol = 1;
    out.time = toc(tstart_pp);
    out.flag = -1;
end


%% TS DIAGRAM and DISPLAY

% Generate the output variable TS 
TS.T = [T_su out.T_ex];
TS.s = [s_su CoolProp.PropsSI('S','H',out.h_ex,'P',P_ex,fluid)];

% If the param.displayResults flag is activated, the results are shown on the
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

