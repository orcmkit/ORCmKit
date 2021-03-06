function [out,TS] = HexModel_2_hConvVar(fluid_h, P_h, in_h_su, m_dot_h, fluid_c, P_c, in_c_su, Q_dot, param)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 01/06/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% "HexModel_2_hConvVar.m" is a matlab code implementing a multi-zone heat exchanger
% model with flow-dependent heat transfer coefficient (see the Documentation/HexModel_hConvVar_MatlabDoc)
% Unlike "HexModel_hConvVar.m" which computes the effective heat transfer
% based on the supply conditions, "HexModel_2_hConvVar.m" derives the cold
% stream flow rate required to meet an imposed heat power Q_dot
%
% The model inputs are:
%       - fluid_h: nature of the hot fluid                        	[-]
%       - P_h: pressure of the hot fluid                            [Pa]
%       - in_h_su: inlet temperature or enthalpy of the hot fluid   [K or J/kg]
%       - m_dot_h: mass flow rate of the hot fluid                  [kg/s]
%       - fluid_c: nature of the cold fluid                        	[-]
%       - P_c: pressure of the cold fluid                           [Pa]
%       - in_c_su: inlet temperature or enthalpy of the cold fluid  [K or J/kg]
%       - Q_dot: heat power to transfer                             [kg/s]
%       - param: structure variable containing the model parameters
%
% The model paramters provided in 'param' should include the following variables:
%
%       - if param.modelType = 'hConvVar':
%           param.m_dot_h_n = HF nominal mass flow rate [kg/s]
%           param.hConv_h_liq_n = HF nominal convective heat transfer coeff for liquid phase [W/m^2.K]
%           param.hConv_h_tp_n = HF nominal convective heat transfer coeff for two-phase [W/m^2.K]
%           param.hConv_h_vap_n = HF nominal convective heat transfer coeff for vapour phase [W/m^2.K]
%           param.m_dot_c_n = CF nominal mass flow rate [kg/s]
%           param.hConv_c_liq_n = CF nominal convective heat transfer coeff for liquid phase [W/m^2.K]
%           param.hConv_c_tp_n = CF nominal convective heat transfer coeff for two-phase [W/m^2.K]
%           param.hConv_c_vap_n = CF nominal convective heat transfer coeff for vapour phase [W/m^2.K]
%           param.A_tot (or param.A_h_tot && param.A_c_tot) = total surfac area of HEX [m^2]
%           param.type_h = type of input for hot fluid ('H' for enthalpy,'T' for temperature)
%           param.type_c = type of input for cold fluid ('H' for enthalpy,'T' for temperature)
%           param.V_h_tot = HEX volume for the hot fluid side [m^3]
%           param.V_c_tot = HEX volume for the cold fluid side [m^3]
%           param.displayTS = flag to display the temperature profiles or not [1/0]
%           param.generateTS = flag to generate the TS output [1/0]
%
% The model outputs are:
%       - out: a structure variable which includes at miniumum the following information:
%               - x_vec =  vector of power fraction in each zone
%               - Qdot_vec =  vector of heat power in each zone [W]
%               - H_h_vec = HF enthalpy vector                  [J/kg]
%               - H_c_vec = CF enthalpy vector                  [J/kg]
%               - T_h_vec = HF temperature vector               [K]
%               - T_c_vec = CF temperature vector               [K]
%               - s_h_vec = HF entropy vector                   [J/kg.K]
%               - s_c_vec = HF entropy vector                   [J/kg.K]
%               - DT_vec = Temperature difference vector        [K]
%               - pinch =  pinch point value                  	[K]
%               - h_h_ex =  HF exhaust enthalpy                 [J/kg]
%               - T_h_ex =  HF exhaust temperature              [K]
%               - h_c_ex =  CF exhaust enthalpy                 [J/kg]
%               - T_c_ex =  CF exhaust temperature              [K]
%               - V_h_vec = HF volume vector                    [m^3]
%               - V_c_vec = CF volume vector                    [m^3]
%               - M_h_vec = HF mass vector                      [kg]
%               - M_c_vec = CF mass vector                      [kg]
%               - M_h = Mass of hot fluid in the HEX            [kg]
%               - M_c = Mass of hot fluid in the HEX            [kg]
%               - time = the code computational time          	[sec]
%               - flag = simulation flag
%
%       - TS : a stucture variable which contains the vectors of temperature
%              and entropy of the fluid (useful to generate a Ts diagram
%              when modelling the entire ORC system
%
% See the documentation for further details or contact rdickes@ulg.ac.be

%% DEMONSTRATION CASE -- COMMENT THIS SECTION IF EXTERNAL CALL FOR SPEED IMPROVEMENT
% if nargin == 0
% % Define a demonstration case if HexModel.mat is not executed externally
% fluid_h = 'R245fa';      % Nature of the hot fluid           [-]
% m_dot_h = 0.1;                  % Mass flow rat of the hot fluid    [kg/s]
% P_h =  3e5;                     % Supply pressure of the hot fluid  [Pa]
% in_h_su =  CoolProp.PropsSI('H', 'T', 35+273.15, 'P', P_h, fluid_h);                 % Supply h or T of the hot fluid  	[J/kg pr K]
% h_ex = CoolProp.PropsSI('H', 'T', 25+273.15, 'P', P_h, fluid_h);
% param.type_h = 'H';             % Type of inputs for the hot fluid  ['H' or 'T']
% param.port_h = 'su';
% fluid_c = 'air';             % Nature of the cold fluid        	[-]
% Q_dot = m_dot_h*(in_h_su-h_ex);                    % Mass flow rat of the cold fluid  	[kg/s]
% P_c = 1e5;                      % Supply pressure of the cold fluid	[Pa]
% in_c_su = CoolProp.PropsSI('H', 'T', 5+273.15, 'P', P_c, fluid_c);           % Supply h or T of the cold fluid  	[J/kg pr K]
% param.type_c = 'H';           % Type of inputs for the cold fluid	['H' or 'T']
% param.port_c = 'su';
% param.displayResults = 1;
% param.displayTS = 1;
% param.m_dot_h_n = 0.618652995944947;
% param.m_dot_c_n = 0.618652995944947;
% param.hConv_c_liq_n = 2.379272937774658e+02;
% param.hConv_c_tp_n = 2.379272937774658e+02;
% param.hConv_c_vap_n = 2.379272937774658e+02;
% param.hConv_h_liq_n = 1.125000000000000e+02;
% param.hConv_h_tp_n = 1.125000000000000e+02;
% param.hConv_h_vap_n = 1.125000000000000e+02;
% param.n = 0.8;
% param.A_h_tot = 5.45000;
% param.A_c_tot = 5.45000;
% param.V_h_tot = 1;
% param.V_c_tot = 1;
% param.generateTS = 1;
% 
% end


%% HEAT EXCHANGER MODELING
tstart_hex = tic;
% Evaluation of the hot fluid (HF) supply temperature
% if strcmp(param.type_h,'H')
%     T_h_su = CoolProp.PropsSI('T','P',P_h,'H',in_h_su, fluid_h);
% elseif strcmp(param.type_h,'T')
%     T_h_su = in_h_su;
% end
% % Evaluation of the cold fluid (CF) supply temperature
% if strcmp(param.type_c,'H')
%     T_c_su = CoolProp.PropsSI('T','P',P_c,'H',in_c_su, fluid_c);
% elseif strcmp(param.type_c,'T')
%     T_c_su = in_c_su;
% end


% Power and enthalpy vectors calculation
m_dot_c_min = HEX_Mdot_c_min(fluid_h, m_dot_h, P_h, in_h_su, Q_dot, fluid_c, P_c, in_c_su, param); %Compute the minimum cold stream flow rate
lb = m_dot_c_min; % Minimum heat power that can be transferred between the two media
ub = 1e2; % Maximum heat power that can be transferred between the two media
f = @(m_dot_c) FCT_Hex_hConvVar_res(fluid_h, m_dot_h, P_h, in_h_su, fluid_c, m_dot_c, P_c, in_c_su,  Q_dot, param); % function to solve in order to find Q_dot_eff in the heat exchanger
if f(lb) > 0
    m_dot_c = m_dot_c_min; % HEX so oversized that the effective heat power is equal to Q_dot_max (pinch equal to 0)
else
    m_dot_c = zeroBrent ( lb, ub, 1e-6, 1e-6, f ); % Solver driving residuals of HEX_hConvVar_res to zero

end


out = FCT_Hex_hConvVar(fluid_h, m_dot_h, P_h, in_h_su, fluid_c, m_dot_c, P_c, in_c_su, Q_dot, param); %Evaluate temperature profile based on Q_dot_eff
out.m_dot_c = m_dot_c;
out.Q_dot_tot = Q_dot;
out.h_h_ex = out.H_h_vec(1);
out.h_c_ex = out.H_c_vec(end);
out.T_h_ex = out.T_h_vec(1);
out.T_c_ex = out.T_c_vec(end);

% Entropy vector calculation
[out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
if param.generateTS
    if strcmp(param.type_h,'H') %if not an incompressible fluid, calculate entropy vector
        for i = 1: length(out.H_h_vec)
            out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h,'H',out.H_h_vec(i),fluid_h);
        end
    end
    if strcmp(param.type_c,'H') %if not an incompressible fluid, calculate entropy vector
        for i = 1: length(out.H_c_vec)
            out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c,'H',out.H_c_vec(i),fluid_c);
        end
    end
end

% Mass calculation
out.V_h_vec = param.V_h_tot*(out.A_h./param.A_h_tot);
out.V_c_vec = param.V_c_tot*(out.A_h./param.A_h_tot);
[out.M_h_vec, out.M_c_vec] = deal(NaN*ones(1, length(out.V_h_vec)));
for i = 1: length(out.V_h_vec)
    if strcmp(param.type_h,'H')
        out.M_h_vec(i) = out.V_h_vec(i)*(CoolProp.PropsSI('D','P',P_h,'H',out.H_h_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h,'H',out.H_h_vec(i+1),fluid_h))/2;
    else
        out.M_h_vec(i) = out.V_h_vec(i)*sf_PropsSI_bar('D', out.T_h_vec(i),  out.T_h_vec(i+1), P_h, fluid_h);
    end
    if strcmp(param.type_c,'H')
        out.M_c_vec(i) = out.V_c_vec(i)*(CoolProp.PropsSI('D','P',P_c,'H',out.H_c_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c,'H',out.H_c_vec(i+1),fluid_c))/2;
    else
        out.M_c_vec(i) = out.V_c_vec(i)*sf_PropsSI_bar('D', out.T_c_vec(i),  out.T_c_vec(i+1), P_c, fluid_c);
    end
end
out.M_h = sum(out.M_h_vec);
out.M_c = sum(out.M_c_vec);

% Flag evaluation
if abs(out.resA) <1e-4 || f(lb) > 0
    out.flag = 1;
else
    out.flag = -2;
end

out.time = toc(tstart_hex);

TS.T_h = out.T_h_vec;
TS.T_c = out.T_c_vec;
TS.s_h = out.s_h_vec;
TS.s_c = out.s_c_vec;
TS.x = out.x_vec;

if param.displayTS
    figure
    hold on
    plot(out.x_vec, TS.T_c-273.15,'s-' ,'linewidth',2)
    plot(out.x_vec, TS.T_h-273.15,'o-' ,'linewidth',2)
    grid on
    xlabel('Heat power fraction [-]','fontsize',14,'fontweight','bold')
    ylabel('Temperature [�C]','fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')
end

end