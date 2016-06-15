function [out,TS] = HexModel(fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 27/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% HexModel is a single matlab code implementing six different modelling
% approaches to simulate counter-flow heat exchangers (see the Documentation/HexModel_MatlabDoc)
%
% The model inputs are:
%       - fluid_h: nature of the hot fluid                        	[-]
%       - P_h_su: inlet pressure of the hot fluid                   [Pa]
%       - in_h_su: inlet temperature or enthalpy of the hot fluid   [K or J/kg]
%       - m_dot_h: mass flow rate of the hot fluid                  [kg/s]
%       - fluid_c: nature of the cold fluid                        	[-]
%       - P_c_su: inlet pressure of the cold fluid                  [Pa]
%       - in_c_su: inlet temperature or enthalpy of the cold fluid  [K or J/kg]
%       - m_dot_c: mass flow rate of the cold fluid                 [kg/s]
%       - param: structure variable containing the model parameters
%
% The model paramters provided in 'param' depends of the type of model selected:
%
%       - if param.modelType = 'CstPinch':
%           param.pinch = pinch point value in the temperature profile [K]
%           param.type_h = type of input for hot fluid ('H' for enthalpy,'T' for temperature)
%           param.type_c = type of input for cold fluid ('H' for enthalpy,'T' for temperature)
%           param.V_h_tot = HEX volume for the hot fluid side [m^3]
%           param.V_c_tot = HEX volume for the cold fluid side [m^3]
%           param.displayResults = flag to display the results or not [1/0]
%           param.displayTS = flag to display the temperature profiles or not [1/0]
%
%       - if param.modelType = 'CstEff':
%           param.epsilon_th = effective thermal efficiency of the heat exchanger [-]
%           param.type_h = type of input for hot fluid ('H' for enthalpy,'T' for temperature)
%           param.type_c = type of input for cold fluid ('H' for enthalpy,'T' for temperature)
%           param.V_h_tot = HEX volume for the hot fluid side [m^3]
%           param.V_c_tot = HEX volume for the cold fluid side [m^3]
%           param.displayResults = flag to display the results or not [1/0]
%           param.displayTS = flag to display the temperature profiles or not [1/0]
%
%       - if param.modelType = 'PolEff':
%           param.CoeffPolEff = polynomial coefficients for calculating the effective thermal efficiency of the heat exchanger
%           param.m_dot_h_n = nominal hot fluid mass flow rate [kg/sec]
%           param.m_dot_c_n = nominal cold fluid mass flow rate [kg/sec]
%           param.type_h = type of input for hot fluid ('H' for enthalpy,'T' for temperature)
%           param.type_c = type of input for cold fluid ('H' for enthalpy,'T' for temperature)
%           param.V_h_tot = HEX volume for the hot fluid side [m^3]
%           param.V_c_tot = HEX volume for the cold fluid side [m^3]
%           param.displayResults = flag to display the results or not [1/0]
%           param.displayTS = flag to display the temperature profiles or not [1/0]
%
%       - if param.modelType = 'hConvCst':
%           param.hConv_h_liq = HF convective heat transfer coeff for liquid phase [W/m^2.K]
%           param.hConv_h_tp = HF convective heat transfer coeff for two-phase [W/m^2.K]
%           param.hConv_h_vap = HF convective heat transfer coeff for vapour phase [W/m^2.K]
%           param.hConv_c_liq = CF convective heat transfer coeff for liquid phase [W/m^2.K]
%           param.hConv_c_tp = CF convective heat transfer coeff for two-phase [W/m^2.K]
%           param.hConv_c_vap = CF convective heat transfer coeff for vapour phase [W/m^2.K]
%           param.A_tot (or param.A_h_tot && param.A_c_tot) = total surfac area of HEX [m^2]
%           param.type_h = type of input for hot fluid ('H' for enthalpy,'T' for temperature)
%           param.type_c = type of input for cold fluid ('H' for enthalpy,'T' for temperature)
%           param.V_h_tot = HEX volume for the hot fluid side [m^3]
%           param.V_c_tot = HEX volume for the cold fluid side [m^3]
%           param.displayResults = flag to display the results or not [1/0]
%           param.displayTS = flag to display the temperature profiles or not [1/0]
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
%           param.displayResults = flag to display the results or not [1/0]
%           param.displayTS = flag to display the temperature profiles or not [1/0]
%
%       - if param.modelType = 'hConvCor':
%           param.C_h_liq / param.n_h_liq / param.m_h_liq = correlation parameters for the hot fluid - liquid zone
%           param.C_h_vap / param.n_h_vap / param.m_h_vap = correlation parameters for the hot fluid - vapour zone
%           param.C_h_tp = coefficient for hot fluid - two phase zone
%           param.C_c_liq / param.n_c_liq / param.m_c_liq = correlation parameters for the cold fluid - liquid zone
%           param.C_c_vap / param.n_c_vap / param.m_c_vap = correlation parameters for the cold fluid - liquid zone
%           param.C_c_tp = coefficient for cold fluid - two phase zone
%           param.CS_h / param.Dh_h = cross-section and hydraulic diameter of the hot fluid
%           param.CS_c / param.Dh_c = cross-section and hydraulic diameter of the cold fluid
%           param.A_tot (or param.A_h_tot && param.A_c_tot) = total surfac area of HEX [m^2]
%           param.type_h = type of input for hot fluid ('H' for enthalpy,'T' for temperature)
%           param.type_c = type of input for cold fluid ('H' for enthalpy,'T' for temperature)
%           param.V_h_tot = HEX volume for the hot fluid side [m^3]
%           param.V_c_tot = HEX volume for the cold fluid side [m^3]
%           param.displayResults = flag to display the results or not [1/0]
%           param.displayTS = flag to display the temperature profiles or not [1/0]
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

%% DEMONSTRATION CASE

if nargin == 0    
    % Define a demonstration case if HexModel.mat is not executed externally  
    
    fluid_h = 'PiroblocBasic';      % Nature of the hot fluid           [-]
    m_dot_h = 0.8; %0.5;                  % Mass flow rat of the hot fluid    [kg/s]
    P_h_su =  3e5;                  % Supply pressure of the hot fluid  [Pa]
    in_h_su =  180+273.15; %400;                 % Supply h or T of the hot fluid  	[J/kg pr K]
    param.type_h = 'T';             % Type of inputs for the hot fluid  ['H' or 'T']
    fluid_c = 'R245fa';             % Nature of the cold fluid        	[-]
    m_dot_c = 0.01;%0.015899945748693; %0.05;                 % Mass flow rat of the cold fluid  	[kg/s]
    P_c_su = 7.57e+05; %4e5;                   % Supply pressure of the cold fluid	[Pa]
    in_c_su = 2.453e+05;           % Supply h or T of the cold fluid  	[J/kg pr K]
    param.type_c = 'H';             % Type of inputs for the cold fluid	['H' or 'T']
    param.displayResults = 1;
    param.displayTS = 1;
    param.modelType = 'hConvVar';
    
%     switch param.modelType
%         
%         case 'CstPinch' % Example of paramters for modelType = CstPinch 
%             param.pinch = 5;
%             param.V_h_tot = 1;
%             param.V_c_tot = 1;
%             
%         case 'CstEff'   % Example of paramters for modelType = CstEff 
%             param.epsilon_th = 1;
%             param.V_h_tot = 1;
%             param.V_c_tot = 1;
%             
%         case 'PolEff'   % Example of paramters for modelType = PolEff
%             param.m_dot_c_n = 0.149;
%             param.m_dot_h_n = 0.149;
%             param.CoeffPolEff = [0.913872425354551 1.601103316261962 -1.161513908064705 0.244892183446199 -2.586460052802316 1.912516752474164];  
%             param.V_h_tot = 1;
%             param.V_c_tot = 1;
%             
%         case 'hConvCst' % Example of paramters for modelType = hConvCst            
%             param.hConv_h_liq = 500;
%             param.hConv_h_tp = 10000;
%             param.hConv_h_vap = 200;
%             param.hConv_c_liq = 500;
%             param.hConv_c_tp = 10000;
%             param.hConv_c_vap = 200;
%             param.A_tot = 5.45000;
%             param.V_h_tot = 1;
%             param.V_c_tot = 1;
%             
%         case 'hConvVar' % Example of paramters for modelType = hConvVar         
%             param.m_dot_h_n = 0.618652995944947;
%             param.m_dot_c_n = 0.618652995944947;
%             param.hConv_c_liq_n = 2.379272937774658e+02;
%             param.hConv_c_tp_n = 2.379272937774658e+02;
%             param.hConv_c_vap_n = 2.379272937774658e+02;
%             param.hConv_h_liq_n = 1.125000000000000e+02;
%             param.hConv_h_tp_n = 1.125000000000000e+02;
%             param.hConv_h_vap_n = 1.125000000000000e+02;
%             param.n = 0.8;
%             param.A_tot = 5.45000;
%             param.V_h_tot = 1;
%             param.V_c_tot = 1;
%             
%         case 'hConvCor' % Example of paramters for modelType = hConvCor           
%             param.m_c_liq = 0.7;
%             param.m_c_vap = 0.7;
%             param.n_c_liq = 0.333;
%             param.n_c_vap = 0.333;
%             param.C_c_liq = 0.308120727539063;
%             param.C_c_tp = 0.308120727539063;
%             param.C_c_vap =0.308120727539063;
%             param.CS_c = 0.020845188;
%             param.Dh_c = 0.356328;            
%             param.m_h_liq = 0.7;
%             param.m_h_vap = 0.7;
%             param.n_h_liq = 0.333;
%             param.n_h_vap = 0.333;            
%             param.C_h_liq = 1.833880424499512;
%             param.C_h_tp = 2.999996185302734;
%             param.C_h_vap = 1.833880424499512;
%             param.CS_h = 0.020845188;
%             param.Dh_h = 0.356328;    
%             param.A_tot = 5.45000;
%             param.V_h_tot = 1;
%             param.V_c_tot = 1;
%     end
%         
path = 'C:\Users\RDickes\Google Drive\PhD\MOR study\ORC\Experimental database\Sun2Power';

    EV_folder = [path '\Evaporator\'];
    load([EV_folder, 'ParametersCalibration_EV.mat'])
    param = EV_hConvVar;
    param.V_h_tot = 0.009;
    param.V_c_tot = 0.009;
    param.displayResults = 0;
    param.displayTS = 1;
    param.generateTS = 1;
end

tstart_hex = tic;


%% HEAT EXCHANGER MODELING
% Modelling section of the code

if not(isfield(param,'displayResults'))
    param.displayResults = 0;
    param.displayTS = 0;
    %if nothing specified by the user, the results are not displayed by
    %default.
end

% Evaluation of the hot fluid (HF) supply temperature
if strcmp(param.type_h,'H')
    T_h_su = CoolProp.PropsSI('T','P',P_h_su,'H',in_h_su, fluid_h);
elseif strcmp(param.type_h,'T')
    T_h_su = in_h_su;
end
% Evaluation of the cold fluid (CF) supply temperature
if strcmp(param.type_c,'H')
    T_c_su = CoolProp.PropsSI('T','P',P_c_su,'H',in_c_su, fluid_c);
elseif strcmp(param.type_c,'T')
    T_c_su = in_c_su;
end

if (T_h_su-T_c_su)>1e-2  && m_dot_h  > 0 && m_dot_c > 0;
    % Check if the operating conditions permit a viable heat transfer
    
    switch param.modelType
        %If yes, select the proper model paradigm chosen by the user
        
        case 'CstPinch' % Model which imposes a constant pinch
            if T_h_su-T_c_su > param.pinch
                
                % Power and enthalpy vectors calculation
                lb = 0; % Minimum heat power that can be transferred between the two media
                ub = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param); %Compute maximum heat power that can be transferred between the two media
                f = @(Q_dot) HEX_CstPinch_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param.pinch, Q_dot, param); % function to solve in order to find Q_dot_eff in the heat exchanger
                Q_dot_eff = zeroBrent ( lb, ub, 1e-6, 1e-6, f ); % Solver driving residuals of HEX_CstPinch_res to zero
                out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, param); %evaluate the temperature profile for a given heat power, cfr documentation of HEX_profile
                out.Q_dot_tot = Q_dot_eff;
                out.h_h_ex = out.H_h_vec(1);
                out.h_c_ex = out.H_c_vec(end);
                out.T_h_ex = out.T_h_vec(1);
                out.T_c_ex = out.T_c_vec(end);
                out.resPinch = abs(1-out.pinch/param.pinch);
                
                % Entropy vector calculation
                [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
                if strcmp(param.type_h,'H') %if not an incompressible fluid, calculate entropy vector
                    for i = 1: length(out.H_h_vec)
                        out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                    end
                end
                if strcmp(param.type_c,'H') %if not an incompressible fluid, calculate entropy vector
                    for i = 1: length(out.H_c_vec)
                        out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                    end
                end
                
                % Mass calculation
                out.V_h_vec = param.V_h_tot*(out.Qdot_vec./out.Q_dot_tot);
                out.V_c_vec = param.V_c_tot*(out.Qdot_vec./out.Q_dot_tot);
                for i = 1: length(out.V_h_vec)
                    if strcmp(param.type_h,'H')
                        out.M_h_vec(i) = out.V_h_vec(i)*(CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i+1),fluid_h))/2;
                    else
                        out.M_h_vec(i) = out.V_h_vec(i)*sf_PropsSI_bar('D', out.T_h_vec(i),  out.T_h_vec(i+1), P_h_su, fluid_h);
                    end
                    if strcmp(param.type_c,'H')
                        out.M_c_vec(i) = out.V_c_vec(i)*(CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i+1),fluid_c))/2;
                    else
                        out.M_c_vec(i) = out.V_c_vec(i)*sf_PropsSI_bar('D', out.T_c_vec(i),  out.T_c_vec(i+1), P_c_su, fluid_c);
                    end
                end
                out.M_h = sum(out.M_h_vec);
                out.M_c = sum(out.M_c_vec);
                
                % Flag evaluation
                if out.resPinch <1e-4
                    out.flag = 1;
                else
                    out.flag = -1;
                end
                
            else
                
                % Power and enthalpy vectors calculation
                Q_dot_eff = 0;
                out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, param); %evaluate the temperature profile for a given heat power, cfr documentation of HEX_profile
                out.h_h_ex = out.H_h_vec(1);
                out.h_c_ex = out.H_c_vec(end);
                out.T_h_ex = out.T_h_vec(1);
                out.T_c_ex = out.T_c_vec(end);
                out.Q_dot_tot = Q_dot_eff;
                
                % Entropy vector calculation
                [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
                if strcmp(param.type_h,'H') %if not an incompressible fluid, calculate entropy vector
                    for i = 1: length(out.H_h_vec)
                        out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                    end
                end
                if strcmp(param.type_c,'H') %if not an incompressible fluid, calculate entropy vector
                    for i = 1: length(out.H_c_vec)
                        out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                    end
                end
                
                % Mass calculation
                out.V_h_vec = param.V_h_tot*(out.Qdot_vec./out.Q_dot_tot);
                out.V_c_vec = param.V_c_tot*(out.Qdot_vec./out.Q_dot_tot);
                for i = 1: length(out.V_h_vec)
                    if strcmp(param.type_h,'H')
                        out.M_h_vec(i) = out.V_h_vec(i)*(CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i+1),fluid_h))/2;
                    else
                        out.M_h_vec(i) = out.V_h_vec(i)*sf_PropsSI_bar('D', out.T_h_vec(i),  out.T_h_vec(i+1), P_h_su, fluid_h);
                    end
                    if strcmp(param.type_c,'H')
                        out.M_c_vec(i) = out.V_c_vec(i)*(CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i+1),fluid_c))/2;
                    else
                        out.M_c_vec(i) = out.V_c_vec(i)*sf_PropsSI_bar('D', out.T_c_vec(i),  out.T_c_vec(i+1), P_c_su, fluid_c);
                    end
                end
                out.M_h = sum(out.M_h_vec);
                out.M_c = sum(out.M_c_vec);
                
                % Flag evaluation                
                out.flag = 2;
            end
            
        case 'CstEff'   % Model which imposes a constant thermal efficiency
                        
            % Power and enthalpy vectors calculation
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param); %Compute the maximum heat power that can be transferred between the two media
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, param); %Evaluate temperature profile based on Q_dot_max
            Q_dot_eff = param.epsilon_th*Q_dot_max; %Effective heat transfer
            out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, param); %Evaluate temperature profile based on Q_dot_eff
            out.Q_dot_tot = Q_dot_eff;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);

            % Entropy vector calculation            
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            if strcmp(param.type_h,'H') %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.H_h_vec)
                    out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                end
            end
            if strcmp(param.type_c,'H') %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.H_c_vec)
                    out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                end
            end
            
            % Mass calculation
            out.V_h_vec = param.V_h_tot*(out.Qdot_vec./out.Q_dot_tot);
            out.V_c_vec = param.V_c_tot*(out.Qdot_vec./out.Q_dot_tot);
            for i = 1: length(out.V_h_vec)
                if strcmp(param.type_h,'H')
                    out.M_h_vec(i) = out.V_h_vec(i)*(CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i+1),fluid_h))/2;
                else
                    out.M_h_vec(i) = out.V_h_vec(i)*sf_PropsSI_bar('D', out.T_h_vec(i),  out.T_h_vec(i+1), P_h_su, fluid_h);
                end
                if strcmp(param.type_c,'H')
                    out.M_c_vec(i) = out.V_c_vec(i)*(CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i+1),fluid_c))/2;
                else
                    out.M_c_vec(i) = out.V_c_vec(i)*sf_PropsSI_bar('D', out.T_c_vec(i),  out.T_c_vec(i+1), P_c_su, fluid_c);
                end
            end
            out.M_h = sum(out.M_h_vec);
            out.M_c = sum(out.M_c_vec);
            
            % Flag evaluation             
            if abs(out_max.pinch) < 1e-2 % Check that Q_dot_max correspond to the situation where the pinch is equal to zero
                out.flag = 1;
            else
                out.flag = -2;
            end        
            
        case 'PolEff'   % Model which computes the thermal efficiency  with a polynomial regressions
            
            % Power and enthalpy vectors calculation
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param); %Compute the maximum heat power that can be transferred between the two media
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, param); %Evaluate temperature profile based on Q_dot_max
            Q_dot_eff = Q_dot_max*max(1e-5,min(1, param.CoeffPolEff(1) + param.CoeffPolEff(2)*(m_dot_h/param.m_dot_h_n) + param.CoeffPolEff(3)*(m_dot_c/param.m_dot_c_n) + param.CoeffPolEff(4)*(m_dot_h/param.m_dot_h_n)^2 + param.CoeffPolEff(5)*(m_dot_h/param.m_dot_h_n)*(m_dot_c/param.m_dot_c_n) + param.CoeffPolEff(6)*(m_dot_c/param.m_dot_c_n)^2));
            out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, param); %Evaluate temperature profile based on Q_dot_eff
            out.Q_dot_tot = Q_dot_eff;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);
            
            % Entropy vector calculation            
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            if strcmp(param.type_h,'H') %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.H_h_vec)
                    out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                end
            end
            if strcmp(param.type_c,'H') %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.H_c_vec)
                    out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                end
            end
            
            % Mass calculation
            out.V_h_vec = param.V_h_tot*(out.Qdot_vec./out.Q_dot_tot);
            out.V_c_vec = param.V_c_tot*(out.Qdot_vec./out.Q_dot_tot);
            for i = 1: length(out.V_h_vec)
                if strcmp(param.type_h,'H')
                    out.M_h_vec(i) = out.V_h_vec(i)*(CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i+1),fluid_h))/2;
                else
                    out.M_h_vec(i) = out.V_h_vec(i)*sf_PropsSI_bar('D', out.T_h_vec(i),  out.T_h_vec(i+1), P_h_su, fluid_h);
                end
                if strcmp(param.type_c,'H')
                    out.M_c_vec(i) = out.V_c_vec(i)*(CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i+1),fluid_c))/2;
                else
                    out.M_c_vec(i) = out.V_c_vec(i)*sf_PropsSI_bar('D', out.T_c_vec(i),  out.T_c_vec(i+1), P_c_su, fluid_c);
                end
            end
            out.M_h = sum(out.M_h_vec);
            out.M_c = sum(out.M_c_vec);
            
            % Flag evaluation 
            if abs(out_max.pinch) < 1e-2 % Check that Q_dot_max correspond to the situation where the pinch is equal to zero
                out.flag = 1;
            else
                out.flag = -2;
            end
            
        case 'hConvCst' % 3-zone moving-boundary model with constant convective heat transfer coefficients
            if isfield(param, 'A_tot')
                % if only one surface area is specified, then it is the
                % same for the hot and the cold fluid (ex: CPHEX)
                param.A_h_tot = param.A_tot;
                param.A_c_tot = param.A_tot;
            end
            
            % Power and enthalpy vectors calculation
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param); %Compute the maximum heat power that can be transferred between the two media
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, param); %Evaluate temperature profile based on Q_dot_max
            lb = 0; % Minimum heat power that can be transferred between the two media
            ub = Q_dot_max;
            f = @(Q_dot) HEX_hConvCst_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su,  Q_dot, param); % function to solve in order to find Q_dot_eff in the heat exchanger
            if f(ub) > 0
                Q_dot_eff = ub; % HEX so oversized that the effective heat power is equal to Q_dot_max
            else
                Q_dot_eff = zeroBrent ( lb, ub, 1e-6, 1e-6, f); % Solver driving residuals of HEX_hConvCst_res to zero
            end
            out = HEX_hConvCst(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, param); %Evaluate temperature profile based on Q_dot_eff
            out.Q_dot_tot = Q_dot_eff;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);
            
            % Entropy vector calculation
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            if strcmp(param.type_h,'H') %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.H_h_vec)
                    out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                end
            end
            if strcmp(param.type_c,'H') %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.H_c_vec)
                    out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                end
            end
            
            % Mass calculation
            out.V_h_vec = param.V_h_tot*(out.A_h./param.A_h_tot);
            out.V_c_vec = param.V_c_tot*(out.A_h./param.A_h_tot);
            for i = 1: length(out.V_h_vec)
                if strcmp(param.type_h,'H')
                    out.M_h_vec(i) = out.V_h_vec(i)*(CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i+1),fluid_h))/2;
                else
                    out.M_h_vec(i) = out.V_h_vec(i)*sf_PropsSI_bar('D', out.T_h_vec(i),  out.T_h_vec(i+1), P_h_su, fluid_h);
                end
                if strcmp(param.type_c,'H')
                    out.M_c_vec(i) = out.V_c_vec(i)*(CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i+1),fluid_c))/2;
                else
                    out.M_c_vec(i) = out.V_c_vec(i)*sf_PropsSI_bar('D',out. T_c_vec(i),  out.T_c_vec(i+1), P_c_su, fluid_c);
                end
            end
            out.M_h = sum(out.M_h_vec);
            out.M_c = sum(out.M_c_vec);
            
            % Flag evaluation              
            if out.resA <1e-4
                out.flag = 1;
            else
                if abs(out_max.pinch) < 1e-2
                    if Q_dot_eff == Q_dot_max
                        out.flag = 2;
                    else
                        out.flag = -1;
                    end
                else
                    out.flag = -2;
                end
            end
            
        case 'hConvVar' % 3-zone moving-boundary model with mass-flow dependent convective heat transfer coefficients
            if isfield(param, 'A_tot')
                % if only one surface area is specified, then it is the
                % same for the hot and the cold fluid (ex: CPHEX)
                param.A_h_tot = param.A_tot;
                param.A_c_tot = param.A_tot;
            end
            
            % Power and enthalpy vectors calculation
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param); %Compute the maximum heat power that can be transferred between the two media
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, param); %Evaluate temperature profile based on Q_dot_max
            lb = 0; % Minimum heat power that can be transferred between the two media
            ub = Q_dot_max; % Maximum heat power that can be transferred between the two media
            f = @(Q_dot) HEX_hConvVar_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su,  Q_dot, param); % function to solve in order to find Q_dot_eff in the heat exchanger
            if f(ub) > 0
                Q_dot_eff = ub; % HEX so oversized that the effective heat power is equal to Q_dot_max
            else
                Q_dot_eff = zeroBrent ( lb, ub, 1e-8, 1e-8, f ); % Solver driving residuals of HEX_hConvVar_res to zero
            end
            out = HEX_hConvVar(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, param); %Evaluate temperature profile based on Q_dot_eff
            out.Q_dot_tot = Q_dot_eff;
            out.epsilon_th = Q_dot_eff/Q_dot_max;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);
            
            % Entropy vector calculation
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            if param.generateTS
                if strcmp(param.type_h,'H') %if not an incompressible fluid, calculate entropy vector
                    for i = 1: length(out.H_h_vec)
                        out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                    end
                end
                if strcmp(param.type_c,'H') %if not an incompressible fluid, calculate entropy vector
                    for i = 1: length(out.H_c_vec)
                        out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                    end
                end
            end
            
            % Mass calculation
            out.V_h_vec = param.V_h_tot*(out.A_h./param.A_h_tot);
            out.V_c_vec = param.V_c_tot*(out.A_h./param.A_h_tot);
            for i = 1: length(out.V_h_vec)
                if strcmp(param.type_h,'H')
                    out.M_h_vec(i) = out.V_h_vec(i)*(CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i+1),fluid_h))/2;
                else
                    out.M_h_vec(i) = out.V_h_vec(i)*sf_PropsSI_bar('D', out.T_h_vec(i),  out.T_h_vec(i+1), P_h_su, fluid_h);
                end
                if strcmp(param.type_c,'H')
                    out.M_c_vec(i) = out.V_c_vec(i)*(CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i+1),fluid_c))/2;
                else
                    out.M_c_vec(i) = out.V_c_vec(i)*sf_PropsSI_bar('D', out.T_c_vec(i),  out.T_c_vec(i+1), P_c_su, fluid_c);
                end
            end
            out.M_h = sum(out.M_h_vec);
            out.M_c = sum(out.M_c_vec);
            
            % Flag evaluation 
            if out.resA <1e-4
                out.flag = 1;
            else
                if abs(out_max.pinch) < 1e-2
                    if Q_dot_eff == Q_dot_max
                        out.flag = 2;
                    else
                        %fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param
                        out.flag = -1;
                    end
                    
                else
                    out.flag = -2;
                end
            end
            
        case 'hConvCor' % 3-zone moving-boundary model with Nusselt-dependant convective heat transfer coefficients (user-defined correlation)
            if isfield(param, 'A_tot')
                % if only one surface area is specified, then it is the
                % same for the hot and the cold fluid (ex: CPHEX)
                param.A_h_tot = param.A_tot;
                param.A_c_tot = param.A_tot;
            end
            
            % Power and enthalpy vectors calculation
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param); %Compute the maximum heat power that can be transferred between the two media
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, param); %Evaluate temperature profile based on Q_dot_max
            lb = 0; % Minimum heat power that can be transferred between the two media
            ub = Q_dot_max; % Maximum  heat power that can be transferred between the two media
            f = @(Q_dot) HEX_hConvCor_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su,  Q_dot, param); % function to solve in order to find Q_dot_eff in the heat exchanger
            if f(ub) > 0
                Q_dot_eff = ub; % HEX so oversized that the effective heat power is equal to Q_dot_max
            else
                Q_dot_eff = zeroBrent ( lb, ub, 1e-6, 1e-6, f ); % Solver driving residuals of HEX_hConvCor_res to zero
            end
            out = HEX_hConvCor(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, param); %Evaluate temperature profile based on Q_dot_eff
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            out.Q_dot_tot = Q_dot_eff;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);
            
            % Entropy vector calculation
            if strcmp(param.type_h,'H')  %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.H_h_vec)
                    out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                end
            end
            if strcmp(param.type_c,'H')  %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.H_c_vec)
                    out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                end
            end
            
            % Mass calculation
            out.V_h_vec = param.V_h_tot*(out.A_h./param.A_h_tot);
            out.V_c_vec = param.V_c_tot*(out.A_h./param.A_h_tot);
            for i = 1: length(out.V_h_vec)
                if strcmp(param.type_h,'H')
                    out.M_h_vec(i) = out.V_h_vec(i)*(CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i+1),fluid_h))/2;
                else
                    out.M_h_vec(i) = out.V_h_vec(i)*sf_PropsSI_bar('D', T_h_vec(i),  T_h_vec(i+1), P_h_su, fluid_h);
                end
                if strcmp(param.type_c,'H')
                    out.M_c_vec(i) = out.V_c_vec(i)*(CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i+1),fluid_c))/2;
                else
                    out.M_c_vec(i) = out.V_c_vec(i)*sf_PropsSI_bar('D', T_c_vec(i),  T_c_vec(i+1), P_c_su, fluid_c);
                end
            end
            out.M_h = sum(out.M_h_vec);
            out.M_c = sum(out.M_c_vec);
            
            % Flag evaluation
            if out.resA <1e-4
                out.flag = 1;
            else
                if abs(out_max.pinch) < 1e-2
                    if Q_dot_eff == Q_dot_max
                        out.flag = 2;
                    else
                        out.flag = -1;
                    end
                else
                    out.flag = -2;
                end
            end
            
        otherwise
            disp('Wrong type of model input')
    end
    
else
    %If no, there is not any heat power transfered
    Q_dot_eff = 0;
    out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, param);
    out.h_h_ex = out.H_h_vec(1);
    out.h_c_ex = out.H_c_vec(end);
    out.T_h_ex = out.T_h_vec(1);
    out.T_c_ex = out.T_c_vec(end);
    out.Q_dot_tot = Q_dot_eff;
    
    % Entropy calculation
    [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
    if param.generateTS
        if strcmp(param.type_h,'H') %if not an incompressible fluid, calculate entropy vector
            for i = 1: length(out.H_h_vec)
                out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
            end
        end
        if strcmp(param.type_c,'H') %if not an incompressible fluid, calculate entropy vector
            for i = 1: length(out.H_c_vec)
                out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
            end
        end
    end
    % Mass calculation
    out.V_h_vec = param.V_h_tot*(out.Qdot_vec./out.Q_dot_tot);
    out.V_c_vec = param.V_c_tot*(out.Qdot_vec./out.Q_dot_tot);
    for i = 1: length(out.V_h_vec)
        if strcmp(param.type_h,'H')
            out.M_h_vec(i) = out.V_h_vec(i)*(CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h_su,'H',out.H_h_vec(i+1),fluid_h))/2;
        else
            out.M_h_vec(i) = out.V_h_vec(i)*sf_PropsSI_bar('D', out.T_h_vec(i),  out.T_h_vec(i+1), P_h_su, fluid_h);
        end
        if strcmp(param.type_c,'H')
            out.M_c_vec(i) = out.V_c_vec(i)*(CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c_su,'H',out.H_c_vec(i+1),fluid_c))/2;
        else
            out.M_c_vec(i) = out.V_c_vec(i)*sf_PropsSI_bar('D', T_c_vec(i),  T_c_vec(i+1), P_c_su, fluid_c);
        end
    end
    out.M_h = sum(out.M_h_vec);
    out.M_c = sum(out.M_c_vec);
    
    if (T_h_su-T_c_su)<1e-2  && (T_h_su-T_c_su >0)  && m_dot_h  > 0 && m_dot_c > 0
        out.flag = 3;
    else
        out.flag = -3;
    end
end
out.time = toc(tstart_hex);

%% TS DIAGRAM and DISPLAY

% Generate the output variable TS
TS.T_h = out.T_h_vec;
TS.T_c = out.T_c_vec;
TS.s_h = out.s_h_vec;
TS.s_c = out.s_c_vec;
TS.x = out.x_vec;

% If the param.displayTS flag is activated (=1), the temperature profile is
% plotted in a new figure
if param.displayTS == 1
    figure
    hold on
    plot(out.x_vec, TS.T_c-273.15,'s-' ,'linewidth',2)
    plot(out.x_vec, TS.T_h-273.15,'o-' ,'linewidth',2)
    grid on
    xlabel('Heat power fraction [-]','fontsize',14,'fontweight','bold')
    ylabel('Temperature [°C]','fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')
end

% If the param.displayResults flag is activated (=1), the results are displayed on the
% command window
if param.displayResults ==1
    in.fluid_h = fluid_h;
    in.m_dot_h = m_dot_h;
    in.in_h_su = in_h_su;
    in.type_h = param.type_h;
    in.P_h_su = P_h_su;
    in.fluid_c = fluid_c;
    in.m_dot_c = m_dot_c;
    in.in_c_su = in_c_su;
    in.type_c = param.type_c;
    in.P_c_su = P_c_su;
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
    disp('Results')
    disp(out)
    
end

end


function res = HEX_CstPinch_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, pinch, Q_dot, info)
% function giving the residual committed on the pinch for a given Q_dot
out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info); %evaluate the temperature profile for a given heat power, cfr documentation of HEX_profile
res = pinch - out.pinch;
end

function res = HEX_hConvCst_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
% function giving the residual committed on the HEX surface area for a given Q_dot
out = HEX_hConvCst(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
res = out.resA;
end

function out = HEX_hConvCst(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info); %evaluate the temperature profile for a given heat power, cfr documentation of HEX_profile
[out.hConv_h, out.hConv_c, out.DTlog, out.eff_h, out.eff_c, out.A_h, out.A_c, out.U] = deal(NaN*ones(1,length(out.T_h_vec)-1));
for j = 1:length(out.T_h_vec)-1
    % Hot side heat transfer coefficient
    if strcmp(info.type_h, 'H')
        if isempty(strfind(fluid_h, 'INCOMP:'))
            if (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) < CoolProp.PropsSI('H','P',P_h_su,'Q',0,fluid_h)
                out.hConv_h(j) = info.hConv_h_liq;
            elseif (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) > CoolProp.PropsSI('H','P',P_h_su,'Q',1,fluid_h)
                out.hConv_h(j) = info.hConv_h_vap;
            else
                out.hConv_h(j) = info.hConv_h_tp;
            end
        else
            out.hConv_h(j) = info.hConv_h_liq;
        end
    elseif strcmp(info.type_h, 'T')
        out.hConv_h(j) = info.hConv_h_liq;
    end
    % Cold side heat transfer coefficient
    if strcmp(info.type_c, 'H')
        if isempty(strfind(fluid_c, 'INCOMP:'))
            if (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) < CoolProp.PropsSI('H','P',P_c_su,'Q',0,fluid_c)
                out.hConv_c(j) = info.hConv_c_liq;
            elseif (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) > CoolProp.PropsSI('H','P',P_c_su,'Q',1,fluid_c)
                out.hConv_c(j) = info.hConv_c_vap;
            else
                out.hConv_c(j) = info.hConv_c_tp;
            end
            out.hConv_c(j) = info.hConv_c_liq;
        else
            out.hConv_c(j) = info.hConv_c_liq;
        end
    elseif strcmp(info.type_c, 'T')
        out.hConv_c(j) = info.hConv_c_liq;
    end
    out.DTlog(j) = deltaT_log(out.T_h_vec(j+1), out.T_h_vec(j),out.T_c_vec(j), out.T_c_vec(j+1));
    
    % Hot side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_h'))
        out.eff_h(j) = 1;
    elseif strcmp(info.fin_h, 'none')
        out.eff_h(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_h(j),info.fin_h.k, info.fin_h.th, info.fin_h.r, info.fin_h.B, info.fin_h.H);
        out.eff_h(j) = 1-info.fin_h.omega_f*(1-eta_eff);
    end
    
    % Cold side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_c'))
        out.eff_c(j) = 1;
    elseif strcmp(info.fin_c, 'none')
        out.eff_c(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_c(j), info.fin_c.k, info.fin_c.th, info.fin_c.r, info.fin_c.B, info.fin_c.H);
        out.eff_c(j) = 1-info.fin_c.omega_f*(1-eta_eff);
    end
    
    % Global heat transfer coefficient and zone surface area
    out.U(j) = (1/out.hConv_h(j)/out.eff_h(j) + 1/out.hConv_c(j)/out.eff_c(j)/(info.A_c_tot/info.A_h_tot))^-1;
    out.A_h(j) = out.Qdot_vec(j)/out.DTlog(j)/out.U(j);
    out.A_c(j) = out.A_h(j)*info.A_c_tot/info.A_h_tot;
end
out.A_h_tot = sum(out.A_h);
out.resA = 1 - out.A_h_tot/info.A_h_tot;
end

function res = HEX_hConvVar_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
% function giving the residual committed on the HEX surface area for a given Q_dot
out = HEX_hConvVar(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
res = out.resA;
end

function out = HEX_hConvVar(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info); %evaluate the temperature profile for a given heat power, cfr documentation of HEX_profile
[out.hConv_h, out.hConv_c, out.DTlog, out.eff_h, out.eff_c, out.A_h, out.A_c, out.U] = deal(NaN*ones(1,length(out.H_h_vec)-1));
for j = 1:length(out.T_h_vec)-1
    % Hot side heat transfer coefficient
    if strcmp(info.type_h, 'H')
        if isempty(strfind(fluid_h, 'INCOMP:'))
            if (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) < CoolProp.PropsSI('H','P',P_h_su,'Q',0,fluid_h)
                out.hConv_h(j) = info.hConv_h_liq_n*(m_dot_h/info.m_dot_h_n)^info.n;
            elseif (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) > CoolProp.PropsSI('H','P',P_h_su,'Q',1,fluid_h)
                out.hConv_h(j) = info.hConv_h_vap_n*(m_dot_h/info.m_dot_h_n)^info.n;
            else
                out.hConv_h(j) = info.hConv_h_tp_n*(m_dot_h/info.m_dot_h_n)^info.n;
            end
        else
            out.hConv_h(j) = info.hConv_h_liq_n*(m_dot_h/info.m_dot_h_n)^info.n;
        end
    elseif strcmp(info.type_h, 'T')
        out.hConv_h(j) = info.hConv_h_liq_n*(m_dot_h/info.m_dot_h_n)^info.n;
    end
    
    % Cold side heat transfer coefficient
    if strcmp(info.type_c, 'H')
        if isempty(strfind(fluid_c, 'INCOMP:'))
            if (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) < CoolProp.PropsSI('H','P',P_c_su,'Q',0,fluid_c)
                out.hConv_c(j) = info.hConv_c_liq_n*(m_dot_c/info.m_dot_c_n)^info.n;
            elseif (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) > CoolProp.PropsSI('H','P',P_c_su,'Q',1,fluid_c)
                out.hConv_c(j) = info.hConv_c_vap_n*(m_dot_c/info.m_dot_c_n)^info.n;
            else
                out.hConv_c(j) = info.hConv_c_tp_n*(m_dot_c/info.m_dot_c_n)^info.n;
            end
        else
            out.hConv_c(j) = info.hConv_c_liq_n*(m_dot_c/info.m_dot_c_n)^info.n;
        end
    elseif strcmp(info.type_c, 'T')
        out.hConv_c(j) = info.hConv_c_liq_n*(m_dot_c/info.m_dot_c_n)^info.n;
    end
    out.DTlog(j) = deltaT_log(out.T_h_vec(j+1), out.T_h_vec(j),out.T_c_vec(j), out.T_c_vec(j+1));
    
    % Hot side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_h'))
        out.eff_h(j) = 1;
    elseif strcmp(info.fin_h, 'none')
        out.eff_h(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_h(j), info.fin_h.k, info.fin_h.th, info.fin_h.r, info.fin_h.B, info.fin_h.H);
        out.eff_h(j) = 1-info.fin_h.omega_f*(1-eta_eff);
    end
    
    % Cold side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_c'))
        out.eff_c(j) = 1;
    elseif strcmp(info.fin_c, 'none')
        out.eff_c(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_c(j), info.fin_c.k, info.fin_c.th, info.fin_c.r, info.fin_c.B, info.fin_c.H);
        out.eff_c(j) = 1-info.fin_c.omega_f*(1-eta_eff);
    end
    
    % Global heat transfer coefficient and zone surface area
    out.U(j) = (1/out.hConv_h(j)/out.eff_h(j) + 1/out.hConv_c(j)/out.eff_c(j)/(info.A_c_tot/info.A_h_tot))^-1;
    out.A_h(j) = out.Qdot_vec(j)/out.DTlog(j)/out.U(j);
    out.A_c(j) = out.A_h(j)*info.A_c_tot/info.A_h_tot;
end
out.A_h_tot = sum(out.A_h);
out.resA = 1 - out.A_h_tot/info.A_h_tot;
end

function res = HEX_hConvCor_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
% function giving the residual committed on the HEX surface area for a given Q_dot
out = HEX_hConvCor(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
res = out.resA;
end

function out = HEX_hConvCor(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info); %evaluate the temperature profile for a given heat power, cfr documentation of HEX_profile
[out.A, out.hConv_h, out.hConv_c, out.DTlog] = deal(NaN*ones(1,length(out.H_h_vec)-1));
for j = 1:length(out.T_h_vec)-1
    
    % Hot side heat transfer coefficient (based on the empirical correlation Nu = C Re^m Pr^n)
    if strcmp(info.type_h, 'H')
        if isempty(strfind(fluid_h, 'INCOMP:'))
            if (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) < CoolProp.PropsSI('H', 'P', P_h_su, 'Q', 0, fluid_h)
                mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                k = CoolProp.PropsSI('L', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                out.hConv_h(j) = hConvCor(mu, Pr, k, m_dot_h, info.CS_h, info.Dh_h, info.C_h_liq, info.m_h_liq, info.n_h_liq);
            elseif (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) > CoolProp.PropsSI('H','P',P_h_su,'Q',1,fluid_h)
                mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                k = CoolProp.PropsSI('L', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                out.hConv_h(j) = hConvCor(mu, Pr, k, m_dot_h, info.CS_h, info.Dh_h, info.C_h_vap, info.m_h_vap, info.n_h_vap);
            else
                mu = CoolProp.PropsSI('V', 'Q', 0, 'P', P_h_su, fluid_h);
                Pr = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_h_su, fluid_h);
                k = CoolProp.PropsSI('L', 'Q', 0, 'P', P_h_su, fluid_h);
                hConv_tp = hConvCor(mu, Pr, k, m_dot_h, info.CS_h, info.Dh_h, info.C_h_liq, info.m_h_liq, info.n_h_liq);
                out.hConv_h(j) = info.C_h_tp*hConv_tp;
            end
        else
            mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
            Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
            k = CoolProp.PropsSI('L', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
            out.hConv_h(j) = hConvCor(mu, Pr, k, m_dot_h, info.CS_h, info.Dh_h, info.C_h_liq, info.m_h_liq, info.n_h_liq);
        end
    elseif strcmp(info.type_h, 'T')
        cp = sf_PropsSI_bar('C', out.T_h_vec(j), out.T_h_vec(j+1), P_h_su, fluid_h);
        k = sf_PropsSI_bar('L', out.T_h_vec(j), out.T_h_vec(j+1), P_h_su, fluid_h);
        mu = sf_PropsSI_bar('V', out.T_h_vec(j), out.T_h_vec(j+1), P_h_su, fluid_h);
        Pr = cp*mu/k;
        out.hConv_h(j) = hConvCor(mu, Pr, k, m_dot_h, info.CS_h, info.Dh_h, info.C_h_liq, info.m_h_liq, info.n_h_liq);
    end
    
    % Cold side heat transfer coefficient (based on the empirical correlation Nu = C Re^m Pr^n)
    if strcmp(info.type_c, 'H')
        if isempty(strfind(fluid_c, 'INCOMP:'))
            if (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) < CoolProp.PropsSI('H', 'P', P_c_su, 'Q', 0, fluid_c)
                mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                k = CoolProp.PropsSI('L', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                out.hConv_c(j) = hConvCor(mu, Pr, k, m_dot_c, info.CS_c, info.Dh_c, info.C_c_liq, info.m_c_liq, info.n_c_liq);
            elseif (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) > CoolProp.PropsSI('H','P',P_c_su,'Q',1,fluid_c)
                mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                k = CoolProp.PropsSI('L', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                out.hConv_c(j) = hConvCor(mu, Pr, k, m_dot_c, info.CS_c, info.Dh_c, info.C_c_vap, info.m_c_vap, info.n_c_vap);
            else
                mu = CoolProp.PropsSI('V', 'Q', 0, 'P', P_c_su, fluid_c);
                Pr = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_c_su, fluid_c);
                k = CoolProp.PropsSI('L', 'Q', 0, 'P', P_c_su, fluid_c);
                hConv_tp = hConvCor(mu, Pr, k, m_dot_c, info.CS_c, info.Dh_c, info.C_c_liq, info.m_c_liq, info.n_c_liq);
                out.hConv_c(j) = info.C_c_tp*hConv_tp;
            end
        else
            mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
            Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
            k = CoolProp.PropsSI('L', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
            out.hConv_c(j) = hConvCor(mu, Pr, k, m_dot_c, info.CS_c, info.Dh_c, info.C_c_liq, info.m_c_liq, info.n_c_liq);
        end
    elseif strcmp(info.type_c, 'T')
        cp = sf_PropsSI_bar('C', T_c_vec(j), T_c_vec(j+1), P_c_su, fluid_c);
        k = sf_PropsSI_bar('L', T_c_vec(j), T_c_vec(j+1), P_c_su, fluid_c);
        mu = sf_PropsSI_bar('V', T_c_vec(j), T_c_vec(j+1), P_c_su, fluid_c);
        Pr = cp*mu/k;
        out.hConv_c(j) = hConvCor(mu, Pr, k, m_dot_c, info.CS_c, info.Dh_c, info.C_c_liq, info.m_c_liq, info.n_c_liq);
    end
    out.DTlog(j) = deltaT_log(out.T_h_vec(j+1), out.T_h_vec(j),out.T_c_vec(j), out.T_c_vec(j+1));
    
    % Hot side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_h'))
        out.eff_h(j) = 1;
    elseif strcmp(info.fin_h, 'none')
        out.eff_h(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_h(j), info.fin_h.k, info.fin_h.th, info.fin_h.r, info.fin_h.B, info.fin_h.H);
        out.eff_h(j) = 1-info.fin_h.omega_f*(1-eta_eff);
    end
    
    % Cold side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_c'))
        out.eff_c(j) = 1;
    elseif strcmp(info.fin_c, 'none')
        out.eff_c(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_c(j), info.fin_c.k, info.fin_c.th, info.fin_c.r, info.fin_c.B, info.fin_c.H);
        out.eff_c(j) = 1-info.fin_c.omega_f*(1-eta_eff);
    end
    
    % Global heat transfer coefficient and zone surface area
    out.U(j) = (1/out.hConv_h(j)/out.eff_h(j) + 1/out.hConv_c(j)/out.eff_c(j)/(info.A_c_tot/info.A_h_tot))^-1;
    out.A_h(j) = out.Qdot_vec(j)/out.DTlog(j)/out.U(j);
    out.A_c(j) = out.A_h(j)*info.A_c_tot/info.A_h_tot;
end
out.A_h_tot = sum(out.A_h);
out.resA = 1 - out.A_h_tot/info.A_h_tot;
end

function [DT_log, DTh, DTc ]= deltaT_log(Th_su, Th_ex, Tc_su, Tc_ex)
% function that provides the mean logarithm temperature difference between two fluids
DTh = max(Th_su-Tc_ex,1e-2);
DTc = max(Th_ex-Tc_su,1e-2);
if DTh ~= DTc;
    DT_log = (DTh-DTc)/log(DTh/DTc);
else
    DT_log = DTh;
end
end

function eta_fin = FinSchmidt(hConv, k, th, r, B, H)
% functions that compute the fin efficiency based on Schmidt's theory and geometrical data of the HEX
m = sqrt(2*hConv/k/th);
phi_f = B/r;
beta_f = H/B;
R_e = r*1.27*phi_f*(beta_f-0.3)^0.5;
phi = (R_e/r - 1)*(1+0.35*log(R_e/r));
eta_fin = tanh(m*R_e*phi)/(m*R_e*phi);
end

function h = hConvCor(mu, Pr, k, m_dot, CS, Dh, C, m, n)
% function that computes the convective heat transfer coefficient of a single-phase flow base on a user-defined Nusselt correlation
G = m_dot/CS;
Re = G*Dh/mu;
Nu = C*Re^m*Pr^n;
h = Nu*k/Dh;
end