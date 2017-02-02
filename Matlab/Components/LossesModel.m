function [out, TS] = LossesModel(fluid, in_P, in_h, M_dot, T_amb, param)
% fluid, in_P, in_h, M_dot, T_amb, param
%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems
%
% Remi Dickes - 27/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes@ulg.ac.be
%
% LossesModel.mat is a single matlab function implementing three different modelling
% approaches to simulate pressure drop and heat losses in the pipelines (see 
% the Documentation/LossesModel_MatlabDoc)
%
% The model inputs are:
%       - in_P: iput pressure of the WF                          	[Pa]
%       - in_h: iput temperature of the WF                        	[J/kg]
%       - M_dot: mass flow rate of the WF                          	[Pa]
%       - fluid: nature of the fluid (string)                       [-]
%       - T_amb : ambien temperature                                [K]
%       - param: structure variable containing the model parameters
%
% The model paramters provided in 'param' depends of the type of model selected:
%       - if param.modelType = 'CstDP':
%           param.dp , cst pressure drop                            [Pa]            
%           param.AU , heat losses coefficient                      [W/K]              
%           param.type_in, 'ex' (if input h/P are the related to the exhaust) or 'su' (if input h/P are related to the supply)
%           param.displayResults, flag to display the results or not[1/0]
%
%       - if param.modelType = 'MdotDP':
%           param.funMdot_dp(M_dot), a user-definedcorrelation                 
%           param.AU , heat losses coefficient                      [W/K]              
%           param.type_in, 'ex' (if input h/P are the related to the exhaust) or 'su' (if input h/P are related to the supply)
%           param.displayResults, flag to display the results or not[1/0]
%
%       - if param.modelType = 'PhiDP':
%           param.funPhi_dp(M_dot), a user-definedcorrelation                 
%           param.AU , heat losses coefficient                      [W/K]              
%           param.type_in, 'ex' (if input h/P are the related to the exhaust) or 'su' (if input h/P are related to the supply)
%           param.type_phi, 'phi_su', 'phi_1' or 'phi_0'
%           param.displayResults, flag to display the results or not[1/0]
%
% The model outputs are:
%       - out: a structure variable which includes
%               - P_ex =  exhaust pressure                      [Pa]
%               - P_su =  supply pressure                       [Pa]
%               - T_ex =  exhaust temperature               	[K]
%               - T_su =  supply temperature                   	[K]
%               - h_ex =  exhaust enthalpy                      [J/kg]
%               - h_su =  supply enthalpy                       [J/kg]
%               - dp = pressure drop                            [Pa]
%               - Q_dot = ambiant losses                        [W]
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
    % Define a demonstration case if LossesModel.mat is not executed externally  
    fluid = 'R245fa';
    in_P = 7.2841e+04; %1.2e+05;
    in_h = 2.0871e+05; %CoolProp.PropsSI('H', 'Q', 0, 'P', in_P, fluid); %2.2121e+05;    
    M_dot = 0.0760; %0.0281;
    T_amb = 278.1500;
    path = 'C:\Users\RDickes\Google Drive\PhD\MOR study\ORC\Experimental database\Sun2Power';
    DP_folder = [path '\PressureDrops\'];
    load([DP_folder, 'ParametersCalibration_DP.mat'])
    param = DPLP_PhiDP;
    param.type_in = 'ex';
    param.displayResults = 1;
end
tstart_dp = tic;


%% LOSSES  MODELING
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
param.P_max = CoolProp.PropsSI('Pcrit','H',in_h,'P',in_P,fluid)-1e3;
switch param.modelType
    case 'CstDP'
        if strcmp(param.type_in, 'su')
            out = CstDP(in_h, in_P, M_dot, T_amb, fluid, param);
            if out.h_ex < param.h_min || out.h_ex > param.h_max
                out.flag = -1;
            else
                out.flag = 1;
            end
        elseif strcmp(param.type_in, 'ex')
                x0 = [in_h in_P];
                ub = 1.5*x0;
                f = @(x) CstDP_res(x, ub, in_h, in_P, M_dot, T_amb, fluid, param);
                options = optimoptions('fsolve','Display','none');
                x = fsolve(f, x0./ub,options);
                out = CstDP(x(1)*ub(1), x(2)*ub(2), M_dot, T_amb, fluid, param);
                if norm(f(x)) < 1e-3 && out.h_ex > param.h_min && out.h_ex < param.h_max
                    out.flag = 1;
                else
                    out.flag = -1;
                end
        end
        
    case 'MdotDP'
        if strcmp(param.type_in, 'su')
            out = MdotDP(in_h, in_P, M_dot, T_amb, fluid, param);
            if out.h_ex < param.h_min || out.h_ex > param.h_max
                out.flag = -1;
            else
                out.flag = 1;
            end
        elseif strcmp(param.type_in, 'ex')
            
            x0 = [in_h in_P];
            ub = 1.5*x0;
            f = @(x) MdotDP_res(x, ub, in_h, in_P, M_dot, T_amb, fluid, param);
            options = optimoptions('fsolve','Display','none');
            x = fsolve(f, x0./ub,options);            
            out = MdotDP(x(1)*ub(1), x(2)*ub(2), M_dot, T_amb, fluid, param);
            out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
            out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
            if norm(f(x)) < 1e-3 && out.h_ex > param.h_min && out.h_ex < param.h_max
                out.flag = 1;
            else
                out.flag = -1;
            end
        end
           
    case 'PhiDP'
        if strcmp(param.type_in, 'su')
            out = PhiDP(in_h, in_P, M_dot, T_amb, fluid, param);
            if out.h_ex < param.h_min || out.h_ex > param.h_max
                out.flag = -1;
            else
                out.flag = 1;
            end
        elseif strcmp(param.type_in, 'ex')
            k_0 = [1.5 1.1 1 2 3 4 5 7 10];
            stop = 0;
            k = 1;
            try
                while not(stop) && k <= length(k_0)
                    x0 = [in_h k_0(k)*in_P];
                    ub = x0;
                    f = @(x) PhiDP_res(x, ub, in_h, in_P, M_dot, T_amb, fluid, param);
                    options = optimoptions('fsolve','Display','none');
                    x = fsolve(f, x0./ub,options);
                    out = PhiDP(x(1)*ub(1), x(2)*ub(2), M_dot, T_amb, fluid, param);
                    if norm(f(x)) < 1e-3 && out.h_ex > param.h_min && out.h_ex < param.h_max
                        out.flag = 1;
                        stop = 1;
                    else
                        out.flag = -1;
                        stop = 0;
                    end
                    k = k+1;
                end
            catch
                out.flag = -1;
            end
        end
          
    otherwise
        disp('Wrong type of model input')
end

if out.flag > 0
    out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
    out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
else
    out.h_ex = in_h;
    out.P_ex = in_P;
    out.h_su = in_h;
    out.P_su = in_P;
    out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
    out.T_ex = out.T_su;
    out.dp = 0;
    out.Q_dot = 0;
end
out.time = toc(tstart_dp);
out = orderfields(out);

%% TS DIAGRAM and DISPLAY
TS.T(1) = out.T_su;
TS.T(2) = out.T_ex;
TS.s(1) = CoolProp.PropsSI('S','P',out.P_su,'H',out.h_su,fluid);
TS.s(2) = CoolProp.PropsSI('S','P',out.P_ex,'H',out.h_ex,fluid);

if param.displayResults ==1
    in.fluid = fluid;
    in.M_dot = M_dot;
    in.P = in_P;
    in.T = CoolProp.PropsSI('T','P',in_P,'H',in_h,fluid);
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

function out = CstDP(h_su, P_su, M_dot, T_amb, fluid, param)
out.dp = param.dp;
out.P_su = P_su;
out.h_su = h_su;
out.P_ex = out.P_su-out.dp;
out.Q_dot = param.AU*(CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_su,fluid)-T_amb);
out.h_ex = h_su-out.Q_dot/M_dot;
end

function res = CstDP_res(x, ub, h_ex, P_ex, M_dot, T_amb, fluid, param)
x = x.*ub;
out = CstDP(x(1), x(2), M_dot, T_amb, fluid, param);
res(1) = h_ex-out.h_ex;
res(2) = P_ex-out.P_ex;
end

function out = MdotDP(h_su, P_su, M_dot, T_amb, fluid, param)
out.P_su = P_su;
out.h_su = h_su;
out.dp = max(0,min(param.funMdot_dp(M_dot), out.P_su - 2*CoolProp.PropsSI('pmin','Q',0,'H',1e5,fluid)));
out.P_ex = out.P_su-out.dp;
out.Q_dot = param.AU*(CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_su,fluid)-T_amb);
out.h_ex = h_su-out.Q_dot/M_dot;

end

function res = MdotDP_res(x, ub, h_ex, P_ex, M_dot, T_amb, fluid, param)
x = x.*ub;
out = MdotDP(x(1), x(2), M_dot, T_amb, fluid, param);
res(1) = h_ex-out.h_ex;
res(2) = P_ex-out.P_ex;
end

function out = PhiDP(h_su, P_su, M_dot, T_amb, fluid, param)
out.P_su = min(param.P_max,P_su);
out.h_su = h_su;
if strcmp(param.type_phi, 'phi_su')
    phi = M_dot^2/CoolProp.PropsSI('D','P',out.P_su,'H',out.h_su,fluid);
elseif strcmp(param.type_phi, 'phi_1')
    phi = M_dot^2/CoolProp.PropsSI('D','P',out.P_su,'Q',1,fluid);
elseif strcmp(param.type_phi, 'phi_0')
    phi = M_dot^2/CoolProp.PropsSI('D','P',out.P_su,'Q',0,fluid);
end
out.dp = max(0,min(param.funPhi_dp(phi), out.P_su - 1e3));
out.P_ex = out.P_su-out.dp;
out.Q_dot = param.AU*(CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_su,fluid)-T_amb);
out.h_ex = h_su-out.Q_dot/M_dot;
end

function res = PhiDP_res(x, ub, h_ex, P_ex, M_dot, T_amb, fluid, param)
x = x.*ub;
out = PhiDP(x(1), x(2), M_dot, T_amb, fluid, param);
res(1) = (h_ex-out.h_ex)/h_ex;
res(2) = (P_ex-out.P_ex)/P_ex;
end

