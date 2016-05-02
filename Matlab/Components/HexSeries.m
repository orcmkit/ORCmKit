function out = HexSeries(fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, in_hex1, in_hex2)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 28/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% HexSeries is a single matlab function used to evaluate the heat transfer 
% occuring in two heat exchangers in series (see the Documentation/HexSeries_MatlabDoc)
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
%       - in_hex1: structure variable describing the first HEX (cfr HexModel.mat)
%       - in_hex2: structure variable describing the second HEX (cfr HexModel.mat)
%
% The model outputs is:
%       - out: a structure variable which includes all the modelling
%       results of the two heat exchanger (see HexMode.mat)
%
% See the documentation for further details or contact rdickes@ulg.ac.be

%% DEMONSTRATION CASE
if nargin == 0
    fluid_h = 'R245fa';
    P_h_su = 3e5;
    in_h_su = CoolProp.PropsSI('H', 'T', 90+273.15, 'P', P_h_su, fluid_h);
    m_dot_h = 0.1;
    fluid_c = 'INCOMP::MPG-30%';
    P_c_su = 2e5;
    in_c_su = CoolProp.PropsSI('H', 'T', 20+273.15, 'P', P_c_su, fluid_c);
    m_dot_c = 2.2;
    in_hex1.modelType = 'CstEff';
    in_hex1.type_h = 'H';
    in_hex1.type_c = 'H';
    in_hex1.epsilon_th = 0.1;
    in_hex1.displayResults = 0;
    in_hex1.displayTS = 0;
    in_hex1.advancedUser = 1;
    in_hex2.modelType = 'CstEff';
    in_hex2.type_h = 'H';
    in_hex2.type_c = 'H';
    in_hex2.epsilon_th = 0.65;
    in_hex2.displayResults = 0;
    in_hex2.displayTS = 1;
    in_hex2.advancedUser = 1;
end

%% MODELING OF 2 HEX IN SERIES 
Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, in_hex1);
lb = 0;
ub = Q_dot_max;
f = @(x) hexseries_res(x, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, in_hex1, in_hex2);
Q_dot_eff = zeroBrent ( lb, ub, 1e-6, 1e-6, f );
out = hexseries(Q_dot_eff, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, in_hex1, in_hex2);
if abs(out.res) < 1e-4
    out.flag = 1;
else
    out.flag = -1;
end

end

%% NESTED FUNCTIONS
function res = hexseries_res(x, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, in_hex1, in_hex2)
var = hexseries(x, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, in_hex1, in_hex2);
res = var.res;
end

function out = hexseries(x, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, in_hex1, in_hex2)
out.in_c_mid = in_c_su + x/m_dot_c;
[out_hex2, TS_hex2] = HexModel(fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, out.in_c_mid, m_dot_c, in_hex2);
out.hex2 = out_hex2;
out.ts2 = TS_hex2;
out.in_h_mid = out_hex2.h_h_ex;
[out_hex1, TS_hex1] = HexModel(fluid_h, P_h_su, out.in_h_mid, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, in_hex1);
out.hex1 = out_hex1;
out.ts1 = TS_hex1;

if x == 0 && out_hex1.Q_dot_tot == 0
    out.res = 0;
elseif x == 0 && out_hex1.Q_dot_tot ~= 0
    out.res = 1;
elseif x ~= 0 && out_hex1.Q_dot_tot == 0
    out.res = 1;
else
    out.res = 1 - x/out_hex1.Q_dot_tot;
end

end