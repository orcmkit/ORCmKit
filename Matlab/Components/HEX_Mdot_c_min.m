function m_dot_c_min = HEX_Mdot_c_min(fluid_h, m_dot_h, P_h, in_h_su, Q_dot, fluid_c, P_c, in_c_su, param)
%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 27/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% HEX_Qdotmax is a single matlab code aiming to calculate the maxium amount of
% heat power that can be transferred between two fluids in a counterflow
% heat exchanger. This maxium is either given by a pinch point of 0K between 
% the temperature profiles, or limited by the maximum and minimum temperatures 
% achievable by the fluids.
% 
% The model inputs are:
%       - fluid_h: nature of the hot fluid                        	[-]
%       - P_h: inlet pressure of the hot fluid                      [Pa]
%       - in_h_su: inlet temperature or enthalpy of the hot fluid   [K or J/kg]
%       - m_dot_h: mass flow rate of the hot fluid                  [kg/s]
%       - fluid_c: nature of the cold fluid                        	[-]
%       - P_c: inlet pressure of the cold fluid                  [Pa]
%       - in_c_su: inlet temperature or enthalpy of the cold fluid  [K or J/kg]
%       - m_dot_c: mass flow rate of the cold fluid                 [kg/s]
%       - param: structure variable containing the model parameters i.e.
%           *param.type_h = type of input for hot fluid ('H' for enthalpy,'T' for temperature)
%           *param.type_c = type of input for cold fluid ('H' for enthalpy,'T' for temperature)
%
% The model outputs are:
%       - Q_dot_max : the maximum amount of power that can be transferred.
%          	
% See the documentation for further details or contact rdickes@ulg.ac.be

%% MODELLING CODE
res_T = 1e-2;
h_c_su = in_c_su;
H_h_vec = H_vec_definition(fluid_h, m_dot_h, P_h, in_h_su, Q_dot, 'hot');
T_h_vec = NaN*ones(1,length(H_h_vec));
for k = 1:length(H_h_vec)
    T_h_vec(k) = CoolProp.PropsSI('T','P', P_h, 'H', H_h_vec(k), fluid_h);
end
j = length(H_h_vec);
stop = 0;
while not(stop) && j >1
    Q_dot_test = m_dot_h*(H_h_vec(j)-H_h_vec(1));
    m_dot_c_test = Q_dot_test/(CoolProp.PropsSI('H','P', P_c, 'T', T_h_vec(j), fluid_c)-h_c_su);
    out_test = HEX_profile_2(fluid_h, m_dot_h, P_h, in_h_su, fluid_c, m_dot_c_test, P_c, h_c_su, Q_dot, param);
    if out_test.pinch < -res_T
        stop = 0;
    else
        stop = 1;
    end
    j = j-1;
end
m_dot_c_min = m_dot_c_test;
end