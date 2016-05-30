function res = FCT_Hex_hConvVar_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, param)

% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems
%
% Remi Dickes - 11/05/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% "FCT_Hex_hConvVar_res.m" is a constitutive subfunction of "HexModel_hConvVar.m"
% It isolates the residual "res" obtained by "FCT_Hex_hConvVar.m" while assuming a
% given heat transfer "Q_dot" in the heat exchanger

out = FCT_Hex_hConvVar(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, param);
res = out.resA;

end

