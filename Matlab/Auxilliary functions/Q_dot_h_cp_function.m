function out = Q_dot_h_cp_function(input_vec, ub, fluid, type)
input_vec = input_vec.*ub;
V_dot = input_vec(1);
T_su = input_vec(2);
T_ex = input_vec(3);

cp = sf_PropsSI_bar('D', T_su+273.15, T_ex+273.15, 1e5, fluid);
rho = sf_PropsSI_bar('D', T_ex+273.15, T_ex+273.15, 1e5, fluid);
m_dot = V_dot*rho;
Q_dot_h = m_dot*cp*(T_su-T_ex);

if strcmp(type, 'Qdot')
    out = Q_dot_h;
elseif strcmp(type, 'inv_Qdot')
    out = 1/Q_dot_h;
end

end

