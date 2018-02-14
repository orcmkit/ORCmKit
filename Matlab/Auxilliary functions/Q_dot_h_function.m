function out = Q_dot_h_function(input_vec, ub, fluid, type)

input_vec = input_vec.*ub;
m_dot = input_vec(1);
T_su = input_vec(2);
P_su = input_vec(3);
T_ex = input_vec(4);
P_ex = input_vec(5);

if abs(T_su - (CoolProp.PropsSI('T', 'Q', 0.5, 'P', P_su, fluid)-273.15))<1e-3
    if strcmp(type, 'Qdot')
        h_su = CoolProp.PropsSI('H', 'Q', 0, 'P', P_su, fluid);
    elseif strcmp(type, 'inv_Qdot')
        h_su = CoolProp.PropsSI('H', 'Q', 1, 'P', P_su, fluid);
    end
else
    h_su = CoolProp.PropsSI('H', 'T', T_su + 273.15, 'P', P_su, fluid);
end

if abs(T_ex - (CoolProp.PropsSI('T', 'Q', 0.5, 'P', P_ex, fluid)-273.15))<1e-3
    if strcmp(type, 'Qdot')
        h_ex = CoolProp.PropsSI('H', 'Q', 1, 'P', P_ex, fluid);
    elseif strcmp(type, 'inv_Qdot')
        h_ex = CoolProp.PropsSI('H', 'Q', 0, 'P', P_ex, fluid);
    end
else
    h_ex = CoolProp.PropsSI('H', 'T', T_ex + 273.15, 'P', P_ex, fluid);
end
Q_dot_h = m_dot*(h_su-h_ex);

if strcmp(type, 'Qdot')
    out = Q_dot_h;
elseif strcmp(type, 'inv_Qdot')
    out = 1/Q_dot_h;
end

end

