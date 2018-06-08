function out = W_dot_exp_hyd_function(input_vec, ub, fluid, type)

input_vec = input_vec.*ub;
m_dot = input_vec(1);
T_su = input_vec(2);
P_su = max(1e4,input_vec(3));
T_ex = input_vec(4);
P_ex = max(1e4,P_su-input_vec(5));

if abs(T_su - (CoolProp.PropsSI('T', 'Q', 0.5, 'P', P_su, fluid)-273.15))<1e-3
    if strcmp(type, 'Wdot')
        h_su = CoolProp.PropsSI('H', 'Q', 0, 'P', P_su, fluid);
    elseif strcmp(type, 'inv_Wdot')
        h_su = CoolProp.PropsSI('H', 'Q', 1, 'P', P_su, fluid);
    end
else
    h_su = CoolProp.PropsSI('H', 'T', T_su + 273.15, 'P', P_su, fluid);
end

if abs(T_ex - (CoolProp.PropsSI('T', 'Q', 0.5, 'P', P_ex, fluid)-273.15))<1e-3
    if strcmp(type, 'Wdot')
        h_ex = CoolProp.PropsSI('H', 'Q', 1, 'P', P_ex, fluid);
    elseif strcmp(type, 'inv_Wdot')
        h_ex = CoolProp.PropsSI('H', 'Q', 0, 'P', P_ex, fluid);
    end
else
    h_ex = CoolProp.PropsSI('H', 'T', T_ex + 273.15, 'P', P_ex, fluid);
end

W_dot_hyd = m_dot*(h_su-h_ex);

if strcmp(type, 'Wdot')
    out = W_dot_hyd;
elseif strcmp(type, 'inv_Wdot')
    out = 1/W_dot_hyd;
end

end

