function out = eta_is_exp(input_vec, ub, fluid, type)

input_vec = input_vec.*ub;
m_dot = input_vec(1);
T_su = input_vec(2);
P_su = input_vec(3);
P_ex = input_vec(4);
Wdot_elec = input_vec(5);

if abs(T_su - (CoolProp.PropsSI('T', 'Q', 0.5, 'P', P_su, fluid)-273.15))<1e-3
    if strcmp(type, 'eta')
        h_su = CoolProp.PropsSI('H', 'Q', 0, 'P', P_su, fluid);
    elseif strcmp(type, 'inv_eta')
        h_su = CoolProp.PropsSI('H', 'Q', 1, 'P', P_su, fluid);
    end
else
    h_su = CoolProp.PropsSI('H', 'T', T_su + 273.15, 'P', P_su, fluid);
end
s_su = CoolProp.PropsSI('S', 'H', h_su, 'P', P_su, fluid);
h_ex_is = CoolProp.PropsSI('H', 'S', s_su, 'P', P_ex, fluid);
W_dot_is = m_dot*(h_su-h_ex_is);
eta_is = Wdot_elec/W_dot_is;

if strcmp(type, 'eta')
    out = eta_is;
elseif strcmp(type, 'inv_eta')
    out = 1/eta_is;
end

end

