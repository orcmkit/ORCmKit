fluid = 'R245fa';
P_crit = 0.99*CoolProp.PropsSI('Pcrit', 'T', 0, 'P', 0, fluid);
P_min = CoolProp.PropsSI('P', 'T', 273.15-10, 'Q', 0, fluid);
P_lin = linspace(P_min, P_crit, 100);
for i=1:100
    s_liq(1,i) = CoolProp.PropsSI('S', 'P', P_lin(i), 'Q', 0, fluid);
    T_liq(1,i) = CoolProp.PropsSI('T', 'P', P_lin(i), 'Q', 0, fluid);
    s_vap(1,i) = CoolProp.PropsSI('S', 'P', P_lin(i), 'Q', 1, fluid);
    T_vap(1,i) = CoolProp.PropsSI('T', 'P', P_lin(i), 'Q', 1, fluid);
end

s_TS_curve = [s_liq flip(s_vap)];
T_TS_curve = [T_liq flip(T_vap)];