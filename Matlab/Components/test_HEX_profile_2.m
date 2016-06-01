fluid_h = 'PiroblocBasic';
fluid_c = 'R245fa';
param1.type_h = 'T';
param1.type_c  ='H';

m_dot_h = 0.1;
m_dot_c = 0.5;
P_h = 1e5;
P_c = 10e5;

Q_dot = 15e3;

%Cas 1:
param1.port_h = 'su';
param1.port_c = 'su';
in_h = 150+273.15; %CoolProp.PropsSI('H', 'T', 100+273.15, 'P', P_h, fluid_h);
in_c = CoolProp.PropsSI('H', 'T', 20+273.15, 'P', P_c, fluid_c);
out1 = HEX_profile_2(fluid_h, m_dot_h, P_h, in_h, fluid_c, m_dot_c, P_c, in_c, Q_dot, param1)

%Cas 2:
param2 = param1;
param2.port_h = 'ex';
param2.port_c = 'su';
out2 = HEX_profile_2(fluid_h, m_dot_h, P_h, out1.T_h_vec(1), fluid_c, m_dot_c, P_c, in_c, Q_dot, param2)

%Cas 3:
param3 = param1;
param3.port_h = 'su';
param3.port_c = 'ex';
out3 = HEX_profile_2(fluid_h, m_dot_h, P_h, in_h, fluid_c, m_dot_c, P_c, out1.H_c_vec(end), Q_dot, param3)

%Cas 4:
param4 = param1;
param4.port_h = 'ex';
param4.port_c = 'ex';
out4 = HEX_profile_2(fluid_h, m_dot_h, P_h, out1.T_h_vec(1), fluid_c, m_dot_c, P_c, out1.H_c_vec(end), Q_dot, param4)


figure;
hold on
plot(out1.x_vec, out1.T_h_vec-273.15, 'r--')
plot(out1.x_vec, out1.T_c_vec-273.15, 'b--')

plot(out2.x_vec, out2.T_h_vec-273.15, 'r:')
plot(out2.x_vec, out2.T_c_vec-273.15, 'b:')

plot(out3.x_vec, out3.T_h_vec-273.15, 'r-.')
plot(out3.x_vec, out3.T_c_vec-273.15, 'b-.')

plot(out4.x_vec, out4.T_h_vec-273.15, 'r-')
plot(out4.x_vec, out4.T_c_vec-273.15, 'b-')

hold off