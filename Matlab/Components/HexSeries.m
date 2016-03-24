function out = HexSeries(fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, in_hex1, in_hex2)

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
    in_hex2.displayTS = 0;
    in_hex2.advancedUser = 1;
end

%% HEX computation
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