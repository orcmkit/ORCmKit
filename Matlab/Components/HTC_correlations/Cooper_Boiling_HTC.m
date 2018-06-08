function h_conv_boiling = Cooper_Boiling_HTC(P_c_su, fluid_c, h_conv_h, info, DT_log, Qdot)
C = 1;
AU_tp = Qdot/DT_log;
p_star = P_c_su/CoolProp.PropsSI('Pcrit', 'Q', 1, 'P', P_c_su, fluid_c);
M = 1e3*CoolProp.PropsSI('M', 'Q', 1, 'P', P_c_su, fluid_c);
Rp = 0.4; % roughness in µm
q = Qdot/info.A_c_tot;
err_q = 1;
k = 0;
while k <= 10 && err_q > 5e-2 %iterate for boiling number
    k = k+1;
    h = C*55*(p_star^(0.12-0.2*log10(Rp)))*((-log10(p_star))^(-0.55))*(q^0.67)*(M^(-0.5));
    U = (1/h +  1/h_conv_h)^-1;
    A_tp = AU_tp/U;
    q_new = Qdot/A_tp;
    err_q = abs(q_new-q)/q;
    q = q_new;
end

if err_q > 5e-2
    display('Cooper boiling: Wrong heat flux')
end
h_conv_boiling = h;

end