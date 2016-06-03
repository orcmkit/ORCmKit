function T_vec = T_vec_definition(fluid, m_dot, P, in_T, Q_dot, port_in)

if strcmp(port_in,'hot')
    T_h = in_T;
    f_T_c = @(x) Tc_def(x, T_h,  P, m_dot, Q_dot, fluid);
    lb =  T_h-Q_dot/m_dot/sf_PropsSI_bar('C', T_h, T_h, P, fluid)-30;
    ub = T_h;
    T_c = zeroBrent ( lb, ub, 1e-6, 1e-6, f_T_c );
elseif strcmp(port_in,'cold')
    T_c = in_T;
    f_T_h = @(x) Th_def(x, T_c,  P, m_dot, Q_dot, fluid);
    lb =  T_c;
    ub = T_c+Q_dot/m_dot/sf_PropsSI_bar('C', T_c, T_c, P, fluid)+30;
    T_h = zeroBrent ( lb, ub, 1e-6, 1e-6, f_T_h );
end

T_vec = [T_c T_h];

end

function err = Th_def(T_h, T_c,  P, Mdot, Qdot, fluid)
err = T_h-(T_c+Qdot/Mdot/sf_PropsSI_bar('C', T_c, T_h, P, fluid));
end

function err = Tc_def(T_c, T_h,  P, Mdot, Qdot, fluid)
err = T_c-(T_h-Qdot/Mdot/sf_PropsSI_bar('C', T_h, T_c, P, fluid));
end