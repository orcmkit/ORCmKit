function [P_su_q_ext, Q_su_q_ext, gamma_q_ext, P_su_T_ext, T_su_T_ext, gamma_T_ext] = gamma_map_generator


% fluid studied
fluid = 'R134a';

% saturated conditions
q_vec = linspace(0.1, 1, 10);
P_su_vec = linspace(1e5, 0.99*CoolProp.PropsSI('Pcrit','Q', 0, 'P',1e5,fluid), 20);
i = 0;
for P = P_su_vec
    for q = q_vec
        i = i + 1;
        P_su_q(i) = P;
        Q_su_q(i) = q;
        s_su_q(i) = CoolProp.PropsSI('S','Q', Q_su_q(i), 'P',P_su_q(i),fluid);
        rho_su_q(i) = CoolProp.PropsSI('D','Q', Q_su_q(i), 'P',P_su_q(i),fluid);
        f_gamma = @(x) gamma_fct(x, P_su_q(i), s_su_q(i), rho_su_q(i), fluid);
        lb = 0.1;
        ub = 2;
        gamma_test = zeroBrent ( lb, ub, 1e-6, 1e-6, f_gamma );
        if abs(f_gamma(gamma_test)) <1e-6
            gamma_q(i) = gamma_test;
        else
            gamma_q(i) = NaN;
        end
    end
end


gamma_q_ext = gamma_q;
Q_su_q_ext = Q_su_q;
P_su_q_ext = P_su_q./1e5;

%superheated conditions
P_su_vec = linspace(1e5, 0.99*CoolProp.PropsSI('Pcrit','Q', 0, 'P',1e5,fluid), 20);
i = 0;
for P = P_su_vec
    T_sat = CoolProp.PropsSI('T','Q', 0, 'P',P,fluid);
    for T = linspace(T_sat+0.1, 210+275.15, 20)
        i = i + 1;
        P_su_T(i) = P;
        T_su_T(i) = T;
        s_su_T(i) = CoolProp.PropsSI('S','T', T_su_T(i), 'P',P_su_T(i),fluid);
        rho_su_T(i) = CoolProp.PropsSI('D','T', T_su_T(i), 'P',P_su_T(i),fluid);
        f_gamma = @(x) gamma_fct(x, P_su_T(i), s_su_T(i), rho_su_T(i), fluid);
        lb = 0.1;
        ub = 2;
        gamma_test = zeroBrent ( lb, ub, 1e-6, 1e-6, f_gamma );
        if abs(f_gamma(gamma_test)) <1e-6
            gamma_T(i) = gamma_test;
        else
            gamma_T(i) = NaN;
        end
    end
end

%super critical conditions
P_su_vec = linspace(1.01*CoolProp.PropsSI('Pcrit','Q', 0, 'P',1e5,fluid), 55e5, 20);
for P = P_su_vec
    T_crit = CoolProp.PropsSI('Tcrit','Q', 0, 'P',P,fluid);
    for T = linspace(T_crit+10, 210+275.15, 20)
        i = i + 1;
        P_su_T(i) = P;
        T_su_T(i) = T;
        s_su_T(i) = CoolProp.PropsSI('S','T', T_su_T(i), 'P',P_su_T(i),fluid);
        rho_su_T(i) = CoolProp.PropsSI('D','T', T_su_T(i), 'P',P_su_T(i),fluid);
        f_gamma = @(x) gamma_fct(x, P_su_T(i), s_su_T(i), rho_su_T(i), fluid);
        lb = 0.1;
        ub = 2;
        gamma_test = zeroBrent ( lb, ub, 1e-6, 1e-6, f_gamma );
        if abs(f_gamma(gamma_test)) <1e-6
            gamma_T(i) = gamma_test;
        else
            gamma_T(i) = NaN;
        end
    end
end


gamma_T_ext = gamma_T;
T_su_T_ext = T_su_T./1e2;
P_su_T_ext = P_su_T./1e5;


end


function res = gamma_fct(x, P_su, s_su, rho_su, fluid)
gamma = max(0.1,min(x,3));
P_thr = P_su *(2/(gamma+1))^(gamma/(gamma-1));
rho_thr = CoolProp.PropsSI('D','P',P_thr,'S',s_su,fluid);
gamma_bis =  log10(P_su/P_thr)/log10(rho_su/rho_thr);
res = gamma_bis-gamma;
end