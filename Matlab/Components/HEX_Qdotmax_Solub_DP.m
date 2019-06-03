function [Q_dot_max, pinch_min] = HEX_Qdotmax_Solub_DP(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, DP_h, DP_c, info)
%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 25/04/2018 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% HEX_Qdotmax_Solub is a single matlab code aiming to calculate the maxium amount of
% heat power that can be transferred between two fluids in a counterflow
% heat exchanger, while accounting for possible lubricant solubility. 
% This maxium is either given by a pinch point of 0K between 
% the temperature profiles, or limited by the maximum and minimum temperatures 
% achievable by the fluids.
% 
% See the documentation for further details or contact rdickes@ulg.ac.be

%% MODELLING CODE
info.n_disc = 2;
if strcmp(info.H.type,'H') && strcmp(info.C.type,'H') % CASE 1 : HOT FLUID AND COLD FLUID MIGHT EXPERIENCE A PHASE CHANGE
    T_h_min = 253.15;
    T_c_max = 463.15;
    h_h_su = in_h_su;
    h_c_su = in_c_su;
    if info.C.solub
        Tsat_pure_c = CoolProp.PropsSI('T', 'P', P_c_su, 'Q', 0, fluid_c);
        Tbubble_min_c = R245fa_POE_Tbubble(1-info.C.C_oil, P_c_su, Tsat_pure_c);
        [T_c_su, ~, ~, ~, ~, ~, ~, ~] = HP_solubMixt(info.C.C_oil, P_c_su, h_c_su, fluid_c, info.C.fluid_lub, Tbubble_min_c, Tsat_pure_c, info.C.fit_DTP_zeta);
    else
        if strcmp(fluid_c(1:3), 'ICP')
            T_c_su = PropsSI_ICP('T', 'H', h_c_su, 'P', P_c_su, fluid_c);
        else
            T_c_su = CoolProp.PropsSI('T', 'H', h_c_su, 'P', P_c_su, fluid_c);
        end
    end % calcul T_c_su
    if info.H.solub
        Tsat_pure_h = CoolProp.PropsSI('T', 'P', P_h_su, 'Q', 0, fluid_h);
        Tbubble_min_h = R245fa_POE_Tbubble(1-info.H.C_oil, P_h_su, Tsat_pure_h);
        [T_h_su, ~, ~, ~, ~, ~, ~, ~] = HP_solubMixt(info.H.C_oil, P_h_su, h_h_su, fluid_h, info.H.fluid_lub, Tbubble_min_h, Tsat_pure_h, info.H.fit_DTP_zeta);
    else
        if strcmp(fluid_h(1:3), 'ICP')
            T_h_su = PropsSI_ICP('T', 'H', h_h_su, 'P', P_h_su, fluid_h);
        else
            T_h_su = CoolProp.PropsSI('T', 'H', h_h_su, 'P', P_h_su, fluid_h);
        end
    end % calcul T_h_su
    if info.C.solub
        [h_c_ex_extmax, ~, ~, ~, ~, ~, ~, ~] = TP_solubMixt(min(T_h_su, T_c_max), P_c_su, info.C.C_oil, fluid_c, info.C.fluid_lub, Tbubble_min_c, Tsat_pure_c, info.C.fit_DTP_zeta);
    else
        if strcmp(fluid_c(1:3), 'ICP')
            h_c_ex_extmax = PropsSI_ICP('H', 'T', min(T_h_su, T_c_max), 'P', P_c_su, fluid_c);
        else
            h_c_ex_extmax = CoolProp.PropsSI('H', 'T', min(T_h_su, T_c_max), 'P', P_c_su, fluid_c);
        end
    end % calcul h_c_ex_extmax
    if info.H.solub
        [h_h_ex_extmax, ~, ~, ~, ~, ~, ~, ~] = TP_solubMixt(max(T_c_su, T_h_min), P_h_su, info.H.C_oil, fluid_h, info.H.fluid_lub, Tbubble_min_h, Tsat_pure_h, info.H.fit_DTP_zeta);
    else
        if strcmp(fluid_h(1:3), 'ICP')
            h_h_ex_extmax = PropsSI_ICP('H', 'T', max(T_c_su, T_h_min), 'P', P_h_su, fluid_h);
        else
            h_h_ex_extmax = CoolProp.PropsSI('H', 'T', max(T_c_su, T_h_min), 'P', P_h_su, fluid_h);
        end
    end % calcul h_h_ex_extmax
    Q_dot_hext_max = m_dot_h*(h_h_su-h_h_ex_extmax);
    Q_dot_cext_max = m_dot_c*(h_c_ex_extmax-h_c_su);
    ub = min(Q_dot_cext_max, Q_dot_hext_max)*1.01;
    lb = 0;
    f = @(x) pinch0_Solub_DP(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, x, DP_h, DP_c,  info);
    [Q_dot_max,pinch_min] = zeroBrent ( lb, ub, 1e-16, 1e-14, f, 1e-8 );
end

end

function err = pinch0_Solub_DP(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, DP_h, DP_c,  info)
out = HEX_profile_Solub_DP(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, DP_h, DP_c , info);
err = out.pinch;
end