function [Q_dot_max, pinch_min] = HEX_Qdotmax_Solub_DP2_Fast(fluid_h, m_dot_h, P_h_su, in_h_su, T_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, T_c_su, info, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, DP_h, DP_c)
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
%info.n_disc = 5;
if strcmp(info.H.type,'H') && strcmp(info.C.type,'H') % CASE 1 : HOT FLUID AND COLD FLUID MIGHT EXPERIENCE A PHASE CHANGE
    T_h_min = 253.15;
    T_c_max = 463.15;
    h_h_su = in_h_su;
    h_c_su = in_c_su;    
    if info.C.solub
        Tsat_pure_c = CoolPropt_T_PQ(info.C.abs_lowLevel,info.C.CP_file, P_c_su, 0); %CoolProp.PropsSI('T', 'P', P_c_su, 'Q', 0, fluid_c);
        Tbubble_min_c = R245fa_POE_Tbubble(1-info.C.C_oil, P_c_su, Tsat_pure_c);        
        [h_c_ex_extmax, ~, ~, ~, ~, ~, ~, ~] = TP_solubMixt_Fast(min(T_h_su, T_c_max), P_c_su, info.C.C_oil, fluid_c, info.C.fluid_lub, Tbubble_min_c, Tsat_pure_c, info.C.fit_DTP_zeta, info.C);
    else
        if strcmp(fluid_c(1:3), 'ICP')
            h_c_ex_extmax = PropsSI_ICP('H', 'T', min(T_h_su, T_c_max), 'P', P_c_su, fluid_c);
        else
            h_c_ex_extmax = CoolPropt_H_PT(info.C.abs_lowLevel,info.C.CP_file, P_c_su, min(T_h_su, T_c_max)); %CoolProp.PropsSI('H', 'T', , 'P', P_c_su, fluid_c);
        end
    end % calcul h_c_ex_extmax
    if info.H.solub
        Tsat_pure_h = CoolPropt_T_PQ(info.H.abs_lowLevel,info.H.CP_file, P_h_su, 0); %CoolProp.PropsSI('T', 'P', P_h_su, 'Q', 0, fluid_h);
        Tbubble_min_h = R245fa_POE_Tbubble(1-info.H.C_oil, P_h_su, Tsat_pure_h);
        [h_h_ex_extmax, ~, ~, ~, ~, ~, ~, ~] = TP_solubMixt_Fast(max(T_c_su, T_h_min), P_h_su, info.H.C_oil, fluid_h, info.H.fluid_lub, Tbubble_min_h, Tsat_pure_h, info.H.fit_DTP_zeta, info.H);
    else
        if strcmp(fluid_h(1:3), 'ICP')
            h_h_ex_extmax = PropsSI_ICP('H', 'T', max(T_c_su, T_h_min), 'P', P_h_su, fluid_h);
        else
            h_h_ex_extmax = CoolPropt_H_PT(info.H.abs_lowLevel,info.H.CP_file, P_h_su, max(T_c_su, T_h_min)); %CoolProp.PropsSI('H', 'T', max(T_c_su, T_h_min), 'P', P_h_su, fluid_h);
        end
    end % calcul h_h_ex_extmax
    Q_dot_hext_max = m_dot_h*(h_h_su-h_h_ex_extmax);
    Q_dot_cext_max = m_dot_c*(h_c_ex_extmax-h_c_su);
    ub = min(Q_dot_cext_max, Q_dot_hext_max)*1.01;
    lb = 0;
    f = @(x) pinch0_Solub_Fast(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, x, info, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, DP_h, DP_c);
    [Q_dot_max, pinch_min ]= zeroBrent ( lb, ub, 1e-12, 1e-8, f, 1e-5);
end

end

function err = pinch0_Solub_Fast(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, DP_h, DP_c)
out = HEX_profile_Solub_DP2_Fast(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, DP_h, DP_c);
err = out.pinch;
end