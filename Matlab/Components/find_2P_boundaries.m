function [h_l, h_v, P_l, P_v, flag_l_bd, flag_v_bd] = find_2P_boundaries(fluid, h_su, h_ex, P_su, P_ex, param)
if strcmp(fluid(1:3), 'ICP')
    h_l = inf;
    h_v = inf;
    P_l = NaN;
    P_v = NaN;
    flag_l_bd = 1;
    flag_v_bd = 1;
else
    if abs(P_su-P_ex)<1e1
        P_l = 0.5*P_su + 0.5*P_ex;
        P_v = 0.5*P_su + 0.5*P_ex;
        flag_l_bd = 1;
        flag_v_bd = 1;
    else
        lb = 0.999*min(P_ex, P_su);
        ub = 1.001*max(P_ex, P_su);
        p_o = 0.5*P_su + 0.5*P_ex;
        
        f_l = @(x) res_find_saturated_liquid_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, param);
        [x_P_l,res_P_l,flag_l_bd] = fminbnd(f_l, lb./ub,ub./ub);
        P_l = x_P_l.*ub;
        
        f_v = @(x) res_find_saturated_vapour_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, param);
        [x_P_v,res_P_v,flag_v_bd] = fminbnd(f_v, lb./ub,ub./ub);
        P_v = x_P_v.*ub;
    end
    h_l = saturated_liquid_enthalpy(P_l, param, fluid);
    h_v = saturated_vapour_enthalpy(P_v, param, fluid);
end
end

function res= res_find_saturated_liquid_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, param)
[res, ~, ~] = find_saturated_liquid_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, param);
end

function [res, P_l_guess, h_l] = find_saturated_liquid_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, param)
P_l_guess = x.*ub;
h_l = saturated_liquid_enthalpy(P_l_guess, param, fluid);
ratio_h = (h_su - h_l)/(h_su - h_ex);
P_l_bis = (1-ratio_h)*P_su + ratio_h*P_ex;
res = abs(P_l_guess - P_l_bis)/P_l_bis;
end

function h_l = saturated_liquid_enthalpy(P_l, param, fluid)
if param.solub
    Tbubble_min_l = R245fa_POE_Tbubble(1-param.C_oil, P_l, CoolProp.PropsSI('T', 'Q', 0, 'P', P_l, fluid));
    h_l = (1-param.C_oil)*CoolProp.PropsSI('H', 'T', Tbubble_min_l, 'Q', 0, fluid) + param.C_oil*PropsSI_ICP('H', 'T', Tbubble_min_l, 'P', P_l, param.fluid_lub);
else
    h_l = CoolProp.PropsSI('H', 'Q', 0, 'P', P_l, fluid);
end
end

function res= res_find_saturated_vapour_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, param)
[res, ~, ~] = find_saturated_vapour_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, param);
end

function [res, P_v_guess, h_v] = find_saturated_vapour_pressure(x ,ub, h_su, h_ex, P_su, P_ex, fluid, param)
P_v_guess = x.*ub;
h_v = saturated_vapour_enthalpy(P_v_guess, param, fluid);
ratio_h = (h_su - h_v)/(h_su - h_ex);
P_v_bis = (1-ratio_h)*P_su + ratio_h*P_ex;
res = abs(P_v_guess - P_v_bis)/P_v_bis;
end

function h_v = saturated_vapour_enthalpy(P_v, param, fluid)
if param.solub
    Tsat_pure_v = CoolProp.PropsSI('T', 'Q', 1, 'P', P_v, fluid);
    Tbubble_min_v = R245fa_POE_Tbubble(1-param.C_oil, P_v, Tsat_pure_v);
    [~, h_v, ~, ~, ~, ~, ~, ~] = xP_solubMixt(param.C_oil, P_v, param.x_lim_vap, fluid, param.fluid_lub, Tbubble_min_v, Tsat_pure_v, param.fit_DTP_zeta);
else
    h_v = CoolProp.PropsSI('H', 'Q', 1, 'P', P_v, fluid);
end
end