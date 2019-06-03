function [h_l, h_v, T_su] = supply_conditions_HEX_Solub(in_su, P_su, fluid, param)

if strcmp(param.type,'H')
    if strcmp(fluid(1:3), 'ICP')
        T_su = PropsSI_ICP('T', 'H', in_su, 'P', P_su, fluid);
        h_l = inf;
        h_v = inf;
    else
        if param.solub
            T_sat_pure = CoolProp.PropsSI('T', 'P', P_su, 'Q', 0, fluid);
            T_min_bubble = R245fa_POE_Tbubble(1-param.C_oil, P_su, T_sat_pure);
            %[T_su, ~, ~, ~, ~, ~, ~, ~] = HP_solubMixt(param.C_oil, P_su, in_su, fluid, param.fluid_lub, T_min_bubble, T_sat_pure);
            [T_su, ~, ~, ~, ~, ~, ~, ~] = HP_solubMixt(param.C_oil, P_su, in_su, fluid, param.fluid_lub, T_min_bubble, T_sat_pure, param.fit_DTP_zeta);

            h_l = (1-param.C_oil)*CoolProp.PropsSI('H', 'T', T_min_bubble, 'Q', 0, fluid) + param.C_oil*PropsSI_ICP('H', 'T', T_min_bubble, 'P', P_su, param.fluid_lub);   
            if R245fa_POE_Tbubble((1- param.x_lim_vap - param.C_oil + param.x_lim_vap*param.C_oil)/(1- param.x_lim_vap + param.x_lim_vap*param.C_oil), P_su, T_sat_pure) > (155+273.15)
                h_v = inf;
            else
                [~, h_v, ~, ~, ~, ~, ~, ~] = xP_solubMixt(param.C_oil, P_su, param.x_lim_vap, fluid, param.fluid_lub, T_min_bubble, T_sat_pure, param.fit_DTP_zeta);
            end
        else
            T_su = CoolProp.PropsSI('T','H', in_su, 'P', P_su, fluid);
            h_l = CoolProp.PropsSI('H','P',P_su,'Q',0,fluid);
            h_v = CoolProp.PropsSI('H','P',P_su,'Q',1,fluid);
        end
    end
elseif strcmp(param.type,'T')
    T_su = in_su;
    h_l = inf;
    h_v = inf;
end

end