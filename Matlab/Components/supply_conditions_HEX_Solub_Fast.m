function [h_l, h_v, T_su] = supply_conditions_HEX_Solub_Fast(in_su, P_su, fluid, param)

if strcmp(param.type,'H')
    if strcmp(fluid(1:3), 'ICP')
        T_su = PropsSI_ICP('T', 'H', in_su, 'P', P_su, fluid);
        h_l = inf;
        h_v = inf;
    else
        if param.solub
            T_sat_pure = CoolPropt_T_PQ(param.abs_lowLevel, param.CP_file, P_su, 0); %CoolProp.PropsSI('T', 'P', P_su, 'Q', 0, fluid);
            T_min_bubble = R245fa_POE_Tbubble(1-param.C_oil, P_su, T_sat_pure);
            [T_su, ~, ~, ~, ~, ~, ~, ~] = HP_solubMixt_Fast(param.C_oil, P_su, in_su, fluid, param.fluid_lub, T_min_bubble, T_sat_pure, param.fit_DTP_zeta, param);
            %h_l = (1-param.C_oil)*CoolProp.PropsSI('H', 'T', T_min_bubble, 'Q', 0, fluid) + param.C_oil*PropsSI_ICP('H', 'T', T_min_bubble, 'P', P_su, param.fluid_lub);   
            h_l = (1-param.C_oil)*CoolPropt_H_QT(param.abs_lowLevel,param.CP_file, 0, T_min_bubble) + param.C_oil*PropsSI_ICP('H', 'T', T_min_bubble, 'P', P_su, param.fluid_lub);               
            if R245fa_POE_Tbubble((1- param.x_lim_vap - param.C_oil + param.x_lim_vap*param.C_oil)/(1- param.x_lim_vap + param.x_lim_vap*param.C_oil), P_su, T_sat_pure) > (155+273.15)
                h_v = inf;
            else
                [~, h_v, ~, ~, ~, ~, ~, ~] = xP_solubMixt_Fast(param.C_oil, P_su, param.x_lim_vap, fluid, param.fluid_lub, T_min_bubble, T_sat_pure, param.fit_DTP_zeta, param);
            end
        else
            T_su = CoolPropt_T_HP(param.abs_lowLevel,param.CP_file, in_su, P_su); %CoolProp.PropsSI('T','H', in_su, 'P', P_su, fluid);
            h_l = CoolPropt_H_PQ(param.abs_lowLevel,param.CP_file, P_su, 0);%CoolProp.PropsSI('H','P',P_su,'Q',0,fluid);
            h_v = CoolPropt_H_PQ(param.abs_lowLevel,param.CP_file, P_su, 1);%CoolProp.PropsSI('H','P',P_su,'Q',1,fluid);
        end
    end
elseif strcmp(param.type,'T')
    T_su = in_su;
    h_l = inf;
    h_v = inf;
end

end