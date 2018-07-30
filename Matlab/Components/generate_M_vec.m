function [M_vec, M_mf_vec, M_lub_vec, alpha_mean_vec, M_tot, M_mf_tot, M_lub_tot] = generate_M_vec(param, out, fluid)

[M_vec, M_mf_vec, M_lub_vec, alpha_mean_vec] = deal(NaN*ones(1, length(out.V_vec)));
    
for i = 1:length(out.V_vec) 
    if param.solub %if solubility
        [~, ~, M_lub_vec(i), M_mf_vec(i), ~, ~, M_vec(i), alpha_mean_vec(i)] = MassMeanInt_solub(out.T_vec(i), out.P_vec(i), out.zeta_r_vec(i), out.C_rv_vec(i), out.Tbubble_min_vec(i), out.Tsat_pure_vec(i), out.T_vec(i+1), out.P_vec(i+1), out.zeta_r_vec(i+1), out.C_rv_vec(i+1), out.Tbubble_min_vec(i+1), out.Tsat_pure_vec(i+1), out.V_vec(i), fluid, param.fluid_lub, param.fit_ratio_rho, param);

    else % if no solubility
        if strcmp(param.type,'T') || strcmp(fluid(1:3),'ICP')
            M_mf_vec(i) = out.V_vec(i)*(0.5*PropsSI_ICP('D', 'T', out.T_vec(i), 'P', out.P_vec(i), fluid) + PropsSI_ICP('D', 'T', out.T_vec(i+1), 'P', out.P_vec(i+1), fluid));
            M_lub_vec(i) = NaN;
            M_vec(i) = M_mf_vec(i);
        else % if compressible
            if strcmp(out.type_zone{i},  'tp')
                q1 = CoolProp.PropsSI('Q','P',out.P_vec(i),'H',out.H_vec(i),fluid);
                q2 = CoolProp.PropsSI('Q','P',out.P_vec(i+1),'H',out.H_vec(i+1),fluid);
                rho_v = CoolProp.PropsSI('D','P',0.5*out.P_vec(i)+0.5*out.P_vec(i+1),'Q',1,fluid);
                rho_l = CoolProp.PropsSI('D','P',0.5*out.P_vec(i)+0.5*out.P_vec(i+1),'Q',0,fluid);
                out.alpha_mean_vec(i) = VoidFraction_Integration(q1, q2, rho_v, rho_l, param);
                out.M_mf_vec(i) =  out.V_vec(i)*(rho_v*out.alpha_mean_vec(i)+ rho_l*(1-out.alpha_mean_vec(i)));
            else
                out.M_mf_vec(i) = out.V_vec(i)*(0.5*CoolProp.PropsSI('D','P',out.P_vec(i),'H',out.H_vec(i),fluid) + 0.5*CoolProp.PropsSI('D','P',out.P_vec(i+1),'H',out.H_vec(i+1),fluid));
            end
            M_lub_vec(i) = NaN;
            M_vec(i) = M_mf_vec(i);
        end
    end
end

M_tot = sum(M_vec);
M_lub_tot = sum(M_lub_vec);
M_mf_tot = sum(M_mf_vec);

end