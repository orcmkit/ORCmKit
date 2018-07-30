function [M_mix_tot, M_lub_tot, M_ref_tot, M_mix_vec, M_lub_vec, M_ref_vec, out] = CycleMassProfiles(T_vec, P_vec, V_vec, C_oil_vec, fluid_ref, fluid_lub, param)

      
for j = 1:length(T_vec)
    
    outb = Mass_function_solub( T_vec(j), P_vec(j),  V_vec(j), C_oil_vec(j), fluid_ref, fluid_lub, param(j));
    out(j,1) = orderfields(outb);
    M_mix_vec(j,1) = out(j,1).M_mix;
    M_ref_vec(j,1) = out(j,1).M_ref;
    M_lub_vec(j,1) = out(j,1).M_lub;
    M_vap_vec(j,1) = out(j,1).M_vap;
    M_liq_vec(j,1) = out(j,1).M_liq;
    alpha_vec(j,1) = out(j,1).alpha;
    rho_liq_vec(j,1) = out(j,1).rho_liq;
    rho_vap_vec(j,1) = out(j,1).rho_vap;
    C_rv_vec(j,1) = out(j,1).C_rv;
    C_rl_vec(j,1) = out(j,1).C_rl;
    zeta_oil_vec(j,1) = out(j,1).zeta_oil;
       
end

M_mix_tot = sum(M_mix_vec);
M_lub_tot = sum(M_lub_vec);
M_ref_tot = sum(M_ref_vec);


end