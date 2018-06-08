function alpha_mean = VoidFraction_Integration(q1, q2, rho_v, rho_l, param)

switch param.correlation.type_void_fraction
    case 'Homogenous'
        f_void = @(q) VoidFraction_homogenous(q, rho_v,  rho_l);
               
    case 'Zivi'
        f_void = @(q) VoidFraction_Zivi(q, rho_v,  rho_l);
                
    case 'SlipRatio'
        f_void = @(q) VoidFraction_SlipRatio(q, rho_v,  rho_l, param.S_ratio);
        
    case 'LockMart'
        f_void = @(q) VoidFraction_Lockhart_Martinelli(q, rho_v,  rho_l, param.mu_v, param.mu_l);
           
end

alpha_mean = 1/(q2-q1)*integral(f_void, q1, q2);
end








