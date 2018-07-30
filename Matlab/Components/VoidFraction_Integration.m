function alpha_mean = VoidFraction_Integration(q1, q2, rho_v, rho_l, param)

switch param.correlation.type_void_fraction
    case 'Homogenous'
        f_void = @(q) VoidFraction_homogenous(q, rho_v,  rho_l);
               
    case 'Zivi'
        f_void = @(q) VoidFraction_Zivi(q, rho_v,  rho_l);
                
    case 'SlipRatio'
        f_void = @(q) VoidFraction_SlipRatio(q, rho_v,  rho_l, param.S_ratio);
        
    case 'LM'
        f_void = @(q) VoidFraction_Lockhart_Martinelli(q, rho_v,  rho_l, param.mu_v, param.mu_l);
           
    case 'Premoli'
        f_void = @(q) VoidFraction_Premoli(q, rho_v, rho_l, param.mu_l, param.sig, param.Dh, param.G);
        
    case 'Hughmark'
        f_void = @(q) VoidFraction_Hughmark(q, rho_v, rho_l, param.mu_v, param.mu_l, param.Dh, param.G);

end
if abs(q2-q1) > 0.01
    alpha_mean = 1/(q2-q1)*integral(f_void, q1, q2);
else
    alpha_mean = f_void(0.5*q1+0.5*q2);
end
end








