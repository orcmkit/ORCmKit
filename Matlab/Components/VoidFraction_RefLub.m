function alpha = VoidFraction_RefLub(q, rho_v, rho_l, param)

switch param.type_void_fraction
    case 'Homogenous'
        alpha = VoidFraction_homogenous(q, rho_v,  rho_l);
        
    case 'Zivi'
        alpha = VoidFraction_Zivi(q, rho_v,  rho_l);
        
    case 'SlipRatio'
        alpha = VoidFraction_SlipRatio(q, rho_v,  rho_l, param.S_ratio);
        
    case 'LockMart'
        alpha = VoidFraction_Lockhart_Martinelli(q, rho_v,  rho_l, param.mu_v, param.mu_l);

    case 'GrahamCondensation'
        alpha = VoidFraction_GrahamCondensation(q, rho_v, param.G, param.D);
    
    case 'Hughmark'
        alpha = VoidFraction_Hughmark(q, rho_v, rho_l, param.mu_v, param.mu_l, param.D, param.G);

    case 'AnnularFlow'
        if q > 0.7
            fit_ratiobis = param.fit_ratio_rho;
            AF = AnnularFlow_VPipe(param.G, param.P_Pa, param.T_K, param.C_oil, param.D, param.fluid_ref, param.fluid_lub, fit_ratiobis, param.fit_DTP_zeta);
            if AF.flag>0
                alpha = AF.alpha;
            else
                alpha = VoidFraction_SlipRatio(q, rho_v,  rho_l, 10);
                disp('Prob Annular flow')
            end
        else
            alpha = VoidFraction_SlipRatio(q, rho_v,  rho_l, 10);
        end
end
end

