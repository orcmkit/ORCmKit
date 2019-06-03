function [mu, Pr, k, mu_rat] = PropsSinglePhase_Fast(T_mean, P_mean, T_wall, T_sat_mean, zone, fluid, param)


if strcmp(fluid(1:3), 'ICP')
    k =  PropsSI_ICP('L', 'T', T_mean, 'P', P_mean, fluid);
    mu = PropsSI_ICP('V', 'T', T_mean, 'P', P_mean, fluid);
    cp = PropsSI_ICP('C', 'T', T_mean, 'P', P_mean, fluid);
    Pr = cp*mu/k;
    mu_wall = PropsSI_ICP('V', 'T', T_wall, 'P', P_mean, fluid);
    mu_rat = mu/mu_wall;
else
    if (strcmp(zone, 'liq')) && (T_mean >= T_sat_mean-5e-2)
        [cp, k, mu, ~] = CoolPropt_Props_QT(param.abs_lowLevel,param.CP_file, 0, T_mean);
        Pr = cp*mu/k;
        mu_rat = 1;
    elseif (strcmp(zone, 'vap') || strcmp(zone, 'vap_wet') || strcmp(zone, 'tp_dryout') ) && (T_mean <= T_sat_mean +5e-2)
        [cp, k, mu, ~] = CoolPropt_Props_QT(param.abs_lowLevel,param.CP_file, 1, T_mean);
        Pr = cp*mu/k;
        mu_rat = 1;       
    else
%         try
            [cp, k, mu, ~] = CoolPropt_Props_PT(param.abs_lowLevel,param.CP_file, P_mean, T_mean);
            [~, ~, mu_wall, ~] = CoolPropt_Props_PT(param.abs_lowLevel,param.CP_file, P_mean, T_wall);
            Pr = cp*mu/k;
            mu_rat = mu/mu_wall;
%         catch
%             mu = refpropm('V',          'T', T_mean,    'P', P_mean/1e3, fluid);
%             cp = refpropm('C',          'T', T_mean,    'P', P_mean/1e3, fluid);
%             k  = refpropm('L',          'T', T_mean,    'P', P_mean/1e3, fluid);
%             mu_wall = refpropm('V', 	'T', T_wall ,   'P', P_mean/1e3, fluid);
%             Pr = cp*mu/k;
%             mu_rat = mu/mu_wall;
%         end
    end
    
end
