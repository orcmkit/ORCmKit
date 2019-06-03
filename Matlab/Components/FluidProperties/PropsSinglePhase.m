function [mu, Pr, k, mu_rat] = PropsSinglePhase(T_mean, P_mean, T_wall, T_sat_mean, zone, fluid)


if strcmp(fluid(1:3), 'ICP')
    k =  PropsSI_ICP('L', 'T', T_mean, 'P', P_mean, fluid);
    mu = PropsSI_ICP('V', 'T', T_mean, 'P', P_mean, fluid);
    cp = PropsSI_ICP('C', 'T', T_mean, 'P', P_mean, fluid);
    Pr = cp*mu/k;
    mu_wall = PropsSI_ICP('V', 'T', T_wall, 'P', P_mean, fluid);
    mu_rat = mu/mu_wall;
else
    if (strcmp(zone, 'liq')) && (T_mean >= T_sat_mean-5e-2)
        mu = CoolProp.PropsSI('V',        'T', T_mean, 'Q', 0, fluid);
        Pr = CoolProp.PropsSI('Prandtl',  'T', T_mean, 'Q', 0, fluid);
        k  = CoolProp.PropsSI('L',        'T', T_mean, 'Q', 0, fluid);
        mu_rat = 1;
    elseif (strcmp(zone, 'vap') || strcmp(zone, 'vap_wet') || strcmp(zone, 'tp_dryout') ) && (T_mean <= T_sat_mean +5e-2)
        mu = CoolProp.PropsSI('V',        'T', T_mean, 'Q', 1, fluid);
        Pr = CoolProp.PropsSI('Prandtl',  'T', T_mean, 'Q', 1, fluid);
        k  = CoolProp.PropsSI('L',        'T', T_mean, 'Q', 1, fluid);
        mu_rat = 1;
    else
        try
            mu = CoolProp.PropsSI('V',        'T', T_mean, 'P', P_mean, fluid);
            Pr = CoolProp.PropsSI('Prandtl',  'T', T_mean, 'P', P_mean, fluid);
            k  = CoolProp.PropsSI('L',        'T', T_mean, 'P', P_mean, fluid);
            mu_wall = CoolProp.PropsSI('V',  	'T', T_wall ,  'P', P_mean, fluid);
            mu_rat = mu/mu_wall;
        catch
            mu = refpropm('V',          'T', T_mean,    'P', P_mean/1e3, fluid);
            cp = refpropm('C',          'T', T_mean,    'P', P_mean/1e3, fluid);
            k  = refpropm('L',          'T', T_mean,    'P', P_mean/1e3, fluid);
            mu_wall = refpropm('V', 	'T', T_wall ,   'P', P_mean/1e3, fluid);
            Pr = cp*mu/k;
            mu_rat = mu/mu_wall;
        end
    end
    
end
