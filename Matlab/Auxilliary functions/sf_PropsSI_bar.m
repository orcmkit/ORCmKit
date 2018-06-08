function [out] = sf_PropsSI_bar(Prop, T1, T2, P, fluid)
T = (T1+T2)/2;

if ~ischar(fluid)
    cp  = cp_bar(fluid.m, fluid.n, fluid.f, T1-273.15, T2-273.15);
    rho = NaN;
    mu = NaN;
    k = NaN;
elseif strcmp(fluid, 'PiroblocBasic')
    cp = max(1e-10,1.93403625e3 + 9.71458203e-1*(T-273.15)+8.24730587e-3*(T-273.15)^2);
    rho = max(1e-10,8.83041605e2 - 6.65325683e-1*(T-273.15));
    nu= max(1e-10,0.007023*(T-273.15)^(-1.585));
    mu = nu*rho;
    k = max(1e-10,9.04426244e-1 - 2.34107821e-4*(T-273.15));
elseif strcmp(fluid, 'PiroblocBasic_2')
    cp = max(1e-10, 699.4141908 + 3.976263771*T );
    rho = max(1e-10,890.0666667 - 0.6565714286*(T-273.15)); %ok
    mu = max(1e-10,51.53590531*(T-273.15)^(-2.087778736)); %ok
    k = max(1e-10,0.1321943333 - 7.310285714e-05*(T-273.15));%ok
    
elseif strcmp(fluid, 'OilLub')
    cp= max(1e-10,1701.03390+4.34455206e-3*(T-273.15));
    rho = max(1e-10,9.64196328e2-6.21731235e-1*(T-273.15));
    mu = NaN;
    k = NaN;
elseif strcmp(fluid, 'RL_32_H')
    cp= NaN;
    rho = 992.74 - 0.7547*(T-273.15);
    mu = NaN;
    k = NaN;
elseif strcmp(fluid, 'RL_32_3MAF')
    cp = 1804.5 + 2.3639*(T-273.15); % T in K
    h = 1804.5*(T-273.15) + 2.3639/2*(T-273.15)^2; % integration of cp
    rho = 981 - 0.7547*(T-273.15 - 20);  %rho0 from datasheet, slope from other POE data
    mu = exp(-14.48 + 3439./T); % Andrate equation fitted on Emkarate brochure
    k = NaN;
elseif strcmp(fluid, 'Planetelf_ACD_100YF')
    cp= NaN;
    rho = 986.54 - 0.7314*(T-273.15);
    mu = NaN;
    k = NaN;
else
    cp = CoolProp.PropsSI('C', 'T', T, 'P',P, fluid);
    rho = CoolProp.PropsSI('D', 'T', T, 'P',P, fluid);
    mu = CoolProp.PropsSI('V', 'T', T, 'P',P, fluid);
    k = CoolProp.PropsSI('L', 'T', T, 'P',P, fluid);
end

switch Prop
    case 'C'
        out = cp;
    case 'D'
        out = rho;
    case 'V'
        out = mu;
    case 'L'
        out = k;
        
    case 'H'
        out = h;
end
end


%POE specific capacity:
% cp = 1804.5 + 2.3639*T ; %Gopfert - POE // T in °C, cp in J/kg.K 
% cp = 1304 + 1.035*T + 2.801e-3*T^2; % Zhelezny - Planetelf ACD 100 FY // T in K, cp in J/kg.K
% cp = 1000*(53.266 + 0.107*T)/sqrt(rho_15degC)