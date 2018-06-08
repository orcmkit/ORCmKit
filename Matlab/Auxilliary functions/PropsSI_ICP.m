function [out] = PropsSI_ICP(Prop, v1, in1, v2, in2, fluid)

if ~ischar(fluid)
    %cp  = cp_bar(fluid.m, fluid.n, fluid.f, T1-273.15, T2-273.15);
    %rho = NaN;
    %mu = NaN;
    %k = NaN;
elseif strcmp(fluid, 'ICP_PiroblocBasic')
    a0 = ((1.93403625e3) - (273.15*9.71458203e-1) + (273.15^2*8.24730587e-3));
    a1 = (9.71458203e-1 - 2*273.15*8.24730587e-3);
    a2 = 8.24730587e-3;
    if strcmp(v1, 'T')
        T = in1;
    elseif strcmp(v1, 'H')
        rt = roots([a2/3 a1/2 a0 -in1]);
        T  = max(rt(imag(rt)==0));
    end
    cp = a0 + a1*T +a2*T^2;
    h = a0*T + a1/2*T^2 + a2/3*T^3;
    rho = max(1e-10,8.83041605e2 - 6.65325683e-1*(T-273.15));
    mu = exp(-15.3 + 3526./T);
    k = max(1e-10,9.04426244e-1 - 2.34107821e-4*(T-273.15));
elseif strcmp(fluid, 'ICP_PiroblocBasic_2')
    a0 = 699.4141908; %ok
    a1 = 3.976263771; %ok
    if strcmp(v1, 'T')
       T = in1;
    elseif strcmp(v1, 'H')
       rt = roots([a1/2 a0 -in1]);
       T  = max(rt(imag(rt)==0));
    end
    cp = a0 + a1*T ; %ok
    h = a0*T + a1/2*T^2; %ok
    rho = max(1e-10,890.0666667 - 0.6565714286*(T-273.15)); %ok
    mu = max(1e-10,51.53590531*(T-273.15)^(-2.087778736)); %ok
    k = max(1e-10,0.1321943333 - 7.310285714e-05*(T-273.15));%ok
elseif strcmp(fluid, 'ICP_RL_32_3MAF')
    a0 = (1804.5 -273.15*2.3639);
    a1 = 2.3639;
    if strcmp(v1, 'T')
        T = in1;
    elseif strcmp(v1, 'H')
        rt = roots([a1/2 a0 -in1]);
        T  = max(rt(imag(rt)==0));
    end

    cp = a0 + a1*T;
    h = a0*T + a1/2*T^2; % integration of cp
    rho = 981 - 0.7547*(T-273.15 - 20);  %rho0 from datasheet, slope from other POE data
    mu = exp(-14.48 + 3439./T); % Andrate equation fitted on Emkarate brochure
    k = NaN;
    
%elseif strcmp(fluid, 'OilLub')
    %cp= max(1e-10,1701.03390+4.34455206e-3*(T-273.15));
    %rho = max(1e-10,9.64196328e2-6.21731235e-1*(T-273.15));
    %mu = NaN;
    %k = NaN;
%elseif strcmp(fluid, 'RL_32_H')
    %cp= NaN;
    %rho = 992.74 - 0.7547*(T-273.15);
    %mu = NaN;
    %k = NaN;

%elseif strcmp(fluid, 'Planetelf_ACD_100YF')
    %cp= NaN;
    %rho = 986.54 - 0.7314*(T-273.15);
    %mu = NaN;
    %k = NaN;
%else
    %cp = CoolProp.PropsSI('C', 'T', T, 'P',P, fluid);
    %rho = CoolProp.PropsSI('D', 'T', T, 'P',P, fluid);
    %mu = CoolProp.PropsSI('V', 'T', T, 'P',P, fluid);
    %k = CoolProp.PropsSI('L', 'T', T, 'P',P, fluid);
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
    case 'T'
        out = T;
end
end


%POE specific capacity:
% cp = 1804.5 + 2.3639*T ; %Gopfert - POE // T in °C, cp in J/kg.K 
% cp = 1304 + 1.035*T + 2.801e-3*T^2; % Zhelezny - Planetelf ACD 100 FY // T in K, cp in J/kg.K
% cp = 1000*(53.266 + 0.107*T)/sqrt(rho_15degC)