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
elseif strcmp(fluid, 'OilLub')
    cp= max(1e-10,1701.03390+4.34455206e-3*(T-273.15));
    rho = max(1e-10,9.64196328e2-6.21731235e-1*(T-273.15));
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
end
end
