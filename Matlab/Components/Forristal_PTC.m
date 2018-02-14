function out = Forristal_PTC(fluid1, T1, P1, m_dot1, T6, D2, D3, D4, D5, L, W_ptc, fluid34, P34, P56, v_wind, eps3, eps4, eps5, alpha3, alpha5, tau5, eta_opt, k45, k23, DNI, cos_theta, GlassEnvelop)

% Forristall's model to simulate parabolic trough collectors
% implement by RDickes - 14/02/2018

if nargin ==0
    fluid1 = 'PiroblocBasic'; %-
    T1 = 150 + 273.15; %K
    P1 = 1e5; %Pa
    m_dot1 = 0.1; %kg/s
    T6 = 25+273.15; %K
    D2 = 0.035;%m
    D3 = 0.042;%m
    D4 = 0.05;%m
    D5 = 0.062;%m
    L = 1; %m
    W_ptc = 2.15;%m
    fluid34 = 'air'; %-
    P34 = 1e5; %Pa
    P56 = 1e5; %Pa
    v_wind = 1; %m/s
    eps3 = 0.3; %-
    eps4 = 0.3;%-
    eps5 = 0.3;%-
    alpha3 = 0.95; %-
    alpha5 = 1e-2;%-
    tau5 = 0.95; %-
    eta_opt = 0.8; %-
    k45 = 1.04; %W/m.K
    k23 = 17; %W/m.K
    DNI = 900; %W/m²
    cos_theta = 0.7; %-
    GlassEnvelop = 1; %-
end

% Resolution of energy balance
out.T_vec = solve_energy_balance(fluid1, T1, P1, m_dot1, T6, D2, D3, D4, D5, L, W_ptc, fluid34, P34, P56, v_wind, eps3, eps4, eps5, alpha3, alpha5, tau5, eta_opt, k45, k23, DNI, cos_theta, GlassEnvelop);


% Detailed information

if GlassEnvelop
    out.T1 = T1;
    out.T2 = out.T_vec(1);
    out.T3 = out.T_vec(2);
    out.T4 = out.T_vec(3);
    out.T5 = out.T_vec(4);
    out.T6 = T6;
    out.T7 = T6-8;
    T2 = out.T_vec(1);
    T3 = out.T_vec(2);
    T4 = out.T_vec(3);
    T5 = out.T_vec(4);
    T7 = T6-8;
    
    eta_opt5 = eta_opt;
    eta_opt3 = eta_opt*tau5;
    
    %Energy balance at point 2 (tube inner wall/HTF interface)
    out_21_conv = convection21(P1,T1,T2,D2,L,m_dot1, fluid1);
    out.Qdot_21_conv = out_21_conv.Qdot;
    out.flag_21_conv = out_21_conv.flag;
    
    out_32_cond = conduction32(T2,T3,D2,D3,L,k23);
    out.Qdot_32_cond = out_32_cond.Qdot;
    out.flag_32_cond = out_32_cond.flag;
    
    out.res1 = out_21_conv.Qdot - out_32_cond.Qdot;
    
    %Energy balance at point 3 (tube outer wall/inner gas interface)
    out_34_conv = convection34(T3,T4,D3,D4,L,P34,fluid34);
    out.Qdot_34_conv = out_34_conv.Qdot;
    out.flag_34_conv = out_34_conv.flag;
    
    out_34_rad = radiation34(T3,T4,D3,D4,L,eps3, eps4);
    out.Qdot_34_rad = out_34_rad.Qdot;
    out.flag_34_rad = out_34_rad.flag;
    
    out_3_solAbs = TubeSolAbsorption(L, W_ptc, DNI, cos_theta, alpha3, eta_opt3);
    out.Qdot_3_solAbs = out_3_solAbs.Qdot;
    out.Qdot_3_solRefl = out_3_solAbs.Qdot_refl;
    out.flag_3_solAbs = out_3_solAbs.flag;    
    
    out.res3 = out_3_solAbs.Qdot - out_32_cond.Qdot - out_34_conv.Qdot - out_34_rad.Qdot;
    
    %Energy balance at point 4 (envelope inner wall/inner gas interface)
    out_45_cond = conduction45(T4,T5,D4,D5,L,k45);
    out.Qdot_45_cond = out_45_cond.Qdot;
    out.flag_45_cond = out_45_cond.flag;   
    
    out.res4 = out_34_conv.Qdot + out_34_rad.Qdot - out_45_cond.Qdot;
    
    %Energy balance at point 5 (envelope outer wall/ambiance interface)
    out_57_rad = radiation57(T5,T7,D5,L,eps5);
    out.Qdot_57_rad = out_57_rad.Qdot;
    out.flag_57_rad = out_57_rad.flag;
    
    out_56_conv = convection56(P56,T5,T6,D5,L,v_wind);
    out.Qdot_56_conv = out_56_conv.Qdot;
    out.flag_56_conv = out_56_conv.flag;
    
    out_5_solAbs = EnvelopeSolAbsorption(L, W_ptc, DNI, cos_theta, alpha5, eta_opt5, tau5);
    out.Qdot_5_solAbs = out_5_solAbs.Qdot;
    out.Qdot_5_solRefl = out_5_solAbs.Qdot_refl;
    out.flag_5_solAbs = out_5_solAbs.flag;
    out.res5 = out_45_cond.Qdot + out_5_solAbs.Qdot - out_56_conv.Qdot - out_57_rad.Qdot;
    
    % Heat losses
    out.Qdot_losses = out.Qdot_56_conv + out.Qdot_57_rad;
    out.eta_ptc = out.Qdot_21_conv/(DNI*cos_theta*L*W_ptc);    
    
    % Not evaluated
    out.Qdot_37_rad = NaN;
    out.flag_37_rad = NaN;
    out.Qdot_36_conv = NaN;
    out.flag_36_conv = NaN;
    
    %residuals
    out.res = [out.res1 out.res3 out.res4 out.res5];
    
    % Flag cheking
    if min([out.flag_21_conv, out.flag_32_cond, out.flag_34_conv, out.flag_34_rad, out.flag_3_solAbs, out.flag_45_cond, out.flag_57_rad, out.flag_56_conv, out.flag_5_solAbs, out.flag_36_conv, out.flag_37_rad]) <0 || norm(out.res) > 1e-4
        display('Check for correlations validity - solution may be wrong')
        out.flag = -1;
    else
        out.flag = 1;
    end
    
elseif not(GlassEnvelop)
    
    out.T1 = T1;
    out.T2 = out.T_vec(1);
    out.T3 = out.T_vec(2);
    out.T4 = NaN;
    out.T5 = NaN;
    out.T6 = T6;
    out.T7 = T6-8;
    T2 = out.T_vec(1);
    T3 = out.T_vec(2);
    T7 = T6-8;
    
    eta_opt3 = eta_opt;
    
    %Energy balance at point 2 (tube inner wall/HTF interface)
    out_21_conv = convection21(P1,T1,T2,D2,L,m_dot1, fluid1);
    out.Qdot_21_conv = out_21_conv.Qdot;
    out.flag_21_conv = out_21_conv.flag;
    
    out_32_cond = conduction32(T2,T3,D2,D3,L,k23);
    out.Qdot_32_cond = out_32_cond.Qdot;
    out.flag_32_cond = out_32_cond.flag;
    
    out.res1 = out_21_conv.Qdot - out_32_cond.Qdot;
    
    %Energy balance at point 3 (tube outer wall/atmosphere interface)
    out_36_conv = convection56(P56,T3,T6,D3,L,v_wind);
    out.Qdot_36_conv = out_36_conv.Qdot;
    out.flag_36_conv = out_36_conv.flag;
    
    out_37_rad = radiation57(T3,T7,D3,L,eps3);
    out.Qdot_37_rad = out_37_rad.Qdot;
    out.flag_37_rad = out_37_rad.flag;
    
    out_3_solAbs = TubeSolAbsorption(L, W_ptc, DNI, cos_theta, alpha3, eta_opt3);
    out.Qdot_3_solAbs = out_3_solAbs.Qdot;
    out.Qdot_3_solRefl = out_3_solAbs.Qdot_refl;
    out.flag_3_solAbs = out_3_solAbs.flag;
    
    out.res3 = out_3_solAbs.Qdot - out_32_cond.Qdot - out_36_conv.Qdot - out_37_rad.Qdot;
    

    % Heat losses
    out.Qdot_losses = out.Qdot_36_conv + out.Qdot_37_rad;
    out.eta_ptc = out.Qdot_21_conv/(DNI*cos_theta*L*W_ptc);
    
    % Not evaluated
    out.Qdot_34_conv = NaN;
    out.flag_34_conv = NaN;
    out.Qdot_34_rad = NaN;
    out.flag_34_rad = NaN;
    out.Qdot_45_cond = NaN;
    out.flag_45_cond = NaN;
    out.Qdot_57_rad = NaN;
    out.flag_57_rad = NaN;
    out.Qdot_56_conv = NaN;
    out.flag_56_conv = NaN;
    out.Qdot_5_solAbs = NaN;
    out.flag_5_solAbs = NaN;
    out.Qdot_5_solRefl = NaN;
    out.res4 = NaN ;
    out.res5 = NaN ;
    
    out.res = [out.res1 out.res3];
    
    % Flag cheking
    if min([out.flag_21_conv, out.flag_32_cond, out.flag_34_conv, out.flag_34_rad, out.flag_3_solAbs, out.flag_45_cond, out.flag_57_rad, out.flag_56_conv, out.flag_5_solAbs, out.flag_36_conv, out.flag_37_rad]) <0 || norm(out.res) > 1e-4
        display('Check for correlations validity - solution may be wrong')
        out.flag = -1;
    else
        out.flag = 1;
    end
    
    
end

out = orderfields(out);

end

%Resolution of energy balance
function T_vec = solve_energy_balance(fluid1, T1, P1, m_dot1, T6, D2, D3, D4, D5, L, W_ptc, fluid34, P34, P56, v_wind, eps3, eps4, eps5, alpha3, alpha5, tau5, eta_opt, k45, k23, DNI, cos_theta, GlassEnvelop)

T7 = T6-8;% sky temperature

if GlassEnvelop
    DT = 20;
    x0 = linspace(T1+DT, T6, 4);
    ub = (T1+200)*ones(1,4);
    f2opt = @(x) energy_balance(x, ub, T1, T6, T7, D2, D3, D4, D5, L, W_ptc, fluid1, P1, m_dot1, fluid34, P34, P56, v_wind, eps3, eps4, eps5, alpha3, alpha5, tau5, eta_opt, k45, k23, DNI, cos_theta, GlassEnvelop);
    options = optimset('Display','none');
    x_opt = fsolve(f2opt, x0./ub, options);%
    T_vec = x_opt.*ub;
    
elseif not(GlassEnvelop)
    DT = 20;
    x0 = linspace(T1+DT, T6, 2);
    ub = (T1+200)*ones(1,2);
    f2opt = @(x) energy_balance(x, ub, T1, T6, T7, D2, D3, D4, D5, L, W_ptc, fluid1, P1, m_dot1, fluid34, P34, P56, v_wind, eps3, eps4, eps5, alpha3, alpha5, tau5, eta_opt, k45, k23, DNI, cos_theta, GlassEnvelop);
    options = optimset('Display','none');
    x_opt = fsolve(f2opt, x0./ub, options);%
    T_vec = x_opt.*ub;
end

end

%Energy balance
function res = energy_balance(x, ub, T1, T6, T7, D2, D3, D4, D5, L, W_ptc, fluid1, P1, m_dot1, fluid34, P34, P56, v_wind, eps3, eps4, eps5, alpha3, alpha5, tau5, eta_opt, k45, k23, DNI, cos_theta, GlassEnvelop)
if GlassEnvelop
    x = x.*ub;
    T2 = x(1);
    T3 = x(2);
    T4 = x(3);
    T5 = x(4);
    
    eta_opt5 = eta_opt;
    eta_opt3 = eta_opt*tau5;
    
    %Energy balance at point 2 (tube inner wall/HTF interface)
    out_21_conv = convection21(P1,T1,T2,D2,L,m_dot1, fluid1);
    out_32_cond = conduction32(T2,T3,D2,D3,L,k23);    
    res1 = out_21_conv.Qdot - out_32_cond.Qdot;
    
    %Energy balance at point 3 (tube outer wall/inner gas interface)
    out_34_conv = convection34(T3,T4,D3,D4,L,P34,fluid34);
    out_34_rad = radiation34(T3,T4,D3,D4,L,eps3, eps4);
    out_3_solAbs = TubeSolAbsorption(L, W_ptc, DNI, cos_theta, alpha3, eta_opt3);        
    res3 = out_3_solAbs.Qdot - out_32_cond.Qdot - out_34_conv.Qdot - out_34_rad.Qdot;
    
    %Energy balance at point 4 (envelope inner wall/inner gas interface)
    out_45_cond = conduction45(T4,T5,D4,D5,L,k45);
    res4 = out_34_conv.Qdot + out_34_rad.Qdot - out_45_cond.Qdot;
    
    %Energy balance at point 5 (envelope outer wall/ambiance interface)
    out_57_rad = radiation57(T5,T7,D5,L,eps5);
    out_56_conv = convection56(P56,T5,T6,D5,L,v_wind);
    out_5_solAbs = EnvelopeSolAbsorption(L, W_ptc, DNI, cos_theta, alpha5, eta_opt5, tau5);    
    res5 = out_45_cond.Qdot + out_5_solAbs.Qdot - out_56_conv.Qdot - out_57_rad.Qdot;
    
    %Residual verctor
    res = [res1 res3 res4 res5];
    
elseif not(GlassEnvelop)
    x = x.*ub;
    T2 = x(1);
    T3 = x(2);    
    
    eta_opt3 = eta_opt;
    
    %Energy balance at point 2 (tube inner wall/HTF interface)
    out_21_conv = convection21(P1,T1,T2,D2,L,m_dot1, fluid1);
    out_32_cond = conduction32(T2,T3,D2,D3,L,k23);
    res1 = out_21_conv.Qdot - out_32_cond.Qdot;
    
    %Energy balance at point 3 (tube outer wall/atmosphere interface)
    out_36_conv = convection56(P56,T3,T6,D3,L,v_wind);
    out_37_rad = radiation57(T3,T7,D3,L,eps3);
    out_3_solAbs = TubeSolAbsorption(L, W_ptc, DNI, cos_theta, alpha3, eta_opt3);    
    res3 = out_3_solAbs.Qdot - out_32_cond.Qdot - out_36_conv.Qdot - out_37_rad.Qdot;
        
    %Residual verctor
    res = [res1 res3 ];
    
end
%res
end

%Convection from absorber inner wall (2) to HTF (1) --> ok
function out_21_conv = convection21(P1,T1,T2,D2,L,m_dot1, fluid1)
flag = 1;
cp1 = sf_PropsSI_bar('C', T1, T1, P1, fluid1);
k1 = sf_PropsSI_bar('L', T1, T1, P1, fluid1);
mu1 = sf_PropsSI_bar('V', T1, T1, P1, fluid1);
Pr1 = cp1*mu1/k1;
cp2 = sf_PropsSI_bar('C', T2, T2, P1, fluid1);
k2 = sf_PropsSI_bar('L', T2, T2, P1, fluid1);
mu2 = sf_PropsSI_bar('V', T2, T2, P1, fluid1);
Pr2 = cp2*mu2/k2;
G1 = m_dot1/(pi*D2^2/4);
Re1 = G1*D2/mu1;
if Re1 > 2300 %turbulent flow 
    f2 = (1.82*log10(Re1)-1.64)^-2;
    Nu1 = ((f2/8)*(Re1-1000)*Pr1)/(1+12.7*sqrt(f2/8)*(Pr1^(2/3)-1))*(Pr1/Pr2)^0.11; % Gnielinski
    if Re1>5e6 || Pr1 < 0.5 || Pr1 > 2000 || Pr2 < 0.5 || Pr2 > 2000
        flag = -1;
    end
else %laminar flow
    Nu_1 = 4.36;
    L_tot = 30;
    Nu_2 = 1.953*(Re1*Pr1*D2/L_tot)^0.33333333333333333333333333333;
    Nu1 = (Nu_1^3 + 0.6^3 + (Nu_2-0.6)^3)^0.3333333333333333333333;
end
h1 = Nu1*k1/D2;
Qdot_21_conv = pi*D2*L*h1*(T2-T1);
out_21_conv.Qdot = Qdot_21_conv;
out_21_conv.flag = flag;
end

%Conduction from absorber outer wall (3) to absorber inner wall (3) --> ok
function out_32_cond = conduction32(T2,T3,D2,D3,L,k23)
Qdot_32_cond = 2*pi*k23*L*(T3-T2)/log(D3/D2);
out_32_cond.Qdot = Qdot_32_cond;
out_32_cond.flag = 1;
end

%Convection from absorber outer wall (3) and glass envelop inner wall (4) -> ok
function out_34_conv = convection34(T3,T4,D3,D4,L,P34,fluid34)
if P34>0
    flag = 1;
    mu34 = CoolProp.PropsSI('V', 'T', (T3+T4)/2, 'P', P34, fluid34);
    Pr34 = CoolProp.PropsSI('Prandtl', 'T', (T3+T4)/2, 'P', P34, fluid34);
    k34 = CoolProp.PropsSI('L', 'T', (T3+T4)/2, 'P', P34, fluid34);
    cp34 = CoolProp.PropsSI('C', 'T', (T3+T4)/2, 'P', P34, fluid34);
    rho34 = CoolProp.PropsSI('D', 'T', (T3+T4)/2, 'P', P34, fluid34);
    alfa34 = k34/(cp34*rho34); %thermal diffusivité [m²/s]
    nu34 = mu34/rho34; %cinematic viscosity [m²/s]
    Ra34 = 9.81*(1/(0.5*T3+0.5*T4))*abs(T3-T4)*D3^3/alfa34/nu34;
    if Ra34<(D4/(D4-D3))^4
        flag = -1;
    end
    Qdot_34_conv = (2.425*k34*L*(T3-T4)*(Pr34*Ra34/(0.861+Pr34))^0.25)/((1+(D3/D4)^0.6)^1.25);
    out_34_conv.Qdot = Qdot_34_conv;
    out_34_conv.flag = flag;
else
    disp('To implement vaccuum case')
end
end

% Radiation from absorber outer wall (3) and glass envelop inner wall (4) -> ok
function out_34_rad = radiation34(T3,T4,D3,D4,L,eps3, eps4)
sigma34 = 5.670367e-8;
Qdot_34_rad = sigma34*pi*D3*L*(T3^4-T4^4)/(1/eps3+(1-eps4)*D3/(eps4*D4));
out_34_rad.Qdot = Qdot_34_rad;
out_34_rad.flag = 1;
end

% Conduction from glass envelop inner wall (4) to glass envelop outer wall (5)-> ok
function out_45_cond = conduction45(T4,T5,D4,D5,L,k45)
Qdot_45_cond = 2*pi*k45*L*(T4-T5)/log(D5/D4);
out_45_cond.Qdot = Qdot_45_cond;
out_45_cond.flag = 1;
end

% Convection from glass envelop outer wall (5) to the atmosphere (6) ->ok
function out_56_conv = convection56(P56,T5,T6,D5,L,v_wind)
flag = 1;
if v_wind < 0.01 % if no wind, natural convection with Churchill & Chu's correlation
    mu56 = CoolProp.PropsSI('V', 'T', (T5+T6)/2, 'P', P56, 'air');
    rho56 = CoolProp.PropsSI('D', 'T', (T5+T6)/2, 'P', P56, 'air');
    Pr56 = CoolProp.PropsSI('Prandtl', 'T', (T5+T6)/2, 'P', P56, 'air');
    k56 = CoolProp.PropsSI('L', 'T', (T5+T6)/2, 'P', P56, 'air');
    cp56 = CoolProp.PropsSI('C', 'T', (T5+T6)/2, 'P', P56, 'air');
    alfa56 = k56/(cp56*rho56); %thermal diffusivité [m²/s]
    nu56 = mu56/rho56; %cinematic viscosity [m²/s]
    Ra56 = 9.81*(1/(0.5*T5 + 0.5*T6))*abs(T5-T6)*D5^3/alfa56/nu56;
    Nu56 = (0.6+(0.387*Ra56^(1/6))/((1+(0.559/Pr56)^(9/16))^(8/27)))^2;
    if Ra56<1e5 || Ra56 >1e12
        flag = -1;
    end
else %if wind, forced convection with Zhukauska's correlation
    Pr6 = CoolProp.PropsSI('Prandtl', 'T', T6, 'P', P56, 'air');
    Pr5 = CoolProp.PropsSI('Prandtl', 'T', T5, 'P', P56, 'air');
    mu6 = CoolProp.PropsSI('V', 'T', T6, 'P', P56, 'air');
    rho6 = CoolProp.PropsSI('V', 'T', T6, 'P', P56, 'air');
    Re6 = rho6*v_wind*D5/mu6;
    k56 = CoolProp.PropsSI('L', 'T', (T5+T6)/2, 'P', P56, 'air');
    if  Re6 <= 40
        C = 0.75;
        m = 0.4;
    elseif Re6 > 40 && Re6 <= 1000
        C = 0.51;
        m = 0.5;
    elseif Re6 > 1000 && Re6 <= 200000
        C = 0.26;
        m = 0.6;
    elseif Re6 > 200000 
        C = 0.076;
        m = 0.7;
    end
        
    if Pr6 <= 10
        n = 0.37;
    else
        n = 0.36;
    end

    if Pr6<0.7 || Pr6>500 || Re6>1e6 %|| Re6<1         
        flag = -1;
    end
    Nu56 = C*(Re6^m)*(Pr6^n)*(Pr6/Pr5)^0.25;
end

h56 = Nu56*k56/D5;
Qdot_56_conv = pi*D5*L*h56*(T5-T6);
out_56_conv.Qdot = Qdot_56_conv;
out_56_conv.flag = flag;

end

% Radiation from glass envelop outer wall (5) to the atmosphere (7) ->ok
function out_57_rad = radiation57(T5,T7,D5,L,eps5)
sigma57 = 5.670367e-8; %Boltzmann constanst [W/m^2.K^4]
Qdot_57_rad = sigma57*pi*D5*L*eps5*(T5^4-T7^4);
out_57_rad.Qdot = Qdot_57_rad;
out_57_rad.flag = 1;
end

% Absorption in the absorber tube
function out_3_solAbs = TubeSolAbsorption(L, W_ptc, DNI, cos_theta, alpha3, eta_opt3)
out_3_solAbs.Qdot = L*W_ptc*DNI*cos_theta*eta_opt3*alpha3;
out_3_solAbs.Qdot_refl = L*W_ptc*DNI*cos_theta*eta_opt3*(1-alpha3);

out_3_solAbs.flag = 1;
end

% Absorption in the glass envelope
function out_5_solAbs = EnvelopeSolAbsorption(L, W_ptc, DNI, cos_theta, alpha5, eta_opt5, tau5)
out_5_solAbs.Qdot = L*W_ptc*DNI*cos_theta*eta_opt5*alpha5;
out_5_solAbs.Qdot_refl = L*W_ptc*DNI*cos_theta*eta_opt5*(1-tau5-alpha5);
out_5_solAbs.flag = 1;
end
