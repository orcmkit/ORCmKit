clear all
close all
clc

fluid = 'R245fa';


%Saturated properties based on pressure
P_max = CoolProp.PropsSI('Pcrit', 'P', 1e5, 'Q',0, fluid)*0.99;
P_min = CoolProp.PropsSI('P', 'T', -40+273.15, 'Q',0, fluid);
N_P = 150;
P_vec = linspace(P_min,P_max,N_P);
m = 0;
for P = P_vec
    m = m +1;
    hLsat_p_vec(m)      = CoolProp.PropsSI('H', 'P', P_vec(m), 'Q',0, fluid);
    hVsat_p_vec(m)      = CoolProp.PropsSI('H', 'P', P_vec(m), 'Q',1, fluid);
    Tsat_p_vec(m)       = CoolProp.PropsSI('T', 'P', P_vec(m), 'Q',0.5, fluid);
    rhoLsat_p_vec(m)    = CoolProp.PropsSI('D', 'P', P_vec(m), 'Q',0, fluid);
    rhoVsat_p_vec(m)    = CoolProp.PropsSI('D', 'P', P_vec(m), 'Q',1, fluid);
    muLsat_p_vec(m)     = CoolProp.PropsSI('V', 'P', P_vec(m), 'Q',0, fluid);
    muVsat_p_vec(m)     = CoolProp.PropsSI('V', 'P', P_vec(m), 'Q',1, fluid);
    sigLsat_p_vec(m)     = CoolProp.PropsSI('I', 'P', P_vec(m), 'Q',0, fluid);
    kLsat_p_vec(m)      = CoolProp.PropsSI('L', 'P', P_vec(m), 'Q',0, fluid);
    PrLsat_p_vec(m)      = CoolProp.PropsSI('Prandtl', 'P', P_vec(m), 'Q',0, fluid);
    iLat_p_vec(m) = hVsat_p_vec(m)-hLsat_p_vec(m);
end
hLsat_p_R245fa  = griddedInterpolant(P_vec, hLsat_p_vec, 'cubic', 'linear');
hVsat_p_R245fa  = griddedInterpolant(P_vec, hVsat_p_vec, 'cubic', 'linear');
Tsat_p_R245fa   = griddedInterpolant(P_vec, Tsat_p_vec , 'cubic', 'linear');
rhoLsat_p_R245fa  = griddedInterpolant(P_vec, rhoLsat_p_vec, 'cubic', 'linear');
rhoVsat_p_R245fa  = griddedInterpolant(P_vec, rhoVsat_p_vec, 'cubic', 'linear');
muLsat_p_R245fa  = griddedInterpolant(P_vec, muLsat_p_vec, 'cubic', 'linear');
muVsat_p_R245fa  = griddedInterpolant(P_vec, muVsat_p_vec, 'cubic', 'linear');
sigLsat_p_R245fa  = griddedInterpolant(P_vec, sigLsat_p_vec, 'cubic', 'linear');
kLsat_p_R245fa  = griddedInterpolant(P_vec, kLsat_p_vec, 'cubic', 'linear');
PrLsat_p_R245fa  = griddedInterpolant(P_vec, PrLsat_p_vec, 'cubic', 'linear');
iLat_p_R245fa = griddedInterpolant(P_vec, iLat_p_vec, 'cubic', 'linear');
clearvars P_vec

%Saturated properties based on temperature
T_min = -40 +273.15;
T_max = 152 +273.15;
N_T = 200;
T_vec = linspace(T_min,T_max,N_T);
m = 0;
for T = T_vec
    m = m +1;
    hLsat_t_vec(m) = CoolProp.PropsSI('H', 'T', T_vec(m), 'Q',0, fluid);
    hVsat_t_vec(m) = CoolProp.PropsSI('H', 'T', T_vec(m), 'Q',1, fluid);
    Psat_t_vec(m) =  CoolProp.PropsSI('P', 'T', T_vec(m), 'Q',0.5, fluid); 
end
hLsat_t_R245fa  = griddedInterpolant(T_vec, hLsat_t_vec, 'cubic', 'linear');
hVsat_t_R245fa  = griddedInterpolant(T_vec, hVsat_t_vec, 'cubic', 'linear');
Psat_t_R245fa   = griddedInterpolant(T_vec, Psat_t_vec , 'cubic', 'linear');
clearvars T_vec

%enthalpy based on temperature/pressure
clearvars P_vec m N_p
N_P = 500;
N_T = 500;
T_min = -50 +273.15;
T_max = 152 +273.15;
P_vec = linspace(P_min,P_max,N_P);
T_vec = linspace(T_min,T_max,N_T);
[P_MX,T_MX] = ndgrid(P_vec,T_vec);
for m = 1:length(P_vec)
    for n = 1:length(T_vec)
        try
            H_MX(m,n) = CoolProp.PropsSI('H', 'P', P_MX(m,n), 'T',T_MX(m,n), fluid);
        catch
            H_MX(m,n) = NaN;
        end
    end
end
h_pt_R245fa  = griddedInterpolant(P_MX,T_MX, H_MX, 'linear', 'linear');
clearvars P_MX H_MX T_MX

%Temperature based on on enthalpy/pressure
N_P = 500;
N_h = 500;
P_vec = linspace(P_min,P_max,N_P);
T_min = -50 +273.15;
T_max = 160 +273.15;
h_min = min([CoolProp.PropsSI('H', 'P', P_min, 'T',T_min, fluid)    CoolProp.PropsSI('H', 'P', P_min, 'T',T_max, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_min, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_max, fluid)]);%CoolProp.PropsSI('H', 'P', P_min, 'Q',0, fluid);
h_max = max([CoolProp.PropsSI('H', 'P', P_min, 'T',T_min, fluid)    CoolProp.PropsSI('H', 'P', P_min, 'T',T_max, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_min, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_max, fluid)]);
h_vec = linspace(h_min,h_max,N_h);
[P_MX,H_MX] = ndgrid(P_vec,h_vec);
for m = 1:length(P_vec)
    for n = 1:length(h_vec)
        T_MX(m,n) = CoolProp.PropsSI('T', 'P', P_MX(m,n), 'H',H_MX(m,n), fluid);
        s_MX(m,n) = CoolProp.PropsSI('S', 'P', P_MX(m,n), 'H',H_MX(m,n), fluid);
        rho_MX(m,n) = CoolProp.PropsSI('D', 'P', P_MX(m,n), 'H',H_MX(m,n), fluid);
    end
end
T_ph_R245fa  = griddedInterpolant(P_MX,H_MX, T_MX, 'cubic', 'linear');
s_ph_R245fa  = griddedInterpolant(P_MX,H_MX, s_MX, 'cubic', 'linear');
rho_ph_R245fa  = griddedInterpolant(P_MX,H_MX, rho_MX, 'cubic', 'linear');
clearvars P_MX H_MX T_MX

save('PropTab_R245fa.mat', 'rho_ph_R245fa', 's_ph_R245fa', 'T_ph_R245fa', 'h_pt_R245fa', 'Psat_t_R245fa', 'hVsat_t_R245fa', 'hLsat_t_R245fa', 'hLsat_p_R245fa', 'hVsat_p_R245fa', 'Tsat_p_R245fa', 'rhoLsat_p_R245fa', 'rhoVsat_p_R245fa', 'muLsat_p_R245fa', 'muVsat_p_R245fa', 'sigLsat_p_R245fa', 'kLsat_p_R245fa', 'PrLsat_p_R245fa', 'iLat_p_R245fa')

if 1
    %Saturation properties check
    P_vec2 = linspace(P_min,P_max,600);
    
    k1 = 0;
    [iLat_p_verif_CP, PrLsat_p_verif_CP, kLsat_p_verif_CP, Psat_t_verif_CP, hVsat_t_verif_CP, hLsat_t_verif_CP, sigLsat_p_verif_CP, Tsat_p_verif_CP, muVsat_p_verif_CP, muLsat_p_verif_CP, rhoVsat_p_verif_CP, rhoLsat_p_verif_CP, hLsat_p_verif_CP, hVsat_p_verif_CP, ] = deal(NaN*ones(length(P_vec2),1));
    tic
    for m = 1:length(P_vec2)
        k1 = k1 + 1;
        hLsat_p_verif_CP(k1)   = CoolProp.PropsSI('H', 'P', P_vec2(m), 'Q',0, fluid);
        hVsat_p_verif_CP(k1)   = CoolProp.PropsSI('H', 'P', P_vec2(m), 'Q',1, fluid);
        rhoLsat_p_verif_CP(k1) = CoolProp.PropsSI('D', 'P', P_vec2(m), 'Q',0, fluid);
        rhoVsat_p_verif_CP(k1) = CoolProp.PropsSI('D', 'P', P_vec2(m), 'Q',1, fluid);
        muLsat_p_verif_CP(k1) = CoolProp.PropsSI('V', 'P', P_vec2(m), 'Q',0, fluid);
        muVsat_p_verif_CP(k1) = CoolProp.PropsSI('V', 'P', P_vec2(m), 'Q',1, fluid);
        Tsat_p_verif_CP(k1) = CoolProp.PropsSI('T', 'P', P_vec2(m), 'Q',1, fluid);
        sigLsat_p_verif_CP(k1) = CoolProp.PropsSI('I', 'P', P_vec2(m), 'Q',0, fluid);
        kLsat_p_verif_CP(k1)   = CoolProp.PropsSI('L', 'P', P_vec2(m), 'Q',0, fluid);
        PrLsat_p_verif_CP(k1)  = CoolProp.PropsSI('Prandtl', 'P', P_vec2(m), 'Q',0, fluid);
        iLat_p_verif_CP(k1) = CoolProp.PropsSI('H', 'P', P_vec2(m), 'Q',1, fluid)-CoolProp.PropsSI('H', 'P', P_vec2(m), 'Q',0, fluid);      
        hLsat_t_verif_CP(k1) = CoolProp.PropsSI('H', 'T', Tsat_p_verif_CP(k1), 'Q',0, fluid);
        hVsat_t_verif_CP(k1) = CoolProp.PropsSI('H', 'T', Tsat_p_verif_CP(k1), 'Q',1, fluid);
        Psat_t_verif_CP(k1) = CoolProp.PropsSI('P', 'T', Tsat_p_verif_CP(k1), 'Q',0, fluid);
    end
    time_CP = toc
    

    T_vec2 = Tsat_p_verif_CP;
    k1 = 0;
    [iLat_p_verif_IT, PrLsat_p_verif_IT, kLsat_p_verif_IT, Psat_t_verif_IT, hVsat_t_verif_IT, hLsat_t_verif_IT, sigLsat_p_verif_IT, Tsat_p_verif_IT, muVsat_p_verif_IT, muLsat_p_verif_IT, rhoVsat_p_verif_IT, rhoLsat_p_verif_IT, hLsat_p_verif_IT, hVsat_p_verif_IT, ] = deal(NaN*ones(length(P_vec2),1));
    tic
    for m = 1:length(P_vec2)
        k1 = k1 + 1;
        hLsat_p_verif_IT(k1) = hLsat_p_R245fa(P_vec2(m));
        hVsat_p_verif_IT(k1) = hVsat_p_R245fa(P_vec2(m));
        rhoLsat_p_verif_IT(k1) = rhoLsat_p_R245fa(P_vec2(m));
        rhoVsat_p_verif_IT(k1) = rhoVsat_p_R245fa(P_vec2(m));
        muLsat_p_verif_IT(k1) = muLsat_p_R245fa(P_vec2(m));
        muVsat_p_verif_IT(k1) = muVsat_p_R245fa(P_vec2(m));
        Tsat_p_verif_IT(k1) = Tsat_p_R245fa(P_vec2(m));
        sigLsat_p_verif_IT(k1) = sigLsat_p_R245fa(P_vec2(m));
        kLsat_p_verif_IT(k1) = kLsat_p_R245fa(P_vec2(m));
        PrLsat_p_verif_IT(k1) = PrLsat_p_R245fa(P_vec2(m));
        iLat_p_verif_IT(k1)  = iLat_p_R245fa(P_vec2(m));
        hLsat_t_verif_IT(k1) = hLsat_t_R245fa(T_vec2(m));
        hVsat_t_verif_IT(k1) = hVsat_t_R245fa(T_vec2(m));
        Psat_t_verif_IT(k1) = Psat_t_R245fa(T_vec2(m));
    end
    time_IT = toc
    
    Rate = time_CP/time_IT
    
    figure;
    subplot(2,2,1)
    hold on
    plot(hLsat_p_verif_CP,hLsat_p_verif_IT,'o')
    plot(hVsat_p_verif_CP,hVsat_p_verif_IT,'s')
    plot(hLsat_t_verif_CP,hLsat_t_verif_IT,'v')
    plot(hVsat_t_verif_CP,hVsat_t_verif_IT,'^')
    plot(iLat_p_verif_CP,iLat_p_verif_IT,'^')    
    grid on
    subplot(2,2,2)
    hold on
    plot(Psat_t_verif_CP/2000,Psat_t_verif_IT/2000,'^')
    plot(rhoLsat_p_verif_CP,rhoLsat_p_verif_IT,'o')
    plot(rhoVsat_p_verif_CP,rhoVsat_p_verif_IT,'s')
    plot(Tsat_p_verif_CP,Tsat_p_verif_IT,'v')
    plot(sigLsat_p_verif_CP*100000,sigLsat_p_verif_IT*100000,'*')
    grid on
    subplot(2,2,3)
    hold on
    plot(kLsat_p_verif_CP*500,kLsat_p_verif_IT*500,'^')
    plot(PrLsat_p_verif_CP,PrLsat_p_verif_IT,'o')
    grid on    
    
    % Other properties check
    T_min = -50 +273.15;
    T_max = 152 +273.15;
    h_min = min([CoolProp.PropsSI('H', 'P', P_min, 'T',T_min, fluid)    CoolProp.PropsSI('H', 'P', P_min, 'T',T_max, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_min, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_max, fluid)]);%CoolProp.PropsSI('H', 'P', P_min, 'Q',0, fluid);
    h_max = max([CoolProp.PropsSI('H', 'P', P_min, 'T',T_min, fluid)    CoolProp.PropsSI('H', 'P', P_min, 'T',T_max, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_min, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_max, fluid)]);
    P_vec3 = linspace(P_min,P_max,182);
    h_vec3 = linspace(h_min+100,h_max-100,223);
    T_vec3 = linspace(T_min, T_max, length(h_vec3));
    [s_IT, rho_IT, rho_CP, s_CP, T_IT, T_CP, H_IT, H_CP,Tsat_CP,Tsat_IT] = deal(NaN*ones(length(P_vec3)*length(h_vec3),1));
    k2 = 0;
    tic
    for m = 1:length(P_vec3)
        for n = 1:length(h_vec3)
            k2 = k2 + 1;
            T_CP(k2) = CoolProp.PropsSI('T', 'P', P_vec3(m), 'H',h_vec3(n), fluid);
            s_CP(k2) = CoolProp.PropsSI('S', 'P', P_vec3(m), 'H',h_vec3(n), fluid);
            h_liqsat = CoolProp.PropsSI('H', 'P', P_vec3(m), 'Q',0, fluid);
            h_vapsat = CoolProp.PropsSI('H', 'P', P_vec3(m), 'Q',1, fluid);
            if h_vec3(n) < 0.999*h_liqsat || h_vec3(n) > 1.001*h_vapsat
                rho_CP(k2) = CoolProp.PropsSI('D', 'P', P_vec3(m), 'H',h_vec3(n), fluid);
            else
                rho_CP(k2) = NaN;
            end
        end
    end
    k2 = 0;
    for m = 1:length(P_vec3)
        for n = 1:length(T_vec3)
            k2 = k2 + 1;
            Tsat_CP(k2) = CoolProp.PropsSI('T', 'P', P_vec3(m), 'Q',0.5, fluid);
            if abs(T_vec3(n)-Tsat_CP(k2))>1
                H_CP(k2) = CoolProp.PropsSI('H', 'P', P_vec3(m), 'T', T_vec3(n), fluid);
            else
                H_CP(k2) = 0;
            end
        end
    end
    time_CP2 = toc
    
    k2 = 0;
    tic
    for m = 1:length(P_vec3)
        for n = 1:length(h_vec3)
            k2 = k2 + 1;
            T_IT(k2) = T_ph_R245fa(P_vec3(m),h_vec3(n));
            s_IT(k2) = s_ph_R245fa(P_vec3(m),h_vec3(n));
            h_liqsat = hLsat_p_R245fa(P_vec3(m));
            h_vapsat = hVsat_p_R245fa(P_vec3(m));
            if h_vec3(n) < 0.999*h_liqsat || h_vec3(n) > 1.001*h_vapsat
                rho_IT(k2) = rho_ph_R245fa(P_vec3(m),h_vec3(n));
            else
                rho_IT(k2) =  NaN;
            end
        end
    end
    k2 = 0;
    for m = 1:length(P_vec3)
        for n = 1:length(T_vec3)
            k2 = k2 + 1;
            Tsat_IT(k2) = Tsat_p_R245fa(P_vec3(m));
            if abs(T_vec3(n)-Tsat_IT(k2))>1
                H_IT(k2) = h_pt_R245fa(P_vec3(m),T_vec3(n));
            else
                H_IT(k2) =  0;
            end
        end
    end
    time_IT2 = toc
    
    Rate2 = time_CP2/time_IT2
    figure
    subplot(2,2,1)
    hold on
    plot(T_IT,T_CP, 'o')
    plot(Tsat_IT,Tsat_CP, 'o')
    grid on
    subplot(2,2,2)
    hold on
    plot(H_IT,H_CP, 'o')    
    grid on
    subplot(2,2,3)
    hold on
    plot(s_IT,s_CP, 'o')    
    subplot(2,2,4)
    hold on
    plot(rho_IT,rho_CP, 'o')    
    grid on
    
end

%save('PropsTab_R245fa.mat','hLsat_p_R245fa','hVsat_p_R245fa','Tsat_p_R245fa','rhoLsat_p_R245fa', 'rhoVsat_p_R245fa','T_ph_R245fa')

% 
% 
% clear all
% close all
% clc
% 
% fluid = 'air';
% 
% P_max = 1.5e5;%CoolProp.PropsSI('Pcrit', 'P', 1e5, 'Q',0, fluid)*0.99;
% P_min = 0.5e5;%CoolProp.PropsSI('P', 'T', -40+273.15, 'Q',0, fluid);
% 
% %Saturated properties based on pressure
% N_P = 3;
% P_vec = linspace(P_min,P_max,N_P);
% m = 0;
% for P = P_vec
%     m = m +1;
%     h_lSat(m) = CoolProp.PropsSI('H', 'P', P_vec(m), 'Q',0, fluid);
%     h_vSat(m) = CoolProp.PropsSI('H', 'P', P_vec(m), 'Q',1, fluid);
%     T_Sat(m) = CoolProp.PropsSI('T', 'P', P_vec(m), 'Q',0, fluid);
%     rho_lSat(m) = CoolProp.PropsSI('D', 'P', P_vec(m), 'Q',0, fluid);
%     rho_vSat(m) = CoolProp.PropsSI('D', 'P', P_vec(m), 'Q',1, fluid);
% end
% hLsat_p_air  = griddedInterpolant(P_vec, h_lSat, 'cubic', 'linear');
% hVsat_p_air  = griddedInterpolant(P_vec, h_vSat, 'cubic', 'linear');
% Tsat_p_air   = griddedInterpolant(P_vec, T_Sat , 'cubic', 'linear');
% rhoLsat_p_air  = griddedInterpolant(P_vec, rho_lSat, 'cubic', 'linear');
% rhoVsat_p_air  = griddedInterpolant(P_vec, rho_vSat, 'cubic', 'linear');
% 
% %Temperature based on on enthalpy/pressure
% clearvars P_vec m N_p
% N_P = 3;
% N_h = 50;
% P_vec = linspace(P_min,P_max,N_P);
% T_min = -50 +273.15;
% T_max = 120 +273.15;
% h_min = min([CoolProp.PropsSI('H', 'P', P_min, 'T',T_min, fluid)    CoolProp.PropsSI('H', 'P', P_min, 'T',T_max, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_min, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_max, fluid)]);
% h_max = max([CoolProp.PropsSI('H', 'P', P_min, 'T',T_min, fluid)    CoolProp.PropsSI('H', 'P', P_min, 'T',T_max, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_min, fluid) CoolProp.PropsSI('H', 'P', P_max, 'T',T_max, fluid)]);
% h_vec = linspace(h_min,h_max,N_h);
% [P_MX,H_MX] = ndgrid(P_vec,h_vec);
% for m = 1:length(P_vec)
%     for n = 1:length(h_vec)
%         T_MX(m,n) = CoolProp.PropsSI('T', 'P', P_MX(m,n), 'H',H_MX(m,n), fluid);
%     end
% end
% T_ph_air  = griddedInterpolant(P_MX,H_MX, T_MX, 'cubic', 'linear');
% 
% if 1
%     % Test properties + speed
%     P_vec2 = linspace(P_min,P_max,3);
%     h_vec2 = linspace(h_min,h_max,500);
%     k1 = 0;
%     k2 = 0;
%     [T_IT] = deal(NaN*ones(length(P_vec2)*length(h_vec2),1));
%     [h_lSat_IT, h_vSat_IT, rho_lSat_IT, rho_vSat_IT, T_Sat_IT, T_IT] = deal(NaN*ones(length(P_vec2),1));
%     tic
%     for m = 1:length(P_vec2)
%         k1 = k1 + 1;
%         h_lSat_IT(k1) = hLsat_p_air(P_vec2(m));
%         h_vSat_IT(k1) = hVsat_p_air(P_vec2(m));
%         rho_lSat_IT(k1) = rhoLsat_p_air(P_vec2(m));
%         rho_vSat_IT(k1) = rhoVsat_p_air(P_vec2(m));
%         T_Sat_IT(k1) = Tsat_p_air(P_vec2(m));
%         for n = 1:length(h_vec2)
%             k2 = k2 + 1;
%             T_IT(k2) = T_ph_air(P_vec2(m),h_vec2(n));
%         end
%     end
%     time_IT = toc
%     
%     [T_CP] = deal(NaN*ones(length(P_vec2)*length(h_vec2),1));
%     [h_lSat_CP, h_vSat_CP, rho_lSat_CP, rho_vSat_CP, T_Sat_CP] = deal(NaN*ones(length(P_vec2),1));
%     tic
%     k1 = 0;
%     k2 = 0;
%     for m = 1:length(P_vec2)
%         k1 = k1 + 1;
%         h_lSat_CP(k1) = CoolProp.PropsSI('H', 'P', P_vec2(m), 'Q',0, fluid);
%         h_vSat_CP(k1) = CoolProp.PropsSI('H', 'P', P_vec2(m), 'Q',1, fluid);
%         rho_lSat_CP(k1) = CoolProp.PropsSI('D', 'P', P_vec2(m), 'Q',0, fluid);
%         rho_vSat_CP(k1) = CoolProp.PropsSI('D', 'P', P_vec2(m), 'Q',1, fluid);
%         T_Sat_CP(k1) = CoolProp.PropsSI('T', 'P', P_vec2(m), 'Q',0, fluid);
%         for n = 1:length(h_vec2)
%             k2 = k2 + 1;
%             T_CP(k2) = CoolProp.PropsSI('T', 'P', P_vec2(m), 'H',h_vec2(n), fluid);
%         end
%     end
%     time_CP = toc
%     
%     rate = time_CP/time_IT
%     figure
%     subplot(2,3,1)
%     plot(h_lSat_CP, h_lSat_IT,'o')
%     title(num2str(MARE(h_lSat_CP, h_lSat_IT)))
%     grid on
%     subplot(2,3,2)
%     plot(h_vSat_CP, h_vSat_IT,'o')
%     title(num2str(MARE(h_vSat_CP, h_vSat_IT)))
%     grid on
%     subplot(2,3,3)
%     plot(rho_lSat_CP, rho_lSat_IT,'o')
%     title(num2str(MARE(rho_lSat_CP, rho_lSat_IT)))
%     grid on
%     subplot(2,3,4)
%     plot(rho_vSat_CP, rho_vSat_IT,'o')
%     title(num2str(MARE(rho_vSat_CP, rho_vSat_IT)))
%     grid on
%     subplot(2,3,5)
%     plot(T_Sat_CP, T_Sat_IT,'o')
%     title(num2str(MAE(T_Sat_CP, T_Sat_IT)))
%     grid on
%     subplot(2,3,6)
%     plot(T_CP, T_IT,'o')
%     title(num2str(MAE(T_CP, T_IT)))
%     grid on
% end
% 
% save('PropsTab_air.mat','hLsat_p_air','hVsat_p_air','Tsat_p_air','rhoLsat_p_air', 'rhoVsat_p_air','T_ph_air')
% 
