clear all
close all
clc

%% RANGE OF REYNOLDS AND PRANDTL NUMBERS
if 1
    W = 0.191;
    L = 0.519-0.06;
    pitch_p = 0.0022;
    th_p = 0.0004;
    b = pitch_p-th_p;
    phi = 1.1407;
    Dh_i = 2*b/phi;
    N_c = 50;
    theta_i = 60*pi/180;
    disp_flag = 0;
    
    % R25FA
    fluid = 'R245fa';
    T_ref_range_vec = linspace(20,50,5)+273.15;
    P_ref_range_vec = linspace(5,12,5)*1e5;
    Mdot_ref_range_vec = linspace(10,150,5)/1000;
    i = 0;
    for jj = 1:length(Mdot_ref_range_vec)
        for kk = 1:length(P_ref_range_vec)
            for ii = 1:length(T_ref_range_vec)
                i = i+1;
                mu_ref_range(i) = CoolProp.PropsSI('V',        'T', T_ref_range_vec(ii), 'P', P_ref_range_vec(kk), fluid);
                Pr_ref_range(i) = CoolProp.PropsSI('Prandtl',  'T', T_ref_range_vec(ii), 'P', P_ref_range_vec(kk), fluid);
                k_ref_range(i)  = CoolProp.PropsSI('L',        'T', T_ref_range_vec(ii), 'P', P_ref_range_vec(kk), fluid);
                Mdot_ref_range(i) = Mdot_ref_range_vec(jj);
                G_ref_range(i) = Mdot_ref_range(i)/N_c/(W*b);
                Re_ref_range(i) = G_ref_range(i)*Dh_i/mu_ref_range(i);
            end
        end
    end
    % Pirobloc basic
    fluid = 'ICP_PiroblocBasic';
    T_htf_range_vec = linspace(50,150,20)+273.15;
    P_htf_range_vec = linspace(1,1,1)*1e5;
    Mdot_htf_range_vec = linspace(100,1000,100)/1000;
    i = 0;
    for jj = 1:length(Mdot_htf_range_vec)
        for kk = 1:length(P_htf_range_vec)
            for ii = 1:length(T_htf_range_vec)
                i = i+1;
                mu_htf_range(i) = PropsSI_ICP('V',        'T', T_htf_range_vec(ii), 'P', P_htf_range_vec(kk), fluid);
                cp_htf_range(i) = PropsSI_ICP('C',  'T', T_htf_range_vec(ii), 'P', P_htf_range_vec(kk), fluid);
                k_htf_range(i)  = PropsSI_ICP('L',        'T', T_htf_range_vec(ii), 'P', P_htf_range_vec(kk), fluid);
                Pr_htf_range(i) =  cp_htf_range(i)*mu_htf_range(i)/k_htf_range(i);
                Mdot_htf_range(i) = Mdot_htf_range_vec(jj);
                G_htf_range(i) = Mdot_htf_range(i)/N_c/(W*b);
                Re_htf_range(i) = G_htf_range(i)*Dh_i/mu_htf_range(i);
            end
        end
    end
    
    figure
    subplot(2,2,1)
    plot(Re_ref_range, 'o')
    subplot(2,2,2)
    plot(Pr_ref_range, 'o')
    subplot(2,2,3)
    plot(Re_htf_range, 'o')
    subplot(2,2,4)
    plot(Pr_htf_range, 'o')
end

%% BPHEX single-phase HTC comparison - R245FA

% Reference conditions
P = 10e5;
T = 70 +273.15;
fluid = 'R245fa';
mu = CoolProp.PropsSI('V',        'T', T, 'P', P, fluid);
mu_rat = 1;
Pr = CoolProp.PropsSI('Prandtl',  'T', T, 'P', P, fluid);
k  = CoolProp.PropsSI('L',        'T', T, 'P', P, fluid);
m_dot = 0.08;
W = 0.191;
L = 0.519-0.06;
pitch_p = 0.0022;
th_p = 0.0004;
b = pitch_p-th_p;
phi = 1.1407;
Dh_i = 2*b/phi;
N_c = 50;
theta_i = 45*pi/180;
disp_flag = 0;
G_i = m_dot/N_c/(W*b);

% influence of mass flux
G_vec = (10:10:150)/1000/N_c/(W*b);
Re_vec = G_vec*Dh_i/mu;
for i_G = 1:length(G_vec)
    [hConv_DittusBoelter_G(i_G) Nu_DittusBoelter_G(i_G) flag_DittusBoelter_G(i_G)]=     DittusBoelter_Pipe_HTC(mu,          Pr, k, G_vec(i_G), Dh_i, 0.4,               disp_flag);
    [hConv_Gnielinski_G(i_G)	Nu_Gnielinski_G(i_G) 	flag_Gnielinski_G(i_G)]=        Gnielinski_Pipe_HTC(mu,             Pr, k, G_vec(i_G), Dh_i, L,                 disp_flag);
    [hConv_Thonon_G(i_G)        Nu_Thonon_G(i_G)        flag_Thonon_G(i_G)]=         	Thonon_BPHEX_HTC(mu,                Pr, k, G_vec(i_G), Dh_i, theta_i,           disp_flag);
    [hConv_Martin1_G(i_G)       Nu_Martin1_G(i_G)       flag_Martin1_G(i_G)] =        	Martin1_BPHEX_HTC(mu,       mu_rat,	Pr, k, G_vec(i_G), Dh_i, theta_i,           disp_flag);
    [hConv_Wanniarachchi_G(i_G) Nu_Wanniarachchi_G(i_G) flag_Wanniarachchi_G(i_G)]=     Wanniarachchi_BPHEX_HTC(mu, mu_rat, Pr, k, G_vec(i_G), Dh_i, theta_i, phi,      disp_flag);
    [hConv_Heavner_G(i_G)       Nu_Heavner_G(i_G)       flag_Heavner_G(i_G)]=        	Heavner_BPHEX_HTC(mu,       mu_rat, Pr, k, G_vec(i_G), Dh_i, theta_i, phi,      disp_flag);
    [hConv_Muley_G(i_G)         Nu_Muley_G(i_G)         flag_Muley_G(i_G)]=         	Muley_BPHEX_HTC(mu,         mu_rat,	Pr, k, G_vec(i_G), Dh_i, theta_i, phi, L,   disp_flag);
    [hConv_Junqi_G(i_G)         Nu_Junqi_G(i_G)         flag_Junqi_G(i_G)]=          	Junqi_BPHEX_HTC(mu,                 Pr, k, G_vec(i_G), Dh_i, theta_i,           disp_flag);
    [hConv_DesideriHTF_G(i_G)	Nu_DesideriHTF_G(i_G)	flag_DesideriHTF_G(i_G)]=       DesideriHTF_BPHEX_HTC(mu,   mu_rat, Pr, k, G_vec(i_G), Dh_i,                    disp_flag);
    [hConv_Kim_G(i_G)           Nu_Kim_G(i_G)           flag_Kim_G(i_G)]=               Kim_BPHEX_HTC(mu,                   Pr, k, G_vec(i_G), Dh_i, theta_i,           disp_flag);
end

% influence of Dh
Dh_vec = 0.002:0.0001:0.005 ;
for i_D = 1:length(Dh_vec)
    [hConv_DittusBoelter_D(i_D) Nu_DittusBoelter_D(i_D) flag_DittusBoelter_D(i_D)]=     DittusBoelter_Pipe_HTC(mu,          Pr, k, G_i, Dh_vec(i_D), 0.4,               disp_flag);
    [hConv_Gnielinski_D(i_D)	Nu_Gnielinski_D(i_D) 	flag_Gnielinski_D(i_D)]=        Gnielinski_Pipe_HTC(mu,             Pr, k, G_i, Dh_vec(i_D), L,                 disp_flag);
    [hConv_Thonon_D(i_D)        Nu_Thonon_D(i_D)        flag_Thonon_D(i_D)]=         	Thonon_BPHEX_HTC(mu,                Pr, k, G_i, Dh_vec(i_D), theta_i,           disp_flag);
    [hConv_Martin1_D(i_D)       Nu_Martin1_D(i_D)       flag_Martin1_D(i_D)] =        	Martin1_BPHEX_HTC(mu,       mu_rat,	Pr, k, G_i, Dh_vec(i_D), theta_i,           disp_flag);
    [hConv_Wanniarachchi_D(i_D) Nu_Wanniarachchi_D(i_D) flag_Wanniarachchi_D(i_D)]=     Wanniarachchi_BPHEX_HTC(mu, mu_rat, Pr, k, G_i, Dh_vec(i_D), theta_i, phi,      disp_flag);
    [hConv_Heavner_D(i_D)       Nu_Heavner_D(i_D)       flag_Heavner_D(i_D)]=        	Heavner_BPHEX_HTC(mu,       mu_rat, Pr, k, G_i, Dh_vec(i_D), theta_i, phi,      disp_flag);
    [hConv_Muley_D(i_D)         Nu_Muley_D(i_D)         flag_Muley_D(i_D)]=         	Muley_BPHEX_HTC(mu,         mu_rat,	Pr, k, G_i, Dh_vec(i_D), theta_i, phi, L,   disp_flag);
    [hConv_Junqi_D(i_D)         Nu_Junqi_D(i_D)         flag_Junqi_D(i_D)]=          	Junqi_BPHEX_HTC(mu,                 Pr, k, G_i, Dh_vec(i_D), theta_i,           disp_flag);
    [hConv_DesideriHTF_D(i_D)	Nu_DesideriHTF_D(i_D)	flag_DesideriHTF_D(i_D)]=       DesideriHTF_BPHEX_HTC(mu,   mu_rat, Pr, k, G_i, Dh_vec(i_D),                    disp_flag);
    [hConv_Kim_D(i_D)           Nu_Kim_D(i_D)           flag_Kim_D(i_D)]=               Kim_BPHEX_HTC(mu,                   Pr, k, G_i, Dh_vec(i_D), theta_i,           disp_flag);
end


% influence of theta
theta_vec = (30:1:60) ;
for i_T = 1:length(theta_vec)
    [hConv_DittusBoelter_T(i_T) Nu_DittusBoelter_T(i_T) flag_DittusBoelter_T(i_T)]=     DittusBoelter_Pipe_HTC(mu,          Pr, k, G_i, Dh_i, 0.4,               disp_flag);
    [hConv_Gnielinski_T(i_T)	Nu_Gnielinski_T(i_T) 	flag_Gnielinski_T(i_T)]=        Gnielinski_Pipe_HTC(mu,             Pr, k, G_i, Dh_i, L,                 disp_flag);
    [hConv_Thonon_T(i_T)        Nu_Thonon_T(i_T)        flag_Thonon_T(i_T)]=         	Thonon_BPHEX_HTC(mu,                Pr, k, G_i, Dh_i, theta_vec(i_T)*pi/180,           disp_flag);
    [hConv_Martin1_T(i_T)       Nu_Martin1_T(i_T)       flag_Martin1_T(i_T)] =        	Martin1_BPHEX_HTC(mu,       mu_rat,	Pr, k, G_i, Dh_i, theta_vec(i_T)*pi/180,           disp_flag);
    [hConv_Wanniarachchi_T(i_T) Nu_Wanniarachchi_T(i_T) flag_Wanniarachchi_T(i_T)]=     Wanniarachchi_BPHEX_HTC(mu, mu_rat, Pr, k, G_i, Dh_i, theta_vec(i_T)*pi/180, phi,      disp_flag);
    [hConv_Heavner_T(i_T)       Nu_Heavner_T(i_T)       flag_Heavner_T(i_T)]=        	Heavner_BPHEX_HTC(mu,       mu_rat, Pr, k, G_i, Dh_i, theta_vec(i_T)*pi/180, phi,      disp_flag);
    [hConv_Muley_T(i_T)         Nu_Muley_T(i_T)         flag_Muley_T(i_T)]=         	Muley_BPHEX_HTC(mu,         mu_rat,	Pr, k, G_i, Dh_i, theta_vec(i_T)*pi/180, phi, L,   disp_flag);
    [hConv_Junqi_T(i_T)         Nu_Junqi_T(i_T)         flag_Junqi_T(i_T)]=          	Junqi_BPHEX_HTC(mu,                 Pr, k, G_i, Dh_i, theta_vec(i_T)*pi/180,           disp_flag);
    [hConv_DesideriHTF_T(i_T)	Nu_DesideriHTF_T(i_T)	flag_DesideriHTF_T(i_T)]=       DesideriHTF_BPHEX_HTC(mu,   mu_rat, Pr, k, G_i, Dh_i,                    disp_flag);
    [hConv_Kim_T(i_T)           Nu_Kim_T(i_T)           flag_Kim_T(i_T)]=               Kim_BPHEX_HTC(mu,                   Pr, k, G_i, Dh_i, theta_vec(i_T)*pi/180,           disp_flag);
end


% Figures
col2eval = { 'Gnielinski', 'Thonon', 'Martin1', 'Wanniarachchi','Muley', 'Junqi','Kim'};%, 'Heavner',  'DesideriHTF',  'DittusBoelter',
figure
subplot(1,3,1)
hold on
ii = 0;
for i = 1:length(col2eval)
    ii = ii+1;
    eval(['LG(ii) = plot(Re_vec, hConv_' col2eval{i} '_G, ''o-'', ''MarkerSize'', 7);'])
    legG{ii} = col2eval{i};
    eval(['plot(Re_vec(flag_' col2eval{i} '_G== 0), hConv_' col2eval{i} '_G(flag_' col2eval{i} '_G== 0), ''*c'');']) % problem with reynolds
    eval(['plot(Re_vec(flag_' col2eval{i} '_G==-1), hConv_' col2eval{i} '_G(flag_' col2eval{i} '_G==-1), ''*m'');']) % problem with Prandtl
    eval(['plot(Re_vec(flag_' col2eval{i} '_G==-2), hConv_' col2eval{i} '_G(flag_' col2eval{i} '_G==-2), ''*r'');']) % problem with reynolds and Prandtl
end
legend(LG, legG, 'location', 'northwest')
grid on
xlabel('Re')
ylabel('hConv')

subplot(1,3,2)
hold on
ii = 0;
for i = 1:length(col2eval)
    ii = ii+1;
    eval(['LD(ii) = plot(Dh_vec, hConv_' col2eval{i} '_D, ''o-'', ''MarkerSize'', 7);'])
    legD{ii} = col2eval{i};
    eval(['plot(Dh_vec(flag_' col2eval{i} '_D== 0), hConv_' col2eval{i} '_D(flag_' col2eval{i} '_D== 0), ''*c'');']) % problem with reynolds
    eval(['plot(Dh_vec(flag_' col2eval{i} '_D==-1), hConv_' col2eval{i} '_D(flag_' col2eval{i} '_D==-1), ''*m'');']) % problem with Prandtl
    eval(['plot(Dh_vec(flag_' col2eval{i} '_D==-2), hConv_' col2eval{i} '_D(flag_' col2eval{i} '_D==-2), ''*r'');']) % problem with reynolds and Prandtl
end
legend(LD, legD, 'location', 'northwest')
grid on
xlabel('Dh')
ylabel('hConv')

subplot(1,3,3)
hold on
ii = 0;
for i = 1:length(col2eval)
    ii = ii+1;
    eval(['LT(ii) = plot(theta_vec, hConv_' col2eval{i} '_T, ''o-'', ''MarkerSize'', 7);'])
    legT{ii} = col2eval{i};
    eval(['plot(theta_vec(flag_' col2eval{i} '_T== 0), hConv_' col2eval{i} '_T(flag_' col2eval{i} '_T== 0), ''*c'');']) % problem with reynolds
    eval(['plot(theta_vec(flag_' col2eval{i} '_T==-1), hConv_' col2eval{i} '_T(flag_' col2eval{i} '_T==-1), ''*m'');']) % problem with Prandtl
    eval(['plot(theta_vec(flag_' col2eval{i} '_T==-2), hConv_' col2eval{i} '_T(flag_' col2eval{i} '_T==-2), ''*r'');']) % problem with reynolds and Prandtl
end
legend(LT, legT, 'location', 'northwest')
grid on
xlabel('theta')
ylabel('hConv')



%% BPHEX boiling HTC comparison