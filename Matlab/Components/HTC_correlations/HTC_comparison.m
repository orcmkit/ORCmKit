clear all
close all
clc

%% RANGE OF REYNOLDS AND PRANDTL NUMBERS
if 0
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

%% BPHEX SINGLE-PHASE HTC comparison - R245FA
if 0
% Reference conditions
P = 10e5;
T = 50 +273.15;
fluid = 'R245fa';
T_sat = CoolProp.PropsSI('T',        'Q', 0.5, 'P', P, fluid)-273.15;
mu = CoolProp.PropsSI('V',        'T', T, 'P', P, fluid);
mu_rat = 1;
Pr = CoolProp.PropsSI('Prandtl',  'T', T, 'P', P, fluid);
k  = CoolProp.PropsSI('L',        'T', T, 'P', P, fluid);
m_dot = 0.08;
 L = 0.519-0.0603;
 W = 0.191;
 pitch_p = 2.2/1000;
 th_p = 0.4/1000;
 b = pitch_p-th_p;
 phi = 1.1414;
 Dh_i = 2*b/phi;
 N_c = 50;
 theta_i = 60*pi/180;
 pitch_co = 7.19360908/1000;
 

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
figure('Name', [num2str(T-273.15) ' - ' num2str(P/1e5)])
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
end

%% BPHEX CONDESNEING HTC comparison
if 1
    P = 4e5;
    fluid = 'R245fa';
    x_i = 0.1;
    T_sat = CoolProp.PropsSI('T',     'Q', 0.5, 'P', P, fluid)-273.15;
    mu_l    = CoolProp.PropsSI('V',	'Q', 0, 'P', P, fluid);
    k_l     = CoolProp.PropsSI('L',	'Q', 0, 'P', P, fluid);
    Pr_l = CoolProp.PropsSI('Prandtl',  'Q', 0, 'P', P, fluid);
    rho_l = CoolProp.PropsSI('D',  'Q', 0, 'P', P, fluid);
    rho_v = CoolProp.PropsSI('D',  'Q', 1, 'P', P, fluid);
    i_fg = CoolProp.PropsSI('H',  'Q', 1, 'P', P, fluid) - CoolProp.PropsSI('H',  'Q', 0, 'P', P, fluid);
    p_star = P/CoolProp.PropsSI('Pcrit',  'Q', 1, 'P', P, fluid);
    m_dot = 0.12;
    
    L = 0.250-0.028757;
    W = 0.1128;
    pitch_p = 0.00231;
    th_p = 0.0003;
    b = pitch_p-th_p;
    phi = 1.1620;
    Dh_i = 2*b/phi;
    N_c = 20;
    theta_i = 60*pi/180;
    pitch_co_i = 7.45870973/1000;
    disp_flag = 0;
    G_i = m_dot/N_c/(W*b);
    
    
    % Impact of quality
    x_vec = 0.00001:0.03:0.999;
    for i_x = 1:length(x_vec)
        [h_Han_X(i_x),      Nu_Han_X(i_x),      flag_Han_X(i_x)]  =     Han_Cond_BPHEX_HTC(x_vec(i_x), mu_l, k_l, Pr_l, rho_l, rho_v, G_i, Dh_i, pitch_co_i, theta_i, disp_flag);
        [h_Longo_X(i_x),    Nu_Longo_X(i_x),    flag_Longo_X(i_x)]  = Longo_Cond_BPHEX_HTC(x_vec(i_x), mu_l, k_l, Pr_l, rho_l, rho_v, i_fg, G_i, Dh_i, 2,    phi, L, disp_flag);
        [h_Shah_X(i_x),     Nu_Shah_X(i_x),     flag_Shah_X(i_x)]   = Shah_Cond_pipe_HTC(x_vec(i_x), mu_l, k_l, Pr_l, p_star, G_i, Dh_i, disp_flag);
    
    end
    
    figure
    hold on
    plot(x_vec, h_Han_X, 'o')
    plot(x_vec, h_Longo_X, 'o')
    plot(x_vec, h_Shah_X, 'o')
    hold off
    legend('Han', 'Longo', 'Shah')
end

%% TUBE CONDESNEING HTC comparison
if 0
    P = 3e5;
    fluid = 'R245fa';
    x_i = 0.1;
    T_sat = CoolProp.PropsSI('T',     'Q', 0.5, 'P', P, fluid)-273.15;
    mu_l    = CoolProp.PropsSI('V',	'Q', 0, 'P', P, fluid);
    mu_v    = CoolProp.PropsSI('V',	'Q', 1, 'P', P, fluid);
    k_l     = CoolProp.PropsSI('L',	'Q', 0, 'P', P, fluid);
    Pr_l = CoolProp.PropsSI('Prandtl',  'Q', 0, 'P', P, fluid);
    rho_l = CoolProp.PropsSI('D',  'Q', 0, 'P', P, fluid);
    rho_v = CoolProp.PropsSI('D',  'Q', 1, 'P', P, fluid);
    i_fg = CoolProp.PropsSI('H',  'Q', 1, 'P', P, fluid) - CoolProp.PropsSI('H',  'Q', 0, 'P', P, fluid);
    p_star = P/CoolProp.PropsSI('Pcrit',  'Q', 1, 'P', P, fluid);
    m_dot = 0.06;
      
    Dh_i = 0.0085;
    N_c = 13;
    G_i = m_dot/N_c/(pi/4*Dh_i^2);
    DT_wall_i = 10;
    disp_flag = 0;
    
    % Impact of quality
    x_vec = 0.00001:0.03:0.999;
    for i_x = 1:length(x_vec)
        [h_Cav_X(i_x),      Nu_Cav_X(i_x),      flag_Cav_X(i_x)]  =     Cavallini_Cond_pipe_HTC(x_vec(i_x), mu_l, mu_v, rho_l, rho_v, k_l, Pr_l, i_fg, DT_wall_i, G_i, Dh_i, disp_flag) ;
        [h_Shah_X(i_x),     Nu_Shah_X(i_x),     flag_Shah_X(i_x)]   =   Shah_Cond_pipe_HTC(x_vec(i_x), mu_l, k_l, Pr_l, p_star, G_i, Dh_i, disp_flag);   
    end
    
    figure
    hold on
    plot(x_vec, h_Cav_X, 'o')
    plot(x_vec, h_Shah_X, 'o')
    hold off
    legend('Cav', 'Shah')
    
    % Impact of wall temperature
    DT_vec = 1:1:10;
    figure
    hold on
    for i_x = 1:length(x_vec)
        [h_Shah_X(i_x),     Nu_Shah_X(i_x),     flag_Shah_X(i_x)]   =   Shah_Cond_pipe_HTC(x_vec(i_x), mu_l, k_l, Pr_l, p_star, G_i, Dh_i, disp_flag);
    end
    plot(x_vec, h_Shah_X, 'o')
    for i_DT = 1:length(DT_vec)
        for i_x = 1:length(x_vec)           
            [h_Cav_X(i_x),      Nu_Cav_X(i_x),      flag_Cav_X(i_x)]  =     Cavallini_Cond_pipe_HTC(x_vec(i_x), mu_l, mu_v, rho_l, rho_v, k_l, Pr_l, i_fg, DT_vec(i_DT), G_i, Dh_i, disp_flag) ;
        end
        plot(x_vec, h_Cav_X, 'o')
        clearvars h_Cav_X Nu_Cav_X flag_Cav_X
    end
    
    
end

%% BPHEX BOILING HTC comparison
if 0
    P = 10e5;
    fluid = 'R245fa';
    x_i = 0.1;
    T_sat = CoolProp.PropsSI('T',     'Q', 0.5, 'P', P, fluid)-273.15;
    mu_l    = CoolProp.PropsSI('V',	'Q', 0, 'P', P, fluid);
    mu_v    = CoolProp.PropsSI('V',	'Q', 1, 'P', P, fluid);
    k_l     = CoolProp.PropsSI('L',	'Q', 0, 'P', P, fluid);
    Pr_l = CoolProp.PropsSI('Prandtl',  'Q', 0, 'P', P, fluid);
    rho_l = CoolProp.PropsSI('D',  'Q', 0, 'P', P, fluid);
    rho_v = CoolProp.PropsSI('D',  'Q', 1, 'P', P, fluid);    
    i_fg = CoolProp.PropsSI('H',  'Q', 1, 'P', P, fluid) - CoolProp.PropsSI('H',  'Q', 0, 'P', P, fluid);
    p_star = P/CoolProp.PropsSI('Pcrit',  'Q', 1, 'P', P, fluid);
    m_dot = 0.09;       
    
    L = 0.519-0.0603;
    W = 0.191;
    pitch_p = 2.2/1000;
    th_p = 0.4/1000;
    b = pitch_p-th_p;
    phi = 1.1414;
    Dh_i = 2*b/phi;
    N_c = 50;
    theta_i = 60*pi/180;
    pitch_co = 7.19360908/1000;
    disp_flag = 0;
    G_i = 55; %m_dot/N_c/(W*b);
    
    Qdot = 1000; 
    DTlog = 15;
    h_conv_h = 400;

    % Impact of quality
    figure
    hold on
    x_vec = 0.00001:0.03:0.999;
    for i_x = 1:length(x_vec)
        rho = CoolProp.PropsSI('D',  'Q', x_vec(i_x), 'P', P, fluid);
        P_crit = CoolProp.PropsSI('Pcrit', 'Q', 1, 'P', P, fluid);
        sigma = CoolProp.PropsSI('I',  'Q', x_vec(i_x), 'P', P, fluid);
        sigma_l = CoolProp.PropsSI('I', 'Q', 0, 'P', 0.1*P_crit, fluid)*1000; %mN/m
        %k = CoolProp.PropsSI('L',	'Q',  x_vec(i_x), 'P', P, fluid);
        p_star = P/P_crit;
        M = 1e3*CoolProp.PropsSI('M', 'Q', 1, 'P', P, fluid);
        [h_HanB_X(i_x),    	Nu_HanB_X(i_x),     flag_HanB_X(i_x)]    = 	Han_Boiling_BPHEX_HTC(x_vec(i_x), mu_l, k_l, Pr_l, rho_l, rho_v,  i_fg, G_i, DTlog, Qdot, h_conv_h, Dh_i, theta_i, pitch_co, disp_flag);
        [h_Amalfi_X(i_x),  	Nu_Amalfi_X(i_x),	flag_Amalfi_X(i_x)]	 = 	Amalfi_Boiling_BPHEX_HTC(x_vec(i_x), rho_l, rho_v, rho, mu_v, mu_l, k_l, sigma, i_fg, G_i, Dh_i, theta_i, Qdot, DTlog, h_conv_h, disp_flag);
        [h_Junqi_X(i_x),   	Nu_Junqi_X(i_x),  	flag_Junqi_X(i_x)]   =	Junqi_Boiling_BPHEX_HTC(x_vec(i_x), mu_l, k_l, Pr_l, rho_l, rho_v,  i_fg, G_i, DTlog, Qdot, h_conv_h, Dh_i, disp_flag);
        [h_Cooper_X(i_x),  	Nu_Cooper_X(i_x),  	flag_Cooper_X(i_x)]  =	Cooper_Boiling_HTC(Dh_i, k_l, P, P_crit, M, h_conv_h, DTlog, Qdot);
        [h_Longo_X(i_x),  	Nu_Longo_X(i_x),  	flag_Longo_X(i_x)]   =  Longo_Boiling_BPHEX_HTC(x_vec(i_x), k_l, rho_l, Pr_l, mu_l, rho_v, P_crit, p_star, sigma_l, G_i, Dh_i, phi, h_conv_h, DTlog, Qdot, fluid);
        %[h_Desideri_X(i_x),	Nu_Desideri_X(i_x),	flag_Desideri_X(i_x)]=  Desideri_Boiling_BPHEX_HTC(x_vec(i_x), rho_l, rho_v, rho, mu_l, k_l, sigma, G_i, Dh_i, disp_flag);
    end
    
    figure
    hold on
    plot(x_vec, h_HanB_X, 'o')
    plot(x_vec, h_Amalfi_X, 'o')
    plot(x_vec, h_Junqi_X, 'o')
    plot(x_vec, h_Cooper_X, 'o')
    plot(x_vec, h_Longo_X, 'o')
    %plot(x_vec, h_Desideri_X, 'o')
    hold off
    legend('Han', 'Amalfi', 'Junqi', 'Cooper', 'Longo')%, 'Desideri')
end
