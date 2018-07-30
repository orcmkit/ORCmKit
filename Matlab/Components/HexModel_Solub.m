function [out,TS] = HexModel_Solub(fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems
% Remi Dickes - 27/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be


%% DEMONSTRATION CASE
if nargin == 0    % Define a demonstration case if HexModel.mat is not executed externally
    clear all
    close all
    clc
    exemple = 'EV';
    if strcmp(exemple, 'EV')
        % Supply fluid conditions
        fluid_h = 'air';                                       	% Nature of the hot fluid           [-]
        m_dot_h = 2;                                          	% Mass flow rat of the hot fluid    [kg/s]
        P_h_su =  1e5;                                              % Supply pressure of the hot fluid  [Pa]
        T_h_su = 220+273.15;
        param.H.type = 'H';
        param.H.solub = 0;
        param.H.fluid_lub = 'ICP_RL_32_3MAF';
        param.H.C_oil = 0.02;
        h_ref = CoolProp.PropsSI('H','P',P_h_su,'T',T_h_su, fluid_h);      % Supply h or T of the cold fluid  	[J/kg pr K]
        h_lub = PropsSI_ICP('H', 'T', T_h_su, 'P', P_h_su, param.H.fluid_lub);
        in_h_su = (1-param.H.C_oil)*h_ref + param.H.C_oil*h_lub;
        load('Ratio_RhoIdeal_R245Fa_POE.mat')
        param.H.fit_ratio_rho = fit_ratio_rho;
        load('DTP_zeta_R245fa_POE.mat')
        param.C.fit_DTP_zeta = fit_DT_Pbar_zeta2;
        param.C_fit = 1;
        
        fluid_c = 'water';                                                     % Nature of the cold fluid        	[-]
        m_dot_c = 0.1;                                                          % Mass flow rat of the cold fluid  	[kg/s]
        T_c_su = 10 + 273.15;
        P_c_su = 5e5;                                                          % Supply pressure of the cold fluid	[Pa]
        param.C.type = 'H';
        param.C.solub = 0;
        param.C.fluid_lub = 'ICP_RL_32_3MAF';
        param.C.C_oil = 0.02;
        h_ref = CoolProp.PropsSI('H','P',P_c_su,'T',T_c_su, fluid_c);      % Supply h or T of the cold fluid  	[J/kg pr K]
        h_lub = PropsSI_ICP('H', 'T', T_c_su, 'P', P_c_su, param.C.fluid_lub);
        in_c_su = (1-param.C.C_oil)*h_ref + param.C.C_oil*h_lub;
        load('Ratio_RhoIdeal_R245Fa_POE.mat')
        param.C.fit_ratio_rho = fit_ratio_rho;
        param.H.fit_DTP_zeta = fit_DT_Pbar_zeta2;
        
        % Geometrical parameters
        param.theta = 55*pi/180;
        param.H.A_tot = 5;
        param.C.A_tot = 5;
        param.H.V_tot = 0.001026;
        param.C.V_tot = 0.00108;
        param.L_hex = 0.221243;
        param.pitch_co = 0.00745870973;
        param.phi =  1.162035402618661;
        param.C.Dh = 0.0032;
        param.H.Dh = 0.003459447096828;
        param.C.CS = 2.267280000000000e-04;
        param.H.CS = 2.267280000000000e-04;
        param.C.n_canals = 20;
        param.H.n_canals = 19;
        param.H.fin = 'none';
        param.C.fin = 'none';
        param.k0 = 13.022;
        param.k1 = 0.01667;
        param.C.R_fooling = 0.20/10000;
        param.H.R_fooling = 0.20/10000;
        param.t = 0.3/1000;
           
        param.typeHEX = 'CounterFlow';
        param.n_disc = 5;
        param.modelType = 'hConvCor';
        param.H.correlation.type_1phase_l = 'Martin1_BPHEX';
        param.H.correlation.type_1phase_v = 'Martin1_BPHEX';
        param.H.correlation.type_2phase_cd = 'Han_cond_BPHEX';
        param.H.correlation.type_2phase_ev = 'Han_boiling_BPHEX';
        param.C.correlation.type_1phase_l = 'Martin1_BPHEX';
        param.C.correlation.type_1phase_v = 'Martin1_BPHEX';
        param.C.correlation.type_2phase_ev = 'Han_boiling_BPHEX';
        param.C.correlation.type_2phase_cd = 'Han_cond_BPHEX';
        
        param.H.correlation.type_void_fraction = 'Homogenous';
        param.C.correlation.type_void_fraction = 'Homogenous';
        
        param.C.correlation.dry_out_incipient = 'UD';
        param.C.x_di = 0.8;
        
        param.hugmark_simplified = 1;
        
        param.displayResults = 1;
        param.displayTS = 1;
        param.generateTS = 1;
    end
    
    if strcmp(exemple, 'REC')
        % Supply fluid conditions
        fluid_h = 'R245fa';                                       	% Nature of the hot fluid           [-]
        m_dot_h = 0.05;                                          	% Mass flow rat of the hot fluid    [kg/s]
        P_h_su =  3e5;                                              % Supply pressure of the hot fluid  [Pa]
        T_h_su = 70+273.15;
        param.H.type = 'H';
        param.H.solub = 0;
        param.H.fluid_lub = 'ICP_RL_32_3MAF';
        param.H.C_oil = 0.02;
        h_ref = CoolProp.PropsSI('H','P',P_h_su,'T',T_h_su, fluid_h);      % Supply h or T of the cold fluid  	[J/kg pr K]
        h_lub = PropsSI_ICP('H', 'T', T_h_su, 'P', P_h_su, param.H.fluid_lub);
        in_h_su = (1-param.H.C_oil)*h_ref + param.H.C_oil*h_lub;
        load('Ratio_RhoIdeal_R245Fa_POE.mat')
        param.H.fit_ratio_rho = fit_ratio_rho;
        load('DTP_zeta_R245fa_POE.mat')
        param.C.fit_DTP_zeta = fit_DT_Pbar_zeta2;
        param.C_fit = 1;
        
        fluid_c = 'R245fa';                                                     % Nature of the cold fluid        	[-]
        m_dot_c = 0.06;                                                          % Mass flow rat of the cold fluid  	[kg/s]
        T_c_su = 30 + 273.15;
        P_c_su = 10e5;                                                          % Supply pressure of the cold fluid	[Pa]
        param.C.type = 'H';
        param.C.solub = 0;
        param.C.fluid_lub = 'ICP_RL_32_3MAF';
        param.C.C_oil = 0.02;
        h_ref = CoolProp.PropsSI('H','P',P_c_su,'T',T_c_su, fluid_c);      % Supply h or T of the cold fluid  	[J/kg pr K]
        h_lub = PropsSI_ICP('H', 'T', T_c_su, 'P', P_c_su, param.C.fluid_lub);
        in_c_su = (1-param.C.C_oil)*h_ref + param.C.C_oil*h_lub;
        load('Ratio_RhoIdeal_R245Fa_POE.mat')
        param.C.fit_ratio_rho = fit_ratio_rho;
        param.H.fit_DTP_zeta = fit_DT_Pbar_zeta2;
        
        % Geometrical parameters
        param.theta = 60*pi/180;
        param.H.A_tot = 1.102;
        param.C.A_tot = 1.102;
        param.H.V_tot = 0.001026;
        param.C.V_tot = 0.00108;
        param.L_hex = 0.221243;
        param.pitch_co = 0.00745870973;
        param.phi =  1.162035402618661;
        param.C.Dh = 0.0032;
        param.H.Dh = 0.003459447096828;
        param.C.CS = 2.267280000000000e-04;
        param.H.CS = 2.267280000000000e-04;
        param.C.n_canals = 20;
        param.H.n_canals = 19;
        param.H.fin = 'none';
        param.C.fin = 'none';
        param.k0 = 13.022;
        param.k1 = 0.01667;
        param.C.R_fooling = 0.20/10000;
        param.H.R_fooling = 0.20/10000;
        param.t = 0.3/1000;
           
        param.typeHEX = 'CounterFlow';
        param.n_disc = 5;
        param.modelType = 'hConvCor';
        param.H.correlation.type_1phase_l = 'Martin1_BPHEX';
        param.H.correlation.type_1phase_v = 'Martin1_BPHEX';
        param.H.correlation.type_2phase_cd = 'Han_cond_BPHEX';
        param.H.correlation.type_2phase_ev = 'Han_boiling_BPHEX';
        param.C.correlation.type_1phase_l = 'Martin1_BPHEX';
        param.C.correlation.type_1phase_v = 'Martin1_BPHEX';
        param.C.correlation.type_2phase_ev = 'Han_boiling_BPHEX';
        param.C.correlation.type_2phase_cd = 'Han_cond_BPHEX';
        
        param.H.correlation.type_void_fraction = 'Homogenous';
        param.C.correlation.type_void_fraction = 'Homogenous';
        
        param.C.correlation.dry_out_incipient = 'UD';
        param.C.x_di = 0.8;
        
        param.hugmark_simplified = 1;
        
        param.displayResults = 1;
        param.displayTS = 1;
        param.generateTS = 1;
    end

    if strcmp(exemple, 'CD')
    end
end

tstart_hex = tic;
 
%% HEAT EXCHANGER MODELLING
param.H.G = m_dot_h/param.H.CS/param.H.n_canals;
param.C.G = m_dot_c/param.C.CS/param.C.n_canals;

% Evaluation of the hot fluid (HF) supply conditions
[h_h_l, h_h_v, T_h_su] = supply_conditions_HEX_Solub(in_h_su, P_h_su, fluid_h, param.H);

% Evaluation of the cold fluid (CF) supply conditions
[h_c_l, h_c_v, T_c_su] = supply_conditions_HEX_Solub(in_c_su, P_c_su, fluid_c, param.C);


% Switch inputs if temperature inversion
flag_reverse = 0;
if T_h_su<T_c_su
    param2 = rmfield(param, {'H', 'C'});
    [fluid_h_int, P_h_su_int, in_h_su_int, m_dot_h_int, T_h_su_int, h_h_l_int, h_h_v_int, H_int] = deal(fluid_h, P_h_su, in_h_su, m_dot_h, T_h_su, h_h_l, h_h_v, param.H);
    [fluid_h, P_h_su, in_h_su, m_dot_h, T_h_su, h_h_l, h_h_v, param2.H] = deal(fluid_c, P_c_su, in_c_su, m_dot_c, T_c_su, h_c_l, h_c_v,  param.C);
    [fluid_c, P_c_su, in_c_su, m_dot_c, T_c_su, h_c_l, h_c_v, param2.C] = deal(fluid_h_int, P_h_su_int, in_h_su_int, m_dot_h_int, T_h_su_int, h_h_l_int, h_h_v_int, H_int);
    clear param; param = param2;
    flag_reverse = 1;
end

if (T_h_su-T_c_su)>1e-2  && m_dot_h  > 0 && m_dot_c > 0; % Check if the operating conditions permit a viable heat transfer
    if isfield(param, 'Q_dot_max')
        Q_dot_max = param.Q_dot_max;
        out_max = HEX_profile_Solub(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, param, h_h_l, h_h_v, h_c_l, h_c_v);
        pinch_Qdot_max = out_max.pinch;
    else
        [Q_dot_max, pinch_Qdot_max] = HEX_Qdotmax_Solub(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param, h_h_l, h_h_v, h_c_l, h_c_v); % compute maximal heat transfer
    end
    lb = 0;
    ub = Q_dot_max;
    f = @(Q_dot) HEX_hConvCorSolub_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su,  Q_dot, param, h_h_l, h_h_v, h_c_l, h_c_v);
    %f(ub)
    %figure
    %Qdot_test = linspace(0, Q_dot_max, 100);
    %hold all
    %for j = 1:length(Qdot_test)
    %   plot(Qdot_test(j), f(Qdot_test(j)), 'o')
    %end
    
    if f(ub) > 0
        Q_dot_eff = ub; % HEX so oversized that the effective heat power is equal to Q_dot_max
    else
        Q_dot_eff = zeroBrent ( lb, ub, eps, 1e-9, f ); % Solver driving residuals of HEX_hConvVar_res to zero
    end
else
    Q_dot_eff = 0;
    Q_dot_max = 0;
    pinch_Qdot_max = -inf;
end
    
out = HEX_hConvCorSolub(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, real(Q_dot_eff), param, h_h_l, h_h_v, h_c_l, h_c_v);
out.H.fluid = fluid_h;
out.C.fluid = fluid_c;
out.Q_dot_tot = real(Q_dot_eff);
out.epsilon_th = real(Q_dot_eff)/Q_dot_max;
out.H.h_ex = out.H.H_vec(1);
out.C.h_ex = out.C.H_vec(end);
out.H.T_ex = out.H.T_vec(1);
out.C.T_ex = out.C.T_vec(end);

% Wall temperatur calculation
for i = 1:length(out.H.T_vec)
    if i == 1
        out.T_wall(i) = out.H.T_vec(i) - out.Qdot_vec(i)/out.H.A_vec(i)/out.H.hConv_vec(i)/out.H.eff_vec(i);
    elseif i==length(out.H.T_vec)
        out.T_wall(i) = out.H.T_vec(i) - out.Qdot_vec(i-1)/out.H.A_vec(i-1)/out.H.hConv_vec(i-1)/out.H.eff_vec(i-1);
    else
        out.T_wall(i) = out.H.T_vec(i) - (0.5*out.Qdot_vec(i)/out.H.A_vec(i)/out.H.hConv_vec(i)/out.H.eff_vec(i) + 0.5*out.Qdot_vec(i-1)/out.H.A_vec(i-1)/out.H.hConv_vec(i-1)/out.H.eff_vec(i-1));
    end
end

% Averaged heat transfer coefficient
out.H.hConv_liq_mean = mean(out.H.hConv_vec(strcmp(out.H.type_zone, 'liq')));
out.H.hConv_tp_mean = mean(out.H.hConv_vec(strcmp(out.H.type_zone, 'tp')));
out.H.hConv_vap_mean = mean(out.H.hConv_vec(strcmp(out.H.type_zone, 'vap')));
out.C.hConv_liq_mean = mean(out.C.hConv_vec(strcmp(out.C.type_zone, 'liq')));
out.C.hConv_tp_mean = mean(out.C.hConv_vec(strcmp(out.C.type_zone, 'tp')));
out.C.hConv_vap_mean = mean(out.C.hConv_vec(strcmp(out.C.type_zone, 'vap')));

% Entropy vector calculation
if param.generateTS
    out.H.s_vec = generate_s_vec(param.H, out.H, fluid_h);
    out.C.s_vec = generate_s_vec(param.C, out.C, fluid_c);
end

% Mass calculation
out.H.V_vec = param.H.V_tot*(out.H.A_vec./sum(out.H.A_vec));
out.C.V_vec = param.C.V_tot*(out.H.A_vec./sum(out.H.A_vec));
[out.H.M_vec, out.H.M_mf_vec, out.H.M_lub_vec, out.H.alpha_mean_vec, out.H.M_tot, out.H.M_mf_tot, out.H.M_lub_tot] = generate_M_vec(param.H, out.H, fluid_h);
[out.C.M_vec, out.C.M_mf_vec, out.C.M_lub_vec, out.C.alpha_mean_vec, out.C.M_tot, out.C.M_mf_tot, out.C.M_lub_tot] = generate_M_vec(param.C, out.C, fluid_c);


% Flag evaluation
if (T_h_su-T_c_su)<1e-2  && (T_h_su-T_c_su >0)  && m_dot_h  > 0 && m_dot_c > 0
    out.flag = 3;
else    
    if out.resA <1e-4
        out.flag = 1;
    else
        if abs(pinch_Qdot_max) < 1e-3
            if Q_dot_eff == Q_dot_max
                out.flag = 2;
            else
                out.flag = -1;
            end
            
        else
            out.flag = -2;
        end
    end
end
out.flag_reverse = flag_reverse;

if flag_reverse
    names_out = fieldnames(rmfield(out,{'H','C'}));
    for i = 1:length(names_out)
        eval(['out2.' names_out{i} ' = flip(out.' names_out{i} ');'])
    end
    names_out_H = fieldnames(out.H);
    for i = 1:length(names_out_H)
        if strcmp(names_out_H{i}, 'fluid')
            eval(['out2.C.' names_out_H{i} ' = (out.H.' names_out_H{i} ');'])
        else
            eval(['out2.C.' names_out_H{i} ' = flip(out.H.' names_out_H{i} ');'])
        end
    end
    names_out_C = fieldnames(out.C);
    for i = 1:length(names_out_C)
        if strcmp(names_out_C{i}, 'fluid')
            eval(['out2.H.' names_out_C{i} ' = (out.C.' names_out_C{i} ');'])
        else
            eval(['out2.H.' names_out_C{i} ' = flip(out.C.' names_out_C{i} ');'])
        end
    end
    out2.x_flux = (out2.H.H_vec-out2.H.H_vec(1))./(out2.H.H_vec(end)-out2.H.H_vec(1));
    clear out; out = out2;
end
out.x_geom = [0 cumsum(out.H.V_vec)./param.H.V_tot];
out.time = toc(tstart_hex);

%% TS DIAGRAM and DISPLAY
% Generate the output variable TS
if param.generateTS
    TS.T_h = out.H.T_vec;
    TS.T_c = out.C.T_vec;
    TS.T_w = out.T_wall;
    %TS.s_h = out.H.s_vec;
    %TS.s_c = out.C.s_vec;
    TS.x_flux = out.x_flux;
    TS.x_geom = out.x_geom;
else
    TS = NaN;
end
% If the param.displayTS flag is activated (=1), the temperature profile is
% plotted in a new figure
if param.displayTS == 1
    figure
    subplot(2,2,1)
    hold on
    ii = 1;
    L2(ii) = plot(TS.x_flux, TS.T_c-273.15,'s-b' ,'linewidth',2);
    lg2{ii} = out.C.fluid;
    ii=ii+1;
    L2(ii) = plot(TS.x_flux, TS.T_w-273.15,'+:k' ,'linewidth',1);
    lg2{ii} = 'wall';
    ii=ii+1;
    L2(ii) = plot(TS.x_flux, TS.T_h-273.15,'o-r' ,'linewidth',2);
    lg2{ii} = out.H.fluid;
    grid on
    xlabel('Heat fraction [-]','fontsize',14,'fontweight','bold')
    ylabel('Temperature [°C]','fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')    %hold on
    legend(L2,lg2,'interpreter', 'none','location', 'southeast')
    
    subplot(2,2,2)
    hold on
    ii = 1;
    L1(ii) = plot(TS.x_geom, TS.T_c-273.15,'s-b' ,'linewidth',2);
    lg1{ii} = out.C.fluid;
    ii=ii+1;
    L1(ii) = plot(TS.x_geom, TS.T_w-273.15,'+:k' ,'linewidth',1);
    lg1{ii} = 'wall';
    ii=ii+1;
    L1(ii) = plot(TS.x_geom, TS.T_h-273.15,'o-r' ,'linewidth',2);
    lg1{ii} = out.H.fluid;
    grid on
    xlabel('Spatial fraction [-]','fontsize',14,'fontweight','bold')
    ylabel('Temperature [°C]','fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')
    legend(L1,lg1,'interpreter', 'none','location', 'southeast')
    

    subplot(2,2,3)
    hold on
    ii = 1;
    L3(ii) = plot(TS.x_flux(1:end-1), out.C.hConv_vec,'s-b' ,'linewidth',2);
    lg3{ii} = out.C.fluid;
    ii=ii+1;
    L3(ii) = plot(TS.x_flux(1:end-1), out.H.hConv_vec,'o-r' ,'linewidth',2);
    lg3{ii} = out.H.fluid;
    ii=ii+1;
    L3(ii) = plot(TS.x_flux(1:end-1), out.U_vec,'v--k' ,'linewidth',2);
    lg3{ii} = 'global U';
    grid on
    xlabel('Heat fraction [-]','fontsize',14,'fontweight','bold')
    ylabel('HTC [W/m².K]','fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')    %hold on
    legend(L3,lg3,'interpreter', 'none','location', 'southeast')
    
    subplot(2,2,4)
    hold on
    ii = 1;
    if param.C.solub
        L3(ii) = plot(TS.x_flux, out.C.x_vec,'^-.b' ,'linewidth',1); lg3{ii} = [out.C.fluid ' - x']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, out.C.zeta_r_vec,'d-b' ,'linewidth',1); lg3{ii} = [out.C.fluid ' - \zeta']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, out.C.C_rv_vec,'s--b' ,'linewidth',1.5); lg3{ii} = [out.C.fluid ' - C_rv']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, out.C.C_rl_vec,'o-b' ,'linewidth',1.5); lg3{ii} = [out.C.fluid ' - C_rl']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, 1-out.C.C_rl_vec - out.C.C_rv_vec,'v:b' ,'linewidth',1.5); lg3{ii} = [out.C.fluid ' - C_lub']; ii = ii+1;
    end
    if param.H.solub
        L3(ii) = plot(TS.x_flux, out.H.x_vec,'^-.r' ,'linewidth',1); lg3{ii} = [out.H.fluid ' - x']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, out.H.zeta_r_vec,'d-r' ,'linewidth',1); lg3{ii} = [out.H.fluid ' - \zeta']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, out.H.C_rv_vec,'s--r' ,'linewidth',1.5); lg3{ii} = [out.H.fluid ' - C_rv']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, out.H.C_rl_vec,'o-r' ,'linewidth',1.5); lg3{ii} = [out.H.fluid ' - C_rl']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, 1-out.H.C_rl_vec - out.H.C_rv_vec,'v:r' ,'linewidth',1.5); lg3{ii} = [out.H.fluid ' - C_lub']; ii = ii+1;
    end
    grid on
    xlabel('Heat fraction [-]','fontsize',14,'fontweight','bold')
    ylabel('Flow composition','fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')    %hold on
    legend(L3,lg3,'interpreter', 'none','location', 'southeast')
    
end

% If the param.displayResults flag is activated (=1), the results are displayed on the
% command window
if param.displayResults ==1
    in.fluid_h = fluid_h;
    in.m_dot_h = m_dot_h;
    in.in_h_su = in_h_su;
    in.type_h = param.H.type;
    in.P_h_su = P_h_su;
    in.fluid_c = fluid_c;
    in.m_dot_c = m_dot_c;
    in.in_c_su = in_c_su;
    in.type_c = param.C.type;
    in.P_c_su = P_c_su;
    in.modelType= param.modelType;
    
    if nargin ==0
        fprintf ( 1, '\n' );
        disp('-------------------------------------------------------')
        disp('--------------------   Demo Code   --------------------')
        disp('-------------------------------------------------------')
        fprintf ( 1, '\n' );
    end
    disp('Working conditions:')
    fprintf ( 1, '\n' );
    disp(in)
    disp('Global results')
    disp(out)
    disp('Hot-side results')
    disp(out.H)
    disp('Cold-side results')
    disp(out.C)
    
end

end





