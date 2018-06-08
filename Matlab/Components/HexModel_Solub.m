function [out,TS] = HexModel_Solub(fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems
% Remi Dickes - 27/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be


%% DEMONSTRATION CASE
if nargin == 0    % Define a demonstration case if HexModel.mat is not executed externally
    clear all
    close all
    %clc
    
    % Supply fluid conditions
    fluid_h = 'ICP_PiroblocBasic';                                       	% Nature of the hot fluid           [-]
    m_dot_h = 0.1;                                                      	% Mass flow rat of the hot fluid    [kg/s]
    P_h_su =  1e5;                                                      	% Supply pressure of the hot fluid  [Pa]
    in_h_su = PropsSI_ICP('H', 'T', 140+273.15, 'P', P_h_su, fluid_h) ;      % Supply h or T of the hot fluid  	[J/kg pr K]
    param.H.type = 'H';
    param.H.solub = 0;
    
    fluid_c = 'R245fa';                                                     % Nature of the cold fluid        	[-]
    m_dot_c = 0.06;                                                          % Mass flow rat of the cold fluid  	[kg/s]
    T_c_su = 50 + 273.15;
    P_c_su = 9e5;                                                          % Supply pressure of the cold fluid	[Pa]
    param.C.type = 'H';
    param.C.solub = 1;
    param.C.fluid_lub = 'ICP_RL_32_3MAF';
    param.C.C_oil = 0.02;
    h_ref = CoolProp.PropsSI('H','P',P_c_su,'T',T_c_su, fluid_c);      % Supply h or T of the cold fluid  	[J/kg pr K]
    h_lub = PropsSI_ICP('H', 'T', T_c_su, 'P', P_c_su, param.C.fluid_lub);
    in_c_su = (1-param.C.C_oil)*h_ref + param.C.C_oil*h_lub;
    load('Ratio_RhoIdeal_R245Fa_POE.mat')
    param.C.fit_ratio_rho = fit_ratio_rho;
    
    % Geometrical parameters 
    param.theta = 60*pi/180;    
    param.H.A_tot = 9.8;
    param.C.A_tot = 9.8;
    param.H.V_tot = 0.009;
    param.C.V_tot = 0.009;
    param.L_hex = 0.459;
    param.pitch_co = 0.0072;
    param.phi =  1.1407;
    param.C.Dh = 0.0032;
    param.H.Dh = 0.0032;
    param.C.CS = 3.4380e-04;
    param.H.CS = 3.4380e-04;
    param.C.n_canals = 50;
    param.H.n_canals = 49;
    param.H.fin = 'none';    
    param.C.fin = 'none';
    param.k0 = 13.022;
    param.k1 = 0.01667;
    param.C.R_fooling = 0.35/1000;
    param.H.R_fooling = 0.35/1000;
    param.t = 0.4/1000;
    
    param.n_disc = 30;
    param.modelType = 'hConvCor';   
    param.H.correlation.type_1phase_l = 'Martin1_BPHEX';
    param.H.correlation.type_1phase_v = 'Martin1_BPHEX';
    param.H.correlation.type_2phase_cd = 'Han_boiling';
    param.H.correlation.type_2phase_ev = 'Han_boiling';
    param.C.correlation.type_1phase_l = 'Martin1_BPHEX';
    param.C.correlation.type_1phase_v = 'Martin1_BPHEX';
    param.C.correlation.type_2phase_ev = 'Han_boiling';
    param.C.correlation.type_2phase_cd = 'Han_boiling';
    
    param.H.correlation.type_void_fraction = 'Homogenous';
    param.C.correlation.type_void_fraction = 'Homogenous';
    
    param.C.correlation.dry_out_incipient = 'Kim_DOI';
    param.C.x_di = 0.8;
    
    param.hugmark_simplified = 1;
    
    param.displayResults = 1;
    param.displayTS = 1;
    param.generateTS = 1;
    
end

tstart_hex = tic;

%% HEAT EXCHANGER MODELLING
C_rl_lim = 0.02;

% Evaluation of the hot fluid (HF) supply conditions
if strcmp(param.H.type,'H')
    if strcmp(fluid_h(1:3), 'ICP')
        T_h_su = PropsSI_ICP('T', 'H', in_h_su, 'P', P_h_su, fluid_h);
        h_h_l = inf;
        h_h_v = inf;
    else
        T_h_su = CoolProp.PropsSI('T','H', in_h_su, 'P', P_h_su, fluid_h);
        if param.H.solub
            T_sat_pure_h = CoolProp.PropsSI('T', 'P', P_h_su, 'Q', 0, fluid_h);
            T_min_bubble_h = R245fa_POE_Tbubble(1-param.H.C_oil, P_h_su);
            h_h_l = (1-param.H.C_oil)*CoolProp.PropsSI('H', 'T', T_min_bubble_h, 'Q', 0, fluid_h) + param.H.C_oil*PropsSI_ICP('H', 'T', T_min_bubble_h, 'P', P_h_su, param.H.fluid_lub);
            x_dryout_cpt = 1-(C_rl_lim/(1-param.H.C_oil)); %limit quality of vapor-like region
            zeta_dryout_cpt = (1- x_dryout_cpt - param.H.C_oil + x_dryout_cpt*param.H.C_oil)/(1- x_dryout_cpt + x_dryout_cpt*param.H.C_oil);
            T_dryout_cpt = R245fa_POE_Tbubble(zeta_dryout_cpt, P_h_su);
            [h_h_v, ~, ~, ~, ~, ~, ~, ~] = TP_solubMixt(T_dryout_cpt, P_h_su, param.H.C_oil, fluid_h, param.H.fluid_lub, T_min_bubble_h, T_sat_pure_h);
        else
            h_h_l = CoolProp.PropsSI('H','P',P_h_su,'Q',0,fluid_h);
            h_h_v = CoolProp.PropsSI('H','P',P_h_su,'Q',1,fluid_h);
        end
    end
elseif strcmp(param.H.type,'T')
    T_h_su = in_h_su;
    h_h_l = inf;
    h_h_v = inf;
end

% Evaluation of the cold fluid (CF) supply conditions
if strcmp(param.C.type,'H')
    if strcmp(fluid_c(1:3), 'ICP')
        T_c_su = PropsSI_ICP('T', 'H', in_c_su, 'P', P_c_su, fluid_c);
        h_c_l = inf;
        h_c_v = inf;
    else
        T_c_su = CoolProp.PropsSI('T','H', in_c_su, 'P', P_c_su, fluid_c);
        if param.C.solub
            T_sat_pure_c = CoolProp.PropsSI('T', 'P', P_c_su, 'Q', 0, fluid_c);
            T_min_bubble_c = R245fa_POE_Tbubble(1-param.C.C_oil, P_c_su);
            h_c_l = (1-param.C.C_oil)*CoolProp.PropsSI('H', 'T', T_min_bubble_c, 'Q', 0, fluid_c) + param.C.C_oil*PropsSI_ICP('H', 'T', T_min_bubble_c, 'P', P_c_su, param.C.fluid_lub);
            x_dryout_cpt = 1-(C_rl_lim/(1-param.C.C_oil)); %limit quality of vapor-like region
            zeta_dryout_cpt = (1- x_dryout_cpt - param.C.C_oil + x_dryout_cpt*param.C.C_oil)/(1- x_dryout_cpt + x_dryout_cpt*param.C.C_oil);
            T_dryout_cpt = R245fa_POE_Tbubble(zeta_dryout_cpt, P_c_su);
            [h_c_v, ~, ~, ~, ~, ~, ~, ~] = TP_solubMixt(T_dryout_cpt, P_c_su, param.C.C_oil, fluid_c, param.C.fluid_lub, T_min_bubble_c, T_sat_pure_c);
        else
            h_c_l = CoolProp.PropsSI('H','P',P_c_su,'Q',0,fluid_c);
            h_c_v = CoolProp.PropsSI('H','P',P_c_su,'Q',1,fluid_c);
        end
    end
elseif strcmp(param.C.type,'T')
    T_c_su = in_c_su;
    h_c_l = inf;
    h_c_v = inf;
end

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
    [Q_dot_max, pinch_Qdot_max] = HEX_Qdotmax_Solub(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param, h_h_l, h_h_v, h_c_l, h_c_v); % compute maximal heat transfer
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
    [out.H.s_vec, out.C.s_vec] = deal(NaN*ones(1, length(out.H.T_vec)));
    if param.generateTS
       if strcmp(param.H.type,'H') && not(strcmp(fluid_h(1:3), 'ICP'))  %if not an incompressible fluid, calculate entropy vector
           for i = 1: length(out.H.H_vec)
               out.H.s_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H.H_vec(i),fluid_h);
           end
       end
       if strcmp(param.C.type,'H') && not(strcmp(fluid_c(1:3), 'ICP')) %if not an incompressible fluid, calculate entropy vector
           for i = 1: length(out.C.H_vec)
               out.C.s_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.C.H_vec(i),fluid_c);
           end
       end
    end
    
    % Mass calculation
    out.H.V_vec = param.H.V_tot*(out.H.A_vec./sum(out.H.A_vec));
    out.C.V_vec = param.C.V_tot*(out.H.A_vec./sum(out.H.A_vec));
    [out.H.M_vec, out.H.M_mf_vec, out.H.M_lub_vec, out.H.alpha_mean_vec, out.C.M_vec, out.C.M_mf_vec, out.C.M_lub_vec, out.C.alpha_mean_vec] = deal(NaN*ones(1, length(out.H.T_vec)-1));
    
    for i = 1:length(out.H.V_vec) % Hot-side fluid        
        if param.H.solub %if solubility    
            [~, ~, out.H.M_lub_vec(i), out.H.M_mf_vec(i), ~, ~, out.H.M_vec(i), out.H.alpha_mean_vec(i)] = MassMeanInt_solub(out.H.T_vec(i), P_h_su, out.H.zeta_r_vec(i), out.H.C_rv_vec(i), out.H.Tbubble_min_vec(i), out.H.Tsat_pure_vec(i), out.H.T_vec(i+1), P_h_su, out.H.zeta_r_vec(i+1), out.H.C_rv_vec(i+1), out.H.Tbubble_min_vec(i+1), out.H.Tsat_pure_vec(i+1), out.H.V_vec(i), fluid_h, param.H.fluid_lub, param.H.fit_ratio_rho, param.H);
        else % if no solubility
            if strcmp(param.H.type,'T') || strcmp(fluid_h(1:3),'ICP')
                out.H.M_mf_vec(i) = out.H.V_vec(i)*(0.5*PropsSI_ICP('D', 'T', out.H.T_vec(i), 'P', P_h_su, fluid_h) + PropsSI_ICP('D', 'T', out.H.T_vec(i+1), 'P', P_h_su, fluid_h));
                out.H.M_lub_vec(i) = NaN;
                out.H.M_vec(i) = out.H.M_mf_vec(i);
            else % if compressible
                if strcmp(out.H.type_zone{i},  'tp')
                    q1 = CoolProp.PropsSI('Q','P',P_h_su,'H',out.H.H_vec(i),fluid_h);
                    q2 = CoolProp.PropsSI('Q','P',P_h_su,'H',out.H.H_vec(i+1),fluid_h);
                    rho_v = CoolProp.PropsSI('D','P',P_h_su,'Q',1,fluid_h);
                    rho_l = CoolProp.PropsSI('D','P',P_h_su,'Q',0,fluid_h);
                    out.H.alpha_mean_vec(i) = VoidFraction_Integration(q1, q2, rho_v, rho_l, param);
                    out.H.M_mf_vec(i) =  out.H.V_vec(i)*(rho_v*out.H.alpha_mean_vec(i)+ rho_l*(1-out.H.alpha_mean_vec(i)));
                else
                    out.H.M_mf_vec(i) = out.H.V_vec(i)*(0.5*CoolProp.PropsSI('D','P',P_h_su,'H',out.H.H_vec(i),fluid_h) + 0.5*CoolProp.PropsSI('D','P',P_h_su,'H',out.H.H_vec(i+1),fluid_h));
                end
                out.H.M_lub_vec(i) = NaN;
                out.H.M_vec(i) = out.H.M_mf_vec(i);
            end
        end        
    end
    out.H.M_tot = sum(out.H.M_vec);
    out.H.M_lub_tot = sum(out.H.M_lub_vec);
    out.H.M_mf_tot = sum(out.H.M_mf_vec);
    
    for i = 1:length(out.C.V_vec) % Cold-side fluid       
        if param.C.solub %if solubility
            [~, ~, out.C.M_lub_vec(i), out.C.M_mf_vec(i), ~, ~, out.C.M_vec(i), out.C.alpha_mean_vec(i)] = MassMeanInt_solub(out.C.T_vec(i), P_c_su, out.C.zeta_r_vec(i), out.C.C_rv_vec(i), out.C.Tbubble_min_vec(i), out.C.Tsat_pure_vec(i), out.C.T_vec(i+1), P_c_su, out.C.zeta_r_vec(i+1), out.C.C_rv_vec(i+1), out.C.Tbubble_min_vec(i+1), out.C.Tsat_pure_vec(i+1), out.C.V_vec(i), fluid_c, param.C.fluid_lub, param.C.fit_ratio_rho, param.C);
        else % if no solubility
            if strcmp(param.C.type,'T') || strcmp(fluid_c(1:3),'ICP')
                out.C.M_mf_vec(i) = out.C.V_vec(i)*(0.5*PropsSI_ICP('D', 'T', out.C.T_vec(i), 'P', P_c_su, fluid_c) + PropsSI_ICP('D', 'T', out.C.T_vec(i+1), 'P', P_c_su, fluid_c));
                out.C.M_lub_vec(i) = NaN;
                out.C.M_vec(i) = out.C.M_mf_vec(i);
            else % if compressible
                if strcmp(out.C.type_zone{i},  'tp')
                    q1 = CoolProp.PropsSI('Q','P',P_c_su,'H',out.C.H_vec(i),fluid_c);
                    q2 = CoolProp.PropsSI('Q','P',P_c_su,'H',out.C.H_vec(i+1),fluid_c);
                    rho_v = CoolProp.PropsSI('D','P',P_c_su,'Q',1,fluid_c);
                    rho_l = CoolProp.PropsSI('D','P',P_c_su,'Q',0,fluid_c);
                    out.C.alpha_mean_vec(i) = VoidFraction_Integration(q1, q2, rho_v, rho_l, param);
                    out.C.M_mf_vec(i) =  out.C.V_vec(i)*(rho_v*out.C.alpha_mean_vec(i)+ rho_l*(1-out.C.alpha_mean_vec(i)));
                else
                    out.C.M_mf_vec(i) = out.C.V_vec(i)*(0.5*CoolProp.PropsSI('D','P',P_c_su,'H',out.C.H_vec(i),fluid_c) + 0.5*CoolProp.PropsSI('D','P',P_c_su,'H',out.C.H_vec(i+1),fluid_c));
                end
                out.C.M_lub_vec(i) = NaN;
                out.C.M_vec(i) = out.C.M_mf_vec(i);
            end
        end
    end
    
    out.C.M_tot = sum(out.C.M_vec); 
    out.C.M_lub_tot = sum(out.C.M_lub_vec);
    out.C.M_mf_tot = sum(out.C.M_mf_vec);  
    
    % Flag evaluation
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
       
    
else
    if 0
        %If no, there is not any heat power transfered
        Q_dot_eff = 0;
        out = HEX_profile_4(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, param, h_h_l, h_h_v, P_h_crit, h_h_crit, h_c_l, h_c_v, P_c_crit, h_c_crit);
        out.H.h_ex = out.H.H_vec(1);
        out.C.h_ex = out.C.H_vec(end);
        out.H.T_ex = out.H.T_vec(1);
        out.C.T_ex = out.C.T_vec(end);
        out.Q_dot_tot = Q_dot_eff;
        
        % Entropy calculation
        [out.H.s_vec, out.C.s_vec] = deal(NaN*ones(1, length(out.H.H_vec)));
        if param.generateTS
            if strcmp(param.H.type,'H') %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.H.H_vec)
                    out.H.s_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H.H_vec(i),fluid_h);
                end
            end
            if strcmp(param.C.type,'H') %if not an incompressible fluid, calculate entropy vector
                for i = 1: length(out.C.H_vec)
                    out.C.s_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.C.H_vec(i),fluid_c);
                end
            end
        end
        % Mass calculation
        out.H.V_vec = param.H.V_tot*(out.Qdot_vec./out.Q_dot_tot);
        out.C.V_vec = param.C.V_tot*(out.Qdot_vec./out.Q_dot_tot);
        for i = 1: length(out.H.V_vec)
            if strcmp(param.H.type,'H')
                out.H.M_vec(i) = out.H.V_vec(i)*(CoolProp.PropsSI('D','P',P_h_su,'H',out.H.H_vec(i),fluid_h)+CoolProp.PropsSI('D','P',P_h_su,'H',out.H.H_vec(i+1),fluid_h))/2;
            else
                out.H.M_vec(i) = out.H.V_vec(i)*sf_PropsSI_bar('D', out.H.T_vec(i),  out.H.T_vec(i+1), P_h_su, fluid_h);
            end
            if strcmp(param.C.type,'H')
                out.C.M_vec(i) = out.C.V_vec(i)*(CoolProp.PropsSI('D','P',P_c_su,'H',out.C.H_vec(i),fluid_c)+CoolProp.PropsSI('D','P',P_c_su,'H',out.C.H_vec(i+1),fluid_c))/2;
            else
                out.C.M_vec(i) = out.C.V_vec(i)*sf_PropsSI_bar('D', out.C.T_vec(i),  out.C.T_vec(i+1), P_c_su, fluid_c);
            end
        end
        out.H.M_tot = sum(out.H.M_vec);
        out.C.M_tot = sum(out.C.M_vec);
        
        if (T_h_su-T_c_su)<1e-2  && (T_h_su-T_c_su >0)  && m_dot_h  > 0 && m_dot_c > 0
            out.flag = 3;
        else
            out.flag = -3;
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
    TS.x_geom = [0 cumsum(out.H.V_vec)./param.H.V_tot];
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
        L3(ii) = plot(TS.x_flux, out.H.C_rv_vec,'s--r' ,'linewidth',2); lg3{ii} = [out.H.fluid ' - C_rv']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, out.H.C_rl_vec,'o-r' ,'linewidth',2); lg3{ii} = [out.H.fluid ' - C_rl']; ii = ii+1;
        L3(ii) = plot(TS.x_flux, 1-out.H.C_rl_vec - out.H.C_rv_vec,'v:r' ,'linewidth',2); lg3{ii} = [out.H.fluid ' - C_lub'];
    end
    grid on
    xlabel('Heat fraction [-]','fontsize',14,'fontweight','bold')
    ylabel('HTC [W/m².K]','fontsize',14,'fontweight','bold')
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





