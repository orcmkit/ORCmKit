function [out, TS] = FCT_ORC_Ext_Npp_Nexp_3_Solub(x, lb, ub, fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param)

% Impose :   All external conditions, N_exp, N_pp, either DT_sc or M_tot
% Guess :    P_pp_ex, P_pp_su, h_ev_ex (and h_pp_su if M_tot is imposed)
% Residuals: Q_dot_rec, N_exp, h_cd_ex (and M_tot, if charge imposed)    
%tic

x = max(x, lb./ub);
x = min(x,ones(1,length(x)));
x = x.*ub;
x(1) = max(x(1), 1.1*x(2));



% --------------------------------------------------------------------------------
% RECEIVER SUPPLY ----------------------------------------------------------------
% --------------------------------------------------------------------------------
out.DT_sc_bis = param.DT_sc;
out.P_lr_su = x(2); 
Tsat_pure_lr_su = CoolProp.PropsSI('T', 'Q', 0, 'P', out.P_lr_su, fluid_wf); 
if C_oil > 0
    Tbubble_min_lr_su = R245fa_POE_Tbubble(1-C_oil, out.P_lr_su, Tsat_pure_lr_su);
    out.T_lr_su = Tbubble_min_lr_su - param.DT_sc;
    [out.h_lr_su,  ~, ~,  ~, out.zeta_lr_su, ~, ~, ~] = TP_solubMixt(out.T_lr_su, out.P_lr_su, C_oil, fluid_wf, fluid_lub, Tbubble_min_lr_su, Tsat_pure_lr_su, param.fit_DTP_zeta);
    out.DT_sc = param.DT_sc;
    [rho_vap_lr_su, ~, ~, rho_liq_lr_su, ~] = R245fa_POE_density(out.T_lr_su , out.P_lr_su, out.zeta_lr_su, fluid_wf, fluid_lub, param.fit_ratio_rho, Tbubble_min_lr_su, Tsat_pure_lr_su);
else
    out.T_lr_su = Tsat_pure_lr_su - param.DT_sc;    
    if param.DT_sc == 0
        out.h_lr_su = CoolProp.PropsSI('H', 'P', out.P_lr_su, 'Q', 0, fluid_wf);
        out.DT_sc = param.DT_sc;
        rho_vap_lr_su = CoolProp.PropsSI('D', 'P', out.P_lr_su, 'Q', 1, fluid_wf);
        rho_liq_lr_su = CoolProp.PropsSI('D', 'P', out.P_lr_su, 'Q', 0, fluid_wf);
    elseif param.DT_sc > 0
        out.T_lr_su = Tsat_pure_lr_su - param.DT_sc;
        out.h_lr_su = CoolProp.PropsSI('H', 'P', out.P_lr_su, 'T', out.T_lr_su, fluid_wf);
        out.DT_sc = param.DT_sc;
        rho_vap_lr_su = 0;
        rho_liq_lr_su = CoolProp.PropsSI('D', 'P', out.P_lr_su, 'H', out.h_lr_su, fluid_wf);
    elseif param.DT_sc < 0
        if param.DT_sc >=-1
            out.T_lr_su = Tsat_pure_lr_su;
            out.h_lr_su = CoolProp.PropsSI('H', 'P', out.P_lr_su, 'Q', abs(param.DT_sc), fluid_wf);
            out.DT_sc = 0;
            rho_vap_lr_su = CoolProp.PropsSI('D', 'P', out.P_lr_su, 'Q', 1, fluid_wf);
            rho_liq_lr_su = CoolProp.PropsSI('D', 'P', out.P_lr_su, 'H', out.h_lr_su, fluid_wf);
        elseif param.DT_sc < -1
            out.T_lr_su = CoolProp.PropsSI('T', 'P', out.P_lr_su, 'Q', 0, fluid_wf) - (param.DT_sc+1);
            out.h_lr_su = CoolProp.PropsSI('H', 'P', out.P_lr_su, 'T', out.T_lr_su, fluid_wf);
            out.DT_sc = param.DT_sc+1;
            rho_vap_lr_su = CoolProp.PropsSI('D', 'P', out.P_lr_su, 'H', out.h_lr_su, fluid_wf);
            rho_liq_lr_su = rho_vap_lr_su;
        end
    end

end

% --------------------------------------------------------------------------------
% PUMP ---------------------------------------------------------------------------
% --------------------------------------------------------------------------------
out.P_pp_ex = x(1);
out.P_pp_su = out.P_lr_su + param.PIPE.H_lr_pp*9.81*rho_liq_lr_su;
out.rp_pp = out.P_pp_ex/out.P_pp_su;
out.h_pp_su = out.h_lr_su - param.PIPE.AU_lr_pp*(out.T_lr_su-T_amb);
[out_PP, ~] = PumpModel3_Solub(out.P_pp_su, out.h_pp_su, out.P_pp_ex, fluid_wf,fluid_lub, C_oil, N_pp, T_amb, param.PP);
out.T_pp_su = out_PP.T_su;
out.m_dot_wf = out_PP.m_dot;
out.W_dot_pp = out_PP.W_dot;
out.T_pp_ex = out_PP.T_ex;
out.h_pp_ex = out_PP.h_ex;
out.eps_is_pp = out_PP.epsilon_is;
out.eps_vol_pp = out_PP.epsilon_vol;
out.time_pp = out_PP.time;
out.flag_pp = out_PP.flag;
i_flag = 1;
out.flag.value(1,i_flag) = out_PP.flag;
out.flag.name{1,i_flag} = 'flag_pp';

out.P_recc_su = out.P_pp_ex - param.PIPE.f_dp_fltd(out.m_dot_wf); 
out.h_recc_su = out.h_pp_ex - param.PIPE.AU_pp_recc*(out.T_pp_ex-T_amb); 


% --------------------------------------------------------------------------------
% EXPANDER -----------------------------------------------------------------------
% --------------------------------------------------------------------------------
out.h_exp_su = x(3); % ok
out.P_exp_ex =  out.P_lr_su + param.CD.H.f_dp(out.m_dot_wf) + param.REC.H.f_dp(out.m_dot_wf); % ok
out.P_exp_su =  out.P_pp_ex - param.PIPE.f_dp_fltd(out.m_dot_wf)- param.EV.C.f_dp(out.m_dot_wf); % ok
out.rp_exp = out.P_exp_su/out.P_exp_ex; % ok
[out_EXP, ~] = ExpanderModel3_Solub(fluid_wf, fluid_lub, C_oil, out.P_exp_su, out.h_exp_su, out.m_dot_wf, out.P_exp_ex, T_amb, param.EXP); % ok
out.T_exp_su = out_EXP.T_su; % ok
out.N_exp_bis = out_EXP.N_exp; % ok
out.W_dot_exp = out_EXP.W_dot; % ok
out.Q_dot_exp = out_EXP.Q_dot_amb; % ok
out.eps_vol_exp = out_EXP.FF; % ok
out.eps_is_exp = out_EXP.epsilon_is; % ok
out.flag_exp = out_EXP.flag; % ok
out.time_exp = out_EXP.time; % ok
i_flag = i_flag+1; % ok
out.flag.value(1,i_flag) = out_EXP.flag; % ok 
out.flag.name{1,i_flag} = 'flag_exp'; % ok

out.T_rech_su = out_EXP.T_ex;  % ok
out.h_rech_su = out_EXP.h_ex;  % ok
out.P_rech_su = out.P_exp_ex;  % ok


% --------------------------------------------------------------------------------
% RECUPERATOR --------------------------------------------------------------------
% --------------------------------------------------------------------------------
param.REC.H.C_oil = C_oil;
param.REC.C.C_oil = C_oil;
[out_REC,~] = HexModel_Solub_DP2(fluid_wf, out.P_rech_su, out.h_rech_su, out.m_dot_wf, fluid_wf, out.P_recc_su, out.h_recc_su, out.m_dot_wf, param.REC); % ok
out.Q_dot_rec = out_REC.Q_dot_tot; % ok
out.pinch_rec = out_REC.pinch; % ok
out.T_recc_ex = out_REC.C.T_ex;
out.P_recc_ex = out_REC.C.P_ex;
out.h_recc_ex = out_REC.C.h_ex;
out.T_rech_ex = out_REC.H.T_ex;
out.P_rech_ex = out_REC.H.P_ex;
out.h_rech_ex = out_REC.H.h_ex;
out.flag_rec = out_REC.flag; % ok
out.time_rec = out_REC.time; % ok
i_flag = i_flag+1; % ok
out.flag.value(1,i_flag) = out_REC.flag; % ok
out.flag.name{1,i_flag} = 'flag_rec'; % ok

out.h_cd_su = out.h_rech_ex; % ok
out.P_cd_su = out.P_rech_ex; % ok
out.T_cd_su = out.T_rech_ex; % ok    
out.h_ev_su = out.h_recc_ex; % ok
out.P_ev_su = out.P_recc_ex; % ok
out.T_ev_su = out.T_recc_ex; % ok

% --------------------------------------------------------------------------------
% EVAPORATOR ---------------------------------------------------------------------
% --------------------------------------------------------------------------------
param.EV.C.C_oil = C_oil;
[out_EV, ~] = HexModel_Solub_DP2(fluid_htf, P_htf_su, h_htf_su, m_dot_htf, fluid_wf, out.P_ev_su, out.h_ev_su, out.m_dot_wf , param.EV);
out.Q_dot_ev = out_EV.Q_dot_tot;
out.T_htf_ev_ex = out_EV.H.T_ex;
out.h_htf_ev_ex = out_EV.H.h_ex;
out.pinch_ev = out_EV.pinch;
out.flag_ev = out_EV.flag;
out.time_ev = out_EV.time;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_EV.flag;
out.flag.name{1,i_flag} = 'flag_ev';

out.h_ev_ex = out_EV.C.h_ex;
out.T_ev_ex = out_EV.C.T_ex;
out.P_ev_ex = out_EV.C.P_ex;


% --------------------------------------------------------------------------------
% CONDENSER ----------------------------------------------------------------------
% --------------------------------------------------------------------------------
param.CD.H.C_oil = C_oil;
[out_CD, ~] = HexModel_Solub_DP2(fluid_wf, out.P_cd_su, out.h_cd_su, out.m_dot_wf, fluid_ctf, P_ctf_su, h_ctf_su, m_dot_ctf , param.CD);
out.Q_dot_cd = out_CD.Q_dot_tot;
out.T_ctf_cd_ex = out_CD.C.T_ex;
out.h_ctf_cd_ex = out_CD.C.h_ex;
out.pinch_cd = out_CD.pinch; 
out.flag_cd = out_CD.flag;
out.time_cd = out_CD.time;
out.W_dot_cd = param.CD.W_dot_aux;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_CD.flag;
out.flag.name{1,i_flag} = 'flag_cd';

out.P_cd_ex = out_CD.H.P_ex;
out.T_cd_ex = out_CD.H.T_ex;
out.h_cd_ex = out_CD.H.h_ex;

% --------------------------------------------------------------------------------
% MASS INVENTORY -----------------------------------------------------------------
% --------------------------------------------------------------------------------
i_mass = 1;
out.Mass.name{1,i_mass} = 'EV';
out.Mass.value_lub(1,i_mass) = out_EV.C.M_lub_tot; 
out.Mass.value_wf(1,i_mass)  = out_EV.C.M_mf_tot;  

i_mass = i_mass + 1;
out.Mass.name{1,i_mass} = 'CD';
out.Mass.value_lub(1,i_mass) = out_CD.H.M_lub_tot; 
out.Mass.value_wf(1,i_mass)  = out_CD.H.M_mf_tot;  

i_mass = i_mass + 1;
out.Mass.name{1,i_mass} = 'RECC';
out.Mass.value_lub(1,i_mass) = out_REC.C.M_lub_tot; 
out.Mass.value_wf(1,i_mass)  = out_REC.C.M_mf_tot;  

i_mass = i_mass + 1;
out.Mass.name{1,i_mass} = 'RECH';
out.Mass.value_lub(1,i_mass) = out_REC.H.M_lub_tot; 
out.Mass.value_wf(1,i_mass)  = out_REC.H.M_mf_tot;  

Zone_mass_pipe =    {'PP_EX_AUX',               'RECC_EX_AUX',              'EV_EX_AUX',                'EXP_EX_AUX',               'RECH_EX_AUX',                  'CD_EX_AUX',                'LR_EX_AUX'                 'PP1',             	'PP2',       	'EXP1',         'EXP2'};
h_mass_pipe =       [out.h_pp_ex,               out.h_recc_ex,              out.h_ev_ex,                out.h_rech_su,              out.h_rech_ex,                  out.h_cd_ex,                out.h_pp_su,                out.h_pp_su,        out.h_pp_ex,    out.h_exp_su,	out.h_rech_su];
T_mass_pipe =       [out.T_pp_ex,               out.T_recc_ex,              out.T_ev_ex,                out.T_rech_su,              out.T_rech_ex,                  out.T_cd_ex,                out.T_pp_su,                out.T_pp_su,        out.T_pp_ex,    out.T_exp_su,	out.T_rech_su];
P_mass_pipe =       [out.P_pp_ex,               out.P_recc_ex,              out.P_ev_ex,                out.P_rech_su,              out.P_rech_ex,               	out.P_cd_ex,                out.P_pp_su,                out.P_pp_su,        out.P_pp_ex,    out.P_exp_su, 	out.P_rech_su];
V_mass_pipe =       [param.PIPE.V_pp_ex_aux     param.PIPE.V_recc_ex_aux 	param.PIPE.V_ev_ex_aux      param.PIPE.V_exp_ex_aux  	param.PIPE.V_rech_ex_aux        param.PIPE.V_cd_ex_aux      param.PIPE.V_lr_ex_aux,     param.PP.V/2,       param.PP.V/2, 	param.EXP.V/2,	param.EXP.V/2];

if C_oil > 0
    for kk = 1:11
        i_mass = i_mass + 1;
        out.Mass.name{1,i_mass} = Zone_mass_pipe{kk};
        if kk<=7
%             eval(['Zone = param.PIPE.geom.' Zone_mass_pipe{kk} ';'])
%             out.Mass.value_wf(1,i_mass) = 0;
%             out.Mass.value_lub(1,i_mass) = 0;
%             for j = 1:length(Zone.D_vec)
%                 param_local.type_void_fraction = 'AnnularFlow';
%                 param_local.G = out.m_dot_wf/(pi*( Zone.D_vec(j)/2)^2);
%                 param_local.P_Pa = P_mass_pipe(kk);
%                 param_local.T_K = T_mass_pipe(kk);
%                 param_local.C_oil = C_oil;
%                 param_local.D = Zone.D_vec(j);
%                 param_local.fluid_ref = fluid_wf;
%                 param_local.fluid_lub = fluid_lub;
%                 param_local.fit_DTP_zeta = [];
%                 param_local.fit_ratio_rho = param.fit_ratio_rho;
%                 [M_local_wf, M_local_lub] = Mass_function_cycle_solub(T_mass_pipe(kk), P_mass_pipe(kk) , Zone.V_vec(j), C_oil, fluid_wf, fluid_lub, param_local);
%                 out.Mass.value_wf(1,i_mass)  = out.Mass.value_wf(1,i_mass)  + M_local_wf;
%                 out.Mass.value_lub(1,i_mass) = out.Mass.value_lub(1,i_mass) + M_local_lub;
%             end
            param_local.type_void_fraction = 'SlipRatio';
            param_local.S_ratio = 60;
            param_local.fit_DTP_zeta = [];
            param_local.fit_ratio_rho = param.fit_ratio_rho;
            [out.Mass.value_wf(1,i_mass), out.Mass.value_lub(1,i_mass)] = Mass_function_cycle_solub(T_mass_pipe(kk), P_mass_pipe(kk) , V_mass_pipe(kk), C_oil, fluid_wf, fluid_lub, param_local);
        else
            param_local.type_void_fraction = 'SlipRatio';
            param_local.S_ratio = 1;
            param_local.fit_DTP_zeta = [];
            param_local.fit_ratio_rho = param.fit_ratio_rho;
            [out.Mass.value_wf(1,i_mass), out.Mass.value_lub(1,i_mass)] = Mass_function_cycle_solub(T_mass_pipe(kk), P_mass_pipe(kk) , V_mass_pipe(kk), C_oil, fluid_wf, fluid_lub, param_local);
        end
        
    end
    
else
    param_local.type_void_fraction = 'SlipRatio';
    param_local.S_ratio = 5;
    for kk = 1:11
        i_mass = i_mass + 1;
        out.Mass.name{1,i_mass} = Zone_mass_pipe{kk};        
        [out.Mass.value_wf(1,i_mass), out.Mass.value_lub(1,i_mass)] = Mass_function_cycle_solub(h_mass_pipe(kk), P_mass_pipe(kk) , V_mass_pipe(kk), C_oil, fluid_wf, fluid_lub, param_local);
    end
end
     
% --------------------------------------------------------------------------------
% LIQUID RECEIVER -----------------------------------------------------------------
% --------------------------------------------------------------------------------
i_mass = i_mass + 1;
out.Mass.name{1,i_mass} = 'LR';
if param.DT_sc > 0   
    out.Mass.value_lub(1,i_mass) = param.LR.V*rho_liq_lr_su*C_oil;
    out.Mass.value_wf(1,i_mass) = param.LR.V*rho_liq_lr_su*(1-C_oil);
    out.level_receiver = 1;
    out.Receiver_type = 1;
elseif param.DT_sc == 0
    Mass_wf_max = param.LR.V*rho_liq_lr_su*(1-C_oil);
    Mass_wf_min = param.LR.V*rho_vap_lr_su;
    out.Mass.value_wf(1,i_mass) = max(Mass_wf_min,min(param.M_wf-sum(out.Mass.value_wf),Mass_wf_max));
    out.Mass.value_lub(1,i_mass) = C_oil/(1-C_oil)*out.Mass.value_wf(1,i_mass);
    out.Receiver_type = 2;
    out.level_receiver = (out.Mass.value_wf(1,i_mass)-Mass_wf_min)/(Mass_wf_max-Mass_wf_min);
else
    out.Mass.value_wf(1,i_mass) = param.LR.V*rho_vap_lr_su;
    out.Mass.value_lub(1,i_mass) = 0;
    out.Receiver_type = 3;
    out.level_receiver = 0;
end

% --------------------------------------------------------------------------------
% GLOBAL CHARGE ------------------------------------------------------------------
% --------------------------------------------------------------------------------
if C_oil >0
delta_Mass_lub = param.M_lub - sum(out.Mass.value_lub); %ok
i_mass_ev = 1;
out.Mass.value_lub(1,i_mass_ev) = out.Mass.value_lub(1,i_mass_ev) + delta_Mass_lub; %ok
out.Mass.value_tot = out.Mass.value_wf + out.Mass.value_lub;
out.M_tot_wf = sum(out.Mass.value_wf); % ok
out.M_tot_lub = sum(out.Mass.value_lub);  % ok
out.M_tot = out.M_tot_wf + out.M_tot_lub;  % ok
else
    out.Mass.value_tot = out.Mass.value_wf;
    out.M_tot_wf = sum(out.Mass.value_wf);
    out.M_tot_lub = 0;
    out.M_tot = out.M_tot_wf;
end
% --------------------------------------------------------------------------------
% ORC PERFORMANCE ----------------------------------------------------------------
% --------------------------------------------------------------------------------
out.W_dot_net = out.W_dot_exp - out.W_dot_pp - out.W_dot_cd;  % ok
out.eff_ORC_gross = out.W_dot_exp/out.Q_dot_ev;  % ok
out.eff_ORC_net = out.W_dot_net/out.Q_dot_ev;  % ok

% --------------------------------------------------------------------------------
% RESIDUALS ----------------------------------------------------------------------
% --------------------------------------------------------------------------------
out.res_ORC_H_cd_ex = (1 - out.h_lr_su/out.h_cd_ex); % ok
out.res_ORC_H_ev_ex = (1 - out.h_exp_su/out.h_ev_ex); %ok
out.res_ORC_N_exp = 1-out.N_exp_bis/N_exp; % ok
out.res_ORC_M = (1 - out.M_tot_wf/param.M_wf); % ok


out.res_vec  = [out.res_ORC_N_exp    out.res_ORC_H_cd_ex     out.res_ORC_H_ev_ex]; % ok
out.res  = norm(out.res_vec); % ok
out.x = x; % ok

if param.display_list
    fprintf('%-15s %-10s %-5s %-50s %-15s %-50s %-60s %-15s %-10s %-100s\n', ' ', ' ', ' ', ['[' num2str(x,'%15.4e') ']'], ' ', ' ', [ num2str(out.res, '%.4g') '  [ ' num2str(out.res_vec,'%15.4e') ' ] '], num2str(out.res_ORC_M,'%15.5e'), ' ', num2str(out.flag.value));
end

if param.eval_complete
    out.PP = out_PP;
    out.EXP = out_EXP;
    out.REC = out_REC;
    out.EV = out_EV;
    out.CD = out_CD;
end
%if strcmp(param.eval_type, 'fast')
    TS = NaN;
% else
%     out = orderfields(out);  
%     TS.EXP = TS_EXP;
%     TS.PP = TS_PP;
%     TS.REC = TS_REC;
%     TS.CD = TS_CD;
%     TS.EV = TS_EV;
%     TS.LossesHP = TS_LossesHP;
%     TS.LossesLP = TS_LossesLP;
%     
%     % Construction of the TS variable
%     TS.cycle.s = TS.PP.s;
%     TS.cycle.T = TS.PP.T;
%     if isfield(param, 'REC')
%         TS.cycle.s = [TS.cycle.s  TS.REC.s_c];
%         TS.cycle.T = [TS.cycle.T  TS.REC.T_c];
%     end
%     if isfield(param, 'PRE')
%         TS.cycle.s = [TS.cycle.s  TS.PRE.s_c];
%         TS.cycle.T = [TS.cycle.T  TS.PRE.T_c];
%     end
%     TS.cycle.s = [TS.cycle.s  TS.EV.s_c];
%     TS.cycle.T = [TS.cycle.T  TS.EV.T_c];
%     if isfield(param, 'LossesHP')
%         TS.cycle.s = [TS.cycle.s TS.LossesHP.s];
%         TS.cycle.T = [TS.cycle.T TS.LossesHP.T];
%     end
%     TS.cycle.s = [TS.cycle.s  TS.EXP.s];
%     TS.cycle.T = [TS.cycle.T  TS.EXP.T];
%     if isfield(param, 'REC')
%         TS.cycle.s = [TS.cycle.s  fliplr(TS.REC.s_h)];
%         TS.cycle.T = [TS.cycle.T  fliplr(TS.REC.T_h)];
%     end
%     TS.cycle.s = [TS.cycle.s fliplr(TS.CD.s_h)];
%     TS.cycle.T = [TS.cycle.T fliplr(TS.CD.T_h)];
%     if isfield(param, 'SUB')
%         TS.cycle.s = [TS.cycle.s  fliplr(TS.SUB.s_h)];
%         TS.cycle.T = [TS.cycle.T  fliplr(TS.SUB.T_h)];
%     end
%     if isfield(param, 'LossesLP')
%         TS.cycle.s = [TS.cycle.s TS.LossesLP.s];
%         TS.cycle.T = [TS.cycle.T TS.LossesLP.T];
%     end
%end
%toc
end
