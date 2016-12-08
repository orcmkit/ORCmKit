function [out, TS] = FCT_ORC_Ext_Npp_Nexp_2(x, lb, ub, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param)

% Impose :   All external conditions, N_exp, N_pp, either DT_sc or M_tot
% Guess :    P_pp_ex, P_pp_su, h_ev_ex (and h_pp_su if M_tot is imposed)
% Residuals: Q_dot_rec, N_exp, h_cd_ex (and M_tot, if charge imposed)    


x = max(x, lb./ub);
x = min(x,ones(1,length(x)));
x = x.*ub;
x(1) = max(x(1), 1.001*x(2));
out.rp_pp = x(1)/x(2);
i_flag = 1;
i_mass = 1;

% --------------------------------------------------------------------------------
% PUMP ---------------------------------------------------------------------------
% --------------------------------------------------------------------------------
out.P_pp_su = x(2);
if strcmp(param.solverType, 'M_imposed')
    %     out.h_pp_su = min(x(4), CoolProp.PropsSI('H', 'P', out.P_pp_su*0.999, 'Q', 0, fluid_wf));
    %     out.T_pp_su = CoolProp.PropsSI('T', 'P', out.P_pp_su, 'H', out.h_pp_su, fluid_wf) ;
    out.T_pp_su = min(x(4), CoolProp.PropsSI('T', 'P', out.P_pp_su, 'Q', 0, fluid_wf)-0.0001);
    out.DT_sc = CoolProp.PropsSI('T', 'P', out.P_pp_su, 'Q', 0, fluid_wf)-out.T_pp_su;
    %if x(4) ~= 0
        out.h_pp_su = CoolProp.PropsSI('H', 'P', out.P_pp_su, 'T', out.T_pp_su, fluid_wf);
    %else
        %out.h_pp_su = CoolProp.PropsSI('H', 'P', out.P_pp_su, 'Q', 0, fluid_wf);
    %end
    out.DT_sc = CoolProp.PropsSI('T', 'P', out.P_pp_su, 'Q', 0, fluid_wf)-out.T_pp_su;
elseif strcmp(param.solverType, 'DTsc_imposed')
    out.T_pp_su = CoolProp.PropsSI('T', 'P', out.P_pp_su, 'Q', 0, fluid_wf) - param.DT_sc;
    if param.DT_sc ~= 0
        out.h_pp_su = CoolProp.PropsSI('H', 'P', out.P_pp_su, 'T', CoolProp.PropsSI('T', 'P', out.P_pp_su, 'Q', 0, fluid_wf) - param.DT_sc, fluid_wf);
    else
        out.h_pp_su = CoolProp.PropsSI('H', 'P', out.P_pp_su, 'Q', 0, fluid_wf);
    end
    out.DT_sc = param.DT_sc;
end
[out_PP, TS_PP] = PumpModel3(out.P_pp_su, out.h_pp_su, x(1), fluid_wf, N_pp, T_amb, param.PP);
out.M_pp = out_PP.M;
out.m_dot_wf = out_PP.m_dot;
out.W_dot_pp = out_PP.W_dot;
out.eps_is_pp = out_PP.epsilon_is;
out.eps_vol_pp = out_PP.epsilon_vol;
out.time_pp = out_PP.time;
out.flag_pp = out_PP.flag;
out.flag.value(1,i_flag) = out_PP.flag;
out.flag.name{1,i_flag} = 'flag_pp';
out.Mass.value(1,i_mass) = out_PP.M;
out.Mass.name{1,i_mass} = 'M_pp';
out.P_recc_su = x(1);
out.T_recc_su = out_PP.T_ex;
out.h_recc_su = out_PP.h_ex;

i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_pp_ex*CoolProp.PropsSI('D','H',out.h_recc_su,'P',out.P_recc_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_pp_ex';

% --------------------------------------------------------------------------------
% LossesHP -----------------------------------------------------------------------
% --------------------------------------------------------------------------------
out.h_dphp_su = max(CoolProp.PropsSI('H', 'P', out.P_recc_su, 'Q', 0, fluid_wf),min(x(3), CoolProp.PropsSI('H', 'P', out.P_recc_su, 'T', T_htf_su, fluid_wf)));
out.P_dphp_su = out.P_recc_su;
param.LossesHP.type_in = 'su';
[out_LossesHP, TS_LossesHP] = LossesModel(fluid_wf, out.P_dphp_su, out.h_dphp_su, out.m_dot_wf, T_amb, param.LossesHP);
out.dphp = out_LossesHP.dp;
out.Q_dot_hp = out_LossesHP.Q_dot;
out.flag_lossesHP = out_LossesHP.flag;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_LossesHP.flag;
out.flag.name{1,i_flag} = 'flag_dphp';
out.T_exp_su = out_LossesHP.T_ex;
out.h_exp_su = out_LossesHP.h_ex;
out.P_exp_su = out_LossesHP.P_ex;

i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_ev_ex*CoolProp.PropsSI('D','H',out.h_exp_su,'P',out.P_exp_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_ev_ex';


% --------------------------------------------------------------------------------
% LossesLP -----------------------------------------------------------------------
% --------------------------------------------------------------------------------
param.LossesLP.type_in = 'ex';
[out_LossesLP, TS_LossesLP] = LossesModel(fluid_wf, out.P_pp_su, out.h_pp_su, out.m_dot_wf, T_amb, param.LossesLP);
out.dplp = out_LossesLP.dp;
out.Q_dot_lp = out_LossesLP.Q_dot;
out.flag_lossesLP = out_LossesLP.flag;
out.P_dplp_su = out_LossesLP.P_su;
out.T_dplp_su = out_LossesLP.T_su;
out.h_dplp_su = out_LossesLP.h_su;
out.P_rech_su = out.P_dplp_su;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_LossesLP.flag;
out.flag.name{1,i_flag} = 'flag_dplp';

% --------------------------------------------------------------------------------
% EXPANDER -----------------------------------------------------------------------
% --------------------------------------------------------------------------------
out.rp_exp = out.P_exp_su/out.P_rech_su;
% [out_EXP, TS_EXP] = ExpModel_2_SemiEmp(fluid_wf, out.P_exp_su, out.h_exp_su, out.m_dot_wf, out.P_rech_su, T_amb, param.EXP);
[out_EXP, TS_EXP] = ExpanderModel2(fluid_wf, out.P_exp_su, out.h_exp_su, out.m_dot_wf, out.P_rech_su, T_amb, param.EXP);
out.N_exp_bis = out_EXP.N_exp;
out.W_dot_exp = out_EXP.W_dot;
out.Q_dot_exp = out_EXP.Q_dot_amb;
out.M_exp = out_EXP.M;
out.eps_vol_exp = out_EXP.FF;
out.eps_is_exp = out_EXP.epsilon_is;
out.flag_exp = out_EXP.flag;
out.time_exp = out_EXP.time;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_EXP.flag;
out.flag.name{1,i_flag} = 'flag_exp';
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = out.M_exp;
out.Mass.name{1,i_mass} = 'M_exp';
out.T_rech_su = out_EXP.T_ex;
out.h_rech_su = out_EXP.h_ex;
out.P_rech_su = out.P_rech_su;

i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_exp_ex*CoolProp.PropsSI('D','H',out.h_rech_su,'P',out.P_rech_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_exp_ex';

% --------------------------------------------------------------------------------
% RECUPERATOR --------------------------------------------------------------------
% --------------------------------------------------------------------------------
param.REC.port_h = 'su';
param.REC.port_c = 'su';
[out_REC, TS_REC] = HexModel(fluid_wf, out.P_rech_su, out.h_rech_su, out.m_dot_wf, fluid_wf, out.P_recc_su, out.h_recc_su, out.m_dot_wf, param.REC);
out.Q_dot_rec = out_REC.Q_dot_tot;
out.M_rech = out_REC.M_h;
out.M_recc = out_REC.M_c;
out.flag_rec = out_REC.flag;
out.time_rec = out_EXP.time;
out.pinch_rec = out_REC.pinch;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_REC.flag;
out.flag.name{1,i_flag} = 'flag_rec';
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = out.M_recc;
out.Mass.name{1,i_mass} = 'M_recc';
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = out.M_rech;
out.Mass.name{1,i_mass} = 'M_rech';
out.h_cd_su = out_REC.h_h_ex;
out.P_cd_su = out.P_rech_su;
out.T_cd_su = out_REC.T_h_ex;    
out.h_ev_su = out_REC.h_c_ex;
out.P_ev_su = out.P_recc_su;
out.T_ev_su = out_REC.T_c_ex;

i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_rech_ex*CoolProp.PropsSI('D','H',out.h_cd_su,'P',out.P_cd_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_rech_ex';
i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_recc_ex*CoolProp.PropsSI('D','H',out.h_ev_su,'P',out.P_ev_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_recc_ex';

% --------------------------------------------------------------------------------
% EVAPORATOR ---------------------------------------------------------------------
% --------------------------------------------------------------------------------
param.EV.port_h = 'su';
param.EV.port_c = 'su';
[out_EV, TS_EV] = HexModel(fluid_htf, P_htf_su, in_htf_su, m_dot_htf, fluid_wf, out.P_ev_su, out.h_ev_su, out.m_dot_wf , param.EV);
out.Q_dot_ev = out_EV.Q_dot_tot;
out.T_htf_ev_ex = out_EV.T_h_ex;
out.h_htf_ev_ex = out_EV.h_h_ex;
out.M_ev = out_EV.M_c;
out.flag_ev = out_EV.flag;
out.time_ev = out_EV.time;
out.pinch_ev = out_EV.pinch;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_EV.flag;
out.flag.name{1,i_flag} = 'flag_ev';
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = out.M_ev;
out.Mass.name{1,i_mass} = 'M_ev';
out.h_ev_ex = out_EV.h_c_ex;
out.T_ev_ex = out_EV.T_c_ex;
out.P_ev_ex = out.P_ev_su;


% --------------------------------------------------------------------------------
% CONDENSER ----------------------------------------------------------------------
% --------------------------------------------------------------------------------
param.CD.port_h = 'su';
param.CD.port_c = 'su';
[out_CD, TS_CD] = HexModel(fluid_wf, out.P_cd_su, out.h_cd_su, out.m_dot_wf, fluid_ctf, P_ctf_su, in_ctf_su, m_dot_ctf , param.CD);
out.Q_dot_cd = out_CD.Q_dot_tot;
out.T_ctf_cd_ex = out_CD.T_c_ex;
out.h_ctf_cd_ex = out_CD.h_c_ex;
out.M_cd = out_CD.M_h;
out.flag_cd = out_CD.flag;
out.time_cd = out_CD.time;
out.pinch_cd = out_CD.pinch; 
out.W_dot_cd = param.CD.W_dot_aux;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_CD.flag;
out.flag.name{1,i_flag} = 'flag_cd';
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = out.M_cd;
out.Mass.name{1,i_mass} = 'M_cd';
out.P_cd_ex = out.P_cd_su;
out.T_cd_ex = out_CD.T_h_ex;
out.h_cd_ex = out_CD.h_h_ex;

i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_cd_ex*CoolProp.PropsSI('D','H',out.h_cd_ex,'P',out.P_cd_ex,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_cd_ex';

out.M_cd_vec = out_CD.M_h_vec;
out.V_cd_vec = out_CD.V_h_vec;

% --------------------------------------------------------------------------------
% LIQUID RECEIVER-----------------------------------------------------------------
% --------------------------------------------------------------------------------
if out.h_pp_su<CoolProp.PropsSI('H','Q',0,'P',out.P_pp_su,fluid_wf)
    i_mass = i_mass + 1;
    out.Mass.value(1,i_mass) = param.V_liq_rec*CoolProp.PropsSI('D','H',out.h_pp_su,'P',out.P_pp_su,fluid_wf);
    out.Mass.name{1,i_mass} = 'M_liq_receiver';
    out.level_receiver = 1;
    out.Receiver_type = 1;
elseif abs((out.h_pp_su-CoolProp.PropsSI('H','Q',0,'P',out.P_pp_su,fluid_wf))/out.h_pp_su)<1e-3   
    i_mass = i_mass + 1;
    out.Mass.value(1,i_mass) = max(param.V_liq_rec*CoolProp.PropsSI('D','Q',1,'P',out.P_pp_su,fluid_wf),min(param.M_tot-sum(out.Mass.value),param.V_liq_rec*CoolProp.PropsSI('D','Q',0,'P',out.P_pp_su,fluid_wf)));
    out.Mass.name{1,i_mass} = 'M_liq_receiver';
    out.Receiver_type = 2;
    out.level_receiver = (out.Mass.value(1,i_mass)-param.V_liq_rec*CoolProp.PropsSI('D','Q',1,'P',out.P_pp_su,fluid_wf))/((param.V_liq_rec*CoolProp.PropsSI('D','Q',0,'P',out.P_pp_su,fluid_wf))-(param.V_liq_rec*CoolProp.PropsSI('D','Q',1,'P',out.P_pp_su,fluid_wf)));
else
    i_mass = i_mass + 1;
    out.Mass.value(1,i_mass) = param.V_liq_rec*CoolProp.PropsSI('D','Q',1,'P',out.P_pp_su,fluid_wf);
    out.Mass.name{1,i_mass} = 'M_liq_receiver';
    out.Receiver_type = 3;
    out.level_receiver = 0;
end
    
out.M_tot = sum(out.Mass.value);

% ORC PERFORAMANCE
out.W_dot_net = out.W_dot_exp - out.W_dot_pp - out.W_dot_cd;
out.eff_ORC_gross = out.W_dot_exp/out.Q_dot_ev;
out.eff_ORC_net = out.W_dot_net/out.Q_dot_ev;

% RESIDUALS and RESULTS
out.res_ORC_H_cd_ex = (1 - out.h_dplp_su/out.h_cd_ex);
out.res_ORC_H_ev_ex = (1 - out.h_dphp_su/out.h_ev_ex);
out.res_ORC_N_exp = 1-out.N_exp_bis/N_exp;
out.res_ORC_M = (1 - out.M_tot/param.M_tot);

if strcmp(param.solverType, 'M_imposed')
    out.res_vec  = [out.res_ORC_N_exp    out.res_ORC_H_cd_ex     out.res_ORC_H_ev_ex    out.res_ORC_M];
elseif strcmp(param.solverType, 'DTsc_imposed')
    out.res_vec  = [out.res_ORC_N_exp    out.res_ORC_H_cd_ex     out.res_ORC_H_ev_ex];
end
out.res  = norm(out.res_vec);
% if out.flag_rec < 0
%     out.res = out.res*1e3;
% end
out.x = x;

if strcmp(param.eval_type, 'fast')
    TS = NaN;
else
    out = orderfields(out);  
    TS.EXP = TS_EXP;
    TS.PP = TS_PP;
    TS.REC = TS_REC;
    TS.CD = TS_CD;
    TS.EV = TS_EV;
    TS.LossesHP = TS_LossesHP;
    TS.LossesLP = TS_LossesLP;
    
    % Construction of the TS variable
    TS.cycle.s = TS.PP.s;
    TS.cycle.T = TS.PP.T;
    if isfield(param, 'REC')
        TS.cycle.s = [TS.cycle.s  TS.REC.s_c];
        TS.cycle.T = [TS.cycle.T  TS.REC.T_c];
    end
    if isfield(param, 'PRE')
        TS.cycle.s = [TS.cycle.s  TS.PRE.s_c];
        TS.cycle.T = [TS.cycle.T  TS.PRE.T_c];
    end
    TS.cycle.s = [TS.cycle.s  TS.EV.s_c];
    TS.cycle.T = [TS.cycle.T  TS.EV.T_c];
    if isfield(param, 'LossesHP')
        TS.cycle.s = [TS.cycle.s TS.LossesHP.s];
        TS.cycle.T = [TS.cycle.T TS.LossesHP.T];
    end
    TS.cycle.s = [TS.cycle.s  TS.EXP.s];
    TS.cycle.T = [TS.cycle.T  TS.EXP.T];
    if isfield(param, 'REC')
        TS.cycle.s = [TS.cycle.s  fliplr(TS.REC.s_h)];
        TS.cycle.T = [TS.cycle.T  fliplr(TS.REC.T_h)];
    end
    TS.cycle.s = [TS.cycle.s fliplr(TS.CD.s_h)];
    TS.cycle.T = [TS.cycle.T fliplr(TS.CD.T_h)];
    if isfield(param, 'SUB')
        TS.cycle.s = [TS.cycle.s  fliplr(TS.SUB.s_h)];
        TS.cycle.T = [TS.cycle.T  fliplr(TS.SUB.T_h)];
    end
    if isfield(param, 'LossesLP')
        TS.cycle.s = [TS.cycle.s TS.LossesLP.s];
        TS.cycle.T = [TS.cycle.T TS.LossesLP.T];
    end
end
end
