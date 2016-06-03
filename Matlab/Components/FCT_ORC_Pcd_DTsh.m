function [out, TS] = FCT_ORC_Pcd_DTsh( x, lb, ub, fluid_wf, fluid_htf, in_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, P_ctf_su, T_amb, N_exp, P_cd, DT_sh, param)
x = max(x, lb./ub);
x = min(x,ones(1,length(x)));
x = x.*ub;

i_flag = 0;
i_mass = 0;

% Expander inlet
out.P_exp_su = x(1);
out.h_exp_su = CoolProp.PropsSI('H', 'T', CoolProp.PropsSI('T', 'P', out.P_exp_su, 'Q', 1, fluid_wf)+ DT_sh, 'P', out.P_exp_su, fluid_wf);
out.P_rech_su = P_cd;

% Expander model
[out_EXP,TS_EXP] = ExpanderModel(fluid_wf, out.P_exp_su, out.h_exp_su, N_exp, out.P_rech_su, T_amb, param.EXP);
out.m_dot_wf = out_EXP.M_dot;
out.h_rech_su = out_EXP.h_ex;
out.W_dot_exp = out_EXP.W_dot;
out.Q_dot_exp = out_EXP.Q_dot_amb;
out.eps_vol_exp = out_EXP.FF;
out.eps_is_exp = out_EXP.epsilon_is;
out.flag_exp = out_EXP.flag;
out.time_exp = out_EXP.time;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_EXP.flag;
out.flag.name{1,i_flag} = 'flag_exp';
i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = out_EXP.M;
out.Mass.name{1,i_mass} = 'M_exp';
TS.EXP = TS_EXP;

% High pressure losses model
param.LossesHP.type_in = 'ex';
[out_LossesHP, TS_LossesHP] = LossesModel(fluid_wf, out.P_exp_su, out.h_exp_su, out.m_dot_wf, T_amb, param.LossesHP);
out.P_dphp_su = out_LossesHP.P_su;
out.h_dphp_su = out_LossesHP.h_su;
out.dphp = out_LossesHP.dp;
out.Q_dot_hp = out_LossesHP.Q_dot;
out.flag_lossesHP = out_LossesHP.flag;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_LossesHP.flag;
out.flag.name{1,i_flag} = 'flag_dphp';
TS.LossesHP = TS_LossesHP;

% Low pressure losses model
param.LossesLP.type_in = 'su';
[out_LossesLP, TS_LossesLP] = LossesModel(fluid_wf, out.P_rech_su, out.h_rech_su, out.m_dot_wf, T_amb, param.LossesLP);
out.P_pp_su = out_LossesLP.P_ex;
out.dplp = out_LossesLP.dp;
out.Q_dot_lp = out_LossesLP.Q_dot;
out.flag_lossesLP = out_LossesLP.flag;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_LossesLP.flag;
out.flag.name{1,i_flag} = 'flag_dplp';
TS.LossesLP = TS_LossesLP;

% Pump model
out.h_pp_su = x(2);
[out_PP,TS_PP] = PpModel_2_SemiEmp(out.P_pp_su, out.h_pp_su, out.P_dphp_su, fluid_wf, out.m_dot_wf, param.PP);
out.h_recc_su = out_PP.h_ex;
out.P_recc_su = out.P_dphp_su;
out.N_pp = out_PP.N_pp;
out.W_dot_pp = out_PP.W_dot;
out.eps_is_pp = out_PP.epsilon_is;
out.eps_vol_pp = out_PP.epsilon_vol;
out.time_pp = out_PP.time;
out.flag_pp = out_PP.flag;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_PP.flag;
out.flag.name{1,i_flag} = 'flag_pp';
i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = out_PP.M;
out.Mass.name{1,i_mass} = 'M_pp';
TS.PP = TS_PP;

% Recuperator model
[out_REC, TS_REC] = HexModel(fluid_wf, out.P_rech_su, out.h_rech_su, out.m_dot_wf, fluid_wf, out.P_recc_su, out.h_recc_su, out.m_dot_wf, param.REC);
out.h_ev_su = out_REC.h_h_ex;
out.P_ev_su = out.P_recc_su;
out.h_cd_su = out_REC.h_c_ex;
out.P_cd_su = out.P_rech_su;
out.Q_dot_rec = out_REC.Q_dot_tot;
out.flag_rec = out_REC.flag;
out.time_rec = out_EXP.time;
out.pinch_rec = out_REC.pinch;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_REC.flag;
out.flag.name{1,i_flag} = 'flag_rec';
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = out_REC.M_c;
out.Mass.name{1,i_mass} = 'M_recc';
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = out_REC.M_h;
out.Mass.name{1,i_mass} = 'M_rech';
TS.REC = TS_REC;

% Evaporator model
[out_EV, TS_EV] = HexModel(fluid_htf, P_htf_su, in_htf_su, m_dot_htf, fluid_wf, out.P_ev_su, out.h_ev_su, out.m_dot_wf , param.EV);
out.h_dphp_su_bis = out_EV.h_c_ex;
out.Q_dot_ev = out_EV.Q_dot_tot;
out.flag_ev = out_EV.flag;
out.time_ev = out_EV.time;
out.pinch_ev = out_EV.pinch;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_EV.flag;
out.flag.name{1,i_flag} = 'flag_ev';
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = out_EV.M_c;
out.Mass.name{1,i_mass} = 'M_ev';
TS.EV = TS_EV;

% Condenser model
out.h_dplp_su = out.h_pp_su;
out.P_dplp_su = out.P_cd_su;
out.Q_dot_cd = out.m_dot_wf*(out.h_cd_su-out.h_dplp_su);
[out_CD,TS_CD] = HexModel_2_hConvVar(fluid_wf, out.P_cd_su, out.h_cd_su, out.m_dot_wf, fluid_ctf, P_ctf_su, in_ctf_su, out.Q_dot_cd, param.CD);
out.Q_dot_cd = out_CD.Q_dot_tot;
out.T_ctf_cd_ex = out_CD.T_c_ex;
out.h_ctf_cd_ex = out_CD.h_c_ex;
out.m_dot_ctf = out_CD.m_dot_c;
out.flag_cd = out_CD.flag;
out.time_cd = out_CD.time;
out.pinch_cd = out_CD.pinch;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_CD.flag;
out.flag.name{1,i_flag} = 'flag_cd';
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = out_CD.M_h;
out.Mass.name{1,i_mass} = 'M_cd';
TS.CD = TS_CD;

% Mass between evaporator and expander
i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_ev_ex*CoolProp.PropsSI('D','H',out.h_exp_su,'P',out.P_exp_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_ev_ex';

% Mass between expander and recuperator
i_mass = i_mass+1;
out.Mass.value(1,i_mass) = param.V_aux_exp_ex*CoolProp.PropsSI('D','H',out.h_rech_su,'P',out.P_rech_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_exp_ex';

% Mass between recuperator and condenser
i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_rech_ex*CoolProp.PropsSI('D','H',out.h_cd_su,'P',out.P_cd_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_rech_ex';

% Mass between condenser and pump
i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_cd_ex*CoolProp.PropsSI('D','H',out.h_pp_su,'P',out.P_pp_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_cd_ex';

% Mass between pump and recuperator
i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_pp_ex*CoolProp.PropsSI('D','H',out.h_recc_su,'P',out.P_recc_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_pp_ex';

% Mass between recuperator and evaporator
i_mass = i_mass + 1;
out.Mass.value(1,i_mass) = param.V_aux_recc_ex*CoolProp.PropsSI('D','H',out.h_ev_su,'P',out.P_ev_su,fluid_wf);
out.Mass.name{1,i_mass} = 'M_aux_recc_ex';


% Liquid receiver
if out.h_pp_su<CoolProp.PropsSI('D','Q',0,'P',out.P_pp_su,fluid_wf)
    i_mass = i_mass + 1;
    out.Mass.value(1,i_mass) = param.V_liq_rec*CoolProp.PropsSI('D','H',out.h_pp_su,'P',out.P_pp_su,fluid_wf);
    out.Mass.name{1,i_mass} = 'M_liq_receiver';
elseif abs((out.h_pp_su-CoolProp.PropsSI('D','Q',0,'P',out.P_pp_su,fluid_wf))/out.h_dplp_su)<1e-5   
    i_mass = i_mass + 1;
    out.Mass.value(1,i_mass) = max(param.V_liq_rec*CoolProp.PropsSI('D','Q',1,'P',out.P_pp_su,fluid_wf),min(param.M_tot-sum(out.Mass.value),param.V_liq_rec*CoolProp.PropsSI('D','Q',0,'P',out.P_pp_su,fluid_wf)));
    out.Mass.name{1,i_mass} = 'M_liq_receiver';
else
    i_mass = i_mass + 1;
    out.Mass.value(1,i_mass) = param.V_liq_rec*CoolProp.PropsSI('D','Q',1,'P',out.P_pp_su,fluid_wf);
    out.Mass.name{1,i_mass} = 'M_liq_receiver';
end

out.M_tot = sum(out.Mass.value);

% RESIDUALS and RESULTS
out.res_ORC_Hev_ex = (1 - out.h_dphp_su/out.h_dphp_su_bis);
out.res_ORC_M = (1 - out.M_tot/param.M_tot);
out.res_vec  = [out.res_ORC_Hev_ex out.res_ORC_M];
out.res  = norm(out.res_vec);
out.x = x;

end
