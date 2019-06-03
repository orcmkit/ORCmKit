function [TS] = generate_TS_solub(in, out, param)
C_oil = in.C_oil;
fluid_wf = 'R245fa';
fluid_lub = 'ICP_RL_32_3MAF';
load(['TS_prop_',fluid_wf , '.mat'])

H_vec = [out.h_pp_su    out.h_pp_ex    out.h_recc_su   out.h_recc_ex    out.h_ev_su   out.h_ev_ex  out.h_exp_su     out.h_rech_su      out.h_rech_ex    out.h_cd_su     out.h_cd_ex   out.h_pp_su];
P_vec = [out.P_pp_su    out.P_pp_ex    out.P_recc_su   out.P_recc_ex    out.P_ev_su   out.P_ev_ex  out.P_exp_su     out.P_rech_su      out.P_rech_ex    out.P_cd_su     out.P_cd_ex   out.P_pp_su];
for ii=1:length(H_vec)
    if C_oil > 0
        Tsat_pure_K = CoolProp.PropsSI('T', 'P', P_vec(ii), 'Q', 0.5, fluid_wf);
        Tbubble_min_K =  R245fa_POE_Tbubble(1-C_oil, P_vec(ii), Tsat_pure_K);
        [T(ii), h_rl(ii), ~, h_rv(ii), ~, x(ii), ~, ~] = HP_solubMixt(C_oil, P_vec(ii), H_vec(ii), fluid_wf, fluid_lub, Tbubble_min_K, Tsat_pure_K, []);
        h_r_eq(ii) = x(ii)*h_rv(ii) + (1-x(ii))*h_rl(ii);
    else
        h_r_eq(ii) = H_vec(ii);
        T(ii) = CoolProp.PropsSI('T', 'P', P_vec(ii), 'H', h_r_eq(ii), fluid_wf);
    end
    s(ii) = CoolProp.PropsSI('S', 'P', P_vec(ii), 'H', h_r_eq(ii), fluid_wf);
    clearvars Tbubble_min_K Tsat_pure_K
end

% EV heat transfer
param.EV.C.C_oil = C_oil;
param.EV.n_disc = 20;
if param.EV.C.solub
    if param.EV.C.C_oil < 1e-4
        param.EV.C.x_lim_vap = 0.999; 
    else
        param.EV.C.x_lim_vap = 0.98;
    end
end
[h_h_l_in, h_h_v_in, ~] = supply_conditions_HEX_Solub(in.h_htf_su, in.P_htf_su, in.fluid_htf, param.EV.H);
[h_c_l_in, h_c_v_in, ~] = supply_conditions_HEX_Solub(out.h_ev_su, out.P_ev_su, fluid_wf, param.EV.C);
EV_profile = HEX_profile_Solub_DP2(in.fluid_htf,  in.m_dot_htf, in.P_htf_su, in.h_htf_su, fluid_wf, out.m_dot_wf, out.P_ev_su, out.h_ev_su, out.Q_dot_ev, param.EV, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, 0, out.P_ev_su-out.P_ev_ex);
for i=1:length(EV_profile.C.T_vec)
    T_ev_vec(i) = EV_profile.C.T_vec(i);
    if C_oil >0
        s_ev_vec(i) = CoolProp.PropsSI('S', 'P', EV_profile.C.P_vec(i), 'H', EV_profile.C.x_vec(i)*EV_profile.C.H_rv_vec(i) + (1-EV_profile.C.x_vec(i))*EV_profile.C.H_rl_vec(i), fluid_wf);
    else
        s_ev_vec(i) = CoolProp.PropsSI('S', 'P', EV_profile.C.P_vec(i), 'H', EV_profile.C.H_vec(i), fluid_wf);
    end
end
clearvars h_h_l_in h_h_v_in h_c_l_in h_c_v_in

% CD heat transfer
param.CD.H.C_oil = C_oil;
param.CD.n_disc = 20;
if param.CD.H.solub
    if param.CD.H.C_oil < 1e-4
        param.CD.H.x_lim_vap = 0.999;
    else
        param.CD.H.x_lim_vap = 0.98;
    end
end
[h_c_l_in, h_c_v_in, ~] = supply_conditions_HEX_Solub(in.h_ctf_su, in.P_ctf_su, in.fluid_ctf, param.CD.C);
[h_h_l_in, h_h_v_in, ~] = supply_conditions_HEX_Solub(out.h_cd_su, out.P_cd_su, fluid_wf, param.CD.H);
CD_profile = HEX_profile_Solub_DP2(fluid_wf, out.m_dot_wf, out.P_cd_su, out.h_cd_su, in.fluid_ctf,  in.m_dot_ctf, in.P_ctf_su, in.h_ctf_su, out.Q_dot_cd, param.CD, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, out.P_cd_su-out.P_cd_ex, 0);
for i=1:length(CD_profile.H.T_vec)
    T_cd_vec(i) = CD_profile.H.T_vec(i);
    if C_oil > 0
        s_cd_vec(i) = CoolProp.PropsSI('S', 'P', CD_profile.H.P_vec(i), 'H', CD_profile.H.x_vec(i)*CD_profile.H.H_rv_vec(i) + (1-CD_profile.H.x_vec(i))*CD_profile.H.H_rl_vec(i), fluid_wf);
    else
        s_cd_vec(i) = CoolProp.PropsSI('S', 'P', CD_profile.H.P_vec(i), 'H', CD_profile.H.H_vec(i), fluid_wf);
    end
end
T_cd_vec = fliplr(T_cd_vec);
s_cd_vec = fliplr(s_cd_vec);
clearvars h_h_l_in h_h_v_in h_c_l_in h_c_v_in

%REC heat transfer
param.REC.H.C_oil = C_oil;
param.REC.C.C_oil = C_oil;
if param.REC.H.solub
    if param.REC.H.C_oil < 1e-4
        param.REC.H.x_lim_vap = 0.999;
    else
        param.REC.H.x_lim_vap = 0.98;
    end
end
if param.REC.C.solub
    if param.REC.C.C_oil < 1e-4
        param.REC.C.x_lim_vap = 0.999;
    else
        param.REC.C.x_lim_vap = 0.98;
    end
end
[h_c_l_in, h_c_v_in, ~] = supply_conditions_HEX_Solub(out.h_recc_su, out.P_recc_su, fluid_wf, param.REC.C);
[h_h_l_in, h_h_v_in, ~] = supply_conditions_HEX_Solub(out.h_rech_su, out.P_rech_su, fluid_wf, param.REC.H);
REC_profile = HEX_profile_Solub_DP2(fluid_wf, out.m_dot_wf, out.P_rech_su, out.h_rech_su, fluid_wf, out.m_dot_wf, out.P_recc_su, out.h_recc_su, out.Q_dot_rec, param.REC, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, out.P_rech_su-out.P_rech_ex, 0);
for i=1:length(REC_profile.H.T_vec)
    T_rech_vec(i) = REC_profile.H.T_vec(i);
    T_recc_vec(i) = REC_profile.C.T_vec(i);
    if C_oil > 0
        s_rech_vec(i) = CoolProp.PropsSI('S', 'P', REC_profile.H.P_vec(i), 'H', REC_profile.H.x_vec(i)*REC_profile.H.H_rv_vec(i) + (1-REC_profile.H.x_vec(i))*REC_profile.H.H_rl_vec(i), fluid_wf);
        s_recc_vec(i) = CoolProp.PropsSI('S', 'P', REC_profile.C.P_vec(i), 'H', REC_profile.C.x_vec(i)*REC_profile.C.H_rv_vec(i) + (1-REC_profile.C.x_vec(i))*REC_profile.C.H_rl_vec(i), fluid_wf);
    else
        s_rech_vec(i) = CoolProp.PropsSI('S', 'P', REC_profile.H.P_vec(i), 'H', REC_profile.H.H_vec(i), fluid_wf);
        s_recc_vec(i) = CoolProp.PropsSI('S', 'P', REC_profile.C.P_vec(i), 'H', REC_profile.C.H_vec(i), fluid_wf);
    end
end
T_rech_vec = fliplr(T_rech_vec);
s_rech_vec = fliplr(s_rech_vec);
clearvars h_h_l_in h_h_v_in h_c_l_in h_c_v_in

TS.T_cycle = [T(1:3) T_recc_vec T(4:5) T_ev_vec T(6:8) T_rech_vec T(9:10) T_cd_vec T(11:12)];
TS.s_cycle = [s(1:3) s_recc_vec s(4:5) s_ev_vec s(6:8) s_rech_vec s(9:10) s_cd_vec s(11:12)];
TS.s_comp_lim = s;
TS.T_com_lim = T;
TS.s_htf_ev = [s(6) s(5)];
TS.T_htf_ev = [in.T_htf_su out.T_htf_ev_ex];
TS.s_ctf_cd = [s(11) s(10)];
TS.T_ctf_cd = [in.T_ctf_su out.T_ctf_cd_ex];

hold on
TS.ts_curve = plot(s_TS_curve, T_TS_curve-273.15, 'k', 'linewidth', 1.5);
TS.ctf_cd =     plot(TS.s_ctf_cd,       TS.T_ctf_cd-273.15, '--');
TS.htf_ev =     plot(TS.s_htf_ev,       TS.T_htf_ev-273.15, '--');
TS.cycle =      plot(TS.s_cycle,        TS.T_cycle-273.15, '-');
TS.comp_lim =   plot(TS.s_comp_lim,     TS.T_com_lim-273.15, 'o', 'linestyle','none');
hold off

end

