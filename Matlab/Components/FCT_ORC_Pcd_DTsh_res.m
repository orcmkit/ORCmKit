function res = FCT_ORC_Pcd_DTsh_res( x, lb, ub, fluid_wf, fluid_htf, in_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, P_ctf_su, T_amb, N_exp, P_cd, DT_sh, param)
out = FCT_ORC_Pcd_DTsh_fast( x, lb, ub, fluid_wf, fluid_htf, in_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, P_ctf_su, T_amb, N_exp, P_cd, DT_sh, param);
res = out.res;
end