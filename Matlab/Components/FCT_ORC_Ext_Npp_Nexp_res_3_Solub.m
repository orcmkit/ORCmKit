function res =  FCT_ORC_Ext_Npp_Nexp_res_3_Solub( x,    lb, ub, fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param)
[out, ~] =      FCT_ORC_Ext_Npp_Nexp_3_Solub(x,         lb, ub, fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param);
res = out.res;
end