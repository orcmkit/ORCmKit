function res = HEX_hConvCorSolub_res_Fast(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info, h_h_l, h_h_v, h_c_l, h_c_v, DP_h, DP_c)
out = HEX_hConvCorSolub_Fast(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, real(Q_dot), info, h_h_l, h_h_v, h_c_l, h_c_v, DP_h, DP_c);
res = out.resA;
end
