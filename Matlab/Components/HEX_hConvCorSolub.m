function out = HEX_hConvCorSolub(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, DP_h, DP_c)
out = HEX_profile_Solub_DP2(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, DP_h, DP_c);

[out.H.A_vec, out.C.A_vec, out.H.hConv_vec, out.C.hConv_vec, out.H.Nu_vec, out.C.Nu_vec, out.H.fConv_vec, out.C.fConv_vec, out.DTlog, out.F, out.H.eff_vec, out.C.eff_vec, out.AU_vec, out.U_vec, out.k] = deal(NaN*ones(1,length(out.H.H_vec)-1));
x_di_c = 1; dry_out_c = 0;
disp_flag = 0;
for j = 1:length(out.H.T_vec)-1
    
    % LMTD for the current cell
    out.DTlog(j) = deltaT_log(out.H.T_vec(j+1), out.H.T_vec(j),out.C.T_vec(j), out.C.T_vec(j+1));
    if strcmp(info.typeHEX, 'CrossFlow')  
        R = (out.H.T_vec(j+1)- out.H.T_vec(j))/(out.C.T_vec(j+1)-out.C.T_vec(j));
        P = (out.C.T_vec(j+1)-out.C.T_vec(j))/(out.H.T_vec(j+1)-out.C.T_vec(j));
        out.F(j) = F_lmtd2(R, P);
    else
        out.F(j) = 1;
    end
    T_wall = (out.H.T_vec(j+1)+ out.H.T_vec(j)+out.C.T_vec(j)+ out.C.T_vec(j+1))/4;
    
    % What type of cells for hot side (1phase/2phase)
    if strcmp(info.H.type, 'H')
        if (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)) < out.H.h_l
            out.H.type_zone{j} = 'liq';
            T_wall_h = T_wall;
        elseif (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)) > out.H.h_v
            if T_wall>(0.5*out.H.Tsat_pure_vec(j)+ 0.5*out.H.Tsat_pure_vec(j+1))
                out.H.type_zone{j} = 'vap';
                T_wall_h = max(T_wall, (0.5*out.H.Tsat_pure_vec(j)+ 0.5*out.H.Tsat_pure_vec(j+1))+5e-2);
            else
                out.H.type_zone{j} = 'vap_wet';
                T_wall_h = max(T_wall, (0.5*out.H.Tsat_pure_vec(j)+ 0.5*out.H.Tsat_pure_vec(j+1))+5e-2);
            end
        else
            out.H.type_zone{j} = 'tp';
            T_wall_h = T_wall;
        end
    elseif strcmp(info.H.type, 'T')
        out.H.type_zone{j} = 'liq';
        T_wall_h = T_wall;
    end
    
    % What type of cells for cold side (1phase/2phase)
    if strcmp(info.C.type, 'H')
        if (0.5*out.C.H_vec(j)+0.5*out.C.H_vec(j+1)) < out.C.h_l
            out.C.type_zone{j} = 'liq';
            T_wall_c = min(T_wall, (0.5*out.C.Tsat_pure_vec(j)+0.5*out.C.Tsat_pure_vec(j+1))-5e-2);
        elseif (0.5*out.C.H_vec(j)+0.5*out.C.H_vec(j+1)) > out.C.h_v
            out.C.type_zone{j} = 'vap';
            T_wall_c = T_wall;
        else
            if (0.5*out.C.x_vec(j)+0.5*out.C.x_vec(j+1)) <= x_di_c
                out.C.type_zone{j} = 'tp';
                T_wall_c = T_wall;
            else
                out.C.type_zone{j} = 'tp_dryout';
                T_wall_c = T_wall;
                dry_out_c = 1;
            end
        end
        
    elseif strcmp(info.C.type, 'T')
        out.C.type_zone{j} = 'liq';
    end
       
    % Hot-side : single phase convective heat transfer coefficient
    if strcmp(out.H.type_zone{j}, 'liq') || strcmp(out.H.type_zone{j}, 'vap') || strcmp(out.H.type_zone{j}, 'vap_wet')
        G_h = m_dot_h/info.H.n_canals/info.H.CS;
        T_h_mean = (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1));
        P_h_mean = (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1));
        T_h_sat_mean = (0.5*out.H.Tsat_pure_vec(j)+0.5*out.H.Tsat_pure_vec(j+1));
        [mu_h, Pr_h, k_h, mu_rat_h] = PropsSinglePhase(T_h_mean, P_h_mean, T_wall_h, T_h_sat_mean, out.H.type_zone{j}, fluid_h);
        

        if strcmp(out.H.type_zone{j}, 'liq')
            type_correlation_h = info.H.correlation.type_1phase_l;
        elseif strcmp(out.H.type_zone{j}, 'vap') || strcmp(out.H.type_zone{j}, 'vap_wet')
            type_correlation_h = info.H.correlation.type_1phase_v;
        end
        
        switch type_correlation_h
            
            case 'UD'
                if strcmp(out.H.type_zone{j}, 'liq')
                    h_nom = info.H.h_nom_liq;
                    m_dot_nom = info.H.m_dot_nom_liq;
                    n_nom = info.H.n_nom_liq;
                elseif strcmp(out.H.type_zone{j}, 'vap') || strcmp(out.H.type_zone{j}, 'vap_wet')
                    h_nom = info.H.h_nom_vap;
                    m_dot_nom = info.H.m_dot_nom_vap;
                    n_nom = info.H.n_nom_vap;
                end
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = UD_HTC(m_dot_h, h_nom, m_dot_nom, n_nom);
                
            case 'S2P_recLiq_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = S2P_recLiq_BPHEX_HTC(mu_h, mu_rat_h, Pr_h, k_h, G_h, info.H.Dh, disp_flag);

            case 'S2P_recVap_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = S2P_recVap_BPHEX_HTC(mu_h, mu_rat_h, Pr_h, k_h, G_h, info.H.Dh, disp_flag);
            
            case 'S2P_rec1p_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = S2P_rec1p_BPHEX_HTC(mu_h, mu_rat_h, Pr_h, k_h, G_h, info.H.Dh, disp_flag);
                
            case 'S2P_rec1phaseOK_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = S2P_rec1phase_BPHEX_HTC(mu_h, mu_rat_h, Pr_h, k_h, G_h, info.H.Dh, disp_flag); 
                
            case 'S2P_evHTF_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = S2P_evHTF_BPHEX_HTC(mu_h, mu_rat_h, Pr_h, k_h, G_h, info.H.Dh, disp_flag);
                
            case 'Martin1_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = Martin1_BPHEX_HTC(mu_h, mu_rat_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, disp_flag);
                
            case 'Martin2_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = Martin2_BPHEX_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, disp_flag);
                
            case 'Wanniarachchi_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = Wanniarachchi_BPHEX_HTC(mu_h, mu_rat_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, info.phi, disp_flag);
                
            case 'Thonon_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = Thonon_BPHEX_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, disp_flag);
                
            case 'Junqi_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = Junqi_BPHEX_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, disp_flag);
                
            case 'Muley_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = Muley_BPHEX_HTC(mu_h, mu_rat_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, info.phi, info.L_hex, disp_flag);
                
            case 'Kim_BPHEX'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = Kim_BPHEX_HTC(mu_h,Pr_h, k_h, G_h, info.H.Dh, info.theta, disp_flag);
                
            case 'DittusBoelter_Pipe'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = DittusBoelter_Pipe_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, 0.3, disp_flag);
                
            case 'Gnielinski_Pipe'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = Gnielinski_Pipe_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.L_hex, disp_flag);
                
            case 'VDI_finnedTubes_staggered'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = VDI_finnedTubes_staggered_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dc, info.H.fin.omega_t, disp_flag);
                
            case 'Wang_finnedTubes_staggered'
                [hConv_1phase_h, Nu_1phase_h, flag_1phase_h] = Wang_finnedTubes_staggered_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.H.Dc, info.H.Pt, info.H.Fp, info.H.Nr, disp_flag);

        end
        if 1
            %%%%
            if strcmp(out.H.type_zone{j}, 'liq')
                hConv_1phase_h = info.C_fit_1pl_h*hConv_1phase_h;
                Nu_1phase_h = info.C_fit_1pl_h*Nu_1phase_h;
            elseif strcmp(out.H.type_zone{j}, 'vap')
                hConv_1phase_h = info.C_fit_1pv_h*hConv_1phase_h;
                Nu_1phase_h = info.C_fit_1pv_h*Nu_1phase_h;
            end
            %%%%
        end
        
        if strcmp(out.H.type_zone{j}, 'liq') || strcmp(out.H.type_zone{j}, 'vap')
            out.H.hConv_vec(j) = hConv_1phase_h;
            out.H.Nu_vec(j) = Nu_1phase_h;
            out.H.fConv_vec(j) = flag_1phase_h;
        end
        
    end
    
    % Hot-side: two-phase convective heat transfer coefficient
    if strcmp(out.H.type_zone{j}, 'tp') || strcmp(out.H.type_zone{j}, 'vap_wet')
        
        G_h = (m_dot_h/info.H.n_canals)/info.H.CS;
        if strcmp(out.H.type_zone{j}, 'tp')
            x_h = min(1,max(0,0.5*out.H.x_vec(j) + 0.5*out.H.x_vec(j+1)));
        else
            x_h = 1;
        end
        mu_h_l = CoolProp.PropsSI('V', 'Q', 0, 'P', (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1)), fluid_h); %to be updated with mixture properties
        k_h_l = CoolProp.PropsSI('L', 'Q', 0, 'P', (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1)), fluid_h); %to be updated with mixture properties
        Pr_h_l = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1)), fluid_h); %to be updated with mixture properties
        rho_h_l = CoolProp.PropsSI('D', 'Q', 0, 'P', (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1)), fluid_h); %to be updated with mixture properties
        rho_h_v = CoolProp.PropsSI('D', 'Q', 1, 'P', (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1)), fluid_h); %to be updated with super heated vapor properties
        i_fg_h = CoolProp.PropsSI('H', 'Q', 1, 'P', (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1)), fluid_h) - CoolProp.PropsSI('H', 'Q', 0, 'P', (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1)), fluid_h);
        
        
        switch info.H.correlation.type_2phase_cd
            
            case 'UD'
                h_nom = info.H.h_nom_tp;
                m_dot_nom = info.H.m_dot_nom_tp;
                n_nom = info.H.n_nom_tp;
                [hConv_2phase_h, Nu_2phase_h, flag_2phase_h] = UD_HTC(m_dot_h, h_nom, m_dot_nom, n_nom);
                
                
            case 'Han_cond_BPHEX'
                [hConv_2phase_h, Nu_2phase_h, flag_2phase_h] = Han_Cond_BPHEX_HTC(x_h, mu_h_l, k_h_l, Pr_h_l, rho_h_l, rho_h_v, G_h, info.H.Dh, info.pitch_co, info.theta, disp_flag);
                
            case 'Longo_cond_BPHEX'
                [hConv_2phase_h, Nu_2phase_h, flag_2phase_h]  = Longo_Cond_BPHEX_HTC(x_h, mu_h_l, k_h_l, Pr_h_l, rho_h_l, rho_h_v, i_fg_h, G_h, info.H.Dh, (0.5*out.H.T_vec(j+1)+ 0.5*out.H.T_vec(j))-T_wall, info.phi, info.L_hex, disp_flag);
                
            case 'Shah_cond_Pipe'
                p_h_star = (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1))/CoolProp.PropsSI('Pcrit', 'Q', 1, 'P', (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1)), fluid_h);
                [hConv_2phase_h, Nu_2phase_h, flag_2phase_h] = Shah_Cond_pipe_HTC(x_h, mu_h_l, k_h_l, Pr_h_l, p_h_star, G_h, info.H.Dh, disp_flag);
                
            case 'Cavallini_cond_Pipe'
                mu_h_v = CoolProp.PropsSI('V', 'Q', 1, 'P', (0.5*out.H.P_vec(j)+0.5*out.H.P_vec(j+1)), fluid_h);                 
                [hConv_2phase_h, Nu_2phase_h, flag_2phase_h] = Cavallini_Cond_pipe_HTC(x_h, mu_h_l, mu_h_v, rho_h_l, rho_h_v, k_h_l, Pr_h_l, i_fg_h, (0.5*out.H.T_vec(j+1)+ 0.5*out.H.T_vec(j))-T_wall, G_h, info.H.Dh, disp_flag);
                

        end
        
        if 1
            %%%%
            hConv_2phase_h = info.C_fit_2p_h*hConv_2phase_h;
            Nu_2phase_h = info.C_fit_2p_h*Nu_2phase_h;
            %%%%
        end
               
        if strcmp(out.H.type_zone{j}, 'tp')
            out.H.hConv_vec(j) = hConv_2phase_h;
            out.H.Nu_vec(j) = Nu_2phase_h;
            out.H.fConv_vec(j) = flag_2phase_h;
        elseif strcmp(out.H.type_zone{j}, 'vap_wet')
            w_vap_wet = ((0.5*out.H.T_vec(j)+ 0.5*out.H.T_vec(j+1)) - (0.5*out.H.Tsat_pure_vec(j)+ 0.5*out.H.Tsat_pure_vec(j+1)))/((0.5*out.H.T_vec(j)+ 0.5*out.H.T_vec(j+1))-T_wall);
            out.H.hConv_vec(j) = hConv_2phase_h - w_vap_wet*(hConv_2phase_h - hConv_1phase_h);
            out.H.Nu_vec(j) = Nu_2phase_h - w_vap_wet*(Nu_2phase_h - Nu_1phase_h);
            out.H.fConv_vec(j) = min(flag_1phase_h,flag_2phase_h);
        end
        
    end
           
    % Cold-side : single phase convective heat transfer coefficient
    if strcmp(out.C.type_zone{j}, 'liq') || strcmp(out.C.type_zone{j}, 'vap') || strcmp(out.C.type_zone{j}, 'tp_dryout')
        G_c = m_dot_c/info.C.n_canals/info.C.CS;
        T_c_mean = (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1));
        T_c_sat_mean = (0.5*out.C.Tsat_pure_vec(j)+0.5*out.C.Tsat_pure_vec(j+1));
        P_c_mean = (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1));
        [mu_c, Pr_c, k_c, mu_rat_c] = PropsSinglePhase(T_c_mean, P_c_mean, T_wall_c, T_c_sat_mean, out.C.type_zone{j}, fluid_c);

        if strcmp(out.C.type_zone{j}, 'liq')
            type_correlation_c = info.C.correlation.type_1phase_l;
        elseif strcmp(out.C.type_zone{j}, 'vap') || strcmp(out.C.type_zone{j}, 'tp_dryout')
            type_correlation_c = info.C.correlation.type_1phase_v;
        end
        
        switch type_correlation_c
            
            case 'UD'
                if strcmp(out.C.type_zone{j}, 'liq')
                    h_nom = info.C.h_nom_liq;
                    m_dot_nom = info.C.m_dot_nom_liq;
                    n_nom = info.C.n_nom_liq;
                elseif strcmp(out.C.type_zone{j}, 'vap')
                    h_nom = info.C.h_nom_vap;
                    m_dot_nom = info.C.m_dot_nom_vap;
                    n_nom = info.C.n_nom_vap;
                end
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = UD_HTC(m_dot_h, h_nom, m_dot_nom, n_nom);
                
            case 'S2P_recLiq_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = S2P_recLiq_BPHEX_HTC(mu_c, mu_rat_c, Pr_c, k_c, G_c, info.C.Dh,  disp_flag);
                
            case 'S2P_recVap_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = S2P_recVap_BPHEX_HTC(mu_c, mu_rat_c, Pr_c, k_c, G_c, info.C.Dh,  disp_flag);
                
            case 'S2P_rec1p_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = S2P_rec1p_BPHEX_HTC(mu_c, mu_rat_c, Pr_c, k_c, G_c, info.C.Dh,  disp_flag);
                
            case 'S2P_rec1phaseOK_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = S2P_rec1phase_BPHEX_HTC(mu_c, mu_rat_c, Pr_c, k_c, G_c, info.C.Dh,  disp_flag);

            case 'S2P_evHTF_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = S2P_evHTF_BPHEX_HTC(mu_c, mu_rat_c, Pr_c, k_c, G_c, info.C.Dh, disp_flag);

            case 'Martin1_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = Martin1_BPHEX_HTC(mu_c, mu_rat_c, Pr_c, k_c, G_c, info.C.Dh, info.theta, disp_flag);
                
            case 'Martin2_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = Martin2_BPHEX_HTC(mu_c, Pr_c, k_c, G_c, info.C.Dh, info.theta, disp_flag);
                
            case 'Wanniarachchi_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = Wanniarachchi_BPHEX_HTC(mu_c, mu_rat_c, Pr_c, k_c, G_c, info.C.Dh, info.theta, info.phi, disp_flag);
                
            case 'Thonon_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = Thonon_BPHEX_HTC(mu_c, Pr_c, k_c, G_c, info.C.Dh, info.theta, disp_flag);
                
            case 'Junqi_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = Junqi_BPHEX_HTC(mu_c, Pr_c, k_c, G_c, info.C.Dh, info.theta, disp_flag);
                
            case 'Muley_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = Muley_BPHEX_HTC(mu_c, mu_rat_c, Pr_c, k_c, G_c, info.C.Dh, info.theta, info.phi, info.L_hex, disp_flag);
                
            case 'Kim_BPHEX'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = Kim_BPHEX_HTC(mu_c,Pr_c, k_c, G_c, info.C.Dh, info.theta, disp_flag);
                
            case 'DittusBoelter_Pipe'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = DittusBoelter_Pipe_HTC(mu_c, Pr_c, k_c, G_c, info.C.Dh, 0.4, disp_flag);
                
            case 'Gnielinski_Pipe'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = Gnielinski_Pipe_HTC(mu_c, Pr_c, k_c, G_c, info.C.Dh, info.L_hex, disp_flag);
                
            case 'VDI_finnedTubes_staggered'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = VDI_finnedTubes_staggered_HTC(mu_c, Pr_c, k_c, G_c, info.C.Dc, info.C.fin.omega_t, disp_flag);
                
            case 'Wang_finnedTubes_staggered'
                [hConv_1phase_c, Nu_1phase_c, flag_1phase_c] = Wang_finnedTubes_staggered_HTC(mu_c, Pr_c, k_c, G_c, info.C.Dh, info.C.Dc, info.C.Pt, info.C.Fp, info.C.Nr, disp_flag);
                
        end
        
        if 1
            %%%%
            if strcmp(out.C.type_zone{j}, 'liq')
                hConv_1phase_c = info.C_fit_1pl_c*hConv_1phase_c;
                Nu_1phase_c = info.C_fit_1pl_c*Nu_1phase_c;
            elseif strcmp(out.C.type_zone{j}, 'vap')
                hConv_1phase_c = info.C_fit_1pv_c*hConv_1phase_c;
                Nu_1phase_c = info.C_fit_1pv_c*Nu_1phase_c;
            end
            %%%%
        end

        
        if strcmp(out.C.type_zone{j}, 'liq') || strcmp(out.C.type_zone{j}, 'vap')
            out.C.hConv_vec(j) = hConv_1phase_c;
            out.C.Nu_vec(j) = Nu_1phase_c;
            out.C.fConv_vec(j) = flag_1phase_c;
        end
    end
    
    % Cold side : two-phase convective heat transfer coefficient
    if strcmp(out.C.type_zone{j}, 'tp') || strcmp(out.C.type_zone{j}, 'tp_dryout')
        G_c = (m_dot_c/info.C.n_canals)/info.C.CS;
        x_c = min(1, max(0, 0.5*out.C.x_vec(j) + 0.5*out.C.x_vec(j+1)));
        mu_c_l = CoolProp.PropsSI('V', 'Q', 0, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c); 
        k_c_l = CoolProp.PropsSI('L', 'Q', 0, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c);
        Pr_c_l = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c); 
        rho_c_l = CoolProp.PropsSI('D', 'Q', 0, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c); 
        rho_c_v = CoolProp.PropsSI('D', 'Q', 1, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c); 
        i_fg_c = CoolProp.PropsSI('H', 'Q', 1, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c) - CoolProp.PropsSI('H', 'Q', 0, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c);
        
        switch info.C.correlation.type_2phase_ev
            
            case 'UD'
                h_nom = info.C.h_nom_tp;
                m_dot_nom = info.C.m_dot_nom_tp;
                n_nom = info.C.n_nom_tp;
                [hConv_2phase_c, Nu_2phase_c, flag_2phase_c] = UD_HTC(m_dot_c, h_nom, m_dot_nom, n_nom);
                                
            case 'Han_boiling_BPHEX'               
                [hConv_2phase_c, Nu_2phase_c, flag_2phase_c] = Han_Boiling_BPHEX_HTC(min(x_c,x_di_c), mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v,  i_fg_c, G_c, out.DTlog(j)*out.F(j), out.Qdot_vec(j), out.H.hConv_vec(j), info.C.Dh, info.theta, info.pitch_co, disp_flag);
               
            case 'Amalfi_boiling_BPHEX'
                rho_c = CoolProp.PropsSI('D', 'Q', x_c, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c);
                sigma_c = CoolProp.PropsSI('I', 'Q', x_c, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c);
                mu_c_v = CoolProp.PropsSI('V', 'Q', 1, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c);
                [hConv_2phase_c, Nu_2phase_c, flag_2phase_c] = Amalfi_Boiling_BPHEX_HTC(min(x_c,x_di_c), rho_c_l, rho_c_v, rho_c, mu_c_v, mu_c_l, k_c_l, sigma_c, i_fg_c, G_c, info.C.Dh, info.theta, out.Qdot_vec(j), out.DTlog(j)*out.F(j), out.H.hConv_vec(j), disp_flag);
                
            case 'Junqi_Boiling_BPHEX'
                [hConv_2phase_c, Nu_2phase_c, flag_2phase_c] = Junqi_Boiling_BPHEX_HTC(min(x_c,x_di_c), mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v,  i_fg_c, G_c, out.DTlog(j)*out.F(j), out.Qdot_vec(j), out.H.hConv_vec(j), info.C.Dh, disp_flag);
        
            case 'Cooper_Boiling'
                P_c_crit = CoolProp.PropsSI('Pcrit', 'Q', 1, 'P', 1e5, fluid_c);
                M_c = 1e3*CoolProp.PropsSI('M', 'Q', 1, 'P', 1e5, fluid_c);
                [hConv_2phase_c, Nu_2phase_c, flag_2phase_c] = Cooper_Boiling_HTC(info.C.Dh, k_c_l, (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), P_c_crit, M_c, out.H.hConv_vec(j), out.DTlog(j)*out.F(j), out.Qdot_vec(j));

        end
        
        if 1
            %%%%
            hConv_2phase_c = info.C_fit_2p_c*hConv_2phase_c;
            Nu_2phase_c = info.C_fit_2p_c*Nu_2phase_c;
            %%%%
        end
        
        if strcmp(out.C.type_zone{j}, 'tp')
            out.C.hConv_vec(j) = hConv_2phase_c;
            out.C.Nu_vec(j) = Nu_2phase_c;
            out.C.fConv_vec(j) = flag_2phase_c;
        elseif strcmp(out.C.type_zone{j}, 'tp_dryout')
            out.C.hConv_vec(j) = hConv_2phase_c - (x_c-x_di_c)/(1-x_di_c)*(hConv_2phase_c - hConv_1phase_c);
            out.C.Nu_vec(j) = Nu_2phase_c - (x_c-x_di_c)/(1-x_di_c)*(Nu_2phase_c - Nu_1phase_c);
            out.C.fConv_vec(j) = min(flag_1phase_c,flag_2phase_c);
        end
        
        
    end
    
    % Hot side heat transfer efficiency (in case of fins)
    if strcmp(info.H.fin, 'none')
        out.H.eff_vec(j) = 1;
    else
        eta_eff = FinSchmidt(out.H.hConv_vec(j), info.H.fin.k, info.H.fin.th, info.H.fin.r, info.H.fin.B, info.H.fin.H);
        out.H.eff_vec(j) = 1-info.H.fin.omega_f*(1-eta_eff);
    end
    
    % Cold side heat transfer efficiency (in case of fins)
    if strcmp(info.C.fin, 'none')
        out.C.eff_vec(j) = 1;
    else
        eta_eff = FinSchmidt(out.C.hConv_vec(j), info.C.fin.k, info.C.fin.th, info.C.fin.r, info.C.fin.B, info.C.fin.H);
        out.C.eff_vec(j) = 1-info.C.fin.omega_f*(1-eta_eff);
    end
    
    % Global heat transfer coefficient and zone surface area
    out.k(j) = info.k0 + info.k1*(((out.H.T_vec(j+1)+out.H.T_vec(j)+out.C.T_vec(j)+out.C.T_vec(j+1))/4)-273.15);
    out.AU_vec(j) = (1/out.H.hConv_vec(j)/out.H.eff_vec(j)/info.H.A_tot + 1/out.C.hConv_vec(j)/out.C.eff_vec(j)/info.C.A_tot + info.C.R_fooling/info.C.A_tot + info.H.R_fooling/info.H.A_tot + info.t/out.k(j)/info.H.A_tot)^-1;
    out.U_vec(j) = out.AU_vec(j)/info.H.A_tot;
    out.H.A_vec(j) = out.Qdot_vec(j)/(out.DTlog(j)*out.F(j))/out.U_vec(j);
    out.C.A_vec(j) = out.H.A_vec(j)*info.C.A_tot/info.H.A_tot;
    
    % Cold-side dryout incipient quality for boiling processes
    if not(dry_out_c) && strcmp(out.C.type_zone{j}, 'tp') %only compute if dry out has not yet started
        switch info.C.correlation.dry_out_incipient
            case 'none'
                x_di_c = inf;
                
            case 'UD'
                x_di_c = info.C.x_di;
                
            case 'Kim_DOI'
                P_star_c = (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1))/CoolProp.PropsSI('Pcrit', 'Q', 1, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c);
                q_c = out.Qdot_vec(j)/out.C.A_vec(j);
                sigma_c_l = CoolProp.PropsSI('I', 'Q', 0, 'P', (0.5*out.C.P_vec(j)+0.5*out.C.P_vec(j+1)), fluid_c);
                x_di_c = Kim_DryOutIncipience(G_c, q_c, info.C.Dh, P_star_c, rho_c_l, rho_c_v, mu_c_l, sigma_c_l, i_fg_c);
        end
    end
    out.C.x_di(j) = x_di_c;
    
end

out.H.A_tot = sum(out.H.A_vec);
out.C.A_tot = sum(out.C.A_vec);
out.resA = 1 - out.H.A_tot/info.H.A_tot;

end
