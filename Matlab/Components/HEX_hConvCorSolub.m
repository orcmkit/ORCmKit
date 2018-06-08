function out = HEX_hConvCorSolub(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info, h_h_l, h_h_v, h_c_l, h_c_v)
out = HEX_profile_Solub(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info, h_h_l, h_h_v, h_c_l, h_c_v);
[out.H.A_vec, out.C.A_vec, out.H.hConv_vec, out.C.hConv_vec, out.H.Nu_vec, out.C.Nu_vec, out.H.fConv_vec, out.C.fConv_vec,out.DTlog, out.H.eff_vec, out.C.eff_vec, out.AU_vec, out.U_vec, out.k] = deal(NaN*ones(1,length(out.H.H_vec)-1));
x_di_c = 1; dry_out_c = 0;
x_di_h = 1;
disp_flag = 0;
for j = 1:length(out.H.T_vec)-1

    % LMTD for the current cell
    out.DTlog(j) = deltaT_log(out.H.T_vec(j+1), out.H.T_vec(j),out.C.T_vec(j), out.C.T_vec(j+1));
    
    T_wall = (out.H.T_vec(j+1)+ out.H.T_vec(j)+out.C.T_vec(j)+ out.C.T_vec(j+1))/4;

    % What type of cells for hot side (1phase/2phase)
    if strcmp(info.H.type, 'H')
        if (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)) < h_h_l
            out.H.type_zone{j} = 'liq';
            T_wall_h = T_wall;
        elseif (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)) > h_h_v
            out.H.type_zone{j} = 'vap';
            T_wall_h = max(T_wall, out.H.Tsat_pure_vec(j)+5e-2);
        else
            if (0.5*out.H.x_vec(j)+0.5*out.H.x_vec(j+1)) <= x_di_h
                out.H.type_zone{j} = 'tp';
                T_wall_h = T_wall;
            else
                out.H.type_zone{j} = 'tp_dryout';
                T_wall_h = max(T_wall, out.H.Tsat_pure_vec(j)+5e-2);
            end
        end
    elseif strcmp(info.H.type, 'T')
        out.H.type_zone{j} = 'liq';
        T_wall_h = T_wall;
    end
    
    % What type of cells for cold side (1phase/2phase)
    if strcmp(info.C.type, 'H')
        if (0.5*out.C.H_vec(j)+0.5*out.C.H_vec(j+1)) < h_c_l
            out.C.type_zone{j} = 'liq';
            T_wall_c = min(T_wall, out.C.Tsat_pure_vec(j)-5e-2);
        elseif (0.5*out.C.H_vec(j)+0.5*out.C.H_vec(j+1)) > h_c_v
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
        
    
    if 0
        % Hot-side convective heat transfer coefficient
        if strcmp(out.H.type_zone{j}, 'liq') || strcmp(out.H.type_zone{j}, 'vap')
            
            % Compute single-phase properties
            G_h = m_dot_h/info.H.n_canals/info.H.CS;
            if info.H.solub
                mu_h = CoolProp.PropsSI('V',        'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h); %to be updated with mixture properties
                Pr_h = CoolProp.PropsSI('Prandtl',  'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h); %to be updated with mixture properties
                k_h  = CoolProp.PropsSI('L',        'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h); %to be updated with mixture properties
            else
                if strcmp(fluid_h(1:3), 'ICP')
                    k_h =  PropsSI_ICP('L', 'T', 0.5*out.H.T_vec(j) + 0.5*out.H.T_vec(j+1), 'P', P_h_su, fluid_h);
                    mu_h = PropsSI_ICP('V', 'T', 0.5*out.H.T_vec(j) + 0.5*out.H.T_vec(j+1), 'P', P_h_su, fluid_h);
                    cp_h = PropsSI_ICP('C', 'T', 0.5*out.H.T_vec(j) + 0.5*out.H.T_vec(j+1), 'P', P_h_su, fluid_h);
                    Pr_h = cp_h*mu_h/k_h;
                else
                    mu_h = CoolProp.PropsSI('V',        'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
                    Pr_h = CoolProp.PropsSI('Prandtl',  'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
                    k_h  = CoolProp.PropsSI('L',        'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
                end
            end
            if strcmp(out.H.type_zone{j}, 'liq')
                type_correlation_h = info.H.correlation.type_1phase_l;
            elseif strcmp(out.H.type_zone{j}, 'vap')
                type_correlation_h = info.H.correlation.type_1phase_v;
            end
            switch type_correlation_h
                
                case 'Martin1_BPHEX'
                    out.H.hConv_vec(j) = Martin1_BPHEX_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, disp_flag);
                    
                case 'Martin2_BPHEX'
                    out.H.hConv_vec(j) = Martin2_BPHEX_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, disp_flag);
                    
                case 'Wanniarachchi_BPHEX'
                    out.H.hConv_vec(j) = Wanniarachchi_BPHEX_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, info.phi, disp_flag);
                    
                case 'Thonon_BPHEX'
                    out.H.hConv_vec(j) = Thonon_BPHEX_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, disp_flag);
                    
                case 'Junqi_BPHEX'
                    out.H.hConv_vec(j) = Junqi_BPHEX_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, disp_flag);
                    
                case 'Muley_BPHEX'
                    out.H.hConv_vec(j) = Muley_BPHEX_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.theta, info.L_hex, info.phi, disp_flag);
                    
                case 'GnielinskiSha_pipe'
                    out.H.hConv_vec(j) = GnielinskiSha_pipe_HTC(mu_h, Pr_h, k_h, G_h, info.H.Dh, info.H.Lt, disp_flag);
            end
            
        elseif strcmp(out.H.type_zone{j}, 'tp')
            
            %         % Compute two-phase properties
            %         mu_h_l = CoolProp.PropsSI('V', 'Q', 0, 'P', P_h_su, fluid_h);
            %         k_h_l = CoolProp.PropsSI('L', 'Q', 0, 'P', P_h_su, fluid_h);
            %         x_h = CoolProp.PropsSI('Q', 'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
            %         Pr_h_l = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_h_su, fluid_h);
            %         rho_h_l = CoolProp.PropsSI('D', 'Q', 0, 'P', P_h_su, fluid_h);
            %         rho_h_v = CoolProp.PropsSI('D', 'Q', 1, 'P', P_h_su, fluid_h);
            %         G_h = (m_dot_h/info.H.n_canals)/info.H.CS;
            %
            %         switch info.H.correlation.type_2phase_cd
            %             case 'Han_cond_BPHEX'
            %                 out.H.hConv_vec(j) = Han_Cond_BPHEX_HTC(x, mu_l, k_l, Pr_l, rho_l, rho_v, G, Dh, pitch_co, theta);
            %
            %             case 'Longo_condensation'
            %                 mu_h_l = CoolProp.PropsSI('V', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 k_h_l = CoolProp.PropsSI('L', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 Pr_h_l = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 x_h = CoolProp.PropsSI('Q', 'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
            %                 rho_h_l = CoolProp.PropsSI('D', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 rho_h_v = CoolProp.PropsSI('D', 'Q', 1, 'P', P_h_su, fluid_h);
            %                 G_h = (m_dot_h/info.H.n_canals)/info.H.CS;
            %                 G_h_eq = G_h * ( (1 - x_h) + x_h * (rho_h_l/rho_h_v)^0.5);
            %                 Re_h_eq = G_h_eq*info.H.Dh/mu_h_l;
            %                 i_fg_h = CoolProp.PropsSI('H', 'Q', 1, 'P', P_h_su, fluid_h) - CoolProp.PropsSI('H', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 g = 9.81;
            %                 if Re_h_eq < 1600
            %                     T_sat = (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1));
            %                     T_wall = (0.25*out.H.T_vec(j)+0.25*out.H.T_vec(j+1) + 0.25*out.C.T_vec(j)+0.25*out.C.T_vec(j+1));
            %                     out.H.hConv_vec(j) = info.H.fact_corr_2p*info.phi*0.943*((k_h_l^3*rho_h_l^2*g*i_fg_h)/(mu_h_l*(T_sat-T_wall)*info.L_hex))^0.25;
            %                 else
            %                     out.H.hConv_vec(j) = info.H.fact_corr_2p*1.875*info.phi*k_h_l/info.H.Dh*Re_h_eq^(info.H.fact2_corr_2p*0.445)*Pr_h_l^0.3333333;
            %                 end
            %
            %             case 'Cavallini_condensation'
            %                 mu_h_l = CoolProp.PropsSI('V', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 mu_h_v = CoolProp.PropsSI('V', 'Q', 1, 'P', P_h_su, fluid_h);
            %                 rho_h_l = CoolProp.PropsSI('D', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 rho_h_v = CoolProp.PropsSI('D', 'Q', 1, 'P', P_h_su, fluid_h);
            %                 Pr_h_l = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 x_h = CoolProp.PropsSI('Q', 'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
            %                 k_h_l = CoolProp.PropsSI('L', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 d_i = info.H.Dh;
            %                 g = 9.81;
            %                 C_T = 2.6; %1.6 for HC or 2.6 for other refrigerant
            %                 X_tt = ((mu_h_l/mu_h_v)^0.1)*((rho_h_v/rho_h_l)^0.5)*((1-x_h)/x_h)^0.9; %Martinelli factor
            %
            %                 G_h = (m_dot_h/info.H.n_canals)/info.H.CS;
            %                 Re_h_l = G_h*info.H.Dh/mu_h_l;
            %                 J_v = x_h*G_h/sqrt(g*d_i*rho_h_v*(rho_h_l-rho_h_v));
            %                 J_v_T = (((7.5/(4.3*X_tt^1.111 + 1))^-3) + ((C_T)^-3) )^-0.333333333333333333;
            %                 h_h_lo = 0.023*Re_h_l^(info.H.fact2_corr_2p*0.8)*Pr_h_l^0.4*k_h_l/d_i;
            %                 h_h_a = h_h_lo*(1 + (1.128*x_h^0.817)*((rho_h_l/rho_h_v)^0.3685)*((mu_h_l/mu_h_v)^0.2363)*((1-mu_h_v/mu_h_l)^2.144)*(Pr_h_l^-0.1));
            %
            %                 if J_v > J_v_T %delta_T-independent flow regime
            %                     out.H.hConv_vec(j) = info.H.fact_corr_2p*h_h_a;
            %                 elseif J_v <= J_v_T %delta_T-dependent flow regime
            %                     i_fg_h = CoolProp.PropsSI('H', 'Q', 1, 'P', P_h_su, fluid_h) - CoolProp.PropsSI('H', 'Q', 0, 'P', P_h_su, fluid_h);
            %                     T_sat = (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1));
            %                     T_wall = (0.25*out.H.T_vec(j)+0.25*out.H.T_vec(j+1) + 0.25*out.C.T_vec(j)+0.25*out.C.T_vec(j+1));
            %                     h_strat = 0.725*((1+0.741*((1-x_h)/x_h)^0.3321)^-1)*(((k_h_l^3*rho_h_l*(rho_h_l-rho_h_v)*g*i_fg_h)/(mu_h_l*d_i*(T_sat-T_wall)))^0.25)+((1-x_h^0.087)*h_h_lo);
            %                     h_h_d = J_v/J_v_T*(h_h_a*(J_v_T/J_v)^0.8 - h_strat) + h_strat;
            %                     out.H.hConv_vec(j) = info.H.fact_corr_2p*h_h_d;
            %                 end
            %
            %             case 'Shah_condensation'
            %                 mu_h_l = CoolProp.PropsSI('V', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 p_h_star = P_c_su/CoolProp.PropsSI('Pcrit', 'Q', 1, 'P', P_c_su, fluid_c);
            %                 x_h = CoolProp.PropsSI('Q', 'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
            %                 k_h_l = CoolProp.PropsSI('L', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 G_h = (m_dot_h/info.H.n_canals)/info.H.CS;
            %                 Re_h_l = G_h*info.H.Dh/mu_h_l;
            %                 Pr_h_l = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_h_su, fluid_h);
            %                 out.H.hConv_vec(j) = info.H.fact_corr_2p*0.023*(k_h_l/info.H.Dh)*(Re_h_l^(info.H.fact2_corr_2p*0.8))*(Pr_h_l^0.4)*(((1-x_h)^0.8)+ ((3.8*(x_h^0.76)*(1-x_h)^0.04)/(p_h_star^0.38)));
            %
            %         end
            
        end
    end
    
    % Hot-side : single phase convective heat transfer coefficient
    if strcmp(out.H.type_zone{j}, 'liq') || strcmp(out.H.type_zone{j}, 'vap') || strcmp(out.H.type_zone{j}, 'tp_dryout')
        G_h = m_dot_h/info.H.n_canals/info.H.CS;
        if info.H.solub
            if (strcmp(out.H.type_zone{j}, 'liq')) && ((0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)) >= out.H.Tsat_pure_vec(j)-5e-2)
                mu_h = CoolProp.PropsSI('V',        'T', (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)), 'Q', 0, fluid_h); %to be updated with mixture properties
                Pr_h = CoolProp.PropsSI('Prandtl',  'T', (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)), 'Q', 0, fluid_h); %to be updated with mixture properties
                k_h  = CoolProp.PropsSI('L',        'T', (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)), 'Q', 0, fluid_h); %to be updated with mixture properties  
                mu_rat_h = 1;
            elseif (strcmp(out.H.type_zone{j}, 'vap') || strcmp(out.H.type_zone{j}, 'tp_dryout') ) && ((0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)) <= out.H.Tsat_pure_vec(j)+5e-2)
                mu_h = CoolProp.PropsSI('V',        'T', (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)), 'Q', 1, fluid_h); %to be updated with mixture properties
                Pr_h = CoolProp.PropsSI('Prandtl',  'T', (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)), 'Q', 1, fluid_h); %to be updated with mixture properties
                k_h  = CoolProp.PropsSI('L',        'T', (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)), 'Q', 1, fluid_h); %to be updated with mixture properties
                mu_rat_h = 1;
            else
                mu_h = CoolProp.PropsSI('V',        'T', (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)), 'P', P_h_su, fluid_h); %to be updated with mixture properties
                Pr_h = CoolProp.PropsSI('Prandtl',  'T', (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)), 'P', P_h_su, fluid_h); %to be updated with mixture properties
                k_h  = CoolProp.PropsSI('L',        'T', (0.5*out.H.T_vec(j)+0.5*out.H.T_vec(j+1)), 'P', P_h_su, fluid_h); %to be updated with mixture properties  
                mu_wall_h = CoolProp.PropsSI('V',  	'T', T_wall_h ,  'P', P_h_su, fluid_h); %to be updated with mixture properties
                mu_rat_h = mu_h/mu_wall_h;
            end
        else
            if strcmp(fluid_h(1:3), 'ICP')
                k_h =  PropsSI_ICP('L', 'T', 0.5*out.H.T_vec(j) + 0.5*out.H.T_vec(j+1), 'P', P_h_su, fluid_h);
                mu_h = PropsSI_ICP('V', 'T', 0.5*out.H.T_vec(j) + 0.5*out.H.T_vec(j+1), 'P', P_h_su, fluid_h);
                cp_h = PropsSI_ICP('C', 'T', 0.5*out.H.T_vec(j) + 0.5*out.H.T_vec(j+1), 'P', P_h_su, fluid_h);
                Pr_h = cp_h*mu_h/k_h;
                mu_wall_h = PropsSI_ICP('V', 'T', T_wall_h, 'P', P_h_su, fluid_h);
                mu_rat_h = mu_h/mu_wall_h;
            else
                mu_h = CoolProp.PropsSI('V',        'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
                Pr_h = CoolProp.PropsSI('Prandtl',  'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
                k_h  = CoolProp.PropsSI('L',        'H', (0.5*out.H.H_vec(j)+0.5*out.H.H_vec(j+1)), 'P', P_h_su, fluid_h);
                mu_wall_h = CoolProp.PropsSI('V',  	'T', T_wall_h ,  'P', P_h_su, fluid_h); 
                mu_rat_h = mu_h/mu_wall_h;
            end
        end
        if strcmp(out.H.type_zone{j}, 'liq')
            type_correlation_h = info.H.correlation.type_1phase_l;
        elseif strcmp(out.H.type_zone{j}, 'vap') || strcmp(out.H.type_zone{j}, 'tp_dryout')
            type_correlation_h = info.H.correlation.type_1phase_v;
        end
            
        switch type_correlation_h
            
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
        end
        
        if strcmp(out.H.type_zone{j}, 'liq') || strcmp(out.H.type_zone{j}, 'vap')
            out.H.hConv_vec(j) = hConv_1phase_h;
            out.H.Nu_vec(j) = Nu_1phase_h;
            out.H.fConv_vec(j) = flag_1phase_h;           
        end
    end
    
    % Cold-side : single phase convective heat transfer coefficient
    if strcmp(out.C.type_zone{j}, 'liq') || strcmp(out.C.type_zone{j}, 'vap') || strcmp(out.C.type_zone{j}, 'tp_dryout')
        G_c = m_dot_c/info.C.n_canals/info.C.CS;
        if info.C.solub
            if (strcmp(out.C.type_zone{j}, 'liq')) && ((0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)) >= out.C.Tsat_pure_vec(j)-5e-2)
                mu_c = CoolProp.PropsSI('V',        'T', (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)), 'Q', 0, fluid_c); %to be updated with mixture properties
                Pr_c = CoolProp.PropsSI('Prandtl',  'T', (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)), 'Q', 0, fluid_c); %to be updated with mixture properties
                k_c  = CoolProp.PropsSI('L',        'T', (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)), 'Q', 0, fluid_c); %to be updated with mixture properties  
                mu_rat_c = 1;
            elseif (strcmp(out.C.type_zone{j}, 'vap') || strcmp(out.C.type_zone{j}, 'tp_dryout') ) && ((0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)) <= out.C.Tsat_pure_vec(j)+5e-2)
                mu_c = CoolProp.PropsSI('V',        'T', (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)), 'Q', 1, fluid_c); %to be updated with mixture properties
                Pr_c = CoolProp.PropsSI('Prandtl',  'T', (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)), 'Q', 1, fluid_c); %to be updated with mixture properties
                k_c  = CoolProp.PropsSI('L',        'T', (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)), 'Q', 1, fluid_c); %to be updated with mixture properties
                mu_rat_c = 1;
            else
                mu_c = CoolProp.PropsSI('V',        'T', (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)), 'P', P_c_su, fluid_c); %to be updated with mixture properties
                Pr_c = CoolProp.PropsSI('Prandtl',  'T', (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)), 'P', P_c_su, fluid_c); %to be updated with mixture properties
                k_c  = CoolProp.PropsSI('L',        'T', (0.5*out.C.T_vec(j)+0.5*out.C.T_vec(j+1)), 'P', P_c_su, fluid_c); %to be updated with mixture properties  
                mu_wall_c = CoolProp.PropsSI('V',  	'T', T_wall_c ,  'P', P_c_su, fluid_c); %to be updated with mixture properties
                mu_rat_c = mu_c/mu_wall_c;
            end
        else
            if strcmp(fluid_c(1:3), 'ICP')
                k_c =  PropsSI_ICP('L', 'T', 0.5*out.C.T_vec(j) + 0.5*out.C.T_vec(j+1), 'P', P_c_su, fluid_c);
                mu_c = PropsSI_ICP('V', 'T', 0.5*out.C.T_vec(j) + 0.5*out.C.T_vec(j+1), 'P', P_c_su, fluid_c);
                cp_c = PropsSI_ICP('C', 'T', 0.5*out.C.T_vec(j) + 0.5*out.C.T_vec(j+1), 'P', P_c_su, fluid_c);
                Pr_c = cp_c*mu_c/k_c;
                mu_wall_c = PropsSI_ICP('V', 'T', T_wall_c, 'P', P_c_su, fluid_c);
                mu_rat_c = mu_c/mu_wall_c;
            else
                mu_c = CoolProp.PropsSI('V',        'H', (0.5*out.C.H_vec(j)+0.5*out.C.H_vec(j+1)), 'P', P_c_su, fluid_c);
                Pr_c = CoolProp.PropsSI('Prandtl',  'H', (0.5*out.C.H_vec(j)+0.5*out.C.H_vec(j+1)), 'P', P_c_su, fluid_c);
                k_c  = CoolProp.PropsSI('L',        'H', (0.5*out.C.H_vec(j)+0.5*out.C.H_vec(j+1)), 'P', P_c_su, fluid_c);
                mu_wall_c = CoolProp.PropsSI('V',  	'T', T_wall_c ,  'P', P_c_su, fluid_c); 
                mu_rat_c = mu_c/mu_wall_c;
            end
        end
        if strcmp(out.C.type_zone{j}, 'liq')
            type_correlation_c = info.C.correlation.type_1phase_l;
        elseif strcmp(out.C.type_zone{j}, 'vap') || strcmp(out.C.type_zone{j}, 'tp_dryout')
            type_correlation_c = info.C.correlation.type_1phase_v;
        end
            
        switch type_correlation_c
            
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
        end
        
        if strcmp(out.C.type_zone{j}, 'liq') || strcmp(out.C.type_zone{j}, 'vap')
            out.C.hConv_vec(j) = hConv_1phase_c;
            out.C.Nu_vec(j) = Nu_1phase_c;
            out.C.fConv_vec(j) = flag_1phase_c;           
        end
    end
    
    % Cold side : two-phase convective heat transfer coefficient
    if strcmp(out.C.type_zone{j}, 'tp') || strcmp(out.C.type_zone{j}, 'tp_dryout')
                
        switch info.C.correlation.type_2phase_ev
            case 'Han_boiling'
                
                G_c = (m_dot_c/info.C.n_canals)/info.C.CS;
                if info.C.solub % lubricant-refrigerant mixture
                    x_c = 0.5*out.C.x_vec(j) + 0.5*out.C.x_vec(j+1);
                    mu_c_l = CoolProp.PropsSI('V', 'Q', 0, 'P', P_c_su, fluid_c); %to be updated with mixture properties
                    k_c_l = CoolProp.PropsSI('L', 'Q', 0, 'P', P_c_su, fluid_c); %to be updated with mixture properties
                    Pr_c_l = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_c_su, fluid_c); %to be updated with mixture properties
                    rho_c_l = CoolProp.PropsSI('D', 'Q', 0, 'P', P_c_su, fluid_c); %to be updated with mixture properties
                    rho_c_v = CoolProp.PropsSI('D', 'Q', 1, 'P', P_c_su, fluid_c); %to be updated with super heated vapor properties
                    i_fg_c = CoolProp.PropsSI('H', 'Q', 1, 'P', P_c_su, fluid_c) - CoolProp.PropsSI('H', 'Q', 0, 'P', P_c_su, fluid_c);                    
                else % pure working fluid
                    x_c = CoolProp.PropsSI('Q', 'H', (0.5*out.C.H_vec(j)+0.5*out.C.H_vec(j+1)), 'P', P_c_su, fluid_c);
                    mu_c_l = CoolProp.PropsSI('V', 'Q', 0, 'P', P_c_su, fluid_c);
                    k_c_l = CoolProp.PropsSI('L', 'Q', 0, 'P', P_c_su, fluid_c);
                    Pr_c_l = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_c_su, fluid_c);
                    rho_c_l = CoolProp.PropsSI('D', 'Q', 0, 'P', P_c_su, fluid_c);
                    rho_c_v = CoolProp.PropsSI('D', 'Q', 1, 'P', P_c_su, fluid_c);
                    i_fg_c = CoolProp.PropsSI('H', 'Q', 1, 'P', P_c_su, fluid_c) - CoolProp.PropsSI('H', 'Q', 0, 'P', P_c_su, fluid_c);                    
                end
                
                [hConv_2phase_c, Nu_2phase_c, flag_2phase_c] = Han_Boiling_BPHEX_HTC(min(x_c,x_di_c), mu_c_l, k_c_l, Pr_c_l, rho_c_l, rho_c_v,  i_fg_c, G_c, out.DTlog(j), out.Qdot_vec(j), out.H.hConv_vec(j), info.C.Dh, info.theta, info.pitch_co, disp_flag);

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
    out.U_vec(j) = out.AU_vec(j)/info.H.A_tot;%(1/out.H.hConv_vec(j)/out.H.eff_vec(j) + 1/out.C.hConv_vec(j)/out.C.eff_vec(j)/(info.C.A_tot/info.H.A_tot))^-1;
    out.H.A_vec(j) = out.Qdot_vec(j)/out.DTlog(j)/out.U_vec(j);
    out.C.A_vec(j) = out.H.A_vec(j)*info.C.A_tot/info.H.A_tot;
    
    % Cold-side dryout incipient quality for boiling processes
    if not(dry_out_c) && strcmp(out.C.type_zone{j}, 'tp') %only compute if dry out has not yet started
        switch info.C.correlation.dry_out_incipient
            case 'none'
                x_di_c = inf;
                
            case 'UD'
                x_di_c = info.C.x_di;
                
            case 'Kim_DOI'
                P_star_c = P_c_su/CoolProp.PropsSI('Pcrit', 'Q', 1, 'P', P_c_su, fluid_c);
                q_c = out.Qdot_vec(j)/out.C.A_vec(j);
                sigma_c_l = CoolProp.PropsSI('I', 'Q', 0, 'P', P_c_su, fluid_c);
                x_di_c = Kim_DryOutIncipience(G_c, q_c, info.C.Dh, P_star_c, rho_c_l, rho_c_v, mu_c_l, sigma_c_l, i_fg_c);
        end
    end
    
    out.C.x_di(j) = x_di_c;
end

out.H.A_tot = sum(out.H.A_vec);
out.C.A_tot = sum(out.C.A_vec);
out.resA = 1 - out.H.A_tot/info.H.A_tot;

end
