function IC = InitialConditions_ORC_Ext_Npp_Nexp_4(fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param)


if param.display
    fprintf('\n');
    dispstat('','init')
end

% Automatic Boundary Conditions

P_pp_su_lb = max(CoolProp.PropsSI('P', 'Q', 0, 'T', T_ctf_su-40, fluid_wf), CoolProp.PropsSI('P_min', 'Q', 0, 'T', 273.15, fluid_wf));
P_pp_ex_ub = CoolProp.PropsSI('P', 'Q', 0, 'T', min(CoolProp.PropsSI('Tcrit', 'Q', 0, 'T',273, fluid_wf)-2, T_htf_su-1), fluid_wf);
rp_max = P_pp_ex_ub/P_pp_su_lb;
rp_min = min(1.01, rp_max);
P_pp_ex_lb = rp_min*P_pp_su_lb;
P_pp_su_ub = P_pp_ex_ub/rp_min;
h_ev_ex_lb = CoolProp.PropsSI('H', 'P', P_pp_ex_lb, 'Q', 0, fluid_wf);
h_ev_ex_ub = CoolProp.PropsSI('H', 'P', P_pp_ex_lb, 'T', T_htf_su, fluid_wf);
P_pp_su_guess0 = linspace(P_pp_su_lb,min(CoolProp.PropsSI('P', 'Q', 0, 'T', T_ctf_su+50, fluid_wf),P_pp_su_ub),param.init(1) );
x_rp_guess0 =  linspace(0.1, 0.8, param.init(2));
x_h_ev_ex_guess0 =  linspace(0, 1, param.init(3));

[res,P_pp_su_guess_vec, P_pp_su_lb_vec, P_pp_su_ub_vec, P_pp_ex_guess_vec, P_pp_ex_lb_vec, P_pp_ex_ub_vec, h_ev_ex_lb_vec, h_ev_ex_guess_vec, h_ev_ex_ub_vec] = deal(NaN*ones(1,length(P_pp_su_guess0)*length(x_rp_guess0)*length(x_h_ev_ex_guess0)));
index = 0;
for i_P_pp_su = 1: length(P_pp_su_guess0)
    for i_rp = 1: length(x_rp_guess0)
        for i_hev_ex = 1: length(x_h_ev_ex_guess0)
            
            index = index+1;
            if param.display
                dispstat(['x0 evaluation: ' num2str(index) '/' num2str(length(res))])
            end
            
            P_pp_su_guess_vec(index) = P_pp_su_guess0(i_P_pp_su);
            P_pp_su_lb_vec(index) = P_pp_su_lb;
            P_pp_su_ub_vec(index) = P_pp_su_ub;
            
            P_pp_ex_guess_vec(index) = x_rp_guess0(i_rp)*P_pp_ex_ub + (1-x_rp_guess0(i_rp))*P_pp_su_guess_vec(index);
            P_pp_ex_lb_vec(index) = P_pp_ex_lb;
            P_pp_ex_ub_vec(index) = P_pp_ex_ub;
            
            h_ev_ex_lb_vec(index)    = h_ev_ex_lb;
            h_ev_ex_guess_vec(index) = (1-x_h_ev_ex_guess0(i_hev_ex))*CoolProp.PropsSI('H', 'P', P_pp_ex_guess_vec(index), 'Q', 0.5, fluid_wf) + x_h_ev_ex_guess0(i_hev_ex)*CoolProp.PropsSI('H', 'P', P_pp_ex_guess_vec(index), 'T', T_htf_su-2, fluid_wf);
            h_ev_ex_ub_vec(index) = h_ev_ex_ub;
            
            param.eval_type = 'fast';
            param.EV.generateTS = 0;
            param.CD.generateTS = 0;
            param.REC.generateTS = 0;
            lb_test = [0 0 0];
            ub_test = [P_pp_ex_guess_vec(index) P_pp_su_guess_vec(index) h_ev_ex_guess_vec(index)];
            x_test = [1 1 1];
            
            [guess, ~] = FCT_ORC_Ext_Npp_Nexp_4(x_test, lb_test, ub_test, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param);

            if any(guess.flag.value < 0)
                res(index) = NaN;
            else
                res(index) = guess.res;
            end
            
        end
    end
end

%Predefined boundary conditions
if isfield(param, 'x0')
        
    param.eval_type = 'fast';
    param.EV.generateTS = 0;
    param.CD.generateTS = 0;
    param.REC.generateTS = 0;
    
    out_x0 = FCT_ORC_Ext_Npp_Nexp_4(max(min(param.x0(1:3), [P_pp_ex_ub P_pp_su_ub h_ev_ex_ub]), [P_pp_ex_lb P_pp_su_lb h_ev_ex_lb])./[P_pp_ex_ub P_pp_su_ub h_ev_ex_ub], [P_pp_ex_lb P_pp_su_lb h_ev_ex_lb], [P_pp_ex_ub P_pp_su_ub h_ev_ex_ub], fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param);  
    
    if any(out_x0.flag.value < 0) 
        res_x0 = NaN;
    else
        res_x0 =  out_x0.res;
    end
        
    res_max = 1e-3;
    res_min = 5e-5;
    
    if param.restrict_bound || res_x0 < res_max
        delta_res = (res_x0-res_min)/(res_max-res_min);       
        delta_P_pp_su = (1-delta_res)*0.5e5 + delta_res*1e5;
        delta_P_pp_ex = (1-delta_res)*0.1e5 + delta_res*0.5e5;
        delta_h_ev_ex = (1-delta_res)*0.01*param.x0(3) + delta_res*0.05*param.x0(3);
        P_pp_su_lb_x0 = max(param.x0(2)-delta_P_pp_su,P_pp_su_lb);
        P_pp_su_ub_x0 = min(param.x0(2)+delta_P_pp_su,P_pp_su_ub);
        P_pp_su_guess_x0 = max(P_pp_ex_lb+1, min(param.x0(2), P_pp_ex_ub-1));       
        P_pp_ex_lb_x0 = max(param.x0(1)-delta_P_pp_ex,P_pp_ex_lb);
        P_pp_ex_ub_x0 = min(param.x0(1)+delta_P_pp_ex,P_pp_ex_ub);
        P_pp_ex_guess_x0 = max(P_pp_ex_lb+1, min(param.x0(1), P_pp_ex_ub-1));        
        h_ev_ex_lb_x0 = max(param.x0(3)-delta_h_ev_ex,h_ev_ex_lb);
        h_ev_ex_ub_x0 = min(param.x0(3)+delta_h_ev_ex,h_ev_ex_ub);
        h_ev_ex_guess_x0 = max(h_ev_ex_lb+1, min(param.x0(3), h_ev_ex_ub-1));       
    else
        P_pp_su_lb_x0 = P_pp_su_lb;
        P_pp_su_ub_x0 = P_pp_su_ub;
        P_pp_su_guess_x0 = max(P_pp_ex_lb+1, min(param.x0(2), P_pp_ex_ub-1));       
        P_pp_ex_lb_x0 = P_pp_ex_lb;
        P_pp_ex_ub_x0 = P_pp_ex_ub;
        P_pp_ex_guess_x0 = max(P_pp_ex_lb+1, min(param.x0(1), P_pp_ex_ub-1));        
        h_ev_ex_lb_x0 = h_ev_ex_lb;
        h_ev_ex_ub_x0 = h_ev_ex_ub;
        h_ev_ex_guess_x0 = max(h_ev_ex_lb+1, min(param.x0(3), h_ev_ex_ub-1));       
    end
    
    
    % Global results
    P_pp_su_guess_vec = [P_pp_su_guess_x0 P_pp_su_guess_vec];
    P_pp_su_lb_vec = [P_pp_su_lb_x0 P_pp_su_lb_vec];
    P_pp_su_ub_vec = [P_pp_su_ub_x0 P_pp_su_ub_vec];
    P_pp_ex_guess_vec = [P_pp_ex_guess_x0 P_pp_ex_guess_vec];
    P_pp_ex_lb_vec = [P_pp_ex_lb_x0 P_pp_ex_lb_vec];
    P_pp_ex_ub_vec = [P_pp_ex_ub_x0 P_pp_ex_ub_vec];
    h_ev_ex_guess_vec = [h_ev_ex_guess_x0 h_ev_ex_guess_vec];
    h_ev_ex_lb_vec = [h_ev_ex_lb_x0 h_ev_ex_lb_vec];
    h_ev_ex_ub_vec = [h_ev_ex_ub_x0 h_ev_ex_ub_vec];

    res = [res_x0 , res];
end

IC.P_pp_su_guess_vec = P_pp_su_guess_vec(not(isnan(res)));
IC.P_pp_su_lb_vec = P_pp_su_lb_vec(not(isnan(res)));
IC.P_pp_su_ub_vec = P_pp_su_ub_vec(not(isnan(res)));
IC.P_pp_ex_guess_vec = P_pp_ex_guess_vec(not(isnan(res)));
IC.P_pp_ex_lb_vec = P_pp_ex_lb_vec(not(isnan(res)));
IC.P_pp_ex_ub_vec = P_pp_ex_ub_vec(not(isnan(res)));
IC.h_ev_ex_guess_vec = h_ev_ex_guess_vec(not(isnan(res)));
IC.h_ev_ex_lb_vec = h_ev_ex_lb_vec(not(isnan(res)));
IC.h_ev_ex_ub_vec = h_ev_ex_ub_vec(not(isnan(res)));
IC.res = res(not(isnan(res)));

end