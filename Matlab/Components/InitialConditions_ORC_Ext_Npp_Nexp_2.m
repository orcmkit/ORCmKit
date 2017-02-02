function IC = InitialConditions_ORC_Ext_Npp_Nexp_2(fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param)


if param.display
    fprintf('\n');
    dispstat('','init')
end

% Automatic Boundary Conditions
%if strcmp(param.solverType, 'M_imposed') && isfield(param, 'x0') && param.limit_P_pp_su
%    P_pp_su_lb = 0.95*param.x0(2);
%else
P_pp_su_lb = max(CoolProp.PropsSI('P', 'Q', 0, 'T', T_ctf_su-20, fluid_wf), CoolProp.PropsSI('P_min', 'Q', 0, 'T', 273.15, fluid_wf));
%end
P_pp_ex_ub = CoolProp.PropsSI('P', 'Q', 0, 'T', min(CoolProp.PropsSI('Tcrit', 'Q', 0, 'T',273, fluid_wf)-2, T_htf_su-1), fluid_wf);
rp_max = P_pp_ex_ub/P_pp_su_lb;
rp_min = min(1.01, rp_max);
P_pp_ex_lb = rp_min*P_pp_su_lb;
P_pp_su_ub = P_pp_ex_ub/rp_min;
h_ev_ex_lb = CoolProp.PropsSI('H', 'P', P_pp_ex_lb, 'Q', 0, fluid_wf);
h_ev_ex_ub = CoolProp.PropsSI('H', 'P', P_pp_ex_lb, 'T', T_htf_su, fluid_wf);

P_pp_su_guess0 = linspace(P_pp_su_lb,min(CoolProp.PropsSI('P', 'Q', 0, 'T', T_ctf_su+50, fluid_wf),P_pp_su_ub),param.init(1) );

%P_pp_su_guess0 = linspace(CoolProp.PropsSI('P', 'Q', 0, 'T', T_pp_su_lb, fluid_wf),CoolProp.PropsSI('P', 'Q', 0, 'T', T_ctf_su+25, fluid_wf),param.init(1) );
x_rp_guess0 =  linspace(0.2, 0.8, param.init(2));
x_h_ev_ex_guess0 =  linspace(0, 1, param.init(3));

if strcmp(param.solverType, 'M_imposed')
    T_pp_su_lb = T_ctf_su-2; %CoolProp.PropsSI('H', 'P', P_pp_su_lb, 'Q', 0, fluid_wf);
    %if strcmp(param.solverType, 'M_imposed') && isfield(param, 'x0') && param.limit_P_pp_su
    %    T_pp_su_ub = param.x0(4)+1;
    %else
    T_pp_su_ub = CoolProp.PropsSI('T', 'P', P_pp_su_ub, 'Q', 0, fluid_wf); %CoolProp.PropsSI('H', 'P', P_pp_su_ub, 'Q', 0, fluid_wf);
    %end
    x_h_pp_su_guess0 = linspace(0, 1, param.init(4));
    [res,P_pp_su_guess_vec, P_pp_su_lb_vec, P_pp_su_ub_vec, P_pp_ex_guess_vec, P_pp_ex_lb_vec, P_pp_ex_ub_vec, h_ev_ex_lb_vec, h_ev_ex_guess_vec, h_ev_ex_ub_vec,T_pp_su_lb_vec,T_pp_su_guess_vec,T_pp_su_ub_vec] = deal(NaN*ones(1,length(P_pp_su_guess0)*length(x_rp_guess0)*length(x_h_ev_ex_guess0)*length(x_h_pp_su_guess0)));
elseif strcmp(param.solverType, 'DTsc_imposed')
    [res,P_pp_su_guess_vec, P_pp_su_lb_vec, P_pp_su_ub_vec, P_pp_ex_guess_vec, P_pp_ex_lb_vec, P_pp_ex_ub_vec, h_ev_ex_lb_vec, h_ev_ex_guess_vec, h_ev_ex_ub_vec] = deal(NaN*ones(1,length(P_pp_su_guess0)*length(x_rp_guess0)*length(x_h_ev_ex_guess0)));
end
index = 0;
for i_P_pp_su = 1: length(P_pp_su_guess0)
    for i_rp = 1: length(x_rp_guess0)
        for i_hev_ex = 1: length(x_h_ev_ex_guess0)
            if strcmp(param.solverType, 'M_imposed')
                
                for i_hpp_su = 1: length(x_h_pp_su_guess0)
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
                    
                    T_pp_su_lb_vec(index)    = T_pp_su_lb;
                    %                     h_pp_su_guess_vec(index) = (1-x_h_pp_su_guess0(i_hpp_su))*CoolProp.PropsSI('H', 'P', P_pp_su_guess0(i_P_pp_su), 'T', T_ctf_su-1, fluid_wf) + x_h_pp_su_guess0(i_hpp_su)*CoolProp.PropsSI('H', 'P', P_pp_su_guess0(i_P_pp_su), 'T', CoolProp.PropsSI('T', 'P', P_pp_su_guess0(i_P_pp_su), 'Q', 0, fluid_wf)-1, fluid_wf);
                    T_pp_su_ub_vec(index) = T_pp_su_ub;
                    %                     h_pp_su_lb_vec(index)    =  T_ctf_su-10; %0.0001; %h_pp_su_lb;
                    
                    T_pp_su_guess_vec(index) = (1-x_h_pp_su_guess0(i_hpp_su))*T_pp_su_lb_vec(index) + x_h_pp_su_guess0(i_hpp_su)*(CoolProp.PropsSI('T', 'P', P_pp_su_guess0(i_P_pp_su), 'Q', 0, fluid_wf)-0.000001);
                    %T_pp_su_guess_vec(index) = (1-x_h_pp_su_guess0(i_hpp_su))*T_pp_su_lb_vec(index) + x_h_pp_su_guess0(i_hpp_su)*T_pp_su_ub;
                    
                    %h_pp_su_ub_vec(index) = CoolProp.PropsSI('T', 'P', P_pp_su_ub, 'Q', 0, fluid_wf); %h_pp_su_ub;
                    param.eval_type = 'fast';
                    param.EV.generateTS = 0;
                    param.CD.generateTS = 0;
                    param.REC.generateTS = 0;
                    lb_test = [0 0 0 0];
                    %                     ub_test = [P_pp_ex_guess_vec(index) P_pp_su_guess_vec(index) h_ev_ex_guess_vec(index) h_pp_su_guess_vec(index)];
                    ub_test = [P_pp_ex_guess_vec(index) P_pp_su_guess_vec(index) h_ev_ex_guess_vec(index) T_pp_su_guess_vec(index)];
                    x_test = [1 1 1 1];
                    
                    [guess, ~] = FCT_ORC_Ext_Npp_Nexp_2(x_test, lb_test, ub_test, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param);
                    
                    if any(guess.flag.value < 0)
                        res(index) = NaN;
                    else
                        res(index) = guess.res;
                    end
                end
                
            elseif strcmp(param.solverType, 'DTsc_imposed')
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
                
                [guess, ~] = FCT_ORC_Ext_Npp_Nexp_2(x_test, lb_test, ub_test, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param);
                
                if any(guess.flag.value < 0)
                    res(index) = NaN;
                else
                    res(index) = guess.res;
                end
            end
        end
    end
end

if strcmp(param.solverType, 'M_imposed') && isfield(param, 'x0_SatLiq')
    
    n_lin = 5;
    x_P_pp_su_x0 = [param.x0_SatLiq(2)+[0.5e4 1e4 3e4 5e4 7e4 1e5 1.5e5] ((1-linspace(0.01, 1,n_lin))*(param.x0_SatLiq(2)+1.5e5) + linspace(0.01, 1,n_lin)*min(param.x0_SatLiq(2)+5e5, P_pp_su_ub))];
    %x_P_pp_su_x0 =  [0.001 0.004 0.008 linspace(0.01, 0.8,10)];
    x_T_pp_su_x0 = linspace(0, 1, 4);
    [h_ev_ex_lb_vec_x0_SatLiq, h_ev_ex_ub_vec_x0_SatLiq, h_ev_ex_guess_vec_x0_SatLiq, P_pp_ex_lb_vec_x0_SatLiq, P_pp_ex_ub_vec_x0_SatLiq, P_pp_ex_guess_vec_x0_SatLiq, P_pp_su_lb_vec_x0_SatLiq, P_pp_su_ub_vec_x0_SatLiq, P_pp_su_guess_vec_x0_SatLiq, T_pp_su_lb_vec_x0_SatLiq, T_pp_su_ub_vec_x0_SatLiq, T_pp_su_guess_vec_x0_SatLiq, res_x0_SatLiq] = deal(NaN*ones(1,length(x_P_pp_su_x0)*length(x_T_pp_su_x0)));
    
    index = 0;
    for i_P_pp_su_x0 = 1: length(x_P_pp_su_x0)
        for i_T_pp_su_x0 = 1: length(x_T_pp_su_x0)
            
            index = index+1;
            if param.display
                dispstat(['x0 evaluation: ' num2str(index+length(res)) '/' num2str(length(res)+length(res_x0_SatLiq))])
            end
            
            h_ev_ex_lb_vec_x0_SatLiq(index) = h_ev_ex_lb;
            h_ev_ex_ub_vec_x0_SatLiq(index) = h_ev_ex_ub;
            h_ev_ex_guess_vec_x0_SatLiq(index) = max(h_ev_ex_lb+1, min(param.x0_SatLiq(3), h_ev_ex_ub-1));
            
            P_pp_ex_lb_vec_x0_SatLiq(index) = P_pp_ex_lb;
            P_pp_ex_ub_vec_x0_SatLiq(index) = P_pp_ex_ub;
            P_pp_ex_guess_vec_x0_SatLiq(index) = max(P_pp_ex_lb+1, min(param.x0_SatLiq(1), P_pp_ex_ub-1));
            
            P_pp_su_lb_vec_x0_SatLiq(index) = P_pp_su_lb;
            P_pp_su_ub_vec_x0_SatLiq(index) = P_pp_su_ub;
            P_pp_su_guess_vec_x0_SatLiq(index) = x_P_pp_su_x0(i_P_pp_su_x0); % (1-x_P_pp_su_x0(i_P_pp_su_x0))*(param.x0_SatLiq(2)+5e3) + x_P_pp_su_x0(i_P_pp_su_x0)*min(param.x0_SatLiq(2)+20e5, P_pp_su_ub);
            
            T_pp_su_lb_vec_x0_SatLiq(index)  = T_pp_su_lb;
            T_pp_su_ub_vec_x0_SatLiq(index)  = T_pp_su_ub;
            T_pp_su_guess_vec_x0_SatLiq(index)  = (1-x_T_pp_su_x0(i_T_pp_su_x0))*T_pp_su_lb + x_T_pp_su_x0(i_T_pp_su_x0)*param.x0_SatLiq(4);
            
            param.eval_type = 'fast';
            param.EV.generateTS = 0;
            param.CD.generateTS = 0;
            param.REC.generateTS = 0;
            lb_test = [0 0 0 0];
            ub_test = [P_pp_ex_guess_vec_x0_SatLiq(index) P_pp_su_guess_vec_x0_SatLiq(index) h_ev_ex_guess_vec_x0_SatLiq(index) T_pp_su_guess_vec_x0_SatLiq(index)];
            x_test = [1 1 1 1];
            
            [guess_x0_SatLiq, ~] = FCT_ORC_Ext_Npp_Nexp_2(x_test, lb_test, ub_test, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param);
            
            if any(guess_x0_SatLiq.flag.value < 0) || guess_x0_SatLiq.res > 1
                res_x0_SatLiq(index) = NaN;
            else
                res_x0_SatLiq(index) = guess_x0_SatLiq.res;
            end
        end
    end
    P_pp_su_guess_vec = [P_pp_su_guess_vec_x0_SatLiq P_pp_su_guess_vec];
    P_pp_su_lb_vec = [P_pp_su_lb_vec_x0_SatLiq P_pp_su_lb_vec];
    P_pp_su_ub_vec = [P_pp_su_ub_vec_x0_SatLiq P_pp_su_ub_vec];
    P_pp_ex_guess_vec = [P_pp_ex_guess_vec_x0_SatLiq P_pp_ex_guess_vec];
    P_pp_ex_lb_vec = [P_pp_ex_lb_vec_x0_SatLiq P_pp_ex_lb_vec];
    P_pp_ex_ub_vec = [P_pp_ex_ub_vec_x0_SatLiq P_pp_ex_ub_vec];
    h_ev_ex_guess_vec = [h_ev_ex_guess_vec_x0_SatLiq h_ev_ex_guess_vec];
    h_ev_ex_lb_vec = [h_ev_ex_lb_vec_x0_SatLiq h_ev_ex_lb_vec];
    h_ev_ex_ub_vec = [h_ev_ex_ub_vec_x0_SatLiq h_ev_ex_ub_vec];
    T_pp_su_guess_vec = [T_pp_su_guess_vec_x0_SatLiq T_pp_su_guess_vec];
    T_pp_su_lb_vec = [T_pp_su_lb_vec_x0_SatLiq T_pp_su_lb_vec];
    T_pp_su_ub_vec = [T_pp_su_ub_vec_x0_SatLiq T_pp_su_ub_vec];
    res = [res_x0_SatLiq , res];
end


%Predefined boundary conditions
if isfield(param, 'x0')
    x0_vec = 1;
    P_pp_su_lb_vec_x0 = P_pp_su_lb*ones(1,length(x0_vec));
    res_x0 = NaN*ones(1,length(x0_vec));
    P_pp_su_ub_vec_x0 = P_pp_su_ub*ones(1,length(x0_vec));
    P_pp_su_guess_vec_x0 = max(P_pp_ex_lb+1, min(param.x0(2), P_pp_ex_ub-1))*ones(1,length(x0_vec));
    P_pp_ex_lb_vec_x0 = P_pp_ex_lb*ones(1,length(x0_vec));
    P_pp_ex_ub_vec_x0 = P_pp_ex_ub*ones(1,length(x0_vec));
    P_pp_ex_guess_vec_x0 = max(P_pp_ex_lb+1, min(param.x0(1), P_pp_ex_ub-1))*ones(1,length(x0_vec));
    h_ev_ex_lb_vec_x0 = h_ev_ex_lb*ones(1,length(x0_vec));
    h_ev_ex_ub_vec_x0 = h_ev_ex_ub*ones(1,length(x0_vec));
    h_ev_ex_guess_vec_x0 = max(h_ev_ex_lb+1, min(param.x0(3), h_ev_ex_ub-1))*ones(1,length(x0_vec));
    if strcmp(param.solverType, 'DTsc_imposed')
        x0_matrix = [P_pp_ex_guess_vec_x0' P_pp_su_guess_vec_x0' h_ev_ex_guess_vec_x0'];
        ub0_matrix = [P_pp_ex_ub_vec_x0' P_pp_su_ub_vec_x0' h_ev_ex_ub_vec_x0'];
        lb0_matrix =[P_pp_ex_lb_vec_x0' P_pp_su_lb_vec_x0' h_ev_ex_lb_vec_x0'];
    elseif strcmp(param.solverType, 'M_imposed')
        T_pp_su_lb_vec_x0 = T_pp_su_lb*ones(1,length(x0_vec));
        T_pp_su_ub_vec_x0 = T_pp_su_ub*ones(1,length(x0_vec));
        T_pp_su_guess_vec_x0 = max(T_pp_su_lb+0.01, min(param.x0(4), T_pp_su_ub-0.01))*ones(1,length(x0_vec));
        x0_matrix = [P_pp_ex_guess_vec_x0' P_pp_su_guess_vec_x0' h_ev_ex_guess_vec_x0' T_pp_su_guess_vec_x0'];
        ub0_matrix = [P_pp_ex_ub_vec_x0' P_pp_su_ub_vec_x0' h_ev_ex_ub_vec_x0' T_pp_su_ub_vec_x0'];
        lb0_matrix =[P_pp_ex_lb_vec_x0' P_pp_su_lb_vec_x0' h_ev_ex_lb_vec_x0' T_pp_su_lb_vec_x0'];
    end
    index = 0;
    for k = 1:length(x0_vec)
        index = index+1;
        if param.display
            dispstat(['x0 evaluation: : ' num2str(index+length(res)) '/' num2str(length(res)+length(x0_vec))])
        end
        param.eval_type = 'fast';
        param.EV.generateTS = 0;
        param.CD.generateTS = 0;
        param.REC.generateTS = 0;
        out_x0 = FCT_ORC_Ext_Npp_Nexp_2(x0_matrix(k,:)./ub0_matrix(k,:), lb0_matrix(k,:), ub0_matrix(k,:), fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param);
        
        if any(out_x0.flag.value < 0) || out_x0.res > 1
            res_x0(k) = NaN;
        else
            res_x0(k) =  out_x0.res;
        end
    end
    
    % Global results
    P_pp_su_guess_vec = [P_pp_su_guess_vec_x0 P_pp_su_guess_vec];
    P_pp_su_lb_vec = [P_pp_su_lb_vec_x0 P_pp_su_lb_vec];
    P_pp_su_ub_vec = [P_pp_su_ub_vec_x0 P_pp_su_ub_vec];
    P_pp_ex_guess_vec = [P_pp_ex_guess_vec_x0 P_pp_ex_guess_vec];
    P_pp_ex_lb_vec = [P_pp_ex_lb_vec_x0 P_pp_ex_lb_vec];
    P_pp_ex_ub_vec = [P_pp_ex_ub_vec_x0 P_pp_ex_ub_vec];
    h_ev_ex_guess_vec = [h_ev_ex_guess_vec_x0 h_ev_ex_guess_vec];
    h_ev_ex_lb_vec = [h_ev_ex_lb_vec_x0 h_ev_ex_lb_vec];
    h_ev_ex_ub_vec = [h_ev_ex_ub_vec_x0 h_ev_ex_ub_vec];
    if strcmp(param.solverType, 'M_imposed')
        T_pp_su_guess_vec = [T_pp_su_guess_vec_x0 T_pp_su_guess_vec];
        T_pp_su_lb_vec = [T_pp_su_lb_vec_x0 T_pp_su_lb_vec];
        T_pp_su_ub_vec = [T_pp_su_ub_vec_x0 T_pp_su_ub_vec];
    end
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

if strcmp(param.solverType, 'M_imposed')
    IC.T_pp_su_guess_vec = T_pp_su_guess_vec(not(isnan(res)));
    IC.T_pp_su_lb_vec = T_pp_su_lb_vec(not(isnan(res)));
    IC.T_pp_su_ub_vec = T_pp_su_ub_vec(not(isnan(res)));
end

IC.res = res(not(isnan(res)));
