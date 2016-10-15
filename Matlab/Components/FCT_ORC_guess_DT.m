function res  = FCT_ORC_guess_DT(dt_x, dt_lb, dt_ub, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param)

load('condition_initial.mat', 'IC')
dt_x = dt_x*dt_ub;
dt_x = min(max(dt_x, 0),dt_ub);

if dt_x < 0.1 && dt_x > 0
    dt_x = 0.1;
end

DT_guess = dt_x;
param.DT_sc =  DT_guess;
param.type_solver = 'DTsc_imposed';
k = 1;
out_ORC_best.res = 1e10;
stop = 0;
options_fmincon = optimset('Disp',param.displayIter,'Algorithm','interior-point','UseParallel',false,'TolX',1e-13,'TolFun',1e-13,'TolCon',1e-6,'MaxIter',1e3,'OutputFcn',@outputfunFS);


while not(stop) && k <= param.Nbr_comb
    
    x0 = [IC.P_pp_ex_guess_vec(k)	IC.P_pp_su_guess_vec(k)    IC.h_ev_ex_guess_vec(k)];
    ub = [IC.P_pp_ex_ub_vec(k)  	IC.P_pp_su_ub_vec(k)       IC.h_ev_ex_ub_vec(k)];
    lb = [IC.P_pp_ex_lb_vec(k)   	IC.P_pp_su_lb_vec(k)       IC.h_ev_ex_lb_vec(k)];
    A_ineq = [-1 1.001 0]; B_ineq = [0];
    param.eval_type = 'fast';
    param.EV.generateTS = 0;
    param.CD.generateTS = 0;
    param.REC.generateTS = 0;
    f = @(x) FCT_ORC_Ext_Npp_Nexp_res_2( x, lb, ub, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param);
    
%     if param.display
        fprintf('%-10s %-10s %-50s %-15s ', num2str(DT_guess), [num2str(k) '/' num2str(param.Nbr_comb)] , ['[' num2str(x0,'%15.4e') ']'] , num2str(f(x0./ub), '%.4g'));
%     end

    x = fmincon(f,x0./ub,A_ineq,B_ineq,[],[],lb./ub,ub./ub,[],options_fmincon);
    [out_ORC, ~] = FCT_ORC_Ext_Npp_Nexp_2(x, lb, ub, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, param);
    
    if any(out_ORC.flag.value<0)
        out_ORC.flag_ORC = - 1;
        stop = 0;
    elseif all(out_ORC.flag.value>0) && any(abs(out_ORC.res_vec) > 5e-5)
        out_ORC.flag_ORC = -2;
        stop = 0;
    else
        out_ORC.flag_ORC = 1;
        stop = 1;
    end
    
%     if param.display
        fprintf('%-50s %-60s %-15s %-10s %-100s \n', ['[' num2str(x.*ub,'%15.4e') ']'], [ num2str(out_ORC.res, '%.4g') '  [ ' num2str(out_ORC.res_vec ,'%15.4e') ' ] '] , num2str(out_ORC.res_ORC_M ,'%15.4e'), num2str(out_ORC.flag_ORC), num2str(out_ORC.flag.value));
%     end
    if out_ORC.res < out_ORC_best.res
        out_ORC_best = out_ORC;
    end
    out_ORC = out_ORC_best;
    k = k+1;
end
if out_ORC.flag_ORC > 0
    res = out_ORC.res_ORC_M;    
    IC.P_pp_ex_guess_vec = [out_ORC.P_recc_su   IC.P_pp_ex_guess_vec];
    IC.P_pp_su_guess_vec = [out_ORC.P_pp_su     IC.P_pp_su_guess_vec];
    IC.h_ev_ex_guess_vec = [out_ORC.h_dphp_su   IC.h_ev_ex_guess_vec];
    IC.P_pp_ex_ub_vec = [IC.P_pp_ex_ub_vec(1) IC.P_pp_ex_ub_vec];
    IC.P_pp_su_ub_vec = [IC.P_pp_su_ub_vec(1) IC.P_pp_su_ub_vec];
    IC.h_ev_ex_ub_vec = [IC.h_ev_ex_ub_vec(1) IC.h_ev_ex_ub_vec];
    IC.P_pp_ex_lb_vec = [IC.P_pp_ex_lb_vec(1) IC.P_pp_ex_lb_vec];
    IC.P_pp_su_lb_vec = [IC.P_pp_su_lb_vec(1) IC.P_pp_su_lb_vec];
    IC.h_ev_ex_lb_vec = [IC.h_ev_ex_lb_vec(1) IC.h_ev_ex_lb_vec];
    save('condition_initial.mat', 'IC')
else
    res = NaN;
end
end

function stop = outputfunFS(x, optimValues, state)
stop = norm(optimValues.fval) < 5e-5;
end

