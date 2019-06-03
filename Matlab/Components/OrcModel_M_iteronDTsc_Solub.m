function [out_ORC, TS_ORC] = OrcModel_M_iteronDTsc_Solub(fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param)

if exist('Guesses_track.mat', 'file')
    delete('Guesses_track.mat')
end

tstart_ORC = tic;
param_bis = param;
param_bis.solverType = 'DTsc_imposed';
param_bis.restrict_bound = 0;

if param_bis.display
    fprintf('\n');
    disp('Start iteration:')
    fprintf('\n');
    fprintf('%-15s %-10s %-5s %-50s %-15s %-50s %-60s %-15s %-10s %-100s\n', 'DT_sc', '#', 'i0', 'x_in', 'res_in', 'x_out', 'res_out', 'res_M', 'flag_ORC', 'flag components');
    fprintf('\n');
    %param_bis.display_list = 1;
    %param_bis.display = 0;
end
f = @(x) ORC_M_res_Solub(x, fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param_bis);
lb = -0.05;
ub = 50;
if isfield(param_bis, 'x0') && length(param_bis.x0) == 4
    DTsc0 = param_bis.x0(4);
else
    DTsc0 = NaN;
end

DT_sc_ok =  zeroBrent_OrcIter_Solub2( DTsc0, lb, ub, 1e-16,1e-10, 5e-4, f ); 
param_bis.DT_sc = DT_sc_ok;

if exist('Guesses_track.mat', 'file')
    load('Guesses_track.mat', 'x0','DT_sc_prev')
    param_bis.x0 = x0;    
    param_bis.init = [1 1 1];
    param_bis.nbr_test = 1;
    if abs(DT_sc_prev-param_bis.DT_sc)<1
        param_bis.restrict_bound = 1;
    else
        param_bis.restrict_bound = 0;
    end
end

[out_ORC, TS_ORC] = OrcModel_DTscImposed_Solub(fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param_bis);

if out_ORC.flag_ORC > 0 && abs(out_ORC.res_ORC_M) > 1e-3
   out_ORC.flag_ORC = -4;
end

%if disp_init
    if out_ORC.flag_ORC > 0
        color = 1;
    else
        color = 2;
    end
    fprintf(color, '\n --> DTsc: %d  (°C) | res_M = %d (-)  | res_ORC = %d °C | flag_ORC = %d (-) \n',  param_bis.DT_sc, out_ORC.res_ORC_M, out_ORC.res, out_ORC.flag_ORC);
%end
out_ORC.time_ORC = toc(tstart_ORC);

if exist('Guesses_track.mat', 'file')
    delete('Guesses_track.mat')
end

end

function res = ORC_M_res_Solub(x, fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param)
param.DT_sc = x;
if exist('Guesses_track.mat', 'file')
    load('Guesses_track.mat', 'x0','DT_sc_prev')
    param.x0 = x0;
    param.init = [2 2 2];
    param.nbr_test = 3;
    if param.DT_sc >0 && abs(DT_sc_prev-param.DT_sc)<1
        param.restrict_bound = 1;
    elseif param.DT_sc < 0 && abs(DT_sc_prev-param.DT_sc)<1e-3
        param.restrict_bound = 1;
    else
        param.restrict_bound = 0;
    end
    
end

[out_DT_Sc_imposed, ~] = OrcModel_DTscImposed_Solub(fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param);

if out_DT_Sc_imposed.flag_ORC > 0
    res = out_DT_Sc_imposed.res_ORC_M;
    x0 = [out_DT_Sc_imposed.P_pp_ex  out_DT_Sc_imposed.P_lr_su out_DT_Sc_imposed.h_exp_su];
    DT_sc_prev = param.DT_sc;
    try
        save('Guesses_track.mat', 'x0','DT_sc_prev')
    catch
    end
    %if param.DT_sc == 0 && res < 0
    %    res = 0;
    %end
else
    res = NaN;
    return
%     if param.DT_sc == 0
%         if out_DT_Sc_imposed.flag_ORC > -3
%             x0 = [out_DT_Sc_imposed.P_pp_ex  out_DT_Sc_imposed.P_lr_su out_DT_Sc_imposed.h_exp_su];
%             DT_sc_prev = param.DT_sc;
%             try
%                 save('Guesses_track.mat', 'x0','DT_sc_prev')
%             catch
%             end
%         end
%         res = 0;
%     elseif param.DT_sc > 0
%         res = out_DT_Sc_imposed.res_ORC_M; %-1;
%     elseif param.DT_sc < 0
%         res = out_DT_Sc_imposed.res_ORC_M; %1;
%     end
end

end