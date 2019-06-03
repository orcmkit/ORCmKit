function [out_ORC, TS_ORC] = OrcModel_DTscImposed_Solub(fluid_wf, fluid_lub, fluid_htf,  h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param)
%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 28/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% OrcModel is a single matlab function developed to calculate the steady-state
% equilibrium conditions of an ORC (see the Documentation/HexModel_MatlabDoc)
%
% The model inputs are:
%       - fluid_wf: nature of the working fluid                     [-]
%       - fluid_htf: nature of the hot source fluid                 [-]
%       - in_htf_su: inlet temperature or enthalpy of the hot source[K or J/kg]
%       - P_htf_su: inlet pressure of the hot source                [Pa]
%       - m_dot_htf: mass flow rate of the hot source               [kg/s]
%       - fluid_ctf: nature of the hot source fluid                 [-]
%       - in_ctf_su: inlet temperature or enthalpy of the hot source[K or J/kg]
%       - P_ctf_su: inlet pressure of the hot source                [Pa]
%       - m_dot_ctf: mass flow rate of the hot source               [kg/s]
%       - T_amb : ambient temperature                               [K]
%       - N_exp : expander speed                                    [rpm]
%       - N_pp  : pump speed                                        [rpm]
%       - param: structure variable containing the model parameters
%
% Please refer to the docucmentation of HexModel, PumpModel, ExpanderModel
% and LossesModel for further details about the component models
%
% The model outputs are:
%       - out: a structure variable which includes at miniumum the following information:
%
%       - TS : a stucture variable which contains the vectors of temperature
%              and entropy of the fluid (useful to generate a Ts diagram
%              when modelling the entire ORC system
%
% See the documentation for further details or contact rdickes@ulg.ac.be

%% DEMONSTRATION CASE
if nargin == 0
    
    fluid_wf = 'R245fa';
    fluid_htf = 'PiroblocBasic';
    P_htf_su = 2e5;
    in_htf_su = 379.4284;
    m_dot_htf = 0.9294;
    fluid_ctf = 'air';
    P_ctf_su = 1e5;
    in_ctf_su = 4.1725e+05;
    m_dot_ctf = 1.4136;
    T_amb = 291.0039;
    N_exp = 5000;
    N_pp = 200;
    param.solverType = 'M_imposed';
    param.DT_sc = 9.0787;
    param.x_cd_ex = 0;
    param.M_tot = 25;
    
    path = 'C:\Users\RDickes\Google Drive\PhD\MOR study\ORC\Experimental database\Sun2Power';
    
    EV_folder = [path '\Evaporator\'];
    load([EV_folder, 'ParametersCalibration_EV.mat'])
    param.EV = EV_hConvVar;
    param.EV.V_h_tot = 0.009;
    param.EV.V_c_tot = 0.009;
    param.EV.displayResults = 0;
    param.EV.displayTS = 0;
    
    CD_folder = [path '\Condenser\'];
    load([CD_folder, 'ParametersCalibration_CD.mat'])
    param.CD = CD_hConvVar;
    param.CD.V_h_tot = 0.014;
    param.CD.V_c_tot = 0.7585;
    param.CD.displayResults = 0;
    param.CD.displayTS = 0;
    param.CD.W_dot_aux = 0;
    
    REC_folder = [path '\Recuperator\'];
    load([REC_folder, 'ParametersCalibration_REC.mat'])
    param.REC = REC_hConvVar;
    param.REC.V_h_tot = 0.001026;
    param.REC.V_c_tot = 0.00108;
    param.REC.displayResults = 0;
    param.REC.displayTS = 0;
    
    PP_folder = [path '\Pump\'];
    load([PP_folder, 'ParametersCalibration_PP.mat'])
    param.PP = PP_SemiEmp;
    param.PP.V = 1.4e-3;
    param.PP.displayResults = 0;
    param.PP.displayTS = 0;
    
    EXP_folder = [path '\Expander\'];
    load([EXP_folder, 'ParametersCalibration_EXP.mat'])
    param.EXP = EXP_SemiEmp;
    param.EXP.modelType = 'SemiEmp';
    param.EXP.V = 1.4e-3;
    param.EXP.displayResults = 0;
    param.EXP.displayTS = 0;
    load('C:\Users\RDickes\Google Drive\PhD\MOR study\ORC\Experimental database\Sun2Power\OffDesign\gamma_R245fa.mat');
    param.EXP.gamma.gamma_PQ_pol = gamma_PQ_R245fa; param.EXP.gamma.gamma_PT_pol = gamma_PT_R245fa;
    
    DP_folder = [path '\PressureDrops\'];
    load([DP_folder, 'ParametersCalibration_DP.mat'])
    param.LossesHP = DPHP_PhiDP;
    param.LossesHP.displayResults = 0;
    param.LossesHP.displayTS = 0;
    
    param.LossesLP = DPLP_PhiDP;
    param.LossesLP.displayResults = 0;
    param.LossesLP.displayTS = 0;
    
    param.V_aux_pp_ex =  2.21035e-4;
    param.V_aux_recc_ex = 9.54259e-5;
    param.V_aux_ev_ex = 9.54259e-5;
    param.V_aux_exp_ex =  7.08995e-4;
    param.V_aux_rech_ex = 7.01469e-4;
    param.V_liq_rec = 5.7e-3;
    param.V_aux_cd_ex = 5.991508e-3 - param.V_liq_rec;
    param.displayTS = 1;
    param.displayResults =0;
    
    param.init = [ 2 2 2 2];
    param.nbr_test = 1;
    param.display = 1;
end

tstart_ORC = tic;


%% ORC MODELING
if not(isfield(param, 'display_list'))
    param.display_list = param.display;
end

% Initial conditions evaluation
IC = InitialConditions_ORC_Ext_Npp_Nexp_3_Solub(fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param);
if max(param.init) <2
    if length(IC.res) < 1 %|| min(IC.res) > 2
        out_ORC.flag_ORC = - 3;
        out_ORC.res = NaN;
        out_ORC.res_ORC_M = NaN;
        out_ORC.res_vec = NaN;
        TS_ORC = NaN;
        return
    end
end
 
% Order guesses results
[~, j_order] = sort(IC.res);

%tuned now
if length(IC.res) > param.nbr_test
    nbr_sort_kept = min(2,length(j_order)-1);
    vec_sort_kept = 1:nbr_sort_kept;
    vec_sort_mixt = (nbr_sort_kept+1):length(j_order);
    j_order = [j_order(vec_sort_kept) j_order(vec_sort_mixt(randperm(length(vec_sort_mixt))) )];
end

res_ordered = IC.res(j_order);

Nbr_comb_x0 = length(j_order);
Nbr_comb_x0_max = param.nbr_test;
if param.display
    disp('x0 residuals and index:')
    fprintf('\n');
    disp(num2str([j_order(1:min(Nbr_comb_x0,Nbr_comb_x0_max)); res_ordered(1:min(Nbr_comb_x0,Nbr_comb_x0_max))]))
    fprintf('\n');
end

% Start evaluation
if not(isempty(res_ordered))
    
    if param.display
        fprintf('\n');
        disp('Start iteration:')
        fprintf('\n');
        fprintf('%-15s %-10s %-5s %-50s %-15s %-50s %-60s %-15s %-10s %-100s\n', 'DT_sc', '#', 'i0', 'x_in', 'res_in', 'x_out', 'res_out', 'res_M', 'flag_ORC', 'flag components');
        fprintf('\n');
    end
    k = 1;
    out_ORC_best.res = 1e10;
    stop = 0;
    options_fmincon = optimset('Disp',param.displayIter,'Algorithm','interior-point','UseParallel',false,'TolX',1e-13,'TolFun',1e-13,'TolCon',1e-6,'MaxIter',1e3,'OutputFcn',@outputfunFS);
    %options_pattsearch = psoptimset('Display','iter','TolX', 1e-8, 'TolFun', 1e-8, 'TolMesh', 1e-8, 'MaxIter', 1e4, 'MaxFunEvals', 1e8, 'OutputFcns',@outputfunPS);

    while not(stop) && k <= min(Nbr_comb_x0,Nbr_comb_x0_max);
        
        x0 = [IC.P_pp_ex_guess_vec(j_order(k))    IC.P_pp_su_guess_vec(j_order(k))    IC.h_ev_ex_guess_vec(j_order(k))];
        ub = [IC.P_pp_ex_ub_vec(j_order(k))       IC.P_pp_su_ub_vec(j_order(k))       IC.h_ev_ex_ub_vec(j_order(k))];
        lb = [IC.P_pp_ex_lb_vec(j_order(k))       IC.P_pp_su_lb_vec(j_order(k))       IC.h_ev_ex_lb_vec(j_order(k))];
        A_ineq = [-1 1.1 0]; B_ineq = [0];
        if param.display_list
            fprintf('%-15s %-10s %-5d %-50s %-15s ', num2str(param.DT_sc,'%15.5e'), [num2str(k) '/' num2str(min(Nbr_comb_x0,Nbr_comb_x0_max))] , j_order(k), ['[' num2str(x0,'%15.4e') ']'] , num2str(IC.res(j_order(k)), '%.4g'));
        end
        param.eval_complete = 0;
        param.EV.generateTS = 0;
        param.CD.generateTS = 0;
        param.REC.generateTS = 0;
        param.EXP.generateTS = 0;
        f = @(x) FCT_ORC_Ext_Npp_Nexp_res_3_Solub( x, lb, ub, fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param);
        x = fmincon(f,x0./ub,A_ineq,B_ineq,[],[],lb./ub,ub./ub,[],options_fmincon);
        %[x, ~, ~,  ] = patternsearch(f,x0./ub,[],[],[],[],lb./ub,ub./ub,[],options_pattsearch);
        param.eval_complete = 1;
        param.EV.generateTS = 0;
        param.CD.generateTS = 0;
        param.REC.generateTS = 0;
        param.EXP.generateTS = 0;
        
        [out_ORC, TS_ORC] = FCT_ORC_Ext_Npp_Nexp_3_Solub(x, lb, ub, fluid_wf, fluid_lub, fluid_htf, h_htf_su, P_htf_su, m_dot_htf, fluid_ctf, h_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, C_oil, param);
        
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
        if param.display_list
            fprintf('%-50s %-60s %-15s %-10s %-100s \n', ['[' num2str(x.*ub,'%15.4e') ']'], [ num2str(out_ORC.res, '%.4g') '  [ ' num2str(out_ORC.res_vec,'%15.4e') ' ] '], num2str(out_ORC.res_ORC_M,'%15.5e'), num2str(out_ORC.flag_ORC), num2str(out_ORC.flag.value));
        end
        if out_ORC.res < out_ORC_best.res
            out_ORC_best = out_ORC;
        end
        out_ORC = out_ORC_best;
        k = k+1;
    end
else
    out_ORC.flag_ORC = - 1;
    TS_ORC = NaN;
end
out_ORC.time_ORC = toc(tstart_ORC);
if param.display
    fprintf('\n')
    dispstat('','keepprev')
end
out_ORC = orderfields(out_ORC);

if param.displayResults ==1 && out_ORC.flag_ORC == - 1
    dispstat('','keepprev')
    dispstat('Error: The model did not converge correctly','keepthis');
end


end

function [stop,options,optchanged] = outputfunPS(optimvalues,options,flag)
stop = optimvalues.fval < 1e-5;
optchanged = 0;
end

function stop = outputfunFS(x, optimValues, state)
%disp(norm(optimValues.fval))
stop = norm(optimValues.fval) < 5e-5;
end