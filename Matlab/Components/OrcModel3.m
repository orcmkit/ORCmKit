function [out_ORC, TS_ORC] = OrcModel3(fluid_wf, fluid_htf, in_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC)
%% DEMONSTRATION CASE
if nargin == 0
    clear all
    close all
    clc
    fluid_wf = 'R245fa';
    fluid_htf = 'water';
    P_htf_su = 6.17e5;
    in_htf_su = CoolProp.PropsSI('H', 'T', 417.6812, 'P', P_htf_su, fluid_htf);
    m_dot_htf = 1.0221;
    fluid_ctf = 'INCOMP::MPG-30%';
    P_ctf_su = 2.61e5;
    in_ctf_su = CoolProp.PropsSI('H', 'T', 312.73, 'P', P_ctf_su, fluid_ctf);
    m_dot_ctf = 1.41;
    T_amb = 20+273.15;
    N_exp = 3050;
    N_pp = 555;
    ORC.DT_sc = 7.9032;
    ORC.PP.modelType = 'CstEff';
    ORC.PRE.modelType = 'CstEff';
    ORC.EV.modelType = 'CstEff';
    ORC.EXP.modelType = 'CstEffVs';
    ORC.CD.modelType = 'CstEff';
    ORC.SUB.modelType = 'CstEff';
    ORC.REC.modelType = 'CstEff';
    ORC.DPHP.modelType = 'MdotDP';
    ORC.DPLP.modelType = 'MdotDP';
    switch ORC.PP.modelType
        case 'CstEff'
            ORC.PP.epsilon_is = 0.8325;
            ORC.PP.epsilon_vol = 0.9270;
            ORC.PP.V_s = 0.056304612/1000;
            ORC.PP.displayResults = 0;
            ORC.PP.advancedUser = 1;
        case 'PolEff'
            ORC.PP.coeffPol_is = [ -0.1364    0.4455   -1.3517   -0.2913    6.3817  -37.8517];
            ORC.PP.coeffPol_vol = [0.9889   -0.0287    0.1525   -0.0002    0.0358   -0.3096];
            ORC.PP.V_s = 0.056304612/1000;
            ORC.PP.displayResults = 0;
            ORC.PP.advancedUser = 1;
        case 'SemiEmp'
            ORC.PP.A_leak = 6.7701e-07;
            ORC.PP.W_dot_0_loss = 133.8142;
            ORC.PP.K_0_loss = 0.9176;
            ORC.PP.V_s = 0.056304612/1000;
            ORC.PP.displayResults = 0;
            ORC.PP.advancedUser = 1;
    end
    switch ORC.EXP.modelType
        case 'CstEffVs'
            ORC.EXP.V_s = 1.3720e-04;
            ORC.EXP.epsilon_s = 0.3968;
            ORC.EXP.FF = 0.9369;
            ORC.EXP.AU_amb = 6.35;
            ORC.EXP.displayResults = 0;
            ORC.EXP.advancedUser = 1;
        case 'Pol2Eff'
            ORC.EXP.V_s = 1.3720e-04;
            ORC.EXP.aPole_s = [  -1.357122817986534   0.812789686165796   0.006294724441691  -0.090801870614303    -0.001277224069567  -0.000011659489633];
            ORC.EXP.aPoleFF = [   0.316346130984361   0.334002213989119   0.003314515707531  -0.065050071542210    0.000928564423544  -0.000042899215077];
            ORC.EXP.AU_amb = 6.35;
            ORC.EXP.displayResults = 0;
            ORC.EXP.advancedUser = 1;
        case 'SemEmp_rv_mc_lk_dp_ht'
            ORC.EXP.V_s = 1.3720e-04;
            ORC.EXP.r_v_in = r_v_in;
            ORC.EXP.A_leak0 = 1.0153e-08;
            ORC.EXP.A_leak_n = 0;
            ORC.EXP.d_su = 0.0131;
            ORC.EXP.alpha = 0.3321;
            ORC.EXP.W_dot_loss_0 = 2.2772e+03;
            ORC.EXP.C_loss = 0;
            ORC.EXP.AU_su_n = 200.1089;
            ORC.EXP.AU_ex_n = 319.1532;
            ORC.EXP.AU_amb = 7.2414;
            ORC.EXP.M_dot_n = 0.6187;
            ORC.EXP.P_su_n = 1.7901e+06;
            ORC.EXP.displayResults = 0;
            ORC.EXP.advancedUser = 1;
    end
    switch ORC.PRE.modelType
        case 'CstEff'
            ORC.PRE.epsilon_th = 0.7069;
            ORC.PRE.type_h = 'H';
            ORC.PRE.type_c = 'H';
            ORC.PRE.displayResults = 0;
            ORC.PRE.displayTS = 0;
            ORC.PRE.advancedUser = 1;
    end
    switch ORC.EV.modelType
        case 'CstEff'
            ORC.EV.epsilon_th = 0.8815;
            ORC.EV.type_h = 'H';
            ORC.EV.type_c = 'H';
            ORC.EV.displayResults = 0;
            ORC.EV.displayTS = 0;
            ORC.EV.advancedUser = 1;
    end
    switch ORC.CD.modelType
        case 'CstEff'
            ORC.CD.epsilon_th = 0.9809;
            ORC.CD.type_h = 'H';
            ORC.CD.type_c = 'H';
            ORC.CD.displayResults = 0;
            ORC.CD.displayTS = 0;
            ORC.CD.advancedUser = 1;
    end
    switch ORC.SUB.modelType
        case 'CstEff'
            ORC.SUB.epsilon_th = 0.7453;
            ORC.SUB.type_h = 'H';
            ORC.SUB.type_c = 'H';
            ORC.SUB.displayResults = 0;
            ORC.SUB.displayTS = 0;
            ORC.SUB.advancedUser = 1;
    end
    switch ORC.REC.modelType
        case 'CstEff'
            ORC.REC.epsilon_th = 0.3996;
            ORC.REC.type_h = 'H';
            ORC.REC.type_c = 'H';
            ORC.REC.displayResults = 0;
            ORC.REC.displayTS = 0;
            ORC.REC.advancedUser = 1;
    end
    switch ORC.DPHP.modelType
        case 'CstDP'
            ORC.DPHP.dp = 0;
            ORC.DPHP.displayResults = 0;
            ORC.DPHP.displayTS = 0;
            ORC.DPHP.advancedUser = 1;
        case 'MdotDP'
            ORC.DPHP.Coeff_dp = [-9.070465145624621e+03 4.382373256958956e+05];
            ORC.DPHP.displayResults = 0;
            ORC.DPHP.displayTS = 0;
            ORC.DPHP.advancedUser = 1;
    end
    switch ORC.DPLP.modelType
        case 'CstDP'
            ORC.DPLP.dp = 0;
            ORC.DPLP.displayResults = 0;
            ORC.DPLP.displayTS = 0;
            ORC.DPLP.advancedUser = 1;
        case 'MdotDP'
            ORC.DPLP.Coeff_dp = [-4.042406334398474e+03 2.032055246425624e+05];
            ORC.DPLP.displayResults = 0;
            ORC.DPLP.displayTS = 0;
            ORC.DPLP.advancedUser = 1;
    end    
    ORC.displayTS = 1;
    ORC.displayResults =0;
end

tstart_ORC = tic;

if strcmp(ORC.CD.type_c, 'T')
    T_ctf_su = in_ctf_su;
else
    T_ctf_su = CoolProp.PropsSI('T', 'H', in_ctf_su, 'P', P_ctf_su, fluid_ctf);
end
if strcmp(ORC.EV.type_h, 'T')
    T_htf_su = in_htf_su;
else
    T_htf_su = CoolProp.PropsSI('T', 'H', in_htf_su, 'P', P_htf_su, fluid_htf);
end
%% INITIAL CONDITIONS
% 1) Automatic intial conditions
fprintf('\n');
dispstat('','init')
P_cd_lb = max(CoolProp.PropsSI('P', 'Q', 0, 'T', T_ctf_su-10, fluid_wf), CoolProp.PropsSI('P_min', 'Q', 0, 'T', 273.15, fluid_wf));
P_ev_ub = CoolProp.PropsSI('P', 'Q', 0, 'T', min(CoolProp.PropsSI('Tcrit', 'Q', 0, 'T',273, fluid_wf)-2, T_htf_su-1), fluid_wf);
rp_max = P_ev_ub/P_cd_lb;
rp_min = min(1.01, rp_max);
P_ev_lb = rp_min*P_cd_lb;
P_cd_ub = P_ev_ub/rp_min;
Q_dot_rec_lb = 0;
P_cd_guess0 = linspace(CoolProp.PropsSI('P', 'Q', 0, 'T', T_ctf_su, fluid_wf),CoolProp.PropsSI('P', 'Q', 0, 'T', T_ctf_su+30, fluid_wf) ,4 );
x_rp_guess0 =  linspace(0.1, 0.9, 4);
x_Q_dot_rec_guess0 = linspace(0.01, 0.9, 4);
index = 0;
[res,P_cd_guess_vec, P_cd_lb_vec, P_cd_ub_vec, P_ev_guess_vec, P_ev_lb_vec, P_ev_ub_vec, Q_dot_rec_guess_vec, Q_dot_rec_lb_vec, Q_dot_rec_max_vec] = deal(NaN*ones(1,length(P_cd_guess0)*length(x_rp_guess0)*length(x_Q_dot_rec_guess0)));
for i_pcd = 1: length(P_cd_guess0)
   for i_rp = 1: length(x_rp_guess0)
       for i_qdot_rec = 1: length(x_Q_dot_rec_guess0)

                   index = index+1;
                   dispstat(['x0 evaluation: ' num2str(index) '/' num2str(length(P_cd_guess_vec))])
                   P_cd_guess_vec(index) = P_cd_guess0(i_pcd);
                   P_cd_lb_vec(index) = P_cd_lb;
                   P_cd_ub_vec(index) = P_cd_ub;
                   
                   P_ev_guess_vec(index) = x_rp_guess0(i_rp)*P_ev_ub + (1-x_rp_guess0(i_rp))*P_cd_guess_vec(index);
                   P_ev_lb_vec(index) = P_ev_lb;
                   P_ev_ub_vec(index) = P_ev_ub;
                
                   z = x_Q_dot_rec_guess0(i_qdot_rec);
                   guess = OrganicRankineCycle_init(P_cd_guess_vec(index), P_ev_guess_vec(index), z, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC);
                   Q_dot_rec_lb_vec(index) = Q_dot_rec_lb;
                   Q_dot_rec_guess_vec(index) = guess.Q_dot_rec_guess;
                   Q_dot_rec_max_vec(index) = guess.Q_dot_rec_max;
                   %dispstat(num2str(guess.flag.value), 'keepthis','keepprev')
                   %guess.flag.name
                   if any(guess.flag.value < 0) 
                       res(index) = NaN;
                   else
                       res(index) = guess.res;
                   end                          
       end
   end
end
Q_dot_rec_ub = 10*max(Q_dot_rec_max_vec);
Q_dot_rec_ub_vec = Q_dot_rec_ub*ones(1, length(Q_dot_rec_guess_vec));
fprintf('\n');

% 2) Guess intial conditions based on domain
if isfield(ORC, 'x0_data') 
    
    %f5lin = @(x,xdata) x(1)*xdata(:,1) + x(2)*xdata(:,2) + x(3)*xdata(:,3) + x(4)*xdata(:,4) + x(5)*xdata(:,5) + x(6);
    f6lin = @(x,xdata) x(1)*xdata(:,1) + x(2)*xdata(:,2) + x(3)*xdata(:,3) + x(4)*xdata(:,4) + x(5)*xdata(:,5) + x(6)*xdata(:,6) + x(7);

    options = optimset('Display','off');
    P_cd_fit = lsqcurvefit(f6lin,[1 1 1 1 1 1 1],[ORC.x0_data.m_dot_ctf', ORC.x0_data.h_ctf_su'./1e5, ORC.x0_data.m_dot_htf',ORC.x0_data.T_htf_su'./1e2, ORC.x0_data.N_pp'./1e2, ORC.x0_data.N_exp'./1e3],ORC.x0_data.P_pp_su'./1e5, [], [], options);
    P_cd_guess = f6lin(P_cd_fit, [m_dot_ctf, in_ctf_su/1e5, m_dot_htf, in_htf_su/1e2, N_pp/1e2, N_exp/1e3])*1e5;
    
    P_ev_fit = lsqcurvefit(f6lin,[1 1 1 1 1 1 1],[ORC.x0_data.m_dot_ctf', ORC.x0_data.h_ctf_su'./1e5, ORC.x0_data.m_dot_htf',ORC.x0_data.T_htf_su'./1e2, ORC.x0_data.N_pp'./1e2, ORC.x0_data.N_exp'./1e3],ORC.x0_data.P_recc_su'./1e5, [], [], options);
    P_ev_guess = f6lin(P_ev_fit, [m_dot_ctf, in_ctf_su/1e5, m_dot_htf, in_htf_su/1e2, N_pp/1e2, N_exp/1e3])*1e5;
    
    Q_dot_rec_fit = lsqcurvefit(f6lin,[1 1 1 1 1 1 1],[ORC.x0_data.m_dot_ctf', ORC.x0_data.h_ctf_su'./1e5, ORC.x0_data.m_dot_htf',ORC.x0_data.T_htf_su'./1e2, ORC.x0_data.N_pp'./1e2, ORC.x0_data.N_exp'./1e3], ORC.x0_data.Q_dot_rec'./1e3, [], [], options);
    Q_dot_rec_guess = f6lin(Q_dot_rec_fit, [m_dot_ctf, in_ctf_su/1e5, m_dot_htf, in_htf_su/1e2, N_pp/1e2, N_exp/1e3])*1e3;
    
    x0_vec = [1 0.9 1.1 0.8 1.2 0.6 1.4 ];
    
    P_cd_guess_vec_x0 = P_cd_guess*x0_vec;
    P_cd_lb_vec_x0 = P_cd_lb*ones(1,length(x0_vec));
    P_cd_ub_vec_x0 = P_cd_ub*ones(1,length(x0_vec));
    
    P_ev_guess_vec_x0 = P_ev_guess*x0_vec;
    P_ev_lb_vec_x0 = P_ev_lb*ones(1,length(x0_vec));
    P_ev_ub_vec_x0 = P_ev_ub*ones(1,length(x0_vec));
    
    Q_dot_rec_guess_vec_x0 = Q_dot_rec_guess*x0_vec;
    Q_dot_rec_lb_vec_x0 = Q_dot_rec_lb*ones(1,length(x0_vec));
    Q_dot_rec_ub_vec_x0 = Q_dot_rec_ub*ones(1,length(x0_vec));
    
    res_x0_data = NaN*ones(1,length(x0_vec));
    for k = 1:length(x0_vec)
        index = index+1;
        dispstat(['x0 evaluation: ' num2str(index) '/' num2str(length(P_cd_guess_vec)+length(x0_vec))])
        out_x0 = OrganicRankineCycle([P_ev_guess_vec_x0(k) P_cd_guess_vec_x0(k) Q_dot_rec_guess_vec_x0(k)]./[P_ev_ub_vec_x0(k) P_cd_ub_vec_x0(k) Q_dot_rec_ub_vec_x0(k)], [P_ev_lb_vec_x0(k) P_cd_lb_vec_x0(k) Q_dot_rec_lb_vec_x0(k)], [P_ev_ub_vec_x0(k) P_cd_ub_vec_x0(k) Q_dot_rec_ub_vec_x0(k)], fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC);     
        if any(out_x0.flag.value < 0)
            res_x0_data(k) = NaN;
        else
            res_x0_data(k) = out_x0.res;
        end      
    end                   
    P_cd_guess_vec = [P_cd_guess_vec_x0 P_cd_guess_vec];
    P_cd_lb_vec = [P_cd_lb_vec_x0 P_cd_lb_vec];
    P_cd_ub_vec = [P_cd_ub_vec_x0 P_cd_ub_vec];    
    P_ev_guess_vec = [P_ev_guess_vec_x0 P_ev_guess_vec];
    P_ev_lb_vec = [P_ev_lb_vec_x0 P_ev_lb_vec];
    P_ev_ub_vec = [P_ev_ub_vec_x0 P_ev_ub_vec];  
    Q_dot_rec_guess_vec = [Q_dot_rec_guess_vec_x0 Q_dot_rec_guess_vec];
    Q_dot_rec_lb_vec = [Q_dot_rec_lb_vec_x0 Q_dot_rec_lb_vec];
    Q_dot_rec_ub_vec = [Q_dot_rec_ub_vec_x0 Q_dot_rec_ub_vec];
    res = [res_x0_data , res];
end


% 3) Guess intial conditions based on a single point
if isfield(ORC, 'x0') 
    x0_vec = [1 1.2 0.6 1.4 0.7 1.3];
    P_cd_guess_vec_x0 = ORC.x0(2)*x0_vec;
    P_cd_lb_vec_x0 = P_cd_lb*ones(1,length(x0_vec));
    P_cd_ub_vec_x0 = P_cd_ub*ones(1,length(x0_vec));    
    P_ev_guess_vec_x0 = ORC.x0(1)*x0_vec;
    P_ev_lb_vec_x0 = P_ev_lb*ones(1,length(x0_vec));
    P_ev_ub_vec_x0 = P_ev_ub*ones(1,length(x0_vec));    
    Q_dot_rec_guess_vec_x0 = ORC.x0(3)*x0_vec;
    Q_dot_rec_lb_vec_x0 = Q_dot_rec_lb*ones(1,length(x0_vec));
    Q_dot_rec_ub_vec_x0 = Q_dot_rec_ub*ones(1,length(x0_vec));   
    res_x0 = NaN*ones(1,length(x0_vec));
    for k = 1:length(x0_vec)
        index = index+1;
        dispstat(['x0 evaluation: : ' num2str(index) '/' num2str(length(P_cd_guess_vec)+length(x0_vec))])
        out_x0 = OrganicRankineCycle([P_ev_guess_vec_x0(k) P_cd_guess_vec_x0(k) Q_dot_rec_guess_vec_x0(k)]./[P_ev_ub_vec_x0(k) P_cd_ub_vec_x0(k) Q_dot_rec_ub_vec_x0(k)], [P_ev_lb_vec_x0(k) P_cd_lb_vec_x0(k) Q_dot_rec_lb_vec_x0(k)], [P_ev_ub_vec_x0(k) P_cd_ub_vec_x0(k) Q_dot_rec_ub_vec_x0(k)], fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC);
        if any(out_x0.flag.value < 0)
            res_x0(k) = NaN;
        else
            res_x0(k) =  out_x0.res;
        end
    end   
    P_cd_guess_vec = [P_cd_guess_vec_x0 P_cd_guess_vec];
    P_cd_lb_vec = [P_cd_lb_vec_x0 P_cd_lb_vec];
    P_cd_ub_vec = [P_cd_ub_vec_x0 P_cd_ub_vec];    
    P_ev_guess_vec = [P_ev_guess_vec_x0 P_ev_guess_vec];
    P_ev_lb_vec = [P_ev_lb_vec_x0 P_ev_lb_vec];
    P_ev_ub_vec = [P_ev_ub_vec_x0 P_ev_ub_vec];  
    Q_dot_rec_guess_vec = [Q_dot_rec_guess_vec_x0 Q_dot_rec_guess_vec];
    Q_dot_rec_lb_vec = [Q_dot_rec_lb_vec_x0 Q_dot_rec_lb_vec];
    Q_dot_rec_ub_vec = [Q_dot_rec_ub_vec_x0 Q_dot_rec_ub_vec];
    res = [res_x0 , res];
end

P_cd_guess_vec = P_cd_guess_vec(not(isnan(res)));
P_cd_lb_vec = P_cd_lb_vec(not(isnan(res)));
P_cd_ub_vec = P_cd_ub_vec(not(isnan(res)));
P_ev_guess_vec = P_ev_guess_vec(not(isnan(res)));
P_ev_lb_vec = P_ev_lb_vec(not(isnan(res)));
P_ev_ub_vec = P_ev_ub_vec(not(isnan(res)));
Q_dot_rec_guess_vec = Q_dot_rec_guess_vec(not(isnan(res)));
Q_dot_rec_lb_vec = Q_dot_rec_lb_vec(not(isnan(res)));
Q_dot_rec_ub_vec = Q_dot_rec_ub_vec(not(isnan(res)));
res = res(not(isnan(res)));
[res_ordered, j_order] = sort(res);
%res_ordered = res_ordered(res_ordered<1.5);
%j_order = j_order(res_ordered<1.5);
Nbr_comb_x0 = length(j_order);
Nbr_comb_x0_max = 3;

disp('x0 residuals and index:')
fprintf('\n');
disp(num2str([j_order(1:min(Nbr_comb_x0,Nbr_comb_x0_max)); res_ordered(1:min(Nbr_comb_x0,Nbr_comb_x0_max))]))
fprintf('\n');

if not(isempty(res_ordered))
    
    fprintf('\n');
    disp('Start iteration:')
    fprintf('\n');
    fprintf('%-10s %-5s %-40s %-15s %-40s %-60s %-10s %-100s\n', '#', 'i0', 'x_in', 'res_in', 'x_out', 'res_out', 'flag_ORC', 'flag');
    fprintf('\n');
    
    k= 1;
    out_ORC_best.res = 1e10;
    stop = 0;
    %options_fsolve = optimoptions('fsolve', 'Algorithm', 'Levenberg-Marquardt', 'Display','iter','TolX', 1e-10, 'TolFun', 1e-10, 'MaxIter', 1e9, 'MaxFunEvals', 1e9,'OutputFcn',@ outputfunFS);
    options_pattsearch = psoptimset('Display','iter','TolX', 1e-10, 'TolFun', 1e-10, 'TolMesh', 1e-10, 'MaxIter', 1e4, 'MaxFunEvals', 1e8, 'OutputFcns',@outputfunPS);
    %options_fminsearch = optimset('Display','iter','TolX', 1e-10, 'TolFun', 1e-10, 'MaxIter', 1e9, 'MaxFunEvals', 1e9,'OutputFcn',@ outputfunFS);

    while not(stop) && k <= min(Nbr_comb_x0,Nbr_comb_x0_max);
      
        x0 = [P_ev_guess_vec(j_order(k))    P_cd_guess_vec(j_order(k))  Q_dot_rec_guess_vec(j_order(k))];
        ub = [P_ev_ub_vec(j_order(k))       P_cd_ub_vec(j_order(k))     Q_dot_rec_ub_vec(j_order(k))];
        lb = [P_ev_lb_vec(j_order(k))       P_cd_lb_vec(j_order(k))     Q_dot_rec_lb_vec(j_order(k))];
        fprintf('%-10s %-5d %-40s %-15s ', [num2str(k) '/' num2str(min(Nbr_comb_x0,Nbr_comb_x0_max))] , j_order(k), ['[' num2str([x0(1)/1e5  x0(2)/1e5 x0(3)/1e3]) ']'] , num2str(res(j_order(k)), '%.4g'));
        f = @(x) OrganicRankineCycle_res( x, lb, ub, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC);
        %[x, ~, ~, output] = fsolve(f,x0./ub, options_fsolve);
        %[x, ~, ~, output] = fminsearch(f,x0./ub, options_fminsearch);
        [x, ~, ~, output] = patternsearch(f,x0./ub,[],[],[],[],lb./ub,ub./ub,[],options_pattsearch);
        [out_ORC, TS_ORC] = OrganicRankineCycle(x, lb, ub, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC);      
        out_ORC.funcCount = NaN; %output.funcCount;
        out_ORC.iterations = NaN; %output.iterations;

        if (abs(out_ORC.res_ORC_Hsu) < 1e-5) && (abs(out_ORC.res_ORC_Mdot) < 1e-5) && (abs(out_ORC.res_ORC_Qdot_rec) < 1e-5) && not(any(out_ORC.flag.value<0))
            out_ORC.flag_ORC = 1;
            stop = 1;
        else
            out_ORC.flag_ORC = - 1;
            stop = 0;
        end
        fprintf('%-40s %-60s %-10s %-100s \n', ['[' num2str([x(1)*ub(1)/1e5  x(2)*ub(2)/1e5 x(3)*ub(3)/1e3]) ']'], [ num2str(out_ORC.res, '%.4g') '  [ ' num2str([out_ORC.res_ORC_Hsu  out_ORC.res_ORC_Mdot   out_ORC.res_ORC_Qdot_rec]) ' ] '], num2str(out_ORC.flag_ORC), num2str(out_ORC.flag.value));
        
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
fprintf('\n')
dispstat('','keepprev')
out_ORC = orderfields(out_ORC);
if ORC.displayResults ==1 && out_ORC.flag_ORC == - 1
    dispstat('','keepprev')
    dispstat('Error: The model did not converge correctly','keepthis');
end


%% TS DIAGRAM and DISPLAY
if ORC.displayTS == 1
    figure
    hold all
    [~,~, ~] = Ts_diagram(TS_ORC);
    hold off
    grid on
    xlabel('Entropy [J/kg.K]','fontsize',14,'fontweight','bold')
    ylabel('Temperature [°C]','fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')  
end

if ORC.displayResults ==1
    in.fluid_wf = fluid_wf;
    in.fluid_htf = fluid_htf;
    in.in_htf_su = in_htf_su;
    in.P_htf_su = P_htf_su;
    in.m_dot_htf = m_dot_htf;
    in.fluid_ctf = fluid_ctf;
    in.in_ctf_su = in_ctf_su;
    in.P_ctf_su = P_ctf_su;
    in.m_dot_ctf = m_dot_ctf;
    in.T_amb = T_amb;
    in.N_exp = N_exp;
    in.N_pp = N_pp;
    in.DT_sc = ORC.DT_sc;
    in.PP_modelType = ORC.PP.modelType;
    in.EV_modelType = ORC.EV.modelType;
    in.EXP_modelType = ORC.EXP.modelType;
    in.CD_modelType = ORC.CD.modelType;
    in.REC_modelType = ORC.REC.modelType;
    in.SUB_modelType = ORC.SUB.modelType;
    in.PRE_modelType = ORC.PRE.modelType;
    if nargin ==0
        fprintf ( 1, '\n' );
        disp('-------------------------------------------------------')
        disp('--------------------   Demo Code   --------------------')
        disp('-------------------------------------------------------')
        fprintf ( 1, '\n' );
    end
    disp('Working conditions:')
    fprintf ( 1, '\n' );
    disp(in)
    disp('Results:')
    disp(out_ORC)
end

end

function res = OrganicRankineCycle_res( x, lb, ub, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC)
[out, ~] = OrganicRankineCycle(x, lb, ub, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC);
res = out.res;
%res = [out.res_ORC_Mdot    out.res_ORC_Hsu     out.res_ORC_Qdot_rec];

end

function [out, TS] = OrganicRankineCycle(x, lb, ub, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC)
x = max(x, lb./ub);
x = min(x,ones(1,length(x)));
x = x.*ub;
x(1) = max(x(1), 1.001*x(2));
out.rp_pp = x(1)/x(2);
i_flag = 1;
out.DT_sc = ORC.DT_sc;

% PUMP
out.P_pp_su = x(2);
out.T_pp_su = CoolProp.PropsSI('T', 'P', out.P_pp_su, 'Q', 0, fluid_wf) - ORC.DT_sc;
out.h_pp_su = CoolProp.PropsSI('H', 'P', out.P_pp_su, 'T', out.T_pp_su, fluid_wf);
[out_PP, TS_PP] = PumpModel(out.P_pp_su, out.h_pp_su, x(1), fluid_wf, N_pp, ORC.PP);
out.m_dot_wf = out_PP.m_dot;
out.W_dot_pp = out_PP.W_dot;
out.eps_is_pp = out_PP.epsilon_is;
out.eps_vol_pp = out_PP.epsilon_vol;
out.time_pp = out_PP.time;
out.flag_pp = out_PP.flag;
T_prev = out_PP.T_ex;
h_prev = out_PP.h_ex;
P_prev = x(1);
out.flag.value(1,i_flag) = out_PP.flag;
out.flag.name{1,i_flag} = 'flag_pp';
TS.PP = TS_PP;

% RECUPERATOR (cold side)
if isfield(ORC, 'REC')
    out.P_recc_su = P_prev;
    out.T_recc_su = out_PP.T_ex;
    out.h_recc_su = out_PP.h_ex;
    out.Q_dot_rec_bis = x(3);
    h_prev = min(out.h_recc_su + x(3)/out.m_dot_wf, CoolProp.PropsSI('H', 'P', out.P_recc_su, 'T', T_htf_su, fluid_wf));
    P_prev = out.P_recc_su;
    T_prev = CoolProp.PropsSI('T', 'P', P_prev, 'H', h_prev, fluid_wf);
end

% PREHEATER &/ EVAPORATOR
if isfield(ORC, 'PRE')
    out.P_pre_su = P_prev;
    out.T_pre_su = T_prev;
    out.h_pre_su = h_prev;
    pre_ev = HexSeries(fluid_htf, P_htf_su, in_htf_su, m_dot_htf, fluid_wf, out.P_pre_su, out.h_pre_su, out.m_dot_wf, ORC.PRE, ORC.EV);
    TS_PRE = pre_ev.ts1;
    out.flag_pre_ev = pre_ev.flag;
    out.Q_dot_pre = pre_ev.hex1.Q_dot_tot;
    out.T_htf_pre_ex = pre_ev.hex1.T_h_ex;
    out.h_htf_pre_ex = pre_ev.hex1.h_h_ex;
    out.flag_pre = pre_ev.hex1.flag;
    out.pinch_pre = pre_ev.hex1.pinch;
    out.time_pre = pre_ev.hex1.time;
    out.h_ev_su = pre_ev.hex1.h_c_ex;
    out.T_ev_su = pre_ev.hex1.T_c_ex;
    out.P_ev_su = out.P_pre_su;
    TS_EV = pre_ev.ts2;
    out.Q_dot_ev = pre_ev.hex2.Q_dot_tot;
    out.T_htf_ev_ex = pre_ev.hex2.T_h_ex;
    out.h_htf_ev_ex = pre_ev.hex2.h_h_ex;
    out.flag_ev = pre_ev.hex2.flag;
    out.pinch_ev = pre_ev.hex2.pinch;
    out.time_ev = pre_ev.hex2.time;
    Q_dot_in = out.Q_dot_ev + out.Q_dot_pre;
    T_prev = pre_ev.hex2.T_c_ex;
    h_prev = pre_ev.hex2.h_c_ex;
    P_prev = out.P_ev_su;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = pre_ev.flag; 
    out.flag.name{1,i_flag} = 'flag_pre_ev';
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = pre_ev.hex1.flag;
    out.flag.name{1,i_flag} = 'flag_pre';
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = pre_ev.hex2.flag;
    out.flag.name{1,i_flag} = 'flag_ev';
    TS.PRE = TS_PRE;
    TS.EV = TS_EV;
else
    out.P_ev_su = P_prev;
    out.T_ev_su = T_prev;
    out.h_ev_su = h_prev;
    [out_EV, TS_EV] = HexModel(fluid_htf, P_htf_su, in_htf_su, m_dot_htf, fluid_wf, out.P_ev_su, out.h_ev_su, out.m_dot_wf , ORC.EV);
    out.Q_dot_ev = out_EV.Q_dot_tot;
    Q_dot_in = out.Q_dot_ev;
    out.T_htf_ev_ex = out_EV.T_h_ex;
    out.h_htf_ev_ex = out_EV.h_h_ex;
    out.flag_ev = out_EV.flag;
    out.time_ev = out_EV.time;
    out.pinch_ev = out_EV.pinch;
    P_prev = out.P_ev_su;
    T_prev = out_EV.T_c_ex;
    h_prev = out_EV.h_c_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_EV.flag;
    out.flag.name{1,i_flag} = 'flag_ev';
    TS.EV = TS_EV;
end

% DPHP
if isfield(ORC, 'DPHP')
    out.P_dphp_su = P_prev;
    out.T_dphp_su = T_prev;
    out.h_dphp_su = h_prev;
    ORC.DPHP.type_in = 'su';
    [out_DPHP, TS_DPHP] = DPModel_2(fluid_wf, out.P_dphp_su, out.h_dphp_su, out.m_dot_wf, T_amb, ORC.DPHP);
    out.dphp = out_DPHP.dp;
    out.Q_dot_hp = out_DPHP.Q_dot;
    T_prev = out_DPHP.T_ex;
    h_prev = out_DPHP.h_ex;
    P_prev = out_DPHP.P_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_DPHP.flag;
    out.flag.name{1,i_flag} = 'flag_dphp';
    TS.DPHP = TS_DPHP;
end

% EXPANDER
out.P_exp_su = P_prev;
out.T_exp_su = T_prev;
out.h_exp_su = h_prev;
if isfield(ORC, 'DPLP')
    ORC.DPLP.type_in = 'ex';
    [out_DPLP_bis, ~] = DPModel_2(fluid_wf, x(2), out.h_pp_su, out.m_dot_wf, T_amb, ORC.DPLP);
    P_exp_ex = x(2)+out_DPLP_bis.dp;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_DPLP_bis.flag;
    out.flag.name{1,i_flag} = 'flag_dplp_bis';
else
    P_exp_ex = x(2);
end
out.rp_exp = out.P_exp_su/P_exp_ex;
[out_EXP, TS_EXP] = ExpanderModel(fluid_wf, out.P_exp_su, out.h_exp_su, N_exp, P_exp_ex, T_amb, ORC.EXP);
out.m_dot_wf_bis = out_EXP.M_dot;
out.W_dot_exp = out_EXP.W_dot;
out.Q_dot_exp = out_EXP.Q_dot_amb;
out.eps_vol_exp = out_EXP.FF;
out.eps_is_exp = out_EXP.epsilon_s;
out.flag_exp = out_EXP.flag;
out.time_exp = out_EXP.time;
T_prev = out_EXP.T_ex;
h_prev = out_EXP.h_ex;
P_prev = P_exp_ex;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_EXP.flag;
out.flag.name{1,i_flag} = 'flag_exp';
TS.EXP = TS_EXP; 

% RECUPERATOR (hot side)
if isfield(ORC, 'REC')
    out.P_rech_su = P_prev;
    out.T_rech_su = T_prev;
    out.h_rech_su = h_prev;
    [out_REC, TS_REC] = HexModel(fluid_wf, out.P_rech_su, out.h_rech_su, out.m_dot_wf, fluid_wf, out.P_recc_su, out.h_recc_su, out.m_dot_wf, ORC.REC);
    out.Q_dot_rec = out_REC.Q_dot_tot;
    out.flag_rec = out_REC.flag;
    out.time_rec = out_EXP.time;
    out.pinch_rec = out_REC.pinch;
    h_prev = out_REC.h_h_ex;
    P_prev = out.P_rech_su;    
    T_prev = out_REC.T_h_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_REC.flag;
    out.flag.name{1,i_flag} = 'flag_rec';
    TS.REC = TS_REC;
end

% CONDENSER &/ SUBCOOLER
if isfield(ORC, 'SUB')
    out.P_cd_su = P_prev;
    out.T_cd_su = T_prev;
    out.h_cd_su = h_prev;
    sub_cd = HexSeries(fluid_wf, out.P_cd_su, out.h_cd_su, out.m_dot_wf, fluid_ctf, P_ctf_su, in_ctf_su, m_dot_ctf, ORC.SUB, ORC.CD);
    out.flag_sub_cd = sub_cd.flag;
    out.Q_dot_cd = sub_cd.hex2.Q_dot_tot;
    out.T_ctf_cd_ex = sub_cd.hex2.T_c_ex;
    out.h_ctf_cd_ex = sub_cd.hex2.h_c_ex;
    out.flag_cd = sub_cd.hex2.flag;
    out.pinch_cd = sub_cd.hex2.pinch;
    out.time_cd = sub_cd.hex2.time;
    out.Q_dot_sub = sub_cd.hex1.Q_dot_tot;
    out.T_ctf_cd_su = sub_cd.hex1.T_c_ex;
    out.h_ctf_cd_su = sub_cd.hex1.h_c_ex;
    out.flag_sub = sub_cd.hex1.flag;
    out.pinch_sub = sub_cd.hex1.pinch;
    out.time_sub = sub_cd.hex1.time;
    out.h_sub_su = sub_cd.hex2.h_h_ex;
    out.T_sub_su = sub_cd.hex2.T_h_ex;
    out.P_sub_su = out.P_cd_su;    
    T_prev = sub_cd.hex1.T_h_ex;
    h_prev = sub_cd.hex1.h_h_ex;
    P_prev = out.P_sub_su;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = sub_cd.flag;
    out.flag.name{1,i_flag} = 'flag_sub_cd';
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = sub_cd.hex1.flag;
    out.flag.name{1,i_flag} = 'flag_sub';
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = sub_cd.hex2.flag;
    out.flag.name{1,i_flag} = 'flag_cd';
    TS.SUB = sub_cd.ts1;
    TS.CD = sub_cd.ts2;
else
    out.P_cd_su = P_prev;
    out.T_cd_su = T_prev;
    out.h_cd_su = h_prev;
    [out_CD, TS_CD] = HexModel(fluid_wf, out.P_cd_su, out.h_cd_su, out.m_dot_wf, fluid_ctf, P_ctf_su, in_ctf_su, m_dot_ctf , ORC.CD);
    out.Q_dot_cd = out_CD.Q_dot_tot;
    out.T_ctf_cd_ex = out_CD.T_c_ex;
    out.h_ctf_cd_ex = out_CD.h_c_ex;
    out.flag_cd = out_CD.flag;
    out.time_cd = out_CD.time;
    out.pinch_cd = out_CD.pinch;
    P_prev = out.P_cd_su;
    T_prev = out_CD.T_h_ex;
    h_prev = out_CD.h_h_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_CD.flag;
    out.flag.name{1,i_flag} = 'flag_cd';
    TS.CD = TS_CD;
end

% DPLP
if isfield(ORC, 'DPLP')
    out.T_dplp_su = T_prev;
    out.h_dplp_su = h_prev;
    out.P_dplp_su = P_prev;
    ORC.DPLP.type_in = 'su';
    [out_DPLP, TS_DPLP] = DPModel_2(fluid_wf, out.P_dplp_su, out.h_dplp_su, out.m_dot_wf, T_amb, ORC.DPLP);
    out.dplp = out_DPLP.dp;
    out.Q_dot_lp = out_DPLP.Q_dot;
    T_prev = out_DPLP.T_ex;
    h_prev = out_DPLP.h_ex;
    P_prev = out_DPLP.P_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_DPLP.flag;
    out.flag.name{1,i_flag} = 'flag_dplp';
    TS.DPLP = TS_DPLP;
end
out.DT_bis =  CoolProp.PropsSI('T', 'P', P_prev, 'Q', 0, fluid_wf)-T_prev;

% ORC PERFORAMANCE
out.W_dot_net = out.W_dot_exp - out.W_dot_pp;
out.eff_ORC_gross = out.W_dot_exp/Q_dot_in;
out.eff_ORC_net = out.W_dot_net/Q_dot_in;

% RESIDUALS and RESULTS
out.res_ORC_Hsu = (1 - out.h_pp_su/h_prev);
%out.res_ORC_Hsu = (1 - out.DT_bis/out.DT_sc);

if out.m_dot_wf == 0 && out.m_dot_wf_bis == 0
    out.res_ORC_Mdot = 0;
elseif out.m_dot_wf == 0 && out.m_dot_wf_bis ~= 0
    out.res_ORC_Mdot = 1;
elseif out.m_dot_wf ~= 0 && out.m_dot_wf_bis == 0
    out.res_ORC_Mdot = 1;
else
    out.res_ORC_Mdot = 1-out.m_dot_wf/out.m_dot_wf_bis;
end

if x(3) == 0 && out.Q_dot_rec == 0
    out.res_ORC_Qdot_rec = 0;
elseif x(3) == 0 && out.Q_dot_rec ~= 0
    out.res_ORC_Qdot_rec = 1;
elseif x(3) ~= 0 && out.Q_dot_rec == 0
    out.res_ORC_Qdot_rec = 1;
else
    out.res_ORC_Qdot_rec = 1 - x(3)/out.Q_dot_rec;
end

% Residuals
if any(out.flag.value < 0)
    out.res_ORC_Mdot = 1e10*out.res_ORC_Mdot;
    out.res_ORC_Hsu = 1e10*out.res_ORC_Hsu;
    out.res_ORC_Qdot_rec = 1e10*out.res_ORC_Qdot_rec;
    out.res  = norm([out.res_ORC_Mdot    out.res_ORC_Hsu     out.res_ORC_Qdot_rec]);
else
    out.res = norm([out.res_ORC_Mdot    out.res_ORC_Hsu     out.res_ORC_Qdot_rec]);
end
out.x = x;


% Construction of the TS variable
TS.cycle.s = TS.PP.s;
TS.cycle.T = TS.PP.T;
if isfield(ORC, 'REC')
    TS.cycle.s = [TS.cycle.s  TS.REC.s_c];
    TS.cycle.T = [TS.cycle.T  TS.REC.T_c];
end
if isfield(ORC, 'PRE')
    TS.cycle.s = [TS.cycle.s  TS.PRE.s_c];
    TS.cycle.T = [TS.cycle.T  TS.PRE.T_c];
end
TS.cycle.s = [TS.cycle.s  TS.EV.s_c];
TS.cycle.T = [TS.cycle.T  TS.EV.T_c];
if isfield(ORC, 'DPHP')
    TS.cycle.s = [TS.cycle.s TS.DPHP.s];
    TS.cycle.T = [TS.cycle.T TS.DPHP.T];
end
TS.cycle.s = [TS.cycle.s  TS.EXP.s];
TS.cycle.T = [TS.cycle.T  TS.EXP.T];
if isfield(ORC, 'REC')
    TS.cycle.s = [TS.cycle.s  fliplr(TS.REC.s_h)];
    TS.cycle.T = [TS.cycle.T  fliplr(TS.REC.T_h)];
end
TS.cycle.s = [TS.cycle.s fliplr(TS.CD.s_h)];
TS.cycle.T = [TS.cycle.T fliplr(TS.CD.T_h)];
if isfield(ORC, 'SUB')
    TS.cycle.s = [TS.cycle.s  fliplr(TS.SUB.s_h)];
    TS.cycle.T = [TS.cycle.T  fliplr(TS.SUB.T_h)];
end
if isfield(ORC, 'DPLP')
    TS.cycle.s = [TS.cycle.s TS.DPLP.s];
    TS.cycle.T = [TS.cycle.T TS.DPLP.T];
end
    
end

function out = OrganicRankineCycle_init(P_low, P_high, z, fluid_wf, fluid_htf, in_htf_su, T_htf_su, P_htf_su, m_dot_htf, fluid_ctf, in_ctf_su, T_ctf_su, P_ctf_su, m_dot_ctf, T_amb, N_exp, N_pp, ORC)
out.rp = P_high/P_low;
i_flag = 1;


% PUMP
out.P_pp_su = P_low;
out.T_pp_su = CoolProp.PropsSI('T', 'P', out.P_pp_su, 'Q', 0, fluid_wf) - ORC.DT_sc;
out.h_pp_su = CoolProp.PropsSI('H', 'P', out.P_pp_su, 'T', out.T_pp_su, fluid_wf);
[out_PP, ~] = PumpModel(out.P_pp_su, out.h_pp_su, P_high, fluid_wf, N_pp, ORC.PP);
out.m_dot_wf = out_PP.m_dot;
out.flag_pp = out_PP.flag;
T_prev = out_PP.T_ex;
h_prev = out_PP.h_ex;
P_prev = P_high;
out.flag.value(1,i_flag) = out_PP.flag;
out.flag.name{1,i_flag} = 'flag_pp';

% DPLP
if isfield(ORC, 'DPLP')
    ORC.DPLP.type_in = 'ex';
    [out_DPLP_bis, ~] = DPModel_2(fluid_wf, P_low, out.h_pp_su, out.m_dot_wf, T_amb, ORC.DPLP);
    P_exp_ex = P_low+out_DPLP_bis.dp;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_DPLP_bis.flag;
    out.flag.name{1,i_flag} = 'flag_dplp_bis';
else
    P_exp_ex = x(2);
end


% RECUPERATOR (cold side)
if isfield(ORC, 'REC')
    out.P_recc_su = P_prev;
    out.T_recc_su = out_PP.T_ex;
    out.h_recc_su = out_PP.h_ex;
    out.Q_dot_rec_max = HEX_Qdotmax(fluid_wf, out.m_dot_wf, P_exp_ex, CoolProp.PropsSI('H', 'P', P_exp_ex, 'T',T_htf_su, fluid_wf), fluid_wf, out.m_dot_wf, out.P_recc_su, out.h_recc_su, ORC.REC); % CoolProp.PropsSI('T', 'P', P_exp_ex, 'Q', 0, fluid_wf) + 100
    out.Q_dot_rec_guess = z(1)*out.Q_dot_rec_max;
    h_prev = min(out.h_recc_su + out.Q_dot_rec_guess/out.m_dot_wf, CoolProp.PropsSI('H', 'P', out.P_recc_su, 'T', T_htf_su, fluid_wf));
    P_prev = out.P_recc_su;
    T_prev = CoolProp.PropsSI('T', 'P', P_prev, 'H', h_prev, fluid_wf);
end

% PREHEATER &/ EVAPORATOR
if isfield(ORC, 'PRE')

    out.P_pre_su = P_prev;
    out.T_pre_su = T_prev;
    out.h_pre_su = h_prev;
    pre_ev = HexSeries(fluid_htf, P_htf_su, in_htf_su, m_dot_htf, fluid_wf, out.P_pre_su, out.h_pre_su, out.m_dot_wf, ORC.PRE, ORC.EV);
    out.flag_pre_ev = pre_ev.flag;
    out.Q_dot_pre = pre_ev.hex1.Q_dot_tot;
    out.T_htf_pre_ex = pre_ev.hex1.T_h_ex;
    out.h_htf_pre_ex = pre_ev.hex1.h_h_ex;
    out.flag_pre = pre_ev.hex1.flag;
    out.pinch_pre = pre_ev.hex1.pinch;
    out.time_pre = pre_ev.hex1.time;
    out.h_ev_su = pre_ev.hex1.h_c_ex;
    out.T_ev_su = pre_ev.hex1.T_c_ex;
    out.P_ev_su = out.P_pre_su;
    out.Q_dot_ev = pre_ev.hex2.Q_dot_tot;
    out.T_htf_ev_ex = pre_ev.hex2.T_h_ex;
    out.h_htf_ev_ex = pre_ev.hex2.h_h_ex;
    out.flag_ev = pre_ev.hex2.flag;
    out.pinch_ev = pre_ev.hex2.pinch;
    out.time_ev = pre_ev.hex2.time;
    T_prev = pre_ev.hex2.T_c_ex;
    h_prev = pre_ev.hex2.h_c_ex;
    P_prev = out.P_ev_su;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = pre_ev.flag; 
    out.flag.name{1,i_flag} = 'flag_pre_ev';
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = pre_ev.hex1.flag;
    out.flag.name{1,i_flag} = 'flag_pre';
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = pre_ev.hex2.flag;
    out.flag.name{1,i_flag} = 'flag_ev';
else
    out.P_ev_su = P_prev;
    out.T_ev_su = T_prev;
    out.h_ev_su = h_prev;
    [out_EV, ~] = HexModel(fluid_htf, P_htf_su, in_htf_su, m_dot_htf, fluid_wf, out.P_ev_su, out.h_ev_su, out.m_dot_wf , ORC.EV);
    out.Q_dot_ev = out_EV.Q_dot_tot;
    out.T_htf_ev_ex = out_EV.T_h_ex;
    out.h_htf_ev_ex = out_EV.h_h_ex;
    out.flag_ev = out_EV.flag;
    out.time_ev = out_EV.time;
    out.pinch_ev = out_EV.pinch;
    P_prev = out.P_ev_su;
    T_prev = out_EV.T_c_ex;
    h_prev = out_EV.h_c_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_EV.flag;
    out.flag.name{1,i_flag} = 'flag_ev';
end


% DPHP
if isfield(ORC, 'DPHP')
    out.P_dphp_su = P_prev;
    out.T_dphp_su = T_prev;
    out.h_dphp_su = h_prev;
    ORC.DPHP.type_in = 'su';
    [out_DPHP, ~] = DPModel_2(fluid_wf, out.P_dphp_su, out.h_dphp_su, out.m_dot_wf, T_amb, ORC.DPHP);
    out.dphp = out_DPHP.dp;
    T_prev = out_DPHP.T_ex;
    h_prev = out_DPHP.h_ex;
    P_prev = out_DPHP.P_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_DPHP.flag;
    out.flag.name{1,i_flag} = 'flag_dphp';
end

% EXPANDER
out.P_exp_su = P_prev;
out.T_exp_su = T_prev;
out.h_exp_su = h_prev;
[out_EXP, ~] = ExpanderModel(fluid_wf, out.P_exp_su, out.h_exp_su, N_exp, P_exp_ex, T_amb, ORC.EXP);
out.m_dot_wf_bis = out_EXP.M_dot;
out.W_dot_exp = out_EXP.W_dot;
out.Q_dot_exp = out_EXP.Q_dot_amb;
out.eps_vol_exp = out_EXP.FF;
out.eps_is_exp = out_EXP.epsilon_s;
out.flag_exp = out_EXP.flag;
out.time_exp = out_EXP.time;
T_prev = out_EXP.T_ex;
h_prev = out_EXP.h_ex;
P_prev = P_exp_ex;
i_flag = i_flag+1;
out.flag.value(1,i_flag) = out_EXP.flag;
out.flag.name{1,i_flag} = 'flag_exp';

% RECUPERATOR (hot side)
if isfield(ORC, 'REC')
    out.P_rech_su = P_prev;
    out.T_rech_su = T_prev;
    out.h_rech_su = h_prev;
    [out_REC, ~] = HexModel(fluid_wf, out.P_rech_su, out.h_rech_su, out.m_dot_wf, fluid_wf, out.P_recc_su, out.h_recc_su, out.m_dot_wf, ORC.REC);
    out.Q_dot_rec = out_REC.Q_dot_tot;
    out.flag_rec = out_REC.flag;
    out.time_rec = out_EXP.time;
    out.pinch_rec = out_REC.pinch;
    h_prev = out_REC.h_h_ex;
    P_prev = out.P_rech_su;    
    T_prev = out_REC.T_h_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_REC.flag;
    out.flag.name{1,i_flag} = 'flag_rec';
end

% CONDENSER &/ SUBCOOLER
if isfield(ORC, 'SUB')
    out.P_cd_su = P_prev;
    out.T_cd_su = T_prev;
    out.h_cd_su = h_prev;
    sub_cd = HexSeries(fluid_wf, out.P_cd_su, out.h_cd_su, out.m_dot_wf, fluid_ctf, P_ctf_su, in_ctf_su, m_dot_ctf, ORC.SUB, ORC.CD);
    out.flag_sub_cd = sub_cd.flag;
    out.Q_dot_cd = sub_cd.hex2.Q_dot_tot;
    out.T_ctf_cd_ex = sub_cd.hex2.T_c_ex;
    out.h_ctf_cd_ex = sub_cd.hex2.h_c_ex;
    out.flag_cd = sub_cd.hex2.flag;
    out.pinch_cd = sub_cd.hex2.pinch;
    out.time_cd = sub_cd.hex2.time;
    out.Q_dot_sub = sub_cd.hex1.Q_dot_tot;
    out.T_ctf_cd_su = sub_cd.hex1.T_c_ex;
    out.h_ctf_cd_su = sub_cd.hex1.h_c_ex;
    out.flag_sub = sub_cd.hex1.flag;
    out.pinch_sub = sub_cd.hex1.pinch;
    out.time_sub = sub_cd.hex1.time;
    out.h_sub_su = sub_cd.hex2.h_h_ex;
    out.T_sub_su = sub_cd.hex2.T_h_ex;
    out.P_sub_su = out.P_cd_su;    
    T_prev = sub_cd.hex1.T_h_ex;
    h_prev = sub_cd.hex1.h_h_ex;
    P_prev = out.P_sub_su;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = sub_cd.flag;
    out.flag.name{1,i_flag} = 'flag_sub_cd';
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = sub_cd.hex1.flag;
    out.flag.name{1,i_flag} = 'flag_sub';
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = sub_cd.hex2.flag;
    out.flag.name{1,i_flag} = 'flag_cd';
else
    out.P_cd_su = P_prev;
    out.T_cd_su = T_prev;
    out.h_cd_su = h_prev;
    [out_CD, ~] = HexModel(fluid_wf, out.P_cd_su, out.h_cd_su, out.m_dot_wf, fluid_ctf, P_ctf_su, in_ctf_su, m_dot_ctf , ORC.CD);
    out.Q_dot_cd = out_CD.Q_dot_tot;
    out.T_ctf_cd_ex = out_CD.T_c_ex;
    out.h_ctf_cd_ex = out_CD.h_c_ex;
    out.flag_cd = out_CD.flag;
    out.time_cd = out_CD.time;
    out.pinch_cd = out_CD.pinch;
    P_prev = out.P_cd_su;
    T_prev = out_CD.T_h_ex;
    h_prev = out_CD.h_h_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_CD.flag;
    out.flag.name{1,i_flag} = 'flag_cd';
end

% DPLP
if isfield(ORC, 'DPLP')
    out.T_dplp_su = T_prev;
    out.h_dplp_su = h_prev;
    out.P_dplp_su = P_prev;
    ORC.DPLP.type_in = 'su';
    [out_DPLP, ~] = DPModel_2(fluid_wf, out.P_dplp_su, out.h_dplp_su, out.m_dot_wf, T_amb, ORC.DPLP);
    out.dplp = out_DPLP.dp;
    out.Q_dot_dplp = out_DPLP.Q_dot;
    T_prev = out_DPLP.T_ex;
    h_prev = out_DPLP.h_ex;
    P_prev = out_DPLP.P_ex;
    i_flag = i_flag+1;
    out.flag.value(1,i_flag) = out_DPLP.flag;
    out.flag.name{1,i_flag} = 'flag_dplp';
end
out.DT_bis =  CoolProp.PropsSI('T', 'P', P_prev, 'Q', 0, fluid_wf)-T_prev;

% RESIDUALS and RESULTS
out.res_ORC_Hsu = (1 - out.h_pp_su/h_prev);

if out.m_dot_wf == 0 && out.m_dot_wf_bis == 0
    out.res_ORC_Mdot = 0;
elseif out.m_dot_wf == 0 && out.m_dot_wf_bis ~= 0
    out.res_ORC_Mdot = 1;
elseif out.m_dot_wf ~= 0 && out.m_dot_wf_bis == 0
    out.res_ORC_Mdot = 1;
else
    out.res_ORC_Mdot = 1-out.m_dot_wf/out.m_dot_wf_bis;
end

if out.Q_dot_rec_guess == 0 && out.Q_dot_rec == 0
    out.res_ORC_Qdot_rec = 0;
elseif out.Q_dot_rec_guess == 0 && out.Q_dot_rec ~= 0
    out.res_ORC_Qdot_rec = 1;
elseif out.Q_dot_rec_guess ~= 0 && out.Q_dot_rec == 0
    out.res_ORC_Qdot_rec = 1;
else
    out.res_ORC_Qdot_rec = 1 - out.Q_dot_rec_guess/out.Q_dot_rec;
end

if any(out.flag.value < 0)
    out.res_ORC_Mdot = 1e2*out.res_ORC_Mdot;
    out.res_ORC_Hsu = 1e2*out.res_ORC_Hsu;
    out.res_ORC_Qdot_rec = 1e2*out.res_ORC_Qdot_rec; 
    out.res  = 1e2*norm([out.res_ORC_Mdot    out.res_ORC_Hsu     out.res_ORC_Qdot_rec]);
else
    out.res = norm([out.res_ORC_Mdot    out.res_ORC_Hsu     out.res_ORC_Qdot_rec]);
end

end

function [stop,options,optchanged] = outputfunPS(optimvalues,options,flag)
stop = optimvalues.fval < 1e-5;
optchanged = 0;
end

function stop = outputfunFS(x, optimValues, state)
%disp(norm(optimValues.fval))
stop = norm(optimValues.fval) < 1e-5;
end