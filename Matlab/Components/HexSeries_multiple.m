function out = HexSeries_multiple(fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param) %in_hex1, in_hex2)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 28/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% HexSeries is a single matlab function used to evaluate the heat transfer 
% occuring in two heat exchangers in series (see the Documentation/HexSeries_MatlabDoc)
%
% The model inputs are:
%       - fluid_h: nature of the hot fluid                        	[-]
%       - P_h_su: inlet pressure of the hot fluid                   [Pa]
%       - in_h_su: inlet temperature or enthalpy of the hot fluid   [K or J/kg]
%       - m_dot_h: mass flow rate of the hot fluid                  [kg/s]
%       - fluid_c: nature of the cold fluid                        	[-]
%       - P_c_su: inlet pressure of the cold fluid                  [Pa]
%       - in_c_su: inlet temperature or enthalpy of the cold fluid  [K or J/kg]
%       - m_dot_c: mass flow rate of the cold fluid                 [kg/s]
%       - in_hex1: structure variable describing the first HEX (cfr HexModel.mat)
%       - in_hex2: structure variable describing the second HEX (cfr HexModel.mat)
%
% The model outputs is:
%       - out: a structure variable which includes all the modelling
%       results of the two heat exchanger (see HexMode.mat)
%
% See the documentation for further details or contact rdickes@ulg.ac.be

%% DEMONSTRATION CASE
if nargin == 0
    fluid_h = 'water';                                               % Nature of the hot fluid           [-]
    m_dot_h = 0.09;                                                          % Mass flow rat of the hot fluid    [kg/s]
    P_h_su =  25e+05;                                                        % Supply pressure of the hot fluid  [Pa]
    in_h_su = CoolProp.PropsSI('H','P',P_h_su,'T',160+273.15, fluid_h);                                                    % Supply h or T of the hot fluid  	[J/kg pr K]
    fluid_c = 'R245fa';                                                     % Nature of the cold fluid        	[-]
    m_dot_c = 0.0252;                                                      % Mass flow rat of the cold fluid  	[kg/s]
    P_c_su = 65.188e5;                                                           % Supply pressure of the cold fluid	[Pa]
    in_c_su = CoolProp.PropsSI('H','P',P_c_su,'T',20+273.15, fluid_c);      % Supply h or T of the cold fluid  	[J/kg pr K]
    L_phex_ev = 0.519-0.06;
    W_phex_ev = 0.191;
    pitch_p_ev = 0.0022;
    th_p_ev = 0.0004;
    b_p_ev = pitch_p_ev-th_p_ev;
    Np_phex_ev = 100;
    Nc_wf_ev = ceil((Np_phex_ev-1)/2);
    Nc_sf_ev = floor((Np_phex_ev-1)/2);
    CS_ev = W_phex_ev*b_p_ev;
    A_eff_ev = 9.8;
    Np_eff_ev = Np_phex_ev-2;
    phi_ev = (A_eff_ev/Np_eff_ev)/(W_phex_ev*L_phex_ev);    
    pitch_co_ev = 0.007213; %computed based on phi_ev and b_phex_ev
    Dh_ev = 2*b_p_ev/phi_ev;
    param.theta = 30*pi/180;

    param.H.A_tot = A_eff_ev;
    param.C.A_tot = A_eff_ev;
    param.H.V_tot = 0.009;
    param.C.V_tot = 0.009;
    param.L_hex = L_phex_ev;
    param.pitch_co = pitch_co_ev;
    param.phi =  phi_ev;
    param.C.Dh = Dh_ev;
    param.H.Dh = Dh_ev;
    param.C.CS = CS_ev;
    param.H.CS = CS_ev;
    param.C.n_canals = Nc_wf_ev;
    param.H.n_canals = Nc_sf_ev;

    param.n_tp_disc = 2;
    param.n_supcrit_disc = 5;
    param.modelType = 'hConvVar';
    param.H.type = 'H';
    param.C.type = 'H';
    param.H.m_dot_n = 1;
    param.H.hConv_vap_n  = 100;
    param.H.n_vap = 0.7;
    param.H.hConv_tp_n  = 2000;
    param.H.n_tp = 0.7;    
    param.H.hConv_liq_n  = 100;
    param.H.n_liq = 0.7;
    param.H.hConv_supcrit_n  = 100;
    param.H.n_supcrit = 0.7; 
    param.C.m_dot_n = 1;
    param.C.hConv_vap_n  = 100;
    param.C.n_vap = 0.7;
    param.C.hConv_tp_n  = 2000;
    param.C.n_tp = 0.7;    
    param.C.hConv_liq_n  = 100;
    param.C.n_liq = 0.7;    
    param.C.hConv_supcrit_n  = 100;
    param.C.n_supcrit = 0.7;    
    
    param.H.correlation.type_1phase = 'Martin';
    param.H.correlation.type_2phase_cd = 'Longo_condensation';
    param.H.correlation.type_2phase_ev = 'Almalfi_boiling';
    param.H.correlation.type_supcrit = 'Griem_SupCrit';
    param.C.correlation.type_1phase = 'Martin';
    param.C.correlation.type_2phase_ev = 'Almalfi_boiling';
    param.C.correlation.type_2phase_cd = 'Longo_condensation';
    
    param.C.correlation.type_supcrit = 'Griem_SupCrit';
    param.H.correlation.void_fraction = 'Hughmark';
    param.C.correlation.void_fraction = 'Hughmark';
    param.hugmark_simplified = 1;
    
    param.H.fact_corr_sp = 0.268088007690337;
    param.H.fact_corr_2p = 1;
    param.C.fact_corr_sp = 0.268088007690337;
    param.C.fact_corr_2p = 1;  
    param.H.fact2_corr_sp = 1;
    param.H.fact2_corr_2p = 1;
    param.C.fact2_corr_sp = 1;
    param.C.fact2_corr_2p = 1;  
    param.displayResults = 0;
    param.displayTS = 0;
    param.generateTS = 1;
    
    param(1) = param;
    param(2) = param(1);
    %param(3) = param(1);

end

%% MODELING OF MULTIPLE HEAT EXCHANGERS IN SERIES
% Evaluation of the hot fluid (HF) supply conditions
if strcmp(param(1).H.type,'H')
    P_h_crit = CoolProp.PropsSI('Pcrit','P',P_h_su,'Q',0,fluid_h);
    h_h_crit = CoolProp.PropsSI('H','P',P_h_su,'T',CoolProp.PropsSI('Tcrit','P',P_h_su,'Q',0,fluid_h)-0.01,fluid_h);
    if P_h_su < P_h_crit
        h_h_l = CoolProp.PropsSI('H','P',P_h_su,'Q',0,fluid_h);
        h_h_v = CoolProp.PropsSI('H','P',P_h_su,'Q',1,fluid_h);
    else
        h_h_l = NaN;
        h_h_v = NaN;
    end
elseif strcmp(param(1).H.type,'T')
    P_h_crit = NaN;
    h_h_crit = NaN;
    h_h_l = NaN;
    h_h_v = NaN;
end

% Evaluation of the cold fluid (CF) supply conditions
if strcmp(param(1).C.type,'H')
    P_c_crit = CoolProp.PropsSI('Pcrit','P',P_c_su,'Q',0,fluid_c);
    h_c_crit = CoolProp.PropsSI('H','P',P_c_su,'T',CoolProp.PropsSI('Tcrit','P',P_c_su,'Q',0,fluid_c)-0.01,fluid_c);
    if P_c_su < P_c_crit
        h_c_l = CoolProp.PropsSI('H','P',P_c_su,'Q', 0, fluid_c);
        h_c_v = CoolProp.PropsSI('H','P',P_c_su,'Q', 1, fluid_c);
    else
        h_c_l = NaN;
        h_c_v = NaN;
    end
elseif strcmp(param(1).C.type,'T')
    P_c_crit = NaN;
    h_c_crit = NaN;
    h_c_l = NaN;
    h_c_v = NaN;
end


N_hex = length(param);
Q_dot_max = HEX_Qdotmax_2(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, param(1), h_h_l, h_h_v, P_h_crit, h_h_crit, h_c_l, h_c_v, P_c_crit, h_c_crit);
lb = 0*ones(1,N_hex);
ub = Q_dot_max*ones(1,N_hex);
f = @(x) hexseries_multiple_res(x, ub, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param, N_hex);

options_fmincon = optimset('Disp','iter','Algorithm','interior-point','UseParallel',false,'TolX',1e-13,'TolFun',1e-13,'TolCon',1e-6,'MaxIter',1e3,'OutputFcn',@outputfunFS);
Q_dot_vec_0 = (ub/N_hex/2);
Q_dot_vec_norm = fmincon(f, Q_dot_vec_0./ub, [], [], [], [], lb./ub, ub./ub, [], options_fmincon);

out = hexseries_multiple(Q_dot_vec_norm, ub, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param, N_hex);
out.Q_dot_vec = Q_dot_vec_norm.*ub;
for i = 1:N_hex
    flag_int(i) = out.int_hex(i).flag;
end


if abs(out.res) < 1e-4 && all(flag_int>0) 
    out.flag = 1;
else
    out.flag = -1;
end

if param(1).displayTS
    figure
    for i = 1:N_hex
        subplot(1,N_hex,i)
        hold on
        plot(out.TS(i).x, out.TS(i).T_c-273.15,'s-' ,'linewidth',2)
        plot(out.TS(i).x, out.TS(i).T_h-273.15,'o-' ,'linewidth',2)
        grid on
        xlabel('Heat fraction [-]','fontsize',14,'fontweight','bold')
        ylabel('Temperature [°C]','fontsize',14,'fontweight','bold')
        set(gca,'fontsize',14,'fontweight','bold')
        hold off
    end
end

end

%% NESTED FUNCTIONS

function res = hexseries_multiple_res(x, ub, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param, N_hex)
var = hexseries_multiple(x, ub, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param, N_hex);
res = norm(var.res);
end

function out = hexseries_multiple(x, ub, fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, param, N_hex)
x = x.*ub;
in_c_int(1) = in_c_su;

for i = 1:N_hex-1
    if strcmp(param(i).C.type, 'H')
        in_c_int(i+1) = in_c_int(i)+x(i)/m_dot_c;
    elseif strcmp(param(i).H.type, 'T')
        options = optimset('Display','off');
        in_c_int(i+1) = fsolve (@(x) Tcex_def(x, in_c_int(i),  P_c_su, m_dot_c, x, fluid_c), in_c_int(i),options);
    end
end

in_h_int(N_hex+1) = in_h_su;
for i = flip(1:N_hex)
    [out.int_hex(i), out.TS(i)] = HexModel2(fluid_h, P_h_su, in_h_int(i+1), m_dot_h, fluid_c, P_c_su, in_c_int(i), m_dot_c, param(i));
    if strcmp(param(i).H.type, 'H')
        in_h_int(i) = out.int_hex(i).H.h_ex;
    elseif strcmp(param(i).H.type, 'T')
        in_h_int(i) = out.int_hex(i).H.T_ex;
    end
end

for i = 1:N_hex-1
    if x(i) == 0 && out.int_hex(i).Q_dot_tot == 0
        out.res(i) = 0;
    elseif x(i) == 0 && out.int_hex(i).Q_dot_tot ~= 0
        out.res(i) = 1;
    elseif x(i) ~= 0 && out.int_hex(i).Q_dot_tot == 0
        out.res(i) = 1;
    else
        out.res(i) = 1 - x(i)/out.int_hex(i).Q_dot_tot;
    end
end
end

function err = Tcex_def(Tex, Tsu,  Psu, Mdot, Qdot, fluid)
err = Tex-(Tsu+Qdot/Mdot/sf_PropsSI_bar('C', Tsu, Tex, Psu, fluid));
end

function stop = outputfunFS(x, optimValues, state)
%disp(norm(optimValues.fval))
stop = norm(optimValues.fval) < 1e-6;
end