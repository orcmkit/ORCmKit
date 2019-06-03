function out = HEX_profile_Solub_DP2_Fast(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, param, h_h_l_in, h_h_v_in, h_c_l_in, h_c_v_in, DP_h_in, DP_c_in)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 27/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% HEX_profile is a single matlab code aiming to calculate the temperature
% profiles occuring between two media if the heat power is provided as
% input. This code has been developed to be as general as possible and can
% be used for multi-phase heat transfer for both hot and cold fluids and can
% also handle incompressible or supercritical medium. The code implements
% in an automatic algorithm the cells division methodology proposed by Bell et al. $
% in:
% "A generalized moving-boundary algorithm to predict the heat transfer rate
% of counterflow heat exchangers for any phase configuration,”
% Appl. Therm. Eng., vol. 79, pp. 192–201, 2015.
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
%       - Q_dot: heat power transferred between the fluids          [W]
%       - param: structure variable containing the model parameters i.e.
%           *param.H.type = type of input for hot fluid ('H' for enthalpy,'T' for temperature)
%           *param.C.type = type of input for cold fluid ('H' for enthalpy,'T' for temperature)

% The model outputs are:
%       - out: a structure variable which includes at miniumum the following information:
%               - x_vec =  vector of power fraction in each zone
%               - Qdot_vec =  vector of heat power in each zone [W]
%               - H_h_vec = HF enthalpy vector                  [J/kg]
%               - H_c_vec = CF enthalpy vector                  [J/kg]
%               - T_h_vec = HF temperature vector               [K]
%               - T_c_vec = CF temperature vector               [K]
%               - DT_vec = Temperature difference vector        [K]
%
% See the documentation for further details or contact rdickes@ulg.ac.be

%% MODELLING CODE
decim = 4; % degree of accuracy used for comparing the entalpies

% Cell division for hot fluid (create a vector of enthalpy the different zones on the hot fluid side)
h_h_su = in_h_su;
h_h_ex = h_h_su - Q_dot/m_dot_h;
if not(param.H.DP_on) || DP_h_in <1e-2
    h_h_l = h_h_l_in;
    h_h_v = h_h_v_in;
    P_h_ex = P_h_su;
else
    P_h_ex = P_h_su - DP_h_in;
    [h_h_l, h_h_v, ~, ~, ~, ~] = find_2P_boundaries_Fast(fluid_h, h_h_su, h_h_ex, P_h_su, P_h_ex, param.H);
end
H_h_vec_disc = sort([linspace(h_h_ex, h_h_su, param.n_disc) h_h_l h_h_v]);
H_h_vec = unique(H_h_vec_disc(not(H_h_vec_disc<h_h_ex | H_h_vec_disc>h_h_su)));
out.H.h_v = h_h_v;
out.H.h_l = h_h_l;

% Cell division for cold fluid (create a vector of enthalpy the different zones on the cold fluid side)
h_c_su = in_c_su;
h_c_ex = h_c_su + Q_dot/m_dot_c;
if not(param.C.DP_on) || DP_c_in <1e-2
    h_c_l = h_c_l_in;
    h_c_v = h_c_v_in;
    P_c_ex = P_c_su;
else
    P_c_ex = P_c_su - DP_c_in;
    [h_c_l, h_c_v, ~, ~, ~, ~] = find_2P_boundaries_Fast(fluid_c, h_c_su, h_c_ex, P_c_su, P_c_ex, param.C);    
end
H_c_vec_disc = sort([linspace(h_c_su, h_c_ex, param.n_disc) h_c_l h_c_v]);
H_c_vec = unique(H_c_vec_disc(not(H_c_vec_disc>h_c_ex | H_c_vec_disc<h_c_su)));
out.C.h_v = h_c_v;
out.C.h_l = h_c_l;

% Cell divitions for entire heat exchanger
j = 1;
while  j < max(length(H_h_vec),length(H_c_vec))-1
    
    Q_dot_h = m_dot_h*(H_h_vec(j+1)-H_h_vec(j));
    Q_dot_c = m_dot_c*(H_c_vec(j+1)-H_c_vec(j));
    if round(Q_dot_h,decim) > round(Q_dot_c,decim)
        H_h_vec = [H_h_vec(1:j), H_h_vec(j)+Q_dot_c/m_dot_h, H_h_vec(j+1:end)];
    elseif round(Q_dot_h,decim) < round(Q_dot_c,decim)
        H_c_vec = [H_c_vec(1:j), H_c_vec(j)+Q_dot_h/m_dot_c, H_c_vec(j+1:end)];
    end
    j = j+1;
end

Q_dot_vec = m_dot_h*diff(H_h_vec);
out.x_flux = (H_h_vec-H_h_vec(1))./(H_h_vec(end)-H_h_vec(1));
out.Qdot_vec = Q_dot_vec;
out.C.H_vec = H_c_vec;
if (H_c_vec(end)- H_c_vec(1))>0
    X_c = (H_c_vec-H_c_vec(1))/(H_c_vec(end)- H_c_vec(1));
    out.C.P_vec = (1-X_c)*P_c_su + X_c*P_c_ex;
else
    out.C.P_vec = linspace(P_c_su, P_c_ex, length(H_c_vec));
end
out.C.T_vec = NaN*ones(1,length(out.C.H_vec));
out.C.Tbubble_min_vec = NaN*ones(1,length(out.C.H_vec));
out.C.Tsat_pure_vec = NaN*ones(1,length(out.C.H_vec));
out.C.H_rl_vec = NaN*ones(1,length(out.C.H_vec));
out.C.H_oil_vec = NaN*ones(1,length(out.C.H_vec));
out.C.H_rv_vec = NaN*ones(1,length(out.C.H_vec));
out.C.zeta_r_vec = NaN*ones(1,length(out.C.H_vec));
out.C.x_vec = NaN*ones(1,length(out.C.H_vec));
out.C.C_rl_vec = NaN*ones(1,length(out.C.H_vec));
out.C.C_rv_vec = NaN*ones(1,length(out.C.H_vec));
for k = 1:length(out.C.H_vec)
    if param.C.solub
        out.C.Tsat_pure_vec(k) = CoolPropt_T_PQ(param.C.abs_lowLevel,param.C.CP_file, out.C.P_vec(k), 0.5);
        out.C.Tbubble_min_vec(k) = R245fa_POE_Tbubble(1-param.C.C_oil, out.C.P_vec(k), out.C.Tsat_pure_vec(k));
        [out.C.T_vec(k), out.C.H_rl_vec(k), out.C.H_oil_vec(k), out.C.H_rv_vec(k), out.C.zeta_r_vec(k), out.C.x_vec(k), out.C.C_rl_vec(k), out.C.C_rv_vec(k)] = HP_solubMixt_Fast(param.C.C_oil, out.C.P_vec(k), out.C.H_vec(k), fluid_c, param.C.fluid_lub, out.C.Tbubble_min_vec(k), out.C.Tsat_pure_vec(k), param.C.fit_DTP_zeta, param.C);
    else
        if strcmp(fluid_c(1:3), 'ICP')
            out.C.T_vec(k) = PropsSI_ICP('T', 'H', out.C.H_vec(k), 'P', out.C.P_vec(k), fluid_c);
        else
            out.C.T_vec(k) = CoolPropt_T_HP(param.C.abs_lowLevel,param.C.CP_file, out.C.H_vec(k), out.C.P_vec(k)); %CoolProp.PropsSI('T', 'H', out.C.H_vec(k), 'P', out.C.P_vec(k), fluid_c);
            out.C.Tsat_pure_vec(k) =  CoolPropt_T_PQ(param.C.abs_lowLevel,param.C.CP_file, out.C.P_vec(k), 0.5);   %CoolProp.PropsSI('T', 'Q', 0.5, 'P', out.C.P_vec(k), fluid_c);
            h_c_l = CoolPropt_H_PQ(param.C.abs_lowLevel,param.C.CP_file, out.C.P_vec(k), 0); %CoolProp.PropsSI('H', 'Q', 0, 'P', out.C.P_vec(k), fluid_c);
            h_c_v = CoolPropt_H_PQ(param.C.abs_lowLevel,param.C.CP_file, out.C.P_vec(k), 1); %CoolProp.PropsSI('H', 'Q', 1, 'P', out.C.P_vec(k), fluid_c);
            if out.C.H_vec(k)<h_c_l
                out.C.x_vec(k) = -2;
            elseif out.C.H_vec(k)>h_c_v
                out.C.x_vec(k) = 2;
            else 
                out.C.x_vec(k) = min(1,max(0,(out.C.H_vec(k) - h_c_l)/(h_c_v - h_c_l)));
            end
        end
    end
end

out.H.H_vec = H_h_vec;
if (H_h_vec(end)- H_h_vec(1))>0
    X_h = (H_h_vec-H_h_vec(1))/(H_h_vec(end)- H_h_vec(1));
    out.H.P_vec = (1-X_h)*P_h_ex + X_h*P_h_su;
else
    out.H.P_vec = linspace(P_h_ex, P_h_su, length(H_h_vec));
end
out.H.T_vec = NaN*ones(1,length(out.H.H_vec));
out.H.Tbubble_min_vec = NaN*ones(1,length(out.H.H_vec));
out.H.Tsat_pure_vec = NaN*ones(1,length(out.H.H_vec));
out.H.H_rl_vec = NaN*ones(1,length(out.H.H_vec));
out.H.H_oil_vec = NaN*ones(1,length(out.H.H_vec));
out.H.H_rv_vec = NaN*ones(1,length(out.H.H_vec));
out.H.zeta_r_vec = NaN*ones(1,length(out.H.H_vec));
out.H.x_vec = NaN*ones(1,length(out.H.H_vec));
out.H.C_rl_vec = NaN*ones(1,length(out.H.H_vec));
out.H.C_rv_vec = NaN*ones(1,length(out.H.H_vec));
for k = 1:length(out.H.H_vec)
    if param.H.solub
        out.H.Tsat_pure_vec(k) = CoolPropt_T_PQ(param.H.abs_lowLevel,param.H.CP_file, out.H.P_vec(k), 0.5);
        out.H.Tbubble_min_vec(k) = R245fa_POE_Tbubble(1-param.H.C_oil, out.H.P_vec(k), out.H.Tsat_pure_vec(k));
        [out.H.T_vec(k), out.H.H_rl_vec(k), out.H.H_oil_vec(k), out.H.H_rv_vec(k), out.H.zeta_r_vec(k), out.H.x_vec(k), out.H.C_rl_vec(k), out.H.C_rv_vec(k)] = HP_solubMixt_Fast(param.H.C_oil, out.H.P_vec(k), out.H.H_vec(k), fluid_h, param.H.fluid_lub, out.H.Tbubble_min_vec(k), out.H.Tsat_pure_vec(k), param.H.fit_DTP_zeta, param.H);
    else
        if strcmp(fluid_h(1:3), 'ICP')
            out.H.T_vec(k) = PropsSI_ICP('T', 'H', out.H.H_vec(k), 'P', out.H.P_vec(k), fluid_h);
        else
            out.H.T_vec(k) = CoolPropt_T_HP(param.H.abs_lowLevel,param.H.CP_file, out.H.H_vec(k), out.H.P_vec(k)); %CoolProp.PropsSI('T', 'H', out.H.H_vec(k), 'P', out.H.P_vec(k), fluid_h);
            out.H.Tsat_pure_vec(k) = CoolPropt_T_PQ(param.H.abs_lowLevel,param.H.CP_file, out.H.P_vec(k), 0.5);   %CoolProp.PropsSI('T', 'Q', 0.5, 'P', out.H.P_vec(k), fluid_h);
            h_h_l = CoolPropt_H_PQ(param.H.abs_lowLevel,param.H.CP_file, out.H.P_vec(k), 0); %CoolProp.PropsSI('H', 'Q', 0, 'P', out.H.P_vec(k), fluid_h);
            h_h_v = CoolPropt_H_PQ(param.H.abs_lowLevel,param.H.CP_file, out.H.P_vec(k), 1); %CoolProp.PropsSI('H', 'Q', 1, 'P', out.H.P_vec(k), fluid_h);
            if out.H.H_vec(k)<h_h_l
                out.H.x_vec(k) = -2;
            elseif out.H.H_vec(k)>h_h_v
                out.H.x_vec(k) = 2;
            else
                out.H.x_vec(k) = min(1,max(0,(out.H.H_vec(k) - h_h_l)/(h_h_v - h_h_l)));
            end
        end
    end
end
out.DT_vec = out.H.T_vec-out.C.T_vec;
out.pinch = min(out.H.T_vec-out.C.T_vec);
out.DT_log_rec = [];
out.A_prim = [];

end
