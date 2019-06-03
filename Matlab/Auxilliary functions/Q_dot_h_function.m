% % CODE 1
% delta_input_vec_h = [d_M_dot_r_ss(i,1)/1000  d_T_cd_su_ss(i,1) d_P_cd_su_ss(i,1)*1e5 d_T_cd_ex_ss(i,1) d_P_cd_ex_ss(i,1)*1e5]; % --> incertitude mesure
% input_vec_h = [m_dot_wf_ss(i,1)  T_cd_su_ss(i,1) P_cd_su_ss(i,1)*1e5 T_cd_ex_ss(i,1) P_cd_ex_ss(i,1)*1e5]; % mesure exp
% ub_h = input_vec_h + delta_input_vec_h; % valeur max possible
% lb_h = input_vec_h - delta_input_vec_h; % valeur min possible
% 
% % calcul Qdot min possible
% f_Qdot_cd_wf_min = @(x) Q_dot_h_function(x, ub_h, fluid_wf, 'Qdot'); % creation fct
% [input_vec_cd_wf_min(i,:), Q_dot_cd_wf_min_ss(i,1), flag_cd_wf_min(i,1)] = patternsearch(f_Qdot_cd_wf_min, input_vec_h./ub_h, [],[],[],[],lb_h./ub_h, ub_h./ub_h,options); % optimiseur
% input_vec_cd_wf_min(i,:) = input_vec_cd_wf_min(i,:).*ub_h;
% delta_input_cd_wf_min(i,:) = input_vec_h-input_vec_cd_wf_min(i,:);
% 
% % calcul Qdot max possible
% f_Qdot_cd_wf_max = @(x) Q_dot_h_function(x, ub_h, fluid_wf, 'inv_Qdot');
% [input_vec_cd_wf_max(i,:), inv_Q_dot_h_max, flag_cd_wf_max(i,1)] = patternsearch(f_Qdot_cd_wf_max, input_vec_h./ub_h, [],[],[],[],lb_h./ub_h, ub_h./ub_h,options);
% Q_dot_cd_wf_max_ss(i,1) = 1/inv_Q_dot_h_max;
% input_vec_cd_wf_max(i,:) = input_vec_cd_wf_max(i,:).*ub_h;
% delta_input_cd_wf_max(i,:) = input_vec_h-input_vec_cd_wf_max(i,:);

% CODE 2
function out = Q_dot_h_function(input_vec, ub, fluid, type)

input_vec = input_vec.*ub;
m_dot = input_vec(1);
T_su = input_vec(2);
P_su = input_vec(3);
T_ex = input_vec(4);
P_ex = input_vec(5);

if abs(T_su - (CoolProp.PropsSI('T', 'Q', 0.5, 'P', P_su, fluid)-273.15))<1e-3
    if strcmp(type, 'Qdot')
        h_su = CoolProp.PropsSI('H', 'Q', 0, 'P', P_su, fluid);
    elseif strcmp(type, 'inv_Qdot')
        h_su = CoolProp.PropsSI('H', 'Q', 1, 'P', P_su, fluid);
    end
else
    h_su = CoolProp.PropsSI('H', 'T', T_su + 273.15, 'P', P_su, fluid);
end

if abs(T_ex - (CoolProp.PropsSI('T', 'Q', 0.5, 'P', P_ex, fluid)-273.15))<1e-3
    if strcmp(type, 'Qdot')
        h_ex = CoolProp.PropsSI('H', 'Q', 1, 'P', P_ex, fluid);
    elseif strcmp(type, 'inv_Qdot')
        h_ex = CoolProp.PropsSI('H', 'Q', 0, 'P', P_ex, fluid);
    end
else
    h_ex = CoolProp.PropsSI('H', 'T', T_ex + 273.15, 'P', P_ex, fluid);
end
Q_dot_h = m_dot*(h_su-h_ex);

if strcmp(type, 'Qdot')
    out = Q_dot_h;
elseif strcmp(type, 'inv_Qdot')
    out = 1/Q_dot_h;
end

end

