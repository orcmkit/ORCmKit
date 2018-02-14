clear all
close all
clc


%% COLD SIDE
T_su_c = 25;
T_ex_c = 40;
P_su_c = 12e5;
P_ex_c = 12e5;
fluid_c = 'water';
m_dot_c = 0.1;

DT_c = 0.1;
DP_c = 0.01e5;
DM_c = 0.01;
delta_input_vec_c = [DM_c DT_c DP_c DT_c DP_c];

input_vec_c = [m_dot_c  T_su_c P_su_c T_ex_c P_ex_c]; 

ub_c = input_vec_c + delta_input_vec_c;
lb_c = input_vec_c - delta_input_vec_c;


% Measured value 
Q_dot_c = Q_dot_c_function(input_vec_c./ub_c, ub_c, fluid_c, 'Qdot');

% Minimum possible value
f_Qdot_c_min = @(x) Q_dot_c_function(x, ub_c, fluid_c, 'Qdot');
[input_vec_c_min, Q_dot_c_min, flag_c_min] = patternsearch(f_Qdot_c_min, input_vec_c./ub_c, [],[],[],[],lb_c./ub_c, ub_c./ub_c);
input_vec_c_min = input_vec_c_min.*ub_c;
delta_input_c_min = input_vec_c-input_vec_c_min;

% Maximum possible value
f_Qdot_c_max = @(x) Q_dot_c_function(x, ub_c, fluid_c, 'inv_Qdot');
[input_vec_c_max, inv_Q_dot_c_max, flag_c_max] = patternsearch(f_Qdot_c_max, input_vec_c./ub_c, [],[],[],[],lb_c./ub_c, ub_c./ub_c);
Q_dot_c_max = 1/inv_Q_dot_c_max;
input_vec_c_max = input_vec_c_max.*ub_c;
delta_input_c_max = input_vec_c-input_vec_c_max;


%% HOT SIDE
T_su_h = 80;
T_ex_h = 50;
P_su_h = 10e5;
P_ex_h = 10e5;
fluid_h = 'water';
m_dot_h = 0.1;

DT_h = 0.1;
DP_h = 0.01e5;
DM_h = 0.01;
delta_input_vec_h = [DM_h DT_h DP_h DT_h DP_h];

input_vec_h = [m_dot_h  T_su_h P_su_h T_ex_h P_ex_h]; 

ub_h = input_vec_h + delta_input_vec_h;
lb_h = input_vec_h - delta_input_vec_h;


% Measured value 
Q_dot_h = Q_dot_h_function(input_vec_h./ub_h, ub_h, fluid_h, 'Qdot');

% Minimum possible value
f_Qdot_h_min = @(x) Q_dot_h_function(x, ub_h, fluid_h, 'Qdot');
[input_vec_h_min, Q_dot_h_min, flag_h_min] = patternsearch(f_Qdot_h_min, input_vec_h./ub_h, [],[],[],[],lb_h./ub_h, ub_h./ub_h);
input_vec_h_min = input_vec_h_min.*ub_h;
delta_input_h_min = input_vec_h-input_vec_h_min;

% Maximum possible value
f_Qdot_h_max = @(x) Q_dot_h_function(x, ub_h, fluid_h, 'inv_Qdot');
[input_vec_h_max, inv_Q_dot_h_max, flag_h_max] = patternsearch(f_Qdot_h_max, input_vec_h./ub_h, [],[],[],[],lb_h./ub_h, ub_h./ub_h);
Q_dot_h_max = 1/inv_Q_dot_h_max;
input_vec_h_max = input_vec_h_max.*ub_h;
delta_input_h_max = input_vec_h-input_vec_h_max;

%% PLOT 

figure;
hold on
% Measured value
plot(Q_dot_c, Q_dot_h,  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'Marker', 's', 'MarkerSize', 10, 'LineStyle', 'none');
vbar = errorbar(Q_dot_c, Q_dot_h, Q_dot_h-Q_dot_h_min, Q_dot_h_max-Q_dot_h);
hbar = herrorbar(Q_dot_c, Q_dot_h, Q_dot_c-Q_dot_c_min, Q_dot_c_max-Q_dot_c);

