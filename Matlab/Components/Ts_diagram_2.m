function [line_htf, line_ctf, line_orc] = Ts_diagram_2(TS, fluid, param)

%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 28/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% Ts_diagram is a single matlab code developed to plot T-s diagrams while 
% controlling the display properties (see the Documentation/Ts_diagram_MatlabDoc)
%
% The model inputs are:
%       - TS: a structure variable which contains the TS vectors of the different components in the ORC
%       - param: a structure variable which contains all the display parameters
%
% The model paramters provided in 'param' must include the following variables:
%       param.color_ctf = line color for CTF
%       param.LineStyle_ctf = line style for CTF
%       param.LineWidth_ctf = line width for CTF
%       param.MarkerType_ctf = marker type for CTF
%       param.alpha_ctf = transparency for CTF
%       param.MarkerSize_ctf = marker size for CTF
%       param.color_htf = line color for HTF
%       param.LineStyle_htf = line style for HTF
%       param.LineWidth_htf = line width for HTF
%       param.MarkerType_htf = marker type for HTF
%       param.alpha_htf = transparency for HTF
%       param.MarkerSize_htf = marker size for HTF
%       param.color_orc = line color for ORC
%       param.LineStyle_orc = line style for ORC
%       param.LineWidth_orc = line width for ORC
%       param.MarkerType_orc = marker type for ORC
%       param.alpha_orc = transparency for ORC
%       param.MarkerSize_orc = marker size for ORC
%
% The model outputs are:
%       - line_htf: graphical variable which includes properties of HTF display (useful for plotting legend)
%       - line_ctf: graphical variable which includes properties of CTF display (useful for plotting legend)
%       - line_orc: graphical variable which includes properties of ORC display (useful for plotting legend)
%
% See the documentation for further details or contact rdickes@ulg.ac.be

%% MODELLING CODE

if nargin == 2   
    % if no display information are probided by the user, default values
    % are given
    param.color_ctf = 'b';
    param.LineStyle_ctf = '--';
    param.LineWidth_ctf = 1;
    param.MarkerType_ctf = 'o';
    param.alpha_ctf = 1;
    param.MarkerSize_ctf = 5;
    param.color_htf = 'r';
    param.LineStyle_htf = '--';
    param.LineWidth_htf = 1;
    param.MarkerType_htf = 'o';
    param.alpha_htf = 1;
    param.MarkerSize_htf = 5;
    param.color_orc = 'k';
    param.LineStyle_orc = '-';
    param.LineWidth_orc = 1;
    param.MarkerType_orc = 'o';
    param.alpha_orc = 1;
    param.MarkerSize_orc = 5;
    
end

hold on
load(['TS_prop_',fluid , '.mat'])
plot(s_TS_curve, T_TS_curve-273.15, 'k', 'linewidth', param.LineWidth_orc+0.5)
% Cold heat transfer fluid lines 
if isfield(TS, 'CD')
    line_ctf =patchline(TS.CD.s_h, TS.CD.T_c-273.15, 'edgecolor',param.color_ctf, 'linestyle',param.LineStyle_ctf, 'linewidth',param.LineWidth_ctf,'edgealpha',param.alpha_ctf);
    patchline([TS.CD.s_h(1) TS.CD.s_h(end)], [TS.CD.T_c(1) TS.CD.T_c(end)]-273.15, 'edgecolor',param.color_ctf, 'linestyle','none', 'Marker', param.MarkerType_ctf,'edgealpha',param.alpha_ctf, 'MarkerSize', param.MarkerSize_ctf);
end
if isfield(TS, 'SUB')
    line_ctf = patchline(TS.SUB.s_h, TS.SUB.T_c-273.15, 'edgecolor',param.color_ctf, 'linestyle',param.LineStyle_ctf, 'linewidth',param.LineWidth_ctf,'edgealpha',param.alpha_ctf);
    patchline([TS.SUB.s_h(1) TS.SUB.s_h(end)], [TS.SUB.T_c(1) TS.SUB.T_c(end)]-273.15, 'edgecolor',param.color_ctf, 'linestyle','none', 'Marker', param.MarkerType_ctf,'edgealpha',param.alpha_ctf, 'MarkerSize', param.MarkerSize_ctf);
end

% Hot heat transfer fluid lines 
if isfield(TS, 'PRE')
    line_htf = patchline(TS.PRE.s_c, TS.PRE.T_h-273.15, 'edgecolor',param.color_htf, 'linestyle',param.LineStyle_htf, 'linewidth',param.LineWidth_htf,'edgealpha',param.alpha_htf);
    patchline([TS.PRE.s_c(1) TS.PRE.s_c(end)], [TS.PRE.T_h(1) TS.PRE.T_h(end)]-273.15, 'edgecolor',param.color_htf, 'linestyle','none', 'Marker', param.MarkerType_htf,'edgealpha',param.alpha_htf, 'MarkerSize', param.MarkerSize_htf);
end
if isfield(TS, 'EV')
    line_htf = patchline(TS.EV.s_c, TS.EV.T_h-273.15, 'edgecolor',param.color_htf, 'linestyle',param.LineStyle_htf, 'linewidth',param.LineWidth_htf,'edgealpha',param.alpha_htf);
    patchline([TS.EV.s_c(1) TS.EV.s_c(end)], [TS.EV.T_h(1) TS.EV.T_h(end)]-273.15, 'edgecolor',param.color_htf, 'linestyle','none', 'Marker', param.MarkerType_htf,'edgealpha',param.alpha_htf, 'MarkerSize', param.MarkerSize_htf);
end

% ORC lines
line_orc = patchline(TS.cycle.s, TS.cycle.T-273.15, 'edgecolor',param.color_orc, 'linestyle',param.LineStyle_orc, 'linewidth',param.LineWidth_orc,'edgealpha',param.alpha_orc);
if isfield(TS, 'PP')
    patchline([TS.PP.s(1) TS.PP.s(end)], [TS.PP.T(1) TS.PP.T(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end
if isfield(TS, 'REC')
    patchline([TS.REC.s_c(1) TS.REC.s_c(end)], [TS.REC.T_c(1) TS.REC.T_c(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end
if isfield(TS, 'PRE')
    patchline([TS.PRE.s_c(1) TS.PRE.s_c(end)], [TS.PRE.T_c(1) TS.PRE.T_c(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end
if isfield(TS, 'EV')
    patchline([TS.EV.s_c(1) TS.EV.s_c(end)], [TS.EV.T_c(1) TS.EV.T_c(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end
if isfield(TS, 'LossesHP')
    patchline([TS.LossesHP.s(1) TS.LossesHP.s(end)], [TS.LossesHP.T(1) TS.LossesHP.T(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end
if isfield(TS, 'EXP')
    patchline([TS.EXP.s(1) TS.EXP.s(end)], [TS.EXP.T(1) TS.EXP.T(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end
if isfield(TS, 'REC')
    if isfield(TS.REC, 's_h')
        patchline([TS.REC.s_h(1) TS.REC.s_h(end)], [TS.REC.T_h(1) TS.REC.T_h(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
    elseif isfield(TS.REC, 's_H')
        patchline([TS.REC.s_H(1) TS.REC.s_H(end)], [TS.REC.T_h(1) TS.REC.T_h(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
    end
end
if isfield(TS, 'CD')
    patchline([TS.CD.s_h(1) TS.CD.s_h(end)], [TS.CD.T_h(1) TS.CD.T_h(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end
if isfield(TS, 'SUB')
    patchline([TS.SUB.s_h(1) TS.SUB.s_h(end)], [TS.SUB.T_h(1) TS.SUB.T_h(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end
if isfield(TS, 'LossesLP')
    patchline([TS.LossesLP.s(1) TS.LossesLP.s(end)], [TS.LossesLP.T(1) TS.LossesLP.T(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end
if isfield(TS, 'LR')
    patchline([TS.LR.s(1) TS.LR.s(end)], [TS.LR.T(1) TS.LR.T(end)]-273.15, 'edgecolor',param.color_orc, 'linestyle','none', 'Marker', param.MarkerType_orc,'edgealpha',param.alpha_orc, 'MarkerSize', param.MarkerSize_orc);
end

hold off
end
