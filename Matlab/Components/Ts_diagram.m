function [line_htf, line_ctf, line_orc] = Ts_diagram(TS, display)
if nargin == 1
    display.color_ctf = 'b';
    display.LineStyle_ctf = '--';
    display.LineWidth_ctf = 1;
    display.MarkerType_ctf = 'o';
    display.alpha_ctf = 1;
    display.MarkerSize_ctf = 5;
    display.color_htf = 'r';
    display.LineStyle_htf = '--';
    display.LineWidth_htf = 1;
    display.MarkerType_htf = 'o';
    display.alpha_htf = 1;
    display.MarkerSize_htf = 5;
    display.color_orc = 'k';
    display.LineStyle_orc = '-';
    display.LineWidth_orc = 1;
    display.MarkerType_orc = 'o';
    display.alpha_orc = 1;
    display.MarkerSize_orc = 5;
    
end


% HCF lines
if isfield(TS, 'CD')
    line_ctf =patchline(TS.CD.s_h, TS.CD.T_c-273.15, 'edgecolor',display.color_ctf, 'linestyle',display.LineStyle_ctf, 'linewidth',display.LineWidth_ctf);
    patchline([TS.CD.s_h(1) TS.CD.s_h(end)], [TS.CD.T_c(1) TS.CD.T_c(end)]-273.15, 'edgecolor',display.color_ctf, 'linestyle','none', 'Marker', display.MarkerType_ctf,'edgealpha',display.alpha_ctf, 'MarkerSize', display.MarkerSize_ctf);
end
if isfield(TS, 'SUB')
    line_ctf = patchline(TS.SUB.s_h, TS.SUB.T_c-273.15, 'edgecolor',display.color_ctf, 'linestyle',display.LineStyle_ctf, 'linewidth',display.LineWidth_ctf,'edgealpha',display.alpha_ctf);
    patchline([TS.SUB.s_h(1) TS.SUB.s_h(end)], [TS.SUB.T_c(1) TS.SUB.T_c(end)]-273.15, 'edgecolor',display.color_ctf, 'linestyle','none', 'Marker', display.MarkerType_ctf,'edgealpha',display.alpha_ctf, 'MarkerSize', display.MarkerSize_ctf);
end

% HTF lines
if isfield(TS, 'PRE')
    line_htf = patchline(TS.PRE.s_c, TS.PRE.T_h-273.15, 'edgecolor',display.color_htf, 'linestyle',display.LineStyle_htf, 'linewidth',display.LineWidth_htf,'edgealpha',display.alpha_htf);
    patchline([TS.PRE.s_c(1) TS.PRE.s_c(end)], [TS.PRE.T_h(1) TS.PRE.T_h(end)]-273.15, 'edgecolor',display.color_htf, 'linestyle','none', 'Marker', display.MarkerType_htf,'edgealpha',display.alpha_htf, 'MarkerSize', display.MarkerSize_htf);
end
if isfield(TS, 'EV')
    line_htf = patchline(TS.EV.s_c, TS.EV.T_h-273.15, 'edgecolor',display.color_htf, 'linestyle',display.LineStyle_htf, 'linewidth',display.LineWidth_htf,'edgealpha',display.alpha_htf);
    patchline([TS.EV.s_c(1) TS.EV.s_c(end)], [TS.EV.T_h(1) TS.EV.T_h(end)]-273.15, 'edgecolor',display.color_htf, 'linestyle','none', 'Marker', display.MarkerType_htf,'edgealpha',display.alpha_htf, 'MarkerSize', display.MarkerSize_htf);
end

% ORC lines
line_orc = patchline(TS.cycle.s, TS.cycle.T-273.15, 'edgecolor',display.color_orc, 'linestyle',display.LineStyle_orc, 'linewidth',display.LineWidth_orc,'edgealpha',display.alpha_orc);
if isfield(TS, 'PP')
    patchline([TS.PP.s(1) TS.PP.s(end)], [TS.PP.T(1) TS.PP.T(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
end
if isfield(TS, 'REC')
    patchline([TS.REC.s_c(1) TS.REC.s_c(end)], [TS.REC.T_c(1) TS.REC.T_c(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
end
if isfield(TS, 'PRE')
    patchline([TS.PRE.s_c(1) TS.PRE.s_c(end)], [TS.PRE.T_c(1) TS.PRE.T_c(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
end
if isfield(TS, 'EV')
    patchline([TS.EV.s_c(1) TS.EV.s_c(end)], [TS.EV.T_c(1) TS.EV.T_c(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
end
if isfield(TS, 'DPHP')
    patchline([TS.DPHP.s(1) TS.DPHP.s(end)], [TS.DPHP.T(1) TS.DPHP.T(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
end
if isfield(TS, 'EXP')
    patchline([TS.EXP.s(1) TS.EXP.s(end)], [TS.EXP.T(1) TS.EXP.T(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
end
if isfield(TS, 'REC')
    if isfield(TS.REC, 's_h')
        patchline([TS.REC.s_h(1) TS.REC.s_h(end)], [TS.REC.T_h(1) TS.REC.T_h(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
    elseif isfield(TS.REC, 's_H')
        patchline([TS.REC.s_H(1) TS.REC.s_H(end)], [TS.REC.T_h(1) TS.REC.T_h(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
    end
end
if isfield(TS, 'CD')
    patchline([TS.CD.s_h(1) TS.CD.s_h(end)], [TS.CD.T_h(1) TS.CD.T_h(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
end
if isfield(TS, 'SUB')
    patchline([TS.SUB.s_h(1) TS.SUB.s_h(end)], [TS.SUB.T_h(1) TS.SUB.T_h(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
end
if isfield(TS, 'DPLP')
    patchline([TS.DPLP.s(1) TS.DPLP.s(end)], [TS.DPLP.T(1) TS.DPLP.T(end)]-273.15, 'edgecolor',display.color_orc, 'linestyle','none', 'Marker', display.MarkerType_orc,'edgealpha',display.alpha_orc, 'MarkerSize', display.MarkerSize_orc);
end

end
