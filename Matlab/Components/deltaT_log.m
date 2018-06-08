function [DT_log, DTh, DTc ]= deltaT_log(Th_su, Th_ex, Tc_su, Tc_ex)
% function that provides the mean logarithm temperature difference between two fluids
DTh = max(Th_su-Tc_ex,1e-2);
DTc = max(Th_ex-Tc_su,1e-2);
if DTh ~= DTc;
    DT_log = (DTh-DTc)/log(DTh/DTc);
else
    DT_log = DTh;
end
end
