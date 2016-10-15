function [DT_log, DTh, DTc ]= deltaT_log(Th_su, Th_ex, Tc_su, Tc_ex)
tol_T = 1e-2;
DTh = Th_su-Tc_ex;
DTc = Th_ex-Tc_su;
if DTh>-tol_T && DTc>-tol_T
    if DTh ~= DTc;
        DT_log = (DTh-DTc)/log(abs(DTh/DTc));
    else
        DT_log = DTh;
    end
else
    DT_log =1e-8;
end
