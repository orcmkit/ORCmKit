function [DT_log, DTh, DTc ]= deltaT_log(Th_su, Th_ex, Tc_su, Tc_ex)

DTh = max(1e-12,Th_su-Tc_ex);
DTc = max(1e-12,Th_ex-Tc_su);
%if DTh>0 && DTc>0
    if DTh ~= DTc;
        DT_log = (DTh-DTc)/log(DTh/DTc);
    else
        DT_log = DTh;
    end
    DT_log = max(1e-12, DT_log);
end
%else
%    DT_log =1e-5;
%end
