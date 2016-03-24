function [DT_log, DTh, DTc ]= HEX_DTlog(Th_su, Th_ex, Tc_su, Tc_ex)
DTh = Th_su-Tc_ex;
DTc = Th_ex-Tc_su;
if DTh>0 && DTc>0
    if DTh ~= DTc;
        DT_log = (DTh-DTc)/log(DTh/DTc);
    else
        DT_log = DTh;
    end
else
    DT_log =1e-8;
end

end
