function [hConv, Nu, flag] = S2P_recVap_BPHEX_HTC(mu, mu_rat, Pr, k, G, Dh, disp_flag) % VERIFIED

% Function implementing the correlation of Dickes for single-phase  flow in the
% recuperator of Sun2Power (BPHEX) --> FOR THE VAPOUR PHASE !!!
% Source: Dickes' PhD thesis

% WARNING : NO INFLUENCE OF CHEVRON ANGLE !!! 


% RDickes - 23/07/2018 (rdickes@ulg.ac.be)
Re_min = 1000; Re_max = 7000;
Pr_min = 0.5; Pr_max = 0.7;

flag = 1;
Re = G*Dh/mu;
C = 0.015585708618164;
m = 0.899999618530274;
Nu = C*Re^m*Pr^(0.499999364217122)*mu_rat^0;

hConv = Nu*k/Dh;



if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['S2P_recVap singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['S2P_recVap singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end
end