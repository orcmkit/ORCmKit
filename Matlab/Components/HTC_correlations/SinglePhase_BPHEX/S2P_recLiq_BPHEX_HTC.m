function [hConv, Nu, flag] = S2P_recLiq_BPHEX_HTC(mu, mu_rat, Pr, k, G, Dh, disp_flag) % VERIFIED

% Function implementing the correlation of Dickes for single-phase  flow in the
% recuperator of Sun2Power (BPHEX) --> FOR THE LIQUID PHASE !!!
% Source: Dickes' PhD thesis

% WARNING : NO INFLUENCE OF CHEVRON ANGLE !!! 


% RDickes - 23/07/2018 (rdickes@ulg.ac.be)
Re_min = 20; Re_max = 500;
Pr_min = 0.6; Pr_max = 6;

flag = 1;
Re = G*Dh/mu;
C = 0.853743743896484;
m = 0.292527580261230;
Nu = C*Re^m*Pr^(1/3)*mu_rat^0;

hConv = Nu*k/Dh;



if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['S2P_recLiq singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['S2P_recLiq singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end
end