function [hConv, Nu, flag] = S2P_rec1phase_BPHEX_HTC(mu, mu_rat, Pr, k, G, Dh, disp_flag) % VERIFIED

% Function implementing the correlation of Dickes for R245fa flow in the
% recuperator of Sun2Power (BPHEX) --> Martin + transition regime at very
% low reynolds
% Source: Dickes' PhD thesis

% WARNING : NO INFLUENCE OF CHEVRON ANGLE !!! 


% RDickes - 06/06/2018 (rdickes@ulg.ac.be)
Re_min = 0; Re_max = 10000;
Pr_min = 0.2; Pr_max = 6;

flag = 1;

Re = G*Dh/mu;

Nu_min = 0.5;
X_cor_emb = [0.199999935913086 -83.100854656982506  -1.253950653076172   0.739998950195312]; %a little less good fit of martin
%X_cor_emb = 1.0e+02 *[0.001853997215271  -4.279656236389161 -0.016638780212402   0.007199992553711]; %really good fit of Martin
C = X_cor_emb(1);
m = max(0, X_cor_emb(2)*Re^X_cor_emb(3)+X_cor_emb(4));

if Pr > 1
    n_murat = 0.16666666667;
else
    n_murat = 0;
end
Nu = max(Nu_min,C*Re^m*Pr^(1/3)*mu_rat^n_murat);
hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Dickes singe-phase BPHEX: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Dickes singe-phase BPHEX: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end
end