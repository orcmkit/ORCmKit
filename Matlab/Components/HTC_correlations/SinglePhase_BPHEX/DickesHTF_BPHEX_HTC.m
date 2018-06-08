function [hConv, Nu, flag] = DickesHTF_BPHEX_HTC(mu, mu_rat, Pr, k, G, Dh, disp_flag) % VERIFIED

% Function implementing the correlation of Dickes for HTFf low in the
% evaporator of Sun2Power (BPHEX)
% Source: Dickes' PhD thesis

% WARNING : NO INFLUENCE OF CHEVRON ANGLE !!! 


% RDickes - 06/06/2018 (rdickes@ulg.ac.be)
Re_min = 0; Re_max = 60;
Pr_min = 40; Pr_max = 500;

flag = 1;

Re = G*Dh/mu;

Nu_lam = 4.5506;
if Re < 5
    Nu = Nu_lam;
else
    X_cor_emb = 1.0e+03 *[   2.236317626953125  -0.004662027740479   0.000216403198242  -0.069780941009521  -0.002640230712891   0.000809321624756];
    C = max(1e-10,  X_cor_emb(1)*Re^X_cor_emb(2)+X_cor_emb(3));
    m = max(0,      X_cor_emb(4)*Re^X_cor_emb(5)+X_cor_emb(6));

%     if Re>= 5 && Re<10
%         C  = 0.641254937648773; m = 0.290755736827850;
%     elseif Re>= 10 && Re<20
%         C  = 0.304503300786018; m = 0.610132354497910;
%     elseif Re>= 20 && Re<40
%         C  = 0.197615483403206; m = 0.819999603033066;
%     elseif Re>= 40
%         C = 0.230211752653122; m = 0.785629661977291;
%     end
    Nu = max(Nu_lam,C*Re^m*Pr^(1/3)*mu_rat^0.14);
end
hConv = Nu*k/Dh;



if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Dickes singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Dickes singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end
end