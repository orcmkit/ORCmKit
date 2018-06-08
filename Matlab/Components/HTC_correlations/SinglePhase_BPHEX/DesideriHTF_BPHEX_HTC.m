function [hConv, Nu, flag] = DesideriHTF_BPHEX_HTC(mu, mu_rat, Pr, k, G, Dh, disp_flag) % verified

% Function implementing the correlation of Desideri for single flow of
% thermal oil in a BPHEX with a chevron angle of 65°
% source: Adriano's thesis page 62


% WARNING : NO INFLUENCE OF CHEVRON ANGLE !!! 

% Outputs:
%   hConv = convective heat transfer coefficient
%   Nu = Nusselt number
%   flag =  +1 if all is ok
%           0  if Re not ok
%           -1 if Pr not ok
%           -2 if Re and Pr not ok

% RDickes - 02/05/2018 (rdickes@ulg.ac.be)

Re_min = 31.6; Re_max = 205;
Pr_min = 76; Pr_max = 117;

flag = 1;

Re = G*Dh/mu;
Nu = 0.283*Re^0.8*Pr^0.3333333333333*mu_rat^0.14;
hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Desideri_HTF singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Desideri_HTF singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end

end