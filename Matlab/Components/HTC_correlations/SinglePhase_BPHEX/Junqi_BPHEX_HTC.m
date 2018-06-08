function [hConv, Nu, flag] = Junqi_BPHEX_HTC(mu, Pr, k, G, Dh, theta, disp_flag) % verified

% Function implementing the correlation of Junqi for single flow in BPHEX
% and validated with R245fa under ORC operation with BPHEX~~ than in S2P
% Source: "Experimental investigation on heat transfer characteristics of
% plat heat exchanger applied in organic Rankine cycle (ORC)", Junqi et
% al., 2018

% Outputs:
%   hConv = convective heat transfer coefficient
%   Nu = Nusselt number
%   flag =  +1 if all is ok
%           0  if Re not ok
%           -1 if Pr not ok
%           -2 if Re and Pr not ok

% RDickes - 02/05/2018 (rdickes@ulg.ac.be)

Re_min = 180; Re_max = 7000;
Pr_min = 2; Pr_max = 12;

flag = 1;

Re = G*Dh/mu;
Nu = 0.964*Re^0.671*Pr^0.32*(theta/pi)^1.022;
hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Junqi singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Junqi singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end

end