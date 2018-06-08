function [hConv, Nu, flag] = Kim_BPHEX_HTC(mu, Pr, k, G, Dh, theta, disp_flag) % verified

% Function implementing the correlation of Kim for single flow of
% water in a BPHEX with impact of chevron angle
% source: Han's paper (both boiling and condensation)


% Outputs:
%   hConv = convective heat transfer coefficient
%   Nu = Nusselt number
%   flag =  +1 if all is ok
%           0  if Re not ok
%           -1 if Pr not ok
%           -2 if Re and Pr not ok

% RDickes - 02/05/2018 (rdickes@ulg.ac.be)

Re_min = -inf; Re_max = inf;
Pr_min = -inf; Pr_max = inf;

flag = 1;

Re = G*Dh/mu;
Nu = 0.295*Re^0.64*Pr^0.32*theta^0.09;
hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Kim singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Kim singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end

end