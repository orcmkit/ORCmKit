function [hConv, Nu, flag] = DittusBoelter_Pipe_HTC(mu, Pr, k, G, Dh, n, disp_flag) % VERIFIED

% Function implementing the correlation of Dittus-Boelter for single flow
% in pipes
% Source: Incroprera and DeWitt,page 412

% Outputs:
%   hConv = convective heat transfer coefficient
%   Nu = Nusselt number
%   flag =  +1 if all is ok
%           0  if Re not ok
%           -1 if Pr not ok
%           -2 if Re and Pr not ok

% RDickes - 26/04/2018 (rdickes@ulg.ac.be)

Re_min = -inf; Re_max = inf;

flag = 1;

Re = G*Dh/mu;

if Re <2300 % laminar flow correlation
    Pr_min = 0.6; Pr_max = inf;
    Nu_lam = 4.36;
    Nu = Nu_lam;
elseif Re >1e4 % turbulent flow correlation --> Dittus Boelter
    Pr_min = 0.6; Pr_max = 160;
    Nu_turb = 0.023*Re^(4/5)*Pr^n;
    Nu = Nu_turb;
else
    Pr_min = 0.6; Pr_max = inf;
    Nu_lam2300 = 4.36;
    Nu_turb10000 = 0.023*(1e4)^(4/5)*Pr^n;
    gamma = (Re-2300)*(1e4-2300);
    Nu = (1-gamma)*Nu_lam2300 + gamma*Nu_turb10000;
end

hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['DittusBoelter singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['DittusBoelter singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end


end