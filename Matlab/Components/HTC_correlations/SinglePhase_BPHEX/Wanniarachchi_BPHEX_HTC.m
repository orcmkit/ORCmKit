function [hConv, Nu, flag] = Wanniarachchi_BPHEX_HTC(mu, mu_rat, Pr, k, G, Dh, theta,phi, disp_flag) % verified

% Function implementing the correlation of Wanniarachchi for single flow in BPHEX
% Source: "Approximate correlations for chevron-type plate heat exchangers", Wanniarachchi et al., 1995

% Outputs:
%   hConv = convective heat transfer coefficient
%   Nu = Nusselt number
%   flag =  +1 if all is ok
%           0  if Re not ok
%           -1 if Pr not ok
%           -2 if Re and Pr not ok

% RDickes - 26/04/2018 (rdickes@ulg.ac.be)

Re_min = 1; Re_max = 1e4;
Pr_min = -inf; Pr_max = inf;

flag = 1;

theta_bis = min(max(20,90-theta*180/pi),62);

Re = G*Dh/mu;
m = 0.646+0.00111*theta_bis;
j_Nu_t = 12.6*theta_bis^(-1.142)*phi^(1-m)*Re^(m);
j_Nu_l = 3.65*theta_bis^(-0.455)*phi^0.661*Re^(0.339);
Nu = (j_Nu_l^3 + j_Nu_t^3)^(1/3)*Pr^(1/3)*mu_rat^(0.17);
hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Wanniarachchi singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Wanniarachchi singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end


end