function [hConv, Nu, flag]  = Martin1_BPHEX_HTC(mu, mu_rat, Pr, k, G, Dh, theta, disp_flag) % Verified

% Function implementing the correlation of Martin for single flow in BPHEX
% Source: "A theoretical approach to predict the performance of
% chevron-type plate heat exchangers", Holger Matrin, 1996

% 1st version: coefficients provided in VDI (N6,page 1516), Martin's original paper, and in Pr. Ngendakumana slides

% no information regarding the validity restiction

% Outputs:
%   hConv = convective heat transfer coefficient
%   Nu = Nusselt number
%   flag =  +1 if all is ok
%           0  if Re not ok
%           -1 if Pr not ok
%           -2 if Re and Pr not ok

% RDickes - 26/04/2018 (rdickes@ulg.ac.be)

Re_min = -inf; Re_max = inf;
Pr_min = -inf; Pr_max = inf;

flag = 1;

                
Re = G*Dh/mu;
if Re < 2000
    f_0 = 64/Re;
    f_90 = 597/Re+3.85;
else
    f_0 = (1.8*log10(Re)-1.5)^(-2);
    f_90 =39/Re^(0.289);
end
f = (((cos(theta))/sqrt(0.18*tan(theta) + 0.36*sin(theta) + f_0/cos(theta)))+((1-cos(theta))/(sqrt(3.8*f_90))))^(-2);
Nu = 0.122*(Pr^0.33333333)*(mu_rat^0.16666666667)*(f*Re^2*sin(2*theta))^(0.374);
hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Martin singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Martin singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end


end
