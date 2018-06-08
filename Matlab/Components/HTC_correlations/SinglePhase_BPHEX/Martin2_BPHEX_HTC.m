function [hConv, Nu, flag] = Martin2_BPHEX_HTC(mu, Pr, k, G, Dh, theta, disp_flag)

% Function implementing the correlaiton of Martin for single flow in BPHEX
% Source: "A theoretical approach to predict the performance of
% chevron-type plate heat exchangers", Holger Matrin, 1996

% 2nd version: coefficients provided in Pr Ngendakuman slides, but not
% found elsewhere :\

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
    f_0 = 16/Re;
    f_90 = 149.25/Re+0.9625;
else
    f_0 = (1.56*log(Re)-3)^(-2);
    f_90 = 9.75/Re^(0.289);
end
f = (((cos(theta))/sqrt(0.045*tan(theta) + 0.09*sin(theta) + f_0/cos(theta)))+((1-cos(theta))/(sqrt(3.8*f_90))))^(-2);
Nu = 0.205*(Pr^0.33333333)*(f*Re^2*sin(2*theta))^(0.374);
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
