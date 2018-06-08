function [hConv, Nu, flag]  = Heavner_BPHEX_HTC(mu, mu_rat, Pr, k, G, Dh, theta, phi, disp_flag) % verified

% Function implementing the correlation of Heavner for single flow in BPHEX
% Source: Ayub, original paper not found

% Outputs:
%   hConv = convective heat transfer coefficient
%   Nu = Nusselt number
%   flag =  +1 if all is ok
%           0  if Re not ok
%           -1 if Pr not ok
%           -2 if Re and Pr not ok

% RDickes - 27/04/2018 (rdickes@ulg.ac.be)

Re_min = 400; Re_max = 10000;
Pr_min = 3.3; Pr_max = 5.9;

flag = 1;

Re = G*Dh/mu;
theta_data = (90-[67     56      33.5      45       22.5])*pi/180;
C_data =    	 [0.089  0.118   0.308     0.195    0.278];
m_data =         [0.718  0.720   0.667     0.692    0.683];
C = interp1(theta_data,C_data,theta);
m = interp1(theta_data,m_data,theta);
Nu = C*phi^(1-m)*Re^m*Pr^0.5*mu_rat^0.166667;
hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Heavner singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Heavner singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end

end