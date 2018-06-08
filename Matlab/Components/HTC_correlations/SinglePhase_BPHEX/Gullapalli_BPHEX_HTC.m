function [hConv, Nu, flag] = Gullapalli_BPHEX_HTC(mu, Pr, k, G, Dh, theta, disp_flag) % verified

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
y = 0.4;
Re = G*Dh/mu;
C0 = 0.236853827;
C1 = -0.429914999;
C2 = 0.194375018;
C3 = 0.146176215;
C4 = -0.147253283; 
C5 = 0.236667683; 
C6 = 0.599681567; 
K0 = -0.284829132;
K1 = 0.696216086;
K2 = 1.120898786;
beta = theta;
Nu = ((C0+C1*beta+C2*beta^2)*Re^(K0 + K1*beta) + (C3+C4*beta+C5*beta^2)*Re^(K2*beta) + C6)*Pr^y;
hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Gullapalli singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Gullapalli singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end

end