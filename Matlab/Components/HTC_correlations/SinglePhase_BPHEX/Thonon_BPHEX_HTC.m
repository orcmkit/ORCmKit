function [hConv, Nu, flag] = Thonon_BPHEX_HTC(mu, Pr, k, G, Dh, theta, disp_flag) % VERIFIED

% Function implementing the correlation of Thonon for single flow in BPHEX
% Source: Ayub, did not managed to get the original paper, no idea about
% restriction on the Prandtl number

% RDickes - 26/04/2018 (rdickes@ulg.ac.be)
Re_min = 50; Re_max = 15000;
Pr_min = -inf; Pr_max = inf;

flag = 1;

Re = G*Dh/mu;

theta_data = ([15     30      45      60])*pi/180;
C_data =    [  0.1    0.2267  0.2998  0.2946];
m_data =    [  0.687  0.631   0.645   0.7];
C = interp1(theta_data,C_data,theta, 'linear', 'extrap');
m = interp1(theta_data,m_data,theta, 'linear', 'extrap');
Nu = C*Re^m*Pr^0.33333333;
hConv = Nu*k/Dh;



if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Thonon singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Thonon singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end
end