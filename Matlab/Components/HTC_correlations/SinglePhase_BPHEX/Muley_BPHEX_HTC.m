function [hConv, Nu, flag] = Muley_BPHEX_HTC(mu, mu_rat, Pr, k, G, Dh, theta, phi, L, disp_flag) % verified

% Function implementing the correlation of Muley and Manglik for single flow 
% in BPHEX:

% Source of laminar flow: "",  
% Muley et al., 1999

% Source of turbulent flow: "Experimental Study of Turbulent Flow Heat 
% Transfer and Pressure Drop in a Plate Heat Exchanger With Chevron Plates",  
% Muley et al., 1999

% Warning, in Muley's correlation, the hydraulic diameter for computing the
% Reynold and the Nusselt numbers is 2*b and not 2*b/pi as generally used!

% Outputs:
%   hConv = convective heat transfer coefficient
%   Nu = Nusselt number
%   flag =  +1 if all is ok
%           0  if Re not ok
%           -1 if Pr not ok
%           -2 if Re and Pr not ok

% RDickes - 27/04/2018 (rdickes@ulg.ac.be)



Re_min = 2; Re_max = inf;

flag = 1;

Re = G*(Dh*phi)/mu;

theta_deg = theta*180/pi;
C1 = (0.2668 - 0.006967*theta_deg + 7.244e-5*theta_deg^2)*(20.78-50.94*phi+41.16*phi^2-10.51*phi^3);
C2 = (0.728+0.0543*sin(4*theta + 3.7));


if Re <400 % laminar flow correlation
    Pr_min = 130; Pr_max = 290;
    Nu = 1.6774*(Dh*phi/L)^0.3333333*(theta_deg/30)^0.38*Re^0.5*Pr^(1/3)*mu_rat^0.14;
elseif Re >1e3 % turbulent flow correlation
    Pr_min = -inf; Pr_max = inf;   
    Nu = C1*Re^C2*Pr^(1/3)*mu_rat^0.14;
else
    Pr_min = 130; Pr_max = 290;
    gamma = (Re-400)/(1e3-400);
    Nu_lam400 =  1.6774*(Dh*phi/L)^0.3333333*(theta_deg/30)^0.38*400^0.5*Pr^(1/3)*mu_rat^0.14;
    Nu_turb1000 = C1*1000^C2*Pr^(1/3)*mu_rat^0.14;
    Nu = (1-gamma)*Nu_lam400 + gamma*Nu_turb1000;
end
hConv = Nu*k/(Dh*phi);
if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Muley singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Muley singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end

end