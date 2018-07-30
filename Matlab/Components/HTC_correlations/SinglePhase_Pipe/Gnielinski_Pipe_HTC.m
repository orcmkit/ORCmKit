function [hConv, Nu, flag] =  Gnielinski_Pipe_HTC(mu, Pr, k, G, Dh, L, disp_flag) % VERIFIED
% Function implementing the correlation of Gnielinski for single flow in a
% pipe
% Reference: VDI, section G1, page 696

% RDickes - 26/04/2018 (rdickes@ulg.ac.be)

Re_min = 0; Re_max = 1e6;

flag = 1;

Re = G*Dh/mu;
if Re > 1e4 %fully turbulent
    Pr_min = 0.1; Pr_max = 1000;
    Nu = Gnielinski_turbulent(Re, Pr);
elseif Re< 2300 % fully laminar
    Pr_min = 0.6; Pr_max = inf;
    Nu = Gnielinski_laminar(Re, Pr, Dh, L);
else % transition region
    Pr_min = 0.1; Pr_max = 1000;
    gamma = (Re-2300)/(1e4-2300);
    Nu_lam2300 =  Gnielinski_laminar(2300, Pr, Dh, L);
    Nu_turb10000 = Gnielinski_turbulent(1e4, Pr);
    Nu = (1-gamma)*Nu_lam2300 + gamma*Nu_turb10000;
end

hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['Gnielinski singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if Pr >= Pr_max || Pr <= Pr_min
    if disp_flag
        display(['Gnielinski singe-phase: Out of validity range --> Pr = ' num2str(Pr) ' is out of [' num2str(Pr_min) ' - ' num2str(Pr_max) '] !!!'])
    end
    flag = flag -2;
end

end

function Nu = Gnielinski_laminar(Re, Pr, Dh, L) %source: VDI section G1 - 3.2.1, page 695, for constant heat flux
Nu_1 = 4.364;
Nu_2 = 1.953*(Re*Pr*Dh/L)^0.33333333333333333333333333333;
Nu = (Nu_1^3 + 0.6^3 + (Nu_2-0.6)^3)^0.3333333333333333333333;
end

function Nu = Gnielinski_turbulent(Re, Pr) %source: VDI section G1 - 4.1, page 696
f = (1.8*log10(Re)-1.5)^-2; %Konakov correlation
Nu = ((f/8)*(Re-1000)*Pr)/(1+12.7*sqrt(f/8)*(Pr^(2/3)-1)); % Gnielinski
end


