function h_boiling = Desideri_Boiling_BPHEX_HTC(x, rho_l, rho_v, rho, mu_l, k_l, sigma, G, Dh)

% Boiling HTC correlation published by Desideri in "An experimental analysis 
% of flow boiling and pressure drop in a brazed plate heat exchanger for 
% organic Rankine cycle power systems"
% in International Journal of Heat and Mass Transfer, 2017

% RDickes - 03/05/2018

% INCONSISTENT FROM A NUMERICAL POINT OF VIEW !!!!!!!!!!!!!!!

g = 9.81; 
Bd = (rho_l-rho_v)*g*Dh^2/sigma; 
rho_star = rho_l/rho_v;
We = (G^2*Dh)/(rho*sigma);
Re_l = G*(1-x)*Dh/mu_l;
    

h_boiling = 1.48e3 * (We^-3.22e-2) * (rho_star^-3.38e-1) * (Re_l^4.51e-1) * (Bd^-4.69e-1) ;
%h_boiling = 2.06e3 * (We^-1.45e-2) * (rho_star^-3.62e-1) * (Re_l^4.14e-1) * (Bd^-4.87e-1) ;
