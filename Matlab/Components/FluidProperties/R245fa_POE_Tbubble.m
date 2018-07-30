function Tbubble_K = R245fa_POE_Tbubble(zeta_r, P_Pa, Tsat_pure_K)
% RDickes - 11/04/2018
% Impletement correlations from Grebner, tuned with Zhelny's data 
theta = R245fa_POE_theta(zeta_r, P_Pa);
Tbubble_K = (1+theta)*Tsat_pure_K;
end
