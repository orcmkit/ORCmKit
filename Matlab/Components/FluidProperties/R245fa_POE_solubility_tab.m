function [ zeta_r, res_zeta_r] = R245fa_POE_solubility_tab(T_K, P_Pa, Tsat_pure_K, fit_DTP_zeta)
% RDickes - 11/04/2018
% Impletement correlations from Grebner, tuned with Zhelny's data
DT = T_K-Tsat_pure_K;
P_bar = P_Pa/1e5;
if not(isempty(fit_DTP_zeta))% && DT > 0 && DT < 140 && P_bar > 1 && P_bar < 20
zeta_r = fit_DTP_zeta([P_bar, DT]);
res_zeta_r = 0;
%    if isnan(zeta_r)
        %disp(['zeta_r = fit_DTP_zeta([' num2str(P_bar) ', ' num2str(DT) '])'])
        %[ zeta_r, res_zeta_r] = R245fa_POE_solubility_1(T_K, P_Pa, Tsat_pure_K);
%    end
else
    %disp(['zeta_r = fit_DTP_zeta([' num2str(P_bar) ', ' num2str(DT) '])'])
    [ zeta_r, res_zeta_r] = R245fa_POE_solubility_1(T_K, P_Pa, Tsat_pure_K);
end
end