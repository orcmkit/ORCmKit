function R2 = R2(Y_data, Y_model)
% R. Dickes - 04/23/2015
% General function computing the coefficient of deterlibation R2 of a dataset

Y_mean_data = sum(Y_data)/length(Y_data);
SS_tot = sum((Y_model-Y_mean_data).^2);
SS_res = sum((Y_model-Y_data).^2);
R2 = 1-SS_res/SS_tot;

end

