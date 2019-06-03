function NRMSE = NRMSE(Y_data, Y_model)
% R. Dickes - 02/19/2016
% General function computing the Root Mean Square Error of a dataset
vec_ok = find(not(isnan(Y_data)) & not(isnan(Y_model)));
NRMSE = sqrt(sum((Y_data(vec_ok)-Y_model(vec_ok)).^2)/length(Y_data(vec_ok)))/(max(Y_data(vec_ok))-min(Y_data(vec_ok)));
end