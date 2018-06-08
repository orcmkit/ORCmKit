function RMSE = RMSE(Y_data, Y_model)
% R. Dickes - 04/23/2015
% General function computing the Root Mean Square Error of a dataset
vec_ok = find(not(isnan(Y_data)) & not(isnan(Y_model)));
RMSE = sqrt(sum((Y_data(vec_ok)-Y_model(vec_ok)).^2)/length(Y_data(vec_ok)));

end