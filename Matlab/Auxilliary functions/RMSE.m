function RMSE = RMSE(Y_data, Y_model)
% R. Dickes - 04/23/2015
% General function computing the Root Mean Square Error of a dataset

RMSE = sqrt(sum((Y_data-Y_model).^2)/length(Y_data));

end