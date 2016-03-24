function NRMSE = NRMSE(Y_data, Y_model)
% R. Dickes - 02/19/2016
% General function computing the Root Mean Square Error of a dataset

NRMSE = sqrt(sum((Y_data-Y_model).^2)/length(Y_data))/(max(Y_data)-min(Y_data))*100;

end