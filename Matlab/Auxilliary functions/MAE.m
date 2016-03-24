function MAE = MAE(Y_data, Y_model)
% R. Dickes - 04/23/2015
% General function computing the Mean Absolute  Error of a dataset

MAE = 1/length(Y_data)*sum(abs(Y_data-Y_model));

end
