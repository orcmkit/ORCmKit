function MARE = MARE(Y_data, Y_model)
% R. Dickes - 04/23/2015
% General function computing the Mean Absolute Relative Error of a dataset
Y_data = Y_data(Y_data~=0);
Y_model = Y_model(Y_data~=0);
MARE = 1/length(Y_data)*sum(abs((Y_data-Y_model)./Y_data));

end