function MARE = MARE(Y_data, Y_model)
% R. Dickes - 04/23/2015
% General function computing the Mean Absolute Relative Error of a dataset
vec_ok = find(not(isnan(Y_data)) & not(isnan(Y_model)) & not(Y_data ==0));
Y_data = Y_data(vec_ok);
Y_model = Y_model(vec_ok);
MARE = 1/length(Y_data)*sum(abs((Y_data-Y_model)./Y_data));

end