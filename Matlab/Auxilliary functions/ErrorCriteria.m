function [ae, are, mae, mare, nrmse, r2] = ErrorCriteria(Y_data, Y_model)
nan = unique([find(isnan(Y_data)),find(isnan(Y_model))]);
Y_data(nan) = [];
Y_model(nan) = [];
ae = AE(Y_data, Y_model);
are = ARE(Y_data, Y_model);
mae = MAE(Y_data, Y_model);
mare = MARE(Y_data, Y_model);
rmse = RMSE(Y_data, Y_model);
nrmse = NRMSE(Y_data, Y_model);

r2 = R2(Y_data, Y_model);
end

