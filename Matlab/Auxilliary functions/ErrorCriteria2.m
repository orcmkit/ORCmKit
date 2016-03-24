function [ae, are, mae, mare, rmse, r2, nrmse] = ErrorCriteria2(Y_data, Y_model)
nan = unique([find(isnan(Y_data)),find(isnan(Y_model))]);
Y_data(nan) = [];
Y_model(nan) = [];
ae = AE(Y_data, Y_model);
are = 100*ARE(Y_data, Y_model);
mae = MAE(Y_data, Y_model);
mare = 100*MARE(Y_data, Y_model);
rmse = RMSE(Y_data, Y_model);
r2 = 100*R2(Y_data, Y_model);
nrmse = NRMSE(Y_data, Y_model);
end

