function ARE = ARE(Y_data, Y_model)

ARE = abs((Y_data-Y_model)./Y_data);

end