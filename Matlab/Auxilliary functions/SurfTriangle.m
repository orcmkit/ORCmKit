function S = SurfTriangle(X,Y)
S = 0.5*abs((X(2)-X(1))*(Y(3)-Y(1))-(X(3)-X(1))*(Y(2)-Y(1)));
end