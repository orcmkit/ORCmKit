function [F, flag] = F_lmtd2(R, P)
% R = (T_h_i - T_h_ex)/(T_c_ex - T_c_i) --> must be positive [0; infinte [
% P = (T_c_ex - T_c_i)/(T_h_i - T_c_i) --> must be within 0 and 1
R_min = 1e-4;%0.2;
R_max = 100;%4;
R = max(R_min,min(R_max,R));
if R < 0.41
    P_max_for_given_R = 0.99999;
else
    P_max_for_given_R = 0.995*((0.966475661367996)*R.^3 + (-1.431274296912407)*R.^2 + (0.247230065033875)*R + (0.309118607270897)) / (R.^4 + (-1.766309686745371)*R.^3 + (1.287179055148762)*R.^2 + (-0.902512766020191)*R + (0.484880205333508));
end

P = max(0,  min(P_max_for_given_R,P));

if R <= 1 % Cdot_min = Cdot_h
    epsilon = P;
    Cr = R;
    Pbis = P;
    Rbis = R;
else % Cdot_max = Cdot_c
    epsilon = P*R;
    Cr = 1/R;
    Pbis = P*R;
    Rbis = 1/R;
end
f = @(x) res_NTU_crossflow(x, epsilon, Cr);
options = optimset('Display','off');
[NTU, res_NTU,flag_ntu] = fzero(f, 1, options);

if abs(res_NTU)< 1e-2  && flag_ntu > 0
    if  R==1
        F=P/NTU/(1-P);
    else
        F = epsilon*log((Pbis-1)/(Rbis*Pbis -1))/NTU/Pbis/(Rbis-1);
    end
    flag = flag_ntu;
else
    %disp(['Prob compute NTU with F_lmtd2(' num2str(R) ',' num2str(P) ')'])
    F = 1;
    flag = flag_ntu;
end
end
function res = res_NTU_crossflow(Ntu, epsilon, C_r)
Ntu = max(1e-4, Ntu);
epsilon_g=1-exp((1/C_r)*Ntu^0.22*(exp(-C_r*Ntu^0.78)-1));
res = (epsilon-epsilon_g);
end