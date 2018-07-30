function F = F_lmtd(LMTD_cc, R1, P1)

f = @(x) res_Flmtd(x, LMTD_cc, R1, P1);
F = fminbnd(f, 1e-5, 0.999999999);

end
function res = res_Flmtd(F_g, LMTD_cc, R1, P1)
LTMD_g = F_g*LMTD_cc;
NTU1_g = P1/LTMD_g;
a = 0.433;
b = 1.6;
c = 0.267;
d = 0.5;
F_n = (1/(1+a*R1^(d*b)*NTU1_g^b)^c);
res = F_g - F_n;
end