function  [hConv, Nu, flag] = Wang_finnedTubes_staggered_HTC(mu, Pr, k, G, Dh, Dc, Pt, Fp, N, disp_flag)

% air-side HTC proposed by Wang et al. in "Heat transfer and friction
% characteristics of plain fin-and- tube heat exchangers , part II : Correlation" 

% Dc is the externel diameter of the tubes forming the bank (tube diameter +2*fin thickness)
% Dh is the actual hydraulic diameter = 4*(A_c/A_o)*L with L : longitudinal depth of the heat exchanger / Ac = minimal flow area/ Ao = total surface area
% Fp is the fin pitch
% Pl is the tube longitudinal pitch

% RDickes - 25/07/2018


Re_Dc = G*Dc/mu;
p1 = -0.361-0.042*N/log(Re_Dc)+0.158*log(N*(Fp/Dc)^0.41);
p2 = -1.224-((0.076*(Pt/Dh)^1.42)/(log(Re_Dc)));
p3 = -0.083+0.058*N/log(Re_Dc);
p4 = -5.735 +1.21*log(Re_Dc/N);
p5 = -0.93;
j = 0.086*(Re_Dc^p1)*(N^p2)*((Fp/Dc)^p3)*((Fp/Dh)^p4)*((Fp/Pt)^p5);
Nu = j*Re_Dc*Pr^0.33333333333333333333;
hConv = Nu*k/Dc;

flag = 1;

end