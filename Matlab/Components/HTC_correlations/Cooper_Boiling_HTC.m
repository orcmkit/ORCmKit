function [h_boiling, Nu, flag] = Cooper_Boiling_HTC(Dh, k_l, P, P_crit, M, h_conv_h, DT_log, Qdot) % ---> VERIFIED

% Boiling HTC correlation published by Cooper in "Saturation Nucleate Pool Boiling - a Simple Correlation"
% in The Institution of Chemical Engineers Symposium Series, 1984

% M is a relative molecular weight, which is relatively the same than *1000

% RDickes - 16/08/2018


C = 1;
flag = 1;
AU_tp = Qdot/DT_log;
p_star = P/P_crit;
Rp = 0.4; % roughness in µm
q_0 = 1;
f_q = @(xx) res_iter_Cooper_boiling(xx, C, p_star, Rp, M, h_conv_h, AU_tp, Qdot);
opts = optimoptions('fsolve', 'display', 'none');
[q,err_q] = fsolve(f_q, q_0, opts);
[~, h, ~, ~, ~] = iter_Cooper_boiling(q, C, p_star, Rp, M, h_conv_h, AU_tp, Qdot);
if err_q > 5e-2
    flag = -5;
    if disp_flag
        display('Cooper boiling: Wrong heat flux')
    end
end
h_boiling = h;
Nu = h*Dh/k_l;

end


function res_q = res_iter_Cooper_boiling(q, C, p_star, Rp, M, h_conv_h, AU_tp, Qdot)
[res_q, ~, ~, ~, ~] = iter_Cooper_boiling(q, C, p_star, Rp, M, h_conv_h, AU_tp, Qdot);
end

function [res_q, h, U, A_tp, q_new] = iter_Cooper_boiling(q, C, p_star, Rp, M, h_conv_h, AU_tp, Qdot)
%h = C*55*(p_star^(0.12-0.2*log10(Rp)))*((-log10(p_star))^(-0.55))*(q^0.67)*(M^(-0.5));
h = C*35*(p_star^(0.12))*((-log10(p_star))^(-0.55))*(q^0.67)*(M^(-0.5));
U = (1/h +  1/h_conv_h)^-1;
A_tp = AU_tp/U;
q_new = Qdot/A_tp;
res_q = abs(q_new-q)/q;
end