function [hConv, Nu, flag] = UD_HTC(m_dot, h_nom, m_dot_nom, n_nom) % VERIFIED
% User-define HTC
% RDickes - 23/07/2018 (rdickes@ulg.ac.be)


hConv = h_nom*(m_dot/m_dot_nom)^n_nom;
Nu = NaN;
flag = 1;

end