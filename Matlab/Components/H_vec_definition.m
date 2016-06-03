function H_vec = H_vec_definition(fluid, m_dot, P, h_in, Q_dot, port_in)

decim = 3; % degree of accuracy used for comparing the entalpies

if strcmp(port_in,'hot')
    h_h = h_in;
    h_c = h_h - Q_dot/m_dot;
elseif strcmp(port_in,'cold')
    h_c = h_in;
    h_h = h_c + Q_dot/m_dot;
end

if not(isempty(strfind(fluid, 'INCOMP:')))
    H_vec = [h_c  h_h];
else
    h_l = CoolProp.PropsSI('H','P',P,'Q',0, fluid);
    h_v = CoolProp.PropsSI('H','P',P,'Q',1, fluid);
    
    if round(h_v,decim) < round(h_c,decim)
        H_vec = [h_c  h_h]; %if vapour-phase only
    elseif round(h_l,decim) > round(h_h,decim)
        H_vec = [h_c  h_h]; % if liquid-phase only
    elseif (round(h_l,decim) < round(h_c,decim)) && (round(h_v,decim) > round(h_h,decim))
        H_vec = [h_c  h_h]; % if two-phase only
    elseif (round(h_l,decim) < round(h_h,decim)) && (round(h_l,decim) > round(h_c,decim))
        if (round(h_v,decim) < round(h_h,decim))
            H_vec = [h_c  h_l  h_v  h_h]; % if liquid, two phase and vapour
        else
            H_vec = [h_c  h_l  h_h]; % if liquid and two phase
        end
    elseif (round(h_v,decim) > round(h_c,decim)) && (round(h_v,decim) < round(h_h,decim))
        if (round(h_l,decim) > round(h_c,decim))
            H_vec = [h_c  h_l  h_v  h_h]; % if liquid, two phase and vapour
        else
            H_vec = [h_c  h_v  h_h]; % of twp phase and vapour
        end
    else
        H_vec = [h_c  h_h]; %in all other cases (should never happen)
    end
end
end