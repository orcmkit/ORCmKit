function out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
%change
decim = 3;
if strcmp(info.type_h,'H') && strcmp(info.type_c,'H')
    h_h_su = in_h_su;
    h_h_ex = h_h_su - Q_dot/m_dot_h;
    if not(isempty(strfind(fluid_h, 'INCOMP:')))
        h_h_l = h_h_su*3;
        h_h_v = h_h_su*3;
    else
        h_h_l = CoolProp.PropsSI('H','P',P_h_su,'Q',0, fluid_h);
        h_h_v = CoolProp.PropsSI('H','P',P_h_su,'Q',1, fluid_h);
    end 
    if round(h_h_v,decim) < round(h_h_ex,decim)
        H_h_vec = [h_h_ex  h_h_su];
    elseif round(h_h_l,decim) > round(h_h_su,decim)
        H_h_vec = [h_h_ex  h_h_su];
    elseif (round(h_h_l,decim) < round(h_h_ex,decim)) && (round(h_h_v,decim) > round(h_h_su,decim))
        H_h_vec = [h_h_ex  h_h_su];
    elseif (round(h_h_l,decim) < round(h_h_su,decim)) && (round(h_h_l,decim) > round(h_h_ex,decim))
        if (round(h_h_v,decim) < round(h_h_su,decim))
            H_h_vec = [h_h_ex  h_h_l  h_h_v  h_h_su];
        else
            H_h_vec = [h_h_ex  h_h_l  h_h_su];
        end
    elseif (round(h_h_v,decim) > round(h_h_ex,decim)) && (round(h_h_v,decim) < round(h_h_su,decim))
        if (round(h_h_l,decim) > round(h_h_ex,decim))
            H_h_vec = [h_h_ex  h_h_l  h_h_v  h_h_su];
        else
            H_h_vec = [h_h_ex  h_h_v  h_h_su];
        end
    else
        H_h_vec = [h_h_ex  h_h_su];
    end
    
    h_c_su = in_c_su;
    h_c_ex = h_c_su + Q_dot/m_dot_c;
    if not(isempty(strfind(fluid_c, 'INCOMP:')))
        h_c_l = h_c_ex*3;
        h_c_v = h_c_ex*3;
    else
        h_c_l = CoolProp.PropsSI('H','P',P_c_su,'Q',0, fluid_c);
        h_c_v = CoolProp.PropsSI('H','P',P_c_su,'Q',1, fluid_c);
    end
    if round(h_c_v,decim) < round(h_c_su,decim)
        H_c_vec = [h_c_su  h_c_ex];
    elseif round(h_c_l,decim) > round(h_c_ex,decim)
        H_c_vec = [h_c_su  h_c_ex];
    elseif (round(h_c_l,decim) < round(h_c_su,decim)) && (round(h_c_v,decim) > round(h_c_ex,decim))
        H_c_vec = [h_c_su  h_c_ex];
    elseif (round(h_c_l,decim) > round(h_c_su,decim)) && (round(h_c_l,decim) < round(h_c_ex,decim))
        if (round(h_c_v,decim) < round(h_c_ex,decim))
            H_c_vec = [h_c_su  h_c_l  h_c_v  h_c_ex];
        else
            H_c_vec = [h_c_su  h_c_l  h_c_ex];
        end
    elseif (round(h_c_v,decim) > round(h_c_su,decim)) && (round(h_c_v,decim) < round(h_c_ex,decim))
        if (round(h_c_l,decim) > round(h_c_su,decim))
            H_c_vec = [h_c_su  h_c_l  h_c_v  h_c_ex];
        else
            H_c_vec = [h_c_su  h_c_v  h_c_ex];
        end
    else
        H_c_vec = [h_c_su  h_c_ex];
    end

    j = 1;
    while  j < max(length(H_h_vec),length(H_c_vec))-1
        Q_dot_h = m_dot_h*(H_h_vec(j+1)-H_h_vec(j));
        Q_dot_c = m_dot_c*(H_c_vec(j+1)-H_c_vec(j));

        if round(Q_dot_h,30) > round(Q_dot_c,30) 
            H_h_vec = [H_h_vec(1:j), H_h_vec(j)+Q_dot_c/m_dot_h, H_h_vec(j+1:end)];
        elseif round(Q_dot_h,30) < round(Q_dot_c,30) 
            H_c_vec = [H_c_vec(1:j), H_c_vec(j)+Q_dot_h/m_dot_c, H_c_vec(j+1:end)];
        end
        j = j+1;
    end 
    Q_dot_vec = m_dot_h*diff(H_h_vec);
    x = (H_h_vec-H_h_vec(1))./(H_h_vec(end)-H_h_vec(1));
    out.x_vec = x;
    out.Qdot_vec = Q_dot_vec;
    out.H_h_vec = H_h_vec;
    out.H_c_vec = H_c_vec;
    out.T_h_vec = NaN*ones(1,length(H_h_vec));
    out.T_c_vec = NaN*ones(1,length(H_c_vec));
    for k = 1:length(H_c_vec)
        out.T_h_vec(k) = CoolProp.PropsSI('T','P', P_h_su, 'H', H_h_vec(k), fluid_h);
        out.T_c_vec(k) = CoolProp.PropsSI('T','P', P_c_su, 'H', H_c_vec(k), fluid_c);        
    end
    out.DT_vec = out.T_h_vec-out.T_c_vec;
    out.pinch = min(out.T_h_vec-out.T_c_vec);
    

elseif strcmp(info.type_h,'T') && strcmp(info.type_c,'H')
    h_c_su = in_c_su;
    h_c_ex = h_c_su + Q_dot/m_dot_c;
    if not(isempty(strfind(fluid_c, 'INCOMP:')))
        h_c_l = h_c_ex*3;
        h_c_v = h_c_ex*3;
    else
        h_c_l = CoolProp.PropsSI('H','P',P_c_su,'Q',0, fluid_c);
        h_c_v = CoolProp.PropsSI('H','P',P_c_su,'Q',1, fluid_c);
    end
    
    if round(h_c_v,decim) < round(h_c_su,decim)
        H_c_vec = [h_c_su  h_c_ex];
    elseif round(h_c_l,decim) > round(h_c_ex,decim)
        H_c_vec = [h_c_su  h_c_ex];
    elseif (round(h_c_l,decim) < round(h_c_su,decim)) && (round(h_c_v,decim) > round(h_c_ex,decim))
        H_c_vec = [h_c_su  h_c_ex];
    elseif (round(h_c_l,decim) > round(h_c_su,decim)) && (round(h_c_l,decim) < round(h_c_ex,decim))
        if (round(h_c_v,decim) < round(h_c_ex,decim))
            H_c_vec = [h_c_su  h_c_l  h_c_v  h_c_ex];
        else
            H_c_vec = [h_c_su  h_c_l  h_c_ex];
        end
    elseif (round(h_c_v,decim) > round(h_c_su,decim)) && (round(h_c_v,decim) < round(h_c_ex,decim))
        if (round(h_c_l,decim) > round(h_c_su,decim))
            H_c_vec = [h_c_su  h_c_l  h_c_v  h_c_ex];
        else
            H_c_vec = [h_c_su  h_c_v  h_c_ex];
        end
    end
      
    T_h_su = in_h_su;
    f_T_h_ex = @(x) Thex_def(x, T_h_su,  P_h_su, m_dot_h, Q_dot, fluid_h);
    lb =  CoolProp.PropsSI('T','P',P_c_su,'H',h_c_su, fluid_c);
    ub = T_h_su;
    T_h_ex = zeroBrent ( lb, ub, 1e-6, 1e-6, f_T_h_ex );
    cp_h = sf_PropsSI_bar('C', T_h_su, T_h_ex, P_h_su, fluid_h);
    Q_dot_vec = NaN*ones(1,length(H_c_vec)-1);
    T_h_vec = NaN*ones(1,length(H_c_vec)-1);
    T_h_vec(1) = T_h_ex;

    j =1;
    while  j <length(H_c_vec)
        Q_dot_vec(j) = m_dot_c*(H_c_vec(j+1)-H_c_vec(j));
        T_h_vec(j+1) = T_h_vec(j) + Q_dot_vec(j)/m_dot_h/cp_h;
        j=j+1;
    end
    x = (H_c_vec-H_c_vec(1))./(H_c_vec(end)-H_c_vec(1));
    out.x_vec = x;
    out.Qdot_vec = Q_dot_vec;
    out.H_h_vec = NaN;
    out.H_c_vec = H_c_vec;
    out.T_h_vec = T_h_vec;
    out.T_c_vec = NaN*ones(1,length(H_c_vec));
    for k = 1:length(H_c_vec)
        out.T_c_vec(k) = CoolProp.PropsSI('T','P', P_c_su, 'H', H_c_vec(k), fluid_c);        
    end
    out.DT_vec = out.T_h_vec-out.T_c_vec;
    out.pinch = min(out.T_h_vec-out.T_c_vec);
    
    
elseif strcmp(info.type_h,'H') && strcmp(info.type_c,'T')
    h_h_su = in_h_su;
    h_h_ex = h_h_su - Q_dot/m_dot_h;
    if not(isempty(strfind(fluid_h, 'INCOMP:')))
        h_h_l = h_h_su*3;
        h_h_v = h_h_su*3;
    else
        h_h_l = CoolProp.PropsSI('H','P',P_h_su,'Q',0, fluid_h);
        h_h_v = CoolProp.PropsSI('H','P',P_h_su,'Q',1, fluid_h);
    end 
    if round(h_h_v,decim) < round(h_h_ex,decim)
        H_h_vec = [h_h_ex  h_h_su];
    elseif round(h_h_l,decim) > round(h_h_su,decim)
        H_h_vec = [h_h_ex  h_h_su];
    elseif (round(h_h_l,decim) < round(h_h_ex,decim)) && (round(h_h_v,decim) > round(h_h_su,decim))
        H_h_vec = [h_h_ex  h_h_su];
    elseif (round(h_h_l,decim) < round(h_h_su,decim)) && (round(h_h_l,decim) > round(h_h_ex,decim))
        if (round(h_h_v,decim) < round(h_h_su,decim))
            H_h_vec = [h_h_ex  h_h_l  h_h_v  h_h_su];
        else
            H_h_vec = [h_h_ex  h_h_l  h_h_su];
        end
    elseif (round(h_h_v,decim) > round(h_h_ex,decim)) && (round(h_h_v,decim) < round(h_h_su,decim))
        if (round(h_h_l,decim) > round(h_h_ex,decim))
            H_h_vec = [h_h_ex  h_h_l  h_h_v  h_h_su];
        else
            H_h_vec = [h_h_ex  h_h_v  h_h_su];
        end
    end

    T_c_su = in_c_su;
    options = optimset('Display','off');
    T_c_ex = fsolve (@(x) Tcex_def(x, T_c_su,  P_c_su, m_dot_c, Q_dot, fluid_c), T_c_su,options);
    cp_c = sf_PropsSI_bar('C', T_c_su, T_c_ex, P_c_su, fluid_c);
    T_c_vec = [T_c_su T_c_ex];
    Q_dot_vec = NaN*ones(1,length(H_h_vec)-1);
    j = 1;
    while  j < length(H_h_vec)-1
        Q_dot_h = m_dot_h*(H_h_vec(j+1)-H_h_vec(j));
        T_c_vec = [T_c_vec(1:j), T_c_vec(j)+Q_dot_h/m_dot_c/cp_c, T_c_vec(j+1:end)];
        Q_dot_vec(j) = Q_dot_h;
        j = j+1;
    end
    x = (H_h_vec-H_h_vec(1))./(H_h_vec(end)-H_h_vec(1));
    out.x_vec = x;
    out.Qdot_vec = Q_dot_vec;
    out.H_h_vec = H_h_vec;
    out.H_c_vec = NaN;
    out.T_c_vec = T_c_vec;
    out.T_h_vec = NaN*ones(1,length(H_h_vec));
    for k = 1:length(H_h_vec)
        out.T_h_vec(k) = CoolProp.PropsSI('T','P', P_h_su, 'H', H_h_vec(k), fluid_h);        
    end
    out.DT_vec = out.T_h_vec-out.T_c_vec;
    out.pinch = min(out.T_h_vec-out.T_c_vec);
    
elseif strcmp(info.type_h,'T') && strcmp(info.type_c,'T')  
    T_c_su = in_c_su;
    options = optimset('Display','off');
    T_c_ex = fsolve (@(x) Tcex_def(x, T_c_su,  P_c_su, m_dot_c, Q_dot, fluid_c), T_c_su,options);
    out.T_c_vec = [T_c_su T_c_ex];
    out.H_c_vec = NaN;   
    T_h_su = in_h_su;
    options = optimset('Display','off');
    T_h_ex = fsolve (@(x) Thex_def(x, T_h_su,  P_h_su, m_dot_h, Q_dot, fluid_h), T_h_su,options);
    out.T_h_vec = [T_h_ex T_h_su]; 
    out.H_h_vec = NaN;
    out.Qdot_vec = Q_dot;
    out.x_vec = [0 1];
    out.DT_vec = out.T_h_vec-out.T_c_vec;
    out.pinch = min(out.T_h_vec-out.T_c_vec);
end
end

function err = Tcex_def(Tex, Tsu,  Psu, Mdot, Qdot, fluid)
err = Tex-(Tsu+Qdot/Mdot/sf_PropsSI_bar('C', Tsu, Tex, Psu, fluid));
end

function err = Thex_def(Tex, Tsu,  Psu, Mdot, Qdot, fluid)
err = Tex-(Tsu-Qdot/Mdot/sf_PropsSI_bar('C', Tsu, Tex, Psu, fluid));
end