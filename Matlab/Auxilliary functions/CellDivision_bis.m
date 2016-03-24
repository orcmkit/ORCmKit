function [H_h_vec, H_c_vec, x] = HEXprofile(fluid_h, m_dot_h, P_h_su, h_h_su, fluid_c, m_dot_c, P_c_su, h_c_su, Q_dot)

h_h_ex = h_h_su - Q_dot/m_dot_h;
h_h_l = CoolProp.PropsSI('H','P',P_h_su,'Q',0, fluid_h);
h_h_v = CoolProp.PropsSI('H','P',P_h_su,'Q',1, fluid_h);
if (h_h_v < h_h_su) && (h_h_v > h_h_ex)
    if (h_h_l < h_h_su) && (h_h_l > h_h_ex)
        H_h_vec = [h_h_ex  h_h_l  h_h_v  h_h_su];
    else
        H_h_vec = [h_h_ex  h_h_v  h_h_su];
    end
else
    H_h_vec = [h_h_ex h_h_su];
end
h_c_ex = h_c_su + Q_dot/m_dot_c;
h_c_l = CoolProp.PropsSI('H','P',P_c_su,'Q',0, fluid_c);
h_c_v = CoolProp.PropsSI('H','P',P_c_su,'Q',1, fluid_c);
if (h_c_l < h_c_ex) && (h_c_l > h_c_su)
    if (h_c_v < h_c_ex) && (h_c_v > h_h_su)
        H_c_vec = [h_c_su  h_c_l  h_c_v  h_c_ex];
    else
        H_c_vec = [h_c_su  h_c_l  h_c_ex];
    end
else
    H_c_vec = [h_c_su  h_c_ex];
end

j = 1;
while  j < max(length(H_h_vec),length(H_c_vec))    
    Q_dot_h = m_dot_h*(H_h_vec(j+1)-H_h_vec(j));
    Q_dot_c = m_dot_c*(H_c_vec(j+1)-H_c_vec(j));
    if round(Q_dot_h-Q_dot_c, 3) > 0
        H_h_vec = [H_h_vec(1:j), H_h_vec(j)+Q_dot_c/m_dot_h, H_h_vec(j+1:end)];      
    elseif round(Q_dot_h - Q_dot_c, 3) < 0
        H_c_vec = [H_c_vec(1:j), H_c_vec(j)+Q_dot_h/m_dot_c, H_c_vec(j+1:end)];   
    end
    j = j+1;
end


x = (H_h_vec-H_h_vec(1))./(H_h_vec(end)-H_h_vec(1));

end