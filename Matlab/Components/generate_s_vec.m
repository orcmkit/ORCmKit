function s_vec = generate_s_vec(param, out, fluid)
s_vec = NaN*ones(1, length(out.T_vec));
if strcmp(param.type,'H') && not(strcmp(fluid(1:3), 'ICP'))  
    if param.solub
        H_ref_s = out.x_vec.*out.H_rv_vec + (1-out.x_vec).*out.H_rl_vec;
    else
        H_ref_s = out.H_vec;
    end
    for i = 1:length(H_ref_s)
        s_vec = CoolProp.PropsSI('S','P',out.P_vec(i),'H',H_ref_s(i),fluid);
    end
end
end