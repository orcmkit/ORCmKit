function out = FCT_Hex_hConvVar(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, param)

% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems
%
% Remi Dickes - 11/05/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% "FCT_Hex_hConvVar.m" is a constitutive subfunction of "HexModel_hConvVar.m"
% It evaluates a multi-zone heat exchanger model with variable heat transfer coefficients
% while assuming a given heat transfer "Q_dot".

out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, param); %evaluate the temperature profile for a given heat power, cfr documentation of HEX_profile
[out.hConv_h, out.hConv_c, out.DTlog, out.eff_h, out.eff_c, out.A_h, out.A_c, out.U] = deal(NaN*ones(1,length(out.H_h_vec)-1));
for j = 1:length(out.T_h_vec)-1
    % Hot side heat transfer coefficient
    if strcmp(param.type_h, 'H')
        if isempty(strfind(fluid_h, 'INCOMP:'))
            if (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) < CoolProp.PropsSI('H','P',P_h_su,'Q',0,fluid_h)
                out.hConv_h(j) = param.hConv_h_liq_n*(m_dot_h/param.m_dot_h_n)^param.n;
            elseif (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) > CoolProp.PropsSI('H','P',P_h_su,'Q',1,fluid_h)
                out.hConv_h(j) = param.hConv_h_vap_n*(m_dot_h/param.m_dot_h_n)^param.n;
            else
                out.hConv_h(j) = param.hConv_h_tp_n*(m_dot_h/param.m_dot_h_n)^param.n;
            end
        else
            out.hConv_h(j) = param.hConv_h_liq_n*(m_dot_h/param.m_dot_h_n)^param.n;
        end
    elseif strcmp(param.type_h, 'T')
        out.hConv_h(j) = param.hConv_h_liq_n*(m_dot_h/param.m_dot_h_n)^param.n;
    end
    
    % Cold side heat transfer coefficient
    if strcmp(param.type_c, 'H')
        if isempty(strfind(fluid_c, 'INCOMP:'))
            if (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) < CoolProp.PropsSI('H','P',P_c_su,'Q',0,fluid_c)
                out.hConv_c(j) = param.hConv_c_liq_n*(m_dot_c/param.m_dot_c_n)^param.n;
            elseif (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) > CoolProp.PropsSI('H','P',P_c_su,'Q',1,fluid_c)
                out.hConv_c(j) = param.hConv_c_vap_n*(m_dot_c/param.m_dot_c_n)^param.n;
            else
                out.hConv_c(j) = param.hConv_c_tp_n*(m_dot_c/param.m_dot_c_n)^param.n;
            end
        else
            out.hConv_c(j) = param.hConv_c_liq_n*(m_dot_c/param.m_dot_c_n)^param.n;
        end
    elseif strcmp(param.type_c, 'T')
        out.hConv_c(j) = param.hConv_c_liq_n*(m_dot_c/param.m_dot_c_n)^param.n;
    end
    out.DTlog(j) = deltaT_log(out.T_h_vec(j+1), out.T_h_vec(j),out.T_c_vec(j), out.T_c_vec(j+1));
    
    % Hot side heat transfer efficiency (in case of fins)
    if not(isfield(param, 'fin_h'))
        out.eff_h(j) = 1;
    elseif strcmp(param.fin_h, 'none')
        out.eff_h(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_h(j), param.fin_h.k, param.fin_h.th, param.fin_h.r, param.fin_h.B, param.fin_h.H);
        out.eff_h(j) = 1-param.fin_h.omega_f*(1-eta_eff);
    end
    
    % Cold side heat transfer efficiency (in case of fins)
    if not(isfield(param, 'fin_c'))
        out.eff_c(j) = 1;
    elseif strcmp(param.fin_c, 'none')
        out.eff_c(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_c(j), param.fin_c.k, param.fin_c.th, param.fin_c.r, param.fin_c.B, param.fin_c.H);
        out.eff_c(j) = 1-param.fin_c.omega_f*(1-eta_eff);
    end
    
    % Global heat transfer coefficient and zone surface area
    out.U(j) = (1/out.hConv_h(j)/out.eff_h(j) + 1/out.hConv_c(j)/out.eff_c(j)/(param.A_c_tot/param.A_h_tot))^-1;
    out.A_h(j) = out.Qdot_vec(j)/out.DTlog(j)/out.U(j);
    out.A_c(j) = out.A_h(j)*param.A_c_tot/param.A_h_tot;
end
out.A_h_tot = sum(out.A_h);
out.resA = 1 - out.A_h_tot/param.A_h_tot;
end
