function Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, info)
%% CODE DESCRIPTION
% ORCmKit - an open-source modelling library for ORC systems

% Remi Dickes - 27/04/2016 (University of Liege, Thermodynamics Laboratory)
% rdickes @ulg.ac.be
%
% HEX_Qdotmax is a single matlab code aiming to calculate the maxium amount of
% heat power that can be transferred between two fluids in a counterflow
% heat exchanger. This maxium is either given by a pinch point of 0K between 
% the temperature profiles, or limited by the maximum and minimum temperatures 
% achievable by the fluids.
% 
% The model inputs are:
%       - fluid_h: nature of the hot fluid                        	[-]
%       - P_h_su: inlet pressure of the hot fluid                   [Pa]
%       - in_h_su: inlet temperature or enthalpy of the hot fluid   [K or J/kg]
%       - m_dot_h: mass flow rate of the hot fluid                  [kg/s]
%       - fluid_c: nature of the cold fluid                        	[-]
%       - P_c_su: inlet pressure of the cold fluid                  [Pa]
%       - in_c_su: inlet temperature or enthalpy of the cold fluid  [K or J/kg]
%       - m_dot_c: mass flow rate of the cold fluid                 [kg/s]
%       - param: structure variable containing the model parameters i.e.
%           *param.type_h = type of input for hot fluid ('H' for enthalpy,'T' for temperature)
%           *param.type_c = type of input for cold fluid ('H' for enthalpy,'T' for temperature)
%
% The model outputs are:
%       - Q_dot_max : the maximum amount of power that can be transferred.
%          	
% See the documentation for further details or contact rdickes@ulg.ac.be

%% MODELLING CODE
res_T = 1e-2;

if strcmp(info.type_h,'H') && strcmp(info.type_c,'H') % CASE 1 : HOT FLUID AND COLD FLUID MIGHT EXPERIENCE A PHASE CHANGE
    
    h_h_su = in_h_su;
    h_c_su = in_c_su;
    T_c_su = CoolProp.PropsSI('T', 'P', P_c_su, 'H', h_c_su, fluid_c);
    T_h_su = CoolProp.PropsSI('T', 'P', P_h_su, 'H', h_h_su, fluid_h);
    T_h_min = CoolProp.PropsSI('Tmin', 'P', P_h_su,'Q', 0, fluid_h)+0.5;  
    T_c_max = CoolProp.PropsSI('Tmax', 'P', P_c_su,'Q', 0, fluid_c)-0.5;  
    
    % Computation of Q_dot_hext_max
    if isempty(strfind(fluid_h, 'INCOMP:')) 
        if T_c_su-T_h_min < res_T;
            h_h_ex_extmax = CoolProp.PropsSI('H', 'P', P_h_su,'T', T_h_min, fluid_h);
        elseif abs(T_c_su-CoolProp.PropsSI('T', 'P', P_h_su, 'Q', 0, fluid_h))<res_T;
            h_h_ex_extmax = CoolProp.PropsSI('H', 'P', P_h_su,'Q', 0, fluid_h);
        else
            h_h_ex_extmax =  CoolProp.PropsSI('H', 'P', P_h_su,'T', T_c_su, fluid_h);
        end
    else
        if T_c_su-T_h_min>res_T;
            h_h_ex_extmax =  CoolProp.PropsSI('H', 'P', P_h_su,'T', T_c_su, fluid_h);
        else
            h_h_ex_extmax = CoolProp.PropsSI('H', 'P', P_h_su,'T', T_h_min, fluid_h);
        end       
    end
    Q_dot_hext_max = m_dot_h*(h_h_su-h_h_ex_extmax);
    
    % Computation of Q_dot_cext_max    
    if isempty(strfind(fluid_c, 'INCOMP:'))
        if T_c_max-T_h_su < res_T;
            h_c_ex_extmax = CoolProp.PropsSI('H', 'P', P_c_su,'T', T_c_max, fluid_c);
        elseif abs(T_h_su-CoolProp.PropsSI('T', 'P', P_c_su, 'Q', 1, fluid_c))<res_T;
            h_c_ex_extmax = CoolProp.PropsSI('H', 'P', P_c_su,'Q', 1, fluid_c);
        else
            h_c_ex_extmax =  CoolProp.PropsSI('H', 'P', P_c_su,'T', T_h_su, fluid_c);
        end
    else
        if T_c_max-T_h_su < res_T;
            h_c_ex_extmax = CoolProp.PropsSI('H', 'P', P_c_su,'T', T_c_max, fluid_c);
        else
            h_c_ex_extmax =  CoolProp.PropsSI('H', 'P', P_c_su,'T', T_h_su, fluid_c);
        end
    end
    Q_dot_cext_max = m_dot_c*(h_c_ex_extmax-h_c_su);
    
    Q_dot_ext_max = min(Q_dot_hext_max,Q_dot_cext_max);
    out_ext_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_ext_max, info);
    if out_ext_max.pinch < -res_T
        [q_c_vec, q_h_vec] = deal(NaN*ones(1,length(out_ext_max.H_c_vec)));
        for k = 1:length(out_ext_max.H_c_vec)           
            if isempty(strfind(fluid_h, 'INCOMP:')) && strcmp(info.type_h,'H')
                q_h_vec(k) = CoolProp.PropsSI('Q','P', P_h_su, 'H', out_ext_max.H_h_vec(k), fluid_h);
            else
                q_h_vec(k) = -1;
            end
            if isempty(strfind(fluid_c, 'INCOMP:')) && strcmp(info.type_c,'H')
                q_c_vec(k) = CoolProp.PropsSI('Q','P', P_c_su, 'H', out_ext_max.H_c_vec(k), fluid_c);
            else
                q_c_vec(k) = -1;
            end
        end
        i_ext_max_rev = find(out_ext_max.T_h_vec-out_ext_max.T_c_vec < -res_T);
        k = 1;
        stop = 0;    
        
        while not(stop) && k <= length(i_ext_max_rev)            
            if q_c_vec(i_ext_max_rev(k)) ~= -1 && q_h_vec(i_ext_max_rev(k)) == -1
                Q_dot_int_max = m_dot_c*(out_ext_max.H_c_vec(i_ext_max_rev(k))-h_c_su) + m_dot_h*(h_h_su-CoolProp.PropsSI('H', 'P', P_h_su, 'T', out_ext_max.T_c_vec(i_ext_max_rev(k)), fluid_h));
            elseif q_c_vec(i_ext_max_rev(k)) == -1 && q_h_vec(i_ext_max_rev(k)) ~= -1
                Q_dot_int_max = m_dot_h*(h_h_su-out_ext_max.H_h_vec(i_ext_max_rev(k))) + m_dot_c*(CoolProp.PropsSI('H', 'P', P_c_su, 'T', out_ext_max.T_h_vec(i_ext_max_rev(k)), fluid_c)-h_c_su);
            else
                f = @(x) pinch0(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, x, info);
                Q_dot_int_max = zeroBrent ( 0, Q_dot_ext_max, 1e-6, 1e-6, f );
            end
            out_int_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_int_max, info);
            if out_int_max.pinch < -res_T
                k = k + 1;
            else
                stop = 1;
            end
        end       
        Q_dot_max = Q_dot_int_max ;
    else
        Q_dot_max = Q_dot_ext_max ;
    end

elseif strcmp(info.type_h,'H') && strcmp(info.type_c,'T')
    h_h_su = in_h_su;
    T_h_su = CoolProp.PropsSI('T', 'P', P_h_su, 'H', h_h_su, fluid_h);
    T_c_su = in_c_su;
    Q_dot_cext_max = m_dot_c*sf_PropsSI_bar('C', T_h_su, T_c_su, P_c_su, fluid_c)*(T_h_su-T_c_su);
    if isempty(strfind(fluid_h, 'INCOMP:'))
        if abs(T_c_su-CoolProp.PropsSI('T', 'P', P_h_su, 'Q', 1, fluid_h))>1e-2;
            if T_c_su>CoolProp.PropsSI('Tmin', 'P', P_h_su,'Q', 0, fluid_h);
                Q_dot_hext_max = m_dot_h*(h_h_su-CoolProp.PropsSI('H', 'P', P_h_su,'T', T_c_su, fluid_h));
            else
                Q_dot_hext_max = m_dot_h*(h_h_su-CoolProp.PropsSI('H', 'P', P_h_su,'T', CoolProp.PropsSI('Tmin', 'P', P_h_su,'Q', 0, fluid_h), fluid_h));
            end
        else
            Q_dot_hext_max = m_dot_h*(h_h_su-CoolProp.PropsSI('H', 'P', P_h_su,'Q', 0, fluid_h));
        end
    else
        if T_c_su>CoolProp.PropsSI('Tmin', 'P', P_h_su,'Q', 0, fluid_h);
            Q_dot_hext_max = m_dot_h*(h_h_su-CoolProp.PropsSI('H', 'P', P_h_su,'T', T_c_su, fluid_h));
        else
            Q_dot_hext_max = m_dot_h*(h_h_su-CoolProp.PropsSI('H', 'P', P_h_su,'T', CoolProp.PropsSI('Tmin', 'P', P_h_su,'Q', 0, fluid_h), fluid_h));
        end
    end
    Q_dot_ext_max = min(Q_dot_hext_max,Q_dot_cext_max);
    out_ext_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_ext_max, info);
    if out_ext_max.pinch < -res_T
        i_ext_max_rev = find(out_ext_max.T_h_vec-out_ext_max.T_c_vec < -res_T);
        k = 1;
        stop = 0;
        while not(stop) && k < length(i_ext_max_rev)+1
            Q_dot_int_max = m_dot_h*(h_h_su-out_ext_max.H_h_vec(i_ext_max_rev(k))) + m_dot_c*sf_PropsSI_bar('C', out_ext_max.T_h_vec(i_ext_max_rev(k)), T_c_su, P_c_su, fluid_c)*(out_ext_max.T_h_vec(i_ext_max_rev(k))-T_c_su);
            out_int_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_int_max, info);
            if out_int_max.pinch < -res_T
                k = k + 1;
            else
                stop = 1;
            end
        end
        Q_dot_max = Q_dot_int_max ;
    else
        Q_dot_max = Q_dot_ext_max ;
    end
    
elseif strcmp(info.type_h,'T') && strcmp(info.type_c,'H')
    
    T_h_su = in_h_su;
    h_c_su = in_c_su;
    T_c_su = CoolProp.PropsSI('T','P',P_c_su,'H',h_c_su,fluid_c);
    T_c_max = CoolProp.PropsSI('Tmax', 'P', P_c_su,'Q', 0, fluid_c)-0.5;

    % Computation of Q_dot_hext_max 
    Q_dot_hext_max = m_dot_h*sf_PropsSI_bar('C', T_h_su, T_c_su, P_h_su, fluid_h)*(T_h_su-T_c_su);
    
    % Computation of Q_dot_cext_max  
    if isempty(strfind(fluid_c, 'INCOMP:'))
        if T_c_max-T_h_su < res_T;
            h_c_ex_extmax = CoolProp.PropsSI('H', 'P', P_c_su,'T', T_c_max, fluid_c);
        elseif abs(T_h_su-CoolProp.PropsSI('T', 'P', P_c_su, 'Q', 1, fluid_c))<res_T;
            h_c_ex_extmax = CoolProp.PropsSI('H', 'P', P_c_su,'Q', 1, fluid_c);
        else
            h_c_ex_extmax =  CoolProp.PropsSI('H', 'P', P_c_su,'T', T_h_su, fluid_c);
        end
    else
        if T_c_max-T_h_su < res_T;
            h_c_ex_extmax = CoolProp.PropsSI('H', 'P', P_c_su,'T', T_c_max, fluid_c);
        else
            h_c_ex_extmax =  CoolProp.PropsSI('H', 'P', P_c_su,'T', T_h_su, fluid_c);
        end
    end
    Q_dot_cext_max = m_dot_c*(h_c_ex_extmax-h_c_su);
        
    Q_dot_ext_max = min(Q_dot_hext_max,Q_dot_cext_max);
    out_ext_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_ext_max, info);
    if out_ext_max.pinch < -res_T
        i_ext_max_rev = find(out_ext_max.T_h_vec-out_ext_max.T_c_vec < -res_T);
        k = 1;
        stop = 0;
        while not(stop) && k <= length(i_ext_max_rev)        
            Q_dot_int_max = m_dot_h*sf_PropsSI_bar('C', T_h_su, out_ext_max.T_c_vec(i_ext_max_rev(k)) - (m_dot_c*(out_ext_max.H_c_vec(i_ext_max_rev(k))-h_c_su))/(m_dot_h*sf_PropsSI_bar('C', out_ext_max.T_c_vec(i_ext_max_rev(k)), out_ext_max.T_c_vec(i_ext_max_rev(k)), P_h_su, fluid_h)), P_h_su, fluid_h)*(T_h_su-out_ext_max.T_c_vec(i_ext_max_rev(k))) + m_dot_c*(out_ext_max.H_c_vec(i_ext_max_rev(k))-h_c_su);
            out_int_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_int_max, info);
            if out_int_max.pinch < -res_T
                k = k + 1;
            else
                stop = 1;
            end
        end
        if k > length(i_ext_max_rev)
            f = @(x) pinch0(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, x, info);
            Q_dot_int_max = zeroBrent ( 0, Q_dot_ext_max, 1e-3, 1e-3, f );
        end
        Q_dot_max = Q_dot_int_max ;
    else
        Q_dot_max = Q_dot_ext_max ;
    end
    
elseif strcmp(info.type_h,'T') && strcmp(info.type_c,'T')
    T_h_su = in_h_su;
    T_c_su = in_c_su;
    Q_dot_hext_max = m_dot_h*sf_PropsSI_bar('C', T_h_su, T_c_su, P_h_su, fluid_h)*(T_h_su-T_c_su);
    Q_dot_cext_max = m_dot_c*sf_PropsSI_bar('C', T_h_su, T_c_su, P_c_su, fluid_c)*(T_h_su-T_c_su);
    Q_dot_max = min(Q_dot_hext_max,Q_dot_cext_max);
end

end

function err = pinch0 (fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
err = out.pinch;
end