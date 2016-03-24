function [out,TS] = HexModel(fluid_h, P_h_su, in_h_su, m_dot_h, fluid_c, P_c_su, in_c_su, m_dot_c, HEX)
%% DEMONSTRATION CASE

if nargin == 0
    clear all
    close all
    clc
    HEX.modelType = 'hConvCst';
    fluid_h = 'R245fa';
    m_dot_h = 0.05;
    P_h_su =  3e5;
    in_h_su =  CoolProp.PropsSI('H','P',P_h_su,'T',100 +273.15,fluid_h);
    HEX.type_h = 'H';
    fluid_c = 'R245fa';
    m_dot_c = 0.05; %0.512129779272778;
    P_c_su = 4e5; %4.251885856949260e+05;
    in_c_su = CoolProp.PropsSI('H','P',P_c_su,'T',20 +273.15,fluid_c);
    HEX.type_c = 'H';

    HEX.A_tot = 5.45000;
    HEX.pinch = 6.406419831493446; 
    HEX.epsilon_th = 1;
    HEX.hConv_h_liq = 500;
    HEX.hConv_h_tp = 10000;
    HEX.hConv_h_vap = 200;
    HEX.hConv_c_liq = 500;
    HEX.hConv_c_tp = 10000;
    HEX.hConv_c_vap = 200;
    HEX.hConv_c_liq_n = 2.379272937774658e+02;
    HEX.hConv_c_tp_n = 2.379272937774658e+02;
    HEX.hConv_c_vap_n = 2.379272937774658e+02;
    HEX.hConv_h_liq_n = 1.125000000000000e+02;
    HEX.hConv_h_tp_n = 1.125000000000000e+02;
    HEX.hConv_h_vap_n = 1.125000000000000e+02;
    HEX.n = 0.8;
    HEX.m = 0.7;
    HEX.C_c_liq = 0.308120727539063;
    HEX.C_c_tp = 0.308120727539063;
    HEX.C_c_vap =0.308120727539063;
    HEX.C_h_liq = 1.833880424499512;
    HEX.C_h_tp = 2.999996185302734;
    HEX.C_h_vap = 1.833880424499512;
    HEX.CS = 0.020845188;
    HEX.Dh = 0.356328;
    HEX.m_dot_h_n = 0.618652995944947;
    HEX.m_dot_c_n = 0.618652995944947;
    HEX.m_dot_n = 0.149;
    HEX.CoeffPolEff = [0.913872425354551 1.601103316261962 -1.161513908064705 0.244892183446199 -2.586460052802316 1.912516752474164];
    HEX.displayResults = 0;
    HEX.displayTS = 1;
    HEX.advancedUser = 1;
end

tstart = tic;

%% INPUTS VERIFICATION

if not(isfield(HEX,'displayResults'))
    HEX.displayResults = 0;
    HEX.displayTS = 0;
end

if not(isfield(HEX,'advancedUser')) || not(HEX.advancedUser)
    if isfield(HEX,'modelType')
        
        switch HEX.modelType
            case 'PinchCst'
                if sum(isfield(HEX,{'pinch', 'type_h', 'type_c'})) ~= 3
                    warning('Error:InputError',['Missing inputs. Check the inputs required by model ' HEX.modelType]);
                    return
                end
            case 'CstEff'
                if  sum(isfield(HEX,{'epsilon_th', 'type_h', 'type_c'})) ~= 3
                    warning('Error:InputError',['Missing inputs. Check the inputs required by model ' HEX.modelType]);
                    return
                end
            case 'PolEff'
                if sum(isfield(HEX,{'m_dot_h_n', 'm_dot_c_n', 'CoeffPolEff', 'type_h', 'type_c'})) ~= 5
                    warning('Error:InputError',['Missing inputs. Check the inputs required by model ' HEX.modelType]);
                    return
                end
            case 'hConvCst'
                if sum(isfield(HEX,{'hConv_h_liq','hConv_h_tp','hConv_h_vap', 'hConv_c_liq','hConv_c_tp','hConv_c_vap','A_tot', 'type_h', 'type_c'})) ~= 9
                    warning('Error:InputError',['Missing inputs. Check the inputs required by model ' HEX.modelType]);
                    return
                end
            case 'hConvVar'
                if sum(isfield(HEX,{'hConv_h_liq_n','hConv_h_tp_n','hConv_h_vap_n', 'hConv_c_liq_n','hConv_c_tp_n','hConv_c_vap_n', 'n', 'm_dot_h_n', 'm_dot_c_n', 'A_tot','type_h', 'type_c'})) ~= 12
                    warning('Error:InputError',['Missing inputs. Check the inputs required by model ' HEX.modelType]);
                    return
                end
            case 'hConvCor'
                if sum(isfield(HEX,{'C_h_liq','C_h_tp','C_h_vap', 'C_c_liq','C_c_tp','C_c_vap', 'n', 'm', 'CS', 'Dh', 'A_tot','type_h', 'type_c'})) ~= 13
                    warning('Error:InputError',['Missing inputs. Check the inputs required by model ' HEX.modelType]);
                    return
                end
            otherwise
                warning('Error:WrongModelType','Wrong modelType. Please choose among the following: {PinchCst|hConvCst|hConvVar|hConvCor}');
                return
        end
    else
        warning('Error: Missing modelType. Please choose among the following: {PinchCst|hConvCst|hConvVar|hConvCor}');
        return
    end
end

%% HEX MODELING
if strcmp(HEX.type_h,'H')
    T_h_su = CoolProp.PropsSI('T','P',P_h_su,'H',in_h_su, fluid_h);
elseif strcmp(HEX.type_h,'T')
    T_h_su = in_h_su;
end
if strcmp(HEX.type_c,'H')
    T_c_su = CoolProp.PropsSI('T','P',P_c_su,'H',in_c_su, fluid_c);
elseif strcmp(HEX.type_c,'T')
    T_c_su = in_c_su;
end

if (T_h_su-T_c_su)>1e-2  && m_dot_h  > 0 && m_dot_c > 0
    switch HEX.modelType
        
        case 'PinchCst'           
            if T_h_su-T_c_su > HEX.pinch
                Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, HEX);
                lb = 0;
                ub = Q_dot_max;
                f=@(Q_dot) HEX_PinchCst_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, HEX.pinch, Q_dot, HEX);
                Q_dot_eff = zeroBrent ( lb, ub, 1e-6, 1e-6, f );
                out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, HEX);
                [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
                out.Q_dot_tot = Q_dot_eff;
                out.h_h_ex = out.H_h_vec(1);
                out.h_c_ex = out.H_c_vec(end);
                out.T_h_ex = out.T_h_vec(1);
                out.T_c_ex = out.T_c_vec(end);
                out.resPinch = abs(1-out.pinch/HEX.pinch);
                if out.resPinch <1e-4
                    out.flag = 1;
                else
                    out.flag = -1;
                end
                if strcmp(HEX.type_h,'H')
                    for i = 1: length(out.H_h_vec)
                        out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                    end
                end
                if strcmp(HEX.type_c,'H')
                    for i = 1: length(out.H_c_vec)
                        out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                    end
                end
            else
                Q_dot_eff = 0;
                out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, HEX);
                out.h_h_ex = out.H_h_vec(1);
                out.h_c_ex = out.H_c_vec(end);
                out.T_h_ex = out.T_h_vec(1);
                out.T_c_ex = out.T_c_vec(end);
                out.Q_dot_tot = Q_dot_eff;
                [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
                if strcmp(HEX.type_h,'H')
                    for i = 1: length(out.H_h_vec)
                        out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                    end
                end
                if strcmp(HEX.type_c,'H')
                    for i = 1: length(out.H_c_vec)
                        out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                    end
                end
                out.flag = 2;
            end
            
        case 'CstEff'
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, HEX);
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, HEX);
            Q_dot_eff = HEX.epsilon_th*Q_dot_max; 
            out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, HEX);
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            out.Q_dot_tot = Q_dot_eff;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);
            if abs(out_max.pinch) < 1e-2
                out.flag = 1;
            else
                out.flag = -2;
            end
            if strcmp(HEX.type_h,'H')
                for i = 1: length(out.H_h_vec)
                    out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                end
            end
            if strcmp(HEX.type_c,'H')
                for i = 1: length(out.H_c_vec)
                    out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                end
            end
            
        case 'PolEff'
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, HEX);
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, HEX);
            Q_dot_eff = Q_dot_max*max(1e-5,min(1, HEX.CoeffPolEff(1) + HEX.CoeffPolEff(2)*(m_dot_h/HEX.m_dot_h_n) + HEX.CoeffPolEff(3)*(m_dot_c/HEX.m_dot_c_n) + HEX.CoeffPolEff(4)*(m_dot_h/HEX.m_dot_h_n)^2 + HEX.CoeffPolEff(5)*(m_dot_h/HEX.m_dot_h_n)*(m_dot_c/HEX.m_dot_c_n) + HEX.CoeffPolEff(6)*(m_dot_c/HEX.m_dot_c_n)^2));
            out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, HEX);
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            out.Q_dot_tot = Q_dot_eff;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);
            if abs(out_max.pinch) < 1e-2
                out.flag = 1;
            else
                out.flag = -2;
            end
            if strcmp(HEX.type_h,'H')
                for i = 1: length(out.H_h_vec)
                    out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                end
            end
            if strcmp(HEX.type_c,'H')
                for i = 1: length(out.H_c_vec)
                    out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                end
            end
            
            
        case 'hConvCst'
            if isfield(HEX, 'A_tot')
                HEX.A_h_tot = HEX.A_tot;
                HEX.A_c_tot = HEX.A_tot;
            end
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, HEX);
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, HEX);
            lb = 1;
            ub = Q_dot_max;
            f = @(Q_dot) HEX_hConvCst_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su,  Q_dot, HEX);
            if f(ub) > 0
                Q_dot_eff = ub;
            else
                Q_dot_eff = zeroBrent ( lb, ub, 1e-6, 1e-6, f);
            end
            out = HEX_hConvCst(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, HEX);
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            out.Q_dot_tot = Q_dot_eff;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);
            if strcmp(HEX.type_h,'H')
                for i = 1: length(out.H_h_vec)
                    out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                end
            end
            if strcmp(HEX.type_c,'H')
                for i = 1: length(out.H_c_vec)
                    out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                end
            end           
            
            if out.resA <1e-4
                out.flag = 1;
            else
                if abs(out_max.pinch) < 1e-2
                    if Q_dot_eff == Q_dot_max
                        out.flag = 2; 
                    else
                        out.flag = -1; 
                    end                 
                else
                    out.flag = -2;
                end 
            end
            
        case 'hConvVar'
            if isfield(HEX, 'A_tot')
                HEX.A_h_tot = HEX.A_tot;
                HEX.A_c_tot = HEX.A_tot;
            end
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, HEX);
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, HEX);
            lb = 0;
            ub = Q_dot_max;
            f = @(Q_dot) HEX_hConvVar_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su,  Q_dot, HEX);
            if f(ub) > 0
                Q_dot_eff = ub;
            else
                Q_dot_eff = zeroBrent ( lb, ub, 1e-6, 1e-6, f );
           end
           out = HEX_hConvVar(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, HEX);
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            out.Q_dot_tot = Q_dot_eff;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);
            if strcmp(HEX.type_h,'H')
                for i = 1: length(out.H_h_vec)
                    out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                end
            end
            if strcmp(HEX.type_c,'H')
                for i = 1: length(out.H_c_vec)
                    out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                end
            end
            if out.resA <1e-4
                out.flag = 1;
            else
                if abs(out_max.pinch) < 1e-2
                    if Q_dot_eff == Q_dot_max
                        out.flag = 2; 
                    else
                        out.flag = -1; 
                    end
                    
                else
                    out.flag = -2;
                end 
            end
            
        case 'hConvCor'
            if isfield(HEX, 'A_tot')
                HEX.A_h_tot = HEX.A_tot;
                HEX.A_c_tot = HEX.A_tot;
            end
            Q_dot_max = HEX_Qdotmax(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, HEX);
            out_max = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_max, HEX);
            lb = 1;
            ub = Q_dot_max;
            f = @(Q_dot) HEX_hConvCor_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su,  Q_dot, HEX);
            if f(ub) > 0
                Q_dot_eff = ub;
            else
                Q_dot_eff = zeroBrent ( lb, ub, 1e-6, 1e-6, f );
            end
            out = HEX_hConvCor(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, HEX);
            [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
            out.Q_dot_tot = Q_dot_eff;
            out.h_h_ex = out.H_h_vec(1);
            out.h_c_ex = out.H_c_vec(end);
            out.T_h_ex = out.T_h_vec(1);
            out.T_c_ex = out.T_c_vec(end);
            if strcmp(HEX.type_h,'H')
                for i = 1: length(out.H_h_vec)
                    out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
                end
            end
            if strcmp(HEX.type_c,'H')
                for i = 1: length(out.H_c_vec)
                    out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
                end
            end
            if out.resA <1e-4
                out.flag = 1;
            else
                if abs(out_max.pinch) < 1e-2
                    if Q_dot_eff == Q_dot_max
                        out.flag = 2; 
                    else
                        out.flag = -1; 
                    end                    
                else
                    out.flag = -2;
                end 
            end
            
        otherwise
            disp('Wrong type of model input')
    end
    
else
    Q_dot_eff = 0;
    out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot_eff, HEX);
    out.h_h_ex = out.H_h_vec(1);
    out.h_c_ex = out.H_c_vec(end);
    out.T_h_ex = out.T_h_vec(1);
    out.T_c_ex = out.T_c_vec(end);
    out.Q_dot_tot = Q_dot_eff;
    [out.s_h_vec, out.s_c_vec] = deal(NaN*ones(1, length(out.H_h_vec)));
    if strcmp(HEX.type_h,'H')
        for i = 1: length(out.H_h_vec)
            out.s_h_vec(i) = CoolProp.PropsSI('S','P',P_h_su,'H',out.H_h_vec(i),fluid_h);
        end
    end
    if strcmp(HEX.type_c,'H')
        for i = 1: length(out.H_c_vec)
            out.s_c_vec(i) = CoolProp.PropsSI('S','P',P_c_su,'H',out.H_c_vec(i),fluid_c);
        end
    end
    
    if (T_h_su-T_c_su)<1e-2  && (T_h_su-T_c_su >0)  && m_dot_h  > 0 && m_dot_c > 0
        out.flag = 3;
    else
        out.flag = -3;
    end
end
out.time = toc(tstart);

%% TS DIAGRAM and DISPLAY
TS.T_h = out.T_h_vec;
TS.T_c = out.T_c_vec;
TS.s_h = out.s_h_vec;
TS.s_c = out.s_c_vec;
TS.x = out.x_vec;
if HEX.displayTS == 1
    figure
    hold on
    plot(out.x_vec, TS.T_c-273.15,'s-' ,'linewidth',2)
    plot(out.x_vec, TS.T_h-273.15,'o-' ,'linewidth',2)
    grid on
    xlabel('Heat power fraction [-]','fontsize',14,'fontweight','bold')
    ylabel('Temperature [°C]','fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14,'fontweight','bold')
end

if HEX.displayResults ==1
    in.fluid_h = fluid_h;
    in.m_dot_h = m_dot_h;
    in.in_h_su = in_h_su;
    in.type_h = HEX.type_h;
    in.P_h_su = P_h_su;
    in.fluid_c = fluid_c;
    in.m_dot_c = m_dot_c;
    in.in_c_su = in_c_su;
    in.type_c = HEX.type_c;
    in.P_c_su = P_c_su;
    in.modelType= HEX.modelType;
    
    if nargin ==0
        fprintf ( 1, '\n' );
        disp('-------------------------------------------------------')
        disp('--------------------   Demo Code   --------------------')
        disp('-------------------------------------------------------')
        fprintf ( 1, '\n' );
    end
    disp('Working conditions:')
    fprintf ( 1, '\n' );
    disp(in)
    disp('Results')
    disp(out)
    
end

end


function res = HEX_PinchCst_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, pinch, Q_dot, info)
out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
res = pinch - out.pinch;
end

function res = HEX_hConvCst_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_hConvCst(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
res = out.resA;
end

function out = HEX_hConvCst(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
[out.hConv_h, out.hConv_c, out.DTlog, out.eff_h, out.eff_c, out.A_h, out.A_c, out.U] = deal(NaN*ones(1,length(out.T_h_vec)-1));
for j = 1:length(out.T_h_vec)-1
    % Hot side heat transfer coefficient
    if strcmp(info.type_h, 'H')
        if isempty(strfind(fluid_h, 'INCOMP:'))
            if (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) < CoolProp.PropsSI('H','P',P_h_su,'Q',0,fluid_h)
                out.hConv_h(j) = info.hConv_h_liq;
            elseif (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) > CoolProp.PropsSI('H','P',P_h_su,'Q',1,fluid_h)
                out.hConv_h(j) = info.hConv_h_vap;
            else
                out.hConv_h(j) = info.hConv_h_tp;
            end
        else
            out.hConv_h(j) = info.hConv_h_liq;
        end
    elseif strcmp(info.type_h, 'T')
        out.hConv_h(j) = info.hConv_h_liq;
    end
    % Cold side heat transfer coefficient
    if strcmp(info.type_c, 'H')
        if isempty(strfind(fluid_c, 'INCOMP:'))
            if (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) < CoolProp.PropsSI('H','P',P_c_su,'Q',0,fluid_c)
                out.hConv_c(j) = info.hConv_c_liq;
            elseif (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) > CoolProp.PropsSI('H','P',P_c_su,'Q',1,fluid_c)
                out.hConv_c(j) = info.hConv_c_vap;
            else
                out.hConv_c(j) = info.hConv_c_tp;
            end
            out.hConv_c(j) = info.hConv_c_liq;
        else
            out.hConv_c(j) = info.hConv_c_liq;
        end
    elseif strcmp(info.type_c, 'T')
        out.hConv_c(j) = info.hConv_c_liq;
    end
    out.DTlog(j) = deltaT_log(out.T_h_vec(j+1), out.T_h_vec(j),out.T_c_vec(j), out.T_c_vec(j+1));
    
    % Hot side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_h'))
        out.eff_h(j) = 1;
    elseif strcmp(info.fin_h, 'none')
        out.eff_h(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_h(j),info.fin_h.k, info.fin_h.th, info.fin_h.r, info.fin_h.B, info.fin_h.H);
        out.eff_h(j) = 1-info.fin_h.omega_f*(1-eta_eff);
    end
    
    % Cold side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_c'))
        out.eff_c(j) = 1;
    elseif strcmp(info.fin_c, 'none')
        out.eff_c(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_c(j), info.fin_c.k, info.fin_c.th, info.fin_c.r, info.fin_c.B, info.fin_c.H);
        out.eff_c(j) = 1-info.fin_c.omega_f*(1-eta_eff);
    end
    
    % Global heat transfer coefficient and zone surface area
    out.U(j) = (1/out.hConv_h(j)/out.eff_h(j) + 1/out.hConv_c(j)/out.eff_c(j)/(info.A_c_tot/info.A_h_tot))^-1;
    out.A_h(j) = out.Qdot_vec(j)/out.DTlog(j)/out.U(j);
    out.A_c(j) = out.A_h(j)*info.A_c_tot/info.A_h_tot;
end
out.A_h_tot = sum(out.A_h);
out.resA = 1 - out.A_h_tot/info.A_h_tot;
end

function res = HEX_hConvVar_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_hConvVar(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
res = out.resA;
end

function out = HEX_hConvVar(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
[out.hConv_h, out.hConv_c, out.DTlog, out.eff_h, out.eff_c, out.A_h, out.A_c, out.U] = deal(NaN*ones(1,length(out.H_h_vec)-1));
for j = 1:length(out.T_h_vec)-1
    % Hot side heat transfer coefficient
    if strcmp(info.type_h, 'H')
        if isempty(strfind(fluid_h, 'INCOMP:'))
            if (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) < CoolProp.PropsSI('H','P',P_h_su,'Q',0,fluid_h)
                out.hConv_h(j) = info.hConv_h_liq_n*(m_dot_h/info.m_dot_h_n)^info.n;
            elseif (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) > CoolProp.PropsSI('H','P',P_h_su,'Q',1,fluid_h)
                out.hConv_h(j) = info.hConv_h_vap_n*(m_dot_h/info.m_dot_h_n)^info.n;
            else
                out.hConv_h(j) = info.hConv_h_tp_n*(m_dot_h/info.m_dot_h_n)^info.n;
            end
        else
            out.hConv_h(j) = info.hConv_h_liq_n*(m_dot_h/info.m_dot_h_n)^info.n;
        end
    elseif strcmp(info.type_h, 'T')
        out.hConv_h(j) = info.hConv_h_liq_n*(m_dot_h/info.m_dot_h_n)^info.n;
    end
    
    % Cold side heat transfer coefficient
    if strcmp(info.type_c, 'H')
        if isempty(strfind(fluid_c, 'INCOMP:'))
            if (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) < CoolProp.PropsSI('H','P',P_c_su,'Q',0,fluid_c)
                out.hConv_c(j) = info.hConv_c_liq_n*(m_dot_c/info.m_dot_c_n)^info.n;
            elseif (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) > CoolProp.PropsSI('H','P',P_c_su,'Q',1,fluid_c)
                out.hConv_c(j) = info.hConv_c_vap_n*(m_dot_c/info.m_dot_c_n)^info.n;
            else
                out.hConv_c(j) = info.hConv_c_tp_n*(m_dot_c/info.m_dot_c_n)^info.n;
            end
        else
            out.hConv_c(j) = info.hConv_c_liq_n*(m_dot_c/info.m_dot_c_n)^info.n;
        end
    elseif strcmp(info.type_c, 'T')
        out.hConv_c(j) = info.hConv_c_liq_n*(m_dot_c/info.m_dot_c_n)^info.n;
    end
    out.DTlog(j) = deltaT_log(out.T_h_vec(j+1), out.T_h_vec(j),out.T_c_vec(j), out.T_c_vec(j+1));
    
    % Hot side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_h'))
        out.eff_h(j) = 1;
    elseif strcmp(info.fin_h, 'none')
        out.eff_h(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_h(j), info.fin_h.k, info.fin_h.th, info.fin_h.r, info.fin_h.B, info.fin_h.H);
        out.eff_h(j) = 1-info.fin_h.omega_f*(1-eta_eff);
    end
    
    % Cold side heat transfer efficiency (in case of fins)
    if not(isfield(info, 'fin_c'))
        out.eff_c(j) = 1;
    elseif strcmp(info.fin_c, 'none')
        out.eff_c(j) = 1;
    else
        eta_eff = FinSchmidt(out.hConv_c(j), info.fin_c.k, info.fin_c.th, info.fin_c.r, info.fin_c.B, info.fin_c.H);
        out.eff_c(j) = 1-info.fin_c.omega_f*(1-eta_eff);
    end
    
    % Global heat transfer coefficient and zone surface area
    out.U(j) = (1/out.hConv_h(j)/out.eff_h(j) + 1/out.hConv_c(j)/out.eff_c(j)/(info.A_c_tot/info.A_h_tot))^-1;
    out.A_h(j) = out.Qdot_vec(j)/out.DTlog(j)/out.U(j);
    out.A_c(j) = out.A_h(j)*info.A_c_tot/info.A_h_tot;
end
out.A_h_tot = sum(out.A_h);
out.resA = 1 - out.A_h_tot/info.A_h_tot;
end

function res = HEX_hConvCor_res(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_hConvCor(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
res = out.resA;
end

function out = HEX_hConvCor(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info)
out = HEX_profile(fluid_h, m_dot_h, P_h_su, in_h_su, fluid_c, m_dot_c, P_c_su, in_c_su, Q_dot, info);
[out.A, out.hConv_h, out.hConv_c, out.DTlog] = deal(NaN*ones(1,length(out.H_h_vec)-1));
for j = 1:length(out.T_h_vec)-1
    if strcmp(info.type_h, 'H')
        if isempty(strfind(fluid_h, 'INCOMP:'))
            if (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) < CoolProp.PropsSI('H', 'P', P_h_su, 'Q', 0, fluid_h)
                mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);              
                Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                k = CoolProp.PropsSI('L', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                out.hConv_h(j) = hConvCor(mu, Pr, k, m_dot_h, info.CS, info.Dh, info.C_h_liq, info.m, info.n);
            elseif (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)) > CoolProp.PropsSI('H','P',P_h_su,'Q',1,fluid_h)
                mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                k = CoolProp.PropsSI('L', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
                out.hConv_h(j) = hConvCor(mu, Pr, k, m_dot_h, info.CS, info.Dh, info.C_h_vap, info.m, info.n);
            else
                mu = CoolProp.PropsSI('V', 'Q', 0, 'P', P_h_su, fluid_h);
                Pr = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_h_su, fluid_h);
                k = CoolProp.PropsSI('L', 'Q', 0, 'P', P_h_su, fluid_h);
                hConv_tp = hConvCor(mu, Pr, k, m_dot_h, info.CS, info.Dh, info.C_h_liq, info.m, info.n);
                out.hConv_h(j) = info.C_h_tp*hConv_tp;
            end
        else
            mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
            Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
            k = CoolProp.PropsSI('L', 'H', (0.5*out.H_h_vec(j)+0.5*out.H_h_vec(j+1)), 'P', P_h_su, fluid_h);
            out.hConv_h(j) = hConvCor(mu, Pr, k, m_dot_h, info.CS, info.Dh, info.C_h_liq, info.m, info.n);
        end
    elseif strcmp(info.type_h, 'T')
        cp = sf_PropsSI_bar('C', out.T_h_vec(j), out.T_h_vec(j+1), P_h_su, fluid_h);
        k = sf_PropsSI_bar('L', out.T_h_vec(j), out.T_h_vec(j+1), P_h_su, fluid_h);
        mu = sf_PropsSI_bar('V', out.T_h_vec(j), out.T_h_vec(j+1), P_h_su, fluid_h);
        Pr = cp*mu/k;
        out.hConv_h(j) = hConvCor(mu, Pr, k, m_dot_h, info.CS, info.Dh, info.C_h_liq, info.m, info.n);
    end
    if strcmp(info.type_c, 'H')
        if isempty(strfind(fluid_c, 'INCOMP:'))
            if (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) < CoolProp.PropsSI('H', 'P', P_c_su, 'Q', 0, fluid_c)
                mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                k = CoolProp.PropsSI('L', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                out.hConv_c(j) = hConvCor(mu, Pr, k, m_dot_c, info.CS, info.Dh, info.C_c_liq, info.m, info.n);
            elseif (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)) > CoolProp.PropsSI('H','P',P_c_su,'Q',1,fluid_c)
                mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                k = CoolProp.PropsSI('L', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
                out.hConv_c(j) = hConvCor(mu, Pr, k, m_dot_c, info.CS, info.Dh, info.C_c_vap, info.m, info.n);
            else
                mu = CoolProp.PropsSI('V', 'Q', 0, 'P', P_c_su, fluid_c);
                Pr = CoolProp.PropsSI('Prandtl', 'Q', 0, 'P', P_c_su, fluid_c);
                k = CoolProp.PropsSI('L', 'Q', 0, 'P', P_c_su, fluid_c);
                hConv_tp = hConvCor(mu, Pr, k, m_dot_c, info.CS, info.Dh, info.C_c_liq, info.m, info.n);
                out.hConv_c(j) = info.C_c_tp*hConv_tp;
            end
        else
            mu = CoolProp.PropsSI('V', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
            Pr = CoolProp.PropsSI('Prandtl', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
            k = CoolProp.PropsSI('L', 'H', (0.5*out.H_c_vec(j)+0.5*out.H_c_vec(j+1)), 'P', P_c_su, fluid_c);
            out.hConv_c(j) = hConvCor(mu, Pr, k, m_dot_c, info.CS, info.Dh, info.C_c_liq, info.m, info.n);
        end
    elseif strcmp(info.type_c, 'T')
        cp = sf_PropsSI_bar('C', T_c_vec(j), T_c_vec(j+1), P_c_su, fluid_c);
        k = sf_PropsSI_bar('L', T_c_vec(j), T_c_vec(j+1), P_c_su, fluid_c);
        mu = sf_PropsSI_bar('V', T_c_vec(j), T_c_vec(j+1), P_c_su, fluid_c);
        Pr = cp*mu/k;
        out.hConv_c(j) = hConvCor(mu, Pr, k, m_dot_c, info.CS, info.Dh, info.C_c_liq, info.m, info.n);
    end
    out.DTlog(j) = deltaT_log(out.T_h_vec(j+1), out.T_h_vec(j),out.T_c_vec(j), out.T_c_vec(j+1));
    out.A(j) = out.Qdot_vec(j)/out.DTlog(j)/(1/(1/out.hConv_h(j) + 1/out.hConv_c(j)));
end
out.A_tot = sum(out.A);
out.resA = 1 - out.A_tot/info.A_tot;
end

function [DT_log, DTh, DTc ]= deltaT_log(Th_su, Th_ex, Tc_su, Tc_ex)
DTh = max(Th_su-Tc_ex,1e-2);
DTc = max(Th_ex-Tc_su,1e-2);
if DTh ~= DTc;
    DT_log = (DTh-DTc)/log(DTh/DTc);
else
    DT_log = DTh;
end
end

function eta_fin = FinSchmidt(hConv, k, th, r, B, H)
m = sqrt(2*hConv/k/th);
phi_f = B/r;
beta_f = H/B;
R_e = r*1.27*phi_f*(beta_f-0.3)^0.5;
phi = (R_e/r - 1)*(1+0.35*log(R_e/r));
eta_fin = tanh(m*R_e*phi)/(m*R_e*phi);
end

function h = hConvCor(mu, Pr, k, m_dot, CS, Dh, C, m, n)
G = m_dot/CS;
Re = G*Dh/mu;
Nu = C*Re^m*Pr^n;
h = Nu*k/Dh;
end