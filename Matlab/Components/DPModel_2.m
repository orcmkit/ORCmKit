function [out, TS] = DPModel_2(fluid, in_P, in_h, M_dot, T_amb, DP)
%fluid, in_P, in_h, M_dot, T_amb, DP
%% DEMONSTRATION CASE
if nargin == 0
    clear all
    fluid = 'R245fa';
    in_P = 1.0113e+05;
    in_h = 2.0076e+05;
    CoolProp.PropsSI('T','P',in_P,'H',in_h,fluid)
    M_dot = 0.079;
    T_amb = 288.15;
    path = 'C:\Users\RDickes\Google Drive\PhD\MOR study\ORC\Experimental database\Sun2Power';
    DP_folder = [path '\PressureDrops\'];
    load([DP_folder, 'ParametersCalibration_DP.mat'])
    DP = DPLP_PhiDP;
    DP.type_in = 'ex';
    DP.displayResults = 1;
    DP.advancedUser = 0;
end
tstart_dp = tic;

%% INPUTS VERIFICATION
if not(isfield(DP,'displayResults'))
    DP.displayResults = 0;
end


%% PRESSURE DROP MODELING
switch DP.modelType
    case 'CstDP'
        if strcmp(DP.type_in, 'su')
            out = CstDP(in_h, in_P, M_dot, T_amb, fluid, DP);
            out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
            out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
            out.flag = 1;
        elseif strcmp(DP.type_in, 'ex')
                x0 = [in_h in_P];
                ub = 1.5*x0;
                f = @(x) CstDP_res(x, ub, in_h, in_P, M_dot, T_amb, fluid, DP);
                options = optimoptions('fsolve','Display','none');
                x = fsolve(f, x0./ub,options);
                out = CstDP(x(1)*ub(1), x(2)*ub(2), M_dot, T_amb, fluid, DP);
                if norm(f(x)) < 1e-3
                    out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
                    out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
                    out.flag = 1;
                else
                    out.h_ex = in_h;
                    out.P_ex = in_P;
                    out.h_su = in_h;
                    out.P_su = in_P;
                    out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
                    out.T_ex = out.T_su;
                    out.dp = 0;
                    out.Q_dot = 0;
                    out.flag = -1;
                end
        end

        
    case 'MdotDP'
        if strcmp(DP.type_in, 'su')
            out = MdotDP(in_h, in_P, M_dot, T_amb, fluid, DP);
            out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
            out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
            out.flag = 1;
        elseif strcmp(DP.type_in, 'ex')
            x0 = [in_h in_P];
            ub = 1.5*x0;
            f = @(x) MdotDP_res(x, ub, in_h, in_P, M_dot, T_amb, fluid, DP);
            options = optimoptions('fsolve','Display','none');
            x = fsolve(f, x0./ub,options);            
            out = MdotDP(x(1)*ub(1), x(2)*ub(2), M_dot, T_amb, fluid, DP);
            out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
            out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
            if norm(f(x)) < 1e-3
                out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
                out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
                out.flag = 1;
            else
                out.h_ex = in_h;
                out.P_ex = in_P;
                out.h_su = in_h;
                out.P_su = in_P;
                out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
                out.T_ex = out.T_su;
                out.dp = 0;
                out.Q_dot = 0;
                out.flag = -1;
            end
        end
        
   
    case 'PhiDP'
        if strcmp(DP.type_in, 'su')
            out = PhiDP(in_h, in_P, M_dot, T_amb, fluid, DP);
            out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
            out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
            out.flag = 1;
        elseif strcmp(DP.type_in, 'ex')
            x0 = [in_h in_P];
            ub = 1.5*x0;
            f = @(x) PhiDP_res(x, ub, in_h, in_P, M_dot, T_amb, fluid, DP);
            options = optimoptions('fsolve','Display','none');
            x = fsolve(f, x0./ub,options);            
            out = PhiDP(x(1)*ub(1), x(2)*ub(2), M_dot, T_amb, fluid, DP);
            if norm(f(x)) < 1e-3
                out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
                out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
                out.flag = 1;
            else
                out.h_ex = in_h;
                out.P_ex = in_P;
                out.h_su = in_h;
                out.P_su = in_P;
                out.T_su = CoolProp.PropsSI('T','P',out.P_su,'H',out.h_su,fluid);
                out.T_ex = out.T_su;
                out.dp = 0;
                out.Q_dot = 0;
                out.flag = -1;
            end
        end
          
    otherwise
        disp('Wrong type of model input')
end
out.time = toc(tstart_dp);
out = orderfields(out);

%% TS DIAGRAM and DISPLAY
TS.T(1) = out.T_su;
TS.T(2) = out.T_ex;
TS.s(1) = CoolProp.PropsSI('S','P',out.P_su,'H',out.h_su,fluid);
TS.s(2) = CoolProp.PropsSI('S','P',out.P_ex,'H',out.h_ex,fluid);

if DP.displayResults ==1
    in.fluid = fluid;
    in.M_dot = M_dot;
    in.P = in_P;
    in.T = CoolProp.PropsSI('T','P',in_P,'H',in_h,fluid);
    in.modelType= DP.modelType;
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
    disp('Results:')
    disp(out)
end

end

function out = CstDP(h_su, P_su, M_dot, T_amb, fluid, DP)
out.dp = DP.dp;
out.P_su = P_su;
out.h_su = h_su;
out.P_ex = out.P_su-out.dp;
out.Q_dot = DP.AU*(CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_su,fluid)-T_amb);
out.h_ex = h_su-out.Q_dot/M_dot;
end

function res = CstDP_res(x, ub, h_ex, P_ex, M_dot, T_amb, fluid, DP)
x = x.*ub;
out = CstDP(x(1), x(2), M_dot, T_amb, fluid, DP);
res(1) = h_ex-out.h_ex;
res(2) = P_ex-out.P_ex;
end

function out = MdotDP(h_su, P_su, M_dot, T_amb, fluid, DP)
out.P_su = P_su;
out.h_su = h_su;
out.dp = max(0,min(DP.funMdot_dp(M_dot), out.P_su - 2*CoolProp.PropsSI('pmin','Q',0,'H',1e5,fluid)));
out.P_ex = out.P_su-out.dp;
out.Q_dot = DP.AU*(CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_su,fluid)-T_amb);
out.h_ex = h_su-out.Q_dot/M_dot;
end

function res = MdotDP_res(x, ub, h_ex, P_ex, M_dot, T_amb, fluid, DP)
x = x.*ub;
out = MdotDP(x(1), x(2), M_dot, T_amb, fluid, DP);
res(1) = h_ex-out.h_ex;
res(2) = P_ex-out.P_ex;
end

function out = PhiDP(h_su, P_su, M_dot, T_amb, fluid, DP)
out.P_su = P_su;
out.h_su = h_su;
if strcmp(DP.type_phi, 'phi_su')
    phi = M_dot^2/CoolProp.PropsSI('D','P',out.P_su,'H',out.h_su,fluid);
elseif strcmp(DP.type_phi, 'phi_1')
    phi = M_dot^2/CoolProp.PropsSI('D','P',out.P_su,'Q',1,fluid);
elseif strcmp(DP.type_phi, 'phi_0')
    phi = M_dot^2/CoolProp.PropsSI('D','P',out.P_su,'Q',0,fluid);
end
out.dp = max(0,min(DP.funPhi_dp(phi), out.P_su - 2*CoolProp.PropsSI('pmin','Q',0,'H',1e5,fluid)));
out.P_ex = out.P_su-out.dp;
out.Q_dot = DP.AU*(CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_su,fluid)-T_amb);
out.h_ex = h_su-out.Q_dot/M_dot;
end

function res = PhiDP_res(x, ub, h_ex, P_ex, M_dot, T_amb, fluid, DP)
x = x.*ub;
out = PhiDP(x(1), x(2), M_dot, T_amb, fluid, DP);
res(1) = (h_ex-out.h_ex)/h_ex;
res(2) = (P_ex-out.P_ex)/P_ex;
end


function res = in_ex(x, h_ex, P_ex, M_dot, fluid, DP)
P_su = x(1);
h_su = x(2);

if strcmp(DP.type_phi, 'phi_su')
    phi = M_dot^2/CoolProp.PropsSI('D','P',P_su,'H',h_ex,fluid);
elseif strcmp(DP.type_phi, 'phi_1')
    phi = M_dot^2/CoolProp.PropsSI('D','P',P_su,'Q',1,fluid);
elseif strcmp(DP.type_phi, 'phi_0')
    phi = M_dot^2/CoolProp.PropsSI('D','P',P_su,'Q',0,fluid);
end
dp = max(0,min(DP.funPhi_dp(phi), P_su - 2*CoolProp.PropsSI('pmin','Q',0,'H',1e5,fluid)));
res = x-dp;
end

function res = PhiDP_ex(x, h_ex, P_ex, M_dot, fluid, DP)
P_su = P_ex + x;
if strcmp(DP.type_phi, 'phi_su')
    phi = M_dot^2/CoolProp.PropsSI('D','P',P_su,'H',h_ex,fluid);
elseif strcmp(DP.type_phi, 'phi_1')
    phi = M_dot^2/CoolProp.PropsSI('D','P',P_su,'Q',1,fluid);
elseif strcmp(DP.type_phi, 'phi_0')
    phi = M_dot^2/CoolProp.PropsSI('D','P',P_su,'Q',0,fluid);
end
dp = max(0,min(DP.funPhi_dp(phi), P_su - 2*CoolProp.PropsSI('pmin','Q',0,'H',1e5,fluid)));
res = x-dp;
end