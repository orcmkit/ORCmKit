function [out, TS] = DPModel(fluid,P_su,h_su,M_dot,DP)

%% DEMONSTRATION CASE
if nargin == 0
    clear all
    DP.modelType = 'LinPhiDP';
    fluid = 'R245fa';
    P_su = 306700;
    h_su = 4.451925057943054e+05;
    M_dot = 0.097130000000000;

    DP.dp = 2.2197e+04;
    DP.CoeffPol_dp = [0.0357   -1.3186    0.0357]*1e5;
    DP.CoeffLin_dp = [-1.259130260099441e+04    7.453088159887394e+07];
    
    DP.displayResults = 1;
    DP.advancedUser = 0;
end
tstart_dp = tic;

%% INPUTS VERIFICATION
if not(isfield(DP,'displayResults'))
    DP.displayResults = 0;
end

if not(isfield(DP,'advancedUser')) || not(DP.advancedUser)
    if isfield(DP,'modelType')
        
        switch DP.modelType
            case 'CstDP'
                if ~isfield(DP,'dp')
                    warning('Error:InputError',['Missing inputs. Check the inputs required by model ' DP.modelType]);
                    return
                end
            case 'MdotDP'
                if ~isfield(DP,'Coeff_dp')
                    warning('Error:InputError',['Missing inputs. Check the inputs required by model ' DP.modelType]);
                    return
                end
            case 'PhiDP'
                if ~isfield(DP,'Coeff_dp')
                    warning('Error:InputError',['Missing inputs. Check the inputs required by model ' DP.modelType]);
                    return
                end
   
            otherwise
                warning('Error:WrongModelType','Wrong modelType. Please choose among the following: {CstDP|PolFlowDP|LinPhiDP}');
                return
        end
    else
        warning('Error: Missing modelType. Please choose among the following: {CstDP|PolFlowDP|LinPhiDP}');
        return
    end
end

%% PRESSURE DROP MODELING
switch DP.modelType
    case 'CstDP'
        out.dp = DP.dp;
        out.P_ex = P_su-out.dp;
        out.h_ex = h_su;
        out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);
    case 'MdotDP'
        if isfield(DP, 'Coeff_dp')
            out.dp = max(0,min(DP.Coeff_dp(1) + DP.Coeff_dp(2)*M_dot, P_su - 2*CoolProp.PropsSI('pmin','Q',0,'H',1e5,fluid)));
        else
            out.dp = max(0,min(DP.funMdot_dp(M_dot), P_su - 2*CoolProp.PropsSI('pmin','Q',0,'H',1e5,fluid)));
        end
        out.P_ex = P_su-out.dp;
        out.h_ex = h_su;
        out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);    
    case 'PhiDP'
        phi = M_dot^2/CoolProp.PropsSI('D','P',P_su,'H',h_su,fluid);
        if isfield(DP, 'Coeff_dp')
        out.dp = max(0,min(DP.Coeff_dp(1) + DP.Coeff_dp(2)*phi, P_su - 2*CoolProp.PropsSI('pmin','Q',0,'H',1e5,fluid)));
        else
            out.dp = max(0,min(DP.funPhi_dp(phi), P_su - 2*CoolProp.PropsSI('pmin','Q',0,'H',1e5,fluid)));
        end
        out.P_ex = P_su-out.dp;
        out.h_ex = h_su;
        out.T_ex = CoolProp.PropsSI('T','P',out.P_ex,'H',out.h_ex,fluid);      
    otherwise
        disp('Wrong type of model input')
end
out.time = toc(tstart_dp);
out = orderfields(out);

%% TS DIAGRAM and DISPLAY
TS.T(1) = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
TS.T(2) = out.T_ex;
TS.s(1) = CoolProp.PropsSI('S','P',P_su,'H',h_su,fluid);
TS.s(2) = CoolProp.PropsSI('S','P',out.P_ex,'H',out.h_ex,fluid);

if DP.displayResults ==1
    in.fluid = fluid;
    in.M_dot = M_dot;
    in.P_su = P_su;
    in.T_su = CoolProp.PropsSI('T','P',P_su,'H',h_su,fluid);
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

