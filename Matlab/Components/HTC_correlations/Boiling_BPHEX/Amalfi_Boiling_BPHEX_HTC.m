function [h_boiling, Nu_boiling, flag] = Amalfi_Boiling_BPHEX_HTC(x, rho_l, rho_v, rho, mu_v, mu_l, k_l, sigma, i_fg, G, Dh, theta, Qdot, DTlog, h_conv_h, disp_flag) %---> VERIFIED

% Boiling HTC correlation published by Amalfi in "Flow boiling and
% frictional pressure gradients in plate heat exchangers. Part 2:
% Comparison of literature methods to database and new prediction methods"
% in International Journal of Refrigeration, 2016

% RDickes - 26/04/2017

% QUEL VALEUR DE K PRENDRE ???? CE SERAIT DU LIQUID?!
flag = 1;

g = 9.81; 
Bd = (rho_l-rho_v)*g*Dh^2/sigma; 
rho_star = rho_l/rho_v;
beta_star = theta/(70*pi/180);
We = (G^2*Dh)/(rho*sigma);
Re_v = G*x*Dh/mu_v;
Re_lo = G*Dh/mu_l;
    
AU_tp = Qdot/DTlog;

if Bd < 4 % Micro-scale flow mechanism
    
    Bo_0 = 1;
    f_Bo = @(xx) res_iter_AmalfiMicro_boiling(xx, beta_star, We, rho_star, h_conv_h, k_l, Dh, AU_tp, Qdot, G, i_fg);
    opts = optimoptions('fsolve', 'display', 'none');
    [Bo,err_Bo] = fsolve(f_Bo, Bo_0, opts);
    [~, Nu_boiling, h, ~, ~, q, ~] = iter_AmalfiMicro_boiling(Bo, beta_star, We, rho_star, h_conv_h, k_l, Dh, AU_tp, Qdot, G, i_fg);
    if disp_flag        
        if err_Bo > 1e-5
            flag = - 5;
            display(['Amalfi boiling: Wrong boiling number --> err_Bo = ' num2str(err_Bo) ' !!!'])
        end
    end
    h_boiling = h;
    
else % Macro-scale flow mechanism
    
    Bo_0 = 1;
    f_Bo = @(xx) res_iter_AmalfiMacro_boiling(xx, beta_star, Re_v, Re_lo, Bd, rho_star, k_l, Dh, h_conv_h, AU_tp, Qdot, G, i_fg);
    opts = optimoptions('fsolve', 'display', 'none');
    [Bo,err_Bo] = fsolve(f_Bo, Bo_0, opts);
    [~, Nu_boiling, h, ~, ~, q, ~] = iter_AmalfiMacro_boiling(Bo, beta_star, Re_v, Re_lo, Bd, rho_star, k_l, Dh, h_conv_h, AU_tp, Qdot, G, i_fg);
    if disp_flag
        if err_Bo > 1e-5
            flag = - 5;
            display(['Amalfi boiling: Wrong boiling number --> err_Bo = ' num2str(err_Bo) ' !!!'])
        end
    end
    h_boiling = h;
   
    
end

end

function res_Bo = res_iter_AmalfiMicro_boiling(Bo_g, beta_star, We, rho_star, h_conv_h, k_l, Dh, AU_tp, Qdot, G, i_fg)
[res_Bo, ~, ~, ~, ~, ~, ~] = iter_AmalfiMicro_boiling(Bo_g, beta_star, We, rho_star, h_conv_h, k_l, Dh, AU_tp, Qdot, G, i_fg);
end

function [res_Bo, Nu, h, U, A_tp, q, Bo] = iter_AmalfiMicro_boiling(Bo_g, beta_star, We, rho_star, h_conv_h, k_l, Dh, AU_tp, Qdot, G, i_fg)
Nu = 982*beta_star^1.101*We^0.315*Bo_g^0.32*rho_star^-0.224;
h = Nu*k_l/Dh;
U = (1/h +  1/h_conv_h)^-1;
A_tp = AU_tp/U;
q = Qdot/A_tp;
Bo = q/(G*i_fg);
res_Bo = abs(Bo-Bo_g)/Bo_g;
end

function res_Bo = res_iter_AmalfiMacro_boiling(Bo_g, beta_star, Re_v, Re_lo, Bd, rho_star, k_l, Dh, h_conv_h, AU_tp, Qdot, G, i_fg)
[res_Bo, ~, ~, ~, ~, ~, ~] = iter_AmalfiMacro_boiling(Bo_g, beta_star, Re_v, Re_lo, Bd, rho_star, k_l, Dh, h_conv_h, AU_tp, Qdot, G, i_fg);
end

function [res_Bo, Nu, h, U, A_tp, q, Bo] = iter_AmalfiMacro_boiling(Bo_g, beta_star, Re_v, Re_lo, Bd, rho_star, k_l, Dh, h_conv_h, AU_tp, Qdot, G, i_fg)
Nu = 18.495*beta_star^0.248*Re_v^0.135*Re_lo^0.351*Bd^0.235*Bo_g^0.198*rho_star^-0.223;
h = Nu*k_l/Dh;
U = (1/h +  1/h_conv_h)^-1;
A_tp = AU_tp/U;
q = Qdot/A_tp;
Bo = q/(G*i_fg);
res_Bo = abs(Bo-Bo_g)/Bo_g;
end
