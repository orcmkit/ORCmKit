
%close all
clear all
clc
fluid = 'R134a';


P_crit = CoolProp.PropsSI('Pcrit','P',1e5,'Q',0,fluid);   
T_crit = CoolProp.PropsSI('Tcrit','P',1e5,'Q',0,fluid);   

T_min = 20 + 273.15;    
T_max = 350 + 273.15;   
T_vec = [linspace(T_min, T_crit-25, 200) linspace(T_crit-25, T_crit+25, 1000) linspace(T_crit+25, T_max, 200)];
P_vec = linspace(P_crit, 3*P_crit, 10);

% Saturated properties
P_sat_vec = linspace(2e5, P_crit, 100);
for i = 1:length(P_sat_vec)
    P_sat(i) = P_sat_vec(i);
    T_sat_l(i) = CoolProp.PropsSI('T','P',P_sat(i),'Q',0,fluid);
    T_sat_v(i) = CoolProp.PropsSI('T','P',P_sat(i),'Q',1,fluid);
    rho_sat_l(i) = CoolProp.PropsSI('D','P',P_sat(i),'Q',0,fluid);
    rho_sat_v(i) = CoolProp.PropsSI('D','P',P_sat(i),'Q',1,fluid);
    cp_sat_l(i) = CoolProp.PropsSI('C','P',P_sat(i),'Q',0,fluid);
    cp_sat_v(i) = CoolProp.PropsSI('C','P',P_sat(i),'Q',1,fluid);
    l_sat_l(i) = CoolProp.PropsSI('L','P',P_sat(i),'Q',0,fluid);
    l_sat_v(i) = CoolProp.PropsSI('L','P',P_sat(i),'Q',1,fluid);
    mu_sat_l(i) = CoolProp.PropsSI('V','P',P_sat(i),'Q',0,fluid);
    mu_sat_v(i) = CoolProp.PropsSI('V','P',P_sat(i),'Q',1,fluid);        
end
T_sat = [T_sat_l flip(T_sat_v)];
rho_sat = [rho_sat_l flip(rho_sat_v)];
cp_sat = [cp_sat_l flip(cp_sat_v)];
l_sat = [l_sat_l flip(l_sat_v)];
mu_sat = [mu_sat_l flip(mu_sat_v)];


fig1 = figure;
fig2 = figure;

for k = 1:length(P_vec)
    for i = 1:length(T_vec)
        T(k,i) = T_vec(i);
        P(k,i) = P_vec(k);
        rho(k,i) = CoolProp.PropsSI('D','P',P(k,i),'T',T(k,i),fluid);
        cp(k,i) = CoolProp.PropsSI('C','P',P(k,i),'T',T(k,i),fluid);
        l(k,i) = CoolProp.PropsSI('L','P',P(k,i),'T',T(k,i),fluid);
        mu(k,i) = CoolProp.PropsSI('V','P',P(k,i),'T',T(k,i),fluid);        
    end
    DcpDt(k,:) = [0 diff(cp(k,:))./diff(T(k,:))]; 
    DrhoDt(k,:) = [0 diff(rho(k,:))./diff(T(k,:))];
    DlDt(k,:) = [0 diff(l(k,:))./diff(T(k,:))];
    DmuDt(k,:) = [0 diff(mu(k,:))./diff(T(k,:))];
    
    [max_cp(k), i_max_cp(k)] = max(cp(k,:));
    i_cp_max = i_max_cp(k)-1;
        
    % critical boundaries defined with cp maximum value
    x_frac = 0.07;
    l_cp(k) = x_frac*max_cp(k)+(1-x_frac)*cp(k,1);
    err_l_cp = abs(cp(k,1:i_max_cp(k)) - l_cp(k)); 
    [min_err_l_cp, i_min_err_l_cp] = min(err_l_cp);
    
    v_cp(k) = x_frac*max_cp(k)+(1-x_frac)*cp(k,end);
    err_v_cp = abs(cp(k,i_max_cp(k):end) - v_cp(k)); 
    [min_err_v_cp, i_min_err_v_cp] = min(err_v_cp);
    i_min_err_v_cp = i_max_cp(k)+i_min_err_v_cp;    
    
    
    % critical boundaries defined with derivative on cp
    [max_DcpDt_pos(k), i_max_DcpDt_pos(k)] = max(DcpDt(k,:));
    [min_DcpDt_neg(k), i_min_DcpDt_neg(k)] = min(DcpDt(k,:));
    
    x_frac = 0.07;
    l_DcpDt(k) = x_frac*max_DcpDt_pos(k)+(1-x_frac)*DcpDt(k,2);
    err_l_DcpDt = abs(DcpDt(k,1:i_max_DcpDt_pos(k)) - l_DcpDt(k)); 
    [min_err_l_DcpDt, i_min_err_l_DcpDt] = min(err_l_DcpDt);
    
    v_DcpDt(k) = x_frac*min_DcpDt_neg(k)+(1-x_frac)*DcpDt(k,end);
    err_v_DcpDt = abs(DcpDt(k,i_min_DcpDt_neg(k):end) - v_DcpDt(k)); 
    [min_err_v_DcpDt, i_min_err_v_DcpDt] = min(err_v_DcpDt);
    i_min_err_v_DcpDt = i_min_DcpDt_neg(k)+i_min_err_v_DcpDt;
    
    % critical boundaries defined with derivative on rho
    [min_DrhoDt(k), i_min_DrhoDt(k)] = min(DrhoDt(k,:));
    
    x_frac = 0.01;
    l_DrhoDt(k) = x_frac*min_DrhoDt(k)+(1-x_frac)*DrhoDt(k,2);
    err_l_DrhoDt = abs(DrhoDt(k,1:i_min_DrhoDt(k)) - l_DrhoDt(k));
    [min_err_l_DrhoDt, i_min_err_l_DrhoDt] = min(err_l_DrhoDt);
    
    v_DrhoDt(k) = x_frac*min_DrhoDt(k)+(1-x_frac)*DrhoDt(k,end);
    err_v_DrhoDt = abs(DrhoDt(k,i_min_DrhoDt(k):end) - v_DrhoDt(k));
    [min_err_v_DrhoDt, i_min_err_v_DrhoDt] = min(err_v_DrhoDt);
    i_min_err_v_DrhoDt = i_min_DrhoDt(k)+i_min_err_v_DrhoDt;
    
    
    % Limits definition
    i_l_limit = i_min_err_l_DrhoDt;%i_min_err_l_DcpDt%i_min_err_l_cp; %
    i_v_limit = i_min_err_v_DrhoDt;%i_min_err_v_DcpDt%i_min_err_v_cp; %
    i_pseudo_crit = i_cp_max;
    T_pseudo_l_limit(k) = T(k,i_l_limit);    
    T_pseudo_v_limit(k) = T(k,i_v_limit);    
    T_pseudo_crit(k) = T(k,i_pseudo_crit);
    h_pseudo_l_limit(k) = CoolProp.PropsSI('D','P',P(k,1),'T',T_pseudo_l_limit(k),fluid); %T(k,i_l_limit);    
    h_pseudo_v_limit(k) = CoolProp.PropsSI('D','P',P(k,1),'T',T_pseudo_v_limit(k),fluid); %T(k,i_v_limit);
    
    figure(fig1)
    subplot(2,2,1)
    hold on
    plot(T_sat, rho_sat, 'k--')
    plot(T(k,:),rho(k,:))
    
    plot(T(k,i_l_limit),rho(k,i_l_limit), 'mv')
    plot(T(k,i_v_limit),rho(k,i_v_limit), 'mv')
    xlabel('T')
    ylabel('rho')
    hold off
    
    subplot(2,2,2)
    hold on
    plot(T(k,:),cp(k,:))
    plot(T(k,i_l_limit),cp(k,i_l_limit), 'mv')
    plot(T(k,i_v_limit),cp(k,i_v_limit), 'mv')    
    xlabel('T')
    ylabel('cp')
    hold off
    
    subplot(2,2,3)
    hold on
    plot(T_sat, l_sat, 'k--')
    plot(T(k,:),l(k,:))
    plot(T(k,i_l_limit),l(k,i_l_limit), 'mv')
    plot(T(k,i_v_limit),l(k,i_v_limit), 'mv')    
    xlabel('T')
    ylabel('L')
    hold off
    
    subplot(2,2,4)
    hold on
    plot(T_sat, mu_sat, 'k--')
    plot(T(k,:),mu(k,:))
    plot(T(k,i_l_limit),mu(k,i_l_limit), 'mv')
    plot(T(k,i_v_limit),mu(k,i_v_limit), 'mv')
    xlabel('T')
    ylabel('mu')
    hold off
    
    if 1
        figure(fig2)
        subplot(2,2,1)
        hold on
        plot(T(k,:),DrhoDt(k,:))
        plot(T(k,i_l_limit),DrhoDt(k,i_l_limit), 'mv')
        plot(T(k,i_v_limit),DrhoDt(k,i_v_limit), 'mv')
        xlabel('T')
        ylabel('DrhoDt')
        hold off
        
        subplot(2,2,2)
        hold on
        plot(T(k,:),DcpDt(k,:))
        plot(T(k,i_l_limit),DcpDt(k,i_l_limit), 'mv')
        plot(T(k,i_v_limit),DcpDt(k,i_v_limit), 'mv')
        xlabel('T')
        ylabel('DcpDt')
        hold off
        
        subplot(2,2,3)
        hold on
        plot(T(k,:),DlDt(k,:))
        xlabel('T')
        ylabel('DlDt')
        hold off
        
        subplot(2,2,4)
        hold on
        plot(T(k,:),DmuDt(k,:))
        xlabel('T')
        ylabel('DmuDt')
        hold off
    end
    
end

figure
subplot(2,3,1)
plot(P_vec, T_pseudo_crit,'--o')
ylabel('T pseudo critical')
subplot(2,3,2)
plot(P_vec, T_pseudo_l_limit,'--o')
ylabel('T l limit')
subplot(2,3,3)
plot(P_vec, T_pseudo_v_limit,'--o')
ylabel('T v limit')
subplot(2,3,5)
plot(P_vec, h_pseudo_l_limit,'--o')
ylabel('h l limit')
subplot(2,3,6)
plot(P_vec, h_pseudo_v_limit,'--o')
ylabel('h v limit')