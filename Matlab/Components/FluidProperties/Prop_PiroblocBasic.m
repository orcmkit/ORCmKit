clear all; close all; clc
load('PiroblocBasic_PropTab.mat')

T_C_tab = PiroblocBasic_PropTab(:,1);
T_K_tab = T_C_tab + 273.15;
nu_tab = PiroblocBasic_PropTab(:,2)*1e-6; % m^2/s
rho_tab = PiroblocBasic_PropTab(:,3)*1e3; %kg/m³
cp_tab = PiroblocBasic_PropTab(:,4)*4186.8; %J/kg.K
k_tab = PiroblocBasic_PropTab(:,5)*1.163; %W/m.°C
Pr_tab = PiroblocBasic_PropTab(:,6);
mu_tab = nu_tab.*rho_tab;
Pr_tab2 = mu_tab.*cp_tab./k_tab;
T_C_extrap = linspace(20,350,100);
T_K_extrap = T_C_extrap+273.15;
for ii = 1:length(T_K_extrap)
    rho_prev(ii,1) = PropsSI_ICP('D', 'T', T_K_extrap(ii), 'P', 1e5, 'ICP_PiroblocBasic');
    cp_prev(ii,1) = PropsSI_ICP('C', 'T', T_K_extrap(ii), 'P', 1e5, 'ICP_PiroblocBasic');
    k_prev(ii,1) = PropsSI_ICP('L', 'T', T_K_extrap(ii), 'P', 1e5, 'ICP_PiroblocBasic');
    mu_prev(ii,1) = PropsSI_ICP('V', 'T', T_K_extrap(ii), 'P', 1e5, 'ICP_PiroblocBasic');
    Pr_prev(ii,1) = mu_prev(ii,1).*cp_prev(ii,1)./k_prev(ii,1);
    
    rho_new(ii,1) = PropsSI_ICP('D', 'T', T_K_extrap(ii), 'P', 1e5, 'ICP_PiroblocBasic_2');
    cp_new(ii,1) = PropsSI_ICP('C', 'T', T_K_extrap(ii), 'P', 1e5, 'ICP_PiroblocBasic_2');
    k_new(ii,1) = PropsSI_ICP('L', 'T', T_K_extrap(ii), 'P', 1e5, 'ICP_PiroblocBasic_2');
    mu_new(ii,1) = PropsSI_ICP('V', 'T', T_K_extrap(ii), 'P', 1e5, 'ICP_PiroblocBasic_2');
    Pr_new(ii,1) = mu_new(ii,1).*cp_new(ii,1)./k_new(ii,1);
end

figure
subplot(2,3,1)
hold on
plot(T_C_tab, rho_tab, 'o')
plot(T_C_extrap, rho_prev, '-')
plot(T_C_extrap, rho_new, '--')
hold off
ylabel('rho')
rho_fit = fit(T_C_tab, rho_tab, 'poly1')
disp(['p1 = ' num2str(rho_fit.p1,10)])
disp(['p2 = ' num2str(rho_fit.p2,10)])


subplot(2,3,2)
hold on
plot(T_C_tab, cp_tab, 'o')
plot(T_C_extrap, cp_prev, '-')
plot(T_C_extrap, cp_new, '--')
hold off
ylabel('cp')
cp_fit = fit(T_K_tab, cp_tab, 'poly1')
disp(['p1 = ' num2str(cp_fit.p1,10)])
disp(['p2 = ' num2str(cp_fit.p2,10)])

subplot(2,3,3)
hold on
plot(T_C_tab, k_tab, 'o')
plot(T_C_extrap, k_prev, '-')
plot(T_C_extrap, k_new, '--')
hold off
ylabel('k')
k_fit = fit(T_C_tab, k_tab, 'poly1')
disp(['p1 = ' num2str(k_fit.p1,10)])
disp(['p2 = ' num2str(k_fit.p2,10)])


subplot(2,3,4)
hold on
plot(T_C_tab, mu_tab, 'o')
plot(T_C_extrap, mu_prev, '-')
plot(T_C_extrap, mu_new, '--')
hold off
ylabel('mu')
mu_fit = fit(T_C_tab, mu_tab, 'power1')
disp(['a = ' num2str(mu_fit.a,10)])
disp(['b = ' num2str(mu_fit.b,10)])


subplot(2,3,5)
hold on
plot(T_C_tab, Pr_tab, 'o')
plot(T_C_tab, Pr_tab2, 'o')
plot(T_C_extrap, Pr_prev, '-')
plot(T_C_extrap, Pr_new, '--')
hold off
ylabel('Pr')