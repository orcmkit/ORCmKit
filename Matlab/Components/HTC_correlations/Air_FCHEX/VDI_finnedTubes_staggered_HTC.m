function  [hConv, Nu, flag] = VDI_finnedTubes_staggered_HTC(mu, Pr, k, G, Dh, omega_t, disp_flag)% verified

% air-side HTC proposed by Schmidt in: VDI Heat Atlas, section M1 - 2.6, page 1275
% Warning, Dh is the externel diameter of the tubes forming the bank !!!!!

% RDickes - 25/07/2018

Re_min = 1000; Re_max = 1e5;
omega_t_min = 5; omega_t_max = 30;
flag = 1;

Re = G*Dh/mu;
Nu = 0.38*Re^0.6*Pr^0.33333333*omega_t^(-0.15);
hConv = Nu*k/Dh;

if Re >= Re_max || Re <= Re_min
    if disp_flag
        display(['VDI FCHEX singe-phase: Out of validity range --> Re = ' num2str(Re) ' is out of [' num2str(Re_min) ' - ' num2str(Re_max) '] !!!'])
    end
    flag = flag -1;
end

if omega_t >= omega_t_max || omega_t <= omega_t_min
    if disp_flag
        display(['VDI FCHEX singe-phase: Out of validity range --> omega_t = ' num2str(omega_t) ' is out of [' num2str(omega_t_min) ' - ' num2str(omega_t_max) '] !!!'])
    end
    flag = flag -2;
end

end
