function [cpbar,Q_4,x,e_min,e] = cp_bar_bis(m,n,f,T_p1,T_p2)
%CP_BAR procedure returns cp_bar (J/kg.K) of the combustion products of a CmHn
%hydrocarbon between two temperatures. The inputs are : the type of fuel (m and n),
%the fuel-air ratio (f) and  the two temperatures (T_p1 and T_p2) in °C.  
%
%In the case of incomplete combustion (negative excess air value), the procedure
%returns also : the energy losses by CO (Q_4 in J/kg of fuel), the fraction (x) of C
%atoms oxidized into CO2 and the mininum excess air (e_min) to get x>0.
%
%Simplified air is used. In the case of incomplete combustion, it is also assumed that
%there is no O2 in the combustion products.

%Default values (in the case of complete combustion)
Q_4 = 0;
x = 1;
e_min = 0;
e = 0;

%"Composition de l'air simplifié"
x_air_N2 = 0.79;
x_air_O2= 0.21;

%"Constante universelle des gaz"
R_u = 8314;

%Masses molaires des différents constituants [kg/mol]
MM_CO = CoolProp.PropsSI('M','T',298,'P',1e5,'CO');
MM_N2 = CoolProp.PropsSI('M','T',298,'P',1e5,'N2');
MM_O2 = CoolProp.PropsSI('M','T',298,'P',1e5,'O2');
MM_CO2 = CoolProp.PropsSI('M','T',298,'P',1e5,'CO2');
MM_H2O = CoolProp.PropsSI('M','T',298,'P',1e5,'water');

if (f>0) %produits de combustion et pas de l'air pur
    %f  stoechiométrique, excès d'air et f
    parfuel = m+n/4;
    MM_f = (12*m+1*n)/1000;
    f_st = MM_f/(parfuel*(MM_O2+(x_air_N2/x_air_O2)*MM_N2));
    e = (f_st/f)-1;
    if f>f_st %défaut d'air
        %Combustion stoechiométrique : CmHn+(m+n/4)*(O2+(79/21)*N2
        e_min = -m/(2*parfuel); %pour avoir x>0
        f_max = f_st/(1+e_min);
        %Hypothèse : suffisamment d'air pour oxydation de H en H2O et pas d'O2 dans les fumées
        %D'où : CmHn+(m+n/4)*(O2+(79/21)*N2)*(1+e) = a*CO2+b*CO+c*H2O+g*N2
        g = parfuel*(1+e)*(79/21);
        c = n/2;
        x=2*(e*parfuel+m/2)/m;
        a = m*x; %x est le nombre d'atomes de C oxydés en CO2
        b = m*(1-x);
        %Fractions molaires des constituants
        nkmoltot = a+ b + c +g;
        x_CO2 = a/nkmoltot;
        x_CO = b/nkmoltot;
        x_H2O = c/nkmoltot;
        x_N2 = g/nkmoltot;
        mmprod = x_CO2*MM_CO2+x_CO*MM_CO+x_H2O*MM_H2O+x_N2*MM_N2;
        %c_p moyen des produits de combustion
        h_CO2 = (CoolProp.PropsSI('H','P',1e5,'T',T_p1+273.15,'CO2')-CoolProp.PropsSI('H','P',1e5,'T',T_p2+273.15,'CO2'))*MM_CO2;
        h_CO = (CoolProp.PropsSI('H','P',1e5,'T',T_p1+273.15,'CO')-CoolProp.PropsSI('H','P',1e5,'T',T_p2+273.15,'CO'))*MM_CO;
        a_H20 = -1.34618e7;
        b_H20 = 1785.27;
        c_H20 = 0.354129;
        d_H20 = -0.00000299621e-9;
        h_H2O = ((a_H20 +b_H20*T_p1 + c_H20*T_p1^2 + d_H20*T_p1^3) - (a_H20 +b_H20*T_p2 + c_H20*T_p2^2 + d_H20*T_p2^3))*MM_H2O;
        h_N2 = (CoolProp.PropsSI('H','P',1e5,'T',T_p1+273.15,'N2')-CoolProp.PropsSI('H','P',1e5,'T',T_p2+273.15,'N2'))*MM_N2;
        cpbar = ((x_CO2*h_CO2+x_CO*h_CO+x_H2O*h_H2O+x_N2*h_N2)/mmprod)/(T_p1-T_p2);
        %Pertes par imbrûlés par kmol de combustible
        T_ref=25;
        LHV_CO_mol = CoolProp.PropsSI('H','P',1e5,'T',T_ref+273.15,'CO')*MM_CO+0.5*CoolProp.PropsSI('H','P',1e5,'T',T_ref+273.15,'O2')*MM_O2-CoolProp.PropsSI('H','P',1e5,'T',T_ref+273.15,'CO2')*MM_CO2;
        LHV_CO_mas = LHV_CO_mol/MM_CO;
        Q_4 = LHV_CO_mol*b;
    else
        %Fractions molaires des constituants
        nkmoltot = m+ n/2 + parfuel*(e+(1+e)*(x_air_N2/x_air_O2));
        x_CO2 = m/nkmoltot;
        x_H2O = (n/2)/nkmoltot;
        x_O2 = parfuel*e/nkmoltot;
        x_N2 = parfuel*(1+e)*(x_air_N2/x_air_O2)/nkmoltot;
        %Masse molaire des produits de combustion
        MM_prod = x_CO2*MM_CO2+x_H2O*MM_H2O+x_O2*MM_O2+x_N2*MM_N2;
        %Enthalpie et chaleur massique des produits de combustion ; /kg de produits
        hp_1 = (CoolProp.PropsSI('H','P',1e5,'T',T_p1+273.15,'CO2')-CoolProp.PropsSI('H','P',1e5,'T',T_p2+273.15,'CO2'))*MM_CO2;
        a_H20 = -1.34618e7;
        b_H20 = 1785.27;
        c_H20 = 0.354129;
        d_H20 = -0.00000299621e-9;
        hp_2 = ((a_H20 +b_H20*T_p1 + c_H20*T_p1^2 + d_H20*T_p1^3) - (a_H20 +b_H20*T_p2 + c_H20*T_p2^2 + d_H20*T_p2^3))*MM_H2O;
        hp_3 = (CoolProp.PropsSI('H','P',1e5,'T',T_p1+273.15,'O2')-CoolProp.PropsSI('H','P',1e5,'T',T_p2+273.15,'O2'))*MM_O2;
        hp_4 = (CoolProp.PropsSI('H','P',1e5,'T',T_p1+273.15,'N2')-CoolProp.PropsSI('H','P',1e5,'T',T_p2+273.15,'N2'))*MM_N2;
        hp = (x_CO2*hp_1+x_H2O*hp_2+x_O2*hp_3+x_N2*hp_4)/MM_prod;
        cpbar = hp/(T_p1-T_p2);
    end
end
end

% FUNCTION gamma (m,n,f,T_p1,T_p2)
% 
% {$GAMMA
% GAMMA  function returns the specific heats ratio of the combustion products of a CmHn
% hydrocarbon between two temperatures. The inputs are : the type of
% fuel (m and n), the fuel-air ratio (f) and  the two temperatures (T_p1 and T_p2)
% in °C. 
% 
% Simplified air is used. In the case of incomplete combustion, it is also assumed that
% there is no O2 in the combustion products.
% 
% Moreover, the function returns the specific heats for pure air if f=0 (whatever m and n).
% 
% When using this function, don't forget to set mass basis calculation in the unit system
% settings.}
% 
% "Check of the units setting"
% If (unitsystem('K')=1)  THEN CALL error('Please, set °C for temperature units!')
% IF (unitsystem('Molar')=1)  then CALL error('Please, set Mass basis!')
% 
% "Composition de l'air simplifié"
% x_air_N2 = 0.79 ; x_air_O2= 0.21 
% 
% "Constante universelle des gaz"
% R_u = 8314 "J/kmol.K"
% 
% "Masses molaires des différents constituants"
% MM_CO = molarmass(CO) ; MM_N2 = molarmass(N2); MM_O2 = molarmass(O2)
% MM_CO2 = molarmass(CO2) ; MM_H2O = molarmass(H2O)
% 
% 
% IF (f>0) Then
% 
% "Combustion stoechiométrique : CmHn+(m+n/4)*(O2+(79/21)*N2"
% parfuel = m+n/4
% MM_f = 12*m+1*n
% f_st = MM_f/(parfuel*(MM_O2+(x_air_N2/x_air_O2)*MM_N2))
% e = (f_st/f)-1
% 
% IF (f>f_st) Then 
% 
% "Hypothèse : suffisamment d'air pour oxydation de H en H2O et pas d'O2 dans les fumées" 
% "D'où : CmHn+(m+n/4)*(O2+(79/21)*N2)*(1+e) = a*CO2+b*CO+c*H2O+g*N2"
% g = parfuel*(1+e)*(79/21)
% c = n/2
% x=2*(e*parfuel+m/2)/m
% a = m*x "x est le nombre d'atomes de C oxydés en CO2"
% b = m*(1-x)
% "Fractions molaires des constituants" 
% nkmoltot = a+ b + c +g
% x_CO2 = a/nkmoltot
% x_CO = b/nkmoltot
% x_H2O = c/nkmoltot
% x_N2 = g/nkmoltot
% MM_prod = x_CO2*MM_CO2+x_CO*MM_CO+x_H2O*MM_H2O+x_N2*MM_N2
% "!c_p moyen des produits de combustion"
% h_CO2 = (enthalpy(CO2,T=T_p1)-enthalpy(CO2,T=T_p2))*MM_CO2
% h_CO = (enthalpy(CO,T=T_p1)-enthalpy(CO,T=T_p2))*MM_CO
% h_H2O =(enthalpy(H2O,T=T_p1)-enthalpy(H2O,T=T_p2))*MM_H2O
% h_N2 = (enthalpy(N2,T=T_p1)-enthalpy(N2,T=T_p2))*MM_N2
% c_bar_p = ((x_CO2*h_CO2+x_CO*h_CO+x_H2O*h_H2O+x_N2*h_N2)/MM_prod)/(T_p1-T_p2)
% gamma = c_bar_p/(c_bar_p-R_u/MM_prod)
% 
% ELSE
% 
% "Fractions molaires des constituants" 
% nkmoltot = m+ n/2 + parfuel*(e+(1+e)*(x_air_N2/x_air_O2))
% x_CO2 = m/nkmoltot
% x_H2O = (n/2)/nkmoltot
% x_O2 = parfuel*e/nkmoltot
% x_N2 = parfuel*(1+e)*(x_air_N2/x_air_O2)/nkmoltot
% 
% "Masse molaire des produits de combustion"
% MM_prod = x_CO2*MM_CO2+x_H2O*MM_H2O+x_O2*MM_O2+x_N2*MM_N2
% 
% "Enthalpie, chaleur massique et gamma des produits de combustion ; /kg de produits"
% hp_1 = (enthalpy(CO2,T=T_p1)-enthalpy(CO2, T=T_p2))*MM_CO2
% hp_2 = (enthalpy(H2O,T=T_p1)-enthalpy(H2O, T=T_p2))*MM_H2O
% hp_3 = (enthalpy(O2,T=T_p1)-enthalpy(O2, T=T_p2))*MM_O2
% hp_4 = (enthalpy(N2,T=T_p1)-enthalpy(N2, T=T_p2))*MM_N2
% hp = (x_CO2*hp_1+x_H2O*hp_2+x_O2*hp_3+x_N2*hp_4)/MM_prod
% cpbar  = hp/(T_p1-T_p2)
% gamma = cpbar/(cpbar-R_u/MM_prod)
% ENDIF
% 
% ELSE
% 
% "Cas de l'air pur"
% MM_a = x_air_O2*MM_O2+x_air_N2*MM_N2
% ha_1 = (enthalpy(O2,T=T_p1)-enthalpy(O2, T=T_p2))*MM_O2
% ha_2 = (enthalpy(N2,T=T_p1)-enthalpy(N2, T=T_p2))*MM_N2
% ha = (x_air_O2*ha_1+x_air_N2*ha_2)/MM_a
% cpbar  = ha/(T_p1-T_p2)
% gamma = cpbar/(cpbar-R_u/MM_a)
% ENDIF
% 
% END
