# -*- coding: utf-8 -*-
"""
2016 Purdue Conferences
Herrick Laboratories, Purdue University, West Lafayette, IN (USA) 

Short Course "Oil Management in Compressors and their Systems"

author: Davide Ziviani
email: dziviani@purdue.edu

v1.0: created July 8th 2016

For illustrational purposes only. If used for research please reference the model as:

Ziviani, D., Woodland, B.J., georges, E., Groll, E.A., Braun, J.E., Horton, W.T., van den Broek, M., De Paepe, M., "Development and a validation of a charge sensitive organic Rankine cycle (ORC) Simulation Tool", Energies, 9(2016), 1-36.
 
"""

#Import CoolProp
from CoolProp.CoolProp import PropsSI
#Python standard imports
import numpy as np
import math
import time
import warnings
import pylab
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.mlab as ml

warnings.simplefilter("ignore",RuntimeWarning)
from random import randint, random


numTol = 1e-5 #numerical tolerance





"""
Expander model based on semi-empirical model with flooding. Assumed homodeneous flow thorugh 


Required Parameters:
    
===========   ==========  ========================================================================
Variable      Units       Description
===========   ==========  ========================================================================
Ref                   
T_amb             varied      
Ref                N/A      A string representing the refrigerant
y               N/A         Flooding fraction
T_amb 
T_suc           [K]
P_suc           [kPa]
P_ex            [kPa]
N_exp           [rpm]
UA_suc_n        [W/K]
UA_ex_n         [W/K]
UA_amb_n        [W/K]
M_dot_n         [kg/s]
V_ratio 
A_suc           [m2]
A_leak          [m2]
x_leak                      Liquid fraction leakages
V_suc_exp       [m3]
T_loss          [Nm]
A_ex            [m2]
===========   ==========  ========================================================================    
All variables are of double-type unless otherwise specified
    
"""


#=============================================================================
"The Oil Properties"

def c_l(T):
    "specific heat [kJ/kg-K] given T in K"
    c = (2.0935*T + 1186.7)/1000 #POE 150 SUS
    return c
    
def u_l(T):
    "internal energy [kJ/kg] given T in K"
    u = (2.0935/2*T**2 + 1186.7*T)/1000 #POE 150 SUS
    return u

def rho_l(T):
    "density [kg/m^3] given T in K"
    rho = -0.7*T+1186 #POE 150 SUS
    return rho

def s_l(T):
    "specific entropy [kJ/kg-K] given T in K"  
    s = (2.0935*(T-298) + 1186.7*np.log(T/298.0))/1000 # POE 150 SUS
    return s

def h_l(T,P):
    "the specific enthalpy [kJ/kg-k]"
    h = u_l(T) + P/rho_l(T)
    return h

def mu_l(T):
    "returns the viscosity given temp in K, [Pa-s]"
    #POE 150 SUS equation only valid from 60 C to 120 C
    mu = 0.000000000517*T**4 - 0.000000795840*T**3 + 0.000460766590*T**2 - 0.118976538068*T + 11.571730524692 
    return mu
    
def k_l(T):
    "Thermal conductivity [W/m-K] for T in Kelvin"
    k = 0.138 #POE 150 SUS
    return k      


def T_l(h,P):
    "find the liquid temperature [K] given h and P"
    #initial guesses, functions fairly linear so guesses not too important
    T = [300, 400]
    h_check = h_l(T[0],P)
    f = [abs(h_check - h)] #function to converge
    
    i=0 #array index
    while abs(f[i])> numTol:
        i=i+1 #update index
        h_check = h_l(T[i],P)
        f = np.append(f,abs(h_check - h))
        T = np.append(T , T[i]-(f[i]*(T[i]-T[i-1]))/(f[i]-f[i-1])) #secant method
        
        if i>100:
            raise Exception("T_l not converging after 100x")
        
    return T[i]


def vol_l(T):
    "return specific volume of liq"
    v = 1.0/rho_l(T)
    return v

def vol_g(gas,T,P):
    "returns the specific volume"
    v = 1.0/(PropsSI('D','T',T,'P',P*1000,gas))
    return v
        
#==========================================================================    
"Mixture Function Calls"

def T_mix_h(gas,h,P,y,Tguess):
    "find the mixture temperature given its specific enthalpy"
    T=[Tguess+5, Tguess-5]
    h_check = (PropsSI('H','T',T[0],'P',P*1000,gas)/1000+ y*h_l(T[0],P))/(1+y)
    f =[abs(h_check - h)]
    
    i=0 
    while abs(f[i])> 1e-5: 
        i=i+1 
        if T[i]< Tguess/2: 
            T[i] = Tguess/2 + random()*(Tguess - Tguess/2)
        h_check = (PropsSI('H','T',T[i],'P',P*1000,gas)/1000+ y*h_l(T[i],P))/(1+y)
        f = np.append(f,abs(h_check - h))
        T = np.append(T , T[i]-(f[i]*(T[i]-T[i-1]))/(f[i]-f[i-1])) #secant method
        
        if i>150:
            raise Exception("T_mix_h not converging after 150x")
    
    return T[i]


def T_mix_s(gas,s,P,y,Tguess):
    "find the mixture temperature given its specific enthalpy"
    T=[Tguess+5, Tguess-5]
    s_check = (PropsSI('S','T',T[0],'P',P*1000,gas)/1000+ y*s_l(T[0]))/(1+y)
    f =[abs(s_check - s)]
    
    i=0 
    while abs(f[i])> numTol:
        i=i+1 
        if T[i]< Tguess/2: 
            T[i] = Tguess/2 + random()*(Tguess - Tguess/2)
            
        s_check = (PropsSI('S','T',T[i],'P',P*1000,gas)/1000+ y*s_l(T[i]))/(1+y)
        f = np.append(f,abs(s_check - s))
        T = np.append(T , T[i]-(f[i]*(T[i]-T[i-1]))/(f[i]-f[i-1])) #secant method
         
        if i>100:
            raise Exception("T_mix_s not converging after 100x")
    
    return T[i]
    

def gamma_mix(gas,T,P,y):
    'find the ratio of specific heats of the oil-refrigerant mixture'
    Cp = (PropsSI('C','T',T,'P',P*1000,gas)/1000+y*c_l(T))/(1+y)
    Cv = (PropsSI('O','T',T,'P',P*1000,gas)/1000+y*c_l(T))/(1+y)
    gamma = Cp/Cv
    return gamma


def Re_mix(gas,T,P,vel,D,y):
    "find the Reynolds number of the oil-refrigerant mixture"
    vol = (vol_g(gas,T,P)+y*vol_l(T))/(1+y)
    rho = 1.0/vol

    "McAdams 1942 eqn"
    muG = PropsSI('V','T',T,'P',P*1000,gas)
    muL = mu_l(T)
    xg = 1.0/(1.0+y)
    xl=y/(1.0+y)    
    mu = muG*muL/(muG*xl+muL*xg)
    Re = rho*vel*D/mu
    
    return Re
    
   
def Pr_mix(gas,T,P,y):
    "find the Prandtl number of the oil-refrigerant mixture"
    cp = (PropsSI('C','T',T,'P',P*1000,gas)+y*c_l(T))/(1+y)*1000 #need [J/kg-K]
   
    "mu = two-phase flow through a frictional pipe (Shannak 2008)"
    muG = PropsSI('V','T',T,'P',P*1000,gas)
    muL = mu_l(T)
    rhoG = PropsSI('D','T',T,'P',P*1000,gas)
    rhoL = rho_l(T)
    x = 1.0/(1.0+y)
    mu = ( muG*x + muL*(1-x)*(rhoG/rhoL) )
    k=k_mix(gas,T,P,y)
    Pr = cp*mu/k
    
    return Pr
    

def k_mix(gas,T,P,y):
    "find the thermal conductivity of the oil-refrigerant mixture"
    xg = 1.0/(1+y)
    xl=y/(1+y)
    
    "Reference: Bell, I. PhD Dissertation 2011 p.134"
    vG = vol_g(gas,T,P)
    vL = vol_l(T)
    void = xg*vG/(xl*vL + xg*vG)
    kG = PropsSI('L','T',T,'P',P*1000,gas)#[W/m-K]
    kL = k_l(T)
    k =(1.0-void)*kL + void*kG
    
    return   k 

#======================================================================   
"Selected Modelling steps in functions for overall clarity"

def gamma_r(gas,P_1,T_1,P_2):
    s_1=PropsSI('S','P',P_1,'T',T_1,gas)
    v_1=1/PropsSI('D','P',P_1,'T',T_1,gas)
    v_2=1/PropsSI('D','P',P_2,'S',s_1,gas)
    "P_1*v_1^gamma=P_2*v_2^gamma"
    r_p=P_1/P_2
    r_v=v_2/v_1
    gamma_r=log10(r_p)/log10(r_v)
    #gamma_tris=ln(r_p)/ln(r_v)
    return gamma_r

def SuctionNozzle(gas,Area,mg,ml,P1,T1):
    y = ml/mg
    h1 = (PropsSI('H','T',T1,'P',P1*1000,gas)/1000+y*h_l(T1,P1))/(1+y)
    s1 = (PropsSI('S','T',T1,'P',P1*1000,gas)/1000+y*s_l(T1))/(1+y)
    
    gamma = gamma_mix(gas,T1,P1,y)
    P_Crit = P1*pow((2.0/(gamma+1)),(gamma/(gamma-1)))  #serve as bound on iterations
    PHigh = P1-1
    PLow = P_Crit
    
    FHigh = Suction_helper(gas,h1,s1,T1,PHigh,Area,mg,ml,y)
    FLow = Suction_helper(gas,h1,s1,T1,PLow,Area,mg,ml,y)  
    
    if FLow>FHigh:
        signChange = -1.0
    else:
        signChange=1.0
    
    FHigh = FHigh*signChange
    FLow = FLow*signChange

    NumSections = 0
    FMid = 10
    
    while abs(FMid) > 1e-4:
       
        NumSections = NumSections+1
        PMid = (PHigh+PLow)/2.0
        FMid = signChange*Suction_helper(gas,h1,s1,T1,PMid,Area,mg,ml,y) 
        
        if FMid>0:
            PHigh = PMid
        else:
            PLow = PMid
                
        if NumSections>100:
            raise Exception("Nozzle not converging after 100x")
    
    P2 = PMid    
    
    "recover static enthalpy in isobaric diffuser"
    T2 = T_mix_h(gas,h1,P2,y,T1)

    return [T2,P2]

def Suction_helper(gas,h1,s1,T1,P2,Area,mg,ml,y):
    
    T2 = T_mix_s(gas,s1,P2,y,T1) #isentropic nozzle
    h2 = (PropsSI('H','T',T2,'P',P2*1000,gas)/1000+y*h_l(T2,P2))/(1+y)
    v2 =(vol_g(gas,T2,P2) +y*vol_l(T2))/(1+y)
    vel2 = (mg+ml)/(Area/v2)
    KE2 = 0.5*vel2**2.0/1000.0
    Ebal = h1 - (h2 + KE2)        
    
    return Ebal 



def IsentropicExpansion(gas,s1,T1,P1,v2,v_ratio,y,P_ex):
    "Perform isentropic expansion along the volume ratio of the machine"
    
    rho2 = 1.0/v2
    "guess outlet pressure"
    P= [P1/v_ratio,0.95*P1/v_ratio]
    
    i=0 
    f = [10.0000] 
            
    while abs(f[i]) > numTol:        
        if f[i] != 10.0000: #only increase index after 1st iteration
            i=i+1
        
        T = T_mix_s(gas,s1,P[i],y,T1)
        
        "Mass (volume) Balance"
        v_out = (vol_g(gas,T,P[i]) +y*vol_l(T))/(1+y)
        rho_out = 1.0/v_out
        Mbal = rho2 - rho_out
        
        if i > 0:
            f = np.append(f, Mbal)              
            P = np.append(P, P[i] - (f[i]*(P[i]-P[i-1]))/(f[i]-f[i-1]))            
        else:                    
            f = [Mbal]         
        if i>50:                  
            raise Exception("IsenExp-P not converging after 50x")
            
    P = P[i]

    return [T,P]

        
def leakage(gas,T_pleleak,P_preleak,P_postleak,y_star,A_leak):

    gamma_leak = gamma_mix(gas,T_pleleak,P_preleak,y_star)
    P_leak_crit = P_preleak*pow((2.0/(gamma_leak+1)),gamma_leak/(gamma_leak-1))
    P_leak_thr = max(P_postleak,P_leak_crit) #choke flow in leakage
    
    s_preleak = (PropsSI('S','T',T_pleleak,'P',P_preleak*1000,gas)/1000+y_star*s_l(T_pleleak))/(1+y_star)
    h_preleak = (PropsSI('H','T',T_pleleak,'P',P_preleak*1000,gas)/1000+y_star*h_l(T_pleleak,P_preleak))/(1+y_star)
    
    T_leak_thr = PropsSI('T','S',s_preleak*1000,'P',P_leak_thr*1000,gas) #isentropic nozzle leaking
    v_leak_thr = 1.0/PropsSI('D','S',s_preleak*1000,'P',P_leak_thr*1000,gas)
    h_leak_thr = PropsSI('H','S',s_preleak*1000,'P',P_leak_thr*1000,gas)/1000    
    vel_leak_thr = math.sqrt(2.0*(h_preleak - h_leak_thr)*1000)  #kJ/kg --> J/kg
    m_leak = 1.0/v_leak_thr*vel_leak_thr*A_leak

    return m_leak


def ExhaustNozzle(gas,Area,mg,ml,P2,T2,P_suc):
    "includes isobaric diffuser"
    
    y = ml/mg
    h2 = (PropsSI('H','T',T2,'P',P2*1000,gas)/1000+y*h_l(T2,P2))/(1+y)
    
    "guess KE_out"
    vel = [10, 12] #m/s
    
    i=0
    f = [10.0000]
    
    while abs(f[i])> numTol:
        if f[i] != 10.0000: 
            i=i+1
        
        KE_guess = 0.5*vel[i]**2.0/1000.0
        
        "isobaric diffuser"
        h2_static = h2 - KE_guess
        P2_static = P2
        T2_static = T_mix_h(gas,h2_static,P2_static,y,T2)
        s2_static = (PropsSI('S','T',T2_static,'P',P2_static*1000,gas)/1000+y*s_l(T2_static))/(1+y)
        
        v2_static = (vol_g(gas,T2_static,P2_static) +y*vol_l(T2_static))/(1+y)
        vel2 = (mg+ml)/(Area/v2_static)
        
        vel_bal = vel2 - vel[i]
        
        if i > 0:
            f = np.append(f, vel_bal)              
            vel = np.append(vel, vel[i] - (f[i]*(vel[i]-vel[i-1]))/(f[i]-f[i-1]))            
        else:                    
            f = [vel_bal]         
        if i>30:                       
            raise Exception("Exhaust Nozzle vel not converging after 30x")   
        
        
    "guess nozzle inlet pressure"
    P1 = [P2+10,P2+50]
    
    j=0
    g=[10.0000]
    
    while abs(g[j])> numTol:
        if g[j] != 10.0000: 
            j=j+1
    
        T1 = T_mix_s(gas,s2_static,P1[j],y,T2) #isentropic nozzle
        h1 = (PropsSI('H','T',T1,'P',P1[j]*1000,gas)/1000+y*h_l(T1,P1[j]))/(1+y)
        
        Ebal = h1 - h2
        
        if j > 0:
            g = np.append(g, Ebal)              
            P1 = np.append(P1, P1[j] - (g[j]*(P1[j]-P1[j-1]))/(g[j]-g[j-1]))            
        else:                    
            g = [Ebal]         
        if j>30:                       
            raise Exception("Exhaust Nozzle P not converging after 30x") 
        
    return [T1,P1[j]]


def Qex_reverse(gas,y,epsilon_ex,C_dot_ex,M_dot,Tw,h_ex,P_ex,T_ex):
    
    "Qex = epsilon*Cdot*(Tex1 - Tw) = mDot*(hex1-hex)"
    "h(Tex1,Pex) - epsilon*Cdot/mDot*Tex1 = -epsilon*Cdot/mDot*Tw + hex"    
    
    balance = -epsilon_ex*C_dot_ex/M_dot*Tw + h_ex
    
    T_ex1 = [T_ex+3, T_ex+5]
    
    i=0
    f = [10.0000]
    
    while abs(f[i]) > numTol:
        
        if f[i] != 10.0000: 
            i=i+1
        
        bal = (PropsSI('H','T',T_ex1[i],'P',P_ex*1000,gas)/1000+y*h_l(T_ex1[i],P_ex))/(1+y) - epsilon_ex*C_dot_ex/M_dot*T_ex1[i]
        
        Ebal = bal - balance
        
        if i > 0:
            f = np.append(f, Ebal)              
            T_ex1 = np.append(T_ex1, T_ex1[i] - (f[i]*(T_ex1[i]-T_ex1[i-1]))/(f[i]-f[i-1]))            
        else:                    
            f = [Ebal]         
        if i>30:
            raise Exception("Qex-reverse converging after 30x")        
    
    return T_ex1[i]


# Main Scroll Compressor Solver ==============================================================

def Expander(gas,T_amb,T_suc,P_suc,P_ex,y,N_exp,UA_suc_n,UA_ex_n,UA_amb_n,M_dot_n,V_ratio,A_suc,A_leak,x_leak,V_suc_comp,T_loss,W_loss_0,A_ex):


    "Solution method for oil-flooded semi-emperical scroll compressor "
    start_time = time.time()
    
    N_exp = N_exp*1.0 #convert int to double

    "Friction Losses"
    W_loss_sh =  W_loss_0/1000 +2.0*math.pi*(N_exp/60)*T_loss/1000 #W --> kW 

    W_loss_el = 0
    W_loss = W_loss_sh + W_loss_el
    
    "generate a starting guess for massflow iterations"
    mguess = (V_suc_comp/V_ratio)*(N_exp/60)/((vol_g(gas,T_suc,P_suc) +y*vol_l(T_suc))/(1+y))
    mLeak_guess = leakage(gas,T_suc,P_suc,P_ex,y,A_leak)
    mguess = mguess + mLeak_guess
    
    M_dot =  [mguess,0.9*mguess]

    i=0 
    f = [10.0000] 
    
    while abs(f[i]) > numTol:        

        if f[i] != 10.0000: 
            i=i+1
        
        "ride herd on massflow"
        if M_dot[i] < 0:
            print "flow herd?"
            M_dot[i] = mguess - random()*(mguess/3)

        "(1) Inlet Conditions"
        h_suc = (PropsSI('H','T',T_suc,'P',P_suc*1000,gas)/1000+y*h_l(T_suc,P_suc))/(1+y)
        s_suc = (PropsSI('S','T',T_suc,'P',P_suc*1000,gas)/1000+y*s_l(T_suc))/(1+y)

        "(2) suction Nozzle pressure drop"        
        nozzle_start_time = time.time()
        nozzlePass=False
        nozReduce=-1
        while nozzlePass == False:   
            nozReduce = nozReduce+1
            "gas and liquid massflow rates"
            m_g = M_dot[i]/(1+y)
            m_l = M_dot[i] - m_g
      
            if time.time() - nozzle_start_time > 5.0: #kill runs taking too long (s)"
                raise Exception("Time's up")            
            
            try:
                "this won't solve if the mass flow rate is too high"
                [T_suc1,P_suc1] = SuctionNozzle(gas,A_suc,m_g,m_l,P_suc,T_suc)
                h_suc1 = (PropsSI('H','T',T_suc1,'P',P_suc1*1000,gas)/1000+y*h_l(T_suc1,P_suc1))/(1+y)
                s_suc1 = (PropsSI('S','H',h_suc1*1000,'P',P_suc1*1000,gas)/1000+y*s_l(T_suc1))/(1+y)
                nozzlePass = True
            except:
                #gradually decrease till a sufficient mass flow range is found
                M_dot[i] = M_dot[i]*(0.95-nozReduce*0.02) 
                
        "Actual UA Values"
        UA_suc = UA_suc_n*pow((M_dot[i]/M_dot_n),0.8)/1000  #W --> kW
        UA_ex = UA_ex_n*pow((M_dot[i]/M_dot_n),0.8)/1000  #W --> kW
        UA_amb = UA_amb_n/1000 #W --> kW
        
        
        "Generate starting guesses for Work"
        try:
            "guess isentropic step"
            V_dot_suc_exp = V_suc_comp/V_ratio*(N_exp/60.0)
            v_ex3 = (V_dot_suc_exp*V_ratio)/(M_dot[i]-mLeak_guess)
            [T_ex3,P_ex3] = IsentropicExpansion(gas,s_suc1,T_suc1,P_suc1,v_ex3,V_ratio,y,P_ex)
            h_ex3 = (PropsSI('H','T',T_ex3,'P',P_ex3*1000,gas)/1000+y*h_l(T_ex3,P_ex3))/(1+y)
            w_exp_int_s = h_suc1 - h_ex3        

            "guess const volume step"
            P_ex2 = P_ex 
            w_exp_v = v_ex3*(P_ex3 - P_ex2)        
            Wexp_guess = (M_dot[i]-mLeak_guess)*(w_exp_int_s+w_exp_v)-W_loss

        except:
            T_ex3 = T_suc - 20
            Wexp_guess = 1.0
           
        W_exp = [Wexp_guess,0.9*Wexp_guess]

        k = 0 
        L = [10.0000]
        
        while abs(L[k]) > numTol:

            if L[k] != 10.0000:
                k=k+1

            "Guess expander body wall temperature"
            TwGuess = W_loss/UA_amb + T_amb #assumes Q_suc ~ Q_ex
            Tw = [TwGuess,TwGuess-5] #2nd guess slightly lower since typically Q_ex > Q_suc (for expander)

            j=0 
            g = [10.0000] 
            
            while abs(g[j]) > 1e-4: 

                if g[j] != 10.0000: 
                    j=j+1

                "Ambient Heat Loss"
                Q_dot_amb = UA_amb*(Tw[j] - T_amb)
                
                "Envelope Energy Balance"
                #hex*m+Qamb+Wexp = hsuc*m
                h_ex = (h_suc*M_dot[i] - W_exp[k] - Q_dot_amb) / M_dot[i]
                T_ex = T_mix_h(gas,h_ex,P_ex,y,T_ex3)
                
                [T_ex0,P_ex0] = ExhaustNozzle(gas,A_ex,m_g,m_l,P_ex,T_ex,P_suc)
                h_ex0 = (PropsSI('H','T',T_ex0,'P',P_ex0*1000,gas)/1000+y*h_l(T_ex0,P_ex0))/(1+y)
                
                P_ex1 = P_ex0 #pre exhaust HX pressure
                P_ex2 = P_ex0 #post const vol expansion pressure
                
                "Reverse solve for exhaust heat transfer"                
                cp_ex = (PropsSI('C','T',T_ex0,'P',P_ex1*1000,gas)/1000+y*c_l(T_ex0))/(1+y) #
                C_dot_ex = M_dot[i]*cp_ex
                NTU_ex = UA_ex/C_dot_ex
                epsilon_ex = 1.0 - math.exp(-NTU_ex)
                
                T_ex1 = Qex_reverse(gas,y,epsilon_ex,C_dot_ex,M_dot[i],Tw[j],h_ex0,P_ex0,T_ex0)
                h_ex1 = (PropsSI('H','T',T_ex1,'P',P_ex1*1000,gas)/1000+y*h_l(T_ex1,P_ex1))/(1+y)
                Q_dot_ex = M_dot[i]*(h_ex - h_ex1)
                
                "(3) Supply Cooling"
                cp_suc = (PropsSI('C','T',T_suc1,'P',P_suc1*1000,gas)/1000+y*c_l(T_suc1))/(1+y)
                C_dot_suc = (M_dot[i])*cp_suc
                NTU_suc = UA_suc/C_dot_suc
                epsilon_suc = 1- math.exp(-NTU_suc)
                Q_dot_suc = epsilon_suc*C_dot_suc*(T_suc1 - Tw[j])
                
                h_suc2 = h_suc1 - Q_dot_suc/M_dot[i]
                P_suc2 = P_suc1

                T_suc2 = T_mix_h(gas,h_suc2,P_suc2,y,T_suc1)

                "Thermal Envelope Energy Balance"
                Ebal =  Q_dot_ex + Q_dot_amb - Q_dot_suc - W_loss
                
                if j > 0:
                    g = np.append(g, Ebal)              
                    Tw = np.append(Tw, Tw[j] - (g[j]*(Tw[j]-Tw[j-1]))/(g[j]-g[j-1]))            
                else:                    
                    g = [Ebal]
                
                if j>12:
                    raise Exception("Twall not converging after 10x")
                    break
                     
                       
            "update mass fraction due to leakage"
            y_star = x_leak*y
            m_Leak = leakage(gas,T_suc2,P_suc2,P_ex2,y_star,A_leak)
            
            if m_Leak > M_dot[i]:   
                m_Leak = 0
            mg_leak = m_Leak/(1+y_star)
            ml_leak = mg_leak*y_star
            mg_int = m_g - mg_leak
            ml_int = m_l - ml_leak
            y_int = ml_int/mg_int
            m_int = mg_int + ml_int
            
            v_suc2 = (vol_g(gas,T_suc2,P_suc2) +y_int*vol_l(T_suc2))/(1+y_int)
            s_suc2 = (PropsSI('S','T',T_suc2,'P',P_suc2*1000,gas)/1000+y_int*s_l(T_suc2))/(1+y_int)
            
            "isentropic step"
            v_ex3 = (V_dot_suc_exp*V_ratio)/(m_int)
            [T_ex3,P_ex3] = IsentropicExpansion(gas,s_suc2,T_suc2,P_suc2,v_ex3,V_ratio,y_int,P_ex)
            h_ex3 = (PropsSI('H','T',T_ex3,'P',P_ex3*1000,gas)/1000+y_int*h_l(T_ex3,P_ex3))/(1+y_int)
            w_exp_int_s = h_suc2 - h_ex3         
            
            "const volume step"
            w_exp_v = v_ex3*(P_ex3 - P_ex2)
            h_ex2 = h_ex3 - w_exp_v
            T_ex2 = T_mix_h(gas,h_ex2,P_ex2,y,T_ex3)
            
            "Work guess Balance"
            Wbal = m_int*(w_exp_int_s+w_exp_v) - W_loss - W_exp[k]

            if k > 0:
                L = np.append(L, Wbal)              
                W_exp = np.append(W_exp, W_exp[k] - (L[k]*(W_exp[k]-W_exp[k-1]))/(L[k]-L[k-1]))            
            else:                    
                L = [Wbal]
            
            if k>30:
                
                raise Exception("Wbal not converging after 30x")
                break

        M_dot_int_calc = V_dot_suc_exp/v_suc2 #determined from machine volume and speed
        
        "Mass Balance"
        Mbal = m_int - M_dot_int_calc

        if i > 0:
            f = np.append(f, Mbal)              
            M_dot = np.append(M_dot, M_dot[i] - (f[i]*(M_dot[i]-M_dot[i-1]))/(f[i]-f[i-1]))            
        else:                    
            f = [Mbal]

        if i>50:
            raise Exception("M_dot not converging after 50x")
    
    "model results"
    W_exp = (m_int*(w_exp_int_s + w_exp_v) - W_loss)*1000 #kW to W
    W_s = m_int*(w_exp_int_s)*1000
    W_v = m_int*(w_exp_v)*1000
    
    h_su_exp = PropsSI('H','P',P_suc*1000,'T',T_suc,gas)
    s_suc = PropsSI('S','P',P_suc*1000,'T',T_suc,gas)
    h_ex_is = PropsSI('H','S',s_suc,'P',P_ex*1000,gas)
    
    M_dot_ref = M_dot[i]/(1+y)
    
    #overall isentropic efficiency
    eta_is = (W_exp)/(M_dot_ref*(h_su_exp - h_ex_is))


    print "Wexp [W]: ",W_s+W_v - W_loss*1000       
    print "mLeak [kg/s]: ",m_Leak
    print "Envelope Bal [W]: ",M_dot[i]*(h_ex - h_suc) + W_exp/1000 + Q_dot_amb
    print "Twall [K]: ",Tw[j]    
    print "oil mass fraction [-]: ",y
    print "Internal pressure ratio Psuc2/P_ex3 [-]: ",P_suc2/P_ex3
    print 'Overall isentropic efficiency -]:', eta_is    
    
    return [W_exp,M_dot_ref,T_ex,eta_is]




def SemiEmp_run():
    
    """
    Example - oil flooded scroll expander
    
    Reference: 
    """

    UA_suc_n = 27.550132090115376 
    UA_ex_n  = 5.03716888709217 
    UA_amb_n = 10.858289415986736 
    M_dot_n = 0.1645732025013935 
    V_ratio = 1.7322395586220594
    A_suc = 0.0003677244637767475 
    A_leak = 1.5601342948837973e-06
    x_leak = 0.0  
    T_loss = 2.8655765611888353
    W_loss_0 = 11.519426758760778 
    A_ex = 0.0014575503229525624 

    Ref = 'R134A'
    T_amb = 25 + 273.15
    V_suc_comp = 104.8e-6 #[m3]
    
    T_suc = 110 + 273.13 #[K]
    P_suc = 2000 #[kPa]
    P_ex = 1000 #[kPa]
    N_exp = 1800 #[rpm]
    y = 0.1 #[-]
    
    "Let's time the simulation"
    time0 = time.time()
    
    [W,M,T,eta] = Expander(Ref,T_amb,T_suc,P_suc,P_ex,y,N_exp,UA_suc_n,UA_ex_n,UA_amb_n,M_dot_n,V_ratio,A_suc,A_leak,x_leak,V_suc_comp,T_loss,W_loss_0,A_ex)

    time1 = time.time() 
    print 'Executed in', (time1 - time0), 'seconds'        

if __name__=='__main__':
    
    SemiEmp_run()






