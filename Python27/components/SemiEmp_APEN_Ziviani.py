from __future__ import division
# -*- coding: utf-8 -*-
"""
Supplemental code for paper:

APEN-D-16-03721
"Characterizing the performance of a single-screw expander in a small-scale organic Rankine cycle for waste heat recovery"

Semi-Emperical Model of single-screw expander
Baseline: no flooding

@author:
Davide Ziviani
davide.ziviani@ugent.be/dziviani@purdue.edu

"""

#Imports
from CoolProp.CoolProp import PropsSI
import numpy as np
import math
import time
import warnings

warnings.simplefilter("ignore",RuntimeWarning)
from random import randint, random

#Numerical tolerance
numTol = 1e-5 



#=============================================================================
"Lubricant oil Properties"

def c_l(T):
    "specific heat [kJ/kg-K] given T in K"
    c = (5.186*T + 337.116)/1000 
    return c
    
def u_l(T):
    "internal energy [kJ/kg] given T in K"
    u = (5.186/2*pow(T,2)+337.116*T)/1000  
    return u

def rho_l(T):
    "density [kg/m^3] given T in K"
    rho = -0.667*T+1050.86 
    return rho

def s_l(T):
    "specific entropy [kJ/kg-K] given T in K"  
    s = (5.186*(T-298) + 337.116*np.log(T/298.0))/1000  
    return s

def h_l(T,P):
    "the specific enthalpy [kJ/kg-k]"
    h = u_l(T) + P/rho_l(T)
    return h

def mu_l(T):
    "returns the viscosity given temp in K, [Pa-s]"
    mu = -0.0001235*T + 0.04808 #Zerol60
    return mu


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
    "returns the specific volume so inverse of density doesn't need to be manually entered"
    v = 1.0/(PropsSI('D','T',T,'P',P*1000+100,gas))
    return v
        
#==========================================================================    

"Mixture Function Calls"

def T_mix_h(gas,h,P,y,Tguess):
    "find the mixture temperature given its specific enthalpy"
    T=[Tguess+5, Tguess-5]
    h_check = (PropsSI('H','T',T[0],'P',P*1000,gas)/1000+ y*h_l(T[0],P))/(1+y)
    f =[abs(h_check - h)]
    
    i=0 #array index
    while abs(f[i])> 1e-5: #numTol:
        i=i+1 #update index
        if T[i]< Tguess/2: #ride herd
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
    
    i=0 #array index
    while abs(f[i])> numTol:
        i=i+1 #update index
        if T[i]< Tguess/2: #ride herd
            T[i] = Tguess/2 + random()*(Tguess - Tguess/2)
            
        s_check = (PropsSI('S','T',T[i],'P',P*1000,gas)/1000+ y*s_l(T[i]))/(1+y)
        f = np.append(f,abs(s_check - s))
        T = np.append(T , T[i]-(f[i]*(T[i]-T[i-1]))/(f[i]-f[i-1])) #secant method
         
        if i>100:
            raise Exception("T_mix_s not converging after 100x")
    
    return T[i]
    
def gamma_mix(gas,T,P,y):
    Cp = (PropsSI('C','T',T,'P',P*1000+100,gas)/1000+y*c_l(T))/(1+y)
    Cv = (PropsSI('O','T',T,'P',P*1000+100,gas)/1000+y*c_l(T))/(1+y)
    gamma = Cp/Cv
    return gamma


#======================================================================   

"Selected Modelling steps in functions for overall clarity"

def SuctionNozzle(gas,Area,mg,ml,P1,T1,h1,s1):
    y = ml/mg

    gamma = gamma_mix(gas,T1,P1,y)
    P_Crit = P1*pow((2/(gamma+1)),(gamma/(gamma-1)))  #serve as bound on iterations
    
    
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
    
#    print "Fhigh: ",FHigh
#    print "FLow: ",FLow
    
    NumSections = 0
    FMid = 10
    
    while abs(FMid) > numTol:
       
        NumSections = NumSections+1
        PMid = (PHigh+PLow)/2.0
        FMid = signChange*Suction_helper(gas,h1,s1,T1,PMid,Area,mg,ml,y) 
        
        if FMid>0:
            PHigh = PMid
        else:
            PLow = PMid
        
        
        if NumSections>100:
#            print "Phigh: ",PHigh
#            print "Plow: ",PLow
#            print "PMid: ",PMid
#            print "Fmid: ",FMid
            raise Exception("Nozzle not converging after 100x")
    
    P2 = PMid    
    
    "recover static enthalpy in isobaric diffuser"
    T2 = PropsSI('T','H',h1*1000,'P',P2*1000,gas)#T_mix_h(gas,h1,P2,y,T1)
#    print "Pnoz numSecs: ",NumSections
    
    return [T2,P2]
    
    
def Suction_helper(gas,h1,s1,T1,P2,Area,mg,ml,y):
    
    T2 = PropsSI('T','S',s1*1000,'P',P2*1000,gas)#T_mix_s(gas,s1,P2,y,T1) #isentropic nozzle
    h2 = PropsSI('H','S',s1*1000,'P',P2*1000,gas)/1000
    v2 = 1.0/(PropsSI('D','H',h2*1000,'P',P2*1000,gas))#(vol_g(gas,T2,P2) +y*vol_l(T2))/(1+y)
    vel2 = (mg+ml)/(Area/v2)
    KE2 = 0.5*vel2**2.0/1000.0
    Ebal = h1 - (h2 + KE2)        

    return Ebal    


def IsentropicExpansion(gas,s1,T1,P1,v2,v_ratio,y,P_ex):
    "Perform isentropic expansion along the volume ratio of the machine"
    
    "guess outlet pressure"
    P= [P1/v_ratio,0.95*P1/v_ratio]
    
    i=0 #secant array index for pressure
    f = [10.0000] #initialize convergence function (unlikely to be a naturally occuring soln)
            
    while abs(f[i]) > numTol:        
        if f[i] != 10.0000: #only increase index after 1st iteration
            i=i+1
        
        "ride herd on Pressure"
        if P[i] > P1:
            P[i] = P_ex - random()*(P1-P_ex)
#        if P[i] < P1/(v_ratio*2.0):
#            P[i] = P_ex + random()*(P1-P_ex)
        
        "guess outlet temperature"
        if y > 1:
            T = [T1+2,T1+5]
        else:
            T= [T1+20,T1+40]
        
        j = 0 #secant array index for temperature
        g = [10.0000] #initialize convergence function (unlikely to be a naturally occuring soln)
        
        
        T[j] = PropsSI('T','S',s1*1000,'P',P[i]*1000+100,gas)
        
#        while abs(g[j]) > numTol:        
#            if g[j] != 10.0000: #only increase index after 1st iteration
#                j=j+1
#                
#            "Entropy balance"
#            s_out = (Props('S','T',T[j],'P',P[i],gas)+y*s_l(T[j]))/(1+y)
#            Sbal = s1 - s_out
#            
#            if j > 0:
#                g = np.append(g, Sbal)              
#                T = np.append(T, T[j] - (g[j]*(T[j]-T[j-1]))/(g[j]-g[j-1]))            
#            else:                    
#                g = [Sbal] 
#            if j>20:
#                raise Exception("IsenExp-T not converging after 20x")

        "Mass (volume) Balance"
        v_out = 1.0/(PropsSI('D','S',s1*1000,'P',P[i]*1000,gas))#(vol_g(gas,T[j],P[i]) +y*vol_l(T[j]))/(1+y)
        Mbal = v2 - v_out
        
        if i > 0:
            f = np.append(f, Mbal)              
            P = np.append(P, P[i] - (f[i]*(P[i]-P[i-1]))/(f[i]-f[i-1]))            
        else:                    
            f = [Mbal]         
        if i>50:  
            print "Mbal: isenExp: ",Mbal                                 
            raise Exception("IsenExp-P not converging after 50x")
            break
            
    P = P[i]
    T = T[j]

    return [T,P]
    
    
    
    
def leakage(gas,T_pleleak,P_preleak,P_postleak,y_star,A_leak):

    gamma_leak = gamma_mix(gas,T_pleleak,P_preleak,y_star)
    P_leak_crit = P_preleak*pow((2.0/(gamma_leak+1)),gamma_leak/(gamma_leak-1))
    P_leak_thr = max(P_postleak,P_leak_crit) #choke flow in leakage
    
    s_preleak = (PropsSI('S','T',T_pleleak,'P',P_preleak*1000+100,gas)/1000+y_star*s_l(T_pleleak))/(1+y_star)
    h_preleak = (PropsSI('H','T',T_pleleak,'P',P_preleak*1000+100,gas)/1000+y_star*h_l(T_pleleak,P_preleak))/(1+y_star)
    
    T_leak_thr = PropsSI('T','S',s_preleak*1000,'P',P_leak_thr*1000,gas)#T_mix_s(gas,s_ex2_preleak,P_leak_thr,y_star,T_ex2) #isentropic nozzle leaking
    v_leak_thr = 1.0/PropsSI('D','S',s_preleak*1000,'P',P_leak_thr*1000,gas)#(vol_g(gas,T_leak_thr,P_leak_thr) +y_star*vol_l(T_leak_thr))/(1+y_star)
    h_leak_thr = PropsSI('H','S',s_preleak*1000,'P',P_leak_thr*1000,gas)/1000#(PropsSI('H','T',T_leak_thr,'P',P_leak_thr,gas)+y_star*h_l(T_leak_thr,P_leak_thr))/(1+y_star)
    
    vel_leak_thr = math.sqrt(2.0*(h_preleak - h_leak_thr)*1000)  #kJ/kg --> J/kg
    
    m_leak = 1.0/v_leak_thr*vel_leak_thr*A_leak

    return m_leak




def ExhaustNozzle(gas,Area,mg,ml,P2,T2,P_suc):

    "Nozzle inlet pressure"
    y = ml/mg
    h2 = (PropsSI('H','T',T2,'P',P2*1000+100,gas)/1000+y*h_l(T2,P2))/(1+y)
    s2 = PropsSI('S','T',T2,'P',P2*1000+100,gas)/1000
    v2 = (vol_g(gas,T2,P2) +y*vol_l(T2))/(1+y)
    vel2 = (mg+ml)/(Area/v2)
    KE = 0.5*vel2**2/1000
    h1 = h2 + KE
    s1=s2
    
    "guess nozzle inlet pressure"
    P1 = [P2+10,P2+50]
    
    j=0
    g=[10.0000]
    
    while abs(g[j])> numTol:
        if g[j] != 10.0000: #only increase index after 1st iteration
            j=j+1   
        
        h1_check = PropsSI('H','S',s1*1000,'P',P1[j]*1000,gas)/1000
        Ebal = h1 - h1_check
        
        if j > 0:
            g = np.append(g, Ebal)              
            P1 = np.append(P1, P1[j] - (g[j]*(P1[j]-P1[j-1]))/(g[j]-g[j-1]))            
        else:                    
            g = [Ebal]         
        if j>30:                       
            raise Exception("Exhaust Nozzle P not converging after 30x") 
    
    T1 = PropsSI('T','P',P1[j]*1000,'S',s1*1000,gas) #isentropic nozzle
        
        
    return [T1,P1[j]]




def Qex_reverse(gas,y,epsilon_ex,C_dot_ex,M_dot,Tw,h_ex,P_ex,T_ex):
    
    "Qex = epsilon*Cdot*(Tex1 - Tw) = mDot*(hex1-hex)"  #compressor
    "h(Tex1,Pex) - epsilon*Cdot/mDot*Tex1 = -epsilon*Cdot/mDot*Tw + hex" #compressor
    
    "Qex = epsilon*Cdot*(Tw - Tex1) = mDot*(hex-hex1)"  #expander
    "h(Tex1,Pex) - epsilon*Cdot/mDot*Tex1 = -epsilon*Cdot/mDot*Tw + hex" #expander
    
    balance = -epsilon_ex*C_dot_ex/M_dot*Tw + h_ex
    
#    T_ex1 = [T_ex+3, T_ex+5]
    h_ex1 = [h_ex, h_ex*0.95]
    
    i=0
    f = [10.0000]
    
    while abs(f[i]) > numTol:
        if f[i] != 10.0000: #only increase index after 1st iteration
            i=i+1
        
#        bal = Props('H','T',T_ex1[i],'P',P_ex,gas) - epsilon_ex*C_dot_ex/M_dot*T_ex1[i]
        T_ex1 = PropsSI('T','H',h_ex1[i]*1000,'P',P_ex*1000,gas)
        bal = h_ex1[i] - epsilon_ex*C_dot_ex/M_dot*T_ex1
        
        Ebal = bal - balance
        
        if i > 0:
            f = np.append(f, Ebal)              
            h_ex1 = np.append(h_ex1, h_ex1[i] - (f[i]*(h_ex1[i]-h_ex1[i-1]))/(f[i]-f[i-1]))            
        else:                    
            f = [Ebal]         
        if i>30: 
            print Ebal                      
            raise Exception("Qex-reverse converging after 30x")        
        
        
    
    
    return h_ex1[i]

# Electro-mechanical losses
def eta_el_inverter(rpm, W_el_exp):
    
    rpm_nom = 2930
    W_el_nom = 11
    ln_N = math.log(rpm/rpm_nom)
    ln_W = math.log(W_el_exp/W_el_nom)
    
    eta_el = 9.55726922e-1 + 2.60983262e-2*ln_N + 2.42349302e-2*ln_N**2 + 1.21191602e-2*ln_N**3 + 4.94828374e-2*ln_W + 3.34143316e-2*ln_W**2 + 2.27446360e-2*ln_W**3
    
    return eta_el


def eta_sh(tau,rpm):
    """
    Shaft Efficiency 11kW ORC Kortrijk
    """
    rpm_nom = 2930                #rpm nominal rotational speed 
    omega_nom = (2*math.pi*rpm_nom)/60                   #Frequency, rad/sec (60Hz)    
    W_dot_nom = 11                #Nominal power output kW
    T_nom = W_dot_nom*1e3 / omega_nom   # Nm  Nominal torque 35Nm
    
    eta = 8.93747915e-1+3.23048796e-2*np.log(rpm/rpm_nom)-1.91761519e-2* (np.log(rpm/rpm_nom))**2 +  1.52204756e-2*(np.log(rpm/rpm_nom))**3 + 7.32867448e-3*np.log(tau/T_nom) - 3.17061820e-2*(np.log(tau/T_nom))**2 + 2.16415080e-2*(np.log(tau/T_nom))**3
    + 1.63125253e-2*np.log(rpm/rpm_nom)*np.log(tau/T_nom) + 4.37556935e-3*np.log(rpm/rpm_nom)*(np.log(tau/T_nom))**2 - 4.11952262e-2*(np.log(tau/T_nom))*(np.log(rpm/rpm_nom))**2
    - 1.62681324e-2*(np.log(tau/T_nom))**2*(np.log(rpm/rpm_nom))**2

    return eta

#==============================================================
# Main Solver 
def Expander(gas,T_amb,T_suc,P_suc,P_ex,y,N_exp,UA_suc_n,UA_ex_n,UA_amb_n,M_dot_n,V_ratio,A_suc,A_leak,x_leak,V_suc_comp,T_loss,A_ex):
    "Solution method for semi-emperical scroll expander"
    start_time = time.time()
    
    
    N_exp = N_exp*1.0 #convert int to double

    "Friction Losses"
    W_loss =  2.0*math.pi*(N_exp/60)*T_loss/1000 #W --> kW 

    "generate a starting guess for massflow iterations"
    V_suc_exp =V_suc_comp/V_ratio
    
    #Mass flow rate guess value
    mguess = 2*6*V_suc_exp*(N_exp/60)/((vol_g(gas,T_suc,P_suc) +y*vol_l(T_suc))/(1+y))
    mLeak_guess = leakage(gas,T_suc,P_suc,P_ex,y,A_leak)

    mguess = mguess + mLeak_guess
    

    M_dot =  [mguess,0.9*mguess]
#    M_dot = [0.0381232839733]
    i=0 #secant array index for massflow rate
    f = [10.0000] #initialize convergence function (unlikely to be a naturally occuring soln)
    
#    iii=0
    while abs(f[i]) > numTol:        
#    while iii < 1:
#        iii=12
        if f[i] != 10.0000: #only increase index after 1st iteration
            i=i+1
        
        "ride herd on massflow"
        if M_dot[i] < 0:
#        if M_dot[i] < mguess/2.0:
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
                [T_suc1,P_suc1] = SuctionNozzle(gas,A_suc,m_g,m_l,P_suc,T_suc,h_suc,s_suc)
                h_suc1 = (PropsSI('H','T',T_suc1,'P',P_suc1*1000,gas)/1000+y*h_l(T_suc1,P_suc1))/(1+y)
                s_suc1 = (PropsSI('S','H',h_suc1*1000,'P',P_suc1*1000,gas)/1000+y*s_l(T_suc1))/(1+y)
                nozzlePass = True
            except:
                M_dot[i] = M_dot[i]*(0.95-nozReduce*0.02) #gradually decrease till a sufficient mass flow range is found
        
    
        "Actual UA Values"
        UA_suc = UA_suc_n*pow((M_dot[i]/M_dot_n),0.8)/1000  #W --> kW
        UA_ex = UA_ex_n*pow((M_dot[i]/M_dot_n),0.8)/1000  #W --> kW
        UA_amb = UA_amb_n/1000 #W --> kW
        
        
        "Generate starting guesses for Work"
        try:
            "guess isentropic step"
            V_dot_suc_exp = 2*6*V_suc_exp*(N_exp/60.0)
            v_ex3 = (V_dot_suc_exp*V_ratio)/(M_dot[i]-mLeak_guess)
            [T_ex3,P_ex3] = IsentropicExpansion(gas,s_suc1,T_suc1,P_suc1,v_ex3,V_ratio,y,P_ex)
            h_ex3 = (PropsSI('H','T',T_ex3,'P',P_ex3*1000,gas)/1000+y*h_l(T_ex3,P_ex3))/(1+y)
            w_exp_int_s = h_suc1 - h_ex3        
        
            "guess const volume step"
            P_ex2 = P_ex 
            w_exp_v = v_ex3*(P_ex3 - P_ex2)
        
            Wexp_guess = (M_dot[i]-mLeak_guess)*(w_exp_int_s+w_exp_v)-W_loss
        
        except:
            Wexp_guess = 1.0
           
        
        W_exp = [Wexp_guess,0.9*Wexp_guess]
#        W_exp = [1050.03133631/1000]
        
        
        k = 0 #index for work guesses
        L = [10.0000]
        
#        kk=0
        while abs(L[k]) > numTol:
#        while kk < 1:
#            kk = 12
        
            if L[k] != 10.0000:
                k=k+1
            
            "ride herd on work"
            if W_exp[k] < 0:
                W_exp[k] = 0.25*Wexp_guess + random()*(0.5*Wexp_guess)
#            print "Wcomp guess: ", W_comp[k]   

    
            "Guess expander body wall temperature"
            try:
                TwGuess = W_loss/UA_amb + T_amb #assumes Q_suc ~ Q_ex
            except ZeroDivisionError:
                TwGuess = T_amb
            Tw = [TwGuess,TwGuess-5] #2nd guess slightly lower since typically Q_ex > Q_suc (for expander)
#            Tw = [381.907097956]
            
            j=0 #secant array index for massflow rate
            g = [10.0000] #initialize convergence function (unlikely to be a naturally occuring soln)
            
#            jjj=0
            while abs(g[j]) > numTol: 
#            while jjj < 1:
#                jjj=12
    #            print "g[j]: ",g[j]
                if g[j] != 10.0000: #only increase index after 1st iteration
                    j=j+1
                    
                "ride herd on wall temp"
                if Tw[j] < 0:
                     Tw[j] = T_suc + random()*(T_amb - T_suc)
#                if Tw[j] > 450:
#                     Tw[j] = T_suc + random()*(T_amb - T_suc)
                    
            
                "Ambient Heat Loss"
                Q_dot_amb = UA_amb*(Tw[j] - T_amb)
                
                "Envelope Energy Balance"
                #hex*m+Qamb+Wexp = hsuc*m
                h_ex = (h_suc*M_dot[i] - W_exp[k] - Q_dot_amb) / M_dot[i]
                
                
                
                T_ex = PropsSI('T','H',h_ex*1000,'P',P_ex*1000,gas)#T_mix_h(gas,h_ex,P_ex,y,T_ex3)
                
                [T_ex0,P_ex0] = ExhaustNozzle(gas,A_ex,m_g,m_l,P_ex,T_ex,P_suc)
#                [T_ex0,P_ex0] = [T_ex,P_ex]
                h_ex0 = (PropsSI('H','T',T_ex0,'P',P_ex0*1000,gas)/1000+y*h_l(T_ex0,P_ex0))/(1+y)
                
                P_ex1 = P_ex0 #pre exhaust HX pressure
                P_ex2 = P_ex0 #post const vol expansion pressure
                
#                P_ex1 = P_ex #assume no outet pressure drop                
                
                "Reverse solve for exhaust heat transfer"                
                cp_ex = (PropsSI('C','T',T_ex0,'P',P_ex1*1000,gas)/1000+y*c_l(T_ex0))/(1+y) #assume cp doesnt change much b/t Tex0 &  Tex1
                C_dot_ex = M_dot[i]*cp_ex
                NTU_ex = UA_ex/C_dot_ex
                epsilon_ex = 1.0 - math.exp(-NTU_ex)
                
                h_ex1 = Qex_reverse(gas,y,epsilon_ex,C_dot_ex,M_dot[i],Tw[j],h_ex0,P_ex0,T_ex0)
                T_ex1 = PropsSI('T','H',h_ex1*1000,'P',P_ex1*1000,gas)
                Q_dot_ex = M_dot[i]*(h_ex - h_ex1)
                                
                
                "(3) Supply Cooling"
                cp_suc = (PropsSI('C','T',T_suc1,'P',P_suc1*1000,gas)/1000+y*c_l(T_suc1))/(1+y)
                C_dot_suc = (M_dot[i])*cp_suc
                NTU_suc = UA_suc/C_dot_suc
                epsilon_suc = 1- math.exp(-NTU_suc)
                Q_dot_suc = epsilon_suc*C_dot_suc*(T_suc1 - Tw[j])
                
                h_suc2 = h_suc1 - Q_dot_suc/M_dot[i]
                P_suc2 = P_suc1
    
#                print "twall: ",Tw[j]
                T_suc2 = PropsSI('T','H',h_suc2*1000,'P',P_suc2*1000,gas) #T_mix_h(gas,h_suc2,P_suc2,y,T_suc1)
#                s_suc2 = (PropsSI('S','T',T_suc2,'P',P_suc2,gas)+y*s_l(T_suc2))/(1+y)
#                v_suc2 = (vol_g(gas,T_suc2,P_suc2) +y*vol_l(T_suc2))/(1+y)
                
                "Thermal Envelope Energy Balance"
                Ebal =  Q_dot_ex + Q_dot_amb - Q_dot_suc - W_loss
                
                if j > 0:
                    g = np.append(g, Ebal)              
                    Tw = np.append(Tw, Tw[j] - (g[j]*(Tw[j]-Tw[j-1]))/(g[j]-g[j-1]))            
                else:                    
                    g = [Ebal]
                
                if j>10:
                    raise Exception("Twall not converging after 10x")
                    break
                
                
#            "Work Balance"
#            T_ex2 = T_ex1
#            P_ex2 = P_ex1
#            h_ex2 = h_ex1#(Props('H','T',T_ex2,'P',P_ex2,gas)+y_int*h_l(T_ex2,P_ex2))/(1+y_int)
            
#            "leakage mixing at outlet"
#            h_leak = (Props('H','T',T_ex2,'P',P_ex2,gas)+y_star*h_l(T_ex2,P_ex2))/(1+y_star)
#            m_int = M_dot[i] + m_Leak
#            h_suc3 = (M_dot[i]*h_suc2 + m_Leak*h_leak) / m_int                        
                       
            "update mass fraction due to leakage"
            y_star = x_leak*y
            m_Leak = leakage(gas,T_suc2,P_suc2,P_ex2,y_star,A_leak)
            if m_Leak > M_dot[i]:   #false value and guess again later
                m_Leak = 0
            mg_leak = m_Leak/(1+y_star)
            ml_leak = mg_leak*y_star
            mg_int = m_g - mg_leak
            ml_int = m_l - ml_leak
            y_int = ml_int/mg_int
            m_int = mg_int + ml_int
            
#            P_suc3 = P_suc2 #isobaric mixing
#            T_suc3 = PropsSI('T','H',h_suc3,'P',P_suc3,gas)#T_mix_h(gas,h_suc3,P_suc3,y_int,T_suc2)
            v_suc2 = (vol_g(gas,T_suc2,P_suc2) +y_int*vol_l(T_suc2))/(1+y_int)
            s_suc2 = (PropsSI('S','T',T_suc2,'P',P_suc2*1000+100,gas)/1000+y_int*s_l(T_suc2))/(1+y_int)
            
            
            "isentropic step"
            v_ex3 = (V_dot_suc_exp*V_ratio)/(m_int)
            [T_ex3,P_ex3] = IsentropicExpansion(gas,s_suc2,T_suc2,P_suc2,v_ex3,V_ratio,y_int,P_ex)
            h_ex3 = (PropsSI('H','T',T_ex3,'P',P_ex3*1000+100,gas)/1000+y_int*h_l(T_ex3,P_ex3))/(1+y_int)
            w_exp_int_s = h_suc2 - h_ex3         
            
            "const volume step"
            w_exp_v = v_ex3*(P_ex3 - P_ex2)
            h_ex2 = h_ex3 - w_exp_v
            T_ex2 = T_mix_h(gas,h_ex2,P_ex2,y,T_ex3)
            
            "Work guess Balance"
            Wbal = m_int*(w_exp_int_s+w_exp_v) - W_loss - W_exp[k]
#            print "Wbal: ",Wbal
#            print "Wcomp[k]: ",W_comp[k]
            
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
#        print "m_int: ",m_int
#        print "mDot calc: ",M_dot_int_calc
        if i > 0:
            f = np.append(f, Mbal)              
            M_dot = np.append(M_dot, M_dot[i] - (f[i]*(M_dot[i]-M_dot[i-1]))/(f[i]-f[i-1]))            
        else:                    
            f = [Mbal] 
        if i>50:
            raise Exception("M_dot not converging after 50x")
    

    
    "model results"
    W_exp_sh = (m_int*(w_exp_int_s + w_exp_v) - W_loss)*1000 #kW to W
#    W_exp = W_exp[0]*1000

    tau_sh = W_exp_sh/(2*math.pi*N_exp/60)  #Nm
    W_exp_el= W_exp_sh*eta_sh(tau_sh,N_exp) #W
    W_exp = W_exp_el*eta_el_inverter(N_exp, W_exp_el/1000) #W

    W_s = m_int*(w_exp_int_s)*1000
    W_v = m_int*(w_exp_v)*1000
    
#    print "Wexp_again: ",W_s+W_v - W_loss*1000
#    print "mLeak: ",m_Leak
#    print "Envelope Bal: ",M_dot[i]*(h_ex - h_suc) + W_exp/1000 + Q_dot_amb

    h_su_exp = PropsSI('H','P',P_suc*1000,'T',T_suc,gas)
    s_suc = PropsSI('S','P',P_suc*1000,'T',T_suc,gas)
    h_ex_is = PropsSI('H','S',s_suc,'P',P_ex*1000,gas)
    
    M_dot_ref = M_dot[i]/(1+y)
    
    #Overall isentropic efficiency
    eta_is = (W_exp)/(M_dot_ref*(h_su_exp - h_ex_is))
    #print eta_is


    return [W_exp,M_dot[i],T_ex,eta_is] 



#======================================================================    

"Create a function to run the model"

def Single_Run(Args):
    
    UA_suc_n = Args[0] 
    UA_ex_n = Args[1]
    UA_amb = Args[2]
    M_dot_n = Args[3]
    V_ratio = Args[4]
    A_suc = Args[5]
    A_leak = Args[6]
    T_loss = Args[7]
    x_leak = 0 #baseline model assumes only gas leakage
    A_ex = Args[8]
    
    #Boundary conditions
    gas = 'R245FA'
    V_suc_comp = 58.47e-6 #[m3]
    T_amb = 25 + 273.15 #[K]
    T_suc = 125 + 273.15 #[K]
    P_suc = 1000 #[kPa]
    P_ex = 150 #[kPa]
    N_exp = 3000  #[rpm]
    y = 0 #[-]
    
    
    "Run expander model"
    time0 = time.time()
    
    [W,M,Tex,eta] = Expander(gas,T_amb,T_suc,P_suc,P_ex,y,N_exp,UA_suc_n,UA_ex_n,UA_amb,M_dot_n,V_ratio,A_suc,A_leak,x_leak,V_suc_comp,T_loss,A_ex)

    time1 = time.time() 
    print 'Executed in', (time1 - time0), 'seconds'        
  
    
    print "Discharge Temperature [C]:",Tex-273.15
    print "Mass flow rate [kg/s]:", M
    print "Electric power output [W]:", W
    print "Overall isentropic efficiency [-]:", eta





if __name__ == '__main__':
    
    #Calibrated coefficients R245fa and single-screw expander
    Args = [46.35324727168612, 0.08090805926654343, 1.2020204580202711, 2.5304350109368547, 5.978417696222445, 0.00010903621290184344, 1.8796046985075052e-05, 3.227410827057272, 0.0005814361248644597]
    Single_Run(Args)

