from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import Correlations
from math import pi
import numpy as np


def ORC_Preconditioner(Cycle,Inputs):
    
    #Assume the heat exchangers are highly effective
        
    def OBJECTIVE(x):
        T_evap=x[0]
        T_cond=x[1]
             
        Ref = Inputs['Ref']
        
        #Pump
        p_cond=PropsSI('P','T',T_cond+273.15,'Q',1,Ref)/1000   #Pa
        p_evap=PropsSI('P','T',T_evap+273.15,'Q',1,Ref)/1000   #Pa

        Cycle.Pump.pin_r=p_cond
        Cycle.Pump.pout_r=p_evap
        Cycle.Pump.Tin=T_cond-Inputs['DT_sc']
        Cycle.Pump.Calculate()
        W_dot_pump=Cycle.Pump.W_dot
        mdot_pump = Cycle.Pump.m_dot    
        h_evap_su_c = Cycle.Pump.hout*1000
        #Use maximum heat transfer to get a guess for the evaporator capacity
        DELTAh=(PropsSI('H','T',Inputs['Tin_h']+273.15,'P',p_evap*1000,Ref)-PropsSI('H','T',Cycle.Pump.Tout,'P',p_evap*1000,Ref))   #J/kg
        epsilon_ev=0.99
        Qevap=epsilon_ev*Cycle.Pump.m_dot*DELTAh
        
        h_evap_ex_c = h_evap_su_c + Qevap/mdot_pump
        
        #use the expander model to find its power and mass flow rate
        params={
                'P_su':p_evap,
                'P_ex':p_cond,
                'h_su':h_evap_ex_c,
                'inputs':'Psu_Pex_h'
                }
        Cycle.Expander.Update(**params)
        Cycle.Expander.Calculate()
        W_dot_exp=Cycle.Expander.W_dot
        #print 'W_dot_exp',W_dot_exp
        
        #Use fixed effectiveness to get a guess for the condenser capacity
        epsilon_cd=0.96
        Qcond=epsilon_cd*Inputs['mdot_c']*4187*(T_cond-Inputs['Tin_c'])

        resids=[W_dot_pump+Qevap-W_dot_exp-Qcond,Cycle.Pump.m_dot-Cycle.Expander.m_dot]
        print 'first guess cycle resid:',resids
        return resids
        
    Tevap_init=Inputs['Tin_h']-35
    Tcond_init=Inputs['Tin_c']+10
    x_0=[Tevap_init,Tcond_init]
    x=fsolve(OBJECTIVE,x_0,xtol=1e-5)
    print 'Exiting Preconditioner....'
    print 'Tev,Tcd:',x
    return x
      


def ORCRegen_Preconditioner(Cycle,Inputs):
    
    def OBJECTIVE_REGEN(x):    
        T_evap=x[0]
        T_cond=x[1]

        Ref = Inputs['Ref']
             
        #Pump
        p_cond=PropsSI('P','T',T_cond+273.15,'Q',1,Ref)/1000   #Pa
        p_evap=PropsSI('P','T',T_evap+273.15,'Q',1,Ref)/1000   #Pa
        Cycle.Pump.pin_r=p_cond
        Cycle.Pump.pout_r=p_evap
        Cycle.Pump.Tin=T_cond-Inputs['DT_sc']
        Cycle.Pump.Calculate()
        W_dot_pump=Cycle.Pump.W_dot
        mdot_pump = Cycle.Pump.m_dot    
        h_evap_su_c = Cycle.Pump.hout*1000

        #Use maximum heat transfer to get a guess for the evaporator capacity
        h_evap_ex_c_max = PropsSI('H','T',Inputs['Tin_h']+273.15,'P',p_evap,Ref)
        DELTAh_evap_max = h_evap_ex_c_max - h_evap_su_c
        Qevap = 0.99*mdot_pump*DELTAh_evap_max
        h_evap_ex_c = h_evap_su_c + Qevap/mdot_pump
        T_evap_ex_c = PropsSI('T','P',p_evap,'H',h_evap_ex_c,Ref)

        #Knowing T_cond and T_evap, the expander map can be used to get the flow rate and expander power
        params={
                'P_su':p_evap,
                'P_ex':p_cond,
                'h_su':h_evap_ex_c,
                'inputs':'Psu_Pex_h'
                }
        Cycle.Expander.Update(**params)
        Cycle.Expander.Calculate()
        W_dot_exp=Cycle.Expander.W_dot
        #print 'Wdot_exp:',W_dot_exp
        
        #Use fixed effectiveness to get a guess for the condenser capacity
        epsilon_cd=0.96
        Qcond=epsilon_cd*Inputs['mdot_c']*4187*(T_cond-Inputs['Tin_c'])

        resids=[W_dot_pump+Qevap-W_dot_exp-Qcond,Cycle.Pump.m_dot-Cycle.Expander.m_dot]
        print 'first guess cycle resid:',resids
        return resids
        
         
    Tevap_init=Inputs['Tin_h']-35
    Tcond_init=Inputs['Tin_c']+10
    x_0=[Tevap_init,Tcond_init]
    x=fsolve(OBJECTIVE_REGEN,x_0,xtol=1e-5)
    print 'Exiting Preconditioner....'
    print 'Tev,Tcd:',x
    
    return [x[0],x[1],Cycle.Expander.hout]
    
