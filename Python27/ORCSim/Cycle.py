from __future__ import division
import sys, copy


from Expander import ExpanderClass
from PHEX_ASME2015 import PHEHXClass
from LineSet import LineSetClass
from Pump import PumpClass
from Pump_Centrifugal_Kortrijk import PumpCentrifugalClass 
from LiquidReceiver import LiquidReceiverClass
from scipy.optimize import brentq, fsolve,newton ,anderson
from math import pi,fabs
from CoolProp.CoolProp import PropsSI           
from first_guess_cycle import ORC_Preconditioner,ORCRegen_Preconditioner
import itertools
from counter import my_counter
from ACHPTools import convert#,Write2CSV

import numpy as np                  
from numpy import *
import matplotlib
from matplotlib.pyplot import plot, show, figure, semilogy, xlim, ylim, title, xlabel, ylabel, legend 
import matplotlib.pyplot as plt
import pylab
from pylab import arange,pi,sin,cos,sqrt,linspace

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset



def F2K(T_F):
    """
    Convert temperature in Fahrenheit to Kelvin
    """                       
    return 5/9*(T_F+459.67)
        
class ORCClass():


    def __init__(self):
        """
        Load up the necessary sub-structures to be filled with
        the code that follows
        """
        self.Evaporator=PHEHXClass()
        self.Expander=ExpanderClass()
        self.Regenerator = PHEHXClass()
        self.Condenser=PHEHXClass()
        self.Pump=PumpCentrifugalClass()
        self.LiquidReceiver=LiquidReceiverClass()
        self.LineSetPumpEx=LineSetClass()
        self.LineSetPumpSu=LineSetClass()
        self.LineSetExpSu=LineSetClass()
        self.LineSetExpEx=LineSetClass()
        self.LineSetCondEx=LineSetClass()
    
    def OutputList(self):
        """
        Return a list of parameters for this component for further output
        It is a list of tuples, and each tuple is formed of items:
        [0] Description of value
        [1] Units of value
        [2] The value itself
        """
        try:
            return [
                    #Residuals
                    ('resid_subcool','K',self.resid_DT_sc),
                    ('resid_charge','kg',self.resid_charge),
                    ('resid_ebal','W',self.resid_EB),
                    ('resid_mdot','kg/s',self.resid_mf),
                    #Hot source
                    ('Th2','C',self.Evaporator.Tout_h-273.15),
                    ('Tc2','C',self.Condenser.Tout_c-273.15),
                    #Pump
                    ('T_pump_su','C',self.Pump.Tin),
                    ('T_pump_ex','C',self.Pump.Tout - 273.15),
                    ('P_pump_su','kPa',self.Pump.pin_r),
                    ('P_pump_ex','kPa',self.Pump.pout_r),
                    ('mdot_pump','kg/s',self.Pump.m_dot),
                    ('Power_pump','W',self.Pump.W_dot),
                    #Expander
                    ('T_exp_su','C',self.Expander.T_su),
                    ('T_exp_ex','C',self.Expander.T_ex-273.15),
                    ('P_exp_su','kPa',self.Expander.P_su),
                    ('P_exp_ex','kPa',self.Expander.P_ex),
                    ('mdot_exp','kg/s',self.Expander.m_dot),
                    ('power_exp','W',self.Expander.W_dot),
                    ('expander isentropic efficiency','-',self.Expander.eta_s),
                    ('delta_T_subcool','K',self.DELTAT_sc_cd),
                    ('delta_T_superheat','K',self.Evaporator.Tout_c-self.Evaporator.Tdew_c),
                    ('eta_cyc','-',self.eta_cycle),
                    ('eta_II_finite','-',self.eta_cycle_II),
                    #Charge
                    ('Charge','kg',self.Charge)
                    ]
        except AttributeError:
            return []
    
    def Calculate(self,Inputs,Tdew_evap,Tdew_cond):
         
        psat_cond=PropsSI('P','T',Tdew_cond+273.15,'Q',1,Inputs['Ref'])/1000
        psat_evap=PropsSI('P','T',Tdew_evap+273.15,'Q',1,Inputs['Ref'])/1000
          
        self.Tdew_cond=Tdew_cond
        self.Tdew_evap=Tdew_evap
        self.p_cd=psat_cond
        self.p_ev=psat_evap
        
        if  self.Tdew_cond>self.Tdew_evap:
            self.Tdew_cond=20
            self.Tdew_evap=70
        else:
            pass 
            
    
        ##Pump
        #M are regression coefficients to compute the mass flow rate.  
        #'eta' are regression coefficients to compute the isentropic efficiency.
        params={
                'Ref': Inputs['Ref'],
                'pin_r': psat_cond,
                'pout_r':psat_evap,
                'f_pp':Inputs['f_pp'],
                'Tin': self.Tdew_cond - Inputs['DT_sc']
              }
        self.Pump.Update(**params)
        self.Pump.Calculate()

        ##PumpEx Line
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_evap,
                'hin':PropsSI('H','T',self.Pump.Tout,'P',psat_evap*1000,Inputs['Ref']),
                'mdot':self.Pump.m_dot,
                #geometric parameters
                'OD':1.5*convert('in','m'),
                'ID':1.3*convert('in','m'),
                'L':3*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetPumpEx.Update(**params)
        self.LineSetPumpEx.Calculate()
        
        ##Evaporapor
        params={
                
                'Ref_c':Inputs['Ref'],
                'mdot_c':self.Pump.m_dot,
                'pin_c':psat_evap,
                'Tin_c': self.Pump.Tout-273.15,
                'hin_c':PropsSI('H','T',self.Pump.Tout,'P',psat_evap*1000,Inputs['Ref']),
                       
                'Ref_h':Inputs['Ref_h'],
                'mdot_h':Inputs['mdot_h'],
                'pin_h':Inputs['pin_h']*1000,
                'hin_h':PropsSI('H','T',Inputs['Tin_h']+273.15,'P',Inputs['pin_h']*1000,Inputs['Ref_h']),
            }
        self.Evaporator.Update(**params)
        self.Evaporator.Calculate()

        ##ExpanderSu Line
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_evap,
                'hin':PropsSI('H','T',self.Evaporator.Tout_c,'P',psat_evap*1000,Inputs['Ref']),
                'mdot':self.Expander.m_dot,
                #geometric parameters
                'OD':2*convert('in','m'),
                'ID':1.8*convert('in','m'),
                'L':3*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetExpSu.Update(**params)
        self.LineSetExpSu.Calculate()


        ##Expander
        params={
                'P_su':psat_evap,
                'P_ex':psat_cond,
                'h_su':self.Evaporator.hout_c,
                'inputs':'Psu_Pex_h'
                }
        self.Expander.Update(**params)
        self.Expander.Calculate()       
        
        ##CondSu Line
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_cond,
                'hin':PropsSI('H','T',self.Expander.T_ex,'P',psat_cond*1000+100,Inputs['Ref']),
                'mdot':self.Pump.m_dot,
                #geometric parameters
                'OD':2*convert('in','m'),
                'ID':1.8*convert('in','m'),
                'L':3*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetExpEx.Update(**params)
        self.LineSetExpEx.Calculate()
        
        ##Condenser
        params={
                'Ref_c':Inputs['Ref_c'],
                'mdot_c':Inputs['mdot_c'],
                'pin_c':Inputs['pin_c'],
                'hin_c':PropsSI('H','T',Inputs['Tin_c']+273.15,'P',Inputs['pin_c']*1000,Inputs['Ref_c']),
                
                'Ref_h':Inputs['Ref'],
                'mdot_h':self.Expander.m_dot,
                'pin_h':psat_cond,

                'hin_h':PropsSI('H','T',self.Expander.T_ex,'P',psat_cond*1000+100,Inputs['Ref']),
            }
        self.Condenser.Update(**params)
        self.Condenser.Calculate()

        ##LiquidReceiverSu Line 
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_cond,
                'hin':PropsSI('H','T',self.Condenser.Tout_h,'P',(1+1e-5)*psat_cond*1000,Inputs['Ref']),
                'mdot':self.Expander.m_dot,
                #geometric parameters
                'OD':2*convert('in','m'),
                'ID':1.8*convert('in','m'),
                'L':3*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetCondEx.Update(**params)
        self.LineSetCondEx.Calculate()

        ##Liquid Receiver
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_cond,
                'Tin': self.Condenser.Tout_h,
                'ID':0.35,
                'h_receiver': 0.9,
                'h_ports': 0.5
              }
        self.LiquidReceiver.Update(**params)
        self.LiquidReceiver.Calculate()    

        ##PumpSu Line
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_cond,
                'hin':PropsSI('H','T',self.Condenser.Tout_h,'P',psat_cond*1000+1,Inputs['Ref']),
                'mdot':self.Expander.m_dot,
                #geometric parameters
                'OD':1.5*convert('in','m'),
                'ID':1.3*convert('in','m'),
                'L':( 11 + 3)*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetPumpSu.Update(**params)
        self.LineSetPumpSu.Calculate()
  
        
        ############################
        #      Cycle Performance
        ############################

        #cycle efficiency
        self.Ref = Inputs['Ref']
        self.T_amb = 25  #Reference ambient temperature[C]
        
        level = (self.Condenser.hsatV_h - self.LiquidReceiver.hin)/(self.Condenser.hsatV_h - self.Condenser.hsatL_h)*(self.LiquidReceiver.rho_in/self.Condenser.rhosatL_h)

        #Back-Work Ratio (BWR)
        self.BWR = self.Pump.W_dot/self.Expander.W_dot
        
        #Heat recovery efficiency
        #From "Dynamic modeling and optimal control strategy of waste heat recovery Organic Rankine Cycles" Quoilin
        self.h_amb = PropsSI('H','T',self.T_amb+273.15,'P',self.Evaporator.pin_h*1000,self.Evaporator.Ref_h)#*1000
        self.Q_maxtoamb = self.Evaporator.mdot_h*(self.Evaporator.hin_h - self.h_amb)

        self.epsilon_hr = self.Evaporator.Q/self.Q_maxtoamb
        
        #first-law efficiency
        self.eta_cycle=(self.Expander.W_dot-self.Pump.W_dot)/self.Evaporator.Q
        self.eta_overall = self.eta_cycle*self.epsilon_hr
        
        
        #Second-law efficiency
        self.Exergy_hs=self.Evaporator.mdot_h*((self.Evaporator.hin_h-self.Condenser.hin_c)-(Inputs['Tin_c']+273.15)*(PropsSI('S','H',self.Evaporator.hin_h,'P',self.Evaporator.pin_h*1000,self.Evaporator.Ref_h)-PropsSI('S','H',self.Condenser.hin_c,'P',self.Condenser.pin_c*1000,self.Condenser.Ref_c)))
 
        self.eta_cycle_II=self.Expander.W_dot/self.Exergy_hs
         
        #residuals    
        resid=np.zeros((2))

        resid[0]=self.Pump.W_dot+self.Evaporator.Q-self.Expander.W_dot-self.Condenser.Q
        resid_mf=self.Pump.m_dot-self.Expander.m_dot
        self.DELTAT_sc_cd=(self.Tdew_cond-(self.Condenser.Tout_h-273.15)) 
        self.DELTAT_sh_ev= self.Evaporator.Tout_c-self.Evaporator.Tdew_c

        self.Charge = self.LineSetPumpEx.Charge + self.LineSetPumpSu.Charge + self.Evaporator.Charge_c + self.Condenser.Charge_h  + self.LineSetExpEx.Charge + self.LineSetExpSu.Charge  +  self.LineSetCondEx.Charge + self.LiquidReceiver.Charge_Tank
        
        e_DT_sc=Inputs['DT_sc']-self.DELTAT_sc_cd

        if Inputs['Solver'] == 'Subcooling':
            resid[1]=Inputs['DT_sc']-self.DELTAT_sc_cd
        
        elif Inputs['Solver'] == 'Charge':
            resid[1]=self.Charge- Inputs['Tot_Charge']



        self.resid_EB=resid[0]
        self.resid_mf=resid_mf
        self.resid_DT_sc=e_DT_sc
        self.resid_charge = self.Charge- Inputs['Tot_Charge']
        
        return resid
    


    def Calculate_Regen(self,Inputs,Tdew_evap,Tdew_cond,hin_reg_h):

        if  Tdew_evap > PropsSI(Inputs['Ref'],'Tcrit') - 273.15:
            Tdew_evap = Inputs['Tin_h'] - 50
            print 'Tdew_evep > Tcrit' 

        elif  Tdew_evap < PropsSI(Inputs['Ref'],'Tmin') - 273.15:
            Tdew_evap = Inputs['Tin_h'] - 50
            print 'Tdew_evep < Tmin'         

        else:
            pass
        
        if Tdew_cond > PropsSI(Inputs['Ref'],'Tcrit') - 273.15:
            Tdew_cond = Inputs['Tin_c'] +5
            print 'Tdew_cond > Tcrit'

        elif Tdew_cond < PropsSI(Inputs['Ref'],'Tmin') - 273.15:
            Tdew_cond = Inputs['Tin_c'] +5
            print 'Tdew_cond < Tmin' 
    
        else:
            pass
        
        if  Tdew_cond>Tdew_evap:
            print "Im here"
            Tdew_cond=20
            Tdew_evap=70
        else:
            pass 

        psat_cond=PropsSI('P','T',Tdew_cond+273.15,'Q',1,Inputs['Ref'])/1000
        psat_evap=PropsSI('P','T',Tdew_evap+273.15,'Q',1,Inputs['Ref'])/1000

        self.Tdew_cond=Tdew_cond
        self.Tdew_evap=Tdew_evap
        self.p_cd=psat_cond
        self.p_ev=psat_evap
        

        
        #Control variable
        self.hin_reg_h = hin_reg_h


        ##Pump
        #M are regression coefficients to compute the mass flow rate.  
        #'eta' are regression coefficients to compute the isentropic efficiency.
        params={
                'Ref': Inputs['Ref'],
                'pin_r': psat_cond,
                'pout_r':psat_evap,
                'f_pp':Inputs['f_pp'],
                'Tin': self.Tdew_cond - Inputs['DT_sc']
              }
        self.Pump.Update(**params)
        self.Pump.Calculate()

        ##PumpEx Line
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_evap,
                'hin':PropsSI('H','T',self.Pump.Tout,'P',psat_evap*1000,Inputs['Ref']),
                'mdot':self.Pump.m_dot,
                #geometric parameters
                'OD':1.5*convert('in','m'),
                'ID':1.3*convert('in','m'),
                'L':3*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetPumpEx.Update(**params)
        self.LineSetPumpEx.Calculate()


        ##Regenerator
        params={
                'Ref_c':Inputs['Ref'],
                'mdot_c':self.Pump.m_dot,
                'pin_c':self.p_ev,
                'Tin_c': self.Pump.Tout-273.15,
                'hin_c':PropsSI('H','T',self.Pump.Tout,'P',psat_evap*1000+100,Inputs['Ref']),

                'Ref_h':Inputs['Ref'],
                'mdot_h':self.Expander.m_dot,
                'pin_h':self.p_cd,
                'hin_h':self.hin_reg_h,   #[J/kg]
            }
        self.Regenerator.Update(**params)
        self.Regenerator.Calculate()

        ##Evaporapor
        params={
                
                'Ref_c':Inputs['Ref'],
                'mdot_c':self.Pump.m_dot,
                'pin_c':psat_evap,
                'Tin_c': self.Regenerator.Tout_c,
                'hin_c':PropsSI('H','T',self.Regenerator.Tout_c,'P',psat_evap*1000+100,Inputs['Ref']),
                       
                'Ref_h':Inputs['Ref_h'],
                'mdot_h':Inputs['mdot_h'],
                'pin_h':Inputs['pin_h'],
                'hin_h':PropsSI('H','T',Inputs['Tin_h']+273.15,'P',Inputs['pin_h']*1000,Inputs['Ref_h']),
                }
        self.Evaporator.Update(**params)
        self.Evaporator.Calculate()
 
    
        ##ExpanderSu Line
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_evap,
                'hin':PropsSI('H','T',self.Evaporator.Tout_c,'P',psat_evap*1000+100,Inputs['Ref']),
                'mdot':self.Pump.m_dot,
                #geometric parameters
                'OD':2*convert('in','m'),
                'ID':1.8*convert('in','m'),
                'L':3*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetExpSu.Update(**params)
        self.LineSetExpSu.Calculate()



        ##Expander
        params={
                'P_su':psat_evap,
                'P_ex':psat_cond,
                'h_su':self.Evaporator.hout_c,
                'inputs':'Psu_Pex_h'
                }
        self.Expander.Update(**params)
        self.Expander.Calculate()  



        ##CondSu Line

        params={
                'Ref':Inputs['Ref'],
                'pin':psat_cond,
                'hin':PropsSI('H','T',self.Regenerator.Tout_h,'P',psat_cond*1000+100,Inputs['Ref']),
                'mdot':self.Pump.m_dot,
                #geometric parameters
                'OD':2*convert('in','m'),
                'ID':1.8*convert('in','m'),
                'L':3*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetExpEx.Update(**params)
        self.LineSetExpEx.Calculate()
        
        ##Condenser
        params={
                'Ref_c':Inputs['Ref_c'],
                'mdot_c':Inputs['mdot_c'],
                'pin_c':Inputs['pin_c'],
                'hin_c':PropsSI('H','T',Inputs['Tin_c']+273.15,'P',Inputs['pin_c']*1000,Inputs['Ref_c']),#*1000,
                
                'Ref_h':Inputs['Ref'],
                'mdot_h':self.Expander.m_dot,
                'pin_h':psat_cond,

                'hin_h':PropsSI('H','T',self.Regenerator.Tout_h,'P',psat_cond*1000-100,Inputs['Ref']),#*1000,
                }
        
        self.Condenser.Update(**params)
        self.Condenser.Calculate()
        
        ##LiquidReceiverSu Line 
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_cond,
                'hin':PropsSI('H','T',self.Condenser.Tout_h,'P',(1+1e-5)*psat_cond*1000,Inputs['Ref']),
                'mdot':self.Expander.m_dot,
                #geometric parameters
                'OD':2*convert('in','m'),
                'ID':1.8*convert('in','m'),
                'L':3*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetCondEx.Update(**params)
        self.LineSetCondEx.Calculate()


        ##Liquid Receiver
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_cond,
                'Tin': self.Condenser.Tout_h,
                'ID':0.30,
                'h_receiver': 0.7,
                'h_ports': 0.5
              }
        self.LiquidReceiver.Update(**params)
        self.LiquidReceiver.Calculate()    




        ##PumpSu Line
        params={
                'Ref':Inputs['Ref'],
                'pin':psat_cond,
                'hin':PropsSI('H','T',self.Condenser.Tout_h,'P',psat_cond*1000+100,Inputs['Ref']),#*1000,
                #'mdot':self.Expander_2.m_dot,
                'mdot':self.Expander.m_dot,
                #geometric parameters
                'OD':1.5*convert('in','m'),
                'ID':1.3*convert('in','m'),
                'L':( 11 + 3)*convert('f','m'),
                'k_tube':300,
                't_insul':1*convert('in','m'),
                'k_insul':0.05,
                'h_air':10,
                'T_air':300
                }
        self.LineSetPumpSu.Update(**params)
        self.LineSetPumpSu.Calculate()



        ############################
        #      ORC with Regen Performance
        ############################
        
        
        #cycle efficiency
        self.Ref = Inputs['Ref']
        self.T_amb = 25  #Reference ambient temperature[C]
        
        level = (self.Condenser.hsatV_h - self.LiquidReceiver.hin)/(self.Condenser.hsatV_h - self.Condenser.hsatL_h)*(self.LiquidReceiver.rho_in/self.Condenser.rhosatL_h)
        

        #Back-Work Ratio (BWR)
        self.BWR = self.Pump.W_dot/self.Expander.W_dot
        
        #Heat recovery efficiency
        #From "Dynamic modeling and optimal control strategy of waste heat recovery Organic Rankine Cycles" Quoilin
        self.h_amb = PropsSI('H','T',self.T_amb+273.15,'P',self.Evaporator.pin_h*1000,self.Evaporator.Ref_h)#*1000
        self.Q_maxtoamb = self.Evaporator.mdot_h*(self.Evaporator.hin_h - self.h_amb)

        self.epsilon_hr = self.Evaporator.Q/self.Q_maxtoamb
        
        #first-law efficiency
        self.eta_cycle=(self.Expander.W_dot-self.Pump.W_dot )/self.Evaporator.Q
        self.eta_overall = self.eta_cycle*self.epsilon_hr
        
        
        #Second-law efficiency
        self.Exergy_hs=self.Evaporator.mdot_h*((self.Evaporator.hin_h-self.Condenser.hin_c)-(Inputs['Tin_c']+273.15)*(PropsSI('S','H',self.Evaporator.hin_h,'P',self.Evaporator.pin_h*1000,self.Evaporator.Ref_h)-PropsSI('S','H',self.Condenser.hin_c,'P',self.Condenser.pin_c*1000,self.Condenser.Ref_c)))
 
        self.eta_cycle_II=self.Expander.W_dot/self.Exergy_hs
         
        #residuals    
        resid=np.zeros((3))

        ##Residuals        
        """
        resid[0]: overall energy balance
        resid[1]: subcooling or refrigerant charge
        resid[2]: enthalpy at expander outlet
        """
        resid[0]=self.Pump.W_dot+self.Evaporator.Q-self.Expander.W_dot-self.Condenser.Q 
        resid_mf=self.Pump.m_dot-self.Expander.m_dot 
        self.DELTAT_sc_cd=(self.Tdew_cond-(self.Condenser.Tout_h-273.15)) 
        self.DELTAT_sh_ev= self.Evaporator.Tout_c-self.Evaporator.Tdew_c
        e_DT_sc=Inputs['DT_sc']-self.DELTAT_sc_cd
        self.Charge = self.LineSetPumpEx.Charge + self.LineSetPumpSu.Charge + self.Evaporator.Charge_c + self.Condenser.Charge_h  + self.LineSetExpEx.Charge + self.LineSetExpSu.Charge  +  self.LineSetCondEx.Charge + self.LiquidReceiver.Charge_Tank + self.Regenerator.Charge_c + self.Regenerator.Charge_h
        
        if Inputs['Solver'] == 'Subcooling':
        
            resid[1]=Inputs['DT_sc']-self.DELTAT_sc_cd
        
        elif Inputs['Solver'] == 'Charge':
        
            resid[1]=self.Charge- Inputs['Tot_Charge']
        
        resid[2]= (self.hin_reg_h-self.Expander.hout) 

        self.resid_EB=resid[0]
        self.resid_mf=resid_mf
        self.resid_DT_sc=e_DT_sc
        self.resid_charge = self.Charge- Inputs['Tot_Charge']
   
        return resid




##########################################################################
#    Solvers
########################################################################              
           
    def PreconditionedSolve(self,Inputs):
        
        iteration_precond=[]
        resid1_precond = []
        resid2_precond = []
        
        def OBJECTIVE_ORC(x):
            

            iteration=my_counter.next()
           # print iteration
            """
            A wrapper function to convert input vector for fsolve to the proper form for the solver
            """
        
            try:
                resids=self.Calculate(Inputs,Tdew_evap=float(x[0]),Tdew_cond=float(x[1]))
                iteration_precond.append(iteration)
                resid1_precond.append(resids[0])
                resid2_precond.append(resids[1])
         
            
                if Inputs['Solver'] == 'Subcooling':
                    print 'EnergyBalance,DT_sc:',resids[0],resids[1]            
                elif Inputs['Solver'] == 'Charge':
                    print 'EnergyBalance,Charge:',resids[0],resids[1]
            
            except ValueError:
                raise
            #print it

            self.iteration_precond = np.array(iteration_precond)
            self.resid1_precond = np.array(resid1_precond)
            self.resid2_precond = np.array(resid2_precond)

            return resids


        
        # Use the preconditioner to determine a reasonably good starting guess
        
        Tdew_cond_init,Tdew_evap_init=ORC_Preconditioner(self,Inputs)

        #Actually run the solver to get the solution
        x=fsolve(OBJECTIVE_ORC,[Tdew_cond_init,Tdew_evap_init],full_output=True,xtol=1e-6)
        #x=MultiDimNewtRaph(OBJECTIVE_ORC,[Tdew_cond_init,Tdew_evap_init],dx=1e-6,args=(),ytol=1e-4,w=1.0,JustOneStep=False)
        
    
    
    
    def PreconditionedSolve_ORCRegen(self,Inputs):
        
        iteration_precond=[]
        resid1_precond = []
        resid2_precond = []
        resid3_precond = []
        
        def OBJECTIVE_ORCRegen(x):
            

            iteration=my_counter.next()
           # print iteration
            """
            A wrapper function to convert input vector for fsolve to the proper form for the solver
            """
        
            try:
                resids=self.Calculate_Regen(Inputs,Tdew_evap=float(x[0]),Tdew_cond=float(x[1]),hin_reg_h=float(x[2]))
                iteration_precond.append(iteration)
                resid1_precond.append(resids[0])
                resid2_precond.append(resids[1])
                resid3_precond.append(resids[2])
                
                if Inputs['Solver'] == 'Subcooling':
                    print 'EnergyBalance,DT_sc,hin_reg_h:',resids[0],resids[1],resids[2]             
                elif Inputs['Solver'] == 'Charge':
                    print 'EnergyBalance,Charge,hin_reg_h:',resids[0],resids[1],resids[2]      
            
            except ValueError:
                raise
            #print it

            self.iteration_precond = np.array(iteration_precond)
            self.resid1_precond = np.array(resid1_precond)
            self.resid2_precond = np.array(resid2_precond)
            self.resid3_precond = np.array(resid3_precond)
            #print 'iteration_precond',self.iteration_precond
            return resids


        
        # Use the preconditioner to determine a reasonably good starting guess
        Tdew_cond_init,Tdew_evap_init,hin_reg_h_init=ORCRegen_Preconditioner(self,Inputs)
        
        #Actually run the solver to get the solution
        x=fsolve(OBJECTIVE_ORCRegen,[Tdew_cond_init,Tdew_evap_init,hin_reg_h_init],full_output=True,xtol=1e-6)






        
    def ConvergencePlot(self):    
        
        matplotlib.rc('text', usetex=True)
        plt.rc('font', family='serif')
        fig=figure(figsize=(8,6))
        ax=fig.add_subplot(1,1,1)

        ax2 = ax.twinx()
        
        plot1 = ax.plot(self.iteration_precond, fabs(self.resid1_precond),'bo-',linewidth=2,markersize = 10,markerfacecolor="blue", markeredgewidth=1, markeredgecolor="black", label = 'Energy Balance')
        plot2 = ax2.plot(self.iteration_precond, fabs(self.resid2_precond),'r-^',linewidth=2,markersize = 10,markerfacecolor="red", markeredgewidth=1, markeredgecolor="black" ,label = 'Ref Charge')
        
        ax.axhline(y=1e-5,color='b',ls='dashed', linewidth=2.5)
        ax2.axhline(y=1e-5,color='r',ls='dashed', linewidth=2.5)
        
        # Only one legend with different axis plots
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()

        ax2.legend(lines + lines2,labels + labels2,loc=0, shadow = True,fancybox=True,prop={'size':15})
        
        ax.set_xlabel('Iteration [-]',fontsize=18)
        ax.set_ylabel(r'$|Resid_{Energy Balance}|$',fontsize=18)
        ax2.set_ylabel(r'$|Resid_{Charge}|$',fontsize=18)
        plt.title('Resids Convergence',fontsize=18)
        
        ax.set_yscale('log')
        ax.tick_params(axis='y',which='minor',bottom='off')
        ax2.set_yscale('log')
        ax2.tick_params(axis='y',which='minor',bottom='off')
        
        
         
        for tickx in ax.xaxis.get_major_ticks():
            tickx.label.set_fontsize(20)
        for ticky in ax.yaxis.get_major_ticks():
            ticky.label.set_fontsize(20)
        for ticky in ax2.yaxis.get_major_ticks():
            ticky.label.set_fontsize(20)        
        
        print 'Saving convergence plot...'
        fig.savefig('ConvergencePlot.png',dpi=600)
                   
              
    def BaselineTs(self,Inputs,**kwargs):

        matplotlib.rc('text', usetex=True)
        plt.rc('font', family='serif')
        fig=figure(figsize=(8,6))
        ax=fig.add_subplot(1,1,1)
    
        if 'Tmin' in kwargs:
            Tmin=float(kwargs['Tmin'])
        else:
            Tmin=220.
    
        Tsat = linspace(Tmin,PropsSI(Inputs['Ref'],'Tcrit')-0.0001,1000)
        (ssatL,psatL,ssatV,psatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
        for i in arange(len(Tsat)):
            ssatL[i] = PropsSI('S','T',Tsat[i],'Q',0,Inputs['Ref'])/1000
            ssatV[i] = PropsSI('S','T',Tsat[i],'Q',1,Inputs['Ref'])/1000
    
        ax.plot(ssatL,Tsat,'k')
        ax.plot(ssatV,Tsat,'k')

        
        epss=0.05
        epsT=5
    
        # Base Cycle
        #Expander
        s_exp=r_[self.Expander.s_su/1000, self.Expander.s_out/1000]
        T_exp=r_[self.Expander.T_su+273.15, self.Expander.T_ex]
        ax.plot(s_exp,T_exp,'k-')
        ax.plot(s_exp,T_exp,'ko-',linewidth=2,markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_exp[0]+epss/2,T_exp[0]+epsT,'3',va='top',ha='left')
    
        s_exp_s=r_[self.Expander.s_su/1000, self.Expander.s_su/1000]
        T_exp_s=r_[self.Expander.T_su+273.15, self.Expander.T_ex_s]
        
        ax.plot(s_exp_s,T_exp_s,'k--',linewidth = 1.5)
        ax.plot(s_exp_s[1],T_exp_s[1],'ko-',linewidth=2,markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_exp_s[1],T_exp_s[1]-epsT,'4s',va='top',ha='left')

        #Condenser
        T_cond = linspace(self.Condenser.Tout_h,self.Expander.T_ex,200)
        s_cond=0.0*T_cond
        p_cond=self.p_cd
        for i in range(len(T_cond)):
            s_cond[i]=PropsSI('S','T',T_cond[i],'P',p_cond*1000+100,self.Ref)/1000
        ax.plot(s_cond,T_cond,'k-',linewidth = 2)
        ax.plot(r_[s_cond[0],s_cond[-1]],r_[T_cond[0],T_cond[-1]],'ko',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_cond[-1]+epss/2,T_cond[-1]+epsT,'4',va='top',ha='left')
        
        s_EXV=r_[self.Condenser.sout_h/1000, self.Pump.s_in]
        T_EXV=r_[self.Condenser.Tout_h, self.Pump.Tin]

        ax.text(s_EXV[0]+epss,T_EXV[0]+epsT,'1',va='top',ha='left')
    
        #Pump
        s_pump=r_[self.Pump.s_in,self.Pump.s_out]
        T_pump=r_[self.Pump.Tin+273.15, self.Pump.Tout]
        ax.plot(s_pump,T_pump,'b-',linewidth=2)
        ax.plot([self.Pump.s_in,self.Pump.s_in],[self.Pump.Tin+273.15,self.Pump.Tout_s],'ko--',linewidth=2,markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.plot([self.Pump.s_in,self.Pump.s_out],[self.Pump.Tout_s,self.Pump.Tout],'ko--',linewidth=2,markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.plot(s_pump[0],T_pump[0],'bo',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.plot(s_pump[1],T_pump[1],'bo',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_pump[1]-epss,T_pump[1]+2*epsT,'2',va='top',ha='left')
         
         #Evaporator
         #Zone1
        s_evap_zone1=r_[self.Pump.s_out,self.Evaporator.s_c_liq]
        T_evap_zone1=r_[self.Pump.Tout_s,self.Evaporator.Tdew_c]
        ax.plot(s_evap_zone1,T_evap_zone1,'k-',linewidth = 2)
        ax.plot(s_evap_zone1[1],T_evap_zone1[1],'ko',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_evap_zone1[1],T_evap_zone1[1]-epsT,'2i',va='top',ha='left')    
        #Zone2
        s_evap_zone2=r_[self.Evaporator.s_c_liq,self.Evaporator.s_c_vap]
        T_evap_zone2=r_[self.Evaporator.Tdew_c,self.Evaporator.Tdew_c]
        ax.plot(s_evap_zone2,T_evap_zone2,'k-',linewidth = 2)
        ax.plot(s_evap_zone2[1],T_evap_zone2[1],'ko',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_evap_zone2[1]-1.5*epss,T_evap_zone2[1]+2*epsT,'2ii',va='top',ha='left') 
        #Zone3
        #s_evap_zone3=r_[self.Evaporator.s_c_vap,self.Expander_2.s_in]
        s_evap_zone3=r_[self.Evaporator.s_c_vap,self.Expander.s_su/1000]
        T_evap_zone3=r_[self.Evaporator.Tdew_c,self.Evaporator.Tout_c]
        ax.plot(s_evap_zone3,T_evap_zone3,'k-',linewidth=2)

        #Regenerator
        T_regen_out = self.Regenerator.Tout_c
        s_regen_out = self.Regenerator.sout_c/1000
        T_regen_out_meas = 74.77 + 273.15
        s_regen_out_meas = PropsSI('S','P',895.072*1000 + 100,'T',T_regen_out,Inputs['Ref'])/1000
        ax.plot(s_regen_out,T_regen_out,'bo',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        
        #Hot source
        T_evap_h = r_[self.Evaporator.Tin_h,self.Evaporator.Tout_h]
        s_evap_h = r_[self.Expander.s_su/1000,self.Regenerator.sout_c/1000]
        ax.plot(s_evap_h,T_evap_h,'ro-',markersize = 7,markerfacecolor="red", markeredgewidth=2, markeredgecolor="black")

        #Cold source
        T_cond_c = r_[self.Condenser.Tin_c,self.Condenser.Tout_c]
        s_cond_c = r_[self.Condenser.sout_h/1000,self.Regenerator.sout_h/1000]
        ax.plot(s_cond_c,T_cond_c,'bo-',markersize = 7,markerfacecolor="blue", markeredgewidth=2, markeredgecolor="black")        
        
        
    #     pylab.xticks([-1.8,-1.4,-1.0,-0.6,-0.2])
        plt.ylim([250,450])
        plt.xlim([0.8,2.0])
        
        plt.xlabel('Entropy [kJ/(kg$\cdot$K)]', fontsize = 17)
        plt.ylabel('Temperature [K]',fontsize = 17)
        plt.title(Inputs['Ref'],fontsize = 17)
        
        for tickx in ax.xaxis.get_major_ticks():
            tickx.label.set_fontsize(20)
        for ticky in ax.yaxis.get_major_ticks():
            ticky.label.set_fontsize(20)     
        
        print 'Saving T-s diagram ...'
        fig.savefig('Ts.png',dpi=600)
        fig.savefig('Ts.pdf')
        
    def BaselineTs_Celsius(self,Inputs,**kwargs):

        #matplotlib.rc('text', usetex=True)
        #plt.rc('font', family='serif')
        fig=figure(figsize=(8,6))
        ax=fig.add_subplot(1,1,1)
    
        if 'Tmin' in kwargs:
            Tmin=float(kwargs['Tmin'])
        else:
            Tmin=220.
    
        Tsat = linspace(Tmin,PropsSI(Inputs['Ref'],'Tcrit')-0.0001,1000)
        (ssatL,psatL,ssatV,psatV)=(0.0*Tsat,0.0*Tsat,0.0*Tsat,0.0*Tsat)
        for i in arange(len(Tsat)):
            ssatL[i] = PropsSI('S','T',Tsat[i],'Q',0,Inputs['Ref'])/1000
            ssatV[i] = PropsSI('S','T',Tsat[i],'Q',1,Inputs['Ref'])/1000
    
        ax.plot(ssatL,Tsat-273.15,'k')
        ax.plot(ssatV,Tsat-273.15,'k')

        
        epss=0.05
        epsT=5
    
        # Base Cycle
        #Expander
        s_exp=r_[self.Expander.s_su/1000, self.Expander.s_out/1000]
        T_exp=r_[self.Expander.T_su, self.Expander.T_ex-273.15]
        ax.plot(s_exp,T_exp,'k-')
        ax.plot(s_exp,T_exp,'ko-',linewidth=2,markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_exp[0]+epss/2,T_exp[0]+epsT,'3',va='top',ha='left')
    
        s_exp_s=r_[self.Expander.s_su/1000, self.Expander.s_su/1000]
        T_exp_s=r_[self.Expander.T_su, self.Expander.T_ex_s-273.15]
        
        ax.plot(s_exp_s,T_exp_s,'k--',linewidth = 1.5)
        ax.plot(s_exp_s[1],T_exp_s[1],'ko-',linewidth=2,markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_exp_s[1],T_exp_s[1]-epsT,'4s',va='top',ha='left')

        #Condenser
        T_cond = linspace(self.Condenser.Tout_h,self.Expander.T_ex,200)-273.15
        s_cond=0.0*T_cond
        p_cond=self.p_cd
        for i in range(len(T_cond)):
            s_cond[i]=PropsSI('S','T',T_cond[i]+273.15,'P',p_cond*1000+100,self.Ref)/1000
        ax.plot(s_cond,T_cond,'k-',linewidth = 2)
        ax.plot(r_[s_cond[0],s_cond[-1]],r_[T_cond[0],T_cond[-1]],'ko',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_cond[-1]+epss/2,T_cond[-1]+epsT,'4',va='top',ha='left')
        
        s_EXV=r_[self.Condenser.sout_h/1000, self.Pump.s_in]
        T_EXV=r_[self.Condenser.Tout_h-273.15, self.Pump.Tin-273.15]

        ax.text(s_EXV[0]+epss,T_EXV[0]+epsT,'1',va='top',ha='left')
    
        #Pump
        s_pump=r_[self.Pump.s_in,self.Pump.s_out]
        T_pump=r_[self.Pump.Tin, self.Pump.Tout-273.15]
        ax.plot(s_pump,T_pump,'b-',linewidth=2)
        ax.plot([self.Pump.s_in,self.Pump.s_in],[self.Pump.Tin,self.Pump.Tout_s-273.15],'ko--',linewidth=2,markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.plot([self.Pump.s_in,self.Pump.s_out],[self.Pump.Tout_s-273.15,self.Pump.Tout-273.15],'ko--',linewidth=2,markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.plot(s_pump[0],T_pump[0],'bo',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.plot(s_pump[1],T_pump[1],'bo',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_pump[1]-epss,T_pump[1]+2*epsT,'2',va='top',ha='left')
         
         #Evaporator
         #Zone1
        s_evap_zone1=r_[self.Pump.s_out,self.Evaporator.s_c_liq]
        T_evap_zone1=r_[self.Pump.Tout_s-273.15,self.Evaporator.Tdew_c-273.15]
        ax.plot(s_evap_zone1,T_evap_zone1,'k-',linewidth = 2)
        ax.plot(s_evap_zone1[1],T_evap_zone1[1],'ko',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_evap_zone1[1],T_evap_zone1[1]-epsT,'2i',va='top',ha='left')    
        #Zone2
        s_evap_zone2=r_[self.Evaporator.s_c_liq,self.Evaporator.s_c_vap]
        T_evap_zone2=r_[self.Evaporator.Tdew_c-273.15,self.Evaporator.Tdew_c-273.15]
        ax.plot(s_evap_zone2,T_evap_zone2,'k-',linewidth = 2)
        ax.plot(s_evap_zone2[1],T_evap_zone2[1],'ko',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        ax.text(s_evap_zone2[1]-1.5*epss,T_evap_zone2[1]+2*epsT,'2ii',va='top',ha='left') 
        #Zone3
        #s_evap_zone3=r_[self.Evaporator.s_c_vap,self.Expander_2.s_in]
        s_evap_zone3=r_[self.Evaporator.s_c_vap,self.Expander.s_su/1000]
        T_evap_zone3=r_[self.Evaporator.Tdew_c-273.15,self.Evaporator.Tout_c-273.15]
        ax.plot(s_evap_zone3,T_evap_zone3,'k-',linewidth=2)

        #Regenerator
        T_regen_out = self.Regenerator.Tout_c
        s_regen_out = self.Regenerator.sout_c/1000
        T_regen_out_meas = 74.77 + 273.15
        s_regen_out_meas = PropsSI('S','P',895.072*1000 + 100,'T',T_regen_out-273.15,Inputs['Ref'])/1000
        ax.plot(s_regen_out,T_regen_out-273.15,'bo',markersize = 7,markerfacecolor="none", markeredgewidth=2, markeredgecolor="black")
        
        #Hot source
        T_evap_h = r_[self.Evaporator.Tin_h-273.15,self.Evaporator.Tout_h-273.15]
        s_evap_h = r_[self.Expander.s_su/1000,self.Regenerator.sout_c/1000]
        ax.plot(s_evap_h,T_evap_h,'ro-',markersize = 7,markerfacecolor="red", markeredgewidth=2, markeredgecolor="black")

        #Cold source
        T_cond_c = r_[self.Condenser.Tin_c-273.15,self.Condenser.Tout_c-273.15]
        s_cond_c = r_[self.Condenser.sout_h/1000,self.Regenerator.sout_h/1000]
        ax.plot(s_cond_c,T_cond_c,'bo-',markersize = 7,markerfacecolor="blue", markeredgewidth=2, markeredgecolor="black")        
        
        
    #     pylab.xticks([-1.8,-1.4,-1.0,-0.6,-0.2])
        plt.ylim([10,160])
        plt.xlim([0.8,2.0])
        
        plt.xlabel('Entropy [kJ/(kg$\cdot$K)]', fontsize = 17)
        plt.ylabel('Temperature [C]',fontsize = 17)
        plt.title(Inputs['Ref'],fontsize = 17)
        
        for tickx in ax.xaxis.get_major_ticks():
            tickx.label.set_fontsize(20)
        for ticky in ax.yaxis.get_major_ticks():
            ticky.label.set_fontsize(20)     
        
        print 'Saving T-s diagram ...'
        fig.savefig('Ts.png',dpi=600)
        fig.savefig('Ts.pdf')        
