from __future__ import division


from CoolProp.CoolProp import PropsSI
import pylab
from ACHPTools import Write2CSV
from matplotlib.pyplot import plot, show, figure, semilogy, xlim, ylim, title, xlabel, ylabel, legend
from math import pi,exp,log,sqrt,tan,cos,sin
from scipy.optimize import brentq
from scipy.constants import g
import numpy as np

from PHEX_ASME2015 import PHEHXClass
from LineSet import LineSetClass

class LiquidReceiverClass():
    "Create Refrigerant buffer tank class"

    def __init__(self,**kwargs):
        #Load up the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    
        self.Condenser=PHEHXClass()
        
    def Update(self,**kwargs):
        #Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

    def OutputList(self):
        """
            Return a list of parameters for this component for further output            
            It is a list of tuples, and each tuple is formed of items with indices:
                [0] Description of value                
                [1] Units of value
                [2] The value itself
        """

        return [
            ('Liquid Receiver Total Volume','m3',self.Volume_tank),
            ('Liquid Receiver Total Charge','Kg',self.Charge_Tank),
            ('Inlet Temperature','K',self.Tin),
            ('Outlet Temperature','K',self.Tout),
            ('Inlet Pressure','kPa',self.pin),
            ('Inlet Density', 'kg/m3',self.rho_in),
            ('Outlet Pressure','kPa',self.pout)
           ]
           
           
                   
    def Calculate(self):
        
        """
        The liquid receiver acts as a damper in the cycle, absorbing the the mass flow rate 
        fluctuations. More concretely, a different explanation can be given.
        When the liquid receiver gets subcooled or saturated liquid at its top, it can be assumed to be
        in thermodynamic equilibrium at each time, because liquid and vapor have the same pressure when
        they enter it (indeed, if the reservoir isn't full, the vapor contained in it must be saturated, as it is in
        presence of liquid). In the inferior part of the tank, the mix of saturated and subcooled liquid (already
        present) allows the working fluid to exit it in a subcooled liquid state. The saturation pressure and
        temperature then reign then in the superior part of the reservoir. Thus, with this component, the
        charge fluctuations are litteraly absorbed, put to an equilibrium value, and the subcooling becomes
        null (this fact can't be stated in the presence of non-condensable gases).
        
        level = (h_v_sat - h)/(h_v_sat - h_l_sat)*(rho/rho_l_sat)
        
        """

        # Density [kg/m^3]

        self.rho_in=PropsSI('D','T',self.Tin, 'P', self.pin*1000+100, self.Ref)
        
        #Static pressure (rho*g*h) between inlet and outlet of the tank"
        
        self.pout=self.pin #+ (self.rho_in*g*self.h_ports)/1000
       # print 'LiquidReceiver.pout', self.pout
        self.Tout = self.Tin #no temperature gradient is observed in the reservoir.
        self.hin = PropsSI('H','T',self.Tin,'P',self.pin*1000+100,self.Ref) #J/kg

        """
        "Calculations"
        "x_ex_tank=0"	"due to the presence of non condensable gas  (air, due to leakage) in the working fluid, 
        "the liquid at the exit of the tank is not saturated..."

        #h_su_tank=h_ex_cd

        #V_ex_tank = m_dot/rho_ex_tank "Check V_dot_su_pump at the beginning of the file!!"
        """
        self.hout = PropsSI('H','T',self.Tout, 'P', self.pout*1000+100, self.Ref) #J/kg
        #print 'LiquidReceiver.hout', self.hout
        self.sout = PropsSI('S','T',self.Tout, 'P', self.pout*1000+100, self.Ref) #J/kg      
        
        #Calculate saturated values
        
        
        
        #Charge of the tank [kg]
        """
        The tank is characterized by an internal diameter and heigth (ID,h)
        and by the maximum level of refrigerant inside
        """
                
        self.Volume_tank = pi*self.ID**2/4.0*self.h_receiver
        self.Charge_Tank = self.Volume_tank * self.rho_in
        #self.Volume_ref = self.Charge_Tank/self.LiquidReceiver.rho_in



if __name__=='__main__':

    pin_list=[527.374817]
    Tin_list=[15.48]
    zip(pin_list,Tin_list)       
    for pin,Tin in zip(pin_list,Tin_list):
        kwds={
              'Ref':'R134A',
              'pin':pin,
              'Tin':Tin+273.15,
              'ID':0.3,
              'h_receiver': 1,
              'h_ports':0.5
              }
        LiquidReceiver=LiquidReceiverClass(**kwds)
        LiquidReceiver.Calculate()
        
        print  'Charge [kg]',LiquidReceiver.Charge_Tank
        print 'pin [kPa]', LiquidReceiver.pin
        print  'pout [kPa]',LiquidReceiver.pout
        print 'Receiver Volume [cm3]', LiquidReceiver.Volume_tank*1e6