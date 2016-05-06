from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from CoolProp.CoolProp import PropsSI
from Correlations import Tsat

class PumpClass():
    """
    Pump Model based on correlations obtained from experimental results
    """
    
    def __init__(self,**kwargs):
        #Load up the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
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
        
        return []
        
    def Calculate(self):
        
        #Local copies of coefficients
        W=self.W
           
        #Power
        
        #Compute the pressure difference between the outlet and the inlet of the pump
        self.DELTAP=self.pout_r-self.pin_r
        #Get the rated power for all pressure differences
        W_dot_rated=W[0]*self.DELTAP+W[1]
        #Speed ratio
        N_ratio=self.N/self.N_rated
        #Get the estimation of the power
        self.W_dot=(W[2]*N_ratio+W[3])*W_dot_rated  
        
        #Mass flow rate
        #Define the slope of the line corresponding to the exhaust temperature as a linear interpolation of the minimum and maximum slope
        slope=(self.slope_max-self.slope_min)/(self.p_max-self.p_min)*(self.pout_r-self.p_min)+self.slope_min
        #Define the intercept of the line corresponding to the exhaust temperature
        intercept=(self.intercept_max-self.intercept_min)/(self.p_max-self.p_min)*(self.pout_r-self.p_min)+self.intercept_min
        
        self.m_dot=slope*self.N+intercept 
        
        #Outlet state
        hin=PropsSI('H','T',self.Tin+273.15,'P',self.pin_r*1000,self.Ref)#*1000
        self.s_in = PropsSI('S','T',self.Tin+273.15,'P',self.pin_r*1000,self.Ref)/1000
        hout=hin+self.W_dot/self.m_dot
        self.Tout=PropsSI('T','H',hout,'P',self.pout_r*1000,self.Ref) #in K
        self.s_out=PropsSI('S','T',self.Tout,'P',self.pout_r*1000 + 100,self.Ref)/1000  
        self.Tout_s = PropsSI('T','S',self.s_in*1000,'P',self.pout_r*1000,self.Ref) 


if __name__=='__main__':

    """
    Example Diaphragm pump WANNER ENGINEERING
    """
    pin_r_list=[794.7276887,780.158035,784.3067128,808.239602,822.8122092,826.29617,887.1980418]
    pout_r_list=[1645.186859,1684.81582,1712.113611,1715.081928,1618.593683,1616.02753,1728.196266]
    N_list=[1099.97098,1099.809986,1099.72049,1099.818785,1099.743137,1099.450796,1099.270196]
    Tin_list=[15.4903837535014,15.3066340782123,15.5798263305322,15.7492877094972,15.5736862745098,15.7364804469274,15.0563305322129]
    W_list_meas = [235.4954587,236.254973,245.3089328,241.3617462,233.9065263,228.6898989,239.6439083]
    zip(pin_r_list,pout_r_list,N_list,Tin_list,W_list_meas)       
    for pin_r,pout_r,N,Tin,Wmeas in zip(pin_r_list,pout_r_list,N_list,Tin_list,W_list_meas):
        kwds={
              'W':[0.1096,114.34,1.0993,-0.0981],
              'Ref':'R134a',
              'pin_r':pin_r,
              'pout_r':pout_r,
              'N':N,
              'N_rated':995,
              'slope_min':0.000133504, #corresponding to the min outlet pressure
              'slope_max':0.000114377,  #corresponding to the max outlet pressure
              'intercept_min':0.004, #corresponding to the min outlet pressure
              'intercept_max':0.025, #corresponding to the max outlet pressure
              'p_min':700.4260866,
              'p_max':2659.623637,
              'Tin':Tin
              }
        Pump=PumpClass(**kwds)
        Pump.Calculate()
        
        print 'Calculated:',Pump.W_dot,'W','Measured:',Wmeas,'W'

