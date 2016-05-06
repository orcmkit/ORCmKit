from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from CoolProp.CoolProp import PropsSI
from Correlations import Tsat
import numpy as np

class PumpCentrifugalClass():
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
        
    
    
    def MassFlowPump(self,xy,pp0,pp1,pp2,pp3,pp4):

        f_pp = xy[0]
        deltaP_pp = xy[1]
        m_dot_pp = pp0 + pp1*deltaP_pp + pp2*deltaP_pp**2 + pp3*f_pp + pp4*f_pp**2
        return m_dot_pp
    
    def IsEffPump(self,xyz,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36):
    
        deltaP_pp_star,f_pp_star,r_pp_star,rho_su_pp_star = xyz[0],xyz[1],xyz[2],xyz[3]
    
        eta_is = c0 + c1*deltaP_pp_star + c2*deltaP_pp_star**2 + c3*deltaP_pp_star**3 + c4*f_pp_star + c5*f_pp_star**2 + c6*f_pp_star**3 + c7*r_pp_star + c8*r_pp_star**2 + c9*r_pp_star**3 + c10*rho_su_pp_star + c11*rho_su_pp_star**2 + c12*rho_su_pp_star**3 + c13*deltaP_pp_star*f_pp_star + c14*deltaP_pp_star*f_pp_star**2 + c15*deltaP_pp_star*r_pp_star + c16*deltaP_pp_star*r_pp_star**2 + c17*deltaP_pp_star*rho_su_pp_star + c18*deltaP_pp_star*rho_su_pp_star**2 + c19*deltaP_pp_star**2*f_pp_star + c20*deltaP_pp_star**2*f_pp_star**2 + c21*deltaP_pp_star**2*r_pp_star + c22*deltaP_pp_star**2*r_pp_star**2 + c23*deltaP_pp_star**2*rho_su_pp_star + c24*deltaP_pp_star**2*rho_su_pp_star**2 + c25*f_pp_star*r_pp_star + c26*f_pp_star*r_pp_star**2 + c27*f_pp_star*rho_su_pp_star + c28*f_pp_star*rho_su_pp_star**2 + c29*f_pp_star*2*r_pp_star + c30*f_pp_star**2*r_pp_star**2 + c31*f_pp_star**2*rho_su_pp_star + c32*f_pp_star**2*rho_su_pp_star**2 + c33*r_pp_star*rho_su_pp_star + c34*r_pp_star*rho_su_pp_star**2 + c35*r_pp_star**2*rho_su_pp_star + c36*r_pp_star**2*rho_su_pp_star**2

        return eta_is
    
    
    
    
    def Calculate(self):
        
        #Local copies of coefficients
        pump_iseff_params=self.eta        
        pump_mdot_params=self.M
        
        #Pump ref parameters
        deltaP_pp_ref = 1000   #kPa
        r_pp_ref = 9
        rho_pp_ref = 1000#Props('D','P',120,'T',15+273.15,'R245FA') #kg/s
        f_pp_ref = 30 #Hz
        
        #Inlet state
        hin=PropsSI('H','T',self.Tin+273.15,'P',self.pin_r*1000+100,self.Ref)/1000
        self.s_in = PropsSI('S','T',self.Tin+273.15,'P',self.pin_r*1000+100,self.Ref)/1000
        
        
        #Compute the pressure difference between the outlet and the inlet of the pump
        deltaP_pp = (self.pout_r-self.pin_r)  #kPa
        deltaP_pp_star = (deltaP_pp-deltaP_pp_ref)/deltaP_pp_ref
        f_pp = self.f_pp
        f_pp_star= (f_pp- f_pp_ref)/f_pp_ref
        
        r_pp = self.pout_r/self.pin_r
        r_pp_star = (r_pp-r_pp_ref)/r_pp_ref
        
        rho_su_pp = PropsSI('D','P',self.pin_r*1000,'T',self.Tin+273.15,self.Ref)
        rho_su_pp_star =  (rho_su_pp-rho_pp_ref)/rho_pp_ref
    
    
        #Mass flow rate
        xy = np.array([f_pp,deltaP_pp])
        self.m_dot = self.MassFlowPump(xy,*pump_mdot_params)

        #Efficiency
        xyz = np.array([deltaP_pp_star,f_pp_star,r_pp_star,rho_su_pp_star])
        self.eta_is_pp = self.IsEffPump(xyz,*pump_iseff_params) 
        #print self.eta_is_pp
        #self.eta_is_pp = m_dot*(h_ex_is_pp-h_su_pp)/(W_dot_pp/1000)
        
        #Outlet state
        self.Tout_s = PropsSI('T','S',self.s_in*1000,'P',self.pout_r*1000+100,self.Ref)
        self.hout_s = PropsSI('H','S',self.s_in*1000,'P',self.pout_r*1000+100,self.Ref)/1000
        #print self.hout_s
        #print hin
        self.hout= hin + (self.hout_s-hin)/self.eta_is_pp   #[kJ/kg]

        self.Tout=PropsSI('T','H',self.hout*1000,'P',self.pout_r*1000+100,self.Ref) #in K
        self.s_out=PropsSI('S','T',self.Tout,'P',self.pout_r*1000+100,self.Ref)/1000   
         
 
        
        #Get the estimation of the power
        self.W_dot=self.m_dot*(self.hout - hin)*1000
        




if __name__=='__main__':
    #Run the code as example
    pin_r_list=[125.006]
    pout_r_list=[713.586]
    fpp_list=[28.16]
    Tin_list=[17.22]
    W_dot_pp_list = [598.3]
    eta_pp_list = [0.117106158]       
    for pin_r,pout_r,fpp,Tin,Wdot,eta_pp in zip(pin_r_list,pout_r_list,fpp_list,Tin_list,W_dot_pp_list,eta_pp_list):
        kwds={
                'eta':[-1.94559985e+03,2.61004607e+03,-2.24868329e+02,-1.31580710e+03,
  -3.41646308e+03,-6.16059114e+02,1.70061718e+03,-6.12847688e+03,-3.08159095e+03,-2.13581139e+02,1.12347642e+04,  -1.91028124e+04,7.30201901e+03,2.57633083e+03,-5.00775452e+03,2.25464045e+02,-3.98897093e+02,-1.43095042e+04,   1.71628057e+04,4.38898501e+03,-3.82497930e+00,3.71939847e+02,1.77508114e+02,-5.73292533e+03,8.45941415e+03, 1.69904379e+03,5.55650945e+02,1.99756918e+04,-2.62485094e+04,-7.59467705e+02,-2.56645354e+02, -2.15969687e+01,-5.61487624e+03,2.84386468e+04,-3.23592011e+04,8.93600159e+03,-4.52678170e+03], 
                'M':[3.05754838e-01,-1.37732305e-03,-4.26385445e-07,-2.68106448e-02,1.98497578e-03],
                'Ref': 'R245fa',
                'pin_r': pin_r,
                'pout_r':pout_r,
                'f_pp':fpp,
                'Tin': Tin
              }
        Pump=PumpCentrifugalClass(**kwds)
        Pump.Calculate()
        
        print 'Measured:',Pump.W_dot,'W','Calculated:',Wdot,'W'
        print Pump.m_dot,'kg/s'
        print Pump.Tout, 'K'
        print 'Measured:',Pump.eta_is_pp, '-','Calculated:',eta_pp,'-'