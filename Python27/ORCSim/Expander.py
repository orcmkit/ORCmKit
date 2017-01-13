from __future__ import division #Make integer 3/2 give 1.5 in python 2.x
from CoolProp.CoolProp import PropsSI
from math import tan,atan,pi,sin,log
import numpy as np
import math
import time
import warnings

warnings.simplefilter("ignore",RuntimeWarning)
from random import randint, random




class ExpanderClass():


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

        return [
                ('ff_exp','-',self.phi),
                ('eta_exp_ad','-',self.eta_s),
                ('mdot_exp','kg/s',self.m_dot),
                ('power_exp','W',self.W_dot),
                ('torque_exp','N-m',self.tau),
                ('P_exp_su','kPa',self.P_su),
                ('T_exp_ex','C',self.T_ex),
                ('P_exp_ex','kPa',self.P_ex),
                ('v_ratio_exp','-',self.r_v)
            ]


    def Calculate(self):

        #We can provide either supply temperature or supply enthalpy to fix the inlet state of the expander
        if 'h' in self.inputs: #the preferred input is enthalpy
            self.v_su = 1/PropsSI('D','H',self.h_su,'P',self.P_su*1000,self.Ref) #expander inlet state properties
        elif 'T' in self.inputs:
            self.v_su = 1/PropsSI('D','T',self.T_su+273.15,'P',self.P_su*1000,self.Ref) #expander inlet state properties
        else:
            raise ValueError('Supply temperature or enthalpy must be specified for this Expander object before the Calculate method can execute.')

        rho_su = 1/self.v_su
        
        "Filling factor"
        self.r_p=self.P_su/self.P_ex

        #compute the fill factor
        p_star = (self.P_su - 1000)/1000
        r_star_p = (self.r_p - self.params.rp_ref)/self.params.rp_ref
        rho_star = (rho_su - 1000)/1000
        self.phi = self.params.k_1 + self.params.k_2*log(self.N/self.params.N_FF_ref) + self.params.k_3*r_star_p + self.params.k_4*p_star + self.params.k_5*rho_star
        #compute the mass flow rate from the fill factor
        V_disp = self.V_disp_compressor/self.r_v_built
        
        m_dot_ideal = V_disp*self.N/60/self.v_su/1e6 #60 s/m and 1e6 cm3/m3
        self.m_dot = self.phi*m_dot_ideal

        "Pacejka's Equation"
        N_star_rot = (self.N - self.params.N_ref)/self.params.N_ref
        p_star = (self.P_su - self.params.P_ref)/self.params.P_ref
        N_star_rot_n = (self.params.N_rot_n - self.params.N_ref)/self.params.N_ref
        r_p_0 = self.params.r_p_0_n + self.params.a_0*N_star_rot
        delta = self.params.delta_n + self.params.a_1*p_star + self.params.a_2*N_star_rot
        r_p_max = self.params.r_p_max_n + self.params.a_3*p_star + self.params.a_4*N_star_rot
        y_max = self.params.y_max_n + self.params.a_5*p_star + self.params.a_6*(N_star_rot - N_star_rot_n)**2
        B = delta/(self.params.xi*y_max)
        E = (B*(r_p_max - r_p_0) - np.tan(pi/(2*self.params.xi)))/(B*(r_p_max - r_p_0) - np.arctan(B*(r_p_max - r_p_0)))
        self.eta_s = y_max*np.sin(self.params.xi*np.arctan(B*(self.r_p - r_p_0) - E*(B*(self.r_p - r_p_0) - np.arctan(B*(self.r_p - r_p_0)))))

        #We can provide either supply temperature or supply enthalpy to fix the inlet state of the expander
        if 'h' in self.inputs: #the preferred input is enthalpy
            self.T_su = PropsSI('T','H',self.h_su,'P',self.P_su*1000,self.Ref) - 273.15
            self.s_su = PropsSI('S','H',self.h_su,'P',self.P_su*1000,self.Ref)
        elif 'T' in self.inputs:
            self.h_su = PropsSI('H','T',self.T_su+273.15,'P',self.P_su*1000,self.Ref)
            self.s_su = PropsSI('S','T',self.T_su+273.15,'P',self.P_su*1000,self.Ref)
        else:
            raise ValueError('Supply temperature or enthalpy must be specified for this Expander object before the Calculate method can execute.')

        #adiabatic reversible expander exit state
        h_ex_s = PropsSI('H','P',self.P_ex*1000,'S',self.s_su,self.Ref)
        self.T_ex_s = PropsSI('T','P',self.P_ex*1000,'H',h_ex_s,self.Ref)
        #Expander power and torque
        self.W_dot = self.eta_s*self.m_dot*(self.h_su - h_ex_s)
        self.tau = self.W_dot/(self.N*2*pi)*60*1000 #2*pi rad/rev, 60 s/min, and 1000 N/kN
        #Expander adiabatic exit state
        h_ex_ad = self.h_su - self.W_dot/self.m_dot
        T_ex_ad = PropsSI('T','P',self.P_ex*1000,'H',h_ex_ad,self.Ref) - 273.15
        v_ex_ad = 1/PropsSI('D','P',self.P_ex*1000,'H',h_ex_ad,self.Ref)
        # s_ex_ad = Props('S','P',self.P_ex,'H',h_ex_ad/1000,self.Ref)
        r_v_ad = v_ex_ad/self.v_su
        #record exit state dependent values as those for the adiabatic case for now
        self.T_ex = T_ex_ad +273.15
        self.h_ex = h_ex_ad
        #print 'expander:',self.h_ex
        self.r_v = r_v_ad
        
        
        self.hout = h_ex_ad
        self.s_out = PropsSI('S','T',self.T_ex,'P',self.P_ex*1000+100,self.Ref)

        #print 'etais:',self.eta_s
        #print 'rp:',self.r_p

if __name__=='__main__':

    def Kortrijk_SSE():

        class struct(): pass
        exp_params = struct()

        exp_params.a_0 =-1.93529023e-03# -2.11131627e-03
        exp_params.a_1 = -1.25256170e-02#3.20637318e-01
        exp_params.a_2 = -1.17292542e-01#8.68498341e-01
        exp_params.a_3 = -4.73874764e-02#4.96963978e+00
        exp_params.a_4 = 7.80410436e+00#-1.23773707e-01
        exp_params.a_5 = 6.36530295e-04#-9.11764461e-02
        exp_params.a_6 = 5.78187315e-01#5.75017501e-01
        exp_params.xi = 1.23877317e+00#1.22425573e+00
        exp_params.r_p_0_n = 3.07600000e+00
        exp_params.delta_n =  7.08651674e-01#3.94590727e-04   
        exp_params.r_p_max_n = 5.99558088e+00    
        exp_params.y_max_n = 5.19000000e-01 
        exp_params.N_rot_n = 3547
        exp_params.N_ref = 3500
        exp_params.N_FF_ref = 3000
        exp_params.P_ref = 1000
        exp_params.rp_ref = 6
        exp_params.k_1 = 1.489
        exp_params.k_2 = 0.3427
        exp_params.k_3 = 0.8187 
        exp_params.k_4 = 0.1435
        exp_params.k_5 = 0.0 
        N = 3000
        P_su = 689
        P_ex = 180
        T_su = 100
        fluid = 'R245fa'
        r_v_builtin = 4.8
        V_disp_compressor = 2*6*57.39 #[cm^3/rev]
    
        kwds={
                'params':exp_params,
                'N':N,
                'P_su':P_su,
                'P_ex':P_ex,
                'T_su':T_su,
                'r_v_built':r_v_builtin,
                'V_disp_compressor':V_disp_compressor,
                'inputs':'Psu_Pex_T',
                'Ref':fluid
                }
        Exp = ExpanderClass(**kwds)
        
        Exp.Calculate()

        print 'eta_s [-]:', Exp.eta_s  
        print 'm_dot [kg/s]:',Exp.m_dot      
        
        
        
        
    #Calculate
    Kortrijk_SSE()