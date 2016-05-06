"""
ASME ORC 2015
3rd International Seminar on ORC Systems
12-14 October 2015, Brussels, Belgium


Moving-boundary PHEX model derived from ACHP (https://github.com/TSTK) and modified 
to include incompressible fluids such as thermal oils and lubricants.
At the moment only INCOMP::T66 is included.


Reference Paper:

D. Ziviani et al. "ORCSIM: a generalized organic Rankine Cycle Simulation Tool". Paper#30

     
For illustrative purposes only!
"""


from __future__ import division

from CoolProp.CoolProp import PropsSI 

from Correlations import ShahEvaporation_Average,PHE_1phase_hdP,Cooper_PoolBoiling,TwoPhaseDensity,TrhoPhase_ph,Phase_ph,LMPressureGradientAvg,KandlikarPHE,Bertsch_MC,AccelPressureDrop,ShahCondensation_Average,LongoCondensation

from math import pi,exp,log,sqrt,tan,cos,sin
from scipy.optimize import brentq
import numpy as np
import pylab
from matplotlib.pyplot import plot, show, figure, semilogy, xlim, ylim, title, xlabel, ylabel, legend


def Tsat(FluidName, p, Q, T_guess):
    if Q<0 or Q>1:
        raise ValueError('Quality must be a value between 0.0 and 1.0')
    return PropsSI('T','P',p*1000,'Q',Q,FluidName)

def IsFluidType(FluidName, Type):
    if Type=='Brine':   #Only CoolProp handles brine so it cannot be possible for the fluid type to be 'Brine' here
        return 0
    else:
        return 1        #I don't have any routines which set the fluid type so I can't really check for any other fluid types


#TODO: Modified for Incompressible as hot fluid
def IsFluidIncompType(FluidName, Type):
    if Type=='INCOMP::T66':
        return 1
    else:
        pass

class PHEHXClass():
    """
    There are a number of possibilities:
        
        Each fluid can:
        a) Not change phase
        c) Evaporate 
        c) Condense
        
        Possibility matrix
        
                                      Hot stream
         Cold Stream  || Subcooled ||  Two-Phase || Superheated ||
                      --------------------------------------------
         Subcooled    ||           ||            ||             ||
                      --------------------------------------------
         Two-Phase    ||           ||            ||             ||
                      --------------------------------------------
         Superheated  ||           ||            ||             ||
                      --------------------------------------------
    
    Hot stream goes to the left in the matrix, cold stream goes down.  If 
    hot stream comes in subcooled, there are only three combinations that 
    might exist.  
    
    Based on inlet states can figure out what states are possible.  
    """
    
    def __init__(self,**kwargs):
        #Load the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
        
    def Update(self,**kwargs):
        #Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    
    def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items:
                [0] Description of value
                [1] Units of value
                [2] The value itself
        """
        return [
            ('Effective Length','m',self.Lp),
            #('Outlet Superheat','K',self.Tout_c-self.Tdew_c),
            ('Total heat transfer area','m2',self.A_c_wetted + self.A_h_wetted),
            ('Q Total','W',self.Q),
            ('Q Superheat Hot','W',self.Q_superheated_h),
            ('Q Two-Phase Hot','W',self.Q_2phase_h),
            ('Q Subcooled Hot','W',self.Q_subcooled_h),
            ('Q Superheat Cold','W',self.Q_superheated_c),
            ('Q Two-Phase Cold','W',self.Q_2phase_c),
            ('Q Subcooled Cold','W',self.Q_subcooled_c),
            ('Inlet hot stream temp','K',self.Tin_h),
            ('Outlet hot stream temp','K',self.Tout_h),
            ('Inlet cold stream temp','K',self.Tin_c),
            ('Outlet cold stream temp','K',self.Tout_c),
            ('Charge Total','kg',self.Charge_h),
            ('Charge Superheat','kg',self.Charge_superheated_h),
            ('Charge Two-Phase','kg',self.Charge_2phase_h),
            ('Charge Subcool','kg',self.Charge_subcooled_h),
            ('Charge Total','kg',self.Charge_c),
            ('Charge Superheat','kg',self.Charge_superheated_c),
            ('Charge Two-Phase','kg',self.Charge_2phase_c),
            ('Charge Subcool','kg',self.Charge_subcooled_c),
            ('Hot HTC Superheat','W/m^2-K',self.h_superheated_h),
            ('Hot HTC Two-Phase','W/m^2-K',self.h_2phase_h),
            ('Hot HTC Subcool','W/m^2-K',self.h_subcooled_h),
            ('Cold Mean HTC Superheat','W/m^2-K',self.h_superheated_c),
            ('Cold Mean HTC Ref. Two-Phase','W/m^2-K',self.h_2phase_c),
            ('Cold Mean HTC Ref. Subcool','W/m^2-K',self.h_subcooled_c),
            ('Pressure Drop Hot','Pa',self.DP_h),
            ('Pressure Drop Hot superheated','Pa',self.DP_superheated_h),
            ('Pressure Drop Hot 2 phase','Pa',self.DP_2phase_h),
            ('Pressure Drop Hot subcooled','Pa',self.DP_subcooled_h),
            ('Pressure Drop Cold','Pa',self.DP_c),
            ('Pressure Drop Cold superheated','Pa',self.DP_superheated_c),
            ('Pressure Drop Cold 2 phase','Pa',self.DP_2phase_c),
            ('Pressure Drop Cold subcooled','Pa',self.DP_subcooled_c),
            ('Area fraction Superheat Hot','-',self.w_superheated_h),
            ('Area fraction Two-Phase Hot','-',self.w_2phase_h),
            ('Area fraction Subcooled Hot','-',self.w_subcooled_h),
            ('Area fraction Superheat Cold','-',self.w_superheated_c),
            ('Area fraction Two-Phase Cold','-',self.w_2phase_c),
            ('Area fraction Subcooled Cold','-',self.w_subcooled_c)
         ]
        

    
    def DetermineHTBounds(self):
        # See if each phase could change phase if it were to reach the
        # inlet temperature of the opposite phase 
        
        #Inlet phases
        #TODO: Modified for Incompressible as hot fluid
        if self.Ref_h =='INCOMP::T66':
            self.Tin_h,rhoin_h,Phasein_h=TrhoPhase_ph(self.Ref_h,self.pin_h,self.hin_h,None,None,None,None)
        else:
            self.Tin_h,rhoin_h,Phasein_h=TrhoPhase_ph(self.Ref_h,self.pin_h,self.hin_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)
        
        
        if self.Ref_c == 'INCOMP::MEG-30%':
            self.Tin_c,rhoin_c,Phasein_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hin_c,None,None,None,None) 
        elif self.Ref_c == 'POE':
            self.Tin_c,rhoin_c,Phasein_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hin_c,None,None,None,None) 
        elif self.Ref_c == 'ACD100FY':
            self.Tin_c,rhoin_c,Phasein_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hin_c,None,None,None,None)
        else: 
            self.Tin_c,rhoin_c,Phasein_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hin_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)
            #print self.hin_c
        # Find the maximum possible rate of heat transfer as the minimum of 
        # taking each stream to the inlet temperature of the other stream
        
        hout_h=PropsSI('H','T',self.Tin_c,'P',self.pin_h*1000,self.Ref_h)#*1000
        #print 'hout_h:',hout_h,
        
        if self.Ref_c == 'POE':
            hout_c = h_l(self.Ref_c,self.Tin_h,self.pin_c)
        elif self.Ref_c == 'ACD100FY':
            hout_c = h_l(self.Ref_c,self.Tin_h,self.pin_c)
        else:
            hout_c=PropsSI('H','T',self.Tin_h,'P',self.pin_c*1000,self.Ref_c)#*1000
            #print 'hout_c', hout_c
        if self.mdot_h==0:  #This is used for calculations with steam as the hot fluid when the mass flow rate is unknown.  It assigns a mass flow rate for the steam that is more than enough to get the cold fluid to reach the steam temperature without condensing all the steam.
            Qmax=self.mdot_c*(hout_c-self.hin_c)
            hin_h_sat=PropsSI('H','T',self.Tin_h,'Q',1.0,self.Ref_h)#*1000
            hout_h=PropsSI('H','T',self.Tin_h,'Q',0.0,self.Ref_h)#*1000
            DELTAh_cond=hin_h_sat-hout_h
            self.mdot_h=1.3*Qmax/(DELTAh_cond)
        else:
            Qmax=min([self.mdot_c*(hout_c-self.hin_c),self.mdot_h*(self.hin_h-hout_h)])
            #print 'Qmax:',Qmax
        
#        Qmax=min([self.mdot_c*(hout_c-self.hin_c),self.mdot_h*(self.hin_h-hout_h)])
        
        # Now we need to check for internal pinch points where the temperature
        # profiles would tend to overlap given the "normal" definitions of 
        # maximum heat transfer of taking each stream to the inlet temperature 
        # of the other stream
        #
        # First we build the same vectors of enthalpies like below
        EnthalpyList_c,EnthalpyList_h=self.BuildEnthalpyLists(Qmax)
        # Then we find the temperature of each stream at each junction
        TList_c=np.zeros_like(EnthalpyList_c)
        TList_h=np.zeros_like(EnthalpyList_h)

        if not len(EnthalpyList_h)==len(EnthalpyList_c):
            raise ValueError('Length of enthalpy lists for both fluids must be the same')
        
        #Make the lists of temperatures of each fluid at each cell boundary
        #TODO: Modified for Incompressible as hot/cold fluid
        for i in range(len(EnthalpyList_h)):
            if self.Ref_c =='INCOMP::MEG-30%':
                TList_c[i] = TrhoPhase_ph(self.Ref_c,self.pin_c,EnthalpyList_c[i],None,None,None,None)[0]            
            elif self.Ref_c =='POE':
                TList_c[i] = TrhoPhase_ph(self.Ref_c,self.pin_c,EnthalpyList_c[i],None,None,None,None)[0] 
            elif self.Ref_c == 'ACD100FY':
                TList_c[i] = TrhoPhase_ph(self.Ref_c,self.pin_c,EnthalpyList_c[i],None,None,None,None)[0] 
            else:
                TList_c[i] = TrhoPhase_ph(self.Ref_c,self.pin_c,EnthalpyList_c[i],self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0]
            if self.Ref_h =='INCOMP::T66':
                TList_h[i] = TrhoPhase_ph(self.Ref_h,self.pin_h,EnthalpyList_h[i],None,None,None,None)[0]
            else:
                TList_h[i] = TrhoPhase_ph(self.Ref_h,self.pin_h,EnthalpyList_h[i],self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0]
                
#         figure()
#         pseudolength = (np.array(EnthalpyList_c) - EnthalpyList_c[0])/(EnthalpyList_c[-1] - EnthalpyList_c[0])
#         plot(pseudolength, TList_c, label=str.split(self.Ref_c,'.')[0]) 
#         plot(pseudolength, TList_h, label=str.split(self.Ref_h,'.')[0])
#         xlabel('Normalized Enthalpy (a pseudolength) [-]')
#         ylabel('Fluid Temperature [K]')
#         title('Theoretical Maximum Heat Transfer')
#         leg=legend(loc='best', fancybox=True)
#         leg.get_frame().set_alpha(0.5)   

#        #Double-check that the edges are not pinched
#        if TList_c[0]-1e-9>TList_h[0] or TList_c[-1]-1e-9>TList_h[-1]:
#            raise ValueError('Outlet or inlet of PHE is pinching.  Why?')
        
        #TODO: could do with more generality if both streams can change phase
        #Check if any internal points are pinched
#        if np.sum(TList_c>TList_h)>0:
#            #Loop over the internal cell boundaries
#            for i in range(1,len(TList_c)-1):
#                #If cold stream is hotter than the hot stream
#                if TList_c[i]-1e-9>TList_h[i]:
#                    #Find new enthalpy of cold stream at the hot stream cell boundary
#                    hpinch=Props('H','T',TList_h[i],'P',self.pin_c,self.Ref_c)*1000
#                    #Find heat transfer of hot stream in right-most cell
#                    Qextra=self.mdot_h*(EnthalpyList_h[i+1]-EnthalpyList_h[i])
#                    Qmax=self.mdot_c*(hpinch-self.hin_c)+Qextra
          
        #Brandon's version for general derating of Qmax (handles hot fluid, cold fluid, or dual pinch points)
        #Replaces above lines up to "if np.sum(TList_c..."
        Qmax_new = None
        
        #TODO: Modified for Incompressible as hot fluid
        if self.Ref_h =='INCOMP::T66':
            
            if self.Ref_c =='ACD100FY':
                
                pass
                
            else:    
                if TList_c[0]+1e-6 < self.Tbubble_c <= TList_c[-1]+1e-6:  
                #then there is a cold fluid phase change boundary at the bubble point 
                #that could result in a hot fluid pinch point 
                    i = 0
                    while TList_c[i] < self.Tbubble_c-1e-9:  #find the boundary where the cold fluid changes phase
                        i+=1
                    if (self.Tbubble_c > TList_h[i]):  #hot fluid pinch point
                        hpinch = PropsSI('H','T',self.Tbubble_c,'P',self.pin_h*1000,self.Ref_h)#*1000
                        Qmax_new = self.mdot_h*(self.hin_h - hpinch) + self.mdot_c*(self.hsatL_c - self.hin_c)
        
        
        
        #TODO: Modified for Incompressible as cold fluid
        elif self.Ref_c =='INCOMP::MEG-30%':
                
            if TList_h[0]-1e-6 <= self.Tdew_h < TList_h[-1]-1e-6:  #then there is a hot fluid phase change boundary at the dew point that could result in a cold fluid pinch point
                i = len(TList_h)-1
                while TList_h[i] > self.Tdew_h+1e-9:
                    i-=1
                if (self.Tdew_h < TList_c[i]): #cold fluid pinch point          
                    hpinch = PropsSI('H','T',self.Tdew_h,'P',self.pin_c,self.Ref_c)*1000
                    Qmax_new = self.mdot_c*(hpinch - self.hin_c) + self.mdot_h*(self.hin_h - self.hsatV_h)
        
        #TODO: Modified to have lubricants as cold fluid; missing thermal oil / lubrican case
        elif self.Ref_c =='POE': 
            if TList_h[0]-1e-6 <= self.Tdew_h < TList_h[-1]-1e-6:  #then there is a hot fluid phase change boundary at the dew point that could result in a cold fluid pinch point
                i = len(TList_h)-1
                while TList_h[i] > self.Tdew_h+1e-9:
                    i-=1
                if (self.Tdew_h < TList_c[i]): #cold fluid pinch point          
                    hpinch = h_l(self.Ref_c,self.Tdew_h,self.pin_c)
                    Qmax_new = self.mdot_c*(hpinch - self.hin_c) + self.mdot_h*(self.hin_h - self.hsatV_h)            
        
        elif self.Ref_c == 'ACD100FY':
            if TList_h[0]-1e-6 <= self.Tdew_h < TList_h[-1]-1e-6:  #then there is a hot fluid phase change boundary at the dew point that could result in a cold fluid pinch point
                i = len(TList_h)-1
                while TList_h[i] > self.Tdew_h+1e-9:
                    i-=1
                if (self.Tdew_h < TList_c[i]): #cold fluid pinch point          
                    hpinch = h_l(self.Ref_c,self.Tdew_h,self.pin_c)
                    Qmax_new = self.mdot_c*(hpinch - self.hin_c) + self.mdot_h*(self.hin_h - self.hsatV_h)         
       
        
        else:    
        
            if (TList_c[0]+1e-6 < self.Tbubble_c <= TList_c[-1] and TList_h[0] <= self.Tdew_h < TList_h[-1]-1e-6) \
            and ((self.Tbubble_c > TList_h[1]) and (self.Tdew_h < TList_c[-2])): #both fluids change phase AND dual pinch points
                hpinch_h = PropsSI('H','T',self.Tbubble_c,'P',self.pin_h*1000,self.Ref_h)#*1000
                hpinch_c = PropsSI('H','T',self.Tdew_h,'P',self.pin_c*1000,self.Ref_c)#*1000
        #            Qmax_h = self.mdot_h*(self.hsatV_h - hpinch_h)
        #            Qmax_c = self.mdot_c*(hpinch_c - self.hsatL_c)
                Qmax_new = min(self.mdot_h*(self.hsatV_h - hpinch_h), self.mdot_c*(hpinch_c - self.hsatL_c)) + self.mdot_c*(self.hsatL_c - self.hin_c) + self.mdot_h*(self.hin_h - self.hsatV_h)
        #            Qmax_new_check = min(self.mdot_h*(self.hin_h - hpinch_h) + self.mdot_c*(self.hsatL_c - self.hin_c), 
        #                                 self.mdot_c*(hpinch_c - self.hin_c) + self.mdot_h*(self.hin_h - self.hsatV_h))
        #            
        #            print 'Qmax_new_check =', Qmax_new_check
                
            else: #I had to add the 1e-6 tolerance because otherwise the test could return false.  Checking quality would be better maybe?
                if TList_c[0]+1e-6 < self.Tbubble_c <= TList_c[-1]+1e-6:  #then there is a cold fluid phase change boundary at the bubble point that could result in a hot fluid pinch point 
                    i = 0
                    while TList_c[i] < self.Tbubble_c-1e-9:  #find the boundary where the cold fluid changes phase
                        i+=1
                    if (self.Tbubble_c > TList_h[i]):  #hot fluid pinch point
                        hpinch = PropsSI('H','T',self.Tbubble_c,'P',self.pin_h*1000,self.Ref_h)#*1000
                        Qmax_new = self.mdot_h*(self.hin_h - hpinch) + self.mdot_c*(self.hsatL_c - self.hin_c)
                    
                if TList_h[0]-1e-6 <= self.Tdew_h < TList_h[-1]-1e-6:  #then there is a hot fluid phase change boundary at the dew point that could result in a cold fluid pinch point
                    i = len(TList_h)-1
                    while TList_h[i] > self.Tdew_h+1e-9:
                        i-=1
                    if (self.Tdew_h < TList_c[i]): #cold fluid pinch point          
                        hpinch = PropsSI('H','T',self.Tdew_h,'P',self.pin_c*1000,self.Ref_c)#*1000
                        Qmax_new = self.mdot_c*(hpinch - self.hin_c) + self.mdot_h*(self.hin_h - self.hsatV_h)
            
#        print 'Qmax =', Qmax
#        print 'Qmax_new =', Qmax_new
        if Qmax_new:
            assert(abs(Qmax_new) <= abs(Qmax))
            Qmax = Qmax_new 
#            EnthalpyList_c,EnthalpyList_h=self.BuildEnthalpyLists(Qmax)  #Rebuild the enthalpy list so we can plot it
#            #Rebuild the temperature list
#            TList_c=np.zeros_like(EnthalpyList_c)
#            TList_h=np.zeros_like(EnthalpyList_h)
#            for i in range(len(EnthalpyList_h)):
#                TList_c[i] = TrhoPhase_ph(self.Ref_c,self.pin_c,EnthalpyList_c[i],self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0]
#                TList_h[i] = TrhoPhase_ph(self.Ref_h,self.pin_h,EnthalpyList_h[i],self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0]
#         
#        #Here we will plot a qualitative look at what the temperature profiles would be for maximum heat transfer
#        #The abscissa is dimensionless enthalpy so it doesn't show us how long each section is in the physical heat exchanger
#        figure()
#        pseudolength = (np.array(EnthalpyList_c) - EnthalpyList_c[0])/(EnthalpyList_c[-1] - EnthalpyList_c[0])
#        plot(pseudolength, TList_c, label=str.split(self.Ref_c,'.')[0]) 
#        plot(pseudolength, TList_h, label=str.split(self.Ref_h,'.')[0])
#        xlabel('Normalized Enthalpy (a pseudolength) [-]')
#        ylabel('Fluid Temperature [K]')
#        title('Theoretical Maximum Heat Transfer (Pinch Points Checked)')
#        leg=legend(loc='best', fancybox=True)
#        leg.get_frame().set_alpha(0.5)   
#        show()
#        end Brandon's version for general derating of Qmax
        
        return Qmax
        
    def PlateHTDP(self,Ref,T,p,mdot_gap):
        """
        For single phase fluids, inputs in K, kPa, outputs in W/m^2-K, J/kg-K
        """
        Inputs={
            'Ref':Ref,
            'T':T,
            'p':p,
            'mdot_gap' : mdot_gap,
            'PlateAmplitude': self.PlateAmplitude,
            'PlateWavelength' : self.PlateWavelength,
            'InclinationAngle': self.InclinationAngle,
            'Bp': self.Bp,
            'Lp': self.Lp
        }
        Outputs=PHE_1phase_hdP(Inputs)
        return Outputs['h'],Outputs['cp'],Outputs
    
    def BuildEnthalpyLists(self,Q):
        #Start the enthalpy lists with inlet and outlet enthalpies
        #Ordered from lowest to highest enthalpies for both streams
        EnthalpyList_h=[self.hin_h-Q/self.mdot_h, self.hin_h]
        EnthalpyList_c=[self.hin_c,self.hin_c+Q/self.mdot_c]
        
        #Save the value of Q and outlet enthalpies
        self.Q=Q
        self.hout_h=EnthalpyList_h[0]
        self.hout_c=EnthalpyList_c[1]
        
        #Find the phase boundaries that exist, and add them to lists
        if IsFluidType(self.Ref_h,'Brine'):
            hsatL_h=1e9
            hsatV_h=1e9
        #TODO: Modified for Incompressible as hot fluid
        elif self.Ref_h =='INCOMP::T66':
            hsatL_h=1e9
            hsatV_h=1e9            
            
        else:
            hsatL_h=PropsSI('H','T',self.Tbubble_h,'D',self.rhosatL_h,self.Ref_h)#*1000
            hsatV_h=PropsSI('H','T',self.Tdew_h,'D',self.rhosatV_h,self.Ref_h)#*1000
        
        if IsFluidType(self.Ref_c,'Brine'):
            hsatL_c=1e9
            hsatV_c=1e9
        #TODO: Modified for Incompressible as cold fluid
        elif self.Ref_c =='INCOMP::MEG-30%':
            hsatL_c=1e9
            hsatV_c=1e9
        elif self.Ref_c =='POE':
            hsatL_c=1e9
            hsatV_c=1e9        
        elif self.Ref_c == 'ACD100FY':
            hsatL_c=1e9
            hsatV_c=1e9
            
        else:
            hsatL_c=PropsSI('H','T',self.Tbubble_c,'D',self.rhosatL_c,self.Ref_c)#*1000
            hsatV_c=PropsSI('H','T',self.Tdew_c,'D',self.rhosatV_c,self.Ref_c)#*1000
        
        # Check whether the enthalpy boundaries are within the bounds set by 
        # the imposed amount of heat transfer
        if hsatV_c<EnthalpyList_c[-1] and hsatV_c>EnthalpyList_c[0]:
            EnthalpyList_c.insert(len(EnthalpyList_c)-1,hsatV_c)
        if hsatL_c<EnthalpyList_c[-1] and hsatL_c>EnthalpyList_c[0]:
            EnthalpyList_c.insert(1,hsatL_c)
            
        if hsatV_h<EnthalpyList_h[-1] and hsatV_h>EnthalpyList_h[0]:
            EnthalpyList_h.insert(len(EnthalpyList_h)-1,hsatV_h)
        if hsatL_h<EnthalpyList_h[-1] and hsatL_h>EnthalpyList_h[0]:
            EnthalpyList_h.insert(1,hsatL_h)
            
        I_h=0
        I_c=0
        while I_h<len(EnthalpyList_h)-1:
            #Try to figure out whether the next phase transition is on the hot or cold side     
            Qbound_h=self.mdot_h*(EnthalpyList_h[I_h+1]-EnthalpyList_h[I_h])
            Qbound_c=self.mdot_c*(EnthalpyList_c[I_c+1]-EnthalpyList_c[I_c])
            if Qbound_h<Qbound_c-1e-9:
                # Minimum amount of heat transfer is on the hot side,
                # add another entry to EnthalpyList_c 
                EnthalpyList_c.insert(I_c+1, EnthalpyList_c[I_c]+Qbound_h/self.mdot_c)
            elif Qbound_h>Qbound_c+1e-9:
                # Minimum amount of heat transfer is on the cold side,
                # add another entry to EnthalpyList_h at the interface
                EnthalpyList_h.insert(I_h+1, EnthalpyList_h[I_h]+Qbound_c/self.mdot_h)
            I_h+=1
            I_c+=1
                    
        self.hsatL_c=hsatL_c
        self.hsatL_h=hsatL_h
        self.hsatV_c=hsatV_c
        self.hsatV_h=hsatV_h

        return EnthalpyList_c,EnthalpyList_h
    
    def PostProcess(self,cellList):
        """
        Combine all the cells to calculate overall parameters like pressure drop
        and fraction of heat exchanger in two-phase on both sides
        """
        def collect(cellList,tag,tagvalue,out):
            collectList=[]
            for cell in cellList:
                if cell[tag]==tagvalue:
                    collectList.append(cell[out])
            return collectList
        self.DP_c=0
        self.DP_c_superheat=0
        self.DP_c_2phase=0
        self.DP_c_subcooled=0
        self.DP_h=0
        self.DP_h_superheat=0
        self.DP_h_2phase=0
        self.DP_h_subcooled=0
        self.Charge_c=0
        self.Charge_h=0
        for cell in cellList:
            self.DP_c+=cell['DP_c']
            self.DP_h+=cell['DP_h']
            self.Charge_c+=cell['Charge_c']
            self.Charge_h+=cell['Charge_h']
        self.w_superheated_h=sum(collect(cellList,'Phase_h','Superheated','w'))
        self.w_2phase_h=sum(collect(cellList,'Phase_h','TwoPhase','w'))
        self.w_subcooled_h=sum(collect(cellList,'Phase_h','Subcooled','w'))
        self.w_superheated_c=sum(collect(cellList,'Phase_c','Superheated','w'))
        self.w_2phase_c=sum(collect(cellList,'Phase_c','TwoPhase','w'))
        self.w_subcooled_c=sum(collect(cellList,'Phase_c','Subcooled','w'))
        
        self.DP_superheated_c=sum(collect(cellList,'Phase_c','Superheated','DP_c'))
        self.DP_2phase_c=sum(collect(cellList,'Phase_c','TwoPhase','DP_c'))
        self.DP_subcooled_c=sum(collect(cellList,'Phase_c','Subcooled','DP_c'))
        self.DP_c=self.DP_superheated_c+self.DP_2phase_c+self.DP_subcooled_c
        
        self.DP_superheated_h=sum(collect(cellList,'Phase_h','Superheated','DP_h'))
        self.DP_2phase_h=sum(collect(cellList,'Phase_h','TwoPhase','DP_h'))
        self.DP_subcooled_h=sum(collect(cellList,'Phase_h','Subcooled','DP_h'))
        self.DP_h=self.DP_superheated_h+self.DP_2phase_h+self.DP_subcooled_h
        
        self.Charge_superheated_c=sum(collect(cellList,'Phase_c','Superheated','Charge_c'))
        self.Charge_2phase_c=sum(collect(cellList,'Phase_c','TwoPhase','Charge_c'))
        self.Charge_subcooled_c=sum(collect(cellList,'Phase_c','Subcooled','Charge_c'))
        self.Charge_c=self.Charge_superheated_c+self.Charge_2phase_c+self.Charge_subcooled_c
        self.Charge_superheated_h=sum(collect(cellList,'Phase_h','Superheated','Charge_h'))
        self.Charge_2phase_h=sum(collect(cellList,'Phase_h','TwoPhase','Charge_h'))
        self.Charge_subcooled_h=sum(collect(cellList,'Phase_h','Subcooled','Charge_h'))
        self.Charge_h=self.Charge_superheated_h+self.Charge_2phase_h+self.Charge_subcooled_h
        
        self.Q_superheated_h=sum(collect(cellList,'Phase_h','Superheated','Q'))
        self.Q_2phase_h=sum(collect(cellList,'Phase_h','TwoPhase','Q'))
        self.Q_subcooled_h=sum(collect(cellList,'Phase_h','Subcooled','Q'))
        self.Q_superheated_c=sum(collect(cellList,'Phase_c','Superheated','Q'))
        self.Q_2phase_c=sum(collect(cellList,'Phase_c','TwoPhase','Q'))
        self.Q_subcooled_c=sum(collect(cellList,'Phase_c','Subcooled','Q'))
        
        w_superheat=collect(cellList,'Phase_c','Superheated','w')
        w_2phase=collect(cellList,'Phase_c','TwoPhase','w')
        h_c_sh=collect(cellList,'Phase_c','Superheated','h_c')
        h_c_2phase=collect(cellList,'Phase_c','TwoPhase','h_c')
        self.xout_h=collect(cellList,'Phase_h','TwoPhase','xout_h')
        
        if len(w_superheat)>0:
            self.h_superheated_c=float(sum(np.array(h_c_sh)*np.array(w_superheat))/sum(w_superheat))
        else:
            self.h_superheated_c=0
            
        if len(w_2phase)>0:
            self.h_2phase_c=float(sum(np.array(h_c_2phase)*np.array(w_2phase))/sum(w_2phase))
        else:
            self.h_2phase_c=0
            
        ### Collect all the cells on the hot side
        w_subcooled_h=collect(cellList,'Phase_h','Subcooled','w')
        w_superheat_h=collect(cellList,'Phase_h','Superheated','w')
        w_2phase_h=collect(cellList,'Phase_h','TwoPhase','w')
        h_h_sh=collect(cellList,'Phase_h','Superheated','h_h')
        h_h_2phase=collect(cellList,'Phase_h','TwoPhase','h_h')
        h_h_subcool=collect(cellList,'Phase_h','Subcooled','h_h')
        
        w_subcooled_c=collect(cellList,'Phase_c','Subcooled','w')
        w_superheat_c=collect(cellList,'Phase_c','Superheated','w')
        w_2phase_c=collect(cellList,'Phase_c','TwoPhase','w')
        h_c_sh=collect(cellList,'Phase_c','Superheated','h_c')
        h_c_2phase=collect(cellList,'Phase_c','TwoPhase','h_c')
        h_c_subcool=collect(cellList,'Phase_c','Subcooled','h_c')
        
        if len(w_subcooled_h)>0:
            self.h_subcooled_h=float(sum(np.array(h_h_subcool)*np.array(w_subcooled_h))/sum(w_subcooled_h))
        else:
            self.h_subcooled_h=0
        
        if len(w_2phase_h)>0:
            self.h_2phase_h=float(sum(np.array(h_h_2phase)*np.array(w_2phase_h))/sum(w_2phase_h))
        else:
            self.h_2phase_h=0
            
        if len(w_superheat_h)>0:
            self.h_superheated_h=float(sum(np.array(h_h_sh)*np.array(w_superheat_h))/sum(w_superheat_h))
        else:
            self.h_superheated_h=0
            
        if len(w_subcooled_c)>0:
            self.h_subcooled_c=float(sum(np.array(h_c_subcool)*np.array(w_subcooled_c))/sum(w_subcooled_c))
        else:
            self.h_subcooled_c=0
        
        if len(w_2phase_c)>0:
            self.h_2phase_c=float(sum(np.array(h_c_2phase)*np.array(w_2phase_c))/sum(w_2phase_c))
        else:
            self.h_2phase_c=0
            
        if len(w_superheat_c)>0:
            self.h_superheated_c=float(sum(np.array(h_c_sh)*np.array(w_superheat_c))/sum(w_superheat_c))
        else:
            self.h_superheated_c=0
            
        
        
        self.q_flux=collect(cellList,'Phase_c','TwoPhase','q_flux')
        #TODO: Modified for Incompressible as hot fluid
        if self.Ref_h =='INCOMP::T66':

            self.Tout_h,self.rhoout_h=TrhoPhase_ph(self.Ref_h,self.pin_h,self.hout_h,None,None,None,None)[0:2]
        else:
            self.Tout_h,self.rhoout_h=TrhoPhase_ph(self.Ref_h,self.pin_h,self.hout_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0:2]
        
        if self.Ref_c =='INCOMP::MEG-30%':
            self.Tout_c,self.rhoout_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hout_c,None,None,None,None)[0:2]        
        elif self.Ref_c =='POE':
            self.Tout_c,self.rhoout_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hout_c,None,None,None,None)[0:2]
        elif self.Ref_c == 'ACD100FY':
            self.Tout_c,self.rhoout_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hout_c,None,None,None,None)[0:2]
        else:    
            self.Tout_c,self.rhoout_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hout_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0:2]
        
        if IsFluidType(self.Ref_c,'Brine'):
            self.sout_c=PropsSI('S','T',self.Tout_c,'P',self.pin_c*1000,self.Ref_c)#*1000
            self.DT_sc_c=1e9
        #TODO: Modified for Incompressible as cold fluid
        elif self.Ref_c =='INCOMP::MEG-30%':
            self.sout_c=PropsSI('S','T',self.Tout_c,'P',self.pin_c*1000,self.Ref_c)#*1000
            self.DT_sc_c=1e9            
        elif self.Ref_c =='POE':
            self.sout_c=s_l(self.Ref_c,self.Tout_c)#*1000
            self.DT_sc_c=1e9
        elif self.Ref_c == 'ACD100FY':
            self.sout_c=s_l(self.Ref_c,self.Tout_c)#*1000
            self.DT_sc_c=1e9            
        else:
            self.sout_c=PropsSI('S','T',self.Tout_c,'D',self.rhoout_c,self.Ref_c)#*1000
            #Effective subcooling for both streams
            hsatL=PropsSI('H','T',self.Tbubble_c,'Q',0,self.Ref_c)#*1000
            cpsatL=PropsSI('C','T',self.Tbubble_c,'Q',0,self.Ref_c)#*1000
            if self.hout_c>hsatL:
                #Outlet is at some quality on cold side
                self.DT_sc_c=-(self.hout_c-hsatL)/cpsatL
            else:
                self.DT_sc_c=self.Tbubble_c-self.Tout_c
        
        if IsFluidType(self.Ref_h,'Brine'):
            self.sout_h=PropsSI('S','T',self.Tout_h,'P',self.pin_h*1000,self.Ref_h)#*1000
            self.DT_sc_h=1e9
        #TODO: Modified for Incompressible as hot fluid
        elif self.Ref_h =='INCOMP::T66':
            self.sout_h=PropsSI('S','T',self.Tout_h,'P',self.pin_h*1000,self.Ref_h)#*1000
            self.DT_sc_h=1e9
        else:
            self.sout_h=PropsSI('S','T',self.Tout_h,'D',self.rhoout_h,self.Ref_h)#*1000
            hsatV=PropsSI('H','T',self.Tdew_h,'Q',0,self.Ref_h)#*1000
            cpsatV=PropsSI('C','T',self.Tdew_h,'Q',0,self.Ref_h)#*1000
            if self.hout_h<hsatV:
                #Outlet is at some quality on hot side
                self.DT_sc_h=-(hsatV-self.hout_h)/cpsatV
            else:
                self.DT_sc_h=self.Tout_h - self.Tbubble_h
        
    def eNTU_CounterFlow(self,Cr,Ntu):
        return ((1 - exp(-Ntu * (1 - Cr))) / 
            (1 - Cr * exp(-Ntu * (1 - Cr))))
    
    def _OnePhaseH_OnePhaseC_Qimposed(self,Inputs):
        """
        Single phase on both sides
        Inputs is a dict of parameters
        """
        
        #Calculate the mean temperature
        Tmean_h=Inputs['Tmean_h']
        Tmean_c=Inputs['Tmean_c']
        #Evaluate heat transfer coefficient for both fluids
        h_h,cp_h,PlateOutput_h=self.PlateHTDP(self.Ref_h, Tmean_h, Inputs['pin_h'],self.mdot_h/self.NgapsHot)
        h_c,cp_c,PlateOutput_c=self.PlateHTDP(self.Ref_c, Tmean_c, Inputs['pin_c'],self.mdot_c/self.NgapsCold)
        
        #Use cp calculated from delta h/delta T
        cp_h=Inputs['cp_h']
        cp_c=Inputs['cp_c']
        #Evaluate UA [W/K] if entire HX was in this section 
        UA_total=1/(1/(h_h*self.A_h_wetted)+1/(h_c*self.A_c_wetted)+self.PlateThickness/(self.PlateConductivity*(self.A_c_wetted+self.A_h_wetted)/2.))
        #Get Ntu [-]
        C=[cp_c*self.mdot_c,cp_h*self.mdot_h]
        Cmin=min(C)
        Cr=Cmin/max(C)
        
        #Effectiveness [-]
        Q=Inputs['Q']
        Qmax=Cmin*(Inputs['Tin_h']-Inputs['Tin_c'])
        epsilon = Q/Qmax
        
        if 1 <= epsilon < 1+1e-6:  #if epsilon is slightly larger than 1
            epsilon = 1-1e-12
#         if 1 <= epsilon < 1+1e-3:  #if epsilon is slightly larger than 1
#             print epsilon
#             epsilon = 1 - 1e-6#-1e-12        
        #Pure counterflow with Cr<1 (Incropera Table 11.4)
        #print 'Cr:',Cr, 'epsilon:',epsilon
        NTU=1/(Cr-1)*log((epsilon-1)/(epsilon*Cr-1))
        
        #Required UA value
        UA_req=Cmin*NTU
        
        #w is required part of heat exchanger for this duty
        w=UA_req/UA_total
        
        #Determine both charge components
        rho_h=PropsSI('D','T',Tmean_h, 'P', self.pin_h*1000, self.Ref_h)
        Charge_h = w * self.V_h * rho_h
        
        #TODO: modified to include lubricants
        if self.Ref_c == 'POE':
            rho_c = rho_l(self.Ref_c,Tmean_c)
        elif self.Ref_c == 'ACD100FY':
            rho_c = rho_l(self.Ref_c,Tmean_c)
        else:
            rho_c=PropsSI('D','T',Tmean_c, 'P', self.pin_c*1000+100, self.Ref_c)
        
        Charge_c = w * self.V_c * rho_c
        
        #Pack outputs
        Outputs={
            'w': w,
            'Tout_h': Inputs['Tin_h']-Q/(self.mdot_h*cp_h),
            'Tout_c': Inputs['Tin_c']+Q/(self.mdot_c*cp_c),
            'Charge_c': Charge_c,
            'Charge_h': Charge_h,
            'DP_h': -PlateOutput_h['DELTAP'],
            'DP_c': -PlateOutput_c['DELTAP'],
            'h_h':h_h,
            'h_c':h_c,
            
        }
        return dict(Inputs.items()+Outputs.items())
    
    def _OnePhaseH_OnePhaseC_wimposed(self,Inputs):
        """
        Single phase on both sides
        Inputs is a dict of parameters that is the return value of _OnePhaseH_OnePhaseC_Qimposed
        """
        w = Inputs['w']
        #Calculate the mean temperature
        Tmean_h=Inputs['Tmean_h']  #We don't need to iterate to find the mean temperature because we already know there is more than enough length for this section's effectiveness to be 1.
        Tmean_c=Inputs['Tmean_c']
        #Evaluate heat transfer coefficient for both fluids
        h_h,cp_h,PlateOutput_h=self.PlateHTDP(self.Ref_h, Tmean_h, Inputs['pin_h'],self.mdot_h/self.NgapsHot)
        h_c,cp_c,PlateOutput_c=self.PlateHTDP(self.Ref_c, Tmean_c, Inputs['pin_c'],self.mdot_c/self.NgapsCold)
        
        #Use cp calculated from delta h/delta T
        cp_h=Inputs['cp_h']
        cp_c=Inputs['cp_c']
        #Evaluate UA [W/K] if entire HX was in this section 
        UA_total=1/(1/(h_h*self.A_h_wetted)+1/(h_c*self.A_c_wetted)+self.PlateThickness/(self.PlateConductivity*(self.A_c_wetted+self.A_h_wetted)/2.))
        UA_actual = UA_total*w  #We impose length fraction here so we already know the actual UA value
        #Get Ntu [-] and Effectiveness [-]
        C=[cp_c*self.mdot_c,cp_h*self.mdot_h]
        Cmin=min(C)
        Cr=Cmin/max(C)
        
        NTU = UA_actual/Cmin
        
        epsilon = (1 - exp(-NTU*(1 - Cr)))/(1 - Cr*exp(-NTU*(1 - Cr)))  #Pure counterflow with Cr<1 (Incropera Table 11.3)
        assert(epsilon <= 1)
        
        Qmax=Cmin*(Inputs['Tin_h']-Inputs['Tin_c'])
        Q = Qmax*epsilon
               
        #Determine both charge components
        rho_h=PropsSI('D','T',Tmean_h, 'P', self.pin_h*1000, self.Ref_h)
        Charge_h = w * self.V_h * rho_h
        
        if self.Ref_c == 'POE':
            rho_c = rho_l(self.Ref_c,Tmean_c)
        elif self.Ref_c == 'ACD100FY':
            rho_c = rho_l(self.Ref_c,Tmean_c)
        else:
            rho_c=PropsSI('D','T',Tmean_c, 'P', self.pin_c*1000, self.Ref_c)
        
        Charge_c = w * self.V_c * rho_c
        
        #Pack outputs
        Outputs={
            'Tout_h': Inputs['Tin_h']-Q/(self.mdot_h*cp_h),
            'Tout_c': Inputs['Tin_c']+Q/(self.mdot_c*cp_c),
            'Charge_c': Charge_c,
            'Charge_h': Charge_h,
            'DP_h': -PlateOutput_h['DELTAP'],
            'DP_c': -PlateOutput_c['DELTAP'],
            'h_h':h_h,
            'h_c':h_c,
            'Q_wimposed':Q
        }
        return dict(Inputs.items()+Outputs.items())  #Overwrites any outputs that were passed back in as inputs with the new outputs

    def _OnePhaseH_TwoPhaseC_Qimposed(self,Inputs):
        """
        The hot stream is all single phase, and the cold stream is evaporating
        """
        #Calculate the mean temperature for the single-phase fluid
        h_h,cp_h,PlateOutput_h=self.PlateHTDP(self.Ref_h, Inputs['Tmean_h'], Inputs['pin_h'],self.mdot_h/self.NgapsHot)
        #Use cp calculated from delta h/delta T
        cp_h=Inputs['cp_h']
        #Mole mass of refrigerant for Cooper correlation
        M=PropsSI('M','T',0,'P',0,self.Ref_c)
        #Reduced pressure for Cooper Correlation
        #pstar=Inputs['pin_c']/(Props('E','T',0,'P',0,self.Ref_c))
        pstar=Inputs['pin_c']/(PropsSI(self.Ref_c,'pcrit')/1000)
        change=999
        w=1
        Q=Inputs['Q']
        """
        The Cooper Pool boiling relationship is a function of the heat flux, 
        therefore the heat flux must be iteratively determined
        """
        while abs(change)>1e-6:
            q_flux=Q/(w*self.A_c_wetted)
            
            #Heat transfer coefficient from Cooper Pool Boiling
            h_c_2phase=Cooper_PoolBoiling(pstar,1.0,q_flux,M) #1.5 correction factor comes from Claesson Thesis on plate HX
            
            G=self.mdot_c/self.A_c_flow
            Dh=self.Dh_c
            x=(Inputs['xin_c']+Inputs['xout_c'])/2

            UA_total=1/(1/(h_h*self.A_h_wetted)+1/(h_c_2phase*self.A_c_wetted)+self.PlateThickness/(self.PlateConductivity*self.A_c_wetted))
            C_h=cp_h*self.mdot_h
            
            Qmax=C_h*(Inputs['Tin_h']-Inputs['Tsat_c'])
            epsilon=Q/Qmax
            
            if 1 <= epsilon < 1+1e-6:  #if epsilon is slightly larger than 1
                epsilon = 1-1e-12
            
            NTU=-log(1-epsilon)
            UA_req=NTU*C_h
            
            change=UA_req/UA_total-w
            w=UA_req/UA_total
        
        #Refrigerant charge
        rho_h=PropsSI('D','T',Inputs['Tmean_h'], 'P', self.pin_h*1000, self.Ref_h)
        Charge_h = w * self.V_h * rho_h
        rho_c=TwoPhaseDensity(self.Ref_c,Inputs['xin_c'],Inputs['xout_c'],self.Tdew_c,self.Tbubble_c,slipModel='Zivi')
        Charge_c = rho_c * w * self.V_c
        
        #Use Lockhart Martinelli to calculate the pressure drop.  Claesson found good agreement using C parameter of 4.67
        DP_frict_c=LMPressureGradientAvg(Inputs['xin_c'],Inputs['xout_c'],self.Ref_c,self.mdot_c/self.A_c_flow,self.Dh_c,self.Tbubble_c,self.Tdew_c,C=4.67)*w*self.Lp
        #Accelerational pressure drop component    
        DP_accel_c=AccelPressureDrop(Inputs['xin_c'],Inputs['xout_c'],self.Ref_c,self.mdot_c/self.A_c_flow,self.Tbubble_c,self.Tdew_c)
        
        #Pack outputs
        Outputs={
            'w':w,
            'Tout_h': Inputs['Tin_h']-Q/(self.mdot_h*cp_h),
            'Tout_c': Inputs['Tsat_c'],
            'Charge_c': Charge_c,
            'Charge_h': Charge_h,
            'DP_h': -PlateOutput_h['DELTAP'],
            'DP_c': DP_frict_c+DP_accel_c,
            'h_h':h_h,
            'h_c':h_c_2phase,
            'q_flux':q_flux
        }
        return dict(Inputs.items()+Outputs.items())
    
    def _OnePhaseH_TwoPhaseC_wimposed(self,Inputs):
        """
        The hot stream is all single phase, and the cold stream is evaporating
        """
        w = Inputs['w']
        #Calculate the mean temperature for the single-phase fluid
        h_h,cp_h,PlateOutput_h=self.PlateHTDP(self.Ref_h, Inputs['Tmean_h'], Inputs['pin_h'],self.mdot_h/self.NgapsHot)
        #Use cp calculated from delta h/delta T
        cp_h=Inputs['cp_h']
        #Mole mass of refrigerant for Cooper correlation
        M=PropsSI('M','T',0,'P',0,self.Ref_c)
        #Reduced pressure for Cooper Correlation
        #pstar=Inputs['pin_c']/Props('E','T',0,'P',0,self.Ref_c)
        pstar=Inputs['pin_c']/(PropsSI(self.Ref_c,'pcrit')/1000)

        C_h=cp_h*self.mdot_h
        Q=C_h*(Inputs['Tin_h']-Inputs['Tsat_c'])  #initial guess for Cooper Pool boiling
        change = 999
        """
        The Cooper Pool boiling relationship is a function of the heat flux, 
        therefore the heat flux must be iteratively determined.
        """
        
        while abs(change) > 1e-6:
            q_flux=Q/(w*self.A_c_wetted)
            
            #Heat transfer coefficient from Cooper Pool Boiling
            
            h_c_2phase=Cooper_PoolBoiling(pstar,1.0,q_flux,M) #1.5 correction factor comes from Claesson Thesis on plate HX
            
            G=self.mdot_c/self.A_c_flow
            Dh=self.Dh_c
            x=(Inputs['xin_c']+Inputs['xout_c'])/2
    
            UA_total=1/(1/(h_h*self.A_h_wetted)+1/(h_c_2phase*self.A_c_wetted)+self.PlateThickness/(self.PlateConductivity*self.A_c_wetted))
            UA_actual = UA_total*w
            
            C_h=cp_h*self.mdot_h
            NTU = UA_actual/C_h
            epsilon = 1 - exp(-NTU)
            
            Qmax=C_h*(Inputs['Tin_h']-Inputs['Tsat_c'])
            
            change=Qmax*epsilon-Q
            Q = Qmax*epsilon
        
        #Refrigerant charge
        rho_h=PropsSI('D','T',Inputs['Tmean_h'], 'P', self.pin_h*1000, self.Ref_h)
        Charge_h = w * self.V_h * rho_h
        rho_c=TwoPhaseDensity(self.Ref_c,Inputs['xin_c'],Inputs['xout_c'],self.Tdew_c,self.Tbubble_c,slipModel='Zivi')
        Charge_c = rho_c * w * self.V_c
        
        #Use Lockhart Martinelli to calculate the pressure drop.  Claesson found good agreement using C parameter of 4.67
        DP_frict_c=LMPressureGradientAvg(Inputs['xin_c'],Inputs['xout_c'],self.Ref_c,self.mdot_c/self.A_c_flow,self.Dh_c,self.Tbubble_c,self.Tdew_c,C=4.67)*w*self.Lp
        #Accelerational pressure drop component    
        DP_accel_c=AccelPressureDrop(Inputs['xin_c'],Inputs['xout_c'],self.Ref_c,self.mdot_c/self.A_c_flow,self.Tbubble_c,self.Tdew_c)
        
        #Pack outputs
        Outputs={
            'Tout_h': Inputs['Tin_h']-Q/(self.mdot_h*cp_h),
            'Tout_c': Inputs['Tsat_c'],
            'Charge_c': Charge_c,
            'Charge_h': Charge_h,
            'DP_h': -PlateOutput_h['DELTAP'],
            'DP_c': DP_frict_c+DP_accel_c,
            'h_h':h_h,
            'h_c':h_c_2phase,
            'q_flux':q_flux,
            'Q_wimposed':Q
        }
        return dict(Inputs.items()+Outputs.items())
    
    
    def _TwoPhaseH_OnePhaseC_Qimposed(self,Inputs):
        """
        Hot stream is condensing, cold stream is single phase  
        """
        #Choose which correlation to use.  Shah seems to work better with steam but underpredicts heat transfer for refrigerant
        #It's possible that Shah only works better for steam because it underpredicts the heat transfer coefficient and the HX has more than enough area for maximum heat transfer
        if 'ater' in self.Ref_h:  #'ater' allows to handle both 'Water' and 'water'
            G_h=self.mdot_h/self.A_h_flow
            h_h_2phase=ShahCondensation_Average(Inputs['xout_h'],Inputs['xin_h'],self.Ref_h,G_h,self.Dh_h,Inputs['pin_h'],self.Tbubble_h,self.Tdew_h)
        else:
            h_h_2phase=LongoCondensation((Inputs['xout_h']+Inputs['xin_h'])/2,self.mdot_h/self.A_h_flow,self.Dh_h,self.Ref_h,self.Tbubble_h,self.Tdew_h)
        
        h_c,cp_c,PlateOutput_c=self.PlateHTDP(self.Ref_c, Inputs['Tmean_c'], Inputs['pin_c'],self.mdot_c/self.NgapsCold)
        #Use cp calculated from delta h/delta T
        cp_c=Inputs['cp_c']
        UA_total=1/(1/(h_c*self.A_c_wetted)+1/(h_h_2phase*self.A_h_wetted)+self.PlateThickness/(self.PlateConductivity*(self.A_c_wetted+self.A_h_wetted)/2.))
        C_c=cp_c*self.mdot_c
        
        Q=Inputs['Q']
        Qmax=C_c*(Inputs['Tsat_h']-Inputs['Tin_c'])
        epsilon = Q/(Qmax+1e-6)
        
        if 1 <= epsilon < 1+1e-6:  #if epsilon is slightly larger than 1
            epsilon = 1-1e-12
        
        #Cr = 0, so NTU is simply
        try:
            NTU=-log(1-epsilon)
        except:
            #pass
            raise Exception('epsilon_two_phase =', epsilon)
        UA_req=NTU*C_c
        w=UA_req/UA_total
        
        if self.Ref_c == 'POE':
            rho_c = rho_l(self.Ref_c,Inputs['Tmean_c'])
        elif self.Ref_c == 'ACD100FY':
            rho_c = rho_l(self.Ref_c,Inputs['Tmean_c'])
        else:
            rho_c=PropsSI('D','T',Inputs['Tmean_c'], 'P', self.pin_c*1000, self.Ref_c)
        
        Charge_c = w * self.V_c * rho_c
        rho_h=TwoPhaseDensity(self.Ref_h,Inputs['xout_h'],Inputs['xin_h'],self.Tdew_h,self.Tbubble_h,slipModel='Zivi')
        Charge_h = w * self.V_h * rho_h
        
        #Use Lockhart Martinelli to calculate the pressure drop.  Claesson found good agreement using C parameter of 4.67
        DP_frict_h=LMPressureGradientAvg(Inputs['xin_h'],Inputs['xout_h'],self.Ref_h,self.mdot_h/self.A_h_flow,self.Dh_h,self.Tbubble_h,self.Tdew_h,C=4.67)*w*self.Lp
        #Accelerational pressure drop component    
        DP_accel_h=-AccelPressureDrop(Inputs['xin_h'],Inputs['xout_h'],self.Ref_h,self.mdot_h/self.A_h_flow,self.Tbubble_h,self.Tdew_h)
        
        #Pack outputs
        Outputs={
            'w': w,
            'Tout_c': Inputs['Tin_c']+Q/(self.mdot_c*cp_c+1.0e-6),
            'Tout_h': Inputs['Tsat_h'],
            'DP_c': -PlateOutput_c['DELTAP'],
            'DP_h': DP_frict_h+DP_frict_h,
            'Charge_c':Charge_c,
            'Charge_h':Charge_h,
            'h_h':h_h_2phase,
            'h_c':h_c,
        }
        
        return dict(Inputs.items()+Outputs.items())
    
    def _TwoPhaseH_OnePhaseC_wimposed(self,Inputs):
        """
        Hot stream is condensing, cold stream is single phase  
        """
        w = Inputs['w']
        #Choose which correlation to use.  Shah seems to work better with steam but underpredicts heat transfer for refrigerant
        #It's possible that Shah only works better for steam because it underpredicts the heat transfer coefficient and the HX has more than enough area for maximum heat transfer
        if 'ater' in self.Ref_h:  #'ater' allows to handle both 'Water' and 'water'
            G_h=self.mdot_h/self.A_h_flow
            h_h_2phase=ShahCondensation_Average(Inputs['xout_h'],Inputs['xin_h'],self.Ref_h,G_h,self.Dh_h,Inputs['pin_h'],self.Tbubble_h,self.Tdew_h)
        else:
            h_h_2phase=LongoCondensation((Inputs['xout_h']+Inputs['xin_h'])/2,self.mdot_h/self.A_h_flow,self.Dh_h,self.Ref_h,self.Tbubble_h,self.Tdew_h)
        
        h_c,cp_c,PlateOutput_c=self.PlateHTDP(self.Ref_c, Inputs['Tmean_c'], Inputs['pin_c'],self.mdot_c/self.NgapsCold)
        #Use cp calculated from delta h/delta T
        cp_c=Inputs['cp_c']
        UA_total=1/(1/(h_c*self.A_c_wetted)+1/(h_h_2phase*self.A_h_wetted)+self.PlateThickness/(self.PlateConductivity*(self.A_c_wetted+self.A_h_wetted)/2.))
        UA_actual = UA_total*w
        
        C_c=cp_c*self.mdot_c
        NTU = UA_actual/C_c
        epsilon = 1 - exp(-NTU)
        
        Qmax=C_c*(Inputs['Tsat_h']-Inputs['Tin_c'])
        Q = Qmax*epsilon
        
        if self.Ref_c == 'POE':
            rho_c = rhp_l(self.Ref_c,Inputs['Tmean_c'])
        elif self.Ref_c == 'ACD100FY':
            rho_c = rhp_l(self.Ref_c,Inputs['Tmean_c'])
        else:
            rho_c=PropsSI('D','T',Inputs['Tmean_c'], 'P', self.pin_c*1000, self.Ref_c)
        
        Charge_c = w * self.V_c * rho_c
        rho_h=TwoPhaseDensity(self.Ref_h,Inputs['xout_h'],Inputs['xin_h'],self.Tdew_h,self.Tbubble_h,slipModel='Zivi')
        Charge_h = w * self.V_h * rho_h
        
        #Use Lockhart Martinelli to calculate the pressure drop.  Claesson found good agreement using C parameter of 4.67
        DP_frict_h=LMPressureGradientAvg(Inputs['xin_h'],Inputs['xout_h'],self.Ref_h,self.mdot_h/self.A_h_flow,self.Dh_h,self.Tbubble_h,self.Tdew_h,C=4.67)*w*self.Lp
        #Accelerational pressure drop component    
        DP_accel_h=-AccelPressureDrop(Inputs['xin_h'],Inputs['xout_h'],self.Ref_h,self.mdot_h/self.A_h_flow,self.Tbubble_h,self.Tdew_h)
        
        #Pack outputs
        Outputs={
            'Tout_c': Inputs['Tin_c']+Q/(self.mdot_c*cp_c),
            'Tout_h': Inputs['Tsat_h'],
            'DP_c': -PlateOutput_c['DELTAP'],
            'DP_h': DP_frict_h+DP_frict_h,
            'Charge_c':Charge_c,
            'Charge_h':Charge_h,
            'h_h':h_h_2phase,
            'h_c':h_c,
            'Q_wimposed':Q
        }
        
        return dict(Inputs.items()+Outputs.items())
    
    def _TwoPhaseH_TwoPhaseC_Qimposed(self,Inputs):
        """
        Hot stream is condensing, cold stream is evaporating 
        """
        #Hot side: Shah correlation for steam condensation
        G_h=self.mdot_h/self.A_h_flow
        h_h_2phase=ShahCondensation_Average(Inputs['xout_h'],Inputs['xin_h'],self.Ref_h,G_h,self.Dh_h,Inputs['pin_h'],self.Tbubble_h,self.Tdew_h)
        #Cold side: Cooper Pool Boiling Correlation
        #Mole mass of refrigerant for Cooper correlation
        M=PropsSI('M','T',0,'P',0,self.Ref_c)
        #Reduced pressure for Cooper Correlation
        #pstar=Inputs['pin_c']/Props('E','T',0,'P',0,self.Ref_c)
        pstar=Inputs['pin_c']/(PropsSI(self.Ref_c,'pcrit')/1000)  
        """
    The Cooper Pool boiling relationship is a function of the heat flux, which is known in this particular case, but w
    has to be determined by iteration
    """
        change=999
        w=1
        while abs(change)>1e-6:
            q_flux=Inputs['Q']/(1.0e-6+w*self.A_c_wetted)
            #Heat transfer coefficient from Cooper Pool Boiling
            h_c_2phase=Cooper_PoolBoiling(pstar,1.0,q_flux,M) #1.5 correction factor comes from Claesson Thesis on plate HX
   
            G=self.mdot_c/self.A_c_flow
            Dh=self.Dh_c
            x=(Inputs['xin_c']+Inputs['xout_c'])/2
            DELTAT=Inputs['Tsat_h']-Inputs['Tsat_c']
            UA_req=Inputs['Q']/DELTAT
            UA_total=1/(1/(1.0e-6+h_h_2phase*self.A_h_wetted)+1/(1.0e-6+h_c_2phase*self.A_c_wetted)+self.PlateThickness/(self.PlateConductivity*self.A_c_wetted))
            change=UA_req/UA_total-w
            w=UA_req/UA_total
        #print UA_req,DELTAT
        #Refrigerant charge
        rho_h=TwoPhaseDensity(self.Ref_h,Inputs['xin_h'],Inputs['xout_h'],self.Tdew_h,self.Tbubble_h,slipModel='Zivi')
        Charge_h = w * self.V_h * rho_h
        rho_c=TwoPhaseDensity(self.Ref_c,Inputs['xin_c'],Inputs['xout_c'],self.Tdew_c,self.Tbubble_c,slipModel='Zivi')
        Charge_c = rho_c * w * self.V_c
     
        #Use Lockhart Martinelli to calculate the pressure drop.  Claesson found good agreement using C parameter of 4.67
        DP_frict_c=LMPressureGradientAvg(Inputs['xin_c'],Inputs['xout_c'],self.Ref_c,self.mdot_c/self.A_c_flow,self.Dh_c,self.Tbubble_c,self.Tdew_c,C=4.67)*w*self.Lp
        #Accelerational pressure drop component    
        DP_accel_c=AccelPressureDrop(Inputs['xin_c'],Inputs['xout_c'],self.Ref_c,self.mdot_c/self.A_c_flow,self.Tbubble_c,self.Tdew_c)
        #Use Lockhart Martinelli to calculate the pressure drop.  Claesson found good agreement using C parameter of 4.67
        DP_frict_h=LMPressureGradientAvg(Inputs['xin_h'],Inputs['xout_h'],self.Ref_h,self.mdot_h/self.A_h_flow,self.Dh_h,self.Tbubble_h,self.Tdew_h,C=4.67)*w*self.Lp
        #Accelerational pressure drop component    
        DP_accel_h=AccelPressureDrop(Inputs['xin_h'],Inputs['xout_h'],self.Ref_h,self.mdot_h/self.A_h_flow,self.Tbubble_h,self.Tdew_h)
    
        #Pack outputs
        Outputs={
            'w':w,
            'Charge_c': Charge_c,
            'Charge_h': Charge_h,
            'Tout_c': Inputs['Tsat_c'],
            'Tout_h': Inputs['Tsat_h'],
            'DP_h': DP_frict_h+DP_accel_h,
            'DP_c': DP_frict_c+DP_accel_c,
            'h_h':h_h_2phase,
            'h_c':h_c_2phase,
            'q_flux':q_flux
        }
        return dict(Inputs.items()+Outputs.items())
        
    def Calculate(self):
        """
        Calculate the PHE
        
        """
        
        # Allocate channels between hot and cold streams
        if not hasattr(self,'MoreChannels') or self.MoreChannels not in ['Hot','Cold']:
            raise KeyError("MoreChannels not found, options are 'Hot' or 'Cold'")
        #There are (Nplates - 1) gaps between the plates
        if self.MoreChannels=='Hot':
            #Hot stream gets the extra channel
            self.NgapsHot=(self.Nplates-1)//2+1
            self.NgapsCold=self.Nplates-1-self.NgapsHot
        else:
            #Cold stream gets the extra channel
            self.NgapsCold=(self.Nplates-1)//2+1
            self.NgapsHot=self.Nplates-1-self.NgapsCold
        
        #Saturation temperatures for cold fluid
        #TODO: adjusting for incompressible liquid
        if self.Ref_c =='INCOMP::MEG-30%':        
            self.Tin_c,self.rhoin_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hin_c,None,None,None,None)[0:2]
        elif self.Ref_c =='POE':
            self.Tin_c,self.rhoin_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hin_c,None,None,None,None)[0:2]
        elif self.Ref_c == 'ACD100FY':
            self.Tin_c,self.rhoin_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hin_c,None,None,None,None)[0:2]
        else:
            self.Tbubble_c=Tsat(self.Ref_c,self.pin_c,0,0)
            self.Tdew_c=Tsat(self.Ref_c,self.pin_c,1,0)
            self.Tsat_c=(self.Tbubble_c+self.Tdew_c)/2.0
            if IsFluidType(self.Ref_c,'Brine'):
                self.rhosatL_c=1
                self.rhosatV_c=1
            else:
                self.rhosatL_c=PropsSI('D','T',self.Tbubble_c,'Q',0,self.Ref_c)
                self.rhosatV_c=PropsSI('D','T',self.Tdew_c,'Q',1,self.Ref_c)
            self.Tin_c,self.rhoin_c=TrhoPhase_ph(self.Ref_c,self.pin_c,self.hin_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0:2]  
                
        #Saturation temperatures for hot fluid
        #TODO: adjusting for incompressible liquid
        if self.Ref_h =='INCOMP::T66':
            self.Tin_h,self.rhoin_h=TrhoPhase_ph(self.Ref_h,self.pin_h,self.hin_h,None,None,None,None)[0:2]

        else:    
            self.Tbubble_h=Tsat(self.Ref_h,self.pin_h,0,0)
            self.Tdew_h=Tsat(self.Ref_h,self.pin_h,1,0)
            self.Tsat_h=(self.Tbubble_h+self.Tdew_h)/2.0
            if IsFluidType(self.Ref_h,'Brine'):
                self.rhosatL_h=1
                self.rhosatV_h=1
            elif self.Ref_h =='INCOMP::T66':
                self.rhosatL_h=1
                self.rhosatV_h=1
            
            else:
                self.rhosatL_h=PropsSI('D','T',self.Tbubble_h,'Q',0,self.Ref_h)
                self.rhosatV_h=PropsSI('D','T',self.Tdew_h,'Q',1,self.Ref_h)
            
            #The rest of the inlet states
            self.Tin_h,self.rhoin_h=TrhoPhase_ph(self.Ref_h,self.pin_h,self.hin_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0:2]
            #print 'self.rhoin_h [kg/m3]:',self.rhoin_h

        #TODO: Adjusted for incompressible cold fluid and lubricants
        if IsFluidType(self.Ref_c,'Brine'):
            self.sin_c=PropsSI('S','T',self.Tin_c,'P',self.pin_c*1000,self.Ref_c)
        elif self.Ref_c =='INCOMP::MEG-30%':
            self.sin_c=PropsSI('S','T',self.Tin_c,'P',self.pin_c*1000,self.Ref_c)
        elif self.Ref_c =='POE':
            self.sin_c = s_l(self.Ref_c,self.Tin_c)
        elif self.Ref_c == 'ACD100FY':
            self.sin_c = s_l(self.Ref_c,self.Tin_c)
        else:
            #self.sin_c=PropsSI('S','T',self.Tin_c,'D',self.rhoin_c,self.Ref_c)
            self.sin_c=PropsSI('S','T',self.Tin_c,'P',self.pin_c*1000+100,self.Ref_c)
        if IsFluidType(self.Ref_h,'Brine'):
            self.sin_h=PropsSI('S','T',self.Tin_h,'P',self.pin_h*1000,self.Ref_h)
        elif self.Ref_h =='INCOMP::T66':
            self.sin_h=PropsSI('S','D',self.rhoin_h,'P',self.pin_h*1000,self.Ref_h)
        else:
            self.sin_h=PropsSI('S','T',self.Tin_h,'P',self.pin_h,self.Ref_h)
            
        #TODO saturation entropy cold and hot sides
        if self.Ref_h =='INCOMP::T66':
            self.s_h_liq = 1
            self.s_h_vap = 1
        else:    
            self.s_h_liq= PropsSI('S','P',self.pin_h*1000,'Q',0,self.Ref_h)/1000
            self.s_h_vap= PropsSI('S','P',self.pin_h*1000,'Q',1,self.Ref_h)/1000
       
        if self.Ref_c =='INCOMP::MEG-30%':
            self.s_c_liq = 1
            self.s_c_vap = 1  
        #TODO:  missing case lubricant/incompressible hot source
        elif self.Ref_c =='POE':
            self.s_c_liq = 1
            self.s_c_vap = 1        
        elif self.Ref_c == 'ACD100FY':
            self.s_c_liq = 1
            self.s_c_vap = 1            
        else:
            self.s_c_liq= PropsSI('S','P',self.pin_c*1000,'Q',0,self.Ref_c)/1000
            self.s_c_vap= PropsSI('S','P',self.pin_c*1000,'Q',1,self.Ref_c)/1000

        
        
        
        # Find HT and Delta P on the hot side
        #---------------
        #Mean values for the hot side based on average of inlet temperatures
        HotPlateInputs={
            'PlateAmplitude': self.PlateAmplitude,
            'PlateWavelength' : self.PlateWavelength,
            'InclinationAngle': self.InclinationAngle,
            'Bp': self.Bp,
            'Lp': self.Lp
        }
        HotPlateOutputs=PHE_1phase_hdP(HotPlateInputs,JustGeo=True)
        #There are (Nplates-2) active plates (outer ones don't do anything)
        self.A_h_wetted=HotPlateOutputs['Ap']*(self.Nplates-2)
        self.V_h=HotPlateOutputs['Vchannel']*self.NgapsHot
        self.A_h_flow=HotPlateOutputs['Aflow']*self.NgapsHot
        self.Dh_h=HotPlateOutputs['Dh']
        
        # Find geometric parameters for cold side of plates
        ColdPlateInputs={
            'PlateAmplitude': self.PlateAmplitude,
            'PlateWavelength' : self.PlateWavelength,
            'InclinationAngle': self.InclinationAngle,
            'Bp': self.Bp,
            'Lp': self.Lp
        }
        ColdPlateOutputs=PHE_1phase_hdP(ColdPlateInputs,JustGeo=True)
        #There are (Nplates-2) active plates (outer ones don't do anything)
        self.A_c_wetted=ColdPlateOutputs['Ap']*(self.Nplates-2)
        self.V_c=ColdPlateOutputs['Vchannel']*self.NgapsCold
        self.A_c_flow=ColdPlateOutputs['Aflow']*self.NgapsCold
        self.Dh_c=HotPlateOutputs['Dh']
        
        #Figure out the limiting rate of heat transfer
        self.Qmax=self.DetermineHTBounds()
        
        def GivenQ(Q):
            """
            In this function, the heat transfer rate is imposed.  Therefore the
            outlet states for both fluids are known, and each element can be solved
            analytically in one shot without any iteration.
            """
            
            EnthalpyList_c,EnthalpyList_h=self.BuildEnthalpyLists(Q)
                
           #Plot temperature v. h profiles
#             for i in range(len(EnthalpyList_c)-1):
#                 hc=np.linspace(EnthalpyList_c[i],EnthalpyList_c[i+1])
#                 Tc=np.zeros_like(hc)
#                 for j in range(len(hc)):
#                     Tc[j],r,Ph=TrhoPhase_ph(self.Ref_c,self.pin_c,hc[j],self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)
#                 pylab.plot(self.mdot_c*(hc-EnthalpyList_c[0])/1000,Tc,'b')
#                
#             for i in range(len(EnthalpyList_h)-1):
#                 hh=np.linspace(EnthalpyList_h[i],EnthalpyList_h[i+1])
#                 Th=np.zeros_like(hh)
#                 for j in range(len(hh)):
#                     Th[j],r,Ph=TrhoPhase_ph(self.Ref_h,self.pin_h,hh[j],self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)
#                 pylab.plot(self.mdot_h*(hh-EnthalpyList_h[0])/1000,Th,'r')
#             pylabshow.()
               
            #Ph(self.Ref_h)
#             pylab.plot(np.array(EnthalpyList_h)/1000,self.pin_h*np.ones_like(EnthalpyList_h))
#             pylab.show()
            
            I_h=0
            I_c=0
            wList=[]
            cellList=[]
            while I_h<len(EnthalpyList_h)-1:
                #Heat Transfer Occuring in Each Cell
                #We already built an enthalpy list above that balances energy so we can use either the hot or cold side to find the heat transfer rate in the cell    
                Qbound_h=self.mdot_h*(EnthalpyList_h[I_h+1]-EnthalpyList_h[I_h])
                Qbound=Qbound_h
                
                #Figure out the inlet and outlet enthalpy for this cell
                hout_h=EnthalpyList_h[I_h]
                hin_h=EnthalpyList_h[I_h+1]
                hin_c=EnthalpyList_c[I_c]
                hout_c=EnthalpyList_c[I_c+1]
                
                # Figure out what combination of phases you have:
                # -------------------------------------------------
                # Hot stream is either single phase or condensing
                # Cold stream is either single phase or evaporating
                
                #Use midpoint enthalpies to figure out the phase in the cell
                if self.Ref_h =='INCOMP::T66':
                    Phase_h=Phase_ph(self.Ref_h,self.pin_h,(hin_h+hout_h)/2,None,None,None,None)
                    #print self.Ref_h,Phase_h
                else:
                    Phase_h=Phase_ph(self.Ref_h,self.pin_h,(hin_h+hout_h)/2,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)
                
                
                if self.Ref_c =='INCOMP::MEG-30%':
                    Phase_c=Phase_ph(self.Ref_c,self.pin_c,(hin_c+hout_c)/2,None,None,None,None)
                    #print self.Ref_c,Phase_c
                elif self.Ref_c == 'POE':
                    Phase_c=Phase_ph(self.Ref_c,self.pin_c,(hin_c+hout_c)/2,None,None,None,None)
                elif self.Ref_c == 'ACD100FY':
                    Phase_c=Phase_ph(self.Ref_c,self.pin_c,(hin_c+hout_c)/2,None,None,None,None)
                else:
                    Phase_c=Phase_ph(self.Ref_c,self.pin_c,(hin_c+hout_c)/2,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)
                
                #Determine inlet and outlet temperatures to the cell ([0] gives the first element of the tuple which is temeperature)
                if self.Ref_h =='INCOMP::T66':
                    Tin_h=TrhoPhase_ph(self.Ref_h,self.pin_h,hin_h,None,None,None,None)[0]
                    Tout_h=TrhoPhase_ph(self.Ref_h,self.pin_h,hout_h,None,None,None,None)[0]
                else:
                
                    Tin_h=TrhoPhase_ph(self.Ref_h,self.pin_h,hin_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0]        
                    Tout_h=TrhoPhase_ph(self.Ref_h,self.pin_h,hout_h,self.Tbubble_h,self.Tdew_h,self.rhosatL_h,self.rhosatV_h)[0]
                
                if self.Ref_c =='INCOMP::MEG-30%':
                    Tin_c=TrhoPhase_ph(self.Ref_c,self.pin_c,hin_c,None,None,None,None)[0]
                    Tout_c=TrhoPhase_ph(self.Ref_c,self.pin_c,hout_c,None,None,None,None)[0]
                elif self.Ref_c =='POE':
                    Tin_c=TrhoPhase_ph(self.Ref_c,self.pin_c,hin_c,None,None,None,None)[0]
                    Tout_c=TrhoPhase_ph(self.Ref_c,self.pin_c,hout_c,None,None,None,None)[0]
                elif self.Ref_c == 'ACD100FY':
                    Tin_c=TrhoPhase_ph(self.Ref_c,self.pin_c,hin_c,None,None,None,None)[0]
                    Tout_c=TrhoPhase_ph(self.Ref_c,self.pin_c,hout_c,None,None,None,None)[0]                    
                else:                                
                    Tin_c=TrhoPhase_ph(self.Ref_c,self.pin_c,hin_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0]
                    Tout_c=TrhoPhase_ph(self.Ref_c,self.pin_c,hout_c,self.Tbubble_c,self.Tdew_c,self.rhosatL_c,self.rhosatV_c)[0]
                    #print 'Tin_c:',Tin_c
                    #print 'Tout_c:',Tout_c
                
                if Phase_h in ['Subcooled','Superheated'] and Phase_c in ['Subcooled','Superheated']:
                    # Both are single-phase
                    Inputs={
                        'Q':Qbound,
                        'cp_h':(hin_h-hout_h)/(1.0e-6+Tin_h-Tout_h),
                        'cp_c':(hin_c-hout_c)/(1.0e-6+Tin_c-Tout_c),
                        'Tmean_h':(Tin_h+Tout_h)/2,
                        'Tmean_c':(Tin_c+Tout_c)/2,
                        'Tin_h':Tin_h,
                        'Tin_c':Tin_c,
                        'pin_h':self.pin_h,
                        'pin_c':self.pin_c,
                        'Phase_c':Phase_c,
                        'Phase_h':Phase_h
                    }
                    Outputs=self._OnePhaseH_OnePhaseC_Qimposed(Inputs)
                    Outputs['Tout_c_balance'] = Tout_c
                    Outputs['Tout_h_balance'] = Tout_h
                    wList.append(Outputs['w'])
                    cellList.append(Outputs)
                    if self.Verbosity>6:
                        print 'w[1-1]: ', Outputs['w']
                elif Phase_h=='TwoPhase' and Phase_c in ['Subcooled','Superheated']:
                    # Hot stream is condensing, and cold stream is single-phase (SH or SC)
                    # TODO: bounding state can be saturated state if hot stream is condensing
                    #Must be two-phase so quality is defined
                    xin_h=(hin_h-self.hsatL_h)/(1.0e-6+self.hsatV_h-self.hsatL_h)
                    xout_h=(hout_h-self.hsatL_h)/(1.0e-6+self.hsatV_h-self.hsatL_h)
                    Inputs={
                        'Q':Qbound,
                        'xin_h':xin_h,
                        'xout_h':xout_h,
                        'Tsat_h':self.Tsat_h,
                        'Tmean_c':(Tin_c+Tout_c)/2,
                        'cp_c':(hin_c-hout_c)/(1.0e-6+Tin_c-Tout_c),
                        'Tin_h':Tin_h,
                        'Tin_c':Tin_c,
                        'pin_h':self.pin_h,
                        'pin_c':self.pin_c,
                        'Phase_c':Phase_c,
                        'Phase_h':Phase_h
                    }
                    Outputs=self._TwoPhaseH_OnePhaseC_Qimposed(Inputs)
                    Outputs['Tout_c_balance'] = Tout_c
                    Outputs['Tout_h_balance'] = Tout_h
                    if self.Verbosity>6:
                        print 'w[2-1]: ', Outputs['w']
                    wList.append(Outputs['w'])
                    cellList.append(Outputs)
                elif Phase_c=='TwoPhase' and Phase_h in ['Subcooled','Superheated']:
                    # Cold stream is evaporating, and hot stream is single-phase (SH or SC)
                    
                    #Must be two-phase so quality is defined
                    xin_c=(hin_c-self.hsatL_c)/(self.hsatV_c-self.hsatL_c)
                    xout_c=(hout_c-self.hsatL_c)/(self.hsatV_c-self.hsatL_c)
                    
                    Inputs={
                        'Q':Qbound,
                        'xin_c':xin_c,
                        'xout_c':xout_c,
                        'Tsat_c':self.Tsat_c,
                        'cp_h':(hin_h-hout_h)/(1.0e-6+Tin_h-Tout_h),
                        'Tmean_h':(Tin_h+Tout_h)/2,
                        'Tin_h':Tin_h,
                        'Tin_c':Tin_c,
                        'pin_h':self.pin_h,
                        'pin_c':self.pin_c,
                        'Phase_c':Phase_c,
                        'Phase_h':Phase_h
                    }
                    Outputs=self._OnePhaseH_TwoPhaseC_Qimposed(Inputs)
                    Outputs['Tout_c_balance'] = Tout_c
                    Outputs['Tout_h_balance'] = Tout_h
                    if self.Verbosity>6:
                        print 'w[1-2]: ', Outputs['w']
                    wList.append(Outputs['w'])
                    cellList.append(Outputs)
                    
                elif Phase_c=='TwoPhase' and Phase_h=='TwoPhase':
                    # Cold stream is evaporating, and hot stream is condensing
                    
                    #Must be two-phase so quality is defined
                    xin_c=(hin_c-self.hsatL_c)/(1.0e-6+self.hsatV_c-self.hsatL_c)
                    xout_c=(hout_c-self.hsatL_c)/(1.0e-6+self.hsatV_c-self.hsatL_c)
                    xin_h=(hin_h-self.hsatL_h)/(1.0e-6+self.hsatV_h-self.hsatL_h)
                    xout_h=(hout_h-self.hsatL_h)/(1.0e-6+self.hsatV_h-self.hsatL_h)
                    
                    Inputs={
                        'Q':Qbound,
                        'xin_c':xin_c,
                        'xout_c':xout_c,
                        'xin_h':xin_h,
                        'xout_h':xout_h,
                        'Tsat_c':self.Tsat_c,
                        'Tsat_h':self.Tsat_h,
                        'Tin_h':Tin_h,
                        'Tin_c':Tin_c,
                        'pin_h':self.pin_h,
                        'pin_c':self.pin_c,
                        'Phase_c':Phase_c,
                        'Phase_h':Phase_h
                    }
                    Outputs=self._TwoPhaseH_TwoPhaseC_Qimposed(Inputs)
                    Outputs['Tout_c_balance'] = Tout_c
                    Outputs['Tout_h_balance'] = Tout_h
                    if self.Verbosity>6:
                        print 'w[2-2]: ', Outputs['w']
                    wList.append(Outputs['w'])
                    cellList.append(Outputs)
                    
                I_h+=1
                I_c+=1
            #end while loop
            
            self.cellList=cellList
            if self.Verbosity>6:
                print 'wsum:', np.sum(wList)
            return np.sum(wList)-1.0
        try:
            brentq(GivenQ,0.01*self.Qmax,self.Qmax)#,xtol=0.000001*self.Qmax)
        except ValueError as e:
            if e.args[0]=='f(a) and f(b) must have different signs':
                #if we get this error, we assume Qmax is actually achieved.  It means there was more than enough
                #area to reach Qmax so np.sum(wList)-1.0 above is still negative even when Qmax is the input  
                print "brentq Exception Occurred!!!"
                GivenQ(self.Qmax)
                
                #check which end of the HX has a DELTA_T of 0
                if self.cellList[0]['Tin_c']-1e-6 < self.cellList[0]['Tout_h'] < self.cellList[0]['Tin_c']+1e-6:  #if outlet of hot stream is equal to inlet of cold stream within a tolerance
                    longcell = 0
                elif self.cellList[-1]['Tin_h']-1e-6 < self.cellList[-1]['Tout_c'] < self.cellList[-1]['Tin_h']+1e-6:
                    longcell = len(self.cellList)-1
                else:
                    raise ValueError('We have a problem. Neither end cell has a DELTA_T of 0!')
                #end if
                
                #find the length fraction of the long cell
                w_totalused = 0
                for dictio in self.cellList:
                    w_totalused += dictio['w']
                #end for
                w_leftover = 1 - w_totalused  #this is the length fraction that was not needed to get maximum heat transfer
                
                self.cellList[longcell]['w'] += w_leftover  #add the unneeded length fraction to the required length fraction for the long cell
                                
                #do something depending on the phases we have in that cell

                if self.cellList[longcell]['Phase_h'] in ['Subcooled','Superheated'] \
                and self.cellList[longcell]['Phase_c'] in ['Subcooled','Superheated']:  #then we have single phase flow only
                    #overwrite the cell with the outputs we get from imposing the oversized length fraction on the cell
                    self.cellList[longcell] = self._OnePhaseH_OnePhaseC_wimposed(self.cellList[longcell])
                    
                elif self.cellList[longcell]['Phase_h']=='TwoPhase' \
                and self.cellList[longcell]['Phase_c'] in ['Subcooled','Superheated']:
                    
                    self.cellList[longcell] = self._TwoPhaseH_OnePhaseC_wimposed(self.cellList[longcell])
                    
                elif self.cellList[longcell]['Phase_h'] in ['Subcooled','Superheated'] \
                and self.cellList[longcell]['Phase_c']=='TwoPhase':
                    
                    self.cellList[longcell] = self._OnePhaseH_TwoPhaseC_wimposed(self.cellList[longcell])
                else:
                    raise ValueError('Phases in long cell are not valid')
                
            else:
                raise  #re-raise the exception if it is a different ValueError than the one above
            #end if
        # Collect parameters from all the pieces
        self.PostProcess(self.cellList)
    #end Calculate()
        
def Example_Evaporator():

    Ref = 'R245fa' 
    Hot = 'INCOMP::T66'
    Tin_h_list= [125]  
    pin_h_list= [200] 
    mdot_h_list=[2]
    Tin_c_list= [50]
    pin_c_list=[1000]
    mdot_c_list=[0.5]
        
    i = 0

    for Tin_h, mdot_c, Tin_c, pin_c, mdot_h in zip(Tin_h_list, mdot_c_list, Tin_c_list, pin_c_list, mdot_h_list):
        
        pin_h = 200
        
        params={
                'Ref_c':Ref,
                'mdot_c':mdot_c,
                'pin_c':pin_c,
                'hin_c':PropsSI('H','T',Tin_c+273.15,'P',pin_c*1000,Ref),

                'Ref_h':Hot,
                'mdot_h':mdot_h,
                'pin_h':pin_h,
                'hin_h':PropsSI('H','T',Tin_h+273.15,'P',pin_h*1000,Hot),
        
        
        #Geometric parameters
        'Bp' : 0.119,
        'Lp' : 0.526, #Center-to-center distance between ports
        'Nplates' : 110,
        'PlateAmplitude' : 0.00102, #[m]
        'PlateThickness' : 0.0003, #[m]
        'PlateWavelength' : 0.0066, #[m]
        'InclinationAngle' : pi/3,#[rad]
        'PlateConductivity' : 15.0, #[W/m-K]
        'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
        'Verbosity':6
        }
        PHE=PHEHXClass(**params)
        PHE.Calculate()
        #print PHE.OutputList()

        wlist = [0]
        Tc_list = [PHE.cellList[0]['Tin_c'] - 273.15]
        Tc_balance_list = [PHE.cellList[0]['Tin_c']]
        
        Th_list = [PHE.cellList[0]['Tout_h'] - 273.15]
        Th_balance_list = [PHE.cellList[0]['Tout_h_balance']]
        for dictio in PHE.cellList:
            Tc_list.append(dictio['Tout_c'] - 273.15)
            Tc_balance_list.append(dictio['Tout_c_balance'])
            
            Th_list.append(dictio['Tin_h'] - 273.15)
            Th_balance_list.append(dictio['Tin_h'])
            
            wlist.append(wlist[-1] + dictio['w'])

        figure()
        plot(wlist, Tc_list,'g-o',lw = 2, ms = 10, label=str.split(PHE.Ref_c,'.')[0])
        plot(wlist, Th_list,'b-o',lw = 2, ms = 10, label=str.split(PHE.Ref_h,'.')[0])
        xlim(0,1)
        xlabel('HX Length Fraction [-]')
        ylabel('Fluid Temperature [K]')
        leg=legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.5)
        
        show()
        
        i += 1
    #end for
    #show()

def Example_Condenser():

    Ref = 'R245fa' 
    Cold = 'water'
    Tin_h_list= [87.99]  
    pin_h_list= [150] 
    mdot_h_list=[0.2946]
    Tin_c_list= [17.5]
    pin_c_list=[200]
    mdot_c_list=[4.166]
       
    i = 0
    zip(Tin_c_list,Tin_h_list,pin_h_list,pin_c_list,mdot_c_list,mdot_h_list)
    for Tin_c,Tin_h,pin_h,pin_c,mdot_c,mdot_h in zip(Tin_c_list,Tin_h_list,pin_h_list,pin_c_list,mdot_c_list,mdot_h_list):
        params={
        'Ref_c':Cold,
        'mdot_c':mdot_c,
        'pin_c':pin_c,
        'hin_c':PropsSI('H','T',Tin_c+273.15,'P',pin_c*1000,Cold),

        'Ref_h':Ref,
        'mdot_h':mdot_h,
        'pin_h':pin_h,
        'hin_h':PropsSI('H','T',Tin_h+273.15,'P',pin_h*1000,Ref),
        
        #Geometric parameters
        'Bp' : 0.119,
        'Lp' : 0.526, #Center-to-center distance between ports
        'Nplates' : 110,
        'PlateAmplitude' : 0.00102, #[m]
        'PlateThickness' : 0.0003, #[m]
        'PlateWavelength' : 0.0066, #[m]
        'InclinationAngle' : pi/3,#[rad]
        'PlateConductivity' : 15.0, #[W/m-K]
        'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
        'Verbosity':6
    }
        PHE=PHEHXClass(**params)
        PHE.Calculate()
        #print PHE.OutputList()

        wlist = [0]
        Tc_list = [PHE.cellList[0]['Tin_c'] - 273.15]
        Tc_balance_list = [PHE.cellList[0]['Tin_c']]
        
        Th_list = [PHE.cellList[0]['Tout_h'] - 273.15]
        Th_balance_list = [PHE.cellList[0]['Tout_h_balance']]
        for dictio in PHE.cellList:
            Tc_list.append(dictio['Tout_c'] - 273.15)
            Tc_balance_list.append(dictio['Tout_c_balance'])
            
            Th_list.append(dictio['Tin_h'] - 273.15)
            Th_balance_list.append(dictio['Tin_h'])
            
            wlist.append(wlist[-1] + dictio['w'])

        figure()
        plot(wlist, Tc_list,'b-o',lw = 2,ms = 10, label=str.split(PHE.Ref_c,'.')[0])
        plot(wlist, Th_list,'g-o',lw = 2, ms = 10, label=str.split(PHE.Ref_h,'.')[0])
        xlim(0,1)
        xlabel('HX Length Fraction [-]')
        ylabel('Fluid Temperature [K]')
        leg=legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.5)
        show()
        
        i += 1
    
if __name__=='__main__':

    Example_Evaporator()
    Example_Condenser()