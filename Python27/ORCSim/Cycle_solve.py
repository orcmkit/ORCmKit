"""
ASME ORC 2015
3rd International Seminar on ORC Systems
12-14 October 2015, Brussels, Belgium


Example of detailed ORC with regenerator  model with subcooling and charge 
based solvers. The linesets are included for the charge estimation, but not 
to evaluate pressure drops between main cycle components. 
Minor adjustments can be done inside Cycle.py to account for lineset pressure drops


Reference Paper:

D. Ziviani et al. "ORCSIM: a generalized organic Rankine Cycle Simulation Tool". Paper#30

     
For illustrative purposes only!
"""

from __future__ import division
from Cycle import ORCClass

from math import pi
from ACHPTools import struct, DataHolderFromFile, Write2CSV, dict_to_file
from CoolProp.CoolProp import PropsSI 
import CoolProp
import numpy as np



def ORCSystem_Kortrijk(Inputs, run):    

    """
    The example refers to the ORC installation at UGent Campus Kortrijk. 
    The ORC system includes:
    - three identical SWEP PHEX (evaporator, condenser, regenerator)
    - one Calpeda multi-stage centrifugal pump
    - a single-screw expander
    - a liquid receiver
    - linesets are used to calculate the refrigerant charge only
    
    """
    
   
    
    #--------------------------------------
    #--------------------------------------
    #       Expander parameters
    #--------------------------------------
    #--------------------------------------
    
    Cycle=ORCClass()
    

    exp_params = struct()
    exp_params.a_0 =-1.93529023e-03
    exp_params.a_1 = -1.25256170e-02
    exp_params.a_2 = -1.17292542e-01
    exp_params.a_3 = -4.73874764e-02
    exp_params.a_4 = 7.80410436e+00
    exp_params.a_5 = 6.36530295e-04
    exp_params.a_6 = 5.78187315e-01
    exp_params.xi = 1.23877317e+00
    exp_params.r_p_0_n = 3.07600000e+00
    exp_params.delta_n =  7.08651674e-01  
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
    Cycle.Expander.params = exp_params
    
    params={
            'N':Inputs['N_exp'],
            'Ref':Inputs['Ref'],
            'r_v_built':Inputs['r_v_built'],
            'V_disp_compressor':Inputs['V_disp_compressor'], 
            }
    Cycle.Expander.Update(**params)
    
    #--------------------------------------
    #--------------------------------------
    #      Regenerator parameters 
    #--------------------------------------
    #--------------------------------------
    params={
            'Ref_c':Inputs['Ref'],                       
            'Ref_h':Inputs['Ref'],
            #Geometric parameters 11 kW Kortrijk:  3xSWEP B200T SC-M
            'Bp' : 0.1635,
            'Lp' : 0.4485, #Center-to-center distance between ports
            'Nplates' : 150,
            'PlateAmplitude' : 0.00102, #[m]
            'PlateThickness' : 0.0003, #[m]
            'PlateWavelength' : 0.0066, #[m]
            'InclinationAngle' : pi/3,#[rad]
            'PlateConductivity' : 15.0, #[W/m-K]
            'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'

            'Verbosity':0 #1
        }
    Cycle.Regenerator.Update(**params)    


    #--------------------------------------
    #--------------------------------------
    #      Condenser parameters 
    #--------------------------------------
    #--------------------------------------
    params={
                'Ref_c':Inputs['Ref_c'],
                'Tin_c': Inputs['Tin_c'],
                'pin_c': Inputs['pin_c'],
                'Ref_h': Inputs['Ref'],
                #Geometric parameters 11 kW Kortrijk:  3xSWEP B200T SC-M
                'Bp' : 0.1635,
                'Lp' : 0.4485, #Center-to-center distance between ports
                'Nplates' : 150,
                'PlateAmplitude' : 0.00102, #[m]
                'PlateThickness' : 0.0003, #[m]
                'PlateWavelength' : 0.0066, #[m]
                'InclinationAngle' : pi/3,#[rad]
                'PlateConductivity' : 15.0, #[W/m-K]
                'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
                'Verbosity':0
            }
    Cycle.Condenser.Update(**params)
        
    #--------------------------------------
    #--------------------------------------
    #           Evaporator
    #--------------------------------------
    #--------------------------------------
    params={
                'Ref_h':Inputs['Ref_h'],
                'mdot_h':Inputs['mdot_h'],
                'Tin_h':Inputs['Tin_h'],
                'pin_h': Inputs['pin_h'],
                'Ref_c':Inputs['Ref'],
                #Geometric parameters 11 kW Kortrijk:  3xSWEP B200T SC-M
                'Bp' : 0.1635,
                'Lp' : 0.4485, #Center-to-center distance between ports
                'Nplates' : 150,
                'PlateAmplitude' : 0.00102, #[m]
                'PlateThickness' : 0.0003, #[m]
                'PlateWavelength' : 0.0066, #[m]
                'InclinationAngle' : pi/3,#[rad]
                'PlateConductivity' : 15.0, #[W/m-K]
                'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
                'Verbosity':0
            }
    Cycle.Evaporator.Update(**params)
    
    #--------------------------------------
    #--------------------------------------
    #           Pump
    #--------------------------------------
    #--------------------------------------
    #M are regression coefficients to compute the mass flow rate.  
    #'eta' are regression coefficients to compute the isentropic efficiency.
    params={
                'eta':[-1.94559985e+03,2.61004607e+03,-2.24868329e+02,-1.31580710e+03,
  -3.41646308e+03,-6.16059114e+02,1.70061718e+03,-6.12847688e+03,-3.08159095e+03,-2.13581139e+02,1.12347642e+04,  -1.91028124e+04,7.30201901e+03,2.57633083e+03,-5.00775452e+03,2.25464045e+02,-3.98897093e+02,-1.43095042e+04,   1.71628057e+04,4.38898501e+03,-3.82497930e+00,3.71939847e+02,1.77508114e+02,-5.73292533e+03,8.45941415e+03, 1.69904379e+03,5.55650945e+02,1.99756918e+04,-2.62485094e+04,-7.59467705e+02,-2.56645354e+02, -2.15969687e+01,-5.61487624e+03,2.84386468e+04,-3.23592011e+04,8.93600159e+03,-4.52678170e+03], 
                'M':[3.05754838e-01,-1.37732305e-03,-4.26385445e-07,-2.68106448e-02,1.98497578e-03],
                'Ref': Inputs['Ref'],
                'f_pp':Inputs['f_pp'],
              }
    Cycle.Pump.Update(**params)    

    #--------------------------------------
    #--------------------------------------
    #           Liquid Receiver
    #--------------------------------------
    #--------------------------------------
    params={
                'Ref': Inputs['Ref'],
              }
    Cycle.LiquidReceiver.Update(**params)

    
    #Select the cycle mode
    if Inputs['Regen'] == 'withoutRegen':
        Cycle.PreconditionedSolve(Inputs)
    elif Inputs['Regen'] == 'withRegen':
        Cycle.PreconditionedSolve_ORCRegen(Inputs)
        Cycle.ConvergencePlot()
        Cycle.BaselineTs_Celsius(Inputs)
        Write2CSV(Cycle,'ORC_Cycle.csv', append=run>0)

    print Cycle.OutputList()





if __name__=='__main__':
    
    # Single run
    """
    Run directly:
    - withRegen
    -Charge or subcooling solvers
    
    To run other combinations, the inputs of the ORC must be adjusted accordingly
    """
    #Pinch right
    Tot_Charge_list=[95]  
    f_pp_list= [35]
    Tin_w_list=[28.59]
    pin_w_list=[200.0]
    mdot_w_list=[4.166]
    mdot_s_list=[3]
    Tin_s_list=[120]
    DT_sc_list= [5]
    
    #Pinch left
#     Tot_Charge_list=[95]  
#     f_pp_list= [35]
#     Tin_w_list=[28.59]
#     pin_w_list=[200.0]
#     mdot_w_list=[4.166]
#     mdot_s_list=[1.5]
#     Tin_s_list=[120]
#     DT_sc_list= [2]
    

    run = 0 #counts the number of run we are on
    
    for f_pp,Tin_c,pin_c,mdot_h,Tin_h,mdot_c,DT_sc,Charge in zip(f_pp_list,Tin_w_list,pin_w_list,mdot_s_list,Tin_s_list,mdot_w_list,DT_sc_list,Tot_Charge_list):
    
        import time
        print "!!!!! Run# !!!!: ",run
        
        Inputs={
                'Regen': 'withRegen',       #'withRegen' or 'withoutRegen' 
                'Solver':'Subcooling',      #'Subcooling' or 'Charge'
                'Ref_c':'Water',             #Cooling medium
                'Ref_h':'INCOMP::T66',       #Hot source
                'Ref':'R245FA',              #working fluid
                'Tot_Charge':Charge,         #Tot_Charge, #working fluid total charge [kg]
                'f_pp':f_pp,                 #pump frequency [Hz]
                'N_exp':3000,                #expander rpm [rev/min]
                'r_v_built':4.8,             #expander built-in volume ratio
                'V_disp_compressor':2*6*57.68, #single-screw expander total displacement volume
                'Tin_c':Tin_c,                #cold source inlet temperature
                'pin_c':pin_c,               #cold source inlet pressure
                'mdot_c':mdot_c,             #cold source mass flow rate
                'Tin_h':Tin_h,               #hot source inlet temperature
                'mdot_h':mdot_h,             #hot source mass flow rate
                'pin_h':140,                #hot source inlet pressure
                'DT_sc':DT_sc               #subcooling at the condenser outlet guess to solve pre-cycle calculation
                }
        time0 = time.time()
        run += 1

        print
        ORCSystem_Kortrijk(Inputs, run)
        time1 = time.time() 
        print 'Run', run,'executed in', (time1 - time0)/60., 'minutes.'
        print           

