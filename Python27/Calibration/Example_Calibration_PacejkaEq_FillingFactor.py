"""
ORCmKit : Organic Rankine Cycle Modelling Kit

Example of calibration of Pacejka equation and filling factor with set of experimental data

Packages required:
- Coolprop


@Author: Davide Ziviani, Ghent University
Version: 1.0
Date: 02/21/2016

"""

from __future__ import division
from math import pi, cos, sin, tan, sqrt, asin, acos, atan, fabs, log,exp
from scipy import integrate
from scipy.optimize import leastsq, curve_fit, minimize
from scipy.stats import chi2
import numpy as np
from numpy import *
import pylab
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

import matplotlib
from matplotlib.pyplot import plot, show, figure, semilogy, xlim, ylim, title, xlabel, ylabel, legend,cm
import matplotlib.mlab as ml

import xlrd

import CoolProp 
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropsPlot



class PacejkaFillingFactorClass():
    """
    Calibration of isentropic efficiency and filling factor
    from experimental data

    """
    
    def __init__(self,**kwargs):
        #Load up the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)
    
    def Update(self,**kwargs):
        #Update the parameters passed in
        # using the dictionary
        self.__dict__.update(kwargs)

        
    def R_squared(self,y_calc,y_meas):
        """
        R^2 = 1 - SS_err/SS_tot
        """
        y_bar = np.mean(y_meas)

        SSTot = np.sum((y_meas - y_bar)**2)
        SSRes= np.sum((y_meas - y_calc)**2)
        Rsquared = 1-(SSRes/SSTot)
        
        return Rsquared
        
    def MAPE(self,y_true,y_pred):
        
        return np.mean(np.abs((y_true - y_pred) / y_true)) * 100
        
        
    def Pacejka(self,data,a0,a1,a2,a3,a4,a5,a6,shape_0,dydx_0,x_0_0,x_m_0,y_m_0,N_exp_m):   #x_0_0,dydx_0,shape_0,x_m_0,y_m_0,N_exp_m):
        """
        data[0,:] : pressure ratio [-]
        data[1,:] : expander rotational speed [rpm]
        data[2,:] : expander inlet pressure [kPa]
        """
        
        rp,Nexp,p= data[0,:],data[1,:],data[2,:]
        N_exp_star_m = (self.N_exp_m -3000)/3000    
        p_star = (p - 10)/10
        N_star = (Nexp - 3000) /3000

        
        x_0=self.x_0_0 + a0*Nexp        
        dydx= dydx_0 + a1*p_star - a2*N_star
        shape=shape_0
        x_m=self.x_m_0 -  a3*p_star + a4*N_star
        y_max=self.y_m_0 +a5*p_star -a6*(N_star - N_exp_star_m)**2
    
        "Equation of Pacejka"
        A = x_0
        C = shape
        D = y_max
        B = dydx/(C*D)        
        """
        Note: Non-numpy functions like math.abs() or math.log10() don't play nicely with numpy arrays
        (similar atan). In order to avoid the following error while working with numpy.array
        
        TypeError: only length-1 arrays can be converted to Python scalars
        
        math.atan  ---> np.arctan
        
        """
        E = (B*(x_m-x_0)-tan(pi/(2*C)))/(B*(x_m-x_0)-np.arctan(B*(x_m-x_0)))

        return D * np.sin(C * np.arctan(B*(rp-A) - E * (B * (rp-A) - np.arctan(B*(rp-A)))))

    def Pacejka_residuals(self,params,y,data):

        
        err = (y - self.Pacejka(data,*params))
        
        return err
    
    def Pacejka_residuals_min(self,params,y,data):

        
        err = ((y - self.Pacejka(data,*params))**2).sum()


        return err
    
    def FillingFactor(self,x,FF1,FF2,FF3,FF4):

        rp_star=self.p_su_exp/self.p_ex_exp
        p_star = (self.p_su_exp - 1000)/1000
        A=FF3*rp_star
        B=FF4*p_star
        
        return  FF1 + FF2*np.log(x/3000)+A+B


    
    def Calculate_etais(self):
        
        #Assign experimental data
        p_su_exp = self.p_su_exp
        p_ex_exp = self.p_ex_exp
        T_su_exp = self.T_su_exp
        T_ex_exp = self.T_ex_exp
        rho_su=PropsSI('D','T',T_su_exp+273.15,'P',p_su_exp,self.Ref)
        rho_ex=PropsSI('D','T',T_ex_exp+273.15,'P',p_ex_exp,self.Ref)
        rp_exp=p_su_exp/p_ex_exp
        eta_is_exp = self.eta_is_exp
        N_exp = self.N_exp
        FF_exp = self.FF_exp

        """
        Pacejka Equation interpolation
        """
        # Initial guess for parameters
        Pacejka_coeff_guess = [self.a0,self.a1,self.a2,self.a3,self.a4,self.a5,self.a6,self.shape_0,self.dydx_0,self.x_0_0,self.x_m_0,self.y_m_0,self.N_exp_m] 
        
        # Fit equation using least squares optimization
        xy = np.array([rp_exp,N_exp,p_su_exp])

        #Set the bounds
        bnds = ((-1, 10), (-10, 2),(-1, 10),(-1, 2),(-1, 10),(-10, 10),(-10, 10),(0,10),(-10,3),(1,10),(1,10),(0.2,0.6),(2000,4000))

        #Run minimization
        optparam = minimize(self.Pacejka_residuals_min,Pacejka_coeff_guess,args=(eta_is_exp,xy), method='SLSQP',bounds = bnds ,tol= 1e-10,options={'maxiter':10000 , 'disp': True})
        
        print 'Opt Parameters SLSQP:', optparam
        
        #Calculate isentropic efficiency
        eta_is_exp_fit_SLSQP = self.Pacejka(xy,*optparam.x)
        
        
        #Coefficient of Determination,RE,MAPE
        Rsquared_eta_is = self.R_squared(eta_is_exp_fit_SLSQP,eta_is_exp)
        
        MAPE_eta_is = self.MAPE(self.eta_is_exp,eta_is_exp_fit_SLSQP)
        RE_eta_is = np.max(np.fabs(self.eta_is_exp - eta_is_exp_fit_SLSQP)/self.eta_is_exp)*100
        print 'MAPE [%]:',MAPE_eta_is
        print 'RE max [%]:',RE_eta_is    
        print 'Rsquared eta_is_exp [SLSQP]:',Rsquared_eta_is
        

        
        """
        Example of Parity plot
        """
        fig = plt.figure(figsize=(9,8))
        ax = fig.add_subplot(111)

        plt1=plt.scatter(self.eta_is_exp,eta_is_exp_fit_SLSQP,s=150,facecolors='blue',alpha= 0.2,marker ='s', edgecolors='k',linewidths=2)       
        plt.plot(np.linspace(0.0,0.9,100),np.linspace(0,0.9,100),'r-',linewidth = 2)
        plt.plot(np.linspace(0.0,0.9,100),(1+RE_eta_is/100 )*np.linspace(0.0,0.9,100),'--b',linewidth = 3)
        plt.plot(np.linspace(0.0,0.9,100),(1-RE_eta_is/100 )*np.linspace(0.0,0.9,100),'--b',linewidth = 3)

        plt.text(0.2,0.50, r'RE$_{max}$ = $\pm$ %1.2f %%' % RE_eta_is,color='blue',fontsize=25)
        plt.text(0.2,0.45,r'MAPE = %1.2f %% ' %MAPE_eta_is,color = 'blue',fontsize = 25)
        plt.xlim(0.15,0.55)
        plt.ylim(0.15,0.55)
        plt.xlabel(r'$\varepsilon_{is,oa,exp,meas} \, [-]$',fontsize=25)
        plt.ylabel(r'$\varepsilon_{is,oa,exp,meas} \, [-]$',fontsize=25)

        for tickx in ax.xaxis.get_major_ticks():
            tickx.label.set_fontsize(25)
        for ticky in ax.yaxis.get_major_ticks():
            ticky.label.set_fontsize(25)

        fig.savefig('ParityPlot_etais.png',dpi=600)




        
    def Calculate_FF(self):
        """
        Fit the filling factor
        """    
        N_exp = self.N_exp
        FF_exp = self.FF_exp
        
        # Initial guess for parameters
        FF_coeff_guess = [self.FF1,self.FF2,self.FF3,self.FF4]

        #Run calibration
        params_FF, params_covariance_FF = curve_fit(self.FillingFactor,N_exp,FF_exp,p0=FF_coeff_guess)

        FF_exp_fit = self.FillingFactor(N_exp,*params_FF)
        
        #Calculate errors
        MAPE_FF = self.MAPE(FF_exp_fit,FF_exp)
        Rsquared_FF = self.R_squared(FF_exp_fit,FF_exp)*100
        RE_FF = np.max(np.fabs(FF_exp - FF_exp_fit)/FF_exp)*100
        print 'MAPE_FF [%]:',MAPE_FF
        print 'Rsquared_FF [%]:',Rsquared_FF
        print 'RE_FF [%]:', RE_FF
        
        
        """
        Example of parity plot
        """
        
        #Expander mass flowrate"
        rho_su_exp=PropsSI('D','P',self.p_su_exp*1000,'T',self.T_su_exp+273.15,self.Ref)
        V_s_exp=self.V_s_comp/self.rv_in
        V_dot_s_exp=2*6*V_s_exp*(N_exp/60)
        m_dot_exp_calc=FF_exp_fit*V_dot_s_exp*rho_su_exp

        figFF=figure(figsize=(8,6))
        ax=figFF.add_subplot(1,1,1)
        plt.plot(self.m_dot_exp,m_dot_exp_calc,'ob',markersize = 15,markerfacecolor="blue", markeredgewidth=2, markeredgecolor="black",label='Fit')
        plt.plot(np.linspace(0.0,0.9,100),np.linspace(0.0,0.9,100),'--k',linewidth = 2)
        
        xlabel(r'$\dot{m}_{exp,data} [kg/s]$',fontsize=18)
        ylabel(r'$\dot{m}_{exp,fit} [kg/s]$',fontsize=18)
        for tickx in ax.xaxis.get_major_ticks():
            tickx.label.set_fontsize(20)
        for ticky in ax.yaxis.get_major_ticks():
            ticky.label.set_fontsize(20)
        plt.xlim(0.10,0.50)
        plt.ylim(0.10,0.50)            
        figFF.savefig('ParityPlot_mdot.png',dpi=600)

        
    
                
        
if __name__=='__main__':        

    def DataIO(start,end,filename,sheet_num):
        file = xlrd.open_workbook(filename)
        sheet = file.sheet_by_index(sheet_num)
        ncol = sheet.ncols 
        i = 0
        data = []
        #print range(ncol)
        for i in range(ncol-1):
            col_values = sheet.col_values(colx=i, start_rowx=start, end_rowx=end)
            data.append(col_values)
        return data
    
    #--------------------------------------------------------------------------
    def Import(start,end,filename,sheet_num):
        
        """
        Data Import from Excel    
        """
        
        data = DataIO(start,end,filename,sheet_num)
        i = 0  
    
        while i < (end - start+1):
            W_dot_exp = np.asarray(data[5])     #[W] 
            rpm_exp = np.asarray(data[4])       #[rpm]
            m_dot = np.asarray(data[7])         #[kg/s]
            p_su_exp = np.asarray(data[1])      #[Pa]
            p_ex_exp = np.asarray(data[2])      #[Pa] 
            T_su_exp= np.asarray(data[8])       #[C]
            T_ex_exp= np.asarray(data[9])       #[C]
            eta_is = np.asarray(data[6])        #[-]

            i=i+1
            Data = [W_dot_exp,rpm_exp,m_dot,p_su_exp,p_ex_exp,T_su_exp,T_ex_exp,eta_is]
        
        return Data
        
    

    #Import Experimental Data"
    start=2
    end=44
    filename = 'Example_ExpData_R245fa.xlsx'
    [W_dot_exp,rpm_exp,m_dot,p_su_exp,p_ex_exp,T_su_exp,T_ex_exp,eta_is_exp] = Import(start,end,filename,sheet_num = 0)
    
    fluid = 'R245FA'
    
    rho_su_exp = PropsSI('D','P',p_su_exp,'T',T_su_exp+273.15,fluid)
    rho_ex_exp = PropsSI('D','P',p_ex_exp,'T',T_ex_exp+273.15,fluid)
    v_su_exp = 1/rho_su_exp
    v_ex_exp = 1/rho_ex_exp    
    r_p_exp = p_su_exp/p_ex_exp
    r_v_exp = v_ex_exp/v_su_exp
    
    V1 = 57.39e-06  #[m3]
    rv_in = 5
    FF = (m_dot* v_su_exp) / (2*6*V1/rv_in*rpm_exp/60)

    kwds={  
            'Ref':'R245fa',
            'T_ex_exp':  T_ex_exp , 
            'T_su_exp': T_su_exp,
            'p_su_exp': p_su_exp/1000, #[kPa] 
            'p_ex_exp': p_ex_exp/1000, #[kPa]
            'eta_is_exp': eta_is_exp, 
            'N_exp':  rpm_exp,
            'FF_exp':FF, 
            'm_dot_exp':m_dot,
            'rv_in': rv_in,
            'V_s_comp':V1,#[m3]
            #"Expander parameters:"
            'x_0_0':3.07600000e+00,     #"Value of the pressure ratio when epsilon=0, in the reference conditions"
            'dydx_0':7.08651674e-01,    #Slope at epsilon=0, in the reference conditions"
            'shape_0':1.23877317e+00,   #"Shape parameter of the epsilon vs. rp curve"
            'x_m_0':5.99558088e+00,     #  "rp corresponding the maximum epsilon. in the reference conditions"
            'y_m_0':5.19000000e-01,     #"Max epsilon in the reference conditions"
            'N_exp_m':3547,             #  "rpm corresponding to epsilon_max"
            'a0':-1.93529023e-03, 
            'a1':-1.25256170e-02, 
            'a2':-1.17292542e-01, 
            'a3':-4.73874764e-02, 
            'a4':7.80410436e+00, 
            'a5':6.36530295e-04, 
            'a6':5.78187315e-01, 
            'FF1':1.489, 
            'FF2':0.3427, 
            'FF3':0.8187 ,
            'FF4':0.1435
            }
    
    Exp=PacejkaFillingFactorClass(**kwds)
    Exp.Calculate_etais()
    Exp.Calculate_FF()




















