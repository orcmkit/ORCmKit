"""
ORCmKit : Organic Rankine Cycle Modelling Kit

Example of calibration of Pacejka equation with set of experimental data

Packages required:
- Coolprop
- Mystic


@Author: Davide Ziviani, Ghent University
Version: 1.0
Date: 02/21/2016

"""

from __future__ import division

#import Coolprop
import CoolProp
from CoolProp.CoolProp import PropsSI
#Import mystic
from mystic.solvers import DifferentialEvolutionSolver 
from mystic.termination import ChangeOverGeneration, VTR
from mystic.strategy import Best1Exp, Rand1Exp
from mystic.monitors import VerboseMonitor
from mystic.tools import random_seed,getch
random_seed(123)
#Python module imports
from math import pi, cos, sin, tan, sqrt, asin, acos, atan, fabs, log,exp
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, figure, semilogy, xlim, ylim, title, xlabel, ylabel, legend,cm
from mpl_toolkits.mplot3d import Axes3D
import pylab
import scipy.optimize
import xlrd
import numpy as np
from mystic.models.abstract_model import AbstractFunction


class CalibrationPacejkaEq(AbstractFunction):
    
    def __init__(self, x, y, nnum=13, nden=1):
        """
        Class initialization:
        
        nnum : number of coefficients
        nden : 
        
        x : array of experimental data 
        y : array of experimental isentropic efficiency
        
        """
        self.rp = x[0,:]
        self.Nexp = x[1,:]
        self.p = x[2,:]
        self.y = y

        self.nnum = nnum
        self.nden = nden-1
        
        AbstractFunction.__init__(self, ndim = self.nnum + self.nden)


    def eval(self, coeffs, x = None):
        if x is None:
            rp = self.rp
            p = self.p
            Nexp = self.Nexp
        else:
            rp = x[0,:]
            Nexp = x[1,:]
            p = x[2,:]
        assert(len(coeffs) == self.nnum+self.nden)

        a0,a1,a2,a3,a4,a5,a6,shape_0,dydx_0,x_0_0,x_m_0,y_m_0,N_exp_m = coeffs 
        
        #Set reference values for mornalized variables
        Nexp_ref = 3500 #[rpm] 
        p_ref = 1000 #[kPa]
        N_exp_star_m = (N_exp_m -Nexp_ref)/Nexp_ref   
        p_star = (p - p_ref)/p_ref               
        Nexp_star = (Nexp - Nexp_ref)/Nexp_ref        

        x_0=x_0_0 + a0*Nexp_star        
        dydx= dydx_0 + a1*p_star - a2*Nexp_star
        shape=shape_0
        x_m=x_m_0 -  a3*p_star + a4*Nexp_star
        y_max=y_m_0 +a5*p_star -a6*(Nexp_star - N_exp_star_m)**2
    
        "Equation of Pacejka"
        A = x_0
        C = shape
        D = y_max
        B = dydx/(C*D)        
        E = (B*(x_m-x_0)-tan(pi/(2*C)))/(B*(x_m-x_0)-np.arctan(B*(x_m-x_0)))

        return D * np.sin(C * np.arctan(B*(rp-A) - E * (B * (rp-A) - np.arctan(B*(rp-A)))))
        
    def function(self, coeffs):
        """
        Residual function
        """
        err = self.eval(coeffs) - self.y
        rsse = np.sqrt(np.sum(np.power(err, 2)))

        return rsse

    #minimizers = [0.] #XXX: there are many periodic local minima
    
def R_squared(y_calc,y_meas):
    """
    R^2 = 1 - SS_err/SS_tot
    """
    y_bar = np.mean(y_meas)

    SSTot = np.sum((y_meas - y_bar)**2)
    SSRes= np.sum((y_meas - y_calc)**2)
    Rsquared = 1-(SSRes/SSTot)
    
    return Rsquared
        
def MAPE(y_true,y_pred):
    """
    Function returning the Mean Absolute Percentage Error
    """

    return np.mean(np.abs((y_true - y_pred) / y_true))*100


def DataIO(start,end,filename,sheet_num):
    """
    Function to read Excel file
    """
    
    file = xlrd.open_workbook(filename)
    sheet = file.sheet_by_index(sheet_num)
    ncol = sheet.ncols 
    i = 0
    data = []

    for i in range(ncol-1):
        col_values = sheet.col_values(colx=i, start_rowx=start, end_rowx=end)
        data.append(col_values)
    return data

#--------------------------------------------------------------------------
def Import(start,end,filename,sheet_num):
    
    "Function to import experimental data"
    
    data = DataIO(start,end,filename,sheet_num)
    i = 0  

    while i < (end - start+1):
        Ref = np.asarray(data[0])
        psu_exp = np.asarray(data[1])/1000    #[kPa]
        rp = np.asarray(data[3])              #[-]
        rpm_exp = np.asarray(data[4])         #[rpm]
        Wdot_exp = np.asarray(data[5])        #[W]            
        eta_s = np.asarray(data[6])           #[-]
        i=i+1
        Data = [Ref,psu_exp,rp,rpm_exp,Wdot_exp,eta_s]
    
    return Data


def parity_plot(eta_is_exp_fit,eta_is_exp,Ref):
    
    RE_eta_is = np.max(np.fabs(eta_is_exp - eta_is_exp_fit)/eta_is_exp)*100
    Rsquared_eta_is = R_squared(eta_is_exp_fit,eta_is_exp)
    MAPE_eta_is = MAPE(eta_is_exp_fit,eta_is_exp)
    
    print 'RE [%]:',RE_eta_is      
    print 'MAPE_eta [%]:', MAPE_eta_is
    print 'Rsquared eta_is_exp:',Rsquared_eta_is 
    
    f1=pylab.figure(figsize=(4,4))
    ax1=f1.add_axes((0.15,0.15,0.8,0.8))

    pylab.plot(eta_is_exp[Ref=='R245FA'],eta_is_exp_fit[Ref== 'R245FA'],'bo', ms = 8,label = 'R245fa')                          
    pylab.plot(np.linspace(0.0,0.9,100),np.linspace(0,0.9,100),'k-',linewidth = 2)        
    pylab.plot(np.linspace(0.0,0.9,100),(1+RE_eta_is/100 )*np.linspace(0.0,0.9,100),'--k',linewidth = 1)
    pylab.plot(np.linspace(0.0,0.9,100),(1-RE_eta_is/100 )*np.linspace(0.0,0.9,100),'--k',linewidth = 1)
    pylab.text(0.45,0.25, r'MAPE = %1.2f %% ' %MAPE_eta_is,color = 'black',fontsize=12)
    pylab.text(0.50,0.4,r'$\pm$ %1.2f %%' % RE_eta_is ,color='black',fontsize=12)
    pylab.xlim(0.2,0.7)
    pylab.ylim(0.2,0.7)
    f1.subplots_adjust(bottom=0.15,left= 0.2)
    ax1.set_xlabel('$\epsilon_{is,oa,exp,meas}$ [-]')
    ax1.set_ylabel('$\epsilon_{is,oa,exp,meas}$ [-]')
    f1.savefig('PacejkaEq_R245fa.png',dpi=600)
    pylab.show()



def main(start,end,filename):
    
    #Import Experimental Data
    [Ref,p_su_exp,rp_exp,N_exp,Wdot_exp,eta_is_exp] = Import(start,end,filename,sheet_num = 0)

    data = np.array([rp_exp,N_exp,p_su_exp])
    
    #Set solver
    ND = 13
    NP = ND*10
    MAX_GENERATIONS = 3000
    
    minrange = [-10,-100,-10,-10,-10,-10,-10,0,0,-100,-10,0,0]
    maxrange = [10,1,1,10,1,1,10,10,10,5,10,0.8,5000]

    solver = DifferentialEvolutionSolver(ND, NP)
    solver.SetRandomInitialPoints(min = [0.1]*ND, max = [5]*ND)
    solver.SetStrictRanges(min=minrange, max=maxrange)
    solver.SetEvaluationLimits(generations=MAX_GENERATIONS)    
    
    pf = CalibrationPacejkaEq(data, eta_is_exp, nnum = 13, nden = 1)

    solver.Solve(pf.function, termination=VTR(1e-8), strategy=Rand1Exp,\
                 CrossProbability=0.9, ScalingFactor=0.9)

    coeff_solution = solver.Solution()
    
    print 'DE coefficients:', coeff_solution
    
    eta_is_exp_fit = pf.eval(coeff_solution)
    parity_plot(eta_is_exp_fit,eta_is_exp,Ref)
    
    return pf.eval(coeff_solution)
    






if __name__ == '__main__':

    #Optimize with DESolver
    start=2
    end=44
    filename = 'Example_ExpData_R245FA.xlsx'
    main(start,end,filename) 
    
