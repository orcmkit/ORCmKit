from __future__ import division
import numpy as np
from scipy.linalg import inv
from numpy import array
import matplotlib
import pylab

from matplotlib.pyplot import plot, show, figure, semilogy, xlim, ylim, title, xlabel, ylabel, legend
import matplotlib.pyplot as plt


from counter import my_counter_2

def MultiDimNewtRaph(f,x0,dx=1e-6,args=(),ytol=1e-4,w=1.0,JustOneStep=False):
    """
    A Newton-Raphson solver where the Jacobian is always re-evaluated rather than
    re-using the information as in the fsolve method of scipy.optimize
    """
    x=np.array(x0)
    error=999
    J=np.zeros((len(x),len(x)))
    
    error_list=[]
    iteration_list = []
    
    #If a float is passed in for dx, convert to a numpy-like list the same shape
    #as x
    if isinstance(dx,float):
        dx=dx*np.ones_like(x)
        
    r0=array(f(x,*args))
    while abs(error)>ytol:
        #Build the Jacobian matrix by columns
        for i in range(len(x)):
            epsilon=np.zeros_like(x)
            epsilon[i]=dx[i]
            J[:,i]=(array(f(x+epsilon,*args))-r0)/epsilon[i]
        v=np.dot(-inv(J),r0)
        x=x+w*v
        #Calculate the residual vector at the new step
        r0=f(x,*args)
        error = np.max(np.abs(r0))
        iteration=my_counter_2.next()
        
        error_list.append(error)
        iteration_list.append(iteration)
        #print error
        #Just do one step and stop
        if JustOneStep==True:
            return x
    
    error_convergence = np.array(error_list)
    iteration_convergence = np.array(iteration_list)
    print '---MultiDimNewtRaph---'
    print iteration_convergence
    print error_convergence
    matplotlib.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig=figure(figsize=(8,6))
    ax=fig.add_subplot(1,1,1)

    plot(iteration_convergence, error_convergence,'bs-',linewidth=2.0,markersize = 10,markerfacecolor="blue", markeredgewidth=2, markeredgecolor="black")
    plt.axhline(y=1e-5,color='r',ls='dashed', linewidth=2.5)
    xlabel('Iteration [-] ',fontsize=18)
    ylabel(r'$|Resid|$',fontsize=18)
    title('Overall Cycle Convergence',fontsize=18)
    
    plt.yscale('log')
    plt.minorticks_off() 
    for tickx in ax.xaxis.get_major_ticks():
        tickx.label.set_fontsize(20)
    for ticky in ax.yaxis.get_major_ticks():
        ticky.label.set_fontsize(20)
    
    show()
    fig.savefig('MultiDimNewtRaph.png',dpi=300)
    return x
        
    