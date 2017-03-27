from __future__ import division
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import plotly.plotly as py
from plotly import tools
from plotly.graph_objs import *
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
init_notebook_mode()

my_dpi = 100

def dataEV(c_L, N, l_x, l_y, itersRinfo = True):
    #import energy data
    yE = np.genfromtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                    str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E.txt')
    avgE = np.mean(yE, axis=0) #average energy data
    varssE = np.var(yE[:,-1], ddof = 0)/(c_L**2) #variance of last column gives variance of steady value
    xE = np.arange(0,len(avgE))
    
    #import vortices data
    yV = np.genfromtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" 
                    + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.txt')
    avgV = np.mean(yV, axis=0)
    xV = np.arange(0,len(avgV))
    
    if itersRinfo:
        print "cL: ", c_L , " E_R, E_iterations: ", yE.shape , " V_R, V_iterations: ", yV.shape
        
    return xE, avgE, xV, avgV, varssE

def plotE(xE, avgE, c_L, N, l_x, l_y, save = False):
    plt.figure()
    plt.plot(xE,avgE)
    #plt.title(r'Average energy density over time for $c_L =  %s$'%(c_L) + '\n' + r'and lattice size $%s^2$' %(N)
              #+ r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
    plt.xlabel(r'$t$', fontsize=16)
    plt.ylabel(r'$\langle U \rangle$', fontsize=16)
    
    if save:
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                        str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E.jpg', bbox_inches='tight', pad_inches=0.1)
    
def plotV(xV, avgV, c_L, N, l_x, l_y, save = False):
    plt.figure()
    plt.plot(xV,avgV)
    #plt.title(r'Average number of vortices over time for $c_L =  %s$'%(c_L) + '\n' + r'and lattice size $%s^2$' %(N)
              #+ r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
    plt.xlabel(r'$t$', fontsize=16)
    plt.ylabel(r'$\langle V \rangle$', fontsize=16)
    
    if save:
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                        str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.jpg', bbox_inches='tight', pad_inches=0.1)

def allEV(c_Lplot, N, l_x, l_y, display = True, save = False, itersRinfo = True, All = False):
    """Function to return average energy and vortices values at steady state for all cL values given lattice size, lambda_x,
    lambda_y with option to plot average energy and vortices against time, to display and/or to save plots"""
    dataE = []
    dataV = []
    avgssE = np.zeros(len(c_Lplot))
    avgssV = np.zeros(len(c_Lplot))
    CVs = np.zeros(len(c_Lplot))
    
    if All:
        fig = plt.figure(figsize=(5,8))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        
        #ax1.set_title(r'Lattice size $%s^2$'%(N) + r' with $\lambda_x = %s$'%(l_x) + r' and $\lambda_y = %s$'%(l_y)
                      #+ '\n' + r'Average energy density over time', fontsize='medium')
        ax1.set_xlabel('$t$', fontsize='small')
        ax1.set_ylabel(r'$\langle U_{ss} \rangle$', fontsize='small')

        #ax2.set_title(r'Average number of vortices over time', fontsize='medium')
        ax2.set_xlabel('$t$', fontsize='small')
        ax2.set_ylabel('$\langle V_{ss} \rangle$', fontsize='small')
        
        color=plt.cm.rainbow(np.linspace(0,1,len(c_Lplot)))
        plt.tight_layout()
       
    for i in range(0,len(c_Lplot)): 
        xE, avgE, xV, avgV, varssE = dataEV(c_Lplot[i],N,l_x,l_y, itersRinfo)
        dataE.append(Scatter(x=xE,y=avgE, name = 'c_L = %s'%(c_Lplot[i])))
        dataV.append(Scatter(x=xV,y=avgV, name = 'c_L = %s'%(c_Lplot[i])))
        avgssE[i] = avgE[-1] #steady state is the last value of average values over realisations
        avgssV[i] = avgV[-1]
        CVs[i] = varssE/(c_Lplot[i]**2) #specific heat for each c_L
        
        if save:
            if All:
                ax1.plot(xE, avgE, label = 'c_L = %s'%(c_Lplot[i]), c=color[i])
                ax2.plot(xV, avgV, label = 'c_L = %s'%(c_Lplot[i]), c=color[i])
                plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5, fontsize = 'x-small')
                plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'allEV.jpg', bbox_inches='tight', pad_inches=0.1)
                
            else:
                plotE(xE, avgE, c_Lplot[i], N, l_x, l_y, save = True)
                plotV(xV, avgV, c_Lplot[i], N, l_x, l_y, save = True)
        
    if display:
        layoutE = Layout(
        title = 'Average energy against time',
        xaxis = dict(title = 'Time'),
        yaxis = dict(title = 'Average energy'),
        width = 700)

        layoutV = Layout(
        title = 'Average number of vortices against time',
        xaxis = dict(title = 'Time'),
        yaxis = dict(title = 'Average number of vortices'),
        width = 700)

        iplot(dict(data=dataE,layout=layoutE))
        iplot(dict(data=dataV,layout=layoutV))

    
    return avgssE, avgssV, CVs

def plotavgE(c_Lplot, avgssE, N, l_x, l_y, display = True, save = False):
    """Plot average steady state energy given c_L and energy values, N, l_x, l_y with option to display and/or to show plots"""
    if save:
        plt.figure(figsize=(550/my_dpi, 330/my_dpi), dpi=my_dpi)
        plt.plot(c_Lplot,avgssE) #"o" for scatter
        #plt.title(r'Average energy density against noise parameter for lattice size $%s^2$' %(N) + '\n' +
                  #r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel(r'$c_L$', fontsize=16)
        plt.ylabel(r'$\langle U_{ss} \rangle$', fontsize=16)
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-"  + str(format(l_x, '.2f')) + "-" +
                    str(format(l_y, '.2f')) + 'avgE.jpg', bbox_inches='tight', pad_inches=0.1)
    
    if display:
        line = Scatter(x=c_Lplot,y=avgssE, mode = "lines+markers")
        data = [line]

        layout = Layout(
            title = 'Average energy density vs noise',
            xaxis = dict(title = 'Noise'),
            yaxis = dict(title = 'Average energy'),
            width = 700)

        iplot(dict(data=data,layout=layout))
        
def plotavgV(c_Lplot, avgssV, N, l_x, l_y, display = True, save = False):
    """Plot average steady state energy given c_L and energy values, N, l_x, l_y with option to display and/or to show plots"""
    if save:
        plt.figure(figsize=(550/my_dpi, 330/my_dpi), dpi=my_dpi)
        plt.plot(c_Lplot,avgssV)
        #plt.title(r'Average number of vortices against noise parameter for lattice size $%s^2$' %(N)+ '\n' + r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel(r'$c_L$', fontsize=16)
        plt.ylabel(r'$\langle V_{ss} \rangle$', fontsize=16)
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-"  + str(format(l_x, '.2f')) + "-" +
                                str(format(l_y, '.2f')) + 'avgV.jpg', bbox_inches='tight', pad_inches=0.1)
        
        plt.figure(figsize=(550/my_dpi, 330/my_dpi), dpi=my_dpi)
        plt.plot(c_Lplot,avgssV/(N**2))
        #plt.title(r'Average number of vortices against noise parameter for lattice size $%s^2$' %(N)+ '\n' + r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel(r'$c_L$', fontsize=16)
        plt.ylabel(r'$\langle V_{ss} \rangle/N^{2}$', fontsize=16)
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-"  + str(format(l_x, '.2f')) + "-" +
                    str(format(l_y, '.2f')) + 'avgVdens.jpg', bbox_inches='tight', pad_inches=0.1)
        
    if display:
        line = Scatter(x=c_Lplot,y=avgssV, mode = "lines+markers")
        data = [line]

        layout = Layout(
            title = 'Average number of vortices vs noise',
            xaxis = dict(title = 'Noise'),
            yaxis = dict(title = 'Average number of vortices'),
            width = 700)

        iplot(dict(data=data,layout=layout))

def plot_specheatdiff(c_Lplot, avgssE, N, l_x, l_y, display = True, save = False):
    """Function to calculate dE/dcL by differentiation returning new cL values and dE/dcL values with option to
    plot and/or save graphs"""
    dcL = np.diff(c_Lplot)
    c_Lnew = c_Lplot[:len(dcL)] + dcL*0.5
    dEdcL = np.diff(avgssE)/dcL
    
    i_max = np.argmax(dEdcL)
    print c_Lnew[i_max]
    
    if save:
        plt.figure(figsize=(550/my_dpi, 330/my_dpi), dpi=my_dpi)
        plt.plot (c_Lnew,dEdcL) #"o" for scatter
        #plt.title('Change in average energy density over change in noise \n' + r'against noise for lattice size $%s^2$' %(N)+
                  #r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel(r'$c_L$', fontsize=16)
        plt.ylabel(r'$d\langle U_{ss} \rangle/dc_L$', fontsize=16)
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-"  + str(format(l_x, '.2f')) + "-" +
                    str(format(l_y, '.2f')) + 'specheatdiff.jpg', bbox_inches='tight', pad_inches=0.1)

    if display:
        line = Scatter(x=c_Lnew,y=dEdcL, mode = "lines+markers")
        data = [line]

        layout = Layout(
            title = 'Specific heat using differentiation',
            xaxis = dict(title = 'Noise'),
            yaxis = dict(title = 'Change in energy/change in noise'),
            width = 700)

        iplot(dict(data=data,layout=layout))

    return c_Lnew, dEdcL


def plotVquench(cL_i, cL_f, Gamma_Qs, N, l_x, l_y, display = False, save = False):
    """Function to plot number of vortices against time for different quench rates with option to display and/or save graphs
    and returns the steady state number of vortices for each quench rate"""
    data = []
    logVss = np.zeros(len(Gamma_Qs))
    logtGamma = np.zeros(len(Gamma_Qs))
    
    plt.figure(figsize=(550/my_dpi, 330/my_dpi), dpi=my_dpi)
    plt.xlabel(r'$t$', fontsize=16)
    plt.ylabel(r'$\langle V \rangle$', fontsize=16)  
    
    for i in range(0,len(Gamma_Qs)): 
        yV = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(cL_i, '.2f')) + "-" + str(format(cL_f, '.2f')) + "-"
        + str(Gamma_Qs[i]) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.txt')
        avgV = np.mean(yV, axis=0) #average number of vortices data
        xV = np.arange(0,len(avgV))
        data.append(Scatter(x=xV,y=avgV, name = 'Gamma_Q = %s'%(Gamma_Qs[i])))
        logtGamma[i] = np.log2(len(avgV)/avgV[-1])
        logVss[i] = np.log2(avgV[-1])
   
        
        if save:
            plt.plot(xV,avgV,label = r'$\tau_Q = %s$'%(Gamma_Qs[i]), linewidth=1)
            plt.legend()
            plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'Vtquench.jpg', bbox_inches='tight', pad_inches=0.1)
        
    if display:
        layout = Layout(
        title = 'Average number of vortices against time',
        xaxis = dict(title = 'Time'),
        yaxis = dict(title = 'Number of vortices'),
        width = 700)

        iplot(dict(data=data,layout=layout))
        
    return logVss, logtGamma

def plotVquenchvscL(cL_i, cL_f, c_Lplot, avgssV, Gamma_Qs, N, l_x, l_y, display = False, save = False, linewidth = 2 ):
    """Function to plot number of vortices against cL for different quench rates with option to display and/or save graphs"""
    data = []
    cLs = []
    
    cLf_index = int(np.where(c_Lplot==cL_f)[0])
    
    if save:
        fig, ax = plt.subplots()
        ax.set_xlabel(r'$c_L$', fontsize=16)
        ax.set_ylabel(r'$\langle V_{ss} \rangle$', fontsize=16)
        color=plt.cm.rainbow(np.linspace(0,1,len(Gamma_Qs)+1))
        plt.tight_layout()
        
    for i in range(0,len(Gamma_Qs)): 
        yV = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(cL_i, '.2f')) + "-" + str(format(cL_f, '.2f')) + "-"
        + str(Gamma_Qs[i]) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.txt')
        avgV = np.mean(yV,axis = 0)
        
        #calculate cL values
        cL = np.zeros(len(avgV))
        cL[0] = cL_i
        dcL = (cL_f - cL_i)/Gamma_Qs[i]
        cL_current = cL_i
        for t in range(1, Gamma_Qs[i]+1):
            cL_current += dcL
            cL[t] = cL_current
            
        cL[Gamma_Qs[i]+1:] = cL_f
        
        cLs.append(cL)
        
        data.append(Scatter(x=cL,y=avgV, name = 'Gamma_Q = %s'%(Gamma_Qs[i]), line = dict(width = linewidth)))
    
        if save:
            ax.plot(cL, avgV, label = r'$\tau_Q$ = %s'%(Gamma_Qs[i]), c = color[i])
    
    if save:
        ax.plot(c_Lplot[cLf_index:],avgssV[cLf_index:], c = color[-1], label = r'non-quenching')
        plt.legend(fontsize = 'small')
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + str(l_x) + str(l_y) +'quenchVvscL.jpg', bbox_inches='tight', pad_inches=0.1)

    data.append(Scatter(x=c_Lplot[cLf_index:],y=avgssV[cLf_index:], name = 'non-quneching'))
    if display:
        layout = Layout(
        title = 'Average number of vortices against noise parameter',
        xaxis = dict(title = 'Noise parameter'),
        yaxis = dict(title = 'Number of vortices'),
        width = 700,
        height = 300)

        iplot(dict(data=data,layout=layout))
        
    #return cLs

def tfromcL(cL, Gamma_Q, cL_f, cL_i):
    #convert noise parameter values to t
    return Gamma_Q*((cL - cL_i)/(cL_f - cL_i))

def plot_t_hat(cLc, cL1s, Gamma_Qs, cL_f, cL_i):
    t1 = np.zeros(len(Gamma_Qs))
    tc = np.zeros(len(Gamma_Qs))
    t_hat = np.zeros(len(Gamma_Qs))
    cLcs = cLc*np.ones(len(Gamma_Qs))
    for i in range(0, len(Gamma_Qs)):
        t1[i] = tfromcL(cL1s[i], Gamma_Qs[i], cL_f, cL_i)
        tc[i] = tfromcL(cLcs[i], Gamma_Qs[i], cL_f, cL_i)
        t_hat[i] = tc[i] - t1[i]

    print 't1: ', t1, 'tc: ', tc, 't_hat: ', t_hat
    
    plt.figure()
    plt.plot(Gamma_Qs, t_hat, "o")
    plt.xlabel(r'$\hat{t}$', fontsize=16)
    plt.ylabel(r'$\Gamma_Q$', fontsize=16)
    #plt.title(r'Quench rate $\Gamma_Q$ against freezeout time $\hat{t}$')
    plt.show()
    
    return t_hat


def tau(tau_guess, c_L, N, l_x, l_y):
    yV = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" 
                    + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.txt')
    avgV = np.mean(yV, axis=0)
    xV = np.arange(0,len(avgV))
    
    f = np.zeros(len(xV) - tau_guess)
    for i in range(0,len(f)):
        f[i] = (avgV[tau_guess+i] - avgV[-1])/avgV[-1]

    data = [Scatter(x= range(tau_guess,tau_guess+len(f)), y=f)]
    layout = Layout(
    xaxis = dict(title = 'time'),
    yaxis = dict(title = 'f', range=[-0.001, 0.001]),
    width = 500)

    iplot(dict(data=data,layout=layout))
    
def uncer_that(cL1all, cLc, cLcstd, Gamma_Qs, cL_f, cL_i):
    cL1 = np.mean(cL1all, axis=0)
    cL1std = np.std(cL1all, axis=0)
    tc = np.zeros(len(Gamma_Qs))
    unt1 = np.zeros(len(Gamma_Qs))
    unthat = np.zeros(len(Gamma_Qs))
    untc = np.zeros(len(Gamma_Qs))
    for i in range (0,len(Gamma_Qs)):
        unt1[i] = Gamma_Qs[i]*cL1std[i]*0.5
        tc[i] = tfromcL(cLc, Gamma_Qs[i], cL_f, cL_i)
        untc[i] = Gamma_Qs[i]*cLcstd*0.5
        unthat[i] = np.sqrt(unt1[i]**2 + untc[i]**2)
    return cL1, unthat

from scipy.optimize import curve_fit
def func(x, a, b, c):
    return a * np.exp(-b*x) + c

def tau(tau_guess, c_L, N, l_x, l_y, fwindow, display = True, save = True):
    yV = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" 
                    + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.txt')
    avgV = np.mean(yV, axis=0)
    xV = np.arange(0,len(avgV))
    
    f = np.zeros(len(xV) - tau_guess)
    for i in range(0,len(f)):
        f[i] = (avgV[tau_guess+i] - avgV[-1])/avgV[-1]
        
    t = range(tau_guess,tau_guess+len(f))
        
    ythat = savgol_filter(f, fwindow, 2)

    
    if save: 
        plt.figure(figsize=(550/my_dpi, 330/my_dpi), dpi=my_dpi)
        plt.plot (x,f, label = 'original data')
        plt.plot (x,ythat, label = 'fitted line')
        plt.xlabel('t', fontsize =16)
        plt.ylabel('f', fontsize =16)
        plt.legend()
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" 
                    + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'ft.jpg', bbox_inches='tight', pad_inches=0.1)
    if display:
        data = [Scatter(x= t, y=f),
               Scatter(x= t, y=ythat)]
        layout = Layout(
        xaxis = dict(title = 'time'),
        yaxis = dict(title = 'f'), #range=[-0.001, 0.001]
        width = 500)

        iplot(dict(data=data,layout=layout))
    
    return t, f

def plottaucLepsi(cL, cLc, Tau, unTau, N, l_x, l_y, ls = False):
    plt.figure(figsize=(550/my_dpi, 330/my_dpi), dpi=my_dpi)
    plt.errorbar(cL, Tau, yerr = unTau, linestyle = "none", capsize=5, elinewidth=1, color = "C0") 
    plt.plot(cL, Tau, label = r'original data', color = "C0")
    if ls:
        popt, pcov = curve_fit(func, cL, Tau)
        plt.plot(cL, func(cL, *popt), label='least-squares', color = "C1")
    plt.xlabel(r'$c_L$')
    plt.ylabel(r'$\tau$')
    plt.legend()
    plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + str(l_x) + str(l_y) +'TaucL.jpg', bbox_inches='tight', pad_inches=0.1)
    
    epsilon = (cL - cLc*np.ones(len(Tau)))/cLc
    plt.figure(figsize=(550/my_dpi, 330/my_dpi), dpi=my_dpi)
    plt.errorbar(epsilon, Tau, yerr = unTau, linestyle = "none", capsize=5, elinewidth=1, color = "C0") 
    plt.plot(epsilon, Tau, label = r'original data', color = "C0")
    if ls:
        popt, pcov = curve_fit(func, epsilon, Tau)
        plt.plot(epsilon, func(epsilon, *popt), label='least-squares', color = "C1")
    plt.xlabel(r'$\epsilon$')
    plt.ylabel(r'$\tau$')
    plt.legend()
    plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + str(l_x) + str(l_y) +'TauEpsilon.jpg', bbox_inches='tight', pad_inches=0.1)

    return epsilon

def plotTaut(c_c, Tau, TauQA, t_i, t_f, epsilon, unTau, save = False, display = True):
    t_c = tfromcL(c_c, TauQA, 2, 7)
    t_const = np.arange(t_i, t_f)
    TauQ = TauQA
    t = epsilon*TauQ
    #TauQ2 = 2*TauQA*c_c/5
    #t2 = 2*eps01*TauQ

    if save:
        plt.figure(figsize=(550/my_dpi, 330/my_dpi), dpi=my_dpi)
        plt.errorbar(t, Tau, yerr = unTau, linestyle = "none", capsize=5, elinewidth=1, color = "C0") 
        plt.plot(t, Tau, label = r'$\tau(t)$ for $\tau_Q = %s$'%(TauQA), color = "C0")
        #plt.errorbar(t2, Tau, yerr = unTau, linestyle = "none", capsize=5, elinewidth=1, color = "C2") 
        #plt.plot(t2, Tau, label = r'$\tau(t)$ for $\tau_Q = %s$'%(TauQA*2), color = "C2")
        plt.plot(t_const,t_const, label = r'$t(t)$', color = "C1")
        plt.xlabel(r'$t$', fontsize =16)
        plt.ylabel(r'$\tau$', fontsize =16)
        plt.legend()
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + str(l_x) + str(l_y) +'Taut.jpg', bbox_inches='tight', pad_inches=0.1)
        
    if display:
        data = [Scatter(x= t, y=Tau, mode = "lines",error_y=dict(type='data',array=unTau, visible=True)),
               Scatter(x= t_const, y=t_const, mode="lines")] 
        layout = Layout(
        xaxis = dict(title = 't'),
        yaxis = dict(title = 'Tau'), #range=[-0.001, 0.001]
        height = 400,
        width = 500)

        iplot(dict(data=data,layout=layout))
    

        
"""    if i <= 9:
        figE = plt.figure()
        figE.add_subplot(5,2,i+1)
        plotE(xE, avgE, c_Lplot[i], N, l_x, l_y)
        plt.savefig('KPZProjectGraphs/Plots'+ str(N) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E1.jpg', bbox_inches='tight', pad_inches=0.1)

        figV = plt.figure()
        figV.add_subplot(5,2,i+1)
        plotE(xV, avgV, c_Lplot[i], N, l_x, l_y)
        plt.savefig('KPZProjectGraphs/Plots'+ str(N) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V1.jpg', bbox_inches='tight', pad_inches=0.1)

    elif i > 9 and i <= 19:
        figE = plt.figure()
        figE.add_subplot(5,2,i+1-10)
        plotE(xE, avgE, c_Lplot[i], N, l_x, l_y)
        plt.savefig('KPZProjectGraphs/Plots'+ str(N) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E2.jpg', bbox_inches='tight', pad_inches=0.1)
        
        figV = plt.figure()
        figV.add_subplot(5,2,i+1-10)
        plotV(xV, avgV, c_Lplot[i], N, l_x, l_y)
        plt.savefig('KPZProjectGraphs/Plots'+ str(N) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V2.jpg', bbox_inches='tight', pad_inches=0.1)

    elif i > 19:
        figE = plt.figure()
        figE.add_subplot(5,2,i+1-20)
        plotE(xE, avgE, c_Lplot[i], N, l_x, l_y)
        plt.savefig('KPZProjectGraphs/Plots'+ str(N) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E3.jpg', bbox_inches='tight', pad_inches=0.1)

        figV = plt.figure()
        figV.add_subplot(5,2,i+1-20)
        plotV(xV, avgV, c_Lplot[i], N, l_x, l_y)
        plt.savefig('KPZProjectGraphs/Plots'+ str(N) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V3.jpg', bbox_inches='tight', pad_inches=0.1)
"""

"""CV = np.zeros(len(c_Lnew))
CV2 = np.zeros(len(c_Lplot)-1)
for i in range(0,len(c_Lplot)):
    if i+1 < len(c_Lplot):
        C = np.var((avgssE64_0_0[i], avgssE64_0_0[i+1]))/(c_Lnew[i]**2)
        print C
        CV[i] = C
        CV2[i] = np.var(avgssE64_0_0)/(c_Lplot[i+1]**2)

line = Scatter(x=c_Lplot[1:],y=CVs64_0_0[1:], mode = "lines+markers")
data = [line]

layout = Layout(
    title = 'Specific heat using equation',
    xaxis = dict(title = 'Noise'),
    yaxis = dict(title = 'C_V'),
    width = 700)

iplot(dict(data=data,layout=layout))"""
    