from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import plotly.plotly as py
from plotly import tools
from plotly.graph_objs import *
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
init_notebook_mode()


def EV(c_L, N, l_x, l_y, plot = False, display = False, save = False):
    """Function to return average energy and vortices value at steady state given c_L, lattice size, lambda_x, lambda_y
    with option to plot average energy and vortices against time, to display and/or to save plots"""
    
    #import energy data
    yE = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                    str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E.txt')
    avgE = np.mean(yE, axis=0) #average energy data
    xE = np.arange(0,len(avgE))
    
    #import vortices data
    yV = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" 
                    + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.txt')
    avgV = np.mean(yV, axis=0)
    xV = np.arange(0,len(avgV))
    
    if plot:
        #plot energy data
        plt.figure()
        plt.plot(xE,avgE)
        #plt.title(r'Average energy density over time for $c_L =  %s$'%(c_L) + '\n' + r'and lattice size $%s^2$' %(N)
                  #+ r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel('time [T]')
        plt.ylabel(r'$<E>/N^2$')
        
        if save:
            plt.savefig('KPZProjectGraphs/Plots'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                        str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E.jpg', bbox_inches='tight', pad_inches=0.1)
        
        if display:
            plt.show()
        

        #plot vortices data
        plt.figure()
        plt.plot(xV,avgV)
        #plt.title(r'Average number of vortices over time for $c_L =  %s$'%(c_L) + '\n' + r'and lattice size $%s^2$' %(N)
                  #+ r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel('time [T]')
        plt.ylabel('number of vortices')
        
        if save:
            plt.savefig('KPZProjectGraphs/Plots'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                        str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.jpg', bbox_inches='tight', pad_inches=0.1)
        
        if display:
            plt.show()
        

    print "cL: ", c_L , " E_R, E_iterations: ", yE.shape , " V_R, V_iterations: ", yV.shape
    
    return avgE[-1], avgV[-1]

def dataEV(c_L, N, l_x, l_y, itersRinfo = True):
    #import energy data
    yE = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                    str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E.txt')
    avgE = np.mean(yE, axis=0) #average energy data
    varssE = np.var(yE[:,-1], ddof = 0)/(c_L**2) #variance of last column gives variance of steady value
    xE = np.arange(0,len(avgE))
    
    #import vortices data
    yV = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" 
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
    plt.xlabel('time [T]')
    plt.ylabel(r'$<E>/N^2$')
    
    if save:
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                        str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E.jpg', bbox_inches='tight', pad_inches=0.1)
    
def plotV(xV, avgV, c_L, N, l_x, l_y, save = False):
    plt.figure()
    plt.plot(xV,avgV)
    #plt.title(r'Average number of vortices over time for $c_L =  %s$'%(c_L) + '\n' + r'and lattice size $%s^2$' %(N)
              #+ r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
    plt.xlabel('time [T]')
    plt.ylabel('number of vortices')
    
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
        ax1.set_xlabel('time [T]', fontsize='small')
        ax1.set_ylabel(r'$<E>/N^2$', fontsize='small')

        #ax2.set_title(r'Average number of vortices over time', fontsize='medium')
        ax2.set_xlabel('time [T]', fontsize='small')
        ax2.set_ylabel('r$<V>$', fontsize='small')
        
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

def avgEVt(c_L, N, l_x, l_y, display = False):
    """Function to return average steady state energy and vortices values and average number of iterations needed for
    steady state energy and vortices values given c_L, lattice size, lambda_x, lambda_y with option to display the whole
    data with energy/vortices values in first column and final iterations in second column"""
    #import energy iterations data
    tE = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                    str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'E_t.txt')
    avgE = np.mean(tE[:,0]) #average of first column gives average steady state energy
    avgEt = np.mean(tE[:,1]) #average of second column gives average final iterations
    stdEt = np.std(tE[:,1]) #variance of final iterations
    tE_indices = np.argsort(tE[:,1]) #indices of data in ascending order
    if display:
        print "cL: ", c_L , " E values, final iterations: ", tE[tE_indices]
        
    #import vortices iterations data
    tV = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(c_L, '.2f')) + "-" + 
                    str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V_t.txt')
    avgV = np.mean(tV[:,0]) #average of first column gives average steady state energy
    avgVt = np.mean(tV[:,1]) #average of second column gives average final iterations
    stdVt = np.std(tV[:,1]) #variance of final iterations
    tV_indices = np.argsort(tV[:,1]) #indices of data in ascending order
    if display:
        print "cL: ", c_L , " V values, final iterations: ", tV[tV_indices]
    
    return format(c_L, '.2f'), format(avgE, '.2f'), format(avgEt, '.0f'), format(stdEt, '.0f'), format(avgV, '.2f'), format(avgVt, '.0f'), format(stdVt, '.0f')

def plotavgE(c_Lplot, avgssE, N, l_x, l_y, display = True, save = False):
    """Plot average steady state energy given c_L and energy values, N, l_x, l_y with option to display and/or to show plots"""
    if save:
        plt.figure()
        plt.plot(c_Lplot,avgssE) #"o" for scatter
        #plt.title(r'Average energy density against noise parameter for lattice size $%s^2$' %(N) + '\n' +
                  #r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel(r'$c_L$')
        plt.ylabel(r'$<E>/N^2$')
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
        plt.figure()
        plt.plot(c_Lplot,avgssV)
        plt.title(r'Average number of vortices against noise parameter for lattice size $%s^2$' %(N)+ '\n' +
                 r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel(r'$c_L$')
        plt.ylabel(r'$number of vortices$')
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-"  + str(format(l_x, '.2f')) + "-" +
                                str(format(l_y, '.2f')) + 'avgV.jpg', bbox_inches='tight', pad_inches=0.1)
        
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
        plt.figure()
        plt.plot (c_Lnew,dEdcL) #"o" for scatter
        #plt.title('Change in average energy density over change in noise \n' + r'against noise for lattice size $%s^2$' %(N)+
                  #r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel(r'$c_L$')
        plt.ylabel(r'$d<E>/dc_L$')
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

def plot_specheateqn(c_Lplot, CVs, N, l_x, l_y, display = True, save = False):
    """Function to plot C_V from using equation with option to plot and/or save graphs"""
    
    if save:
        plt.figure()
        plt.plot (c_Lplot[8:23],CVs[8:23]) #"o" for scatter
        #plt.title(r'Specific heat using equation for lattice size $%s^2$'%(N) + '\n' +
                  #r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        plt.xlabel(r'$c_L$')
        plt.ylabel(r'$C_V$')
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-"  + str(format(l_x, '.2f')) + "-" +
                    str(format(l_y, '.2f')) + 'specheateqn.jpg', bbox_inches='tight', pad_inches=0.1)

    if display:
        line = Scatter(x=c_Lplot[8:23],y=CVs[8:23], mode = "lines+markers")
        data = [line]

        layout = Layout(
            title = 'Specific heat using equation',
            xaxis = dict(title = 'Noise'),
            yaxis = dict(title = 'C_V'),
            width = 700)

        iplot(dict(data=data,layout=layout))


def plotquenchV(xV, avgV, Gamma_Q, cL_i, cL_f, N, l_x, l_y, save = False):
    plt.figure()
    plt.plot(xV,avgV)
    #plt.title(r'Average number of vortices over time for $\Gamma_Q =  %s$'%(Gamma_Q) + r', $c_{L_i} = %s$'%(cL_i) +  r' and $c_{L_f} = %s$'%(cL_f) + '\n' + r'with lattice size $%s^2$' %(N) + r', $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
    plt.xlabel('time [T]')
    plt.ylabel('number of vortices')
    
    if save:
        plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + "-" + str(format(cL_i, '.2f')) + "-" + str(format(cL_f, '.2f')) + "-"
        + str(Gamma_Q) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.jpg', bbox_inches='tight', pad_inches=0.1)

def plotVquench(cL_i, cL_f, Gamma_Qs, N, l_x, l_y, display = False, save = False):
    """Function to plot number of vortices against time for different quench rates with option to display and/or save graphs
    and returns the steady state number of vortices for each quench rate"""
    data = []
    logVss = np.zeros(len(Gamma_Qs))
    logtGamma = np.zeros(len(Gamma_Qs))
    

    for i in range(0,len(Gamma_Qs)): 
        yV = np.loadtxt('KPZProjectGraphs/'+ str(N) + "-" + str(format(cL_i, '.2f')) + "-" + str(format(cL_f, '.2f')) + "-"
        + str(Gamma_Qs[i]) + "-" + str(format(l_x, '.2f')) + "-" + str(format(l_y, '.2f')) + 'V.txt')
        avgV = np.mean(yV, axis=0) #average number of vortices data
        xV = np.arange(0,len(avgV))
        data.append(Scatter(x=xV,y=avgV, name = 'Gamma_Q = %s'%(Gamma_Qs[i])))
        logtGamma[i] = np.log2(len(avgV)/avgV[-1])
        logVss[i] = np.log2(avgV[-1])
        
        
        if save:
            plotquenchV(xV, avgV, Gamma_Qs[i], cL_i, cL_f, N, l_x, l_y, save = True)
        
    if display:
        layout = Layout(
        title = 'Average number of vortices against time',
        xaxis = dict(title = 'Time'),
        yaxis = dict(title = 'Number of vortices'),
        width = 700)

        iplot(dict(data=data,layout=layout))
        
    return logVss, logtGamma

def plotVquenchvscL(cL_i, cL_f, c_Lplot, avgssV, Gamma_Qs, N, l_x, l_y, display = False, save = False):
    """Function to plot number of vortices against cL for different quench rates with option to display and/or save graphs"""
    data = []
    cLs = []
    
       
    if save:
        fig, ax = plt.subplots()
        ax.set_title(r'Average number of vortices against noise parameter for lattice size $%s^2$' %(N)+ '\n' +
                     r' with $\lambda_x = %s$ '%(l_x) + r'and $\lambda_y = %s$'%(l_y))
        ax.set_xlabel(r'$c_L$')
        ax.set_ylabel(r'average number of vortices')
        color=plt.cm.rainbow(np.linspace(0,1,len(Gamma_Qs)))
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
        
        data.append(Scatter(x=cL,y=avgV, name = 'Gamma_Q = %s'%(Gamma_Qs[i]), mode = "lines+markers"))
        
        cLf_index = int(np.where(c_Lplot==cL_f)[0])
        
        if save:
            ax.plot(cL, avgV, label = r'$\Gamma_Q$ = %s'%(Gamma_Qs[i]), c = color[i])
            ax.plot(c_Lplot[cLf_index:],avgssV[cLf_index:], label = r'steady state')
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 'small')
            plt.savefig('KPZProjectGraphs/Plots/'+ str(N) + 'quenchVvscL.jpg', bbox_inches='tight', pad_inches=0.1)
            
    data.append(Scatter(x=c_Lplot[cLf_index:],y=avgssV[cLf_index:], name = 'steady state'))
    if display:
        layout = Layout(
        title = 'Average number of vortices against noise parameter',
        xaxis = dict(title = 'Noise parameter'),
        yaxis = dict(title = 'Number of vortices'),
        width = 700)

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

    print 't1: ', t1, 'tc: ', tc[0], 't_hat: ', t_hat
    
    plt.figure()
    plt.plot(Gamma_Qs, t_hat, "o")
    plt.xlabel(r'$\hat{t}$')
    plt.ylabel(r'$\Gamma_Q$')
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
    