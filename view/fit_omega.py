from array import array
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit

pi = np.pi
clight = 299792458

def get_the_turn_number(c,value):
    print(len(c))
    for i in range(len(c)):
        if c[i]>=value:
            print(i,",",c[i])
            return i
    return i

def func_exp_fit(x, a, b, c):
    return a*np.exp(b*x)+c

def fit_growth(N_mode,nTurns,thresh,nBunch,Gamma0,eta,Qs,T0,nDynamic,M1_1,M1_2,temppath):
    ids = np.array([i for i in range(nBunch)])
    qs = np.ndarray([N_mode])

    # calculate the components
    sincomp = []
    coscomp = []
    a1 = np.ndarray([N_mode,nTurns])
    b1 = np.ndarray([N_mode,nTurns])
    a2 = np.ndarray([N_mode,nTurns])
    b2 = np.ndarray([N_mode,nTurns])

    c = np.ndarray([N_mode,nTurns])
    d = np.ndarray([N_mode,nTurns])
    e = np.ndarray([N_mode,nTurns])

    value = thresh # the value where we stop the plotting

    idx = np.zeros(N_mode)
    qidx = 0

    fig1,axes1 = plt.subplots(N_mode+1,1)
    fig1.set_figheight(10*(N_mode))
    fig1.set_figwidth(30)

    for mu in range(N_mode):
        #sincomp.append(np.sin((nBunch-mu)*2*pi*f0*np.array(centroids)))
        #coscomp.append(np.cos((nBunch-mu)*2*pi*f0*np.array(centroids)))
        sincomp.append(np.sin((nBunch+mu)*2*pi/nBunch*ids))
        coscomp.append(np.cos((nBunch+mu)*2*pi/nBunch*ids))
        for i in range(nTurns):
            a1[mu][i] = np.sum(sincomp[mu]*M1_2[i*nBunch:(i+1)*nBunch])/len(sincomp[mu])*2
            b1[mu][i] = np.sum(coscomp[mu]*M1_2[i*nBunch:(i+1)*nBunch])/len(coscomp[mu])*2
            a2[mu][i] = np.sum(sincomp[mu]*M1_1[i*nBunch:(i+1)*nBunch])/len(sincomp[mu])*2
            b2[mu][i] = np.sum(coscomp[mu]*M1_1[i*nBunch:(i+1)*nBunch])/len(coscomp[mu])*2
        c[mu] = np.sqrt(np.array(a1[mu])**2+np.array(b1[mu])**2) #gamma
        d[mu] = np.sqrt(np.array(a2[mu])**2+np.array(b2[mu])**2)*Gamma0/eta*Qs*2*pi/T0 #time
        e[mu] = np.sqrt(c[mu]**2+d[mu]**2)

        idx[mu] = get_the_turn_number(c[mu],value)
        print("The turn number when amplitude reach ",value,":",idx[mu])

        # fit the growth rate
        rng1 = nDynamic
        rng2 = int(idx[mu])
        if rng1>rng2:
            rng1 = 0
        if c[mu][rng2]<c[mu][0]:
            rng2 = nTurns
        #rng2 = rng1+7000
        cNew = c[mu][rng1:rng2]
        turn_temp = range(rng2-rng1)
        # the guessed value of growth rate
        q_guess = 1
        popt, pcov = curve_fit(func_exp_fit, turn_temp, \
                               cNew,bounds=([-cNew[0]*20-1,-q_guess-1,-cNew[0]*0.5-1], [cNew[0]*20+1, q_guess+1, cNew[0]*0.5+1]),\
                               maxfev=200000000)
        print(popt)
        print(pcov)
        qs[mu] = popt[1]
        c_fit =  popt[0]*np.exp(popt[1]*turn_temp)+popt[2]
        axes1[mu].plot(c[mu][rng1:rng2],'gx')
        axes1[mu].plot(c_fit,'r-')
        axes1[mu].tick_params(labelsize=50)

        axes1[mu].set_xlabel('# of Turn',fontsize=30)
        axes1[mu].set_ylabel('Amp mu = '+str(-mu),fontsize=30)

        axes1[mu].tick_params(labelsize=30)
        y_labels = axes1[mu].get_yticks()
        axes1[mu].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.6f'))
    qidx+=1 # for each parameter sample point
    fn_after = os.path.join(temppath,'CB_Modes.png')
    plt.savefig(fn_after,bbox_inches='tight')
    plt.show()
    return qs