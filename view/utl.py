from array import array
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
pi = np.pi
clight = 299792458
# Get the input parameters
def get_input_para(inputfile):
    tempinput = {}
    with open(inputfile) as inputfile:
        for line in inputfile:
            if len(line.split())>1:
                tempinput[line.split()[0]] = line.split()[1:]
        for i in tempinput:
            for j in range(len(tempinput[i])):
                tempinput[i][j] = float(tempinput[i][j])
    return tempinput
def read_M1(cwd,nTurns,nBeam,nBunch,Gamma0):
    
    M1_all = array('d')
    M1_fn = 'M1_all.bin'
    datafile = os.path.join(cwd,M1_fn)    
    with open(datafile, mode='rb') as file: # b is important -> binary
        M1_all.fromfile(file,2*nTurns*nBeam*nBunch)
    M1_1 = []
    for i in range(nTurns):
        for j in range(nBunch):
            M1_1.append(M1_all[i*nBunch*2+j])
    M1_2 = []
    for i in range(nTurns):
        for j in range(nBunch):
            M1_2.append(M1_all[i*nBunch*2+nBunch+j]-Gamma0)
    return M1_all,M1_1,M1_2

def get_first_bunch_M1(nTurns,nBunch,Gamma0,M1_all):
    M1_1_0 = []
    M1_2_0 = []
    for i in range(nTurns):
        M1_1_0.append(M1_all[i*nBunch*2])
        M1_2_0.append(M1_all[i*nBunch*2+nBunch]-Gamma0)
    return M1_1_0,M1_2_0

def plot_firt_bunch_centroids(M1_1_0,M1_2_0,T0,h,temppath,nTurns,turn,n_turn_disp):
    fig1,axes1 = plt.subplots(1,1)
    fig1.set_figheight(10)
    fig1.set_figwidth(30)
    
    rng1 = turn
    rng2 = rng1+n_turn_disp
    axes1.plot(np.array(M1_1_0[rng1:rng2])/T0*h[0]*360,np.array(M1_2_0[rng1:rng2]),'rx-')
    axes1.plot(np.array(M1_1_0[rng1:rng1+10])/T0*h[0]*360,np.array(M1_2_0[rng1:rng1+10]),'gx-')
    axes1.plot(np.array(M1_1_0[rng2-10:rng2])/T0*h[0]*360,np.array(M1_2_0[rng2-10:rng2]),'bx-')
    #axes1.plot(np.array(M1_2_0[rng1:rng2]),'rx-')
    #axes1.axvline(x = 0)
    #axes1.axvline(x = 360)
    #axes1.axvline(x = 180)
    axes1.set_ylabel('first_centroid(Gamma)',fontsize=30)
    fn_after = os.path.join(temppath,'Centroids_first.png')
    plt.savefig(fn_after,bbox_inches='tight')
    plt.show()
    return
def plot_centroids_gamma(M1_1,M1_2,T0,h,temppath,nTurns,turn,n_turn_disp,nBunch):
    
    fig1,axes1 = plt.subplots(1,1)
    fig1.set_figheight(10)
    fig1.set_figwidth(30)
    
    rng1 = turn*nBunch
    rng2 = rng1+nBunch*n_turn_disp
    
    axes1.plot((np.array(M1_2[rng1:rng2])),'r.')
    
    axes1.set_ylabel('centroid (d_gamma)',fontsize=30)
    axes1.set_xlabel('Bunch # ',fontsize=30)
    axes1.tick_params(labelsize=30)
    fn_after = os.path.join(temppath,'Centroids_gamma_'+str(turn)+'_'+str((rng2-rng1)/nBunch)+'.pdf')
    plt.savefig(fn_after,bbox_inches='tight')
    plt.show()
    return

def get_centroids(Pattern,Trf,nBunch,t0,fill_step):
    centroids = np.ndarray(nBunch)
    bunch_idx = 0
    bucket_idx = 0
    nTrain = int(len(Pattern)/3)
    print("nTrain: ",nTrain)
    for i in range(nTrain):
        for j in range(Pattern[i*3]):
            #print("Bucket_ID: ",bucket_idx/424)
            centroids[bunch_idx]=Trf*bucket_idx+t0
            bunch_idx+=1
            bucket_idx+=Pattern[i*3+1]
            #print(bunch_idx,bucket_idx)
        bucket_idx+=Pattern[i*3+2]
    return centroids

def get_plot_idx(Pattern,turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain):
    bunch_idx = 0
    for i in range(train_start):
        bunch_idx += Pattern[i*3]
    bunch_idx+=bunch_start
    start = turn_start*nBunch+bunch_idx
    for i in range(train_end-train_start):
        bunch_idx += Pattern[(train_start+i)*3]
    bunch_idx+=bunch_end-bunch_start
    end = (turn_end)*nBunch+bunch_idx
    return start,end
    
def plot_centroids_t(M1_1,centroids,T0,h,temppath,nTurns,nTrain, turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,Pattern, nBunch):
    
    fig1,axes1 = plt.subplots(1,1)
    fig1.set_figheight(10)
    fig1.set_figwidth(30)
    
    
    rng1,rng2 = get_plot_idx(Pattern,turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
    rng3,rng4 = get_plot_idx(Pattern,0,0,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
    axes1.plot((np.array(M1_1[rng1:rng2])-centroids[rng3:rng4]),'r.')
    #axes1.plot(np.array(centroids),'r.')
    
    #axes1.axvline(x = nfill*nBunch)
    #axes1.axvline(x = (nfill+n_q_ramp)*nBunch)
    
    axes1.set_ylabel('centroid (d_t)',fontsize=30)
    axes1.set_xlabel('Bunch # ',fontsize=30)
    axes1.tick_params(labelsize=30)
    fn_after = os.path.join(temppath,'Centroids_t_'+str(turn_start)+'_'+str((rng2-rng1)/nBunch)+'.pdf')
    plt.savefig(fn_after,bbox_inches='tight')
    plt.show()
    return
def plot_centroids_phi(M1_1,centroids,T0,h,f,temppath,nTurns,nTrain, turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,Pattern, nBunch):
    
    fig1,axes1 = plt.subplots(1,1)
    fig1.set_figheight(10)
    fig1.set_figwidth(30)
    
    
    rng1,rng2 = get_plot_idx(Pattern,turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
    rng3,rng4 = get_plot_idx(Pattern,0,0,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
    axes1.plot((np.array(M1_1[rng1:rng2])-centroids[rng3:rng4])*f[0]*360,'r.')
    #axes1.plot(np.array(centroids),'r.')
    
    #axes1.axvline(x = nfill*nBunch)
    #axes1.axvline(x = (nfill+n_q_ramp)*nBunch)
    
    axes1.set_ylabel('centroid (d_phi [degree])',fontsize=30)
    axes1.set_xlabel('Bunch # ',fontsize=30)
    axes1.set_title("Phase slip along the train",fontsize=30)
    axes1.tick_params(labelsize=30)
    fn_after = os.path.join(temppath,'Centroids_phi_'+str(turn_start)+'_'+str((rng2-rng1)/nBunch)+'.png')
    return

def plot_dPhi(M1_1,centroids,T0,h,f,Trf, temppath,nTurns,nTrain, turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,Pattern, nBunch):
    fig1,axes1 = plt.subplots(1,1)
    fig1.set_figheight(10)
    fig1.set_figwidth(30)
    
    rng1,rng2 = get_plot_idx(Pattern,turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
    rng3,rng4 = get_plot_idx(Pattern,0,0,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
    
    Phases = (np.array(M1_1[rng1:rng2])-centroids[rng3:rng4])*f[0]*360
    Phases1 = Phases[:int(nBunch/2)]
    Phases2 = Phases[int(nBunch/2):]
    dPhases = Phases1-Phases2
    dTime = dPhases/2/pi*Trf[0]

    axes1.plot(dTime,'r.')

    axes1.set_ylabel('Time difference [degree]',fontsize=30)
    axes1.set_xlabel('Bunch # ',fontsize=30)
    axes1.set_title("Time difference between colliding bunches",fontsize=30)
    axes1.tick_params(labelsize=30)
    fn_after = os.path.join(temppath,'Time_Difference'+str(turn_start)+'.png')
    plt.savefig(fn_after,bbox_inches='tight')
    plt.show()
    
    fig1,axes1 = plt.subplots(1,1)
    fig1.set_figheight(10)
    fig1.set_figwidth(30)
    axes1.plot(dPhases,'r.')
    axes1.set_ylabel('Phase difference [degree]',fontsize=30)
    axes1.set_xlabel('Bunch # ',fontsize=30)
    axes1.set_title("Phase difference between colliding bunches",fontsize=30)
    axes1.tick_params(labelsize=30)
    fn_after = os.path.join(temppath,'Phase_Difference_'+str(turn_start)+'.png')
    plt.savefig(fn_after,bbox_inches='tight')
    plt.show()
    
    omegarf = 2*pi/Trf
    Lambda = clight/(omegarf/2/pi)
    fig1,axes1 = plt.subplots(1,1)
    fig1.set_figheight(10)
    fig1.set_figwidth(30)
    axes1.plot(dPhases/360*Lambda[0]*1000,'r.')
    dz = dPhases[7:]/360*Lambda[0]*1000
    dz_abs_mean = np.mean(np.abs(dz))
    axes1.set_ylabel('Position difference [mm]',fontsize=30)
    axes1.set_xlabel('Bunch # ',fontsize=30)
    axes1.set_title("Position difference between colliding bunches",fontsize=30)
    axes1.tick_params(labelsize=30)
    #axes1.text(100,dz_abs_mean,r'$<\|dz\|>=$'+'{:.2e}'.format(dz_abs_mean)+'[mm]', fontsize=15)
    fn_after = os.path.join(temppath,'Position_Difference_'+str(turn_start)+'.png')
    plt.savefig(fn_after,bbox_inches='tight')
    plt.show()
    return
def plot_anim(n_start,n_end,n_step,omegarf,for_gif1,for_gif2,cal_dphase,nBunch,nTrain,Pattern,M1_1,centroids,f,temppath):
    
    Lambda = clight/(omegarf/2/pi)
    dz_abs_mean_record = []
    dz_abs_record = []
    
    phi_v_train = []
    for n in range(n_start,n_end,n_step):

        turn_start = n
        turn_end = n
        train_start = 0
        train_end = 0
        bunch_start = 0
        bunch_end = nBunch

        rng1,rng2 = get_plot_idx(Pattern,turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
        rng3,rng4 = get_plot_idx(Pattern,0,0,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
        phi_v_train.append((np.array(M1_1[rng1:rng2])-centroids[rng3:rng4])*f[0]*360)

    ymin = np.min(phi_v_train)
    ymax = np.max(phi_v_train)
    if for_gif1:
        for i in range(len(phi_v_train)):
            fig1,axes1 = plt.subplots(1,1)
            fig1.set_figheight(23.4/2)
            fig1.set_figwidth(63.1/2)
            axes1.plot(phi_v_train[i],'r.')
            axes1.set_ylabel('centroid (d_phi [degree])',fontsize=30)
            axes1.set_xlabel('Bunch # ',fontsize=30)
            axes1.set_title("Phase slip along the train",fontsize=30)
            #axes1.set_ylim([ymin,ymax])
            axes1.tick_params(labelsize=30)
            axes1.text(-20,np.max(phi_v_train[i]),"Turn number: "+str(i*n_step+n_start), fontsize=15)
            fn_after = os.path.join(temppath,'Centroids_phi_'+'{:0>4d}'.format(n_start)+'_'+'{:0>4d}'.format(i*n_step+n_start)+'.png')
            plt.savefig(fn_after,bbox_inches='tight')
            plt.show()

    if cal_dphase:
        for i in range(len(phi_v_train)):
            Phases1 = phi_v_train[i][:int(nBunch/2)]
            Phases2 = phi_v_train[i][int(nBunch/2):]
            dPhases = Phases1-Phases2
            dz = dPhases[:]/360*Lambda[0]*1000
            #print(len(dz_abs_record))
            #print(len(dz))
            dz_abs_mean = np.mean(np.abs(dz))
            dz_abs_mean_record.append(dz_abs_mean)
            dz_abs_record.append(dz)
        fig1,axes1 = plt.subplots(1,1)
        fig1.set_figheight(10)
        fig1.set_figwidth(15)
        axes1.plot(range(n_start,n_end,n_step),dz_abs_mean_record,'r.-')
        axes1.set_ylabel('k [mm]',fontsize=30)
        axes1.set_xlabel('Turn # ',fontsize=30)
        axes1.set_title("k v.s. turns",fontsize=30)
        axes1.text(800,dz_abs_mean_record[-1]*1.25,"Final mismatch: "+'{:.2e}'.format(dz_abs_mean_record[-1])+" [mm]", fontsize=15)
        axes1.tick_params(labelsize=30)
        fn_after = os.path.join(temppath,'Position_Difference_vs_Turn.png')
        plt.savefig(fn_after,bbox_inches='tight')
        plt.show()

# Get the VI data
def get_VI(datafile,n_stride,nTurns,store_step,NpRF,h,nRF,nBeam,omegarf,V0,V0Q,II,IQ):
    test = array('d')
    with open(datafile, mode='rb') as file: # b is important -> binary
        test.fromfile(file,int((n_stride)*nTurns/store_step*NpRF*h[0]))
    time = np.array(test[0::n_stride],dtype='float64')
    Ibi = np.array(test[1::n_stride],dtype='float32')
    #Ibi_1 = np.array(test[2::n_stride])
    Ibi2 = np.array(test[3::n_stride],dtype='float32')

    V = []
    Ig = []
    #U = []
    Vref = []
    Iref = []
    for i in range(nRF):
        V.append(np.array(test[1+2*nBeam+0+i*3::n_stride],dtype='float32'))
        Ig.append(np.array(test[1+2*nBeam+1+i*3::n_stride],dtype='float32'))

        Vref.append(V0[i]*np.sin(omegarf[i]*time)+V0Q[i]*np.cos(omegarf[i]*time))
        Iref.append(II[i]*np.sin(omegarf[i]*time)+IQ[i]*np.cos(omegarf[i]*time))

    return V,Ig,Ibi,Ibi2,Vref,Iref

# Plot VI
def plot_VI(startTurn,startRF,nRF,nRFsamp,NpRF,nBeam,h,V,Vref,Iref,Ig,Ibi,Ibi2,cwd):
    print(nRF)
    rng1 = NpRF*(h[0]*startTurn+startRF)
    rng2 = NpRF*(h[0]*startTurn+startRF+nRFsamp)
    rng3 = 50
    step1 = 1
    step2 = 1

    if nRF>1:
        fig1,axes1 = plt.subplots(4,nRF+1)
        for i in range(nRF):
            axes1[0][i].plot(V[i][rng1:rng2:step1],'r.-',ms=10)
            axes1[0][i].plot((Vref[i])[rng1:rng2:step1],'g.-',ms=10)
            axes1[0][i].axhline(y=0)
            #axes1[0][i].set_xlim([int(NpRF/2)-rng3,int(NpRF/2)+rng3])
            #axes1[0][i].set_ylim([-4e6,4e6])
            axes1[0][0].axhline(y=3.7e6/14*9)
            axes1[0][1].axhline(y=3.7e6/14*5)

            #axes1[0][1].axhline(y=3.7e6/14*6)
            #axes1[0][1].axhline(y=7.02e7)
            axes1[0][i].axvline(x=int(NpRF/2))
            axes1[1][i].axvline(x=int(NpRF/2))
            axes1[1][i].plot(Ibi[rng1:rng2:step1],'r.-',ms=10)
            if nBeam == 2 :
                axes1[1][i].plot(Ibi2[rng1:rng2:step1],'g.-',ms=10)
            #axes1[1][i].set_xlim([int(NpRF/2)-rng3,int(NpRF/2)+rng3])
            #axes1[2][i].set_xlim([int(NpRF/2)-rng3,int(NpRF/2)+rng3])
            axes1[2][i].axhline(y=0)

            axes1[2][i].plot((V[i]-Vref[i])[rng1:rng2:step1],'r.-',ms=10)
            axes1[2][i].axvline(x=int(NpRF/2))
            #axes1[3][i].plot((Ig[i])[rng1:rng2:step1],'r.-',ms=10)
            #axes1[3][i].plot((Iref[i])[rng1:rng2:step1],'g-',ms=10)
            axes1[3][i].plot((Iref[i]-Ig[i])[rng1:rng2:step1],'g-',ms=10)
            axes1[3][i].axvline(x=int(NpRF/2))
            #axes1[0][i].set_ylim([-3e6,3e6])
        Vsum = np.sum(V,0)
        Vref_sum = np.sum(Vref,0)

        axes1[0][nRF].plot(Vsum[rng1:rng2:step1],'r.-',ms=10)
        axes1[0][nRF].plot(Vref_sum[rng1:rng2:step1],'g.-',ms=10)
        #axes1[0][nRF].set_xlim([int(NpRF/2)-rng3,int(NpRF/2)+rng3])
        #axes1[0][nRF].set_ylim([-4e6,4e6])
        axes1[0][nRF].axhline(y=0)
        axes1[0][nRF].axhline(y=3.7e6)
        axes1[1][nRF].axvline(x=int(NpRF/2))
        axes1[0][nRF].axvline(x=int(NpRF/2))
        #axes1[1][nRF].set_xlim([int(NpRF/2)-rng3,int(NpRF/2)+rng3])
        axes1[1][nRF].plot(Ibi[rng1:rng2:step1],'r.-',ms=10)
        if nBeam == 2 :
            axes1[1][nRF].plot(Ibi2[rng1:rng2:step1],'g.-',ms=10)
        axes1[2][nRF].plot((Vsum-Vref_sum)[rng1:rng2:step1],'r.-',ms=10)
        axes1[2][nRF].axvline(x=int(NpRF/2))
        axes1[3][nRF].axvline(x=int(NpRF/2))

    else:
        fig1,axes1 = plt.subplots(4,1)
        for i in range(nRF):
            axes1[0].plot(V[i][rng1:rng2:step1],'rx-',ms=10)
            axes1[0].plot((Vref[i])[rng1:rng2:step1],'g.-',ms=10)
            axes1[0].axhline(y=0)
            
            axes1[0].axvline(x=int(NpRF/2))
            axes1[1].plot(Ibi[rng1:rng2:step1],'r.-',ms=10)
            axes1[1].axvline(x=int(NpRF/2))
            axes1[2].plot((V[i]-Vref[i])[rng1:rng2:step1],'r.-',ms=10)
            axes1[2].axhline(y=0)
            axes1[2].axvline(x = int(NpRF/2))
            axes1[3].plot((Ig[i])[rng1:rng2:step1],'r.-',ms=10)
            axes1[3].plot((Iref[i])[rng1:rng2:step1],'g.-',ms=10)
            axes1[3].axhline(y=0)
            axes1[3].axvline(x=int(NpRF/2))
    fig1.set_figheight(30)
    fig1.set_figwidth(30)
    #print((V[0]-Vref[0])[rng1:rng2:step1][int(NpRF/2)]/1e6)
    fn_VI = os.path.join(cwd,'VI2.pdf')
    plt.savefig(fn_VI,bbox_inches='tight')
    plt.show()