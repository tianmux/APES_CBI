import numpy as np
import copy
import sys
from array import array
import itertools
import os
import sys
from scipy import signal
from scipy.optimize import curve_fit
import subprocess
import shutil

pi = np.pi
clight = 299792458
E0Au = 196.9665687*931.5e6
E0Elec = 0.51099895000e6

# this is the function generates the beam pattern
def get_beam_pattern(nTrain, nBunch, gap, h, fill_step):
    pattern = np.ndarray(nTrain*2)
    nBunch_avai = h/fill_step # the total available postions for bunch
    nBunch_avai_pTrain = int(nBunch_avai/nTrain)

    nBunch_pTrain = int(nBunch_avai_pTrain*(1-gap)) # the number of bunches per train
    nGap_pTrain = nBunch_avai_pTrain-nBunch_pTrain # the number of bunches per train 
    for i in range(nTrain):
        pattern[i*2] = nBunch_pTrain
        pattern[i*2+1] = nGap_pTrain
    for i in range(int(nBunch_avai-nBunch_avai_pTrain*nTrain)):
        pattern[i*2]+=1
    temp1 = 0
    temp2 = 0
    
    for i in range(nTrain):
        temp1 += pattern[i*2]
        temp2 += pattern[i*2+1]
    print(temp1)
    print(temp2)
    
    return temp1, pattern

# The function to read in the default input parameters
def get_input(tempinput,input_fn):
    with open(input_fn) as inputfile:
        for line in inputfile:
            if len(line.split())>1:
                tempinput[line.split()[0]] = line.split()[1:]
    for i in tempinput:
        for j in range(len(tempinput[i])):
            tempinput[i][j] = float(tempinput[i][j])
    return tempinput

# The function to generate the whole parameter space
# The current upper limit of dimensions is 2, so we can only scan 2 different parameters. 
def get_scan_para_space(para_to_scan,input_original,input_space,working_dir):
    nd = para_to_scan["nVar"] # number of dimension
    samp_idx = 0
    
    if nd==1:
        for i in range(para_to_scan["nSamp"][0]):
            tempinput = copy.deepcopy(input_original)
            tempinput["Sample_idx"] = samp_idx
            tempinput = generate_input(tempinput,para_to_scan,i,samp_idx)
            samp_idx+=1
            input_space.append(tempinput)
    if nd==2:
        for i in  range(para_to_scan["nSamp"][0]):
            for j in  range(para_to_scan["nSamp"][1]):
                tempinput = copy.deepcopy(input_original)
                tempinput = generate_input(tempinput,para_to_scan,i,samp_idx)
                samp_idx+=1
                input_space.append(tempinput)
    return input_space

# Generate the input parameters based on the parameter which we want to modify.
def generate_input(tempinput,para_to_scan,i,samp_idx):
    # Default input
    type_of_particle =  1.0 
    csv_out = 0
    dynamicOn = 0
    n_per_show = 10
    turn_start = 0
    
    mainRF = 0
    main_detune = 0
    detune_slow_factor = 1.0 
    R = 100e3/2/pi
    GMTSQ = 1/1.48e-5#90089.0664
    Gamma = 45.5e9/E0Elec; 
    print("Gamma0 = ",Gamma)
    nBeam = 1
    beam_shift = 94
    
    #----------------------------#
    #important inputs
    n_turns = 100
    n_dynamicOn = 40
    n_fill = 10
    n_q_ramp = 20
    n_detune_start = 10
    n_detune_ramp = 30
    n_detune_ramp_tot = 30
    n_I_ramp_start = 10
    n_I_ramp_end = 30
    step_store = 10
    
    n_bunches = 12196
    Prad = 30e6
    t_rad_long=  0.425 
    Ek_damp = 1
    Ek_damp_always = 1
    Npar = 1440
    
    beta = np.sqrt(1-1/Gamma**2)
    T0 = 2*np.pi*R/(clight*beta)
    f0 = 1/T0
    IbDC = 0.838
    NperBunch = IbDC/n_bunches/f0/1.6e-19
    N_bins = 65 
    fill_step = 16
    siglong = 0.038e-2*Gamma #11.368 
    A = 1.0 
    n_ini_CBI = 1
    mu = np.array([0])
    CBI_ini_amp = np.array([0])
    
    nRF = 3
    nRF1 = 3
    nRFc = 0
    nRF2 = 0
    nHOM = 0
    # the max Eacc = 20MV/m, so max Vcav = 20*0.23=4.6, 
    # this will determine the ratio between focusing and defocusing cavities.40_20
    nCav = np.array([45,6,9])
    h = np.array([216816,216816,216816])
    f = h*f0
    # beam pattern related
    nTrain = 150
    gap = 0.1 # percentage of gap
    n_bunches, pattern = get_beam_pattern(nTrain, n_bunches, gap,h[0],fill_step)
    NperBunch = IbDC/n_bunches/f0/1.6e-19
    
    # convert RoQ from total to per cavity
    
    RoQ = np.array([213/4,213/4,213/4])*nCav
    RoQ = RoQ/nCav
    RoQacc = RoQ*2
    delay = np.array([0,0,0])
    n_fb_on = np.array([50.0,50,50])
    gII = np.array([10.0,10,10])
    gQQ = np.array([10.0,10,10])
    gIQ = np.array([10.0,10,10])
    gQI = np.array([10.0,10,10])
    gIIi = np.array([0.0,0,0])
    gQQi = np.array([0.0,0,0])
    gIQi = np.array([0.0,0,0])
    gQIi = np.array([0.0,0,0])
    PA_cap = np.array([1.0,1,1])
    
    #----------------------------#
    # the following parameters need to be derived from the input parameters
    # here just show some place holder.
    QL = np.array([1,0,1])
    Vref_I = np.array([1,1,1])
    Vref_Q = np.array([1,1,1])
    Iref_I = np.array([1,1,1])
    Iref_Q = np.array([1,1,1])
    I_I_ref_ini = np.array([1,1,1])
    I_I_ref_final = np.array([1,1,1])
    I_Q_ref_ini = np.array([1,1,1])
    I_Q_ref_final = np.array([1,1,1])
    detune = np.array([0.0,0,1])
    detune_ini = np.array([0.0,0,1])
    detune_mid = np.array([1,1,1])
    detune_final = np.array([1,1,1])
    
    
    
    if int(type_of_particle==2):
        atomicZ = 79
        Ek = Gamma*E0Au
    else:
        atomicZ =1
    if int(type_of_particle==1):  
        Ek = Gamma*E0Elec
    
    eta = 1/GMTSQ-1/Gamma**2
    print("eta = ",eta)
    if nRF == 1:
        Qs = np.sqrt(h[int(mainRF)]*atomicZ*np.abs(Vref_I[int(mainRF)])*eta/(2*np.pi*Ek))
    elif nRF != 1 :
        Qs = np.sqrt(h[int(mainRF)]*atomicZ*np.abs(np.sum(Vref_I))*eta/(2*np.pi*Ek))
    
    omegarf = 2*np.pi*(np.array(h)*f0)
    omegac = 2*np.pi*(np.array(h)*f0+detune_final)
    Trf = 2*np.pi/omegarf
    Rsh = [RoQ[i]*QL[i] for i in range(int(nRF))]
    Th = 2*np.pi/omegarf[0]
    dthat =Th/N_bins
    bucket_height = 2*Qs/(h[mainRF]*eta)*Gamma
    
    print("Bucket height = ",bucket_height)
    print("Ek = ",Ek)
    print("Qs = ",Qs)
    
    
    nd = para_to_scan["nVar"] # number of dimension
    idx0 = int(para_to_scan["Vars"][0][para_to_scan["Vars"][0].find("_")+1:])
    if nd==2:
        idx1 = int(para_to_scan["Vars"][1][para_to_scan["Vars"][1].find("_")+1:])
    # what really matters to the APES are the I_I and I_Q, namely the two orthogonal driving currents from the generator.
    
    
    
    thetaL = np.zeros(nRF)
    Vs = np.zeros(nRF)
    Vq = np.zeros(nRF)
    Phis = np.zeros(nRF)
    PhisPhasor = np.zeros(nRF)
    PhiIQ = np.zeros(nRF)
    PhiIQIg = np.zeros(nRF)
    VrefI = np.zeros(nRF)
    VrefQ = np.zeros(nRF)
    Vreftot = np.zeros(nRF)
    IrefI = np.zeros(nRF)
    IrefQ = np.zeros(nRF)
    QL = np.zeros(nRF)
    Rsh = np.zeros(nRF)
    
    # Need to calculate the required voltage and phase
    # then calculate the inputs for the code, namely VrefI, VrefQ, IrefI, IrefQ.
    
    # for fundamental, 
    Vtot = 0.1e9 # total voltage, just for calcualating the required voltages.
    Urad0 = Prad/IbDC # radiation caused Voltage total, depends on beam kinetic energy
    U_loss = Urad0/nCav 
    print("Urad0 : ",Urad0)
    V0 = Vtot/nCav # old voltage per cavity.
    
    # this is the required real voltage on beam, to compensate radiation loss
    Vsynch_need = U_loss
    # this is the required imaginary voltage, to provide bucket. per cavity
    Vquard_need = V0*np.sin(np.arccos(U_loss/V0))*h[0]/h[0] # calculating it from known parameters
    
    
    # new cavity voltage per cavity, 
    # once the total voltage and phasor phase of the focusing cavity is chosen, 
    # the result should be decided.
    
    if nRF == 3 :
        
        Vs[0] = 707930.6240073383
        Vq[0] = 3519260.0939708874
        Vs[1] = 303840.04162928666 
        Vq[1] = -4556823.877762136
        Vs[2] = 202560.0277528578
        Vq[2] = -3037882.5851747575
        
    if nRF == 1:
        Qs = np.sqrt(h[int(mainRF)]*atomicZ*np.abs(Vq[0]*nCav)*eta/(2*np.pi*Ek))
    elif nRF != 1 :
        Qs = np.sqrt(h[int(mainRF)]*atomicZ*np.abs(np.sum(Vq*nCav))*eta/(2*np.pi*Ek))
    
    print("Qs = ",Qs)
    PhisPhasor = np.arctan(Vq/Vs)
    
    Pbeam = IbDC*Vs
    print("Beam power per cavity: ",Pbeam)

    
    Vreftot = np.sqrt(Vs**2+Vq**2)
    Qbeam0 = Vreftot**2/(RoQacc*Pbeam)
    Qbeam = Qbeam0
    QL = Qbeam
    Rsh = RoQ*QL
    dv0 = (para_to_scan["max"][0]-para_to_scan["min"][0])/para_to_scan["nSamp"][0]
    
    ###################################################################################
    print(idx0,idx1,i,samp_idx,thetaL)
    thetaL[idx0]=para_to_scan["min"][0]+dv0*i
    print(thetaL[idx0])
    if nd==2:
        stride = para_to_scan["nSamp"][0]
        dv1 = (para_to_scan["max"][1]-para_to_scan["min"][1])/para_to_scan["nSamp"][1]
        thetaL[idx1]=para_to_scan["min"][1]+dv1*int(samp_idx%stride)
    thetaL[0]=2
    # convert thetaL to radian
    thetaL = thetaL/180*pi

    print(idx0,idx1,i,samp_idx,thetaL)
    Vbr = 2*IbDC*Rsh
    Vgr = Vreftot/np.cos(thetaL)*(1+Vbr/Vreftot*np.cos(PhisPhasor))
    ###################################################################################
    print("Vbr = ",Vbr)
    print("Vgr = ",Vgr)
    
    tgPhi = -(Vbr*np.sin(PhisPhasor)/Vreftot+(1+Vbr*np.cos(PhisPhasor)/Vreftot)*np.tan(thetaL))
    tgPhi_ini = -np.tan(thetaL)
    delta_f_ini = f*(tgPhi_ini/2/QL+np.sqrt((tgPhi_ini/2/QL)**2+1))-f
    delta_f = f*(tgPhi/2/QL+np.sqrt((tgPhi/2/QL)**2+1))-f
    
    VrefI = Vreftot*np.sin(PhisPhasor)
    VrefQ = -Vreftot*np.cos(PhisPhasor)

    I_I = Vgr/Rsh*np.sin(PhisPhasor+thetaL) 
    I_Q = -Vgr/Rsh*np.cos(PhisPhasor+thetaL)
    
    I_I_ini = Vreftot/(Rsh)/np.cos(thetaL)*np.sin(PhisPhasor+thetaL)
    I_Q_ini = -Vreftot/(Rsh)/np.cos(thetaL)*np.cos(PhisPhasor+thetaL)
    
    print("Vnew : ",Vreftot)
    print("QL : ",QL)
    print("ThetaL : ",thetaL/pi*180, " [degree]")
    
    print("Tan(PhisPhasor) : ",np.tan(PhisPhasor))
    print("PhisPhasor : ",PhisPhasor/pi*180)

    print("detune tan : ", tgPhi)
    print("detune angle : ", np.arctan(tgPhi)/pi*180, " [degree]")
    print("delta_f : ",delta_f, " [Hz]")
    print("VrefI : ",VrefI)
    print("VrefQ : ",VrefQ)
    
    print("II : ",I_I)
    print("IQ : ",I_Q)
    print("II_ini : ",I_I_ini)
    print("IQ_ini : ",I_Q_ini)

    print("VrefTot : ",Vreftot)
    print("IbDC : ", IbDC)
    
    # now write to the input file
    tempinput['type'] = np.array([type_of_particle])
    tempinput['csv_out'] = np.array([csv_out])
    tempinput['dynamicOn'] = np.array([dynamicOn])
    tempinput['n_per_show'] = np.array([n_per_show])
    tempinput['turn_start'] = np.array([turn_start])
    tempinput['mainRF'] = np.array([mainRF])
    tempinput['main_detune'] = np.array([main_detune])
    tempinput['detune_slow_factor'] = np.array([detune_slow_factor])
    tempinput['R'] = np.array([R])
    tempinput['GMTSQ'] = np.array([GMTSQ])
    tempinput['Gamma'] = np.array([Gamma])
    tempinput['nBeam'] = np.array([nBeam])
    tempinput['beam_shift'] = np.array([beam_shift])
    tempinput['nTrain'] = np.array([nTrain])
    tempinput['Pattern'] = pattern.astype(int)

    tempinput['n_turns'] = np.array([n_turns])
    tempinput['n_dynamicOn'] = np.array([n_dynamicOn])
    tempinput['n_bunches'] = np.array([n_bunches])
    tempinput['n_fill'] = np.array([n_fill])
    tempinput['n_q_ramp'] = np.array([n_q_ramp])
    tempinput['n_detune_start'] = np.array([n_detune_start])
    tempinput['n_detune_ramp'] = np.array([n_detune_ramp])
    tempinput['n_detune_ramp_tot'] = np.array([n_detune_ramp_tot])
    tempinput['n_I_ramp_start'] = np.array([n_I_ramp_start])
    tempinput['n_I_ramp_end'] = np.array([n_I_ramp_end])
    tempinput['step_store'] = np.array([step_store])
    tempinput['Prad'] = np.array([Prad])
    tempinput['t_rad_long'] = np.array([t_rad_long])
    tempinput['Ek_damp'] = np.array([Ek_damp])
    tempinput['Ek_damp_always'] = np.array([Ek_damp_always])
    tempinput['Npar'] = np.array([Npar])
    tempinput['NperBunch'] = np.array([NperBunch])
    tempinput['N_bins'] = np.array([N_bins])
    tempinput['fill_step'] = np.array([fill_step])
    tempinput['siglong'] = np.array([siglong])
    tempinput['A'] = np.array([A])
    tempinput['n_ini_CBI'] = np.array([n_ini_CBI])
    tempinput['mu'] = np.array(mu)
    tempinput['CBI_ini_amp'] = np.array(CBI_ini_amp)

    tempinput['nRF'] = np.array([nRF])
    tempinput['nRF1'] = np.array([nRF1])
    tempinput['nRFc'] = np.array([nRFc])
    tempinput['nRF2'] = np.array([nRF2])
    tempinput['nHOM'] = np.array([nHOM])
    tempinput['nCav'] = np.array(nCav)
    tempinput['h'] = np.array(h)
    tempinput['I_Q_ref_ini'] = np.array(I_Q_ref_ini)
    tempinput['I_Q_ref_final'] = np.array(I_Q_ref_final)
    tempinput['RoQ'] = np.array(RoQ)*nCav
    tempinput['delay'] = np.array(delay)
    tempinput['n_fb_on'] = np.array(n_fb_on)
    tempinput['gII'] = np.array(gII)
    tempinput['gQQ'] = np.array(gQQ)
    tempinput['gIQ'] = np.array(gIQ)
    tempinput['gQI'] = np.array(gQI)
    tempinput['gIIi'] = np.array(gIIi)
    tempinput['gQQi'] = np.array(gQQi)
    tempinput['gIQi'] = np.array(gIQi)
    tempinput['gQIi'] = np.array(gQIi)
    tempinput['PA_cap'] = np.array(PA_cap)
    
    tempinput['QL'] = np.array(QL)
    tempinput['Vref_I'] = np.array(VrefI)*nCav
    tempinput['Vref_Q'] = np.array(VrefQ)*nCav
    tempinput['Iref_I'] = np.array(I_I)
    tempinput['Iref_Q'] = np.array(I_Q)
    tempinput['I_I_ref_ini'] = np.array(I_I_ini)
    tempinput['I_I_ref_final'] = np.array(I_I)
    tempinput['I_Q_ref_ini'] = np.array(I_Q_ini)
    tempinput['I_Q_ref_final'] = np.array(I_Q)
    tempinput['detune'] = np.array(detune)
    tempinput['detune_ini'] = np.array(delta_f_ini)
    tempinput['detune_mid'] = np.array((delta_f-delta_f_ini)/2)
    tempinput['detune_final'] = np.array(delta_f)
    
    return tempinput

def write_input(input_space,home_dir):
    for i in range(len(input_space)):
        work_dir = os.path.join(home_dir,str(i))
        try:
            os.mkdir(work_dir)
        except OSError:
            print ("Creation of the directory %s failed" % work_dir)
        else:
            print ("Successfully created the directory %s" % work_dir)
            os.chdir(work_dir)
            fn1 = 'input.txt'
            inputfile1 = os.path.join(work_dir,fn1)
            with open(inputfile1,'w') as wrt_to_input:
                for j in input_space[i]:
                    wrt_to_input.write(j+' ')
                    for k in range(len(input_space[i][j])):
                        wrt_to_input.write(str(input_space[i][j][k])+' ')
                        #print(tempinput[i][j])
                    wrt_to_input.write('\n')
                wrt_to_input.close()
                print("Generated input file.")
    os.chdir(home_dir)
def run_APES(input_space,home_dir):
    args = ("../runavx2")
    for i in range(len(input_space)):
        work_dir = os.path.join(home_dir,str(i))
        try:
            os.chdir(work_dir)
        except OSError:
            print ("Changing of the directory %s failed" % work_dir)
        else:
            print ("Successfully changed to the directory %s" % work_dir)
            popen = subprocess.Popen(args, stdout=subprocess.PIPE,cwd=work_dir)
            print("Simulation started...")
            err = popen.wait()
            output = popen.stdout.read()
            print(output.decode("utf-8"))
            
# 
input_fn = "/media/tianmu/data/temp/APES/Counter_Phasing/CEPC_d1d2_g10_Train150/input.txt"
tempinput = {}
tempinput = get_input(tempinput,input_fn)

home_dir = "/media/tianmu/data/temp/APES/Counter_Phasing/CEPC_d1d2_g10_Train150/"

para_to_scan = {
    "nVar" : 2, # number of variables
    "Vars" : ["thetaL_1","thetaL_2"],
    "nSamp" : [1,1],
    "min" : [-10,-10],
    "max" : [10,10]
    }

input_space = []
input_space = get_scan_para_space(para_to_scan,tempinput,input_space,home_dir)
write_input(input_space,home_dir)
run_APES(input_space,home_dir)
