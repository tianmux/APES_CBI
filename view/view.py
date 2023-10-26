import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import os
import numpy as np
from array import array
from numpy import pi,arange,cos,sin,sqrt,exp,log10
from matplotlib.figure import Figure

clight = 299792458
E0Au = 196.9665687*931.5e6
E0Elec = 0.51099895000e6
clight = 299792458

def set_working_directory():
    directory = entry.get()  # get text from entry field
    os.chdir(directory)  # set working directory
    print("Current working directory:", os.getcwd())

def open_file():
    inputfile = 'input.txt'
    
    inputfile = os.path.join(os.getcwd(),inputfile)
    tempinput = {}
    try:
        with open(inputfile) as fn:
            for line in fn:
                if len(line.split())>1:
                    tempinput[line.split()[0]] = line.split()[1:]
            for i in tempinput:
                for j in range(len(tempinput[i])):
                    tempinput[i][j] = float(tempinput[i][j])
        fn.close()
        nRF = int(tempinput['nRF'][0])
        nRF1 = int(tempinput['nRF1'][0])
        nRF2 = int(tempinput['nRF2'][0])
        nRFc = int(tempinput['nRFc'][0])

        E0Au = 196.9665687*931.5e6
        E0Elec = 0.51099895000e6
        nTurns = int(tempinput['n_turns'][0])
        nDynamic = int(tempinput['n_dynamicOn'][0])
        nfill = int(tempinput['n_fill'][0])
        n_q_ramp = int(tempinput['n_q_ramp'][0])
        NpRF = int(tempinput['N_bins'][0])
        h = [int(i) for i in tempinput['h']]
        detune_ini = np.array([i for i in tempinput['detune_ini']])
        detune_final = np.array([i for i in tempinput['detune_final']])

        step = int(tempinput['step_store'][0])
        fill_step = int(tempinput['fill_step'][0])
        nBeam = int(tempinput['nBeam'][0])
        beam_shift = int(tempinput['beam_shift'][0])
        nBunch = int(tempinput['n_bunches'][0])
        nGap = h[0]/fill_step-nBunch
        nPar = int(tempinput['Npar'][0])
        NperBunch = int(tempinput['NperBunch'][0])
        nTot = nBunch*nPar*nBeam
        Gamma0 = tempinput['Gamma'][0]
        Rring = tempinput['R'][0]
        n_record = nTurns/step
        clight = 299792458
        beta = np.sqrt(1-1/Gamma0**2)
        T0 = 2*np.pi*Rring/(clight*beta)
        f0 = 1/T0
        frf = f0*np.array(h)
        V0 = [i for i in tempinput['Vref_I']]
        V0Q = [i for i in tempinput['Vref_Q']]
        II = [i for i in tempinput['Iref_I']]
        IQ = [i for i in tempinput['Iref_Q']]
        mainRF = int(tempinput['mainRF'][0])
        t_rad_long = tempinput['t_rad_long'][0]
        nTrain = int(tempinput['nTrain'][0])
        Pattern = np.array(tempinput['Pattern']).astype(int)

        if int(tempinput['type'][0]==2):
            atomicZ = 79
            Ek = Gamma0*E0Au
        else:
            atomicZ =1
        if int(tempinput['type'][0]==1):  
            Ek = Gamma0*E0Elec
            
        GMTSQ = tempinput['GMTSQ'][0]
        Ek_damp = tempinput['Ek_damp'][0]

        eta = 1/GMTSQ-1/Gamma0**2
        if nRF ==1:
            Qs = np.sqrt(h[mainRF]*atomicZ*np.abs(V0[mainRF])*eta/(2*np.pi*Ek))
        elif nRF == 2 :
            Qs = np.sqrt(h[mainRF]*atomicZ*np.abs(V0[mainRF]+V0[1])*eta/(2*np.pi*Ek))
            
        elif nRF == 3 :
            Qs = np.sqrt(h[mainRF]*atomicZ*np.abs(V0[mainRF]+V0[1]+V0[2]*h[2]/h[0])*eta/(2*np.pi*Ek))
        omega0 = f0*2*pi
        omegarf = 2*np.pi*(np.array(h)*f0)
        omegac = 2*np.pi*(np.array(h)*f0+detune_final)
        Trf = 2*np.pi/omegarf

        RoQ0 = [i for i in tempinput['RoQ']]

        RoQ = np.array(RoQ0)/1
        QL0 = [i for i in tempinput['QL']]
        QL = np.array(QL0)*1
        R = [RoQ[i]*QL[i] for i in range(nRF)]

        Th = 2*np.pi/omegarf[0]
        dthat =Th/NpRF

        pattern = 'd'+'dd'*nBeam+3*nRF*'d'
        n_stride = 1+2*nBeam+3*nRF
        stride = len(pattern)*8
        test = array('d')
        bucket_height = 2*Qs/(h[mainRF]*eta)


        print(Qs)
        print(bucket_height)
        return tempinput
    except: 
        print('input file is not found')
    
    
def plot_first_bunch_centroid():
    global M1_1,M1_2
    tempinput = open_file()
    nTurns = int(tempinput['n_turns'][0])
    
    nBunch = int(tempinput['n_bunches'][0])
    nBeam = int(tempinput['nBeam'][0])
    Gamma0 = tempinput['Gamma'][0]
    h = [int(i) for i in tempinput['h']]
    Rring = tempinput['R'][0]
    beta = np.sqrt(1-1/Gamma0**2)
    T0 = 2*np.pi*Rring/(clight*beta)

    
    M1_1_0 = []
    M1_2_0 = []
    for i in range(nTurns):
        M1_1_0.append(M1_all[i*nBunch*2])
        M1_2_0.append(M1_all[i*nBunch*2+nBunch]-Gamma0)
    
    global canvas
    if canvas is not None:
        canvas.get_tk_widget().destroy()

    # Create a new Figure and Canvas
    fig = Figure(figsize=(6,5), dpi=100)
    ax = fig.add_subplot(111)

    turn_start = int(entry_start_n.get())  # get text from entry field
    n_turn_disp = int(entry_n_plot.get())
    #turns = 15000
    #n_turn_disp = 500#int(1/Qs)*7000
    rng1 = turn_start
    rng2 = rng1+n_turn_disp
    
    ax.plot(np.array(M1_1_0[rng1:rng2])/T0*h[0]*360,np.array(M1_2_0[rng1:rng2]),'rx-')
    ax.plot(np.array(M1_1_0[rng1:rng1+10])/T0*h[0]*360,np.array(M1_2_0[rng1:rng1+10]),'gx-')
    ax.plot(np.array(M1_1_0[rng2-10:rng2])/T0*h[0]*360,np.array(M1_2_0[rng2-10:rng2]),'bx-')
    ax.set_ylabel('first_centroid(Gamma)',fontsize=10)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack()

def get_plot_idx(Pattern,turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain):
    bunch_idx = 0
    for i in range(train_start):
        bunch_idx += Pattern[i*2]
    bunch_idx+=bunch_start
    start = turn_start*nBunch+bunch_idx
    for i in range(train_end-train_start):
        bunch_idx += Pattern[(train_start+i)*2]
    bunch_idx+=bunch_end-bunch_start
    end = (turn_end)*nBunch+bunch_idx
    return start,end
def get_centroids(Pattern,Trf,nBunch,t0,fill_step):
    centroids = np.ndarray(nBunch)
    bunch_idx = 0
    bucket_idx = 0
    nTrain = int(len(Pattern)/2)
    for i in range(nTrain):
        for j in range(Pattern[i*2]):
            centroids[bunch_idx]=Trf*fill_step*bucket_idx+t0
            bunch_idx+=1
            bucket_idx+=1
        bucket_idx+=Pattern[i*2+1]
    print(bunch_idx,bucket_idx)
    return centroids

def plot_dphi():
    tempinput = open_file()
    nTurns = int(tempinput['n_turns'][0])
    mainRF = int(tempinput['mainRF'][0])
    fill_step = int(tempinput['fill_step'][0])

    nBunch = int(tempinput['n_bunches'][0])
    nBeam = int(tempinput['nBeam'][0])
    nTrain = int(tempinput['nTrain'][0])

    Gamma0 = tempinput['Gamma'][0]
    h = [int(i) for i in tempinput['h']]
    Rring = tempinput['R'][0]
    beta = np.sqrt(1-1/Gamma0**2)
    T0 = 2*np.pi*Rring/(clight*beta)
    Pattern = np.array(tempinput['Pattern']).astype(int)
    f0 = 1/T0
    detune_final = np.array([i for i in tempinput['detune_final']])
    f = f0*np.array(h)
    omega0 = f0*2*pi
    omegarf = 2*np.pi*(np.array(h)*f0)
    omegac = 2*np.pi*(np.array(h)*f0+detune_final)
    Trf = 2*np.pi/omegarf

    turns = np.array([i for i in range(nTurns*nBeam*nBunch)])

    turn_start = int(entry_start_n_phi.get())  # get text from entry field
    turn_end = int(entry_n_plot_phi.get())-1+turn_start
    #turn_start = 15000
    #turn_end = 15000
    train_start = 0
    train_end = 0
    bunch_start = 0
    bunch_end = nBunch
    bunch_idx = 0

    centroids = get_centroids(Pattern,Trf[mainRF],nBunch,Trf[mainRF]/2,fill_step)#[(i*fill_step+0.5)*Trf[mainRF] for i in range(nBunch)]
    rng1,rng2 = get_plot_idx(Pattern,turn_start,turn_end,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
    rng3,rng4 = get_plot_idx(Pattern,0,0,train_start,train_end,bunch_start,bunch_end,nBunch,nTrain)
    print("Read M1_all.bin done!")
    print(rng1,rng2,rng3,rng4)
    global canvas
    if canvas is not None:
        canvas.get_tk_widget().destroy()

    # Create a new Figure and Canvas
    fig = Figure(figsize=(6,5), dpi=100)
    ax = fig.add_subplot(111)

    ax.plot((np.array(M1_1[rng1:rng2])-centroids[rng3:rng4])*f[0]*360,'r.')

    ax.set_ylabel('centroid (d_t)',fontsize=10)
    ax.set_xlabel('Bunch # ',fontsize=10)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack()

def read_M1s():
    try:
        tempinput = open_file()
        nTurns = int(tempinput['n_turns'][0])
        mainRF = int(tempinput['mainRF'][0])
        fill_step = int(tempinput['fill_step'][0])

        nBunch = int(tempinput['n_bunches'][0])
        nBeam = int(tempinput['nBeam'][0])
        nTrain = int(tempinput['nTrain'][0])

        Gamma0 = tempinput['Gamma'][0]
        h = [int(i) for i in tempinput['h']]
        Rring = tempinput['R'][0]
        beta = np.sqrt(1-1/Gamma0**2)
        T0 = 2*np.pi*Rring/(clight*beta)
        Pattern = np.array(tempinput['Pattern']).astype(int)
        f0 = 1/T0
        detune_final = np.array([i for i in tempinput['detune_final']])
        f = f0*np.array(h)
        omega0 = f0*2*pi
        omegarf = 2*np.pi*(np.array(h)*f0)
        omegac = 2*np.pi*(np.array(h)*f0+detune_final)
        Trf = 2*np.pi/omegarf
        
        global M1_all
        M1_fn = 'M1_all.bin'
        cwd = os.getcwd()
        datafile = os.path.join(cwd,M1_fn)    
        with open(datafile, mode='rb') as file: # b is important -> binary
            M1_all.fromfile(file,2*nTurns*nBeam*nBunch)
        global M1_1,M1_2
        for i in range(nTurns):
            for j in range(nBunch):
                M1_1.append(M1_all[i*nBunch*2+j])
                M1_2.append(M1_all[i*nBunch*2+nBunch+j]-Gamma0)
        
        print("Read M1_all.bin done!")
    except:
        print("Error: Please open file first!")
        return

root = tk.Tk()

# Create text entry field
default_value = tk.StringVar()
default_value.set("/home/txin/Documents/Work_Folder/Code/APES/unit_Test/BEPCII/20230726_recon")
entry = tk.Entry(root,textvariable=default_value, width=50)
entry.pack()

# Create button to set the current working directory
button_setDir = tk.Button(root, text="Set Working Directory", command=set_working_directory)
button_setDir.pack()

# Create button to open file
button_open = tk.Button(root, text="Open File", command=open_file)
button_open.pack()

# Create button to read data
M1_all = array('d')
M1_1 = []
M1_2 = []
button_read = tk.Button(root, text="Read M1s", command=read_M1s)
button_read.pack()

# Create button to plot M1s

label = tk.Label(root, text="Starting Turn Number")  # Create a new Label widget
label.pack()
entry_start_n = tk.Entry(root,)  # Create a new Entry widget
entry_start_n.insert(0, "15000")
entry_start_n.pack()

label = tk.Label(root, text="Number of Turns to plot")  # Create a new Label widget
label.pack()
entry_n_plot = tk.Entry(root)  # Create a new Entry widget
entry_n_plot.insert(0, "100")
entry_n_plot.pack()

button = tk.Button(root, text="Plot first bunch centroids", command=plot_first_bunch_centroid)
button.pack()

label = tk.Label(root, text="Starting Turn Number")  # Create a new Label widget
label.pack()
entry_start_n_phi = tk.Entry(root,)  # Create a new Entry widget
entry_start_n_phi.insert(0, "15000")
entry_start_n_phi.pack()

label = tk.Label(root, text="Number of Turns to plot")  # Create a new Label widget
label.pack()
entry_n_plot_phi = tk.Entry(root)  # Create a new Entry widget
entry_n_plot_phi.insert(0, "1")
entry_n_plot_phi.pack()

button = tk.Button(root, text="Plot d_phis", command=plot_dphi)
button.pack()

# Declare the canvas as a global variable
canvas = None
fig = Figure(figsize=(6,5), dpi=100)
ax = fig.add_subplot(111)
canvas = FigureCanvasTkAgg(fig, master=root)
#canvas.get_tk_widget().pack()

# Get screen width and height
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

# Calculate width and height for Tk root window
window_width = int(screen_width / 2)
window_height = int(screen_height / 2)

# Set the dimension of the screen and where it is placed
root.geometry(f'{window_width}x{window_height}+{window_width//2}+{window_height//2}')

root.mainloop()

