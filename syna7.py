#1 checks for diff a7 density/neuron
#2 checks for all R density combos
#3 checks for diff nt conc combos
# Import important libraries
from neuron import h,gui
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import itertools
import os
import math
from scipy.interpolate import interp1d

# Load model hoc files
h.load_file("stdlib.hoc")
h.load_file("stdrun.hoc")
h.load_file("driver.hoc")
h.load_file("cholinergic_neuron.hoc")
h.load_file("GABAfsi.hoc")
h.load_file("GABAnfsi.hoc")
h.load_file("Pyramidal_cell.hoc")

#Accessing different neurons defined in hoc 

gabafsi = h.gabafsi
gabanfsi = h.gabanfsi
chol = h.chol
pyrcell = h.pyrcell
pcell= h.pyrcell

# Accessing different ion channels defined in mod files
Ca12L = h.Ca12L
Ca13L = h.Ca13L
Ca21PQ = h.Ca21PQ
Ca22N = h.Ca22N
Ca31T = h.Ca31T
Ca33T = h.Ca33T

KV11HH = h.KV11
KV12HH = h.KV12
KV13HH = h.KV13
KV14HH = h.KV14
KV31HH = h.KV31
KV32HH = h.KV32
KV72HH = h.KV72

NaV11= h.NaV11
NaV12= h.NaV12
NaV16= h.NaV16

HCN1 = h.HCN1
HCN2 = h.HCN2
leak= h.leak
NaKpump= h.NaKpump
capump= h.capump
achrelease = h.achrelease
gabafrelease = h.gabafrelease
gabanrelease = h.gabanrelease
glurelease = h.glurelease

ampa= h.AMPA
nmda = h.NMDA
gabaa = h.GABAA
gabab = h.GABAB
#a4b2 = h.a4b2
#a3b4 = h.a3b4
a7 = h.a7

# Setting simulation parameters
h.tstop = 300
h.dt = 0.025
v_init = -65
thresholdpyrcell = -20
thresholdgabafsi = -20
thresholdgabanfsi = -25

df = pd.DataFrame(columns =['Number of receptors'])                 # creating a pandas dataframe with a column of input parameters

rec_names = ["a7_gabafsi", "a7_gabanfsi", "a7_pyrcell"]

# Defining the number of receptors for each type of receptors

num_a7gf = [0, 10000]
num_a7gn = [0, 10000]
num_a7pyr = [0, 10000]

folder_name = 'vara7density4'
os.makedirs(folder_name, exist_ok = True)
simulation = 1

num_rec = [num_a7gf, num_a7gn, num_a7pyr]
num_rec_comb = list(itertools.product(*num_rec))                         # making all the possible combination of number of different receptors

#for pyrcell
num_dend_pyr = 86                                                              
num_axon_pyr = 47

pyrcell_soma_nseg = pyrcell.soma.nseg                                           # number of segments in pyrcell soma

pyrcell_dend_nseg = []                                                       # Creating list containing number of segments in each dendrite  
for i in range(0,num_dend_pyr):
    pyrcell_dend_nseg.append(pyrcell.dend[i].nseg)    

pyrcell_axon_nseg = []                                                       # Creating list containing number of segments in each dendrite  
for i in range(0,num_axon_pyr):
    pyrcell_axon_nseg.append(pyrcell.axon[i].nseg)    

#for gabafsi
num_dend_gabafsi = 84                                                              
num_axon_gabafsi = 457

gabafsi_soma_nseg = gabafsi.soma.nseg                                           # number of segments in pyrcell soma

gabafsi_dend_nseg = []                                                       # Creating list containing number of segments in each dendrite  
for i in range(0,num_dend_gabafsi):
    gabafsi_dend_nseg.append(gabafsi.dend[i].nseg)    

gabafsi_axon_nseg = []                                                       # Creating list containing number of segments in each dendrite  
for i in range(0,num_axon_gabafsi):
    gabafsi_axon_nseg.append(gabafsi.axon[i].nseg)    

#for gabanfsi
num_dend_gabanfsi = 69                                                              
num_axon_gabanfsi = 99

gabanfsi_soma_nseg = gabanfsi.soma.nseg                                           # number of segments in pyrcell soma

gabanfsi_dend_nseg = []                                                       # Creating list containing number of segments in each dendrite  
for i in range(0,num_dend_gabanfsi):
    gabanfsi_dend_nseg.append(gabanfsi.dend[i].nseg)    

gabanfsi_axon_nseg = []                                                       # Creating list containing number of segments in each dendrite  
for i in range(0,num_axon_gabanfsi):
    gabanfsi_axon_nseg.append(gabanfsi.axon[i].nseg)    

for num_recs in num_rec_comb:                                             # Taking a particular combination of the number of receptors for each kind of receptors at one time
    num_a7gf= num_recs[0]
    num_a7gn= num_recs[1]
    num_a7pyr= num_recs[2]
    '''
    num_gabaa = num_recs[0]
    num_gabab = num_recs[1]
    num_ampa = num_recs[2]
    num_nmda = num_recs[3]
    num_a7 = num_recs[4]

#why in loop, match w original
    #a7 on pyr
    for i in range(0, num_dend_pyr):                                           # Loop over each dendrite
        num_a7pyr_dec =[]                                                  # Creating list containing how the no. of a7 receptors decreases down the length of pyrcell dendrite 
        for j in range (0, pyrcell_dend_nseg[i]):                             # Loop over every segment to implement a7 receptor
            num_a7pyr_dec.append(int(num_a7pyr - (num_a7pyr/pyrcell_dend_nseg[i])*j))
            recep_list1 =[]
            for k in range(0, num_a7pyr_dec[j]):
                recep =  h.a7()
                recep.loc(pyrcell.dend[i]((2*(j+1)-1)/(2*pyrcell_dend_nseg[i]))) # formula in bracket calculates the middle node of every segment
                h.setpointer(chol.soma(0.5)._ref_T_achrelease, 'C', recep)     # setting the pointer of a7 receptor to reference to the presynaptic transmitter concentration
                recep_list1.append(recep)

    #a7 on p
    for i in range(0, num_dend_pyr):                                           # Loop over each dendrite
        num_a7p_dec =[]                                                  # Creating list containing how the no. of a7 receptors decreases down the length of pyrcell dendrite 
        for j in range (0, pyrcell_dend_nseg[i]):                             # Loop over every segment to implement a7 receptor
            num_a7p_dec.append(int(num_a7pyr - (num_a7pyr/pyrcell_dend_nseg[i])*j))
            recep_list2 =[]
            for k in range(0, num_a7p_dec[j]):
                recep =  h.a7()
                recep.loc(pyrcell.dend[i]((2*(j+1)-1)/(2*pyrcell_dend_nseg[i]))) # formula in bracket calculates the middle node of every segment
                h.setpointer(chol.soma(0.5)._ref_T_achrelease, 'C', recep)     # setting the pointer of a7 receptor to reference to the presynaptic transmitter concentration
                recep_list2.append(recep)

    #a7 on gabafsi
    for i in range(0, num_dend_gabafsi):                                           # Loop over each dendrite
        num_a7gf_dec =[]                                                  # Creating list containing how the no. of a7 receptors decreases down the length of pyrcell dendrite 
        for j in range (0, gabafsi_dend_nseg[i]):                             # Loop over every segment to implement a7 receptor
            num_a7gf_dec.append(int(num_a7gf - (num_a7gf/gabafsi_dend_nseg[i])*j))
            recep_list3 =[]
            for k in range(0, num_a7gf_dec[j]):
                recep =  h.a7()
                recep.loc(gabafsi.dend[i]((2*(j+1)-1)/(2*gabafsi_dend_nseg[i]))) # formula in bracket calculates the middle node of every segment
                h.setpointer(chol.soma(0.5)._ref_T_achrelease, 'C', recep)     # setting the pointer of a7 receptor to reference to the presynaptic transmitter concentration
                recep_list3.append(recep)

    #a7 on gabanfsi
    for i in range(0, num_dend_gabanfsi):                                           # Loop over each dendrite
        num_a7gn_dec =[]                                                  # Creating list containing how the no. of a7 receptors decreases down the length of pyrcell dendrite 
        for j in range (0, gabanfsi_dend_nseg[i]):                             # Loop over every segment to implement a7 receptor
            num_a7gn_dec.append(int(num_a7gn - (num_a7gn/gabanfsi_dend_nseg[i])*j))
            recep_list4 =[]
            for k in range(0, num_a7gn_dec[j]):
                recep =  h.a7()
                recep.loc(gabanfsi.dend[i]((2*(j+1)-1)/(2*gabanfsi_dend_nseg[i]))) # formula in bracket calculates the middle node of every segment
                h.setpointer(chol.soma(0.5)._ref_T_achrelease, 'C', recep)     # setting the pointer of a7 receptor to reference to the presynaptic transmitter concentration
                recep_list4.append(recep)'''

    #implementing Rs on GABAfsi   
    recep_list6 = []
    for i in range(0, 100):                                           
        recep = h.GABAA()                                                    
        recep.loc(gabafsi.soma(0.5))
        h.setpointer(gabanfsi.soma(0.5)._ref_T_gabanrelease, 'C', recep)
        recep_list6.append(recep)  

    recep_list7 = []
    for i in range(0, 100):                                           
        recep = h.GABAB()                                                   
        recep.loc(gabafsi.axon[456](0.9))
        h.setpointer(gabafsi.soma(0.5)._ref_T_gabafrelease, 'C', recep)
        recep_list7.append(recep)  

    recep_list8 = []
    for i in range(0, 100):                                           
        recep = h.AMPA()                                                   
        recep.loc(gabafsi.dend[83](0.5))
        h.setpointer(pyrcell.soma(0.5)._ref_T_glurelease, 'C', recep)
        recep_list8.append(recep)  

    recep_list9 = []
    for i in range(0, 100):                                           
        recep = h.NMDA()                                                   
        recep.loc(gabafsi.dend[83](0.5))
        h.setpointer(pyrcell.soma(0.5)._ref_T_glurelease, 'C', recep)
        recep_list9.append(recep)

    recep_list1 = []
    for i in range(0, num_a7gf):                                           
        recep = h.a7()                                                   
        recep.loc(gabafsi.soma(0.5))
        h.setpointer(chol.soma(0.5)._ref_T_achrelease, 'C', recep)
        recep_list1.append(recep)      

    #implementing Rs on GABAnfsi   
    recep_list10 = []
    for i in range(0, 100):                                           
        recep = h.GABAA()                                                    
        recep.loc(gabanfsi.soma(0.5))
        h.setpointer(gabafsi.soma(0.5)._ref_T_gabafrelease, 'C', recep)
        recep_list10.append(recep)  

    recep_list11 = []
    for i in range(0, 100):                                           
        recep = h.GABAB()                                                   
        recep.loc(gabanfsi.axon[98](0.9))
        h.setpointer(gabanfsi.soma(0.5)._ref_T_gabanrelease, 'C', recep)
        recep_list11.append(recep) 

    recep_list12 = []
    for i in range(0, 100):                                           
        recep = h.AMPA()                                                   
        recep.loc(gabanfsi.dend[68](0.5))
        h.setpointer(pyrcell.soma(0.5)._ref_T_glurelease, 'C', recep)
        recep_list12.append(recep)   

    recep_list13 = []
    for i in range(0, 100):                                           
        recep = h.NMDA()                                                    
        recep.loc(gabanfsi.dend[68](0.5))
        h.setpointer(pyrcell.soma(0.5)._ref_T_glurelease, 'C', recep)
        recep_list13.append(recep) 

    recep_list2 = []
    for i in range(0, num_a7gn):                                           
        recep = h.a7()                                                   
        recep.loc(pyrcell.soma(0.5))
        h.setpointer(chol.soma(0.5)._ref_T_achrelease, 'C', recep)
        recep_list2.append(recep)       
    
    #implementing Rs on Pyrcell  
    recep_list14 = []
    #gap= num_gabaa/2
    #gbp= num_gabab/2

    for i in range(0, 100):                                           
        recep = h.GABAA()                                                   
        recep.loc(pyrcell.soma(0.5))
        h.setpointer(gabafsi.soma(0.5)._ref_T_gabafrelease, 'C', recep)
        recep_list14.append(recep)  

    recep_list15 = []
    for i in range(0, 100):                                           
        recep = h.GABAA()                                                   
        recep.loc(pyrcell.soma(0.5))
        h.setpointer(gabanfsi.soma(0.5)._ref_T_gabanrelease, 'C', recep)
        recep_list15.append(recep)  

    recep_list16 = []
    for i in range(0, 100):                                           
        recep = h.GABAB()                                                   
        recep.loc(pyrcell.dend[85](0.5))
        h.setpointer(gabafsi.soma(0.5)._ref_T_gabafrelease, 'C', recep)
        recep_list16.append(recep) 

    recep_list17 = []
    for i in range(0, 100):                                           
        recep = h.GABAB()                                                  
        recep.loc(pyrcell.dend[85](0.5))
        h.setpointer(gabanfsi.soma(0.5)._ref_T_gabanrelease, 'C', recep)
        recep_list17.append(recep) 

    recep_list18 = []
    for i in range(0, 100):                                           
        recep = h.AMPA()                                                   
        recep.loc(pyrcell.dend[83](0.5))
        h.setpointer(pcell.soma(0.5)._ref_T_glurelease, 'C', recep)
        recep_list18.append(recep)   

    recep_list19 = []
    for i in range(0, 100):                                           
        recep = h.NMDA()                                                   
        recep.loc(pyrcell.dend[83](0.5))
        h.setpointer(pcell.soma(0.5)._ref_T_glurelease, 'C', recep)
        recep_list19.append(recep) 
    
    recep_list3 = []
    for i in range(0, num_a7pyr):                                           
        recep = h.a7()                                                   
        recep.loc(pyrcell.soma(0.5))
        h.setpointer(chol.soma(0.5)._ref_T_achrelease, 'C', recep)
        recep_list3.append(recep)      

    #implementing Rs on Pcell  
    recep_list20 = []
    for i in range(0, 100):                                           
        recep = h.GABAA()                                                   
        recep.loc(pcell.soma(0.5))
        h.setpointer(gabafsi.soma(0.5)._ref_T_gabafrelease, 'C', recep)
        recep_list20.append(recep)  

    recep_list21 = []
    for i in range(0, 100):                                           
        recep = h.GABAA()                                                   
        recep.loc(pcell.soma(0.5))
        h.setpointer(gabanfsi.soma(0.5)._ref_T_gabanrelease, 'C', recep)
        recep_list21.append(recep)  

    recep_list22 = []
    for i in range(0, 100):                                           
        recep = h.GABAB()                                                   
        recep.loc(pcell.dend[85](0.5))
        h.setpointer(gabafsi.soma(0.5)._ref_T_gabafrelease, 'C', recep)
        recep_list22.append(recep) 

    recep_list23 = []
    for i in range(0, 100):                                           
        recep = h.GABAB()                                                   
        recep.loc(pcell.dend[85](0.5))
        h.setpointer(gabanfsi.soma(0.5)._ref_T_gabanrelease, 'C', recep)
        recep_list23.append(recep) 

    recep_list24 = []
    for i in range(0, 100):                                           
        recep = h.AMPA()                                                   
        recep.loc(pcell.dend[83](0.5))
        h.setpointer(pyrcell.soma(0.5)._ref_T_glurelease, 'C', recep)
        recep_list24.append(recep)   

    recep_list25 = []
    for i in range(0, 100):                                           
        recep = h.NMDA()                                                    
        recep.loc(pcell.dend[83](0.5))
        h.setpointer(pyrcell.soma(0.5)._ref_T_glurelease, 'C', recep)
        recep_list25.append(recep)  

    recep_list4 = []
    for i in range(0, num_a7pyr):                                           
        recep = h.a7()                                                   
        recep.loc(pcell.soma(0.5))
        h.setpointer(chol.soma(0.5)._ref_T_achrelease, 'C', recep)
        recep_list4.append(recep)      
     

    #print("No. of gabaa, gabab, ampa, nmda, a7 receptors in the zero end are: ", num_gabaa, num_gabab, num_ampa, num_nmda, num_a7)
    print("No. of a7gf, afgn, a7pyr receptors in the zero end are: ", num_a7gf, num_a7gn, num_a7pyr)
    #measure from where?
    voltage_vecpyrcell = h.Vector()
    voltage_vecpyrcell.record(pyrcell.axon[46](0.5)._ref_v)
    voltage_vecgabafsi = h.Vector()
    voltage_vecgabafsi.record(gabafsi.axon[456](0.5)._ref_v)
    voltage_vecgabanfsi = h.Vector()
    voltage_vecgabanfsi.record(gabanfsi.axon[98](0.5)._ref_v)
    
    
    cai_vecpyrcell = h.Vector()
    cai_vecpyrcell.record(pyrcell.soma(0.5)._ref_cai)
    '''    cai_vecgabafsi = h.Vector()
    cai_vecgabafsi.record(gabafsi.soma(0.5)._ref_cai)
    cai_vecgabanfsi = h.Vector()
    cai_vecgabanfsi.record(gabanfsi.soma(0.5)._ref_cai)'''

    spike_times_vecpyrcell = h.Vector()
    spike_times_vecgabafsi = h.Vector()
    spike_times_vecgabanfsi = h.Vector()
    
    spike_counterpyrcell = h.APCount(pyrcell.axon[46](0.5)) # Calculation of spike frequency 
    spike_counterpyrcell.thresh = thresholdpyrcell
    spike_counterpyrcell.record(spike_times_vecpyrcell)

    spike_countergabafsi = h.APCount(gabafsi.axon[456](0.5)) # Calculation of spike frequency 
    spike_countergabafsi.thresh = thresholdgabafsi
    spike_countergabafsi.record(spike_times_vecgabafsi)

    spike_countergabanfsi = h.APCount(gabanfsi.axon[98](0.5)) # Calculation of spike frequency 
    spike_countergabanfsi.thresh = thresholdgabanfsi
    spike_countergabanfsi.record(spike_times_vecgabanfsi)
    
    time_vec = []
    for i in range(1, 301):
        time_vec.append(i)
    
    h.run()
    
    spike_freqpyrcell = spike_counterpyrcell.n             # Retrieval of spike frequency
    vol_peakpyrcell = voltage_vecpyrcell.max()
    ahp_peakpyrcell = voltage_vecpyrcell.min()    
    
    voltage_valuespyrcell = list(voltage_vecpyrcell)  # Convert the recorded voltage to a Python list

    spike_timespyrcell = list(spike_times_vecpyrcell)    
    isipyrcell = np.diff(spike_timespyrcell)

    spike_freqgabafsi = spike_countergabafsi.n             # Retrieval of spike frequency
    vol_peakgabafsi = voltage_vecgabafsi.max()
    ahp_peakgabafsi = voltage_vecgabafsi.min()    
    
    voltage_valuesgabafsi = list(voltage_vecgabafsi)  # Convert the recorded voltage to a Python list

    spike_timesgabafsi = list(spike_times_vecgabafsi)    
    isigabafsi = np.diff(spike_timesgabafsi)

    spike_freqgabanfsi = spike_countergabanfsi.n             # Retrieval of spike frequency
    vol_peakgabanfsi = voltage_vecgabanfsi.max()
    ahp_peakgabanfsi = voltage_vecgabanfsi.min()    
    
    voltage_valuesgabanfsi = list(voltage_vecgabanfsi)  # Convert the recorded voltage to a Python list

    spike_timesgabanfsi = list(spike_times_vecgabanfsi)    
    isigabanfsi = np.diff(spike_timesgabanfsi)
    
    low_lim = 0
    up_lim = 300
    vol_avg = []
   
    
    while(up_lim<=h.tstop):
        vol_add = 0
        for j in range(low_lim, up_lim): 
            vol_add = vol_add + voltage_vecpyrcell.x[j]
        vol_avg.append(vol_add/100)        
        low_lim = low_lim + 100
        up_lim = up_lim + 100
        
    vol_points = []
    for i in range(100): 
        vol_points.append(voltage_vecpyrcell.x[i])
        
        
    cai_points = []
    for i in range(100): 
        cai_points.append(cai_vecpyrcell.x[i])
          

    df = df._append({'Number of receptors': list(zip(rec_names, num_recs)), 
                    'Spike frequency pyrcell': spike_freqpyrcell,
                    'Spike frequency pyrcell': spike_freqgabafsi,
                    'Spike frequency pyrcell': spike_freqgabanfsi,
                    #'Avg plots per 100 data points': vol_avg,
                    'V pyrcell': voltage_valuespyrcell,
                    'V gabafsi': voltage_valuesgabafsi,
                    'V gabanfsi': voltage_valuesgabanfsi,
                    'Vmax pyrcell': vol_peakpyrcell,
                    'Vmax gabafsi': vol_peakgabafsi,
                    'Vmax gabanfsi': vol_peakgabanfsi,
                    'AHP peak pyrcell': ahp_peakpyrcell,
                    'AHP peak gabafsi': ahp_peakgabafsi,
                    'AHP peak gabanfsi': ahp_peakgabafsi,
                    'Spike times pyrcell': spike_timespyrcell,
                    'Spike times gabafsi': spike_timesgabafsi,
                    'Spike timesgabanfsi': spike_timesgabanfsi,
                    'ISI pyrecll': isipyrcell,
                    'ISI gbabafsi': isigabafsi,
                    'ISI gabanfsi': isigabanfsi,
                    'Voltage vector': vol_points,
                    'cai points' : cai_points}, ignore_index = True)  
                    
    plt.plot(voltage_vecpyrcell)
    plt.xlabel('Time (ms)')
    plt.ylabel('Pyramidal cell voltage (mV) at axon 0.5 location')
    plt.title('Membrane potential at simulation {}'.format(simulation))
    
    plt.plot(voltage_vecgabafsi)
    plt.xlabel('Time (ms)')
    plt.ylabel('GABAfsi voltage (mV) at axon 0.5 location')
    #plt.title('GABAfsi Membrane potential at simulation {}'.format(simulation))

    plt.plot(voltage_vecgabanfsi)
    plt.xlabel('Time (ms)')
    plt.ylabel('GABAnfsi voltage (mV) at axon 0.5 location')
    #plt.title('GABAnfsi cell Membrane potential at simulation {}'.format(simulation))

    file_name = f'{folder_name}/membrane potential in simulation{simulation}.jpg'
    
    plt.savefig(file_name)
    plt.clf()
    
    simulation += 1  

df.to_excel('outputa73.xlsx')              # Converting the pandas dataframe to an excel file   	
    
