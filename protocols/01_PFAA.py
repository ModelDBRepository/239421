# 01 - Excitatory and inhibitory activity on Purkinje cell model
# Protocols to reproduce all the images based on the excitation/inhibition burst/pause behaviors
import matplotlib as mpl
mpl.use('tkagg')  
from Purkinje_py3 import Purkinje_py3
from neuron import h
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

#fixed time step only
Fixed_step = h.CVode()
Fixed_step.active(0) #the model does not work with the variable time step!

#Instantiation of the cell template
cell = Purkinje_py3()

#this code discover the number of cores available in a CPU and activate the multisplit to use them all.
cores = multiprocessing.cpu_count()
h.load_file("parcom.hoc")
p = h.ParallelComputeTool()
p.change_nthread(cores,1)
p.multisplit(1)
print ('cores', cores)

stimdata = dict()
stimdata['timeglobal'] =  5000

synapsesdata = dict()
#number of parallel fibers synapses, 1 for each section. MAX 1111
synapsesdata['npf'] = 100

#number of ascending axon synapses, 1 for each section. MAX 383
synapsesdata['naa'] = 0

#number of stellate cell synapses, 1 for each section. MAX 1332
synapsesdata['nstl'] = 25

#Excitation
#parallel fiber
synapsesdata['syninterval'] = 2
synapsesdata['synnumber'] = 10
synapsesdata['synstart'] = 600
synapsesdata['synnoise'] = 0

#ascending axon
synapsesdata['synaainterval'] = 2
synapsesdata['synaanumber'] = 10
synapsesdata['synaastart'] = 600
synapsesdata['synaanoise'] = 0

#Inhibition
#Stellate on parallel fiber
synapsesdata['synstlinterval'] = 7
synapsesdata['synstlnumber'] = 3
synapsesdata['synstlstart'] = 600
synapsesdata['synstlnoise'] = 0

# delay factor
synapsesdata['synpfdelay'] = 0
synapsesdata['synaadelay'] = 0
synapsesdata['synstldelay'] = 4

#Neuron control menu
h.nrncontrolmenu()

#Voltage graph
h('load_file("vm.ses")')




cell.createsyn((int(synapsesdata['npf'])), (int(synapsesdata['naa'])), (int(synapsesdata['nstl'])))

#PF bursts
spk_stimpf = []
totalstim_pf = int(stimdata['timeglobal']/  synapsesdata['synstart'])

for j in range(int(totalstim_pf)):
    spk_stim_pf = h.NetStim()
    spk_stim_pf.interval=synapsesdata['syninterval']
    spk_stim_pf.number=synapsesdata['synnumber']
    spk_stim_pf.noise=synapsesdata['synnoise']
    spk_stim_pf.start=(synapsesdata['synstart'] * (totalstim_pf - j)) + synapsesdata['synpfdelay']
    
    spk_stimpf.append(spk_stim_pf)
    spk_nc_pfsyn = []
    j = j-1

print('len pf', len(cell.PF_syn))

for m in range(int(totalstim_pf)):	
    spk_nc_pfsyn.append([h.NetCon(spk_stimpf[m],PF.input,0,0.1,1) for PF in cell.PF_syn])



#Ascending Axon.
spk_stimaa = []
totalstim_aa = int(stimdata['timeglobal']/  synapsesdata['synaastart'])

for j in range(int(totalstim_aa)):
    spk_stim_aa = h.NetStim()
    spk_stim_aa.interval=synapsesdata['synaainterval']
    spk_stim_aa.number=synapsesdata['synaanumber']
    spk_stim_aa.noise=synapsesdata['synaanoise']
    spk_stim_aa.start=(synapsesdata['synaastart'] * (totalstim_aa - j)) + synapsesdata['synaadelay']
    
    spk_stimaa.append(spk_stim_aa)
    spk_nc_aasyn = []
    j = j-1

print('len aa', len(cell.AA_syn))

for m in range(int(totalstim_aa)):	
    spk_nc_aasyn.append([h.NetCon(spk_stimaa[m],AA.input,0,0.1,1) for AA in cell.AA_syn])
    
    
#Stellate cell
spk_stim_stl = []
totalstim_sc = int(stimdata['timeglobal']/  synapsesdata['synstlstart'])

for j in range(int(totalstim_sc)):
    spk_stim_sc = h.NetStim()
    spk_stim_sc.interval=synapsesdata['synstlinterval']
    spk_stim_sc.number=synapsesdata['synstlnumber']
    spk_stim_sc.noise=synapsesdata['synstlnoise']
    spk_stim_sc.start=(synapsesdata['synstlstart'] * (totalstim_sc - j)) + synapsesdata['synstldelay']
    
    spk_stim_stl.append(spk_stim_sc)
    spk_nc_stlsyn = []
    j = j-1

print('len stl', len(cell.stl_syn))

for m in range(int(totalstim_sc)):	
    spk_nc_stlsyn.append([h.NetCon(spk_stim_stl[m],stl.input,0,0.1,1) for stl in cell.stl_syn])
      

#Basic properties of the simulation. dt, temperature, sim duration and initial voltage
h.dt = 0.025
h.celsius = 37
h.tstop = stimdata['timeglobal']
h.v_init = -65

#initialization and run.    
def initialize():
    h.finitialize()
    h.run()
    
initialize()

#save files
np.savetxt('01_vm_soma.txt', np.column_stack((np.array(cell.rec_t), np.array(cell.vm_soma))), delimiter = ' ')

img = plt.plot(np.array(cell.rec_t), np.array(cell.vm_soma))
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.savefig('01_vm_soma.eps')
plt.close()
