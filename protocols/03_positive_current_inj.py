# 03 - Positive current injections from 0.1 to 1.5nA.

from Purkinje import Purkinje
from neuron import h
import multiprocessing
import numpy as np


#fixed time step only
Fixed_step = h.CVode()
Fixed_step.active(0) #the model does not work with the variable time step!

#Instantiation of the cell template
cell = Purkinje()

#Dictionary with the parameters
stimdata = dict()
stimdata['stim0del'] = 300
stimdata['stim0dur'] = 1000
stimdata['stim0amp'] = 0.1

stimdata['stim1del'] = 1300
stimdata['stim1dur'] = 1000
stimdata['stim1amp'] = 0.2

stimdata['stim2del'] = 2300
stimdata['stim2dur'] = 1000
stimdata['stim2amp'] = 0.5

stimdata['stim3del'] = 3300
stimdata['stim3dur'] = 1000
stimdata['stim3amp'] = 1

stimdata['stim4del'] = 4300
stimdata['stim4dur'] = 1000
stimdata['stim4amp'] = 1.5

stimdata['timeglobal'] = 5500

#this code discover the number of cores available in a CPU and activate the multisplit to use them all.
cores = multiprocessing.cpu_count()
h.load_file("parcom.hoc")
p = h.ParallelComputeTool()
p.change_nthread(cores,1)
p.multisplit(1)
print 'cores', cores

#Neuron control menu
h.nrncontrolmenu()

#Voltage graph
h('load_file("vm.ses")')

stim = [h.IClamp(0.5,sec=cell.soma), h.IClamp(0.5,sec=cell.soma), h.IClamp(0.5,sec=cell.soma), h.IClamp(0.5,sec=cell.soma), h.IClamp(0.5,sec=cell.soma)] 

stim[0].delay = stimdata['stim0del']
stim[0].dur = stimdata['stim0dur']
stim[0].amp = stimdata['stim0amp'] 

stim[1].delay = stimdata['stim1del']
stim[1].dur = stimdata['stim1dur']
stim[1].amp = stimdata['stim1amp'] 

stim[2].delay = stimdata['stim2del']
stim[2].dur = stimdata['stim2dur']
stim[2].amp = stimdata['stim2amp'] 

stim[3].delay = stimdata['stim3del']
stim[3].dur = stimdata['stim3dur']
stim[3].amp = stimdata['stim3amp'] 

stim[4].delay = stimdata['stim4del']
stim[4].dur = stimdata['stim4dur']
stim[4].amp = stimdata['stim4amp'] 

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

#save vectors for time and voltage
time = np.array(cell.rec_t)
vm_soma = np.array(cell.vm_soma)
vm_NOR3 = np.array(cell.vm_NOR3)

timevolt_soma = np.column_stack((time,vm_soma))
timevolt_NOR3 = np.column_stack((time,vm_NOR3))

#save files
np.savetxt('03_vm_soma.txt', timevolt_soma, delimiter = ' ')
np.savetxt('03_vm_NOR3.txt', timevolt_NOR3, delimiter = ' ')