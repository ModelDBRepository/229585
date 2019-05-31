# 05 - Calcium spikes and sodium bursts.

from Purkinje import Purkinje
from neuron import h
import multiprocessing
import numpy as np

#fixed time step only
Fixed_step = h.CVode()
Fixed_step.active(0) #the model does not work with the variable time step!

#Instantiation of the cell template
cell = Purkinje()

stimdata = dict()
stimdata['stim0del'] = 1000
stimdata['stim0dur'] = 4000
stimdata['stim0amp'] = 2

stimdata['timeglobal'] = 4000

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

stim = [h.IClamp(0.5,sec=cell.soma)] 

stim[0].delay = stimdata['stim0del']
stim[0].dur = stimdata['stim0dur']
stim[0].amp = stimdata['stim0amp']  

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
np.savetxt('05_vm_soma.txt', timevolt_soma, delimiter = ' ')
np.savetxt('05_vm_NOR3.txt', timevolt_NOR3, delimiter = ' ')
