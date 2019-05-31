# 01 - No sodium and calcium channels in the AIS.
# This protocol tests the absence of channels in the AIS and it is critical to validate the Purkinje cell. 

from Purkinje import Purkinje
from neuron import h
import multiprocessing
import numpy as np

#fixed time step only
Fixed_step = h.CVode()
Fixed_step.active(0) #the model does not work with the variable time step!

#Instantiation of the cell template
cell = Purkinje()

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

#To turn off the AIS CHANNELS
cell.axonAIS.pcabar_Cav3_1 = 0
cell.axonAIS.gbar_Nav1_6 = 0
cell.axonAIS.pcabar_Cav2_1 = 0

#Basic properties of the simulation. dt, temperature, sim duration and initial voltage
h.dt = 0.025
h.celsius = 37
h.tstop = 1000
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
np.savetxt('01_vm_soma.txt', timevolt_soma, delimiter = ' ')
np.savetxt('01_vm_NOR3.txt', timevolt_NOR3, delimiter = ' ')