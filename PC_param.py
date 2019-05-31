pc_param = dict()

#Conductances for all the channels with the same order as in the template

#SOMA
pc_param['eleak'] = -63
pc_param['LeakSoma'] = 1.1E-3
pc_param['Cav3.1Soma'] = 7e-6
pc_param['Cav2.1Soma'] = 2.2e-4 
pc_param['HCNSoma'] = 0.0004
pc_param['Nav1.6Soma'] = 0.214
pc_param['Kv3.4Soma'] = 0.05
pc_param['Kv1.1Soma'] = 0.002 
pc_param['Cav3.2Soma'] = 0.0008 
pc_param['Kca3.1Soma'] = 0.01 
pc_param['Cav3.3Soma'] = 0.0001 
pc_param['PC_KirSoma'] = 0.00003 
pc_param['Kca1.1Soma'] = 0.01 
pc_param['Kca2.2Soma'] = 1e-3 


#DEND
pc_param['Cav2.1Dend'] = 1e-3 
pc_param['Kca1.1Dend'] = 3.5e-2
pc_param['Kv4.3Dend'] = 0.001
pc_param['Kv1.1Dend'] = 0.0012 
pc_param['Kv1.5Dend'] = 0.13195e-3
pc_param['Kv3.3Dend'] = 0.01 
pc_param['Cav3.3Dend'] = 0.0001 
pc_param['Cav3.2Dend'] = 0.0012 
pc_param['Kca3.1Dend'] = 0.002
pc_param['Cav3.1Dend'] = 5e-6 
pc_param['Kca2.2Dend'] = 1e-3
pc_param['PC_KirDend'] = 0.00001
pc_param['Nav1.6Dend'] = 0.016
pc_param['HCNDend'] = 0.000004 

#AIS
pc_param['Cav3.1Ais'] = 8.2e-6
pc_param['Nav1.6AIS'] = 0.50 
pc_param['Cav2.1AIS'] = 2.2e-4 
pc_param['Kv3.4AIS'] = 0.01 

#AISK
pc_param['Kv1.1AisK'] = 0.01 

#First node Of Ranvier
pc_param['Nav1.6Nor'] = 0.03 
pc_param['Kv3.4Nor'] = 0.02  
pc_param['Cav3.1Nor'] = 1e-5 
pc_param['Cav2.1Nor'] = 2.2e-4 

#Second node Of Ranvier
pc_param['Nav1.6Nor2'] = 0.03 
pc_param['Kv3.4Nor2'] = 0.02  
pc_param['Cav3.1Nor2'] = 1e-5 
pc_param['Cav2.1Nor2'] = 2.2e-4  

#Third node Of Ranvier
pc_param['Nav1.6Nor3'] = 0.03
pc_param['Kv3.4Nor3'] = 0.02  
pc_param['Cav3.1Nor3'] = 1e-5 
pc_param['Cav2.1Nor3'] = 2.2e-4 

#Axon collateral
pc_param['Nav1.6Axoncoll'] = 0.03 
pc_param['Kv3.4Axoncoll'] = 0.02
pc_param['Cav3.1Axoncoll'] = 1e-5 
pc_param['Cav2.1Axoncoll'] = 2.2e-4 
