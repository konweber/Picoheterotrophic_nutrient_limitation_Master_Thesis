##################################################
### NOT EXTENSIVELY TESTED, MIGHT CONTAIN BUGS ###
##################################################

### Import modules
import numpy as np
import time

start_time = time.time()

### import all variables from namelist (init) file
from Column_model_v2_init import *

from Column_model_v2_prognostics import *


B_data = np.zeros((len(t),len(z))) # biomass (molC/m3)
LDOC_data = np.zeros((len(t),len(z))) # LDOC (molC/m3)
POC_data = np.zeros((len(t),len(z))) # POC (molC/m3)
Fe_data = np.zeros((len(t),len(z))) # Fe (molFe/m3)
Lig_data = np.zeros((len(t),len(z))) # Lig (mol Lig/m3)
Feprime_data = np.zeros(len(z))
### First time step
B_data[0,:] = B_t0
LDOC_data[0,:] = LDOC_t0
POC_data[0,:] = POC_t0
Fe_data[0,:] = Fe_t0
Lig_data[0,:] = Lig_t0



###############################################################
### Time stepping #############################################
###############################################################

for t_step in range(0,len(t)-1):
    ### Prokaryotic prognostic variables
    Feprime_data[:] = calc_Feprime(Lig_data[t_step,:], Fe_data[t_step,:], beta, z)
    B_data[t_step+1,:] = B_prog(B_data[t_step,:], LDOC_data[t_step,:], Fe_data[t_step,:],
                                 z, dt, dz, µ0, m_l, m_q_mol, K_LDOC, K_Fe, w_POC, f_B_part)
    LDOC_data[t_step+1,:] = LDOC_prog(B_data[t_step,:], LDOC_data[t_step,:], Fe_data[t_step,:], POC_data[t_step,:],
                                       z, dt, µ0, m_l, m_q_mol, K_LDOC, K_Fe, kappa, v_remin, phi, CUE_t0)
    POC_data[t_step+1,:] = POC_prog(POC_data[t_step,:], POC_z0, z, dt, dz, w_POC, v_remin, t_step)
    Fe_data[t_step+1,:] = Fe_prog(Fe_data[t_step,:], B_data[t_step,:], LDOC_data[t_step,:], POC_data[t_step,:], Feprime_data[:], K_LDOC, K_Fe, µ0, m_l, m_q_mol, rB_Fe_C, kscav, rPOC_Fe_C, z, dt, v_remin, Fe_z0, t_step)
    Lig_data[t_step+1,:] = Lig_prog(Lig_data, B_data, LDOC_data, POC_data, Fe_data, gamma, pi, rCLig, flinldoc, µ0, m_l, m_q_mol, K_LDOC, K_Fe, philig, v_remin, CUE_t0, dt, z)
    
print("Time stepping finished")

Prokar_abundance_data = B_data * 12 * 10**(-3) / q_c

MSE = np.sum((Prokar_abundance_data[-1,:] - herndl_B_abdunance_fit)**2) / len(z)
ME = np.sum(np.abs(Prokar_abundance_data[-1,:] - herndl_B_abdunance_fit)) / len(z)
print(f"MSE: {MSE:.6f}")
print(f"ME: {ME:.6f}")
### Until here the results are the same as for the automatic parameter estimation in the loop

### Save data
np.save(f'Output_data/B_data_{file_name_ext}', B_data)
np.save(f'Output_data/LDOC_data_{file_name_ext}', LDOC_data)
np.save(f'Output_data/POC_data_{file_name_ext}', POC_data)
np.save(f'Output_data/Fe_data_{file_name_ext}', Fe_data)
np.save(f'Output_data/Lig_data_{file_name_ext}', Lig_data)

print("Data saved")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.6f} seconds")
