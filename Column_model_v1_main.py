### Import modules
import numpy as np
import time

start_time = time.time()

### import all variables from namelist (init) file
from Column_model_v1_init import *

from Column_model_v1_prognostics import *


B_data = np.zeros((len(t),len(z))) # biomass (molC/m3)
LDOC_data = np.zeros((len(t),len(z))) # LDOC (molC/m3)
POC_data = np.zeros((len(t),len(z))) # POC (molC/m3)
### First time step
B_data[0,:] = B_t0
LDOC_data[0,:] = LDOC_t0
POC_data[0,:] = POC_t0


###############################################################
### Time stepping #############################################
###############################################################

for t_step in range(0,len(t)-1):
    ### Prokaryotic prognostic variables
    B_data[t_step+1,:] = B_prog(B_data[t_step,:], LDOC_data[t_step,:], z, dt, dz, µ0, m_l, m_q_mol, K_LDOC, w_POC, f_B_part)
    LDOC_data[t_step+1,:] = LDOC_prog(B_data[t_step,:], LDOC_data[t_step,:], POC_data[t_step,:], z, dt, µ0, m_l, m_q_mol, K_LDOC, kappa, w_POC, v_remin, phi, CUE_t0)
    POC_data[t_step+1,:] = POC_prog(POC_data[t_step,:], POC_z0, z, dt, dz, w_POC, v_remin, t_step)

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

print("Data saved")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Elapsed time: {elapsed_time:.6f} seconds")
