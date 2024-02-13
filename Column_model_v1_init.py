import numpy as np
import json

### File name extension
file_name_ext = 'best_param_100mg_updated_run'

### Define constants
### Numerical setup
dt = 0.1 # time step (days)
t_max = 2500 # maximum time (days)
dz = 50 # depth step (m)
z_max = 5000 # maximum depth (m)
z0 = 200 # depth of the euphotic zone (m), used to calculate intial POC profile
### Set up vertical grid
### z denotes the upper boundary of each grid box
z = np.arange(z0, z_max, dz) # depth grid (m)
t = np.arange(0,t_max,dt)

herndl_B_abdunance_fit = (z + 25)**(-0.680975835) * 4.97015566E11 # fit determined from Herndl et al. (2023)
# + 200 from z0 (depth of the euphotic zone) and +25 for the middle of the grid box

### Prokaryotic parameters
n_species = 1 # number of prokaryotic species
µ0 = 10**(-4) * 3600 * 24 # maximum growth rate (1/day) # a little more than 10^-5 /s
m_l = 0.000471199305939 # 10**(-8) * 3600 * 24 # linear mortality rate (1/day) # Zakem et al. (2021) (log10(m_l) ~ -3.5 - -1.5)
m_q = 1.70850330829918E-11 # 10**(-16) * 3600 * 24 / 1000 # 5 * 10**(-15) # quadratic mortality rate (m3/cell/day) # 5.8 * 10^-17 L/cell/s (log10(m_q) ~ -15 - -10)
q_c = 10 * 10**(-18) # carbon quota (kg C/cell)
### f_B_part is problematic, because there is no clear upper boundary condition!
f_B_part = 0 # fraction of biomass associated with particles (dimensionless)
CUE = 0.1 # Carbon use efficiency (dimensionless)
# CUE_depth_dependent = False # If True, CUE is a function of depth
# CUE_min = 0.01 # Minimum CUE
# CUE_max = 0.15 # Maximum CUE
K_LDOC = 10**(-5) * 1000 # half-saturation constant for LDOC (mol C/m3)
kappa = 0.5 # fraction of dead prokaryotic biomass that yields LDOC (dimensionless)

### POC parameters
v_remin = 0.087454901051186 # remineralization rate (1/day); need to verify with literature; in the range of Uchimiya et al. (2018)
w_POC_z0 = 5.96238766293015 # sinking rate of POC (m/day); this is pretty low to be honest
dw_POC_dz = 0.057083546535192# gradient of sinking rate of POC (m/day/m); need to verify with literature
phi = 0.528359089919721 # fraction of POC that yields LDOC (dimensionless)
b_t0 = 0.7 # Martin curve parameter for initial POC profile (dimensionless)
###################################################################
F_z0_t0 = 80 * 10**(-3) / 12 # surface POC flux (mol C/m2/day)
#################################################################
### AR(1) process parameters
# np.random.seed(123)
# phi = 0.99
# sigma = 10**(-3) / 12 # standard deviation of the noise (mol C/m2/day)

CUE_t0 = np.ones_like(z) * CUE
### POC input from mixed layer
POC_z0 = np.ones_like(t) * F_z0_t0 * dt/dz # amount of POC in the upper most grid box that is added in each time step (mol C/m3)

### Calculated parameters
m_q_mol = m_q/q_c * 12 * 10**(-3) # quadratic mortality rate (m3/mol C/day)


### Calculate vertical POC velocity for each level
w_POC = np.ones_like(z) * w_POC_z0 + dw_POC_dz * z # sinking velocity of POC (m/day)

### Exponentially decaying biomass
B_t0 = ((z + 25)**(-0.680975835) * 4.97015566E11) * q_c / (12 * 10**(-3)) # initial biomass concentration (molC/m3)
LDOC_t0 = 5 * 10**(-5) * np.ones_like(z) # np.minimum(10**(-6) * np.ones_like(z), np.ones_like(z) * F_z0_t0 * (b_t0)/z0 * ((z + dz/2 + z0)/z0)**(-b_t0 - 1) * phi) # initial LDOC concentration (molC/m3), Martin curve profile for the flux, differentiated
POC_t0 = np.ones_like(z) * F_z0_t0 * ((z + dz/2 + z0)/z0)**(-b_t0) / w_POC_z0 # initial POC concentration (molC/m3), Martin curve profile for the flux, differentiated

# if CUE_depth_dependent == False:
#     CUE_t0 = np.ones_like(z) * CUE 
# elif CUE_depth_dependent == True:
#     CUE_t0 = CUE_max * 10**(z/z_max * 1/np.log10(CUE_min/CUE_max)) # CUE decrease exponentially with depth



# t_eq = t_max // dt #1000 // dt # time for system to equilibrate (days)
# ar1_t = np.linspace(0, int(len(t) - t_eq), int(len(t) - t_eq))
# for i in range(1,len(ar1_t)):
#     noise = np.random.normal(0,sigma)
#     ar1_t[i] = phi * ar1_t[i-1] + noise

# season_t = (F_z0_t0 + 3 * F_z0_t0/4 * np.sin(2 * np.pi * dt * np.linspace(0, int(len(t) - t_eq), int(len(t) - t_eq)) / 365))
# forcing_t = (np.maximum(0, ar1_t + season_t)).tolist()
# POC_z0 = np.array(int(t_eq) * [F_z0_t0] + forcing_t) * dt / dz # Flux times time step times area divided by volume

param_dict = { 
    'file_name_ext': file_name_ext,
    'dt': dt,
    't_max': t_max,
    'dz': dz,
    'z_max': z_max,
    'µ0': µ0,
    'm_l': m_l,
    'm_q': m_q,
    'q_c': q_c,
    'f_B_part': f_B_part,
    'CUE': CUE,
    # 'CUE_depth_dependent': CUE_depth_dependent,
    # 'CUE_min': CUE_min,
    # 'CUE_max': CUE_max,
    'K_LDOC': K_LDOC,
    'kappa': kappa,
    'v_remin': v_remin,
    'w_POC_z0': w_POC_z0,
    'dw_POC_dz': dw_POC_dz,
    'phi': phi,
    'z0': z0,
    'b_t0': b_t0,
    'F_z0_t0': F_z0_t0,
    'B_t0': B_t0.tolist(), # convert to list to be able to save as json
    'LDOC_t0': LDOC_t0.tolist(),
    'POC_t0': POC_t0.tolist(),
    'CUE_t0': CUE_t0.tolist(),
    'POC_z0': POC_z0.tolist(),
    'w_POC': w_POC.tolist(),
    't': t.tolist(),
    'z': z.tolist()
    # 't_eq': t_eq,
    # 'ar1_t': ar1_t.tolist(),
    # 'season_t': season_t.tolist(),
    # 'forcing_t': forcing_t
}

with open(f'Parameter_dicts/param_dict_{file_name_ext}.json', 'w') as fp:
    json.dump(param_dict, fp)


