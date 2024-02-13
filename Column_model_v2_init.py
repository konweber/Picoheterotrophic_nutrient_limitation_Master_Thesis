import numpy as np
import json

### File name extension
file_name_ext = 'test_version2'

### Define constants
### Numerical setup
dt = 0.1 # time step (days)
t_max = 1000 # maximum time (days)
dz = 50 # depth step (m)
z_max = 5000 # maximum depth (m)
### Set up vertical grid
### z denotes the upper boundary of each grid box
z = np.arange(0, z_max, dz) # depth grid (m)
t = np.arange(0,t_max,dt)

herndl_B_abdunance_fit = (z + 225)**(-0.680975835) * 4.97015566E11 # fit determined from Herndl et al. (2023)
# + 200 from z0 (depth of the euphotic zone) and +25 for the middle of the grid box

### Prokaryotic parameters
µ0 = 10**(-4) * 3600 * 24 # maximum growth rate (1/day) # a little more than 10^-5 /s
m_l = 0.002194278163009 # 10**(-8) * 3600 * 24 # linear mortality rate (1/day) # Zakem et al. (2021)
m_q = 1.95860584347654E-12 # 10**(-16) * 3600 * 24 / 1000 # 5 * 10**(-15) # quadratic mortality rate (m3/cell/day) # 5.8 * 10^-17 L/cell/s
q_c = 10 * 10**(-18) # carbon quota (kg C/cell)
### f_B_part is problematic, because there is no clear upper boundary condition!
f_B_part = 0 # fraction of biomass associated with particles (dimensionless)
CUE = 0.1 # Carbon use efficiency (dimensionless)
# CUE_depth_dependent = False # If True, CUE is a function of depth
# CUE_min = 0.01 # Minimum CUE
# CUE_max = 0.15 # Maximum CUE
K_LDOC = 10**(-5) * 1000 # half-saturation constant for LDOC (mol C/m3)
kappa = 0.5 # fraction of dead prokaryotic biomass that yields LDOC (dimensionless)
K_Fe = 10**(-9) * 1000 # half-saturation constant for Fe (mol Fe/m3)

### Fe parameters
rB_Fe_C = 4E-5 # ratio of Fe to C in prokaryotic biomass (mol Fe/mol C)
kscav = 1E-3 # scavenging rate (1/day)
rPOC_Fe_C = 1E-5 # ratio of Fe to C in POC (mol Fe/mol C)

### Ligand parameters
rCLig = 30 # number of carbon atoms per ligand molecule (dimensionless)
beta = 10E9 # binding constant (m3/mol)
gamma = kappa / rCLig * 1E-2 # fraction of dead prokaryotic biomass that yields LDOC (mol Lig/mol C)
pi = 0

flinldoc = 0.5 # fraction of ligand that is LDOC (dimensionless)

### POC parameters
v_remin = 0.04 # remineralization rate (1/day); need to verify with literature; in the range of Uchimiya et al. (2018)
w_POC_z0 = 154 # sinking rate of POC (m/day); this is pretty low to be honest
dw_POC_dz = 0.025 # gradient of sinking rate of POC (m/day/m); need to verify with literature
phi = 0.5 # fraction of POC that yields LDOC (dimensionless)
philig = phi / rCLig * 1E-2 # fraction of POC that yields Ligands (mol Lig/mol C)
z0 = 200 # depth of the euphotic zone (m), used to calculate intial POC profile
b_t0 = 0.7 # Martin curve parameter for initial POC profile (dimensionless)
F_z0_t0 = 20 * 10**(-3) / 12 # surface POC flux (mol C/m2/day)
### AR(1) process parameters
# np.random.seed(123)
# phi = 0.99
# sigma = 10**(-3) / 12 # standard deviation of the noise (mol C/m2/day)

CUE_t0 = np.ones_like(z) * CUE
### POC input from mixed layer
POC_z0 = np.ones_like(t) * F_z0_t0 * dt/dz # amount of POC in the upper most grid box that is added in each time step (mol C/m3)
Fe_z0 = 0 * np.ones_like(t) # amount of Fe in the upper most grid box that is added in each time step (mol Fe/m3)
### Calculated parameters
m_q_mol = m_q/q_c * 12 * 10**(-3) # quadratic mortality rate (m3/mol C/day)


### Calculate vertical POC velocity for each level
w_POC = np.ones_like(np.arange(0,z_max,dz)) * w_POC_z0 + dw_POC_dz * z # sinking velocity of POC (m/day)

### Exponentially decaying biomass
B_t0 = ((z + 25 + z0)**(-0.680975835) * 4.97015566E11) * q_c / (12 * 10**(-3)) # initial biomass concentration (molC/m3)
LDOC_t0 = 5 * 10**(-5) * np.ones_like(z) # np.minimum(10**(-6) * np.ones_like(z), np.ones_like(z) * F_z0_t0 * (b_t0)/z0 * ((z + dz/2 + z0)/z0)**(-b_t0 - 1) * phi) # initial LDOC concentration (molC/m3), Martin curve profile for the flux, differentiated
POC_t0 = np.ones_like(z) * F_z0_t0 * ((z + dz/2 + z0)/z0)**(-b_t0) / w_POC_z0 # initial POC concentration (molC/m3), Martin curve profile for the flux, differentiated
Fe_t0 = np.ones_like(z) * 0.2 * 10E-9 * 1000 # initial Fe concentration (molFe/m3)
Lig_t0 = np.ones_like(z) * 1 * 10E-9 * 1000 # initial Ligand concentration (mol Lig/m3)

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
    'K_Fe': K_Fe,
    'kappa': kappa,
    'rB_Fe_C': rB_Fe_C,
    'kscav': kscav,
    'rPOC_Fe_C': rPOC_Fe_C,
    'rCLig': rCLig,
    'beta': beta,
    'gamma': gamma,
    'pi': pi,
    'flinldoc': flinldoc,
    'v_remin': v_remin,
    'w_POC_z0': w_POC_z0,
    'dw_POC_dz': dw_POC_dz,
    'phi': phi,
    'philig': philig,
    'z0': z0,
    'b_t0': b_t0,
    'F_z0_t0': F_z0_t0,
    'B_t0': B_t0.tolist(), # convert to list to be able to save as json
    'LDOC_t0': LDOC_t0.tolist(),
    'POC_t0': POC_t0.tolist(),
    'Fe_t0': Fe_t0.tolist(),
    'Lig_t0': Lig_t0.tolist(),
    'CUE_t0': CUE_t0.tolist(),
    'POC_z0': POC_z0.tolist(),
    'Fe_z0': Fe_z0.tolist(),
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


