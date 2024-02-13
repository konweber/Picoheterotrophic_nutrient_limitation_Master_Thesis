import numpy as np

def B_prog(B_data, LDOC_data, Fe_data, z, dt, dz, µ0, m_l, m_q_mol, K_LDOC, K_Fe, w_POC, f_B_part):
    B_data_new = np.zeros_like(B_data)
    i = np.arange(1,len(z))
    if len(B_data.shape) == 1:
        B_data_new[i] = ( # biomass concentration in molC/m3
            B_data[i] # previous biomass
            + dt * (µ0 * B_data[i] * np.minimum(LDOC_data[i] / (LDOC_data[i] + K_LDOC), Fe_data[i] / (Fe_data[i] + K_Fe)) # growth # Add limitation by Fe
                                                    - m_l * B_data[i] - m_q_mol * B_data[i]**2) # mortality
                                                    + dt/dz * w_POC[i] * f_B_part * (B_data[i-1] - B_data[i]) # sinking of biomass associated with particles
                                                    )
        B_data_new[0] = (
            B_data[0] # previous biomass
            + dt * (µ0 * B_data[0] * (LDOC_data[0] / (LDOC_data[0] + K_LDOC)) # growth # Add limitation by Fe
                                                    - m_l * B_data[0] - m_q_mol * B_data[0]**2) # mortality
                                                    - dt/dz * w_POC[0] * f_B_part * B_data[0] # sinking of biomass associated with particles
        )

    # else:
    #     j = np.arange(0, len(B_data.shape[1])) # index for the different prokaryotic species
    #     B_data_new[i, j] = (
    #         B_data[i, j] # previous biomass
    #         + dt * (µ0[j] * B_data[i, j] * (LDOC_data[i, j] / (LDOC_data[i, j] + K_LDOC[j])) # growth
    #                                                 - m_l[j] * B_data[i, j] - m_q_mol[j] * B_data[i, j]**2) # mortality
    #                                                 + dt/dz * w_POC[i] * f_B_part[j] * (B_data[i-1, j] - B_data[i, j]) # sinking of biomass associated with particles
    #     )

    #     B_data_new[0, j] = (
    #         B_data[0, j] # previous biomass
    #         + dt * (µ0[j] * B_data[0, j] * (LDOC_data[0, j] / (LDOC_data[0, j] + K_LDOC[j])) # growth
    #                                                 - m_l[j] * B_data[0, j] - m_q_mol[j] * B_data[0, j]**2) # mortality
    #                                                 - dt/dz * w_POC[i] * f_B_part[j] * B_data[0, j] # sinking of biomass associated with particles

    #     )
    
    return B_data_new

### Just adapted for one prokaryotic species!
def LDOC_prog(B_data, LDOC_data, Fe_data, POC_data, z, dt, µ0, m_l, m_q_mol, K_LDOC, K_Fe, kappa, v_remin, phi, CUE_t0):
    LDOC_data_new = np.zeros_like(LDOC_data)
    i = np.arange(0,len(z))
    if len(B_data.shape) == 1:
        LDOC_data_new[i] = ( # LDOC concentration in molC/m3
            LDOC_data[i] + dt * 
            (-µ0 * B_data[i] * np.minimum((LDOC_data[i] / (LDOC_data[i] + K_LDOC)), (Fe_data[i] / (Fe_data[i] + K_Fe))) * 1/CUE_t0[i] # consumtion of LDOC by growth
                                                            + kappa * ( m_q_mol * B_data[i]**2 + m_l * B_data[i]) # production by bacterial mortality
                                                            + v_remin * POC_data[i] * phi) # production by POC remineralization
                                                            )
    
    # else:
    #     LDOC_data_new[i] = ( # LDOC concentration in molC/m3
    #         LDOC_data[i] + dt * 
    #         (np.sum([µ0[j] * B_data[i, j] * (LDOC_data[i, j] / (LDOC_data[i, j] + K_LDOC[j])) * 1/CUE_t0[i, j] for j in len(B_data.shape[1])]) # consumtion of LDOC by growth
    #         + kappa * ( np.sum([m_q_mol[j] * B_data[i, j]**2 for j in len(B_data.shape[1])]) + np.sum([m_l[j] * B_data[i, j] for j in len(B_data.shape[1])]) ) # production by bacterial mortality
    #         + v_remin * POC_data[i] * phi) # production by POC remineralization
    #         )

    return LDOC_data_new

def POC_prog(POC_data, POC_z0, z, dt, dz, w_POC, v_remin, t_step):
    POC_data_new = np.zeros_like(POC_data)
    i = np.arange(1,len(z))
    POC_data_new[i] = ( # POC concentration in molC/m3
        POC_data[i] # previous POC
          + dt/dz * (w_POC[i] * (POC_data[i-1] - POC_data[i])) # sinking of POC
            - dt * v_remin * POC_data[i] # remineralization of POC
        )
    # upper boundary condition, POC coming from the mixed layer
    POC_data_new[0] = (
        POC_data[0] # previous POC
          + POC_z0[t_step] # POC input from mixed layer
            - dt/dz * w_POC[0] * POC_data[0] # sinking of POC
              - dt * v_remin * POC_data[0] # remineralization of POC
        )
    
    return POC_data_new


def Fe_prog(Fe_data, B_data, LDOC_data, POC_data, Lig_data, Feprime_data, tau_Fe, K_LDOC, K_Fe, µ0, m_l, m_q_mol, rB_Fe_C, kscav, rPOC_Fe_C, z, dt, v_remin, Fe_z0, t_step):
    Fe_data_new = np.zeros_like(Fe_data)
    i = np.arange(1,len(z))
    Fe_data_new[i] = ( # Fe concentration in molFe/m3
        Fe_data[i] + dt * # previous Fe
            (-µ0 * B_data[i] * np.minimum((LDOC_data[i] / (LDOC_data[i] + K_LDOC)), (Fe_data[i] / (Fe_data[i] + K_Fe))) * rB_Fe_C
             + (m_l * B_data[i] + m_q_mol * B_data[i]**2) * rB_Fe_C
             - kscav * Feprime_data[i]
             + v_remin * POC_data[i] * rPOC_Fe_C
             - (1/tau_Fe) * np.maximum(0, Fe_data[i] - Lig_data[i])
            )
    )
    # upper boundary condition, Fe coming from the mixed layer
    Fe_data_new[0] = (
        Fe_data[0] + dt * # previous Fe
            (-µ0 * B_data[0] * np.minimum((LDOC_data[0] / (LDOC_data[0] + K_LDOC)), (Fe_data[0] / (Fe_data[0] + K_Fe))) * rB_Fe_C +
             (m_l * B_data[0] + m_q_mol * B_data[0]**2) * rB_Fe_C
             - kscav * Feprime_data[0]
             + v_remin * POC_data[0] * rPOC_Fe_C
             + Fe_z0[t_step]
            - (1/tau_Fe) * np.maximum(0, Fe_data[0] - Lig_data[0])
            )
    )

    return Fe_data_new

def Lig_prog(Lig_data, B_data, LDOC_data, POC_data, Fe_data, gamma, pi, rCLig, flinldoc, µ0, m_l, m_q_mol, K_LDOC, K_Fe, philig, v_remin, CUE_t0, dt, z):
    Lig_data_new = np.zeros_like(Lig_data)
    i = np.arange(0,len(z))
    Lig_data_new[i] = ( # Lig concentration in molC/m3
        Lig_data[i] + dt * # previous Lig
            (gamma * (m_l * B_data[i] + m_q_mol * B_data[i]**2)
            + pi * B_data[i]
            - µ0 * B_data[i] * np.minimum((LDOC_data[i] / (LDOC_data[i] + K_LDOC)), (Fe_data[i] / (Fe_data[i] + K_Fe))) * 1/CUE_t0[i] * Lig_data[i] * rCLig * flinldoc / LDOC_data[i]
            + v_remin * POC_data[i] * philig
            )
    )

    return Lig_data_new

def calc_Feprime(Lig_data, Fe_data, beta, z):
    Feprime_data = np.zeros_like(Fe_data)
    i = np.arange(0,len(z))
    Feprime_data[i] = ( -( 1 + (Lig_data[i] - Fe_data[i]) * beta) + np.sqrt((1 + (Lig_data[i] - Fe_data[i]) * beta)**2 + 4 * Fe_data[i] * beta) ) / (2 * beta)

    return Feprime_data

