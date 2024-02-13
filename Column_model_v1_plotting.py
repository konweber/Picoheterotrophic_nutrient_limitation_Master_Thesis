### Plotting
from tkinter import font
from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from matplotlib.animation import PillowWriter
import ffmpeg
import time

### Only necessary for animations
### Can be removed if not used
ffmpeg.FFMPEG_BINARY = 'C:/Users/konst/Downloads/ffmpeg-master-latest-win64-gpl/ffmpeg-master-latest-win64-gpl/bin'

# Enable the use of LaTeX for rendering text labels
mpl.rcParams['text.usetex'] = True
# Set the font family to sans-serif (Helvetica, Arial, etc.)
mpl.rcParams['font.family'] = 'sans-serif'
# Set the font style for the labels (optional)
mpl.rcParams['font.style'] = 'normal'
# Add a LaTeX preamble to change font family for numbers to sans-serif
mpl.rcParams['text.latex.preamble'] = r'\usepackage{sfmath}'

plt.rcParams['animation.writer'] = 'ffmpeg'

### Import data
from Column_model_v1_init import *
#from Column_model_v1_main import * # need to change that, only import variables!
from Column_model_v1_main import Prokar_abundance_data, LDOC_data, POC_data



### Intial conditions
fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True, dpi=300)
axs = axs.flatten()
axs[0].invert_yaxis()
axs[0].set_ylabel(r'Depth [m]')
axs[0].set_xlabel(r'Biomass [mol C $\cdot$ m$^{-3}$]')
axs[0].plot(B_t0, z)
axs[1].set_xlabel(r'LDOC [mol C $\cdot$ m$^{-3}$]')
axs[1].plot(LDOC_t0, z)
axs[2].set_xlabel(r'POC [mol C $\cdot$ m$^{-3}$]')
axs[2].plot(POC_t0, z)

plt.savefig(f'Plots/Initial_conditions/initial_conditions_{file_name_ext}.png', dpi=300)


### CUE profile
fig = plt.figure(figsize=(10,10), dpi=300)
plt.plot(CUE_t0, z)
plt.gca().invert_yaxis()
plt.xlabel(r'CUE [-]')
plt.ylabel(r'Depth [m]')
plt.savefig(f'Plots/Initial_conditions/initial_CUE_profile_{file_name_ext}.png', dpi=300)

### POC input from mixed layer
fig = plt.figure(figsize=(10,10), dpi=300)
plt.plot(t, POC_z0)
plt.xlabel(r'Time [days]')
plt.ylabel(r'POC additional concentration [mol C $\cdot$ m$^{-3}$]') 
plt.savefig(f'Plots/Boundary_conditions/upperlevel_POC_z0_{file_name_ext}.png', dpi=300)


###########################################################################################
### Time series animation #################################################################
###########################################################################################

start_time = time.time()
red_factor = 50
fig, axs = plt.subplots(figsize=(12, 8), ncols=4, sharey=True, dpi=300)

# Customize each subplot as needed

axs[0].invert_yaxis()

axs[0].set_xlabel(r'Prokaryotic abundance' + '\n' + r'[cells $\cdot$ m$^{-3}$]')
axs[1].set_xlabel(r'LDOC' + '\n' + r'[µmol C $\cdot$ m$^{-3}$]')
axs[2].set_xlabel(r'POC flux' + '\n' + r'[mg C $\cdot$ m$^{-2}$ $\cdot$ day$^{-1}$]')
axs[3].set_xlabel(r'POC' + '\n' + r'[mg C $\cdot$ m$^{-3}$]')

axs[0].set_ylabel(r'Depth [m]')

time_text = fig.text(0.5, 0.92, '', transform=fig.transFigure, ha="center", fontsize=12)

def init():
    line1, = axs[0].plot([], [], lw=2)
    line2, = axs[1].plot([], [], lw=2)
    line3, = axs[2].plot([], [], lw=2)
    line4, = axs[3].plot([], [], lw=2)
    return line1, line2, line3, line4

def update(frame):
    axs[0].clear()
    axs[1].clear()
    axs[2].clear()
    axs[3].clear()
    
    axs[0].invert_yaxis()

    axs[0].set_xlabel(r'Prokaryotic abundance' + '\n' + r'[cells $\cdot$ m$^{-3}$]')
    axs[1].set_xlabel(r'LDOC' + '\n' + r'[µmol C $\cdot$ m$^{-3}$]')
    axs[2].set_xlabel(r'POC flux' + '\n' + r'[mg C $\cdot$ m$^{-2}$ $\cdot$ day$^{-1}$]')
    axs[3].set_xlabel(r'POC' + '\n' + r'[mg C $\cdot$ m$^{-3}$]')

    axs[0].set_ylabel(r'Depth [m]')
    time_text.set_text(f'Time Step: {frame * red_factor}')

    lines = []
    
    lines.append(axs[0].plot(Prokar_abundance_data[0, :], z + z0, color='gray', alpha=0.8)[0])
    lines.append(axs[0].plot(Prokar_abundance_data[frame * red_factor, :], z + z0, color='limegreen')[0])
    
    #lines.append(axs[1].plot(LDOC_data[0, :], z + z0, color='gray')[0])
    lines.append(axs[1].plot(LDOC_data[frame * red_factor, :] * 1E6, z + z0, color='orangered')[0])

    #lines.append(axs[2].plot(POC_data[0, :] * w_POC, z + z0, color='gray')[0])
    lines.append(axs[2].plot(POC_data[frame * red_factor, :] * w_POC * 10E3 * 12, z + z0, color='royalblue')[0])

    #lines.append(axs[3].plot(POC_data[0, :], z + z0, color='gray')[0])
    lines.append(axs[3].plot(POC_data[frame * red_factor, :] * 1E3 * 12, z + z0, color='purple')[0])
    
    return lines


num_frames = len(t) // red_factor

animation = FuncAnimation(fig, update, frames=np.arange(0, num_frames), init_func=init, blit=True)

print("Animation finished")

animation.save(f'Plots/Animations/animation_{file_name_ext}.mp4', writer="ffmpeg", fps=10)

print("Animation saved")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")

#HTML(animation.to_jshtml())
#HTML(animation.to_html5_video())