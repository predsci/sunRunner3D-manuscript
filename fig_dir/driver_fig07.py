import os
import sys
import pysunrunner
import pysunrunner.pload as pp
import pysunrunner.io as io
import pysunrunner.pviz as pviz

import numpy as np
import matplotlib.pyplot as plt

#run name
myrun = 'run09f-2'

# list of variables that will be plotted
var_list = ['vx1', 'Bx1', 'rho', 'prs']

# variables that will be plotted with R^2 scaling
var_scaled = ['Bx1', 'rho']

# variables that will be plotted on log10 scale
var_log10 = 'prs'

# path to location of PLUTO output files
local_path = "/Volumes/work1/michal/runs/cr2261-high/hmi_mast_mas_std_0201/"

pluto_path = local_path+myrun+'/output/'

# directory for saving plots
png_dir = local_path+myrun+'/png_files/'

# check if png_dir exists and if not create it
if not os.path.exists(png_dir):
    # If it doesn't exist, create it
    os.makedirs(directory)
    print(f"Directory '{png_dir}' created.")
else:
    print(f"Directory '{png_dir}' already exists.")


# read the dbl.out file and retrieve time information from it

wdir = pluto_path

nlinf = io.nlast_info(w_dir=wdir, datatype='dbl')

time = io.read_time(w_dir=wdir,datatype='dbl')
time = np.array(time, dtype=np.float32)

# Time in hours - conversion factor for time from PLUTO units to hours

time_fac_pluto = 1.49597871e+08/3600 
time_h = time * time_fac_pluto

# set the zero of time to the end of the relaxation == CME initiation

time_h = time_h - time_h[1]

# For figure 7 plots are made 20 hours after CME initiation

time_fig = 20.0

# find the PLUTO dump index at time_fig
time_idx = np.argmin(np.abs(time_h- time_fig))

# Load PLUTO results for this time point

D = pp.pload(time_idx,w_dir=wdir,datatype='dbl')

subplot_kw = {'projection': "polar"}

cmap_list = ['rainbow', 'coolwarm', 'rainbow', 'rainbow']

title_list = ['(a) Radial Velocity', '(b) Scaled Radial Magnetic Field', '(c) Scaled Density', '(d) Log$_{10}$(Pressure)']

zmin_list = [200, -100, 0, 2]

zmax_list = [2200, 100, 150, 8]

# Code for the radial rays, including periodic boundary conditions
r_coords = np.array(D.x1)
p_coords = np.array(D.x3)

i_phi_slice = 250 + 5*np.arange(18)

n_slice = len(i_phi_slice)
indices = np.linspace(0, 1, n_slice)

for ii in range(0, n_slice):
    if (i_phi_slice[ii] >= len(p_coords)):
        i_phi_slice[ii] = i_phi_slice[ii] - len(p_coords)

fi_slice_all = p_coords[i_phi_slice]
# End of radial rays code

# 
# Create a figure and a 2x2 grid of subplots
fig, axs = plt.subplots(2, 2, subplot_kw=subplot_kw, figsize=(10, 10))

# Flatten the 2x2 array of axes for easier iteration
axs = axs.flatten()

# Iterate over the variable list and create a polar plot for each item
for ii, var_name in enumerate(var_list):
    
    if var_name in var_scaled:
        r_scale = True
    else: 
        r_scale = False
    if var_name in var_log10:
        log_scale = True
    else:
        log_scale = False
    axs[ii] = pviz.plot_equatorial_cut(D = D, var_name = var_name, ax=axs[ii], cmap = cmap_list[ii], title = title_list[ii],
        r_scale = r_scale, log_scale = log_scale, zmin = zmin_list[ii], zmax = zmax_list[ii])

    # Start adding radial rays
    cmap = plt.colormaps[cmap_list[ii]]
    # cmap = plt.cm.get_cmap(cmap_list[ii])
    colors = cmap(indices)
    jcount = 0
    for fi_slice in fi_slice_all:
        axs[ii].plot([fi_slice, fi_slice], [np.min(r_coords), np.max(r_coords)], color=colors[jcount], linestyle='-', linewidth=1)
        jcount = jcount + 1
    # End of adding radial rays

plt.tight_layout()
plt.subplots_adjust(bottom = 0)
file_name = png_dir +'figure07_'+myrun+'.png'
plt.savefig(file_name, dpi=150)
print("Saving Plot to: ", file_name)

# to plot to screen uncomment the line below
#plt.show()




