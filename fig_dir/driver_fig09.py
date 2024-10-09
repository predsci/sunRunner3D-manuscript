import os
import sys
import pysunrunner
import pysunrunner.pload as pp
import pysunrunner.io as io
import pysunrunner.pviz as pviz

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#run name
myrun = 'run09f-2'

# list of variables that will be plotted
var_list = ['vx1', 'Bx1', 'rho']

# variables that will be plotted with R^2 scaling
var_scaled = ['Bx1', 'rho']

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

# Filename for plot

file_name = png_dir +'figure09_'+myrun+'.png'

# read the dbl.out file and retrieve time information from it

wdir = pluto_path

nlinf = io.nlast_info(w_dir=wdir, datatype='dbl')

time = io.read_time(w_dir=wdir,datatype='dbl')
time = np.array(time, dtype=np.float32)

# Time in hours - conversion factor for time from PLUTO units to hours
time_fac_pluto = 1.49597871e+08/3600 #time_unit_str = '(hours)'; %1AU/1km/s
time_h = time * time_fac_pluto

# conversion factor from PLUTO magnetic field units to nT

b_fac_pluto = 0.0458505

# units are [km/s], N[cm^-3], and will not be converted 

v_fac_pluto, rho_fac_pluto = 1.0, 1.0

# set the zero of time to the end of the relaxation == CME initiation
time_h = time_h - time_h[1]

# For figure 09 plots are made 20 hours after CME initiation

time_fig = 20.0

# phi cut is at 295 which is center of CME

phi_cut = np.deg2rad(295.0)

# find the PLUTO dump index at time_fig
time_idx = np.argmin(np.abs(time_h- time_fig))

# Load PLUTO results for this time point
D = pp.pload(time_idx,w_dir=wdir,datatype='dbl')

subplot_kw = {'projection': "polar"}

cmap_list = ['rainbow', 'coolwarm', 'rainbow']

title_list = ['(a) Radial Velocity', '(b) Scaled Radial Magnetic Field', '(c) Scaled Density']

units_list = [v_fac_pluto, b_fac_pluto, rho_fac_pluto]

zmin_list = [200, -100, 0]

zmax_list = [2300, 100, 150]


# For the radial rays 
#

# Code for the radial rays
r_coords = np.array(D.x1)
t_coords = np.array(D.x2)
# Convert from co-latitude to latitude

t_coords = np.pi / 2 - t_coords

i_th_slice = 85 + 5*np.arange(18)

th_slice_all = t_coords[i_th_slice]

n_slice = len(th_slice_all)

indices = np.linspace(0, 1, n_slice)

# End of radial rays section

# read space craft location file 

spacecraft_file = 'results/Solar-MACH_2022-09-06_12-00-00.csv'
df_sc = pd.read_csv(spacecraft_file)

spacecraft_list = df_sc['Spacecraft/Body'].values
spacecraft_lat  = df_sc['Carrington latitude (°)'].values
spacecraft_dist = df_sc['Heliocentric distance (AU)'].values

abbreviations = {
    'STEREO A': 'STA',
    'Earth': 'Earth',
    'BepiColombo': 'Bepi',
    'Parker Solar Probe': 'PSP',
    'Solar Orbiter': 'SolO'
}

sc_colors = {
    'STA': 'red',
    'Earth': 'green',
    'Bepi': 'orange',
    'PSP': 'purple',
    'SolO': 'dodgerblue' #'cornflowerblue'   
}

# 
# 
# Create a figure and a 1x3 grid of subplots
fig, axs = plt.subplots(1,3, subplot_kw=subplot_kw, figsize=(10, 5))

# Iterate over the variable list and create a polar plot for each item
for ii, var_name in enumerate(var_list):
    
    if var_name in var_scaled:
        r_scale = True
    else: 
        r_scale = False

    axs[ii] = pviz.plot_phi_cut(D=D, var_name = var_name,
        phi_cut = phi_cut, ax=axs[ii],cmap = cmap_list[ii], title = title_list[ii],
        r_scale = r_scale, zmin = zmin_list[ii], zmax = zmax_list[ii], conversion_units = units_list[ii])

    # Start adding radial rays
    cmap = plt.colormaps[cmap_list[ii]]
    colors = cmap(indices)
    jcount = 0
    for th_slice in th_slice_all:
        axs[ii].plot([th_slice, th_slice], [np.min(r_coords), np.max(r_coords)], color=colors[jcount], linestyle='-', linewidth=1)
        jcount = jcount + 1
    # End of adding radial rays

    # adding spacecraft locations with markers and labels
    icount = 0
    for spacecraft, abbrev in abbreviations.items():
        sc_t = np.deg2rad(spacecraft_lat[icount])
        sc_r   = spacecraft_dist[icount]
        
        # Plot a marker at the specified (sc_phi, sc_r)
        axs[ii].plot(sc_t, sc_r, marker = 'o', color = sc_colors[abbrev], markersize=7)  # 'ko' means a black marker (circle)

        # Add text label next to the marker
        axs[ii].text(sc_t, sc_r, abbrev, fontsize=10, ha='right', va='bottom', color=sc_colors[abbrev])

        icount = icount + 1

plt.tight_layout()
plt.subplots_adjust(bottom = -0.1)

plt.savefig(file_name, dpi=150)
print("Saving Plot to: ", file_name)

# to plot to screen uncomment the line below
# plt.show()



