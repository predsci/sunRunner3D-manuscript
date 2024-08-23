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

# variable that will be plotted
var_name = 'vx1'

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

# Panels (a), (b) and (c) of figure 8 are made 10, 15 and 20 hours after CME initiation
# You only need to change the time_fig and panel_fig variables

time_fig = 10.0 #10.0 15.0 20.0

panel_name = 'a' #'b' 'c'

##
# File name for figure
file_name = png_dir +'figure08_'+panel_name+'.png'

# read the dbl.out file and retrieve time information from it
wdir = pluto_path

time = io.read_time(w_dir=wdir,datatype='dbl')
time = np.array(time, dtype=np.float32)

# Time in hours - conversion factor for time from PLUTO units to hours
time_fac_pluto = 1.49597871e+08/3600 #time_unit_str = '(hours)'; %1AU/1km/s
time_h = time * time_fac_pluto

# set the zero of time to the end of the relaxation == CME initiation
time_h = time_h - time_h[1]
time_idx = np.argmin(np.abs(time_h- time_fig))

# Load PLUTO results for this time point
D = pp.pload(time_idx,w_dir=wdir,datatype='dbl')


# Plot is made as a function of radial distance, hence r_cut is set to None
r_cut = None

# Plot is made at theta == equatorial
theta_slice = 0.0

# The slices are at a series of phi values.  Same as the radial rays in Figure 7 

# For the radial rays, including periodic boundary conditions
# the coordinates are D.x1 D.x2 and D.x3 (r, theta, phi)

p_coords = np.array(D.x3)

i_phi_slice = 250 + 5*np.arange(18)

n_slice = len(i_phi_slice)

for ii in range(0, n_slice):
    if (i_phi_slice[ii] >= len(p_coords)):
        i_phi_slice[ii] = i_phi_slice[ii] - len(p_coords)

phi_slice = p_coords[i_phi_slice]

# End of phi slice calculation

# Plot color scheme

cmap = 'rainbow'

# Plot title
title = None

# Min and Max y-axis values

ymin = 0.0
ymax = 2500.00 

# X and Y axes labels

xlabel = 'R (AU)'
ylabel = 'V (km s$^{-1}$)'

# 
# 
# Create a figure
fig, axs = plt.subplots(1,1, figsize=(10, 6))

axs = pviz.plot_slice(D=D, var_name = var_name,
        r_cut=r_cut, theta_cut=theta_slice, phi_cut=phi_slice, ax=axs, cmap = cmap, 
        title = title, xlabel = xlabel, ylabel = ylabel, ymin = ymin, ymax = ymax)


plt.savefig(file_name, dpi=150)
print("Saving Plot to: ", file_name)

# to plot to screen uncomment the line below
#plt.show()


