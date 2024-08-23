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

# name of variable that will be plotted
var_name = 'vx1'

local_path = "/Volumes/work1/michal/runs/cr2261-high/hmi_mast_mas_std_0201/"

pluto_path = local_path+myrun+'/output/'
# the 'vx1' need to be moved to a 'tmp' directory. Also copy to this directory the dbl.out and grid.out
# edit the first line of dble.out deleting ALL output variables except for 'vx1'
# This 'hack' is needed to expediate reading all the PLUTO output files for a single variable

pluto_path = local_path+myrun+'/tmp/'

# directory for saving plots
png_dir = local_path+myrun+'/png_files/'

# r_val holds the location of distance for the cuts
# for figure 15: r_val = 1.0 
r_val = 1.0 

# set the file name
file_name = png_dir +'figure11_'+myrun+'.png'

# Meridonial location for phi
p_val = np.deg2rad(290)

# the stack plot is for 1-D cuts in the theta_dimension
stack_dim = 't'

# read the dbl.out file and retrieve time information from it
wdir = pluto_path

time = io.read_time(w_dir=wdir,datatype='dbl')
time = np.array(time, dtype=np.float32)

# Time in hours - conversion factor for time from PLUTO units to hours
time_fac_pluto = 1.49597871e+08/3600 #time_unit_str = '(hours)'; %1AU/1km/s
time_h = time * time_fac_pluto

# set the zero of time to the end of the relaxation == CME initiation
time_h = time_h - time_h[1]

# retrieve variable
print('Retrieving variable. This may take some time..')

nlinf = io.nlast_info(w_dir=wdir, datatype='dbl')

nlast = nlinf['nlast']

pluto_list = []

for ii in range(0, nlast):
	D = pp.pload(ii,w_dir=wdir,datatype='dbl')
	pluto_list.append(D)


#####################################
# Calculate  t_val for the stack-plot

# the coordinates are D.x1 D.x2 and D.x3 (r, theta, phi)
r_coords = np.array(D.x1)
t_coords = np.array(D.x2)
p_coords = np.array(D.x3)

# location is the same as radial rays in figure11

# Convert from co-latitude to latitude

t_coords = np.pi / 2 - t_coords

i_th_slice = 85 + 5*np.arange(18)

th_slice_all = t_coords[i_th_slice]

n_slice = len(th_slice_all)

indices = np.linspace(0, 1, n_slice)

##
#####################################

# X and Y axes labels

xlabel = 'Time (hours since launch)'
ylabel = ''

# Plot title
title = None

# Plot color scheme
cmap = 'rainbow'

# value of y-shift for stack plot
yshift = 2000

# 
# Create a figure
fig, axs = plt.subplots(1,1, figsize=(10, 6))

axs = pviz.plot_stack(pluto_list=pluto_list, var_name = var_name,
        r_val = r_val, t_val=th_slice_all, p_val=p_val, stack_dim = stack_dim, time = time_h, ax=axs, cmap = cmap,
        title = title, xlabel = xlabel, ylabel = ylabel, log_scale = False, yshift = yshift)

plt.savefig(file_name, dpi=150)
print("Saving Plot to: ", file_name)

# to plot to screen uncomment the line below
# plt.show()

