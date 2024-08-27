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

# Panels (a), (b) and (c) of Figure 10 are plotted at three different radial distances
# The first is in units of solar-radii and the last two are in AU
# PLUTO works in AU 

R_sun_to_au =  0.004650467260962157

# r_val holds the location of distance and latitude for the stack plots
# for panel (a): r_val = 31.0 * R_sun_to_au 
# for panel (b): r_val = 0.7
# for panel (c): r_val = 1.0 

## set r_val and panel name 
r_val = 1.0 * R_sun_to_au 

panel_name = 'a' #'b' 'c'

# set figure filename
file_name = png_dir +'figure10'+panel_name+'_'+myrun+'.png'


# For latitude: equatorial location
t_val = 0.0 

# the stack plot is for 1-D cuts in the phi_dimension
stack_dim = 'p'

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
# Calculate  p_val for the stack-plot

# the coordinates are D.x1 D.x2 and D.x3 (r, theta, phi)
r_coords = np.array(D.x1)
t_coords = np.array(D.x2)
p_coords = np.array(D.x3)

# location is the same as radial rays in figure07 
# including periodic boundary conditions

i_phi_slice = 250 + 5*np.arange(18)

n_slice = len(i_phi_slice)
for ii in range(0, n_slice):
	if (i_phi_slice[ii] >= len(p_coords)):
		i_phi_slice[ii] = i_phi_slice[ii] - len(p_coords)

fi_slice_all = p_coords[i_phi_slice]

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
yshift = 1200

# 
# Create a figure
fig, axs = plt.subplots(1,1, figsize=(10, 6))

axs = pviz.plot_stack(pluto_list=pluto_list, var_name = var_name,
        r_val = r_val, t_val=t_val, p_val=fi_slice_all, stack_dim = stack_dim, time = time_h, ax=axs, cmap = cmap,
        title = title, xlabel = xlabel, ylabel = ylabel, log_scale = False, yshift = yshift)

plt.savefig(file_name, dpi=150)
print("Saving Plot to: ", file_name)

# to plot to screen uncomment the line below
#plt.show()

