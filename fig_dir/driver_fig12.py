import os
import numpy as np
import matplotlib.pyplot as plt

import pysunrunner
import pysunrunner.pload as pp
import pysunrunner.io as io

from scipy import interpolate

import pyspedas
from pyspedas import time_string
from pyspedas import time_double

import pandas as pd

###############################################################################
# Load SOLO data

# change this path to match your location of SOLO cdv data file
csv_file_path = 'results/orbiter.csv'
orbiter = pd.read_csv(csv_file_path)
orbiter.index = pd.to_datetime(orbiter['Time'])

###############################################################################

# Define the run identifier and related parameters
myrun = 'run09f-2'

cme_rad = '44'
my_var = 'vx1'

# Path to PLUTO results 
local_path = "/Volumes/work1/michal/runs/cr2261-high/hmi_mast_mas_std_0201/"

pluto_path = local_path+myrun+'/output/'

# Directory to save plot
png_dir=local_path+'/'+myrun+'/png_files'

wdir = pluto_path

# Retrieve the total number of PLUTO dumps
nlinf = io.nlast_info(w_dir=wdir, datatype='dbl')

# read the dbl.out file and retrieve time information from it
time = io.read_time(w_dir=wdir,datatype='dbl')
time = np.array(time, dtype=np.float32)

# Remove the first (zero) element corresponding to the start of relaxation
time = time[1:]

nlast = nlinf['nlast']
nlast1 = nlast+1

# Load PLUTO results starting from frame 1 (end of relaxation, CME initiation)
pluto_var_list = []

for ii in range(1, nlast1):
	D = pp.pload(ii,w_dir=wdir,datatype='dbl')
	pluto_var = getattr(D, my_var)
	pluto_var_list.append(pluto_var)

# Number of time steps
n_timesteps = len(pluto_var_list)

# Extract the grid in r, theta, and phi coordinates
r_coords = np.array(D.x1)
t_coords = np.array(D.x2)
p_coords = np.array(D.x3)

# Convert from co-latitude to latitude
t_coords = np.pi/2 - t_coords

# Conversion factor from astronomical units (AU) to solar radii (Rs)
r_fac_pluto=2.149e+02

# Conversion factor from PLUTO time units to hours
time_fac_pluto = 1.49597871e+08/3600 
time_h = time * time_fac_pluto

# Convert time to seconds, with time_s[1] being zero (eruption time)
# We use time_h[1] instead of time_h[0] because the latter corresponds to equilibration

time_s = (time_h - time_h[0]) * 60 * 60 

# Define the time of the CME (in UT)
cme_time = '2022-09-05/13:10:00'
cme_time_int = time_double(cme_time)

# Add the CME time to time_s to get the date/time for PLUTO output

time_pluto_int = np.zeros(len(time_s))
for ii in range(n_timesteps):
	time_pluto_int[ii] = time_s[ii] + cme_time_int

# Convert PLUTO times to string format
time_pluto_string = time_string(time_pluto_int,fmt='%Y-%m-%d/%H:%M:%S')



# Convert to datetime objects for comparison with SOLO data
pluto_times = pyspedas.time_datetime(time=time_pluto_string, tz='UTC')


# find nearest - we will update this soon to use the OMNI location 

# =====================================================
# define find_nearest function
# =====================================================
def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx, array[idx]

# =====================================================
# Perform Time Shift Based on Solar Rotation
# =====================================================

# Solar rotation in pluto time units

omega_corotate = 428.7616

twopi = 2.0 * np.pi

# Initialize necessary variables
ntime = len(time_h)
t0 = time[1]
t0_h = t0 * time_fac_pluto
n_phi     = len(p_coords)
n_r       = len(r_coords)
n_theta   = len(t_coords)

print("Doing Phi Shift")

for ii in range(n_timesteps):

	phi_shift = omega_corotate * (time[ii] - t0)
	pluto_var = pluto_var_list[ii]
	tmp = pluto_var * 0.0
	
	for k in range(n_phi):
		pv = p_coords[k] - phi_shift
		pv = pv % twopi
		if (pv < 0.0):
			pv = pv + twopi
	    
		ip1 = np.argmax(p_coords > pv)
		ip = ip1 - 1
		if (ip < 0):
			ip = ip1
	    
		if (p_coords[ip] == p_coords[ip1]):
			ap = 0.0
		else:
			ap = (pv - p_coords[ip])/(p_coords[ip1] - p_coords[ip])
		tmp[:,:,k] = (1.0-ap) * pluto_var[:, :, ip] + ap * pluto_var[:, :, ip1]
	
	# Update the shifted data 
	pluto_var_list[ii] = tmp


# =====================================================
# Time-Series Analysis with SOLO Data
# =====================================================

# Convert SOLO timestamps to datetime and ensure proper format
ts_orbiter = pd.DatetimeIndex(orbiter['Time'])
ts_model = pd.DatetimeIndex(pluto_times)

# Remove timezone information to align time zones
ts_model = ts_model.tz_localize(None)
ts_orbiter = ts_orbiter.tz_localize(None)


# Function to find the closest index in ts_orbiter for each timestamp in ts_model
def find_closest_indices(ts_model, ts_orbiter):
    indices = []
    for timestamp in ts_model:
        # Find the absolute difference between the current ts_model timestamp and all ts_orbiter timestamps
        differences = np.abs(ts_orbiter - timestamp)
        # Find the index of the minimum difference
        closest_index = differences.argmin()
        indices.append(closest_index)
    return indices

# Get the indices of ts_orbiter that match ts_model
matching_indices = find_closest_indices(ts_model, ts_orbiter)

# Function to create an index array for latitude and longitude sampling
def create_index_array(nlat, idx):
    if nlat % 2 == 0:
        raise ValueError("nlat must be odd.")
    
    half_nlat = (nlat - 1) // 2
    start_idx = idx - half_nlat
    end_idx = idx + half_nlat + 1  # +1 because the end index is exclusive in numpy.arange
    
    return np.arange(start_idx, end_idx)

# =====================================================
# Create Time-Series Plot
# =====================================================

# Meshgrid for latitude and longitude
y, x = np.meshgrid(np.rad2deg(t_coords), np.rad2deg(p_coords))

# Define minimum and maximum values for the color scale
zmin, zmax = 200, 2200

# Define CME longitude and latitude (in degrees)
cme_lon = 295.0
cme_lat = -40.0

# Number of latitude and longitude samples (grid resolution ~1 degree)
nlat = 7
nlon = 7

# Initialize arrays for storing model data
icount = 0

model_data = np.zeros((nlat, nlon, nlast))
model_data_1d = np.zeros(nlast)

nstart = 0

# Loop through each time step to extract data for the time-series plot
for ii in range(nstart, nlast):
	
	pluto_var = pluto_var_list[ii]

	index = matching_indices[ii]

	# Get orbiter longitude, latitude, and radial distance
	orbiter_lon = orbiter.iloc[index, orbiter.columns.get_loc('lon')]
	orbiter_lat = orbiter.iloc[index, orbiter.columns.get_loc('lat')]	
	orbiter_r = orbiter.iloc[index, orbiter.columns.get_loc('rAU')]

	# Find the closest grid points in the PLUTO data
	r_idx = np.argmin(np.abs(r_coords - orbiter_r))
	lon_idx = np.argmin(np.abs(p_coords - np.deg2rad(orbiter_lon)))
	lat_idx = np.argmin(np.abs(t_coords - np.deg2rad(orbiter_lat)))

	# Sample a range of latitudes and longitudes around the selected indices
	lat_idx_arr = create_index_array(nlat, lat_idx)
	lon_idx_arr = create_index_array(nlon, lon_idx)

	# Store the sampled data in the model_data array
	model_data[:,:,ii] = pluto_var[r_idx,lat_idx_arr,lon_idx_arr]
	model_data_1d[ii] = pluto_var[r_idx,lat_idx,lon_idx]

## Reshape the model data to 2D for plotting
model_data_2d = model_data.reshape(nlat*nlon, nlast)

# Create a DataFrame for the model data (transposed for proper orientation)

model_data_2d = model_data_2d.T
model_df = pd.DataFrame(model_data_2d, index = pd.to_datetime(pluto_times))

# Create a DataFrame for the 1D model data (representing a single point over time)
model_df_1d = pd.DataFrame({'values': model_data_1d}, index=pd.to_datetime(pluto_times))

# Define the SOLO CME time and add UTC timezone information
solo_cme_time = pd.to_datetime('2022-09-06/10:00:00')
solo_cme_time = solo_cme_time.tz_localize('UTC')

# Find the time in the model where the maximum velocity occurs
model_max_time = model_df_1d['values'].idxmax()

# Sort the model values in descending order to find the second highest peak
sorted_values = model_df_1d['values'].sort_values(ascending=False)

# Get the timestamp of the second highest value in the model
model_second_max_time = sorted_values.index[1]

# Calculate the time shift needed to align the model with the observed SOLO CME time
if (myrun == 'run09f'):
	shift_periods = (solo_cme_time - model_second_max_time)
else:
	shift_periods = (solo_cme_time - model_max_time)

# Apply the calculated time shift to the model DataFrame
shifted_model_df_1d = model_df_1d.copy()
shifted_model_df_1d.index = shifted_model_df_1d.index + shift_periods

# =====================================================
# Display Results
# =====================================================

# Calculate the minimum and maximum values for each timestep in the model

model_mnmx = np.zeros((nlast, 2))
for ii in range(0, nlast):
	model_mnmx[ii, 0] = np.min(model_data_2d[ii,:])
	model_mnmx[ii, 1] = np.max(model_data_2d[ii,:])

# Create a DataFrame for the minimum and maximum values over time
df = pd.DataFrame(model_mnmx, index = pd.to_datetime(pluto_times), columns=["model_min", "model_max"])

# Apply the same time shift to the min/max DataFrame
shifted_df = df.copy()
shifted_df.index = shifted_df.index + shift_periods

# =====================================================
# Plot the results
# =====================================================

# Create a figure and axis for the plot
fig, ax = plt.subplots(figsize=(10, 5), gridspec_kw={'hspace':0})

# Plot the observed SOLO data
ax.plot(orbiter.index, orbiter['vr'], color = 'orange', zorder = 1)
# Plot the shifted model data
ax.plot(shifted_model_df_1d.index, shifted_model_df_1d['values'], color = 'black', zorder = 4)
# Add a vertical line indicating the SOLO CME time
ax.axvline(pd.Timestamp('09/06/2022 10:00:00'), linestyle='dotted', color='grey', zorder = 3)
# Fill between the min and max model values over time
ax.fill_between(shifted_df.index, shifted_df["model_min"], df["model_max"], color='cornflowerblue', alpha=0.4, zorder = 2)
# Set the labels for the x and y axes
ax.set_ylabel(r'$\rm V \; [km  s^{-1}]$')
ax.set_xlabel("Date")
# Set the y-axis limits based on the min/max model values
ax.set_ylim([np.min(df.values), np.max(df.values)])
# Save the plot to a file
save_name = 'figure12_'+myrun+'.png'
file_name=png_dir+'/'+save_name
plt.savefig(file_name, dpi=150)
print("Saving Plot to: ", file_name)

# Uncomment the line below to display the plot on the screen
#plt.show()









