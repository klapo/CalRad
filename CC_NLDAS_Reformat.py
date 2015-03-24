# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# Compile and write NLDAS SWdwn over the CA study area for Jan 2005 (example)
# netcdf
from netCDF4 import Dataset,num2date, date2num
import numpy as np
from datetime import datetime, timedelta

# OS interaction
import sys
import os

# Graphing
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap

# <codecell>

# Output file name
fname = '../../NLDAS_Jan_SWdwn.nc'

# Directory of .nc files
wd = '/Users/karllapo/Desktop/UW/DATA/NLDAS_Jan2005/test'
os.chdir(wd)
content = os.listdir(os.getcwd())
num_files = len([name for name in os.listdir('.') if os.path.isfile(name)])

# Initialize SWdwn
SWdwn = np.empty((50,50,num_files))
SWdwn[:] = np.nan
# Set iteration counter
iter = 0
# Initialize date variable
dates = []
# Bounding box
LL = [120.5,35]
LR = [115,35]
UR = [115,39]
UL = [122.5,39]

# <codecell>

# NLDAS is stored with the CONUS as a single time slice
# e.g., for Jan 21st, 2005 @ 15UTG: NLDAS_FORA0125_H.A20050121.1500.002.2014325204303.pss.nc
# Iterate over nc files and store in a single netcdf structure
for files in content:
    # Hokey catch for NLDAS files
    if files[-2:] == 'nc' and len(files) > 30:
        f = Dataset(files,'r')
        print(files)

        #Bounding box for CA study area
        lat = f.variables['lat_110'][:]
        lon = f.variables['lon_110'][:]

        # Extract only within cartesian box
        ind_lat = np.flatnonzero( (lat > LL[1]) & (lat < UL[1]) )
        ind_lon = np.flatnonzero( (-lon > UR[0]) & (-lon < LL[0]) )
        SWdwn_this_time_step = f.variables['DSWRF_110_SFC'][ind_lat,ind_lon]
        nlat = ind_lat.size 
        nlon = ind_lon.size
        lat = lat[ind_lat]
        lon = lon[ind_lon]
        
        # Fill this time slice to 3D array
        SWdwn[0:nlat,0:nlon,iter] = SWdwn_this_time_step
        
        # Time stamp
        YY = int(files[-20:-16])
        MM = int(files[-16:-14])
        DD = int(files[-14:-12])
        HH = int(files[-11:-9])
        dates.append(datetime(YY,MM,DD,HH,0,0))
        
        # Increment
        iter = iter + 1

# <codecell>

# Trim
SWdwn = SWdwn[0:nlat,0:nlon,0:iter]

# Open netcdf for writing, remove and re-write if previous file version
if os.path.isfile(fname):
    os.remove(fname)
ncfile = Dataset(fname,'w')

# create the lat, lon, time dimensions.
ncfile.createDimension('latitude',nlat)
ncfile.createDimension('longitude',nlon)
ncfile.createDimension('time',None)

# Define the coordinate variables (not needed for time)
lat_nc_out = ncfile.createVariable('latitude',np.float32,('latitude',))
lon_nc_out = ncfile.createVariable('longitude',np.float32,('longitude',))
t_out = ncfile.createVariable('time',np.float64, ('time',))

# Assign units attributes to coordinate vars
lat_nc_out.units = 'degrees_north'
lon_nc_out.units = 'degrees_east'
t_out.units = 'hours since 0001-01-01 00:00:00.0'
t_out.calendar = 'gregorian'

# write data to coordinate vars.
lat_nc_out[:] = lat
lon_nc_out[:] = lon
t_out[:] = date2num(dates,
            units=t_out.units,
            calendar=t_out.calendar)

# create output variable and units
SWdwn_nc_out = ncfile.createVariable('DSWRF_110_SFC',np.float64,('time','latitude','longitude'))
SWdwn_nc_out.units =  'W/m^2'

# write data to variables along time dimension.
for nt in range(iter):
    SWdwn_nc_out[nt,:,:] = SWdwn[:,:,nt]
    
ncfile.close()

# <codecell>

f = Dataset(fname,'r')
print(f.variables['latitude'][:])

# <codecell>


# <codecell>


