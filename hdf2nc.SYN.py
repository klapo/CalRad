# Compile and write NLDAS SWdwn over the CA study area for Jan 2005 (example)
# Memory management

# use this script to conver to netcdf, extract relevant info, close the netcdf file, delete it, ad nauseum, until I have a small manageable netcdf file
import gc

# netcdf
from netCDF4 import Dataset,num2date, date2num
import numpy as np
from datetime import datetime, timedelta

# OS interaction
import sys
import os

# Directories
hdfdir= '/home/dingo/2/ceres.syn/ed3a/orig.data'

# Change this directory name
ncdir = '/home/dingo/2/snowmelt/Lapo/CloudClimatology/SYN.netcdf.2002.2014'
years = np.arange(2002,2014)

# Output file names
fname = 'SYN.irrad.nc'
fname_CA_domain = 'CA.SYN.irrad.nc'

# Set initilization flag for reading variables
init_flag = 1
# Set iteration counter
iter = 0
# Initialize date variable
dates = []
# Bounding box - rectangular domain
LL_rect = [-123.5,34.5]
LR_rect = [-115.5,34.5]
UR_rect = [-115.5,41]
UL_rect = [-123.5,41]
# Bounding box - ragged domain
LL_rag = [-120,34.5]
LR_rag = [-115,34.5]
UR_rag = [-118.5,41]
UL_rag = [-123.5,41]

# Directory contents
os.chdir(hdfdir)
content = os.listdir(os.getcwd())
for files in content:
	if files[0:-8] == 'CER_SYN1deg-3Hour_Terra-Aqua-MODIS_Edition3A_':
		print(files) 
		os.system('ncl_convert2nc '+files+' -B -e hdf'+' -o '+ncdir)
		netcdf_file = files+".nc"

		if os.path.isfile(netcdf_file) 
		## Initialize - Grab first file, get details for pre-allocation
		print('Pre-allocating....   \n')
	
			# Need to somehow read the created netcdf file
			f = Dataset(netcdf_file,'r')

			## Bounding box for CA study area
			lat = f.variables['Colatitude__1371'][:]
			lon = f.variables['Longitude__1371'][:]

			# Extract only within rectangular domain
			ind_lat_rect = np.flatnonzero( (lat > LR_rect[1]) & (lat < UL_rect[1]) )
			ind_lon_rect = np.flatnonzero( (lon < LR_rect[0]) & (lon > LL_rect[0]) )
			
			# 2D array of lat-lon for operands below
			lat_m,lon_m = np.meshgrid(lat[ind_lat_rect],lon[ind_lon_rect])

			## Ragged domain, CA study area - indices for NaN'ing data outside 
			line_west_m = (UL_rag[1]-LL_rag[1])/(UL_rag[0]-LL_rag[0])
			line_west_b = LL_rag[1]-line_west_m*LL_rag[0]
			line_east_m = (UR_rag[1]-LR_rag[1])/(UR_rag[0]-LR_rag[0])
			line_east_b = LR_rag[1]-line_east_m*LR_rag[0]	
			ind_lon_rag,ind_lat_rag = np.nonzero( (lon_m < (lat_m-line_west_b)/line_west_m) | (lon_m > (lat_m-line_east_b)/line_east_m) \
							 | (lat_m < LR_rag[1]) | (lat_m > UL_rag[1]))
			break

	# Pre-allocate - global
	nlat = ind_lat_rect.size 
	nlon = ind_lon_rect.size
	num_files = len([name for name in os.listdir('.') if os.path.isfile(name)])
	SWdwn = np.empty((nlat,nlon,num_files))
	SWdwn[:,:,:] = np.nan
	LWdwn = np.empty_like(SWdwn)
	LWdwn[:,:,:] = np.nan

	# Pre-allocate - CA domain
	nlat = ind_lat_rect.size 
	nlon = ind_lon_rect.size
	num_files = len([name for name in os.listdir('.') if os.path.isfile(name)])
	SWdwn = np.empty((nlat,nlon,num_files))
	SWdwn[:,:,:] = np.nan
	LWdwn = np.empty_like(SWdwn)
	LWdwn[:,:,:] = np.nan


	## Iterate over nc files and store in a single netcdf structure
	for files in content:
		# NLDAS is stored with the CONUS as a single time slice with each year in a separate directory
		# e.g., for Jan 21st, 2005 @ 15UTG: 2005/NLDAS_FORA0125_H.A20050121.1500.002.2014325204303.pss.nc
		# Catch for NLDAS files
		if files[-2:] == 'nc' and len(files) > 30:
		
			f = Dataset(files,'r')
			print(files)

			# Fill this time slice to 3D array
			SWdwn_this_time_step = f.variables['DSWRF_110_SFC'][ind_lat_rect,ind_lon_rect]
			LWdwn_this_time_step = f.variables['DLWRF_110_SFC'][ind_lat_rect,ind_lon_rect]
			SWdwn[:,:,iter] = SWdwn_this_time_step
			LWdwn[:,:,iter] = LWdwn_this_time_step
			
			# Time stamp
			YY = int(files[-20:-16])
			MM = int(files[-16:-14])
			DD = int(files[-14:-12])
			HH = int(files[-11:-9])
			dates.append(datetime(YY,MM,DD,HH,0,0))
			
			# Increment
			iter = iter + 1
			collect = gc.collect()

## NaN data outside CA domain
SWdwn[ind_lat_rag,ind_lon_rag,:] = np.nan
LWdwn[ind_lat_rag,ind_lon_rag,:] = np.nan

## Sort by date
ind = dates.argsort()
dates = sorted(dates)
SWdwn = SWdwn[:,:,ind]
LWdwn = LWdwn[:,:,ind]

## Open netcdf for writing, remove and re-write if previous file version
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
lat_nc_out[:] = lat[ind_lat_rect]
lon_nc_out[:] = lon[ind_lon_rect]
t_out[:] = date2num(dates,
            units=t_out.units,
            calendar=t_out.calendar)

# create output variable and units
SWdwn_nc_out = ncfile.createVariable('DSWRF_110_SFC',np.float64,('time','latitude','longitude'))
SWdwn_nc_out.units =  'W/m^2'
LWdwn_nc_out = ncfile.createVariable('DLWRF_110_SFC',np.float64,('time','latitude','longitude'))
LWdwn_nc_out.units =  'W/m^2'

# write data to variables along time dimension.
for nt in range(iter):
    SWdwn_nc_out[nt,:,:] = SWdwn[:,:,nt]
    LWdwn_nc_out[nt,:,:] = LWdwn[:,:,nt]
    
ncfile.close()

