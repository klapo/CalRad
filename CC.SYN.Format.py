## Import statements
import gc

# netcdf
from netCDF4 import Dataset,num2date, date2num
import numpy as np
from datetime import datetime, timedelta
import pandas as pd

# OS interaction
import sys
import os

# OS interaction
import sys
import os

## Directory of .nc files
# Change these directory names
hdfdir= '/home/dingo/2/ceres.syn/ed3a/orig.data'
ncdir = '/home/dingo/2/snowmelt/Lapo/CloudClimatology/SYN.netcdf.2002.2014'

## netcdf parameters
# Input file name structure
hdf_name = 'CER_SYN1deg-3Hour_Terra-Aqua-MODIS_Edition3A_'
# Output file names
fname_CA_domain = 'CA.SYN.irrad.nc'

# Set initilization flag for reading variables
init_flag = 1
# Set iteration counter
iter = 0
# Initialize date variable
dates = []
# Bounding box - rectangular domain
LL_rect = [-124+360,34]
LR_rect = [-115+360,34]
UR_rect = [-115+360,41.5]
UL_rect = [-124+360,41.5]

# Larger sub-domain to test
# LL_rect = [-150+360,30]
# LR_rect = [-100+360,30]
# UR_rect = [-150+360,55]
# UL_rect = [-100+360,55]

# Bounding box - ragged domain
LL_rag = [-120.3+360,34.4]
LR_rag = [-114.7+360,34.4]
UR_rag = [-118.2+360,41.6]
UL_rag = [-123.8+360,41.6]

# Fields to read - arbitrarily I am using the 1386 downwelling version. 
# May need to fix this in the future
SYN_netcdf_fields = ['Total_Sky_SW_flux_Direct__1391', \
    'Total_Sky_SW_flux_Diffuse__1391', \
    'Untuned_Total_Sky_LW_Surface_Down__1386', \
    'Untuned_Total_Sky_SW_Surface_Down__1386', \
    'Surface_altitude_above_sea_level__1371']

# Directory contents
os.chdir(hdfdir)
content = os.listdir(os.getcwd())
for files in content:
	if files[0:-8] == hdf_name:

		## Initialize - Grab first file, get details for pre-allocation
		print('Pre-allocating....   \n')
		os.system('ncl_convert2nc '+files+' -B -e hdf'+' -o '+ncdir)
		netcdf_file = ncdir+"/"+files+".nc"

		if not os.path.isfile(netcdf_file):
			print("Skipping found file: " +files)
			continue

		# Need to somehow read the created netcdf file
		f = Dataset(netcdf_file,'r')

		## Rectangular bounding box for CA study area
		# SYN lat/lon dim: NOT UNIFORM, LOTS OF PROBLEMS
		lat = 90-f.variables['Colatitude__1371'][:]
		lon = f.variables['Longitude__1371'][:]

		# Extract only within rectangular domain - use this index variable for data extraction from netcdf
		ind_lat_rect,ind_lon_rect = np.nonzero( (lat > LR_rect[1]) & (lat < UL_rect[1]) & \
										   (lon < LR_rect[0]) & (lon > LL_rect[0]) )        
		# lat/lon as series
		lat_ser = lat[:,0]
		lon_ser = lon[125,:]

		## lat/lon for output, in mesh format
		ind_lat_s_rect = np.flatnonzero( (lat_ser > LR_rect[1]) & (lat_ser < UL_rect[1]) )
		ind_lon_s_rect = np.flatnonzero( (lon_ser < LR_rect[0]) & (lon_ser > LL_rect[0]) )
		lat_m,lon_m = np.meshgrid(lat_ser[ind_lat_s_rect],lon_ser[ind_lon_s_rect])

		## Ragged domain, CA study area - indices for NaN'ing data outside 
		line_west_m = (UL_rag[1]-LL_rag[1])/(UL_rag[0]-LL_rag[0])
		line_west_b = LL_rag[1]-line_west_m*LL_rag[0]
		line_east_m = (UR_rag[1]-LR_rag[1])/(UR_rag[0]-LR_rag[0])
		line_east_b = LR_rag[1]-line_east_m*LR_rag[0]
		ind_lon_rag,ind_lat_rag = np.nonzero((lon_m < (lat_m-line_west_b)/line_west_m) | \
										 (lon_m > (lat_m-line_east_b)/line_east_m) | \
										 (lat_m < LR_rag[1]) | (lat_m > UL_rag[1]))
		f.close()
		os.system('rm '+netcdf_file)
		collect = gc.collect()
		break
        
# Pre-allocate
os.chdir(hdfdir)
nlat = ind_lat_s_rect.size
nlon = ind_lon_s_rect.size
num_files = len([name for name in os.listdir('.') if name[0:-8] == hdf_name])
SWdwn_Direct = np.empty((num_files*8,nlat,nlon))
SWdwn_Diffuse = np.empty((num_files*8,nlat,nlon))
LWdwn = np.empty((num_files*8,nlat,nlon))
SWdwn = np.empty((num_files*8,nlat,nlon))
Elev = np.empty((num_files*8,nlat,nlon))

## Iterate over hdf files, convert to netcdf, read/extract, delete netcdf, and store read data in a single netcdf structure
for files in content:
	if files[0:-8] == hdf_name:
		# Convert and process each file
		print("Working on: "+files)
		os.system('ncl_convert2nc '+files+' -B -e hdf'+' -o '+ncdir)
		netcdf_file = ncdir+"/"+files+".nc"
		if not os.path.isfile(netcdf_file):
			print("Skipping found file: " +files)
			continue
		f = Dataset(netcdf_file,'r')

		# Fill this time slice to 3D array
		SWdwn_Direct_dt = f.variables['Total_Sky_SW_flux_Direct__1391'][:,ind_lat_s_rect,ind_lon_s_rect]
		SWdwn_Diffuse_dt = f.variables['Total_Sky_SW_flux_Diffuse__1391'][:,ind_lat_s_rect,ind_lon_s_rect]
		LWdwn_dt = f.variables['Untuned_Total_Sky_LW_Surface_Down__1386'][:,ind_lat_s_rect,ind_lon_s_rect]
		SWdwn_dt = f.variables['Untuned_Total_Sky_SW_Surface_Down__1386'][:,ind_lat_s_rect,ind_lon_s_rect]
		Elev_dt = f.variables['Surface_altitude_above_sea_level__1371'][:,ind_lat_s_rect,ind_lon_s_rect]

		dt_ind = np.arange(iter*8,iter*8+8)
		SWdwn_Direct[dt_ind,:,:] = SWdwn_Direct_dt
		SWdwn_Diffuse[dt_ind,:,:] = SWdwn_Diffuse_dt
		LWdwn[dt_ind,:,:] = LWdwn_dt
		SWdwn[dt_ind,:,:] = SWdwn_dt
		Elev[dt_ind,:,:] = Elev_dt

		# Time stamp
		YY = int(netcdf_file[-11:-7])
		MM = int(netcdf_file[-7:-5])
		DD = int(netcdf_file[-5:-3])
		HH = [1, 4, 7, 10, 13, 16, 19, 22] # UTC time stamp center points
		mm = int(30)
		[dates.append(datetime(YY,MM,DD,int(HOUR),mm,0)) for HOUR in HH]

		# Increment: release memory, delete netcdf
		iter = iter + 1
		f.close()
		os.system('rm '+netcdf_file)
		collect = gc.collect()

## NaN data outside CA domain
SWdwn[:,ind_lat_rag,ind_lon_rag] = np.nan
LWdwn[:,ind_lat_rag,ind_lon_rag] = np.nan
SWdwn_Diffuse[:,ind_lat_rag,ind_lon_rag] = np.nan
SWdwn_Direct[:,ind_lat_rag,ind_lon_rag] = np.nan
Elev[:,ind_lat_rag,ind_lon_rag] = np.nan

## Sort by date
ind = np.argsort(dates)
dates = sorted(dates)
SWdwn_Direct = SWdwn_Direct[ind,:,:]
SWdwn_Diffuse = SWdwn_Diffuse[ind,:,:]
LWdwn = LWdwn[ind,:,:]
SWdwn = SWdwn[ind,:,:]
Elev = Elev[ind,:,:]

## Open netcdf for writing, remove and re-write if previous file version
os.chdir(ncdir)
if os.path.isfile(fname_CA_domain):
    os.remove(fname_CA_domain)
ncfile = Dataset(fname_CA_domain,'w')

# create the lat, lon, time dimensions.
ncfile.createDimension('latitude',nlat)
ncfile.createDimension('longitude',nlon)
ncfile.createDimension('time',None)

# Define the coordinate variables (not needed for time)
lat_nc_out = ncfile.createVariable('latitude',np.float32,('latitude',))
lon_nc_out = ncfile.createVariable('longitude',np.float32,('longitude',))
t_out = ncfile.createVariable('time',np.float64, ('time',))

# Assign units attributes to coordinate vars
lat_nc_out.units = 'degrees'
lon_nc_out.units = 'degrees'
t_out.units = 'hours since 0001-01-01 00:00:00.0'
t_out.calendar = 'gregorian'

# write data to coordinate vars.
lat_nc_out[:] = lat_ser[ind_lat_s_rect]
lon_nc_out[:] = lon_ser[ind_lon_s_rect]
t_out[:] = date2num(dates,
            units=t_out.units,
            calendar=t_out.calendar)

# create output variable and units
SWdwn_nc_out = ncfile.createVariable('SWdwn',np.float64,('time','latitude','longitude'))
SWdwn_nc_out.units =  'W/m^2'
SWdwn_nc_out[:,:,:] = SWdwn[:,:,:]

LWdwn_nc_out = ncfile.createVariable('LWdwn',np.float64,('time','latitude','longitude'))
LWdwn_nc_out.units =  'W/m^2'
LWdwn_nc_out[:,:,:] = LWdwn[:,:,:]

SWdwn_Direct_nc_out = ncfile.createVariable('SWdwn_Direct',np.float64,('time','latitude','longitude'))
SWdwn_Direct_nc_out.units =  'W/m^2'
SWdwn_Direct_nc_out[:,:,:] = SWdwn_Direct[:,:,:]

SWdwn_Diffuse_nc_out = ncfile.createVariable('SWdwn_Diffuse',np.float64,('time','latitude','longitude'))
SWdwn_Diffuse_nc_out.units =  'W/m^2'
SWdwn_Diffuse_nc_out[:,:,:] = SWdwn_Diffuse[:,:,:]

Elev_nc_out = ncfile.createVariable('Elev',np.float64,('time','latitude','longitude'))
Elev_nc_out.units =  'm asl'
Elev_nc_out[:,:,:] = Elev[:,:,:]

ncfile.close()

