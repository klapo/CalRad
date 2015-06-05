## Import statements
# netcdf
from netCDF4 import Dataset,num2date, date2num
import numpy as np
from datetime import datetime, timedelta
import pandas as pd

# OS interaction
import sys
import os
import subprocess

## Directory of .nc files
datadir = '/home/dingo/2/snowmelt/Lapo/CloudClimatology/VIC_MTCLIM/VIC.Disagg.ForcingFluxes'
#datadir = '/home/lapok/Livneh_VIC/VIC.Disagg.ForcingFluxes/'
#datadir = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/VIC_MTCLIM'
outdir = '/home/dingo/2/snowmelt/Lapo/CloudClimatology/VIC_MTCLIM/'
#outdir = '/home/lapok/Livneh_VIC/'
#outdir = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/VIC_MTCLIM'

## File names
fname = 'CA.MTCLIM.irrad.nc'
fin_name = 'mtclim_irrad_'

## netcdf parameters
# Set initilization flag for reading variables
init_flag = 1
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

## Directory contents
os.chdir(datadir)
p = subprocess.Popen('ls -l *.gz | wc -l', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
for line in p.stdout.readlines():
    num_files = int(line)
content = os.listdir(os.getcwd())

## Dimensions
ndates = 105192 # Found using line counts from vim
lat_fname = []
lon_fname = []

print("Preparing dimensions")
for files in content:
    if os.path.isfile(files) and files[0:13] == fin_name:        
        # Strip out lat and lon from name
        lat_fname.append(float(files[13:21]))
        lon_fname.append(float(files[22:32]))
 
## Create spatial variables
lat = np.unique(lat_fname)
lon = np.unique(lon_fname)
nlat = lat.shape[0]
nlon = lon.shape[0]
lat_m,lon_m = np.meshgrid(lat,lon)

## Create time variable
t_beg = datetime(2000,1,1,0)
t_end = datetime(2011,12,31,23)
dates = pd.date_range(start=t_beg,end=t_end,freq='1H')

## Pre-allocation
SWdwn = np.empty((ndates,nlat,nlon))
LWdwn = np.empty_like(SWdwn)

## Read each file
content = os.listdir(os.getcwd())
for files in content:
    if os.path.isfile(files) and files[0:13] == fin_name:
		print("Processing: "+files)
	
		# Unzip the file - use Popen to wait until the action executes
		p = subprocess.Popen('gunzip '+files, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p.wait()
		#		os.system('gunzip '+files)
		files_unzip = files[0:-3]

		# Read the text file (2 columns, each the time vector in length, for each spatial location
		f = pd.read_csv(files_unzip,delimiter='\t',names=['SWdwn','LWdwn'])

		# Strip out lat and lon from name
		lat_fname = float(files_unzip[13:21])
		lon_fname = float(files_unzip[22:32])

		# Find where the files belongs in the numpy arrays and assign
		ind_lon,ind_lat = np.nonzero( (lat_fname == lat_m) & (lon_fname == lon_m) )
		LWdwn[:,ind_lat[0],ind_lon[0]] = f.LWdwn[:]
		SWdwn[:,ind_lat[0],ind_lon[0]] = f.SWdwn[:]

		# Re zip the file
		#        os.system('gzip '+files_unzip)
		p = subprocess.Popen('gxip '+files_unzip, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		p.wait()

# Extract only within rectangular domain - use this index variable for data extraction from netcdf
ind_lat_s_rect = np.flatnonzero( (lat > LR_rect[1]) & (lat < UL_rect[1]) )
ind_lon_s_rect = np.flatnonzero( (lon < LR_rect[0]) & (lon > LL_rect[0]) )

# 2D array of lat-lon for operands below
lat_m,lon_m = np.meshgrid(lat[ind_lat_s_rect],lon[ind_lon_s_rect])

# Indices for extracting from numpy array
ind_lat_rect, ind_lon_rect = np.nonzero( (lat_m > LR_rect[1]) & (lat_m < UL_rect[1]) \
                                            & (lon_m < LR_rect[0]) & (lon_m > LL_rect[0]) )

## Ragged domain, CA study area - indices for NaN'ing data outside 
line_west_m = (UL_rag[1]-LL_rag[1])/(UL_rag[0]-LL_rag[0])
line_west_b = LL_rag[1]-line_west_m*LL_rag[0]
line_east_m = (UR_rag[1]-LR_rag[1])/(UR_rag[0]-LR_rag[0])
line_east_b = LR_rag[1]-line_east_m*LR_rag[0]
#ind = np.flatnonzero( (lon_m > (lat_m-line_west_b)/line_west_m) & (lon_m < (lat_m-line_east_b)/line_east_m) \
#                         & (lat_m > LR[1]) & (lat_m < UL[1]))
ind_lon_rag,ind_lat_rag = np.nonzero( (lon_m < (lat_m-line_west_b)/line_west_m) | (lon_m > (lat_m-line_east_b)/line_east_m) \
                 | (lat_m < LR_rag[1]) | (lat_m > UL_rag[1]))


## CA domain
# Let's hack python cause it is totally mental
# Trim to rectangular domain
SWdwn_out = SWdwn[:,ind_lat_s_rect.min():ind_lat_s_rect.max()+1,ind_lon_s_rect.min():ind_lon_s_rect.max()+1]
# NaN data outside 
SWdwn_out[:,ind_lat_rag,ind_lon_rag] = np.nan

## No sorting by date--output is already sorted

## Output variables/parameters
lat_out = lat[ind_lat_s_rect]
lon_out = lon[ind_lon_s_rect]
nlat = lat_out.shape[0]
nlon = lon_out.shape[0]

## Open netcdf for writing, remove and re-write if previous file version
os.chdir(outdir)
if os.path.isfile(fname):
    os.remove(fname)
ncfile = Dataset(fname,'w')

## create the lat, lon, time dimensions.
ncfile.createDimension('latitude',nlat)
ncfile.createDimension('longitude',nlon)
ncfile.createDimension('time',ndates)

## Define the coordinate variables
lat_nc_out = ncfile.createVariable('latitude',np.float32,('latitude',))
lon_nc_out = ncfile.createVariable('longitude',np.float32,('longitude',))
t_out = ncfile.createVariable('time',np.float64, ('time',))

## Assign units attributes to coordinate vars
lat_nc_out.units = 'degrees_north'
lon_nc_out.units = 'degrees_east'
t_out.units = 'hours since 0001-01-01 00:00:00.0'
t_out.calendar = 'gregorian'

## write data to coordinate vars.
lat_nc_out[:] = lat_out
lon_nc_out[:] = lon_out
t_out[:] = date2num(dates.to_pydatetime(),
            units=t_out.units,
            calendar=t_out.calendar)

## create output variable and units
SWdwn_nc_out = ncfile.createVariable('SWdwn',np.float64,('time','latitude','longitude'))
SWdwn_nc_out.units =  'W/m^2'
SWdwn_nc_out = SWdwn

LWdwn_nc_out = ncfile.createVariable('LWdwn',np.float64,('time','latitude','longitude'))
LWdwn_nc_out.units =  'W/m^2'
LWdwn_nc_out = LWdwn

ncfile.close()

