# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

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
# datadir = '/home/dingo/2/snowmelt/MODIS/gridded.08222013/orig.data'
datadir = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/MODIS.IRRAD'
# outdir = '/home/dingo/2/snowmelt/Lapo/'
outdir = '/Users/karllapo/gdrive/SnowHydrology/proj/CloudClimatology/data/MODIS.IRRAD'

## OUTPUT
fname = 'CA.MODIS.irrad.nc'

# <codecell>

## Functions/classes for reading binary files in python
## From: http://code.activestate.com/recipes/577610-decoding-binary-files/
import struct

class BinaryReaderEOFException(Exception):
    def __init__(self):
        pass
    def __str__(self):
        return 'Not enough bytes in file to satisfy read request'

class BinaryReader:
    # Map well-known type names into struct format characters.
    typeNames = {
        'int8'   :'b',
        'uint8'  :'B',
        'int16'  :'h',
        'uint16' :'H',
        'int32'  :'i',
        'uint32' :'I',
        'int64'  :'q',
        'uint64' :'Q',
        'float'  :'f',
        'double' :'d',
        'char'   :'s'}

    def __init__(self, fileName):
        self.file = open(fileName, 'rb')
        
    def read(self, typeName):
        typeFormat = BinaryReader.typeNames[typeName.lower()]
        typeSize = struct.calcsize(typeFormat)
        value = self.file.read(typeSize)
        if typeSize != len(value):
            raise BinaryReaderEOFException
        return struct.unpack(typeFormat, value)[0]
    
    def __del__(self):
        self.file.close()

# <codecell>

## Data directory
os.chdir(datadir)
content = os.listdir(os.getcwd())

## Initiliaze empty numpy object
init = 1

## Initialize - Grab first file, get details for pre-allocation
for files in content:

    # See README file for details - MODIS files names: ssda_grid_0.05deg_ins_(DATES)
    if files[-3:] == 'dat' and len(files) > 20:
        print('Pre-allocating....')
        print('Found file: '+files)
        
        binaryReader = BinaryReader(files)
        try:
            lon1 = binaryReader.read('float')
            lon2 = binaryReader.read('float')
            lat1 = binaryReader.read('float')
            lat2 = binaryReader.read('float')
            resolution = binaryReader.read('float')
            nlat = binaryReader.read('int32')
            nlon = binaryReader.read('int32')
            ntimes = binaryReader.read('int32')
            misv = binaryReader.read('int32')
            
            num_dates = ntimes
            Mdates = np.empty((num_dates,6))
            while num_dates > 0:
                date_ind = ntimes-num_dates
                for k in np.arange(0,6):
                    Mdates[date_ind,k] = binaryReader.read('int32')
                num_dates = num_dates-1
                
            number_binary_items = ntimes*nlat*nlon
            num_bin = number_binary_items
            SW = np.empty(number_binary_items)
            while num_bin > 0:
                binary_ind = number_binary_items-num_bin
                SW[binary_ind] = binaryReader.read('float')
                num_bin = num_bin - 1
                
        except BinaryReaderEOFException:
            # One of our attempts to read a field went beyond the end of the file.
            print "Error: File seems to be corrupted."

        # Re-shape into an array of space and time
        SW = SW.reshape(nlat,nlon,ntimes,order='F')
    
        if init == 1:
            modis_sw = SW
            MODISdates = []
            [MODISdates.append(datetime(int(d[0]),int(d[1]),int(d[2]),int(d[3]),int(d[4]),int(d[5]))) for d in Mdates]
            init = 0
        else:
            modis_sw = np.concatenate((modis_sw,SW),axis=2)
            [MODISdates.append(datetime(int(d[0]),int(d[1]),int(d[2]),int(d[3]),int(d[4]),int(d[5]))) for d in Mdates]
            

# <codecell>

# ###### Graphing ######
# ## Check that I'm extracting the appropriate location
# from matplotlib.pyplot import subplots
# import matplotlib.pyplot as plt
# from matplotlib import cm
# from mpl_toolkits.basemap import Basemap, addcyclic

# # must insert this statement to render the plots within the notebook
# # this is specific to the ipython notebook
# %matplotlib inline

# fig, axs = subplots(1,1)
# fig.set_size_inches(10,5)
# SW_reshape = SW.reshape(nlat,nlon,ntimes,order='F')
# # axs.imshow(SW_reshape[0:120,31:190,62], origin="lower",vmin=50,vmax=1000)
# # axs.imshow(modis_sw[ind_lat_rect.min():ind_lat_rect.max()+1,
#     #\ind_lon_rect.min():ind_lon_rect.max()+1,62], origin="lower",vmin=50,vmax=1000)
# # axs.imshow(modis_sw[:,:,62], origin="lower",vmin=50,vmax=1000)
# axs.imshow(modis_sw_out[:,:,62], origin="lower",vmin=50,vmax=1000)

# <codecell>

## Study domain indices
# NOTE: MODIS domain only goes to 35N -- the rectangle etc have been adjusted to account for this.
# Bounding box - rectangular domain
LL_rect = [-123.5,35.1]
LR_rect = [-115.5,35.1]
UR_rect = [-115.5,41]
UL_rect = [-123.5,41]
# Bounding box - ragged domain
LL_rag = [-120,34.5]
LR_rag = [-115,34.5]
UR_rag = [-118.5,41]
UL_rag = [-123.5,41]

# MODIS grid -- add .05 to accomodate python indexing slicing off the last entry
lat = np.arange(lat1,lat2+.05,.05)
lon = np.arange(lon1,lon2+.05,.05)

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
# For now I'm assuming this worked...
ind_lon_rag,ind_lat_rag = np.nonzero( (lon_m < (lat_m-line_west_b)/line_west_m) | \
                (lon_m > (lat_m-line_east_b)/line_east_m) | (lat_m < LR_rag[1]) | (lat_m > UL_rag[1]))


# Bizarre behavior here: python won't broadcast the modis_sw to modis_sw_out using the indices:
#    - ind_lat_rect and ind_lon_rect
#    - instead construct the indices as ind_lat_rect.min():ind_lat_rect.max(),ind_lon_rect.min():ind_lon_rect.max()

# <codecell>

## NaN data outside CA domain
# modis data w/in study region (rectangle)
# See above note about indexing
modis_sw_out = modis_sw[ind_lat_rect.min():ind_lat_rect.max()+1,ind_lon_rect.min():ind_lon_rect.max()+1,:]
# modis data w/in study region (ragged domain)
modis_sw_out[ind_lat_rag,ind_lon_rag,:] = np.nan

## No sorting by date--MODIS data should already be sorted

## Output variables/parameters
lat_out = lat[ind_lat_rect]
lon_out = lon[ind_lon_rect]
nlat = lat_out.shape[0]
nlon = lon_out.shape[0]
ndate = np.shape(MODISdates)[0]

## Open netcdf for writing, remove and re-write if previous file version
os.chdir(outdir)
if os.path.isfile(fname):
    os.remove(fname)
ncfile = Dataset(fname,'w')

# create the lat, lon, time dimensions.
ncfile.createDimension('latitude',nlat)
ncfile.createDimension('longitude',nlon)
ncfile.createDimension('time',ndate)

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
lat_nc_out[:] = lat[ind_lat_rect]
lon_nc_out[:] = lon[ind_lon_rect]
t_out[:] = date2num(dates,
            units=t_out.units,
            calendar=t_out.calendar)

## create output variable and units
SWdwn_nc_out = ncfile.createVariable('SWdwn',np.float64,('time','latitude','longitude'))
SWdwn_nc_out.units =  'W/m^2'

## write data to variables along time dimension.
for nt in range(count):
    SWdwn_nc_out[nt,:,:] = modis_sw_out[:,:,nt]

ncfile.close()

