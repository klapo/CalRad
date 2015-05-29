import os
import sys
import numpy as np

####################
## Download Flags ##
####################
SURF_FLAG = 0		# Download surface met data
RAD_FLAG = 1		# Download surface rad data
dl_url = list()

if SURF_FLAG:
	dl_url.append('ftp://ftp.etl.noaa.gov/users/mhughes/Sierra_precip/WRFmicroandLBCoutput/surfoutd02_CONUSnl_WRF_ERAI_')
	dl_url.append('ftp://ftp.etl.noaa.gov/users/mhughes/Sierra_precip/WRFmicroandLBCoutput/surfoutd02_CONUSn1_WRF_NARR_')
	
if RAD_FLAG:
	dl_url.append('ftp://ftp.etl.noaa.gov/users/mhughes/Sierra_precip/WRFmicroandLBCoutput/radoutd02_CONUSnl_WRF_ERAI_')
	dl_url.append('ftp://ftp.etl.noaa.gov/users/mhughes/Sierra_precip/WRFmicroandLBCoutput/radoutd02_CONUSn1_WRF_NARR_')

microphysics_list = ['Morr', 'Thom', 'WSM6']
dates_year = [2007,2008]
dates_month = np.arange(1,13)

for url in dl_url:
	for microphysics in microphysics_list:
		for years in dates_year:
			for months in dates_month:
				if months < 10:
					filename = url+microphysics+"_"+str(years)+"-0"+str(months)+".nc"
				elif months >= 10:
					filename = url+microphysics+"_"+str(years)+"-"+str(months)+".nc"
				
				if months < 7 and years == 2007:
					continue

				# Grab the file!
				print(filename)
				os.system('curl -O '+filename)
