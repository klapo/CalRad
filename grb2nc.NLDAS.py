import os

grbdir= '/home/dingo/2/snowmelt/NLDAS/data.2012/'
ncdir = '/home/dingo/2/snowmelt/NLDAS/data.2012/netcdf/'
years = np.arange(2004,2011)

# Directory contents
for dir_year in years:
	os.chdir(grbdir+dir_year)
	content = os.listdir(os.getcwd())
    for files in content:
		if files[-3:] == 'grb':
			print(files)
			os.system('ncl_convert2nc '+files+' -B'+' -o '+ncdir)

