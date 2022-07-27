#
#  File:
#    Ion_pairs_regrid_RBSP-FB.py
#
#  Synopsis:
#    Read in netCDF files for 
#       Prod[time,pressure]...currently constant modeled as input 
#           to mo_solarproton.F90 following Jackman spedata.F90  
#            (vertical pressure levels are 0-400 spaced by 2 km)
#       and then regrid to waccm lev grid and time
#
#  Author:
#    Katharine Duderstadt 
#  
#  Date:
#    9 June, 2017
#    4 Dec 17 -- modify for RBSP-FB runs
#    18 Feb 2020 - modify for cesm2
#
#  Description:
#    Also use this routine to set the time interval over which the REP 
#     event occurs
#
#  Input:
#      Ion_pairs_RBSP-FU_Mar2013_L4_50per.nc
#      Ion_pairs_RBSP-FU_Mar2013_L4-5_50per.nc
#      Ion_pairs_RBSP-FU_Mar2013_L5_50per.nc
#      Ion_pairs_RBSP-FU_Mar2013_L5-5_50per.nc
#      jackman_pressure.csv   (to read in pressure levels)
#
#  Output:
#    cesm2_REP_ion_pairs_RBSP-FB.nc
#
#  Notes:
#    Need to then run "Add_dates_REP_ion_pairs.ncl" to get it in the 
#    right format for WACCM -- following jackman SPE input files
#
#    WACCM is in days since 1850-01-01 00:00:00  (note, used to be 1900)
#    Use Julian dates
#         Modify to begin on Feb 1, 2013 -- 3350+2452975 = 2456325
#         End May 31, 2013 -- 3469+2452975 = 2456444
#
#         "Event" begins Feb 14, 2013 -- 2456325
#          Event ends Mar 14, 2013  --  2456366
#
# Event begins Mar 4, 2013 -- 3381+2452975 = 2456355.5
# Event ends Mar 14, 2013 - 3391+2452975 = 2456366.5
#
# ADD A DAY
# Event begins Mar 4, 2013 -- 3381+2452975 = 2456355+1 = 2456356.5
# Event ends Mar 14, 2013 - 3391+2452975 = 2456366+1 = 2456376.5
#
# Choose a profile and create Prod array
#
# Set time index for WACCM to 48*(2453064-2452975) = 4272
#
# MODIFIED 4Dec17 to read in 11 days...and place Jan 1 - Jan 11
#      Julien Jan 1 2004 noon -- 2453006.0
#      Julien Jan 11 2004 noon -- 2453017.0
#
#   Right now we are calculating geomagnetic latitude from L-shell using simple dipole.
#      We may need to use IRBM or something more rigorous (or SSCweb with a POES satellite)
#
#        r = L cos^2(lambda)   where we assume r = 1 Earth radius
#
#            L = 3       glat = 55      (55.9)
#            L = 4       glat = 60      (58.75)
#            L = 4.5     glat = 62      (61.25)
#            L = 5       glat = 63      (63.75)
#            L = 6       glat = 66      (66.25)
#            L = 7       glat = 68      (68.5)
#
# --------------- import modules and packages -----------

import numpy, os
import Ngl
import Nio


from scipy import interpolate
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d

from math import exp, fabs

filename_in_L4 = "Ion_pairs_RBSP-FU_Mar2013_L5_50per.nc"
filename_in_L4_5 = "Ion_pairs_RBSP-FU_Mar2013_L4-5_50per.nc"
filename_in_L5 = "Ion_pairs_RBSP-FU_Mar2013_L5_50per.nc"
filename_in_L5_5 = "Ion_pairs_RBSP-FU_Mar2013_L5-5_50per.nc"

filename_pressure = "Jackman_Pressure.csv"
filename_out = "cesm2_REP_ion_pairs_RBSP-FB.nc"

# ---------------- Create Variables ----------------------
#  time_WACCM_index is the number of timesteps produced in output file
#  currently 30 min timesteps, so days*48
# --------------------------------------------------------

time_orig_index = 29
lev_orig_index = 201

#time_WACCM_index = 5712
time_WACCM_index = (2456459-2456277)*48
number_L = 4

lev_WACCM_index = 58
glat_WACCM_index = 32
chosen_date_index = 23

REP_prod_orig = numpy.zeros((time_orig_index, lev_orig_index, number_L),'f')
cesm2_REP_prod_orig = numpy.zeros((time_orig_index, lev_orig_index, glat_WACCM_index),'f')
REP_lev_orig = numpy.zeros(lev_orig_index,'f')
#REP_time_orig = numpy.zeros(time_orig_index, 'f') 
REP_julian_orig = numpy.zeros(time_orig_index,'float64') 
profile_chosen = numpy.zeros(lev_orig_index,'float64') 
REP_prod_origlev_time = numpy.zeros(lev_orig_index,'f')
 
REP_prod_newtime = numpy.zeros((time_WACCM_index, lev_orig_index, number_L),'float64')
REP_prod_WACCM = numpy.zeros((time_WACCM_index,lev_WACCM_index, number_L),'f')
cesm2_REP_prod_WACCM = numpy.zeros((time_WACCM_index,lev_WACCM_index, glat_WACCM_index),'f')
REP_lev_WACCM = numpy.zeros(lev_WACCM_index,'f')
REP_glat_WACCM = numpy.zeros(glat_WACCM_index,'f')
REP_time_WACCM = numpy.zeros(time_WACCM_index,'float64')
REP_julian_WACCM = numpy.zeros(time_WACCM_index,'float64')

# -------------------------------------------------------
# Read in netCDF file for REP ion pairs 
#
#      Ion_pairs_FIRE_FU-3_losscone.nc
# -------------------------------------------------------

nc_file_L4 = Nio.open_file(filename_in_L4,"r")

REP_file_ion_pair_L4 = nc_file_L4.variables["Prod"]

REP_file_julian = nc_file_L4.variables["julian"]
REP_file_pressure = nc_file_L4.variables["pressure"]

REP_prod_orig[:,:,0] = REP_file_ion_pair_L4[:,:]

REP_julian_orig[:] = REP_file_julian[:]
REP_lev_orig[:] = REP_file_pressure[:]

nc_file_L4.close()

# -----

nc_file_L4_5 = Nio.open_file(filename_in_L4_5,"r")

REP_file_ion_pair_L4_5 = nc_file_L4_5.variables["Prod"]
REP_prod_orig[:,:,1] = REP_file_ion_pair_L4_5[:,:]

nc_file_L4_5.close()

# -----

nc_file_L5 = Nio.open_file(filename_in_L5,"r")

REP_file_ion_pair_L5 = nc_file_L5.variables["Prod"]
REP_prod_orig[:,:,2] = REP_file_ion_pair_L5[:,:]

nc_file_L5.close()

# -----

nc_file_L5_5 = Nio.open_file(filename_in_L5_5,"r")

REP_file_ion_pair_L5_5 = nc_file_L5_5.variables["Prod"]
REP_prod_orig[:,:,3] = REP_file_ion_pair_L5_5[:,:]

nc_file_L5_5.close()


# --------------- Set Time Interval ----------------------
#
# ---------------------------------------------------------

#start_julian_int = 2456325
start_julian_int = 2456277
start_julian = numpy.float64(start_julian_int)
#end_julian_int = 2456444
end_julian_int = 2456459
end_julian = numpy.float64(end_julian_int)
start_julian_event = 2456338.5
end_julian_event = 2456366.5
event_duration = (end_julian_event - start_julian_event)*48
julian_event_index =int(48*(start_julian_event - start_julian_int))

interval = numpy.float64(1./48.) 
print(interval) 

#REP_julian_WACCM[:] = numpy.arange(start_julian,end_julian,interval)
REP_julian_WACCM = numpy.linspace(start_julian,end_julian,num=time_WACCM_index,endpoint=False,dtype='float64')


print('REP_julian_WACCM ',REP_julian_WACCM[:])

#REP_time_WACCM[:] = REP_julian_WACCM[:]-2415020.500000
REP_time_WACCM[:] = REP_julian_WACCM[:]-2415020.500000+18262.

ntime_event_start = 1
ntime_event_end = 1
#REP_prod_newtime[julian_event_index:julian_event_index_end,:] = REP_prod_orig[23,:] 
for nday in range (0,time_orig_index):
#   for ntimestep in range (0,48*time_orig_index):
   ntime_event_start = julian_event_index+48*(nday)
   ntime_event_end = ntime_event_start + 48
   REP_prod_newtime[ntime_event_start:ntime_event_end,:,:] = REP_prod_orig[nday,:,:]

#profile_chosen[:] = REP_prod_orig[chosen_date_index,:] 
#print('chosen date index', chosen_date_index)
#print('julian chosen date', REP_julian_WACCM[chosen_date_index])

# ---------------------------------------------------------
# Interpolate from Fang levels (0-400 km in 2 km layers)
# to WACCM pressure levels
# 
# Get pressure levels from jackman_pressure.csv file  
# --------------------------------------------------------

#data2 = Ngl.asciiread(filename_pressure,lev_WACCM_index,type="float",sep=",")
data2 = Ngl.asciiread(filename_pressure,(1,58),type="float",sep=",")

print(data2)
REP_lev_WACCM[:] = data2[:] 

#print('REP_lev_orig ',REP_lev_orig)
#print('REP_lev_WACCM ',REP_lev_WACCM)
REP_lev_WACCM[57] = 1.09e-8

for nLshell in range (0,number_L):

   for ntime in range (0,time_WACCM_index):
      REP_prod_origlev_time[:] = REP_prod_newtime[ntime,:,nLshell]
      REP_interp_lev =interpolate.interp1d(REP_lev_orig, REP_prod_origlev_time)
      REP_prod_WACCM[ntime,:,nLshell] = REP_interp_lev(REP_lev_WACCM)

# -------------------------------------------------------
# Set glat 
# -------------------------------------------------------

glat_mee = [-82.5, -72.9375, -70.1875, -68.5, -66.25, -63.75, -61.25, -58.75,\
    -55.9375, -52.1875, -47.5, -42.5, -37.5, -32.5, -25, -10, 10, 25, 32.5, \
    37.5, 42.5, 47.5, 52.1875, 55.9375, 58.75, 61.25, 63.75, 66.25, 68.5, \
    70.1875, 72.9375, 82.5]

absglat = 0.
for ntime in range (0,time_WACCM_index):
   for nlev in range (0,lev_WACCM_index):
      for nglat in range (0, glat_WACCM_index):
         absglat = fabs(glat_mee[nglat])
         if absglat >= 55. and absglat < 57.5:  #Lshells 3-4
            cesm2_REP_prod_WACCM[ntime,nlev,nglat] = REP_prod_WACCM[ntime,nlev,0]
         elif absglat >= 57.5 and absglat < 61.:  #Lshells 4
            cesm2_REP_prod_WACCM[ntime,nlev,nglat] = REP_prod_WACCM[ntime,nlev,0]
         elif absglat >= 61. and absglat < 62.:  #Lshells 4.5
            cesm2_REP_prod_WACCM[ntime,nlev,nglat] = REP_prod_WACCM[ntime,nlev,1]
         elif absglat >= 62. and absglat < 64.0:  #Lshells 5.
            cesm2_REP_prod_WACCM[ntime,nlev,nglat] = REP_prod_WACCM[ntime,nlev,2]
         elif absglat >= 64. and absglat < 69.0:  #Lshells 5.5
            cesm2_REP_prod_WACCM[ntime,nlev,nglat] = REP_prod_WACCM[ntime,nlev,3]
         else:
            cesm2_REP_prod_WACCM[ntime,nlev,nglat] = 0.
          

# -------------------------------------------------------
# Create netCDF file
# -------------------------------------------------------

ncfile2 = Nio.open_file(filename_out,'w') 

print('Write netCDF file to ',filename_out)

# create dimensions

ntime = time_WACCM_index 
nlevel = lev_WACCM_index
nglat = glat_WACCM_index 

#ncfile2.create_dimension('pressure',nlevel)
ncfile2.create_dimension('time',None)
ncfile2.create_dimension('plev',nlevel)
ncfile2.create_dimension('glat',nglat)

#    define variables
time = ncfile2.create_variable('time','d',('time',))
plev = ncfile2.create_variable('plev','f',('plev',))
glat = ncfile2.create_variable('glat','f',('glat',))
Prod_REP = ncfile2.create_variable('rep_ion','f',('time','plev','glat'))

ncfile2.variables['time'][:] = REP_time_WACCM 
ncfile2.variables['plev'][:] = REP_lev_WACCM 
ncfile2.variables['glat'][:] = glat_mee 
ncfile2.variables['rep_ion'][:,:,:] = cesm2_REP_prod_WACCM 

#   close ncfile
ncfile2.close()

Ngl.end()



