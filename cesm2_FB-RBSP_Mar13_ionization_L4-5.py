#
#  File:
#    cesm2_FB-RBSP_Mar13_ionization.py
#      based on Fang_parameterization_extend.py
#    Modified to read in RBSP data
#      and then scale according to FB conjunctions
#      fit/extrapolate using chosen method
#      calculate ionization rates
#      and save formatted for WACCM input
#
#   This version extends to 400 km in order to use in WACCM-X

#  Synopsis:
#     Calculates ion pair production rate profiles for electrons
#     with energies 100eV - 1MeV using Xiaohua Fang's parameterization
#
#     < 50 keV  multi-stream model of Lummerzheim et al., 1989 
#     and Lummerzheim and Lilenstein, 1994 
#
#     > 50 keV      two stream model (Solomon et al, 1988 and 
#     Solomon and Abreau, 1989) multiplied by a scaling factor 
#     (1.06 to 1.16 depending on energy) so that mean energy loss
#     per ion pair is 35 eV in accordance with laboratory measurements 
#     (Rees., 1989]
# 
#     Transform first-principle model results into normalized quantities 
#     and then fit a curve following the procedures initially used by 
#     Roble and Ridley, 1987
#
#  Author:
#    Katharine Duderstadt 
#  
#  Initial draft:
#    17 Feb, 2016

#  Most Recent modification:
#    18 Feb 2020
#
#  Input Files:
#
#   RBSPa_FEDU_2013Feb15-Mar15_5L.csv
#       (values from Chia-Lin 17 Feb 2020 -- RBSPa_FEDU_2013Feb15-Mar15.txt)
#       reformatted with julian date and separated into four L shell files 
#
#   FB-RBSP_ratios.csv (currently 1 MLT and 1 Lshell conjunctions, 50, 75, 100 percentile)
#    (need to get newest ratios_for .5 MLT and .5 Lshell)
#
#   MSIS90_1Jan04.csv  (change for correct month and year)
#
#   Fang_Table1.csv
#
#
#  Output:
#
#      Ion_pairs_RBSP-FU_Mar2013_3-7L_50per.nc
#
#      Fang_Fig1_ext.ps
#      Fang_Fig2.ps
#      Diff_flux_FIRE_FB-RBSP_5L_50per.ps
#      IonPairs_5L_Mar2013_50per.ps
#      NO_Fig3b_RBSP-FU_L5_50per.ps
#
#  Notes:
#   Modified on 1Apr17 to calculate ion pairs and NO estimate from
#   FU-3 examples.  Use FU3_flux_loss_cone.csv input file
#   Modifies on 9Jul17 to look at minimum outside of loss cone region
#   18 Nov, 2017  RBSP-FB conjunctions
#   Modified on 10May17 to consider different ratios based on energy level
#
#   Modified to extend ratios outside of FIREBIRD range -- currently extrapolating exponentially
#   at higher energies, but keeping constant at < 50 keV thinking that auroral processes will take over
#
#   Loss cone precipitation - integrate over steradians - multiply the flux by
#   2pi*(1-cos(66.3)) = 3.758 to take into account steradians of the loss cone
#   (Dan Marsh uses a value of 66.3 from Marshall and Bortnik
#   (https://doi.org/10.1002/2017JA024873)
# ---------------------------------

import numpy as np, os

#
#  Import Ngl support functions.
#
import Ngl
import numpy, os 
import Nio

from scipy import interpolate
from scipy.interpolate import interp1d

from math import exp,sqrt,cos,floor

# -------------------------------------------------------------------
# Set energy levels and altitude levels. 
# Energy levels from 100ev to 1 MeV  (0.1-1000. keV)
#  
# In CEDAR proposal, we propose to study 30 keV  to 4 MeV
# We might want to change
#     imax(1) energy_min, energy_max
#  energy max has been extended to 5MeV
#
# Question: How can we extrapolate to include higher energies?
# Fang only considers to 1 MeV
#
# datetime refers to row of .csv file
#
# FB_ratio is set to 10 energy levels   
# Read in .csv file with 10 energies (rows) and 4 columns 
# (energy, % FB/RB 50, %FB/RB 75, %FB/RB 100)
# Note that the energies are based on FIREBIRD flux energies.
#
# Set multiply to enhance precipitation by 10, 100, etc.
# -------------------------------------------------------------------
 
filenameoutnc = "Ion_pairs_RBSP-FU_Mar2013_L4-5_50per.nc"

#ts = 11
ts = 29
nplev = 201
enrbsp = 14
enrbsp1 = 15

datetime = 4 

#FB_ratio = 50. #min
#FB_ratio = 1. #none 

#multiply = 100.
multiply = 1.

losscone = 3.758

#imax = 60
#imax1 = 59

imax = 120 
imax1 = 119 

#energy_min = 0.1
energy_min = 2.
#energy_max = 1000.
energy_max = 2000.

energy_lev = np.zeros((imax),'f')
energy_lev[0]=energy_min
energy_lev[imax1]=energy_max
Prod_all = np.zeros((ts,nplev),'f')
pressure_new = np.zeros((nplev),'f') 

flaglo = 0
flaghi = 0
for nlev in range (0,imax1):
   ALO=np.log10(energy_lev[0])
   AHI=np.log10(energy_lev[imax1])
   SP=(AHI-ALO)/imax
   energy_lev[nlev]=10.**(ALO+SP*nlev)
   if energy_lev[nlev] >= 265.4 and flaglo == 0:
     indexFBlow = nlev
     flaglo = 1
   if energy_lev[nlev] >= 913 and flaghi == 0:
     indexFBhigh = nlev
     flaghi = 1
   nlev=nlev+1

print("Energy Levels ", energy_lev[:])
print("index FB low for E=265.4 keV", indexFBlow, energy_lev[indexFBlow])
print("index FB high for E=913 keV ", indexFBhigh, energy_lev[indexFBhigh])

#energy_rbsp = [56.7, 80.4, 110.9, 145.9, 184.9, 221.1, 335.5, 458.2, 592.1, 737.4, 893.7, 1058.9, 1556., 1697.]
energy_rbsp = [58., 82., 110., 145., 182., 221., 338., 460., 597., 741., 879., 1088., 1650., 1768.]

# --------------------------------------------------------------
# Read in FB-RBSP ratio file
# FB-RBSP_ratios.csv
# This file includes both FU3 and FU4 ratios
#
# Remember that we need to divide ratios by 100.
#
# Interpolate ratios to RBSP energy levels
# Consant extrapolation beyond FB-RBSP ratios (extending first and last values)
# -------------------------------------------------------------

data_ratio = np.zeros((10,4),dtype=np.float64)
ratio_energies = np.zeros((10),dtype=np.float64)
FB_ratios = np.zeros((10,3),'f')
FB_ratios_enlev = np.zeros((imax,3),'f')
FB_ratios_enrbsp = np.zeros((14,3),'f')
FB_ratio_50 = np.zeros((10),'f')
FB_ratio_75 = np.zeros((10),'f')
FB_ratio_100 = np.zeros((10),'f')
FB_ratio_50_enlev = np.zeros((imax),'f')
FB_ratio_75_enlev = np.zeros((imax),'f')
FB_ratio_100_enlev = np.zeros((imax),'f')
FB_ratio_50_enrbsp = np.zeros((14),'f')
FB_ratio_75_enrbsp = np.zeros((14),'f')
FB_ratio_100_enrbsp = np.zeros((14),'f')

filename_3 = "FB-RBSP_ratios.csv"

data_3     = Ngl.asciiread(filename_3,(10,4),"float",sep=",")

data_ratio[:,:] = data_3[:,:]

ratio_energies[:] = data_3[:,0]
FB_ratios[:,:] = data_3[:,1:]/100.
FB_ratio_50[:] = FB_ratios[:,0]
FB_ratio_75[:] = FB_ratios[:,1]
FB_ratio_100[:] = FB_ratios[:,2]

FB_enrbsp_interp = interp1d(ratio_energies, FB_ratios,kind='linear',axis=0,bounds_error=False,fill_value=0.)

FB_ratios_enlev[:,:] = FB_enrbsp_interp(energy_lev)
FB_ratio_50_enlev[:] = FB_ratios_enlev[:,0] 
FB_ratio_75_enlev[:] = FB_ratios_enlev[:,1] 
FB_ratio_100_enlev[:] = FB_ratios_enlev[:,2] 

FB_ratios_enrbsp[:,:] = FB_enrbsp_interp(energy_rbsp)
FB_ratio_50_enrbsp[:] = FB_ratios_enrbsp[:,0] 
FB_ratio_75_enrbsp[:] = FB_ratios_enrbsp[:,1] 
FB_ratio_100_enrbsp[:] = FB_ratios_enrbsp[:,2] 

# Extrapolate above and below FIREBIRD energies 
FB_ratio_50_enlev[0:indexFBlow] = FB_ratio_50[0]
FB_ratio_50_enlev[indexFBhigh:imax] = FB_ratio_50[9]
FB_ratio_75_enlev[0:indexFBlow] = FB_ratio_75[0]
FB_ratio_75_enlev[indexFBhigh:imax] = FB_ratio_75[9]
FB_ratio_100_enlev[0:indexFBlow] = FB_ratio_100[0]
FB_ratio_100_enlev[indexFBhigh:imax] = FB_ratio_100[9]

print("energy_lev index_high", energy_lev[indexFBhigh])
print("FB_ratio_50_enrbsp", FB_ratio_50_enrbsp[:])
print("FB_ratio50[0]", FB_ratio_50[0])
print("FB_ratio50[9]", FB_ratio_50[9])
print("FB_ratio_50_enlev", FB_ratio_50_enlev[:])

# -------------------------------------------------------------
# Read in netCDF file of differential energy flux from FIREBIRD data. 
#       Electron Energy (E) in keV
#       Qmono(E) in cm-2 s-1 sr-1 keV-1
# from csv files saved from excel files firebird_flux.csv
#
# Interpolate to energy levels
# Multiply by loss cone solid angle
#
# Multiply by any hypothetical multiplication factors
#
# -------------------------------------------------------------------

energy_fire = energy_rbsp

data_fire = np.zeros((ts,enrbsp+1),dtype=np.float64)
julian = np.zeros((ts),dtype=np.float64)
Qmono_all = np.zeros((ts,enrbsp),'f')

filename_fire = "RBSP_FEDU_files/RBSPa_FEDU_2013Feb15-Mar15_L4-5.csv"

data1     = Ngl.asciiread(filename_fire,(ts,enrbsp1),"double",sep=",")

data_fire[:,:] = data1[:,:]

julian[:] = data_fire[:,0]
#Qmono_all[:,:] = data_fire[:,1:]
Qmono_all[:,:] = data_fire[:,1:]*losscone*multiply

#print("date", julian[:])
#print("ts=0", Qmono_all[0,:])

# -----------------------------------------------------------
# Read in Table 1 of Fang et al., 2010 of Parameterization 
# Coefficients for Isotropically Incident Monoenergetic 
# 100 eV to 1 MeV electrons 
# 
# Currently in Fang_Table1.csv
# ----------------------------------------------------------

Pij = np.zeros((8,4),dtype=np.float64)

filename_table = "Fang_Table1.csv"
 
data2 = Ngl.asciiread(filename_table,(8,4),"double",sep=",")

Pij[:,:] = data2[:,:]
 
#print("Fang_Table1 i = 1 ", Pij[0,:])
#print("Fang_Table1 j = 0 ", Pij[:,0])

# ------------------------------------------------------------
# Read in MSIS90 temperature and density as a function of height
# from 0-400 km in intervals of 2 km
# 
# http://omniweb.gsfc.nasa.gov/vitmo/msis_vitmo.html
#
# MSIS90.csv  height(km), density(g cm-3), temp(k)
#
# Calculate pressure = rho*rd*temp
#    which should be a reasonable value throughout the altitudes of interest
#    (in homosphere) 
#    gas constant for dry air -- rd = 287 
#    need to multiply density by 1000 to get kg m-3
#
#  Rerun MSIS90 to get reasonable date for each simulation.
#  Question: Will we need to do this monthly? Daily?
# ------------------------------------------------------------

rd = 287.

data_MSIS90 = np.zeros((nplev,3),dtype=np.float64)

filename_MSIS90 = "MSIS90_4Mar13.csv"
 
data3 = Ngl.asciiread(filename_MSIS90,(nplev,3),"double",sep=",")

data_MSIS90[:,:] = data3[:,:]

z = np.zeros(nplev)
rho = np.zeros(nplev)
temp= np.zeros(nplev)

z[:] = data_MSIS90[:,0]
rho[:] = data_MSIS90[:,1]
temp[:] = data_MSIS90[:,2]

pressure_new[:] = rd*1000.*rho[:]*temp[:]
 

# -----------------------------------------------------------
# Note that height is in km, density is in g cm-3 and 
# temp is in K
# -----------------------------------------------------------

#print("height ", z[:])
#print("density ", rho[:])
#print("temp ", temp[:])
#print("pressure ", pressure_new[:])

# -----------------------------------------------------
# Set acceleration of gravity and average molecular weight profiles
# and calculate scale height profile
#
# Molecular weight profiles from Handbook of Geophysics 1985
# Assume exospheric temperature = 900 K (from MSIS90 9Dec13)
#
# k - Boltzman constant = 1.38e-16 cm2 g s-2 K-1
# 
# -----------------------------------------------------

k = 1.381e-16
Av = 6.022e23
g0 = 9.807
re = 6371.
gmks = np.zeros(nplev)
g = np.zeros(nplev)
gmks[:] = g0*(re/(re+z[:]))**2.
g[:] = gmks[:]*100.

#print("acceleration of gravity profile ", g[:])

mwm = np.zeros(nplev)
mw = np.zeros(nplev)

#mwm_z = [75.,80.,85.,90.,95.,100.,105.,110.,115.,120.]
#mwm_data = [28.96, 28.95, 28.93, 28.89, 28.81, 28.37, 27.51, 26.73, 26.05, 25.45] 
mwm_z = [75.,80.,85.,90.,95.,100.,105.,110.,115.,120.,130.,140.,150.,160.,180.,200.,250.,300.,350.,400.]
mwm_data = [28.96, 28.95, 28.93, 28.89, 28.81, 28.37, 27.51, 26.73, 26.05, 25.45, 24.39, 23.51, 22.73, 22.03, 20.81, 19.79, 17.97, 16.91, 16.26, 15.74] 

mwm_interp = interp1d(mwm_z, mwm_data,kind='linear',bounds_error=False,fill_value=0.)

mw[:] = mwm_interp(z[:])/Av   
mw[0:38] = 28.97/Av

#print("height", z[:])
#print("molecular weight", mw[:])

H = np.zeros(nplev)

H[:] = k*temp[:]/(mw[:]*g[:])

#print("scale height ", H[:])

# Scale height is in cm
 
# ---------------------------------------------------------------------
# Loop through energy levels and calculate Ci coefficents
#
# Emono -- energy_lev[60]
#
# Parameterization Coefficiences Pij  --   Pij[8,4] 
# 
# --------------------------------------------------------------------- 

coeff = np.zeros((8,imax))
en_ln = np.zeros(imax)
en_ln_j = np.zeros(4)
en_ln_j_P = np.zeros(4)

en_ln[:] = np.log(energy_lev[:]) 
#print("Pij...j values",Pij[1,:])

for n_en in range (0,imax):
#   for i in range (1,8):
   for i in range (0,8):
      sum_en_ln_j_P = 0.
      for j in range (0,4):
         en_ln_j[j] = en_ln[n_en]**j
         en_ln_j_P[j] = en_ln_j[j]*Pij[i,j]
         sum_en_ln_j_P = sum_en_ln_j_P+en_ln_j_P[j]
      coeff[i,n_en] = np.exp(sum_en_ln_j_P)
           
#print("coeff at E = 1keV", coeff[:,15])
# --------------------------------------------------------------------
# Calculate normalized atmospheric column mass (y) as a function of 
# height (z)
#
# use rho(z), H(z), energy_lev(imax)
# --------------------------------------------------------------------

y = np.zeros((nplev,imax))
y_fact1 = np.zeros(nplev)

y_fact1[:] = (rho[:]*H[:]/6.e-6)**0.7
for n_en in range (0,imax):
   y[:,n_en] = 2.*y_fact1[:]/energy_lev[n_en]

#print("y at E = 1MeV", y[:,imax1])
#print("y at E = 1keV", y[:,15])
# --------------------------------------------------------------------
# Calculate normalized energy dissipation f(height,energy_lev)
#
# coeff(8,imax), y(imax,energy_lev), energy_lev(imax) 
# -------------------------------------------------------------------

f = np.zeros((nplev,imax))
f_part1 = np.zeros(nplev)
f_part2 = np.zeros(nplev)
f_part3 = np.zeros(nplev)
f_part4 = np.zeros(nplev)

f_part1[:] = 0.
f_part2[:] = 0.
f_part3[:] = 0.
f_part4[:] = 0.

for n_en in range (0,imax):
   f_part1[:] = coeff[0,n_en]*(y[:,n_en]**coeff[1,n_en])
   f_part2[:] = np.exp(-coeff[2,n_en]*(y[:,n_en]**coeff[3,n_en]))
   f_part3[:] = coeff[4,n_en]*(y[:,n_en]**coeff[5,n_en])
   f_part4[:] = np.exp(-coeff[6,n_en]*(y[:,n_en]**coeff[7,n_en]))
   f[:,n_en] = f_part1[:]*f_part2[:] + f_part3[:]*f_part4[:]

#print("f at E = 1MeV", f[:,imax1])
#print("f at E = 1MeV", f[:,15])


# ---------------------------------------------------------------
# Plot Figure 1 of Fang
#
# 100 eV -> 0
# 1 keV -> 15
# 10 keV -> 30
# 100 keV -> 45
# 1 MeV ->59 
#
# horizontal f(:,en) vs. vertical y(:,en)
# --------------------------------------------------------------- 

# ---------------- Begin plots -------------------------------

wks_type = "ps"
wks = Ngl.open_wks(wks_type,"Fang_Fig1_ext")  # Open a workstation.

#
#  Set resources for titling.
#
resources = Ngl.Resources()

#
resources = Ngl.Resources()
resources.tiXAxisString = "Normalized Energy Deposition, f"
resources.tiYAxisString = "Normalized Atmospheric Column Mass, y"
resources.trYAxisType = "LogAxis"
resources.trYMaxF = 10
resources.trYMinF = 0.1 
resources.trXMaxF = 0.8 
resources.trXMinF = 0 

plot = Ngl.xy(wks,f[:,45],y[:,45],resources)   # Draw an XY plot.
#plot = Ngl.xy(wks,f[:,15],y[:,15],resources)   # Draw an XY plot.

# ---------------- End Plot Figure 1 ------------------------
# Begin plot Figure 2
# -----------------------------------------------------------

Qtot_fig2 = np.zeros(nplev)
Prod_fig2 = np.zeros(nplev)

Qmono_fig2 = 0.

# ------- set Qmono_fig2 for energy level --------------------
# Energy should equal 1 erg, recalling that 1 keV = 1.6e-9 erg
# ------------------------------------------------------------
n_en = 0 
Qmono_fig2 = 1.6e9
 
for nz in range (0,nplev):
   Qtot_sum_fig2 = 0.
   Qtot_fig2[nz] = f[nz,n_en]*Qmono_fig2/(0.035*H[nz])
   Qtot_sum_fig2 = Qtot_sum_fig2+Qtot_fig2[nz]
   Prod_fig2[nz] = Qtot_sum_fig2
   if Prod_fig2[nz] <= 1.e-10:
      Prod_fig2[nz] = 1.e-10

#print("Prod_fig2 = ", Prod_fig2[:])

 
# ---------------- plot for figure 2 -----------------------------

wks_type = "ps"
wks_fig2 = Ngl.open_wks(wks_type,"Fang_Fig2")  # Open a workstation.

#
#  Set resources for titling.
#
resources2 = Ngl.Resources()

#
resources2 = Ngl.Resources()
resources2.tiXAxisString = "Ion pair production rate (cm-3 s-1)"
resources2.tiYAxisString = "Altitude (km)"
resources2.trXAxisType = "LogAxis"
resources2.trYMaxF = 400. 
resources2.trYMinF = 50. 
resources2.trXMaxF = 1.e5 
#resources2.trXMinF = 1.e-1 
resources2.trXMinF = 10. 

plot2 = Ngl.xy(wks_fig2,Prod_fig2[:],z[:],resources2)   # Draw an XY plot.

# ----------------  End plot ---------------------------------

diff_flux = np.zeros((ts,imax),'f')

# ------------------------------------------
# Begin looping through the timesteps
# Interpolate to new energy levels -- exponential 
# ------------------------------------------

ntime = 0
for ntime in range (0,ts):
   Qmono_fire = np.zeros((enrbsp),'f')
   Qmono_fire_log = np.zeros((enrbsp),'f')
   Qmono = np.zeros((imax),'f')
   Qmono_new = np.zeros((imax),'f')
   Qmono_log = np.zeros((imax),'f')

   Qmono_fire[:] = Qmono_all[ntime,:] 

#   print("FIREBIRD energy", energy_fire[:])
#   print("energy_lev", energy_lev[:])
 
#   print("FIREBIRD flux", Qmono_fire[:])

#   Qmono_fire_log[:] = np.log10(Qmono_fire[:])
   Qmono_fire_log[:] = np.log(Qmono_fire[:])

# Linear least squares fit to log data

   energy_array = np.zeros((enrbsp,2),'f')
#   energy_array[:,1] = np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]) 
   energy_array[:,1] = np.array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]) 
   energy_array[:,0] = energy_fire[:]
 

   m, yint = np.linalg.lstsq(energy_array,Qmono_fire_log)[0]
#   print("ts, slope, int")
#   print(ntime,m,yint)
  
#   a = 10**m
#   C = 10**yint 

   a = exp(m)
   C = exp(yint) 
  
#   print("ts, a, C")
#   print(ntime,a,C)

   Qmono_new[:] = C*(a**energy_lev[:])

# Create tail continuing exp fit to last two data points.

#   slope = (Qmono_fire_log[4]-Qmono_fire_log[3])/(energy_fire[4]-energy_fire[3])
   slope = (Qmono_fire_log[enrbsp-2]-Qmono_fire_log[enrbsp-1])/(energy_fire[enrbsp-2]-energy_fire[enrbsp-1])
   if slope <= 0.:
      energy_fire_fin = energy_max
#      Qmono_fire_log_fin = Qmono_fire_log[4]+slope*(energy_fire_fin-energy_fire[4])
      Qmono_fire_log_fin = Qmono_fire_log[enrbsp-2]+slope*(energy_fire_fin-energy_fire[enrbsp-2])
   else:
      energy_fire_fin = energy_max
      Qmono_fire_log_fin = Qmono_fire_log[enrbsp-2]

# Set lower energies as constant
   energy_fire_start = energy_min
   Qmono_fire_log_start = Qmono_fire_log[0]

# Extrapolate exponentially for lower energies

#   slope2 = (Qmono_fire_log[1]-Qmono_fire_log[0])/(energy_fire[1]-energy_fire[0])
#   energy_fire_start = energy_min
#   Qmono_fire_log_start = Qmono_fire_log[0]-slope2*(energy_fire[0]-energy_fire_start)

   energy_fire_ext = np.zeros((enrbsp+2),'f')
   energy_fire_ext[0] = energy_fire_start
   energy_fire_ext[1:enrbsp+1] = energy_fire[:]
   energy_fire_ext[enrbsp+1] = energy_fire_fin

   Qmono_fire_log_ext = np.zeros((enrbsp+2),'f')
   Qmono_fire_log_ext[0] = Qmono_fire_log_start
   Qmono_fire_log_ext[1:enrbsp+1] = Qmono_fire_log[:]
   Qmono_fire_log_ext[enrbsp+1] = Qmono_fire_log_fin

#   print("Qmono_fire_log_ext", Qmono_fire_log_ext[:])

#   en_interp = interp1d(energy_fire, Qmono_fire_log,kind='linear',bounds_error=False,fill_value=-6.)
   en_interp = interp1d(energy_fire_ext, Qmono_fire_log_ext,kind='linear',bounds_error=False,fill_value=-6.)

   Qmono_log[:] = en_interp(energy_lev[:])   

#   print("Qmono_log", Qmono_log[:])

   for n_en in range (0,imax):
#      if Qmono_log[n_en] > 0.:
      if Qmono_log[n_en] > -6.:
#         Qmono[n_en] = 10.**Qmono_log[n_en]
         Qmono[n_en] = np.exp(Qmono_log[n_en])
      else:
         Qmono[n_en] = 0.   

   
#   print("Qmono", Qmono[:])

#   diff_flux[ntime,:] = Qmono[:]
#   diff_flux[ntime,:] = Qmono_new[:]

# -----------------------------------------------------------------
# ************  HERE NEED TO INTERPOLATE FB-RATIO TO ENERGY LEVELS

#   diff_flux[ntime,:] = Qmono[:]/FB_ratio
#   diff_flux[ntime,:] = Qmono_new[:]/FB_ratio
   diff_flux[ntime,:] = Qmono[:]*FB_ratio_50_enlev[:]


#   print("DoubleMAX interpolated", Qmono)

# ---------------------------------------------------------------
# Calculate ionization profile Qtot(height,energy_lev)
# and Prod(height) -- total ion pair production 
#
# Prod_fig2(height, energy_lev), Q_mono(imax), H(height), eps = 0.035 
#
# delta epsilon (eps) is 0.035 keV required to produceone ion pair. 
# ---------------------------------------------------------------

   Qenergy = np.zeros((nplev,imax))
   Prod = np.zeros(nplev)
   Qenergy[:,:] = 0. 

   for nz in range (0,nplev):
      for n_en in range (0,imax):
         Qenergy[nz,n_en] = f[nz,n_en]*1.0/(0.035*H[nz])
      

   for nz in range (0,nplev):
      Qtot_sum = 0.
      for n_en in range (0,imax-1):
#      for n_en in range (0,imax):
         Qtot = 0.
#         Qtot[n_en] = f[nz,n_en]*Qmono[n_en]/(0.035*H[nz])
#         Qtot = Qenergy[nz,n_en]*Qmono[n_en]*energy_lev[n_en]*(energy_lev[n_en+1]-energy_lev[n_en])
         Qtot = Qenergy[nz,n_en]*Qmono[n_en]*energy_lev[n_en]*(energy_lev[n_en+1]-energy_lev[n_en])*FB_ratio_50_enlev[n_en]
#         Qtot = Qenergy[nz,n_en]*Qmono_new[n_en]*energy_lev[n_en]*(energy_lev[n_en+1]-energy_lev[n_en])
         Qtot_sum = Qtot_sum+Qtot
         Prod[nz] = Qtot_sum

#   print("Prod = ", Prod[:])

   Prod_all[ntime,:] = Prod[:]

 
# ---------------- plot for figure 3 -----------------------------
# Specify timestep in Prod_all[ts,:]
# Not sure if this is right... 
# ----------------------------------------------------------------

wks_type = "ps"
wks_fig3 = Ngl.open_wks(wks_type,"IonPairs_L4-5_Mar2013_50per")  # Open a workstation.

Prod_fig3 = np.zeros(nplev)
Prod_fig3_all = np.zeros([ts,nplev],'f')
#Prod_fig3[:] = Prod_all[0,:]
Prod_fig3[:] = Prod_all[datetime,:]
Prod_fig3_all[:,:] = Prod_all[0:ts,:]

#
#  Set resources for titling.
#
resources3 = Ngl.Resources()
resleg_res3 = Ngl.Resources()

#
resources3 = Ngl.Resources()
resources3.tiXAxisString = "Ion pair production rate (cm-3 s-1)"
resources3.tiYAxisString = "Altitude (km)"
#resources3.trXAxisType = "LogAxis"
resources3.xyLineColors = ["red","brown","hot pink","orange","goldenrod","dark green","green","aquamarine","blue","purple","black"]
cmap_res3 = ["white","black","purple","blue","aquamarine","green","dark green","goldenrod","orange","hot pink","brown","red"]
Ngl.define_colormap(wks_fig3,cmap_res3) 
#resources3.xyLineColors = ["CornFlowerBlue", "DeepPink4", "DarkSeaGreen","dark green","cadet blue","khaki","goldenrod","Chartreuse2","salmon","DarkSlateBlue","DodgerBlue","DarkOrchid3","hot pink","Firebrick","lavender","Burlywood", "Coral3","honeydew","navy","LightBlue","turquoise","aquamarine","red","blue","green","orange","brown","purple","cyan","black","yellow","gray"]
#resources3.trYMaxF = 200. 
#resources3.trYMaxF = 120. 
#resources3.trYMaxF = 105. 
resources3.trYMaxF = 105. 
resources3.trYMinF = 45. 
#resources3.trXMaxF = 80. 
#resources3.trXMaxF = 1.e4
resources3.trXMaxF = 200.
#resources3.trXMinF = .001 
resources3.trXMinF = 0 
resources3.xyLineThicknessF = 3.0

resleg_res3.lgLabelFontHeightF = 0.014
resleg_res3.vpWidthF          = 0.2
resleg_res3.vpHeightF         = 0.35
resleg_res3.lgLineThicknessF  = 4.0

#for its in range (0,ts):
for nz in range (0,nplev):
   if Prod_fig3[nz] <= 1.e-10: 
      Prod_fig3[nz] = 1.e-10

#print("Prod = ", Prod_fig3[:])

labels_res3 = ["Mar 14","Mar 13","Mar 12","Mar 11","Mar 10","Mar 9","Mar 8","Mar 7","Mar 6","Mar 5","Mar 4"]

#plot3 = Ngl.xy(wks_fig3,Prod_fig3[:],z[:],resources3)   # Draw an XY plot.
#lg = Ngl.legend_ndc(wks_fig3,11,labels_res3,0.3,0.9,resleg_res3)
#lg = Ngl.legend_ndc(wks_fig3,11,labels_res3,0.7,0.93,resleg_res3)
#plot3 = Ngl.xy(wks_fig3,Prod_fig3[:],z[:],resources3)   # Draw an XY plot.
plot3 = Ngl.xy(wks_fig3,Prod_fig3_all[:,:],z[:],resources3)   # Draw an XY plot.


# ----------------  End plot ---------------------------------

# ---------------- plot for figure 3b -----------------------------
# Specify timestep in Prod_all[ts,:]
# Not sure if this is right... 
# ----------------------------------------------------------------

wks_type = "ps"
wks_fig3b = Ngl.open_wks(wks_type,"NO_Fig3b_RBSP-FU_L4-5_50per")  # Open a workstation.

NO_fig3b = np.zeros(nplev)
NO_fig3b_all = np.zeros([ts,nplev],'f')
NO_fig3b[:] = 0.69*Prod_fig3[:]
for nz in range (0,nplev):
   for nts in range (0,ts):
      NO_fig3b_all[nts,nz] = 0.69*Prod_fig3_all[nts,nz]

#
#  Set resources for titling.
#
resources3b = Ngl.Resources()
resources3b.tiXAxisString = "NO production rate (cm-3 s-1)"
resources3b.tiYAxisString = "Altitude (km)"
#resources3b.trXAxisType = "LogAxis"
resources3b.trYMaxF = 105. 
#resources3b.trYMinF = 55. 
resources3b.trYMinF = 45. 
resources3b.trXMaxF = 200. 
#resources3b.trXMaxF = 5000.
#resources3b.trXMinF = .00001 
resources3b.trXMinF = 0. 
resources3b.xyLineThicknessF = 3.0
resources3b.xyLineColors = ["red","brown","hot pink","orange","goldenrod","dark green","green","aquamarine","blue","purple","black"]
cmap = ["white","black","purple","blue","aquamarine","green","dark green","goldenrod","orange","hot pink","brown","red"]
Ngl.define_colormap(wks_fig3b,cmap) 

#print("NO = ", NO_fig3b[:])
#lg = Ngl.legend_ndc(wks_fig3b,11,labels_res3,0.70,0.93,resleg_res3)

#plot3b = Ngl.xy(wks_fig3b,NO_fig3b[:],z[:],resources3b)   # Draw an XY plot.
plot3b = Ngl.xy(wks_fig3b,NO_fig3b_all[:,:],z[:],resources3b)   # Draw an XY plot.

# ----------------  End plot ---------------------------------


# --------------- plot for figure 4 -------------------------------

Prod_plot = np.zeros((nplev,ts),'f')
Prod_plot = np.rollaxis(Prod_all,1)
julian_plot = np.zeros((ts),'f')
#julian_plot[:] = 100000.*(julian[:]-2457050.)
julian_plot[:] = 100000.*(julian[:]-2457055.)

#print(julian_plot)

wks_type = "ps"
wks_fig4 = Ngl.open_wks(wks_type,"RBSP-FU_ions_Mar2013_L4-5_50per")  # Open a workstation.

resources4 = Ngl.Resources()

resources4.cnFillOn = True
resources4.tiXAxisString = "Time"
resources4.tiYAxisString = "Altitude (km)"
resources4.trXTensionF = 1.
#resources4.sfXArray = julian_plot
resources4.sfYArray = z
resources4.trYMaxF = 120.
resources4.trYMinF = 0.
#resources4.trXMinF = 575780.
#resources4.trXMaxF = 575870.
resources4.cnLevelSelectionMode = "ExplicitLevels"
#resources4.cnLevels = [1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3]
#resources4.cnLevels = [20.,40.,60.,80.,100.,120.,140.,160.,180.]
resources4.cnLevels = [10.,20.,50.,100.,200.,500.,1000.,2000.,5000.,10000.]
resources4.cnLineLabelsOn = False

#plot4 = Ngl.contour(wks_fig4,Prod_plot[:,0:4],resources4)   # Draw an XY plot.
plot4 = Ngl.contour(wks_fig4,Prod_plot[:,0:ts],resources4)   # Draw an XY plot.

# ---------------  End plot --------------------------------


# --------------- plot for figure 5 -------------------------------

#Prod_lshell = np.zeros((nplev,ts),'f')
#julian_lshell = np.zeros((ts),'f')

#nstep = 0
#newstep = 0
#for nstep in range (0,ts):
#   if lshell[nstep] >= 4.5 and lshell[nstep] < 5.5:
#      Prod_lshell[:,newstep] = Prod_plot[:,nstep]
#      julian_lshell[newstep] = julian_plot[nstep]
#      newstep = newstep+1

#print(julian_lshell)

#Prod_lshell_new = np.zeros((nplev,newstep),'f')
#julian_lshell_new = np.zeros((newstep),'f')
#Prod_lshell_new[:,:] = Prod_lshell[:,0:newstep]
#for nstep in range (0,newstep):
#  julian_lshell_new[nstep] = round(julian_lshell[nstep]/100.,0)

#print(julian_lshell_new)

#wks_type = "ps"

#wks_fig5 = Ngl.open_wks(wks_type,"Firebird_lshell")  # Open a workstation.

#resources5 = Ngl.Resources()

#resources5.cnFillOn = True
#resources5.tiXAxisString = "Time"
#resources5.tiYAxisString = "Altitude (km)"
#resources5.trXTensionF = .001 
#resources5.sfXArray = julian_lshell_new
#resources5.sfXArray = [0,1,3,4,6,7,9,10]
#resources5.sfYArray = z
#resources5.trYMaxF = 120.
#resources5.trYMinF = 0.
#resources5.trXMinF = 691.
#resources5.trXMaxF = 1291.

#plot5 = Ngl.contour(wks_fig5,Prod_lshell_new,resources5)   # Draw an XY plot.

# ---------------  End plot --------------------------------

# ---------------- plot for figure 6 -----------------------------
# specify timestep diff_flux_plot[:] = diff_flux[1,:]
# firebird_flux[:] = Qmono_all[1,:]
# For loss cone 2 Feb 15 use 2 and 5
# ------------------------------------------------------------------

wks_type = "ps"
wks_fig6 = Ngl.open_wks(wks_type,"Diff_flux_FIRE_FB-RBSP_L4-5_50per")  # Open a workstation.

resources6 = Ngl.Resources()
resources6a = Ngl.Resources()
resleg = Ngl.Resources()

resources6.tiYAxisString = "Differential Energy Flux (cm-2 s-1 sr-1 keV-1)"
resources6.tiXAxisString = "Energy (keV)"
resources6.trXAxisType = "LogAxis"
resources6.trYAxisType = "LogAxis"
resources6.trYMaxF = 1.e4 
#resources6.trYMaxF = 1.e3 
#resources6.trYMaxF = 1.e5 
resources6.trYMinF = 1.e-3 
#resources6.trXMaxF = 1.e5 
resources6.trXMaxF = 1.e4 
#resources6.trXMinF = 10. 
resources6.trXMinF = 50. 
#resources6.trYMaxF = 1.e3 
#resources6.trYMinF = 1.e-5 
#resources6.trXMaxF = 1.e4 
#resources6.trXMinF = 100. 
resources6.nglMaximize = False
resources6.vpHeightF   = 0.75
resources6.vpWidthF    = 0.65
resources6.vpXF        = 0.20
resources6.vpYF        = 0.90
#resources6.xyMarkLineModes = "Markers"
#resources6.xyMarkerColor     = "red"
#resources6.xyMarkers = 3 
resources6.xyLineThicknessF = 3.0

resources6.xyLineColors = ["red","brown","hot pink","orange","goldenrod","dark green","green","aquamarine","blue","purple","black"]
cmap = ["white","black","purple","blue","aquamarine","green","dark green","goldenrod","orange","hot pink","brown","red"]
Ngl.define_colormap(wks_fig6,cmap) 

resources6a.trXAxisType = "LogAxis"
resources6a.trYAxisType = "LogAxis"
resources6a.trYMaxF = 1.e4 
#resources6a.trYMaxF = 1.e3 
#resources6a.trYMaxF = 1.e5 
resources6a.trYMinF = 1.e-3 
#resources6a.trXMaxF = 1.e5 
resources6a.trXMaxF = 1.e4 
#resources6a.trXMinF = 10. 
resources6a.trXMinF = 50. 
#resources6a.trYMaxF = 1.e3 
#resources6a.trYMinF = 1.e-5 
#resources6a.trXMaxF = 1.e4 
#resources6a.trXMinF = 100. 
resources6a.xyMarkLineModes = "Markers"
resources6a.xyMarkerColor     = "black"
resources6a.xyMarkers = 3 
resources6a.nglMaximize = False
resources6a.vpHeightF   = 0.75
resources6a.vpWidthF    = 0.65
resources6a.vpXF        = 0.20
resources6a.vpYF        = 0.90


diff_flux_plot = np.zeros((imax),'f')
diff_flux_plot_all = np.zeros((ts,imax),'f')
#diff_flux_plot = np.zeros((32,imax),'f')
#diff_flux_plot[:] = diff_flux[8,:]
diff_flux_plot[:] = diff_flux[datetime,:]
diff_flux_plot_all[:,:] = diff_flux[:,:]

for its in range (0,ts):
#   print("its = ", its)
   for en in range (0,imax):
      if diff_flux_plot_all[its,en] <= 0.:
         diff_flux_plot_all[its,en] = 1.e-6

for en in range (0,imax):
   if diff_flux_plot[en] <= 0.:
      diff_flux_plot[en] = 1.e-6

firebird_flux = np.zeros((enrbsp),'f') 

#firebird_flux[:] = Qmono_all[23,:]
#firebird_flux[:] = Qmono_all[1,:]
#firebird_flux[:] = Qmono_all[datetime,:]
firebird_flux[:] = Qmono_all[datetime,:]*FB_ratio_50_enrbsp[:] 

for en in range (0,14):
   if firebird_flux[en] <= 0.:
      firebird_flux[en] = 1.e-6


print("ENERGY_LEV",energy_lev[21:])
print("FLUX",diff_flux_plot_all[datetime,21:])

resources6.nglDraw  = False
resources6.nglFrame  = False
resources6a.nglDraw  = False
resources6a.nglFrame  = False

resleg.lgLabelFontHeightF = 0.014
resleg.vpWidthF          = 0.13
resleg.vpHeightF         = 0.35
resleg.lgLineThicknessF  = 4.0

plot6 = Ngl.xy(wks_fig6,energy_lev[21:],diff_flux_plot_all[:,21:],resources6)   # Draw an XY plot.
plot7 = Ngl.xy(wks_fig6,energy_fire,firebird_flux,resources6a)   # Draw an XY plot.

#labels = ["Mar 4","Mar 5","Mar 6","Mar 7","Mar 8","Mar 9","Mar 10","Mar 11","Mar 12","Mar 13","Mar 14"]
labels = ["Mar 14","Mar 13","Mar 12","Mar 11","Mar 10","Mar 9","Mar 8","Mar 7","Mar 6","Mar 5","Mar 4"]

Ngl.draw(plot6)
Ngl.draw(plot7)
#lg = Ngl.legend_ndc(wks_fig6,11,labels,0.68,0.9,resleg)
Ngl.frame(wks_fig6)


# --------------------------------------------------------------
# Save netCDF file
# 
#dimensions:
#        time = UNLIMITED ; // (1488 currently)
#        pressure = nplev ;
#variables:
#        float pressure(pressure) ;
#        float Prod(time, pressure) ;
#        float time(time)  ;
#        double julian(time) ;
#
# Note: numpy and netcdf have conflicts in variable type
#
# convert pressure to hPa (currently in Pa) -- divide by 100 
# --------------------------------------------------------------

timeint = ts 

Prodnew = np.zeros((timeint,nplev),'f')
timewaccm = np.zeros((timeint),'float64')
julian_new = np.zeros((timeint),'float64')
lshell_new = np.zeros((timeint),'f')

timewaccm[:] = julian[:]-2415020.5
julian_new = julian
#lshell_new = lshell

if os.path.exists(filenameoutnc):
   os.remove(filenameoutnc)

ncfile = Nio.open_file(filenameoutnc,'w')

#print('Write netCDF file to ',filenameoutnc)
print("julian dates")
print(julian)

# create dimensions
nlevel = nplev 

ncfile.create_dimension('pressure', nlevel)
ncfile.create_dimension('time', None)

#    define variables
time = ncfile.create_variable('time','d',('time',))
julian = ncfile.create_variable('julian','d',('time',))
pressure = ncfile.create_variable('pressure','f',('pressure',))
Prod = ncfile.create_variable('Prod','f',('time','pressure'))
#lshell = ncfile.create_variable('lshell','f',('time',))

Prodnew[:,:]=Prod_all[:,:]
pressure_new = pressure_new/100.

ncfile.variables['time'][:] = timewaccm 
ncfile.variables['pressure'][:] = pressure_new[:]
ncfile.variables['Prod'][:,:] = Prodnew[:,:] 
ncfile.variables['julian'][:] = julian_new[:] 
#ncfile.variables['lshell'][:] = lshell_new[:] 


Ngl.end()
