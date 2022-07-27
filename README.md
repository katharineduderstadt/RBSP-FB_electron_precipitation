# RBSP-FB_electron_precipitation

These files were written for NSF grant #1650738 
https://nsf.gov/awardsearch/showAward?AWD_ID=1650738&HistoricalAwards=false 

The goal is to 

cesm2_FB-RBSP_Mar13_ionization_L4_corrected.py (and for L4-5, L5, L5-5)

  1) read in electron flux data in the lowest pitch angle bin from the Van Allen Probes RBSP-ECT datasets 
        INPUT: RBSP_FEDU_files/RBSPa_FEDU_2013Feb15-Mar15_L5.csv
        
  2) scale this flux to ratios from FIREBIRD to RBSP during conjunctions (given that FIREBIRD is in 
     low Earth orbit near the top of Earth's atmosphere (400-600) km for each L-shell bin 
     [From L-shells 3-7 currently centered around 4 (3-4.25), 4.5 (4.25-4.75), 5 (4.75-4.25), and 5.5 (5.25-7)] 
        INPUT: FB-RBSP_ratios.csv  (.5 MLT and .5 Lshell conjunctions in 50, 75, 100 percentiles)
        
  3) multiply this flux by the loss cone solid angle to estimate precipitation into the atmosphere 
        2pi*(1-cos(66.3)) = 3.758
        
  4) Calculate atmospheric ionization profiles based on the Fang et al. (2010) methods. 
      Currently using MSIS90 atmosphere 
        INPUT: MSIS90_4Mar13.csv
               Fang_Table1.csv
               
  OUTPUT: Ion_pairs_RBSP-FU_Mar2013_L4_50per.nc (and for L4-5, L5, L5-5)

cesm2_Ion_pairs_regrid_RBSP-FB_3-7L.py

  5) Compile L-shell ionization files 
      Ion_pairs_RBSP-FU_Mar2013_L4_50per.nc (and for other L-shells) 
      Jackman_Pressure.csv (pressure levels from Charley Jackman's code -- 
      every 2 km...same as traditionally used in WACCM solar proton input files) 
      
      INPUT: Ion_pairs_RBSP-FU_Mar2013_L4_50per.nc (and for L4-5, L5, L5-5)
      
      OUTPUT: cesm2_REP_ion_pairs_RBSP-FB.nc 
      
cesm2_Add_dates_ion_pairs_RBSP.ncl
  6) Configure format (especially dates) for input file to read into WACCM 
  
      INPUT: cesm2_REP_ion_pairs_RBSP-FB.nc
      
      OUTPUT: cesm2_MEE_RBSP_WACCM.nc
  
  
  
  
