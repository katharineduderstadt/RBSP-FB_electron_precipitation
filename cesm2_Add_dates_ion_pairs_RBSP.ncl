; ===============================================
; cesm2_Add_dates_ion_pairs.ncl
;
; Last modified 29Dec20 
; 
; SYNOPSIS
; Used to add dates in mmddyyyy and datesec to 
; REP_ion_pairs.nc file so that the format is consistent
; with the ion pair files
;
; INPUT FILES
; cesm2_REP_ion_pairs_RBSP-FB.nc
;
; OUTPUT FILES
; cesm2_MEE_RBSP_WACCM.nc
;
; NOTES
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; ================================================;
begin

  filenameREPions = "cesm2_REP_ion_pairs_RBSP-FB.nc" 

  filenameWACCM = "cesm2_MEE_RBSP_WACCM.nc"
  
;  filenameSPEall = "spes_1963-2014_c150717.nc"

; --------------------------------------------------
; Read in netCDF file for differential energy flux 
; --------------------------------------------------

  f = addfile(filenameREPions,"r")
 
  rep_time = f->time
  rep_plev = f->plev
  rep_glat = f->glat
  rep_ion_orig = f->rep_ion 
 
; --------------------------------------------------
; Read in netCDF file for differential energy flux 
; --------------------------------------------------

;  fall = addfile(filenameSPEall,"r")

;  daily_date = fall->date 
;  daily_datesec = fall->datesec 

;  dailytimeint = 18994 
;  timeallint = 20449 

; ------------------------------------------
; Create output file as input to spedata.F90 WACCM
;	pressure = 58 ;
;	time = UNLIMITED ; // 
;variables:
;	float pressure(pressure) ;
;		pressure:long_name = "pressure" ;
;		pressure:units = "mbar" ;
;		pressure:positive = "down" ;
;	float rep_ion(time, pressure) ;
;		Prod:long_name = "REP ionpair production rate" ;
;		Prod:units = "/cm3/sec" ;
;	float time(time) ;
;		time:long_name = "time" ;
;		time:units = "days since 1900-01-01 00:00:00" ;
;		time:calendar = "gregorian" ;
;	int date(time) ;
;		date:long_name = "current date (YYYYMMDD)" ;
;	int datesec(time) ;
;		datesec:long_name = "current seconds of current date" ;
;	int day(time) ;
; -------------------------------------------

; ------------------------------------------
; Convert julian to WACCM time
; ------------------------------------------

timeint = dimsizes(rep_time)
npresslev = dimsizes(rep_plev)
nglat = dimsizes(rep_glat)

plev = new((/npresslev/),float)
date = new((/timeint/),integer)
time = new((/timeint/),double)
date = new((/timeint/),integer)
datestring = new((/timeint/),string)
datesec = new((/timeint/),integer)
day = new((/timeint/),integer)
day_all = new((/timeint/),integer)
daily_day = new((/timeint/),integer)
datej2g = new((/timeint,4/),integer)
dateyear= new((/timeint/),string)
year = new((/timeint/),integer)
datemonth = new((/timeint/),string)
month = new((/timeint/),integer)
dateday = new((/timeint/),string)
dayinteger = new((/timeint/),integer)
dateround = new((/timeint/),double)
datesecfrac = new((/timeint/),double)
datesecreal = new((/timeint/),double)
datesec = new((/timeint/),integer)
datesec_all = new((/timeint/),integer)
rep_ion = new((/timeint,npresslev,nglat/),float)
julian = new((/timeint/),double)
julian_all = new((/timeint/),double)
daily_julian = new((/timeint/),double)
timewaccm = new((/timeint/),double)

timewaccm(:) = rep_time(:)
print("timewaccm "+timewaccm)
;julian(:) = rep_time(:)+2415020.5
julian(:) = rep_time(:)+2396758.5
time(:) = timewaccm(:) 
print("time "+time)

plev = rep_plev
glat = rep_glat


print("glat"+ glat)

year@calendar = "gregorian"

datej2g = jul2greg(julian)

;printVarSummary(datej2g)
dateyear(:)= sprinti("%0.4i",datej2g(:,0))
datemonth(:) = sprinti("%0.2i",datej2g(:,1))
dateday(:) = sprinti("%0.2i",datej2g(:,2))
datestring(:) = dateyear(:)+datemonth(:)+dateday(:)
date = stringtointeger(datestring)
year = stringtointeger(dateyear)
month = stringtointeger(datemonth)
dayinteger = stringtointeger(dateday)


day(:) = day_of_year(year(:),month(:),dayinteger(:))

dateround = floor(julian)
datesecfrac = julian-dateround
datesecreal = datesecfrac*86400.
datesectemp = doubletointeger(datesecreal)  ; truncated 
   jb=0
   do while(jb.lt.timeint)
       if(datesectemp(jb).ge.43200.and.datesectemp(jb).lt.86400) then 
            datesec(jb) = datesectemp(jb)-43200
       else
            datesec(jb) =datesectemp(jb)+43200 
       end if
       jb=jb+1
   end do 

; ---------------------------------------------
; Define netCDF variables
; ---------------------------------------------

plev!0 = "plev"
plev&plev = plev
plev@long_name = "Pressure level"
plev@axis = "Z"
plev@_CoordinateAxisType = "Pressure"
;plev@units = "hPa"
plev@units = "mbar"
plev@positive = "down"

glat!0 = "glat"
glat&glat = glat 
glat@long_name = "geomagnetic latitude"
glat@axis = "Y"
glat@_CoordinateAxisType = "Lat"
glat@units = "degrees_north"

; days since 1900-01-01 00:00:00 -- Julian date 2415020.5
; days since 1850-01-01 00:00:00 -- Julian date 2396758.5
time!0 = "time"
time&time = time
time@axis = "T"
time@long_name = "time"
time@standard_name = "time"
time@units = "days since 1850-01-01 00:00:00"
time@calendar = "gregorian"
time@_CoordinateAxisType = "Time"
;print(time) 

date!0 = "time"
date&time = time
date@long_name = "current date (YYYYMMDD)"

datesec!0 = "time"
datesec&time = time
datesec@long_name = "current seconds of current date"

;day!0 = "time"
;day&time = time

;rep_ion = rep_ion_orig
mee_ion_rates = rep_ion_orig

mee_ion_rates!0 = "time"
mee_ion_rates&time = time
mee_ion_rates@cell_methods = "time: mean"
mee_ion_rates@long_name = "Ion pair production rate by solar protons"
mee_ion_rates@units = "/cm3/sec"

mee_ion_rates!1 = "plev"
mee_ion_rates&plev = plev

mee_ion_rates!2 = "glat"
mee_ion_rates&glat = glat 

julian!0 = "time"
julian&time = time

; -------------------------------------------------

system("/bin/rm -f "+filenameWACCM)
ncdf = addfile(filenameWACCM,"c")

filedimdef(ncdf,"time",-1,True)  ; time unlimited

;ncdf->plev = plev
;ncdf->glat = glat
;ncdf->julian = julian
;ncdf->time = time
ncdf->mee_ion_rates = mee_ion_rates 
ncdf->date = date
ncdf->datesec = datesec
;ncdf->day = day


end


