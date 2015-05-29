PRO latitude_average_and_remove_120days,input_file,varname,ndays_per_year,output_file,ensemble=ensemble,missing_value=missing_value,$
                                        analysis_file=analysis_file,analysis_offset=analysis_offset,forecast_offset=forecast_offset,$
                                        ndays_to_remove=ndays_to_remove,analysis_varname=analysis_varname,ntimes_per_day=ntimes_per_day,$
                                        analysis_mult=analysis_mult

; Procedure to latitude average between 15S-15N and remove the
; previous 120 days, in accordance with WH04 MJO indices.
;
; Specifying the /ENSEMBLE keyword removes the first 120 days from
; *every year*, to allow this procedure to be used with a single file
; containing a timeseries from each ensemble member.
 
;longitude=OPEN_AND_EXTRACT(input_file,'longitude')
;latitude=OPEN_AND_EXTRACT(input_file,'latitude')
;box=[-15.1,MIN(longitude),15.1,MAX(longitude)]
;DEFINE_BOUNDARIES,box,latitude,longitude,box_tx,/LIMIT
;n_lon=N_ELEMENTS(longitude)
;n_lat=N_ELEMENTS(latitude)

  IF KEYWORD_SET(analysis_varname) THEN BEGIN
     our_analysis_varname=analysis_varname
  ENDIF ELSE $
     our_analysis_varname=varname
  print,our_analysis_varname

  IF KEYWORD_SET(ntimes_per_day) THEN BEGIN
     our_ntpd=ntimes_per_day
  ENDIF ELSE $
     our_ntpd=1
  
  IF KEYWORD_SET(analysis_mult) THEN BEGIN
     our_analysis_mult=analysis_mult
  ENDIF ELSE $
     our_analysis_mult=1

ncid_in=NCDF_OPEN(input_file)
varid=NCDF_VARID(ncid_in,varname)
varstruct=NCDF_VARINQ(ncid_in,varid)

NCDF_DIMINQ,ncid_in,varstruct.dim(0),longitude_name,n_lon
NCDF_DIMINQ,ncid_in,varstruct.dim(1),latitude_name,n_lat
longitude=OPEN_AND_EXTRACT(input_file,longitude_name)
latitude=OPEN_AND_EXTRACT(input_file,latitude_name)
box=[-15.1,MIN(longitude),15.1,MAX(longitude)]
DEFINE_BOUNDARIES,box,latitude,longitude,box_tx,/LIMIT
n_lon=N_ELEMENTS(longitude)
n_lat=N_ELEMENTS(latitude)

IF varstruct.ndims eq 3 THEN BEGIN
   NCDF_DIMINQ,ncid_in,varstruct.dim(2),time_name,n_time
   time=OPEN_AND_EXTRACT(input_file,time_name)
   n_year=time/ndays_per_year/our_ntpd
ENDIF
IF varstruct.ndims eq 4 THEN BEGIN
   NCDF_DIMINQ,ncid_in,varstruct.dim(3),time_name,n_time
   NCDF_DIMINQ,ncid_in,varstruct.dim(2),year_name,n_year
   time=OPEN_AND_EXTRACT(input_file,time_name)
   year=OPEN_AND_EXTRACT(input_file,year_name)
ENDIF

IF KEYWORD_SET(analysis_file) THEN BEGIN
   ncid_analysis=NCDF_OPEN(analysis_file)
   varid_analysis=NCDF_VARID(ncid_analysis,our_analysis_varname)
   varstruct_analysis=NCDF_VARINQ(ncid_analysis,varid_analysis)
   IF varstruct_analysis.ndims eq 3 THEN BEGIN
      NCDF_DIMINQ,ncid_analysis,varstruct_analysis.dim(2),time_name_analysis,n_time_analysis
      time_analysis=OPEN_AND_EXTRACT(analysis_file,time_name_analysis)
      n_year_analysis=time_analysis/ndays_per_year
   ENDIF
   IF varstruct_analysis.ndims eq 4 THEN BEGIN
      NCDF_DIMINQ,ncid_analysis,varstruct_analysis.dim(3),time_name_analysis,n_time_analysis
      NCDF_DIMINQ,ncid_analysis,varstruct_analysis.dim(2),year_name_analysis,n_year_analysis
      time_analysis=OPEN_AND_EXTRACT(analysis_file,time_name_analysis)
      year_analysis=OPEN_AND_EXTRACT(analysis_file,year_name_analysis)
   ENDIF
ENDIF

IF KEYWORD_SET(ndays_to_remove) THEN BEGIN
   our_ndays_remove=ndays_to_remove
ENDIF ELSE $
   our_ndays_remove=120
   
weights=fltarr(n_lat)
FOR i=0,n_lat-1 DO $
   weights(i)=COS(latitude(i)*3.14159/180.)
weights=weights/TOTAL(weights)

IF varstruct.ndims eq 4 THEN BEGIN
   variable=OPEN_AND_EXTRACT(input_file,varname,$
                             offset=[box_tx(1),box_tx(0),0,0],$
                             count=[n_lon,n_lat,n_year,ndays_per_year*our_ntpd])
   IF TOTAL(where(ABS(variable) gt 1000)) ge 0 THEN $
      variable[where(ABS(variable) gt 1000)]=!Values.F_NaN
   IF KEYWORD_SET(missing_value) THEN BEGIN
      IF TOTAL(where(ABS(variable) eq missing_value)) ge 0 THEN $
         variable[where(ABS(variable) eq missing_value)]=!Values.F_NaN
   ENDIF
   variable_minus_120d=fltarr(n_lon,n_lat,n_year,ndays_per_year*our_ntpd)
   IF KEYWORD_SET(ensemble) THEN BEGIN
      variable_minus_120d=variable
   ENDIF ELSE BEGIN
      FOR i=0,n_lon-1 DO BEGIN
         FOR j=0,n_lat-1 DO BEGIN
            FOR k=0,n_year-1 DO BEGIN
               IF k eq 0 THEN BEGIN
                  year_start=our_ndays_remove
                  variable_minus_120d(i,j,k,0:our_ndays_remove*our_ntpd-1)=!Values.F_NaN
               ENDIF ELSE $
                  year_start=0            
               FOR m=year_start*our_ntpd,ndays_per_year*our_ntpd-1 DO BEGIN
                  last_120d=fltarr(our_ndays_remove*our_ntpd)
                  IF m ge 120*our_ntpd THEN BEGIN
;                     print,m
                     last_120d=variable(i,j,k,m-our_ndays_remove*our_ntpd:m-our_ntpd)
                  ENDIF ELSE IF m gt 0 THEN BEGIN
                     last_120d[our_ndays_remove*our_ntpd-m:our_ndays_remove*our_ntpd-1]=variable(i,j,k,0:m-1)
                     last_120d[0:(our_ndays_remove*our_ntpd-m-1)]=variable(i,j,k-1,ndays_per_year*our_ntpd-(our_ndays_remove*our_ntpd-m):ndays_per_year*our_ntpd-1)
                  ENDIF
                  IF TOTAL(where(ABS(last_120d) ge 10000)) gt 0 THEN BEGIN
                     bad=where(ABS(last_120d) ge 10000)
                     IF N_ELEMENTS(bad) eq 120 THEN BEGIN
                        mean_last_120d=0.
                     ENDIF ELSE BEGIN
                        good=where(last_120d lt 10000)
                        mean_last_120d=MEAN(last_120d[good])
                     ENDELSE
                  ENDIF ELSE $
                     mean_last_120d=MEAN(last_120d)                     
                  variable_minus_120d(i,j,k,m)=variable(i,j,k,m)-mean_last_120d
                  IF ABS(variable_minus_120d(i,j,k,m)) ge 10000 THEN $
                     variable_minus_120d(i,j,k,m)=!Values.F_NaN
               ENDFOR
            ENDFOR
         ENDFOR
      ENDFOR  
   ENDELSE               
   variable_latavg=fltarr(n_lon,n_year,ndays_per_year*our_ntpd)
   FOR i=0,n_lon-1 DO $
      FOR j=0,n_year-1 DO $
         FOR k=0,ndays_per_year*our_ntpd-1 DO $
            variable_latavg(i,j,k)=TOTAL(REFORM(variable_minus_120d(i,*,j,k))*weights,/NaN)
ENDIF ELSE IF varstruct.ndims eq 3 THEN BEGIN
   variable=OPEN_AND_EXTRACT(input_file,varname,$
                             offset=[box_tx(1),box_tx(0),0],$
                             count=[n_lon,n_lat,n_time])
   IF TOTAL(where(ABS(variable) gt 1000)) ge 0 THEN $
      variable[where(ABS(variable) gt 1000)]=!Values.F_NaN
   IF KEYWORD_SET(missing_value) THEN BEGIN
      IF TOTAL(where(ABS(variable) ge ABS(missing_value))) ge 0 THEN $
         variable[where(ABS(variable) ge ABS(missing_value))]=!Values.F_NaN      
   ENDIF
   variable_minus_120d=fltarr(n_lon,n_lat,n_time)   
   IF KEYWORD_SET(analysis_file) THEN BEGIN
      IF KEYWORD_SET(analysis_offset) THEN BEGIN
         our_analysis_offset=analysis_offset
      ENDIF ELSE $
         our_analysis_offset=0
      IF KEYWORD_SET(forecast_offset) THEN BEGIN
         our_forecast_offset=forecast_offset
      ENDIF ELSE $
         our_forecast_offset=0      
      variable_analysis=OPEN_AND_EXTRACT(analysis_file,our_analysis_varname,$
                                         offset=[box_tx(1),box_tx(0),our_analysis_offset],$
                                         count=[n_lon,n_lat,our_ndays_remove])*our_analysis_mult
      IF TOTAL(where(ABS(variable_analysis) ge 10000)) gt 0 THEN $
         variable_analysis[where(ABS(variable_analysis) ge 10000)]=!Values.F_NaN
      FOR i=0,n_lon-1 DO BEGIN
         FOR j=0,n_lat-1 DO BEGIN
            mean_last_120d=fltarr(n_time)
            FOR k=0,n_time-1 DO BEGIN
               ;print,k
               IF k ge our_ntpd and k lt our_ndays_remove*our_ntpd THEN BEGIN
;                  IF i eq 0 and j eq 0 THEN $
;                     print,'analysis contributes '+STRTRIM(STRING(our_ndays_remove-k/our_ntpd),1)+' forecast contributes: '+STRTRIM(STRING(k/our_ntpd*our_ntpd),1)
                  mean_last_120d(k)=MEAN(REFORM(variable_analysis(i,j,k/our_ntpd:our_ndays_remove-1)),/NaN)*((our_ndays_remove-k/our_ntpd)/FLOAT(our_ndays_remove))+$
                                    MEAN(REFORM(variable(i,j,0:k/our_ntpd*our_ntpd-1)),/NaN)*((k/our_ntpd)/FLOAT(our_ndays_remove))
               ENDIF ELSE IF k le our_ntpd THEN BEGIN
                  mean_last_120d(k)=MEAN(variable_analysis(i,j,k/our_ntpd:our_ndays_remove-1),/NaN)
               ENDIF ELSE IF k ge our_ndays_remove*our_ntpd THEN $
                  mean_last_120d(k)=MEAN(variable(i,j,k/our_ntpd*our_ntpd-our_ndays_remove:k/our_ntpd*our_ntpd-1),/NaN)
            ENDFOR
            IF our_forecast_offset gt 0 THEN BEGIN
               variable_minus_120d(i,j,0:forecast_offset-1)=!Values.F_NaN
               variable_minus_120d(i,j,forecast_offset:n_time-1)=variable(i,j,forecast_offset:n_time-1)-mean_last_120d(*)
            ENDIF ELSE $
               variable_minus_120d(i,j,*)=variable(i,j,*)-mean_last_120d(*)
         ENDFOR
      ENDFOR
   ENDIF ELSE BEGIN
      FOR i=0,n_lon-1 DO BEGIN
         FOR j=0,n_lat-1 DO BEGIN
            variable_minus_120d(i,j,0:our_ndays_remove-1)=!Values.F_NaN
            FOR k=our_ndays_remove,n_time-1 DO BEGIN
               last_120d=REFORM(variable(i,j,k-our_ndays_remove:k-1))
               IF TOTAL(where(ABS(last_120d) ge 10000)) gt 0 THEN BEGIN
                  bad=where(ABS(last_120d) ge 10000)
                  IF N_ELEMENTS(bad) eq our_ndays_remove THEN BEGIN
                     mean_last_120d=0.
                  ENDIF ELSE BEGIN
                     good=where(ABS(last_120d) lt 10000)
                     mean_last_120d=MEAN(last_120d[good])
                  ENDELSE
               ENDIF ELSE $
                  mean_last_120d=MEAN(last_120d,/NaN)
               variable_minus_120d(i,j,k)=variable(i,j,k)-mean_last_120d
            ENDFOR
         ENDFOR
      ENDFOR
   ENDELSE
   variable_latavg=fltarr(n_lon,n_time)
   FOR i=0,n_lon-1 DO $
      FOR j=0,n_time-1 DO $
         variable_latavg(i,j)=TOTAL(REFORM(variable_minus_120d(i,*,j))*weights,/NaN)
ENDIF

out_id=NCDF_CREATE(output_file,/CLOBBER)
out_dimids=intarr(varstruct.ndims)
out_varids=intarr(varstruct.ndims+1)
out_dimids(0)=NCDF_DIMDEF(out_id,'longitude',n_lon)
out_dimids(2)=NCDF_DIMDEF(out_id,time_name,n_time)
IF varstruct.ndims eq 4 THEN $
   out_dimids(3)=NCDF_DIMDEF(out_id,year_name,n_year)
out_varids(0)=NCDF_VARDEF(out_id,'longitude',[out_dimids(0)])
out_varids(2)=NCDF_VARDEF(out_id,time_name,[out_dimids(2)])
IF varstruct.ndims eq 4 THEN BEGIN
   out_varids(3)=NCDF_VARDEF(out_id,year_name,[out_dimids(3)])
   out_varids(4)=NCDF_VARDEF(out_id,varname,[out_dimids(0),out_dimids(3),out_dimids(2)])
ENDIF ELSE $
   out_varids(3)=NCDF_VARDEF(out_id,varname,[out_dimids(0),out_dimids(2)])

NCDF_CONTROL,out_id,/ENDEF

NCDF_VARPUT,out_id,out_varids(0),longitude
NCDF_VARPUT,out_id,out_varids(2),time
IF varstruct.ndims eq 4 THEN BEGIN
   NCDF_VARPUT,out_id,out_varids(3),year
   NCDF_VARPUT,out_id,out_varids(4),variable_latavg
ENDIF ELSE $
   NCDF_VARPUT,out_id,out_varids(3),variable_latavg

NCDF_CLOSE,out_id

END
