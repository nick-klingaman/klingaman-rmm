PRO remove_anncycle_three_harmonics,input_file,clim_input_file,varname,ndays_per_year,output_file,missing_value=missing_value,$
                                    clim_varname=clim_varname,start_date=start_date,ntimes_per_day=ntimes_per_day,obs_mult=obs_mult,$
                                    output_anncycle=output_anncycle,output_file_anncycle=output_file_anncycle

; Generalized procedure to remove a climatology and the first three
; harmonics of the annual cycle from an input file
;
; NPK 12/2/10
;

;longitude=OPEN_AND_EXTRACT(input_file,'longitude')
;latitude=OPEN_AND_EXTRACT(input_file,'latitude')
;n_lon=N_ELEMENTS(longitude)
;n_lat=N_ELEMENTS(latitude)

IF KEYWORD_SET(ntimes_per_day) THEN BEGIN
   our_ntpd=ntimes_per_day
ENDIF ELSE $
   our_ntpd=1
IF KEYWORD_SET(obs_mult) THEN BEGIN
   our_obs_mult=obs_mult
ENDIF ELSE $
   our_obs_mult=1

ncid_in=NCDF_OPEN(input_file)
varid=NCDF_VARID(ncid_in,varname)
varstruct=NCDF_VARINQ(ncid_in,varid)

NCDF_DIMINQ,ncid_in,varstruct.dim(0),longitude_name,n_lon
NCDF_DIMINQ,ncid_in,varstruct.dim(1),latitude_name,n_lat
longitude=OPEN_AND_EXTRACT(input_file,longitude_name)
latitude=OPEN_AND_EXTRACT(input_file,latitude_name)
n_lon=N_ELEMENTS(longitude)
n_lat=N_ELEMENTS(latitude)

IF varstruct.ndims eq 3 THEN BEGIN
   NCDF_DIMINQ,ncid_in,varstruct.dim(2),time_name,n_time
   time=OPEN_AND_EXTRACT(input_file,time_name)
   n_time=N_ELEMENTS(time)
   n_year=n_time/ndays_per_year
   IF n_time/our_ntpd le ndays_per_year THEN $
      n_year=1
ENDIF
IF varstruct.ndims eq 4 THEN BEGIN
   NCDF_DIMINQ,ncid_in,varstruct.dim(3),time_name,n_time
   NCDF_DIMINQ,ncid_in,varstruct.dim(2),year_name,n_year
   time=OPEN_AND_EXTRACT(input_file,time_name)
   year=OPEN_AND_EXTRACT(input_file,year_name)
ENDIF

IF KEYWORD_SET(clim_varname) THEN BEGIN
   clim_varname=clim_varname
ENDIF ELSE $
   clim_varname=varname

IF varstruct.ndims eq 4 THEN BEGIN
   variable=OPEN_AND_EXTRACT(input_file,varname)
   variable_clim=OPEN_AND_EXTRACT(clim_input_file,clim_varname)*our_obs_mult
   variable_first_three_harm=fltarr(n_lon,n_lat,ndays_per_year)
   variable_anomalies=fltarr(n_lon,n_lat,n_year,ndays_per_year*our_ntpd)
   FOR i=0,n_lon-1 DO BEGIN
      FOR j=0,n_lat-1 DO BEGIN
         variable_this_gridpt=REFORM(variable(i,j,*,*))
         variable_clim_this_gridpt=REFORM(variable_clim(i,j,*))
;         IF i eq 30 and j eq 30 THEN $
;            STOP
         clim_fft_freq=FFT(variable_clim_this_gridpt,-1,DIMENSION=1,/DOUBLE)
         clim_fft_freq(4:N_ELEMENTS(clim_fft_freq)-4)=0.
         variable_first_three_harm(i,j,*)=FFT(clim_fft_freq,1,DIMENSION=1,/DOUBLE)
         FOR k=0,n_year-1 DO BEGIN
            FOR m=0,ndays_per_year-1 DO $
               variable_anomalies(i,j,k,m*our_ntpd:(m+1)*our_ntpd-1)=variable_this_gridpt(k,m*our_ntpd:(m+1)*our_ntpd-1)-$
               variable_first_three_harm(i,j,m)            
            IF KEYWORD_SET(missing_value) THEN BEGIN
               temp_this_gridpt=REFORM(variable_this_gridpt(k,*))
               temp_anomalies=REFORM(variable_anomalies(i,j,k,*))
               IF TOTAL(where(temp_this_gridpt eq missing_value)) gt 0 THEN $
                  temp_anomalies[where(temp_this_gridpt eq missing_value)]=missing_value
               variable_anomalies(i,j,k,*)=temp_anomalies
            ENDIF
         ENDFOR
      ENDFOR
   ENDFOR
ENDIF ELSE IF varstruct.ndims eq 3 THEN BEGIN
   variable=OPEN_AND_EXTRACT(input_file,varname)
   variable_clim=OPEN_AND_EXTRACT(clim_input_file,clim_varname)*our_obs_mult
   variable_first_three_harm=fltarr(n_lon,n_lat,ndays_per_year)
   variable_anomalies=fltarr(n_lon,n_lat,n_time)
   IF n_year eq 1 THEN BEGIN
      variable_this_gridpt=fltarr(n_time)
   ENDIF ELSE $
      variable_this_gridpt=fltarr(n_year,n_time)
   FOR i=0,n_lon-1 DO BEGIN
      FOR j=0,n_lat-1 DO BEGIN
         IF n_year eq 1 THEN BEGIN
            variable_this_gridpt(*)=REFORM(variable(i,j,*))
         ENDIF ELSE BEGIN
            FOR k=0,n_year-1 DO $
               variable_this_gridpt(k,*)=REFORM(variable(i,j,k*ndays_per_year*our_ntpd:(k+1)*ndays_per_year*our_ntpd-1))
         ENDELSE
         variable_clim_this_gridpt=REFORM(variable_clim(i,j,*))
         clim_fft_freq=FFT(variable_clim_this_gridpt,-1,DIMENSION=1,/DOUBLE)
         clim_fft_freq(4:N_ELEMENTS(clim_fft_freq)-4)=0.
         variable_first_three_harm(i,j,*)=FFT(clim_fft_freq,1,DIMENSION=1,/DOUBLE)
         IF n_year eq 1 THEN BEGIN
            IF KEYWORD_SET(start_date) THEN BEGIN
               IF start_date+(n_time/our_ntpd) lt ndays_per_year THEN BEGIN                  
                  FOR m=0,n_time-1 DO BEGIN
                     variable_anomalies(i,j,m)=variable_this_gridpt(m)-$
                        variable_first_three_harm(i,j,start_date+m/our_ntpd)
                     ;IF longitude(i) eq 91.25 and latitude(j) eq 0.0 THEN $
                     ;   print,start_date+m/our_ntpd,variable_first_three_harm(i,j,start_date+m/our_ntpd)
                  ENDFOR
               ENDIF ELSE BEGIN
                  variable_first_three_harm_cross=fltarr(n_time)
                  FOR m=0,(ndays_per_year-start_date)*our_ntpd-1 DO BEGIN
                     ;print,m,start_date+m/our_ntpd
                     variable_first_three_harm_cross(m)=variable_first_three_harm(i,j,start_date+m/our_ntpd)
                  ENDFOR
                  FOR m=(ndays_per_year-start_date)*our_ntpd,n_time-1 DO BEGIN
                     ;print,m,m/our_ntpd-(ndays_per_year-start_date)
                     variable_first_three_harm_cross(m)=variable_first_three_harm(i,j,m/our_ntpd-(ndays_per_year-start_date))
                  ENDFOR
                  variable_anomalies(i,j,*)=variable_this_gridpt(*)-variable_first_three_harm_cross(*)
               ENDELSE
            ENDIF ELSE BEGIN
               FOR m=0,n_time-1 DO $
                  variable_anomalies(i,j,m)=variable_this_gridpt(m)-variable_first_three_harm(i,j,m/our_ntpd)
            ENDELSE
         ENDIF ELSE BEGIN
            FOR k=0,n_year-1 DO $               
               variable_anomalies(i,j,k*ndays_per_year:(k+1)*ndays_per_year-1)=variable_this_gridpt(k,*)-variable_first_three_harm(i,j,*)
         ENDELSE
            IF KEYWORD_SET(missing_value) THEN BEGIN
               temp_this_gridpt=REFORM(variable_this_gridpt(k,*))
               temp_anomalies=REFORM(variable_anomalies(i,j,k,*))
               IF TOTAL(where(temp_this_gridpt eq missing_value)) gt 0 THEN $
                  temp_anomalies[where(temp_this_gridpt eq missing_value)]=missing_value
               variable_anomalies(i,j,k,*)=temp_anomalies
            ENDIF
         ENDFOR
   ENDFOR
ENDIF

out_id=NCDF_CREATE(output_file,/CLOBBER)
out_dimids=intarr(varstruct.ndims)
out_varids=intarr(varstruct.ndims+1)
out_dimids(0)=NCDF_DIMDEF(out_id,'longitude',n_lon)
out_dimids(1)=NCDF_DIMDEF(out_id,'latitude',n_lat)
out_dimids(2)=NCDF_DIMDEF(out_id,time_name,n_time)
IF varstruct.ndims eq 4 THEN $
   out_dimids(3)=NCDF_DIMDEF(out_id,year_name,n_year)
out_varids(0)=NCDF_VARDEF(out_id,'longitude',[out_dimids(0)])
out_varids(1)=NCDF_VARDEF(out_id,'latitude',[out_dimids(1)])
out_varids(2)=NCDF_VARDEF(out_id,time_name,[out_dimids(2)])
IF varstruct.ndims eq 4 THEN BEGIN
   out_varids(3)=NCDF_VARDEF(out_id,year_name,[out_dimids(3)])
   out_varids(4)=NCDF_VARDEF(out_id,varname,[out_dimids(0),out_dimids(1),out_dimids(3),out_dimids(2)])
ENDIF ELSE $
   out_varids(3)=NCDF_VARDEF(out_id,varname,[out_dimids(0),out_dimids(1),out_dimids(2)])

IF KEYWORD_SET(missing_value) THEN BEGIN
   IF varstruct.ndims eq 4 THEN BEGIN
      NCDF_ATTPUT,out_id,out_varids(4),'missing_value',missing_value
   ENDIF ELSE $
      NCDF_ATTPUT,out_id,out_varids(3),'missing_value',missing_value
ENDIF
NCDF_CONTROL,out_id,/ENDEF

NCDF_VARPUT,out_id,out_varids(0),longitude
NCDF_VARPUT,out_id,out_varids(1),latitude
NCDF_VARPUT,out_id,out_varids(2),time
IF varstruct.ndims eq 4 THEN BEGIN
   NCDF_VARPUT,out_id,out_varids(3),year
   NCDF_VARPUT,out_id,out_varids(4),variable_anomalies
ENDIF ELSE $
   NCDF_VARPUT,out_id,out_varids(3),variable_anomalies

NCDF_CLOSE,out_id

IF KEYWORD_SET(output_anncycle) THEN BEGIN
   out_id=NCDF_CREATE(output_file_anncycle,/CLOBBER)
   out_dimids=intarr(3)
   out_varids=intarr(4)
   out_dimids(0)=NCDF_DIMDEF(out_id,'longitude',n_lon)
   out_dimids(1)=NCDF_DIMDEF(out_id,'latitude',n_lat)
   out_dimids(2)=NCDF_DIMDEF(out_id,time_name,ndays_per_year*our_ntpd)
   out_varids(0)=NCDF_VARDEF(out_id,'longitude',[out_dimids(0)])
   out_varids(1)=NCDF_VARDEF(out_id,'latitude',[out_dimids(1)])
   out_varids(2)=NCDF_VARDEF(out_id,time_name,[out_dimids(2)])
   out_varids(3)=NCDF_VARDEF(out_id,varname,[out_dimids(0),out_dimids(1),out_dimids(2)])
   NCDF_CONTROL,out_id,/ENDEF

   NCDF_VARPUT,out_id,out_varids(0),longitude
   NCDF_VARPUT,out_id,out_varids(1),latitude
   NCDF_VARPUT,out_id,out_varids(2),findgen(ndays_per_year*our_ntpd)/our_ntpd+our_ntpd/2.
   NCDF_VARPUT,out_id,out_varids(3),variable_first_three_harm

   NCDF_CLOSE,out_id
ENDIF   

;STOP

END
