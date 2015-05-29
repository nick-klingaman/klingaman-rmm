PRO lanczos_filter_file,infile,varname,upper_limit,lower_limit,outfile,remove_mean=remove_mean,continuous=continuous,field_multiplier=field_multiplier

; Accept an input file, a variable name, an upper and lower limit for the lanczos
; bandpass filter, and an output file.

; Filter the input file at every gridpoint, then output the result.

; The resulting output file will have a number of time points equal to
; the number of time points in the input file minus the lower limit
; divided by two.
;
; out_time = in_time-(lower_limit/2)
;
; You can use the REMOVE_MEAN keyword to remove the ensemble-mean
; daily climatology.  This works for four-dimensional (lon x lat x
; members x time) files ONLY.
;
; NPK 24/2/09

; Protect inputs
my_infile=infile
my_varname=varname
my_upper_limit=upper_limit
my_lower_limit=lower_limit
my_outfile=outfile

id=NCDF_OPEN(my_infile,/NOWRITE)
varid=NCDF_VARID(id,my_varname)
var_struct=NCDF_VARINQ(id,varid)
nlon=0
nlat=0
nmem=0
ntime=0
FOR i=0,var_struct.ndims-1 DO BEGIN
    NCDF_DIMINQ,id,var_struct.dim(i),dim_name,dim_size
    IF dim_name EQ 'longitude' OR dim_name EQ 'lon' OR dim_name EQ 'x' THEN BEGIN
       lon_dim=var_struct.dim(i)
       lon_name=dim_name
       nlon=dim_size
    ENDIF ELSE IF dim_name EQ 'latitude' OR dim_name EQ 'lat' OR dim_name EQ 'y' THEN BEGIN
       lat_dim=var_struct.dim(i)
       lat_name=dim_name
       nlat=dim_size
    ENDIF ELSE IF dim_name EQ 'time' OR dim_name EQ 't' THEN BEGIN
       time_dim=var_struct.dim(i)
       time_name=dim_name
       ntime=dim_size
    ENDIF ELSE IF dim_name EQ 'member' OR dim_name EQ 'record' OR dim_name EQ 'year' THEN BEGIN
       mem_dim=var_struct.dim(i)
       mem_name=dim_name
       nmem=dim_size
    ENDIF
 ENDFOR
IF (nlat EQ 0) THEN BEGIN
   print,'LANCZOS_FILTER_FILE : Could not find latitude dimension (tried latitude,lat,y).'
   STOP
ENDIF ELSE IF (nlon EQ 0) THEN BEGIN
    print,'LANCZOS_FILTER_FILE : Could not find longitude dimension (tried longitude,lon,x).'
    STOP
ENDIF ELSE IF (ntime EQ 0) THEN BEGIN
    print,'LANCZOS_FILTER_FILE : Could not find time dimension (tried time,t).'
    STOP
ENDIF ELSE IF (ntime LE lower_limit) THEN BEGIN
    print,'LANCZOS_FILTER_FILE : There are not enough time points in the input file to perform the desired filter.'
    STOP
ENDIF
longitude=OPEN_AND_EXTRACT(my_infile,lon_name)
latitude=OPEN_AND_EXTRACT(my_infile,lat_name)
time=OPEN_AND_EXTRACT(my_infile,time_name)

; Create output file
outid=NCDF_CREATE(my_outfile,/CLOBBER)
dimids=intarr(4)
dimids(0)=NCDF_DIMDEF(outid,'longitude',nlon)
dimids(1)=NCDF_DIMDEF(outid,'latitude',nlat)
IF (nmem GT 0) THEN $
  dimids(3)=NCDF_DIMDEF(outid,'member',nmem)
IF KEYWORD_SET(continuous) and nmem gt 0 THEN BEGIN
   dimids(2)=NCDF_DIMDEF(outid,'time',ntime)
ENDIF ELSE $
   dimids(2)=NCDF_DIMDEF(outid,'time',ntime-lower_limit)
varids=intarr(5)
varids(0)=NCDF_VARDEF(outid,'longitude',[dimids(0)])
varids(1)=NCDF_VARDEF(outid,'latitude',[dimids(1)])
IF (nmem GT 0) THEN $
  varids(2)=NCDF_VARDEF(outid,'member',[dimids(3)])
varids(3)=NCDF_VARDEF(outid,'time',[dimids(2)])
IF (nmem GT 0) THEN BEGIN
   varids(4)=NCDF_VARDEF(outid,varname,[dimids(2),dimids(3),dimids(1),dimids(0)])
ENDIF ELSE $
   varids(4)=NCDF_VARDEF(outid,varname,[dimids(0),dimids(1),dimids(2)])

NCDF_ATTPUT,outid,varids(4),'missing_value',2E20

NCDF_CONTROL,outid,/ENDEF

NCDF_VARPUT,outid,varids(0),longitude
NCDF_VARPUT,outid,varids(1),latitude
IF (nmem GT 0) THEN $
  NCDF_VARPUT,outid,varids(2),indgen(nmem)
IF KEYWORD_SET(continuous) and nmem gt 0 THEN BEGIN
   NCDF_VARPUT,outid,varids(3),time
ENDIF ELSE $
   NCDF_VARPUT,outid,varids(3),time(lower_limit/2:ntime-lower_limit/2-1)

; Read data
IF KEYWORD_SET(field_multiplier) THEN BEGIN
   our_field_multiplier=field_multiplier
ENDIF ELSE $
   our_field_multiplier=1.
variable=OPEN_AND_EXTRACT(my_infile,my_varname)*our_field_multiplier
variable_filtered=fltarr(nlon,nlat,ntime-lower_limit)
;IF KEYWORD_SET(continuous) and nmem GT 0 THEN $
;   variable_filtered=fltarr(nlon,nlat,ntime*nmem)

IF TOTAL(where(variable ge 1E20)) gt 0 THEN $
  variable[where(variable ge 1E20)]=!Values.F_NaN

IF KEYWORD_SET(remove_mean) AND nmem GT 0 THEN BEGIN
    print,'Now removing ensemble mean ...'
    FOR i=0,nlon-1 DO BEGIN
        FOR j=0,nlat-1 DO BEGIN
            FOR k=0,ntime-1 DO BEGIN
               variable_ensmean=MEAN(variable(i,j,k,*),/NaN)
                IF (variable_ensmean gt 1E10) THEN BEGIN
                    print,'Suspect missing data: i=',i,'j=',j,'k=',k
                    STOP
                ENDIF                
                variable(i,j,k,*)=variable(i,j,k,*)-variable_ensmean
            ENDFOR
        ENDFOR
    ENDFOR
ENDIF
                
; Perform filtering
IF (nmem GT 0) THEN BEGIN
   IF KEYWORD_SET(continuous) THEN BEGIN
      FOR i=0,nlon-1 DO BEGIN
         print,'---> Filtering for longitude band '+STRTRIM(STRING(i+1),1)+' of '+STRTRIM(STRING(nlon),1)+' ...'
         FOR j=0,nlat-1 DO BEGIN
            ;print,'-------> Filtering for latitude '+STRTRIM(STRING(j+1),1)+' of '+STRTRIM(STRING(nlat),1)+' ...'
            this_point_ts=fltarr(ntime*nmem)
            this_point_ts_filtered=fltarr(ntime*nmem)
            FOR m=0,nmem-1 DO BEGIN
               time_start=m*ntime
               time_end=(m+1)*ntime-1
               this_point_ts(time_start:time_end)=variable(i,j,m,*)
            ENDFOR
            IF N_ELEMENTS(where(FINITE(this_point_ts) eq 1)) le lower_limit THEN BEGIN
               this_point_ts_filtered(*)=2E20
            ENDIF ELSE IF TOTAL(this_point_ts) eq 0 and this_point_ts(0) eq 0 and this_point_ts(ntime-1) eq 0 THEN BEGIN    
               this_point_ts_filtered(*)=2E20
               print,'Suspect no data (total=0) at '+STRTRIM(STRING(i),1)+STRTRIM(STRING(j),1)
;               STOP
            ENDIF ELSE BEGIN
               ;print,'---> Filtering for latitude point '+STRTRIM(STRING(j+1),1)+' of '+STRTRIM(STRING(nlat),1)+' ...'
               good=where(FINITE(this_point_ts) eq 1)
               bad=where(FINITE(this_point_ts) eq 0)               
               this_point_good=this_point_ts[good]
                                ; this_point_good=this_point_good-MEAN(this_point_good)
               ;print,'-------> Filtering ...'
               this_point_good=LANCZOS_BANDPASS(this_point_good,1./FLOAT(lower_limit),1./FLOAT(upper_limit),$
                                                2*lower_limit+1,/DETREND)
               ;print,'-------> .... complete.'
               this_point_ts[good]=this_point_good
               IF TOTAL(bad) gt 0 THEN $
                  this_point_ts[bad]=2E20
               this_point_ts_filtered(*)=this_point_ts
            ENDELSE
            ;print,'-------> Writing to disk ...'
            FOR m=0,nmem-1 DO BEGIN
               put_to_netcdf=fltarr(ntime)
               IF m eq 0 THEN BEGIN
                  time_start=m*ntime+lower_limit/2
                  put_start=lower_limit/2
                  put_to_netcdf(0:lower_limit/2)=2E20
               ENDIF ELSE BEGIN
                  time_start=m*ntime
                  put_start=0
               ENDELSE
               IF m eq nmem-1 THEN BEGIN
                  time_end=(m+1)*ntime-lower_limit/2-1
                  put_end=ntime-lower_limit/2-1
                  put_to_netcdf(lower_limit/2:ntime-1)=2E20
               ENDIF ELSE BEGIN
                  time_end=(m+1)*ntime-1
                  put_end=ntime-1
               ENDELSE
               put_to_netcdf(put_start:put_end)=this_point_ts_filtered(time_start:time_end)
               NCDF_VARPUT,outid,varids(4),put_to_netcdf,offset=[0,m,j,i],$
                           count=[ntime,1,1,1]
            ENDFOR
            ;print,'--------> ... complete.'
         ENDFOR
      ENDFOR
   ENDIF ELSE BEGIN
      FOR k=0,nmem-1 DO BEGIN
         print,'Now filtering for member '+STRTRIM(STRING(k+1),1)+' ...'
         FOR i=0,nlon-1 DO BEGIN
            FOR j=0,nlat-1 DO BEGIN
               this_point=REFORM(variable(i,j,k,*))
               this_point_allmem=REFORM(variable(i,j,*,*))
;                print,i,j,N_ELEMENTS(where(FINITE(this_point) eq 1))
               IF N_ELEMENTS(where(FINITE(this_point) eq 1)) le lower_limit THEN BEGIN
                  variable_filtered(i,j,*)=2E20                  
               ENDIF ELSE BEGIN
                  good=where(FINITE(this_point) eq 1)
                  bad=where(FINITE(this_point) eq 0)
                  good_allmem=where(FINITE(this_point_allmem) eq 1)
                  this_point_good=this_point[good]
                                ;this_point_good=this_point_good-MEAN(this_point_allmem[good_allmem],/NaN)
                  this_point_good=LANCZOS_BANDPASS(this_point_good,1./FLOAT(lower_limit),1./FLOAT(upper_limit),$
                                                   2*lower_limit+1,/DETREND)           
                  this_point[good]=this_point_good
                  IF TOTAL(bad) gt 0 THEN $
                     this_point[bad]=2E20
                  variable_filtered(i,j,*)=this_point(lower_limit/2:ntime-lower_limit/2-1)                  
               ENDELSE
            ENDFOR
         ENDFOR
         NCDF_VARPUT,outid,varids(4),variable_filtered,offset=[0,0,k,0],$
                     count=[nlon,nlat,1,ntime-lower_limit]
      ENDFOR
   ENDELSE
ENDIF ELSE BEGIN
   FOR i=0,nlon-1 DO BEGIN
	print,'Now filtering for longitude band '+STRTRIM(STRING(i+1),1)+' of '+STRTRIM(STRING(nlon),1)
        FOR j=0,nlat-1 DO BEGIN
            this_point=REFORM(variable(i,j,*))
            ;print,i,j,N_ELEMENTS(where(FINITE(this_point) eq 1))
            IF N_ELEMENTS(where(FINITE(this_point) eq 1)) le lower_limit THEN BEGIN
                variable_filtered(i,j,*)=2E20
            ENDIF ELSE BEGIN
                good=where(FINITE(this_point) eq 1)
                bad=where(FINITE(this_point) eq 0)
                this_point_good=this_point[good]
                this_point_good=this_point_good-MEAN(this_point_good,/NaN)
                this_point_good=LANCZOS_BANDPASS(this_point_good,1./FLOAT(lower_limit),1./FLOAT(upper_limit),$
                                                 2*lower_limit+1,/DETREND)
                this_point[good]=this_point_good
                IF TOTAL(bad) gt 0 THEN $
                  this_point[bad]=2E20
                variable_filtered(i,j,*)=this_point(lower_limit/2:ntime-lower_limit/2-1)            
            ENDELSE
        ENDFOR
    ENDFOR
    NCDF_VARPUT,outid,varids(4),variable_filtered
ENDELSE

NCDF_CLOSE,outid

END
