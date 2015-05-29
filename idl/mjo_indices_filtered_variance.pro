PRO mjo_indices_filtered_variance,infile,runid,varname,upper_limit,lower_limit,n_years,n_days_per_year,box=box,mask_file=mask_file,mask_plot=mask_plot
  
; Procedure to plot filtered variance for a given variable in a netCDF
; file.  Use pre-defined contour ranges based on the variable name.

print,infile
ncid_in=NCDF_OPEN(infile)
varid=NCDF_VARID(ncid_in,varname)
varstruct=NCDF_VARINQ(ncid_in,varid)

CASE varname OF 
   'olr' : BEGIN
      mylevs=['2.0','4.0','6.0','8.0','10.0','12.0','14.0','16.0','20.0','24.0','28.0','32.0']
      units='W m!U-2!N'   
      multiplier=1.
   END
   'U' : BEGIN
      mylevs=['0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5']
      units='m s!U-1!N'
      multiplier=1.
   END
   'u850' : BEGIN
      mylevs=['0.25','0.5','1.0','1.25','1.50','1.75','2.00','2.25','2.50','2.75','3.0','3.25','3.5','3.75']
      units='m s!U-1!N' 
      multiplier=1.      
   END
   'u200' : BEGIN
      mylevs=['2.5','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.0','7.5','8.0','8.5','9.0']
      units='m s!U-1!N'
      multiplier=1.
   END
   'temp_1' : BEGIN
      mylevs=['0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.50','0.60','0.70','0.80']
      units='!Uo!NC'
      multiplier=1.
   END
   'sst' : BEGIN
      mylevs=['0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.50','0.60','0.70','0.80']
      units='!Uo!NC'
      multiplier=1.
   END     
   'precip' : BEGIN
      mylevs=['0.30','0.90','1.50','1.80','2.40','3.00','3.60','4.20','4.80','5.40','6.00','6.60']
      units='mm day!U-1!N'
      multiplier=1.
   END
ENDCASE

NCDF_DIMINQ,ncid_in,varstruct.dim(0),longitude_name,n_lon
NCDF_DIMINQ,ncid_in,varstruct.dim(1),latitude_name,n_lat
longitude=OPEN_AND_EXTRACT(infile,longitude_name)
latitude=OPEN_AND_EXTRACT(infile,latitude_name)
IF KEYWORD_SET(box) THEN $
   DEFINE_BOUNDARIES,box,latitude,longitude,box_tx,/LIMIT
n_lon=N_ELEMENTS(longitude)
n_lat=N_ELEMENTS(latitude)

IF varstruct.ndims eq 3 THEN BEGIN
   NCDF_DIMINQ,ncid_in,varstruct.dim(2),time_name,n_time
   time=OPEN_AND_EXTRACT(infile,time_name)
   n_year=time/ndays_per_year
ENDIF
IF varstruct.ndims eq 4 THEN BEGIN
   NCDF_DIMINQ,ncid_in,varstruct.dim(3),time_name,n_time
   NCDF_DIMINQ,ncid_in,varstruct.dim(2),year_name,n_year
   time=OPEN_AND_EXTRACT(infile,time_name)
   year=OPEN_AND_EXTRACT(infile,year_name)
ENDIF
NCDF_ATTGET,ncid_in,varid,'missing_value',missing

IF KEYWORD_SET(mask_plot) THEN BEGIN
   mask_longitude=OPEN_AND_EXTRACT(mask_file,'longitude')
   mask_latitude=OPEN_AND_EXTRACT(mask_file,'latitude')
   IF KEYWORD_SET(box) THEN $
      DEFINE_BOUNDARIES,box,mask_latitude,mask_longitude,mask_box_tx,/LIMIT
   mask_nlon=N_ELEMENTS(mask_longitude)
   mask_nlat=N_ELEMENTS(mask_latitude)
   lsm=REFORM(OPEN_AND_EXTRACT(mask_file,'lsm',$
                                offset=[mask_box_tx(1),mask_box_tx(0),0,0],$
                                count=[mask_nlon,mask_nlat,1,1]))
ENDIF

IF varstruct.ndims eq 4 THEN BEGIN
   IF KEYWORD_SET(box) THEN BEGIN
      variable=REFORM(OPEN_AND_EXTRACT(infile,varname,$
                                       offset=[box_tx(1),box_tx(0),0,0],$
                                       count=[n_lon,n_lat,n_years,n_days_per_year]))*multiplier
   ENDIF ELSE $
      variable=REFORM(OPEN_AND_EXTRACT(infile,varname))*multiplier
   IF TOTAL(where(variable ge missing)) gt -1 THEN $
      variable[where(variable ge missing)]=!Values.F_NaN
   variable_stddev=fltarr(n_lon,n_lat)
   FOR i=0,n_lon-1 DO BEGIN
      FOR j=0,n_lat-1 DO BEGIN
         IF (KEYWORD_SET(mask_plot)) THEN BEGIN
            IF (lsm(i,j) ne 1) THEN BEGIN
               IF N_ELEMENTS(where(FINITE(variable(i,j,*,*)) eq 1)) ge 2 THEN BEGIN
                  thispt_ts=REFORM(variable(i,j,*,*))
                  variable_stddev(i,j)=STDDEV(thispt_ts,/NaN)
               ENDIF ELSE $
                  variable_stddev(i,j)=!Values.F_NaN         
            ENDIF ELSE $
               variable_stddev(i,j)=!Values.F_NaN
         ENDIF ELSE BEGIN
            IF N_ELEMENTS(where(FINITE(variable(i,j,*,*)) eq 1)) ge 2 THEN BEGIN
               thispt_ts=REFORM(variable(i,j,*,*))
               variable_stddev(i,j)=STDDEV(thispt_ts,/NaN)
            ENDIF ELSE $
               variable_stddev(i,j)=!Values.F_NaN         
         ENDELSE
      ENDFOR
   ENDFOR         
ENDIF

psfile='/home/ss901165/idl/mjo_indices/mjo_indices_filtered_variance.'+runid+'_'+varname+'_filter'+$
       STRTRIM(STRING(upper_limit),1)+STRTRIM(STRING(lower_limit),1)+'.ps'
IF KEYWORD_SET(box) THEN BEGIN
   PSOPEN,file=psfile,FONT=2,CHARSIZE=130,MARGIN=1200,SPACE2=1500,XOFFSET=500,YOFFSET=2000,TFONT=2,TCHARSIZE=100,SPACE3=600,XSIZE=26000,YSIZE=(ABS(box(2))+ABS(box(0)))/180.*21000
ENDIF ELSE $
   PSOPEN,file=psfile,FONT=2,CHARSIZE=130,MARGIN=1200,SPACE2=500,XOFFSET=1000,YOFFSET=1000,TFONT=2,TCHARSIZE=100,SPACE3=600
CS,SCALE=26,NCOLS=N_ELEMENTS(mylevs)+1
LEVS,MANUAL=mylevs
IF KEYWORD_SET(box) THEN BEGIN
   MAP,LONMIN=box(1),LONMAX=box(3),LATMIN=box(0),LATMAX=box(2)
ENDIF ELSE $
   MAP
pink=FSC_COLOR("white",50)
CON,X=longitude,Y=latitude,FIELD=variable_stddev,CB_TITLE='Standard deviation in filtered '+varname+' ('+units+')',$
    TITLE='Standard deviation in '+STRTRIM(STRING(upper_limit),1)+'-'+STRTRIM(STRING(lower_limit),1)+' day filtered '+varname+$
    ' for '+STRTRIM(STRING(n_years),1)+' years of '+runid,/NOLINELABELS,THICK=80,COL=50
PSCLOSE,/NOVIEW

END
