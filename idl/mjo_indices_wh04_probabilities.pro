PRO mjo_indices_wh04_probabilities,input_file,runid,start_year=start_year,nyear=nyear,$
                                   obs_input_file=obs_input_file,obs_start_year=obs_start_year,obs_nyear=obs_nyear,$
                                   prev_day=prev_day,amp_threshold=amp_threshold
  
; Procedure to plot probability of occurrence of each phase, using the
; WH04 diagram and EOFs.

ncid_in=NCDF_OPEN(input_file)
varid=NCDF_VARID(ncid_in,'rmm1')
varstruct=NCDF_VARINQ(ncid_in,varid)

NCDF_DIMINQ,ncid_in,varstruct.dim(1),time_name,ndays_per_year
IF KEYWORD_SET(nyear) THEN BEGIN
   our_nyear=nyear
ENDIF ELSE $
   NCDF_DIMINQ,ncid_in,varstruct.dim(0),year_name,our_nyear
IF KEYWORD_SET(start_year) THEN BEGIN
   our_start_year=start_year
ENDIF ELSE $
   our_start_year=0

; Keyword PREV_DAY switches this program to compute probabilities
; of persistence/next/previous from the *day before*
; strong MJO activity in each phase, not the day after.
; Probability of "death" becomes probability of "birth".
; Probability of "illness" becomes probability of "recovery".
; Probability of "genesis" becomes probability of "death".
IF KEYWORD_SET(prev_day) THEN BEGIN
   lag_offset=-1
ENDIF ELSE $
   lag_offset=1
   
; Threshold amplitude for "strong" MJO activity.
IF KEYWORD_SET(amp_threshold) THEN BEGIN
   our_amp_threshold=amp_threshold 
ENDIF ELSE $
   our_amp_threshold=[1,1,1,1,1,1,1,1,1]

; Read phase and amplitude values for each year and bin
freq_phase=fltarr(9)
freq_phase(*)=0.
freq_phase_persist=fltarr(9)
freq_phase_persist(*)=0.
freq_phase_prev=fltarr(9)
freq_phase_prev(*)=0.
freq_phase_next=fltarr(9)
freq_phase_next(*)=0.
freq_phase_death=fltarr(9)
freq_phase_death(*)=0.
freq_phase_ill=fltarr(9)
freq_phase_ill(*)=0.
freq_phase_gen=fltarr(9)
freq_phase_gen(*)=0.
ndays=ndays_per_year*our_nyear

; Read missing value
NCDF_ATTGET,ncid_in,varid,'missing_value',missing

phase=REFORM(OPEN_AND_EXTRACT(input_file,'phase',offset=[our_start_year,0],count=[our_nyear,ndays_per_year]))
amplitude=REFORM(OPEN_AND_EXTRACT(input_file,'amplitude',$
                                  offset=[our_start_year,0],count=[our_nyear,ndays_per_year]))
rmm1=REFORM(OPEN_AND_EXTRACT(input_file,'rmm1',$
                             offset=[our_start_year,0],count=[our_nyear,ndays_per_year]))
rmm2=REFORM(OPEN_AND_EXTRACT(input_file,'rmm2',$
                             offset=[our_start_year,0],count=[our_nyear,ndays_per_year]))

phase_ts=intarr(ndays)
amplitude_ts=fltarr(ndays)
rmm1_ts=fltarr(ndays)
rmm2_ts=fltarr(ndays)
FOR i=0,our_nyear-1 DO BEGIN
   start=(i*ndays_per_year)
   stop=(i+1)*ndays_per_year-1
   phase_ts(start:stop)=phase(i,*)
   amplitude_ts(start:stop)=amplitude(i,*)
   rmm1_ts(start:stop)=rmm1(i,*)
   rmm2_ts(start:stop)=rmm2(i,*)
ENDFOR

IF TOTAL(where(rmm1_ts eq missing)) ge 1 THEN $
   rmm1_ts[where(rmm1_ts eq missing)]=!Values.F_NaN
IF TOTAL(where(rmm2_ts eq missing)) ge 1 THEN $
   rmm2_ts[where(rmm2_ts eq missing)]=!Values.F_NaN

n_lag_days=30
lag_composites=fltarr(9,n_lag_days,2)
lag_composites_genesis=fltarr(9,n_lag_days,2)
lag_composites(*,*,*)=0.
lag_composites_genesis(*,*,*)=0.

FOR j=LONG(1),LONG(ndays-2) DO BEGIN
   IF amplitude_ts(j) eq missing THEN BEGIN
      ndays=ndays-1
   ENDIF ELSE BEGIN
      IF amplitude_ts(j) lt our_amp_threshold(phase_ts(j)) THEN BEGIN
         freq_phase(0)=freq_phase(0)+1
         IF amplitude_ts(j+lag_offset) lt our_amp_threshold(phase_ts(j)) THEN $
            freq_phase_persist(0)=freq_phase_persist(0)+1
      ENDIF ELSE BEGIN
         freq_phase(phase_ts(j))=freq_phase(phase_ts(j))+1
         IF j lt ndays-n_lag_days-1 THEN BEGIN
                                ; print,'Adding to lag
                                ; composite',lag_composites(phase_ts(j),0,0)
            FOR k=0,n_lag_days-1 DO BEGIN
               IF FINITE(rmm1_ts(j+k)) eq 1 THEN $
                  lag_composites(phase_ts(j),k,0)=lag_composites(phase_ts(j),k,0)+rmm1_ts(j+k)
               IF FINITE(rmm2_ts(j+k)) eq 1 THEN $
                  lag_composites(phase_ts(j),k,1)=lag_composites(phase_ts(j),k,1)+rmm2_ts(j+k)
            ENDFOR
         ENDIF
         IF amplitude_ts(j+lag_offset) gt our_amp_threshold(phase_ts(j)) and phase_ts(j+lag_offset) eq phase_ts(j) THEN $
            freq_phase_persist(phase_ts(j))=freq_phase_persist(phase_ts(j))+1            
         IF phase_ts(j) ne 8 THEN BEGIN
            IF amplitude_ts(j+lag_offset) gt our_amp_threshold(phase_ts(j)) and phase_ts(j+lag_offset)-phase_ts(j) eq 1 THEN $
               freq_phase_next(phase_ts(j))=freq_phase_next(phase_ts(j))+1
         ENDIF ELSE $
            IF amplitude_ts(j+lag_offset) gt our_amp_threshold(phase_ts(j)) and phase_ts(j+lag_offset) eq 1 THEN $
               freq_phase_next(8)=freq_phase_next(8)+1
         IF phase_ts(j) ne 1 THEN BEGIN
            IF amplitude_ts(j+lag_offset) gt our_amp_threshold(phase_ts(j)) and phase_ts(j+lag_offset)-phase_ts(j) eq -1 THEN $
               freq_phase_prev(phase_ts(j))=freq_phase_prev(phase_ts(j))+1
         ENDIF ELSE $
            IF amplitude_ts(j+lag_offset) gt our_amp_threshold(phase_ts(j)) and phase_ts(j+lag_offset) eq 8 THEN $
               freq_phase_prev(1)=freq_phase_prev(1)+1         
         IF amplitude_ts(j+lag_offset) lt our_amp_threshold(phase_ts(j)) and j le ndays-12 THEN BEGIN
            count=0
            FOR k=j+lag_offset,j+10*(SIGN(lag_offset)),SIGN(lag_offset) DO BEGIN
               IF amplitude_ts(k) lt our_amp_threshold(phase_ts(j)) THEN $
                  count=count+1
            ENDFOR
            IF count eq 10 THEN BEGIN
               freq_phase_death(phase_ts(j))=freq_phase_death(phase_ts(j))+1
               this_amplitude=0.
               k=j+lag_offset
               WHILE this_amplitude lt our_amp_threshold(phase_ts(j)) DO BEGIN
                  k=k+lag_offset
                  this_amplitude=amplitude_ts(k)                  
                  code=0
                  IF k eq ndays-2 or k eq 1 THEN BEGIN
                     code=-99
                     this_amplitude=99
                  ENDIF
               ENDWHILE
               IF code eq 0 THEN BEGIN
                  freq_phase_gen(phase_ts(k))=freq_phase_gen(phase_ts(k))+1
                  IF k lt ndays-n_lag_days-1 THEN BEGIN
                     lag_composites_genesis(phase_ts(k),*,0)=lag_composites_genesis(phase_ts(k),*,0)+rmm1_ts(k:k+n_lag_days-1)
                     lag_composites_genesis(phase_ts(k),*,1)=lag_composites_genesis(phase_ts(k),*,1)+rmm2_ts(k:k+n_lag_days-1)
                  ENDIF
               ENDIF
            ENDIF ELSE IF count lt 10 THEN $
               freq_phase_ill(phase_ts(j))=freq_phase_ill(phase_ts(j))+1
         ENDIF
      ENDELSE
   ENDELSE
ENDFOR
freq_phase_persist=freq_phase_persist/FLOAT(freq_phase)
freq_phase_death=freq_phase_death/FLOAT(freq_phase)
freq_phase_ill=freq_phase_ill/FLOAT(freq_phase)
freq_phase_prev=freq_phase_prev/FLOAT(freq_phase)
freq_phase_next=freq_phase_next/FLOAT(freq_phase)
FOR i=1,8 DO BEGIN
   lag_composites(i,*,0)=lag_composites(i,*,0)/FLOAT(freq_phase(i))
   lag_composites(i,*,1)=lag_composites(i,*,1)/FLOAT(freq_phase(i))
   lag_composites_genesis(i,*,0)=lag_composites_genesis(i,*,0)/FLOAT(freq_phase_gen(i))
   lag_composites_genesis(i,*,1)=lag_composites_genesis(i,*,1)/FLOAT(freq_phase_gen(i))
ENDFOR
freq_phase_gen=freq_phase_gen/TOTAL(freq_phase_gen)
freq_phase=freq_phase/FLOAT(ndays)

;print,freq_phase_persist(3)+freq_phase_death(3)+freq_phase_ill(3)+freq_phase_prev(3)+freq_phase_next(3)

points=(2*!PI/99.0)*findgen(100)
x=COS(points)
y=SIN(points)
mylevs=[0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09]
;mylevs=[0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,0.10,0.105]
u=fltarr(1,1)
v=fltarr(1,1)

IF KEYWORD_SET(prev_day) THEN BEGIN
   psfile='/home/ss901165/idl/mjo_indices/mjo_indices_wh04_probabilities.freq_phase.'+runid+'.prev_day.ps'
ENDIF ELSE $
   psfile='/home/ss901165/idl/mjo_indices/mjo_indices_wh04_probabilities.freq_phase.'+runid+'.ps'
PSOPEN,file=psfile,FONT=2,CHARSIZE=120,MARGIN=1500,SPACE3=1500,XOFFSET=500,YOFFSET=500,TFONT=2,TCHARSIZE=100,CB_WIDTH=110,$
       XSIZE=17000,YSIZE=17000
GSET,XMIN=-4,XMAX=4,YMIN=-4,YMAX=4,TITLE='Fractional frequency of occurrence for Wheeler and Hendon phases - '+runid
GPLOT,X=REPLICATE(0,4),Y=[-4,-3,-2,-1],STYLE=0,THICK=80
GPLOT,X=REPLICATE(0,4),Y=[1,2,3,4],STYLE=0,THICK=80
GPLOT,X=[1,2,3,4],Y=REPLICATE(0,4),STYLE=0,THICK=80
GPLOT,X=[-1,-2,-3,-4],Y=REPLICATE(0,4),STYLE=0,THICK=80
GPLOT,X=[SQRT(2)/2.,1,2,3,4],Y=[SQRT(2)/2.,1,2,3,4],STYLE=0,THICK=80
GPLOT,X=[-SQRT(2)/2.,-1,-2,-3,-4],Y=[SQRT(2)/2.,1,2,3,4],STYLE=0,THICK=80
GPLOT,X=[-SQRT(2)/2.,-1,-2,-3,-4],Y=[-SQRT(2)/2.,-1,-2,-3,-4],STYLE=0,THICK=80
GPLOT,X=[SQRT(2)/2.,1,2,3,4],Y=[-SQRT(2)/2.,-1,-2,-3,-4],STYLE=0,THICK=80

CS,SCALE=26,NCOLS=N_ELEMENTS(mylevs)+2
GPLOT,X=[1,4,4,3,2,1,SQRT(2)/2.],Y=[0,0,-4,-3,-2,-1,-SQRT(2)/2.],FILLCOL=NEAREST_LOW(mylevs,freq_phase(4))+2,THICK=0
GPLOT,X=[0,0,0,0,4,3,2,1,SQRT(2)/2.],Y=[-1,-2,-3,-4,-4,-3,-2,-1,-SQRT(2)/2.],FILLCOL=NEAREST_LOW(mylevs,freq_phase(3))+2,THICK=0
GPLOT,X=[0,0,0,0,-4,-3,-2,-1,-SQRT(2)/2.],Y=[-1,-2,-3,-4,-4,-3,-2,-1,-SQRT(2)/2.],FILLCOL=NEAREST_LOW(mylevs,freq_phase(2))+2,THICK=0
GPLOT,X=[-4,-3,-2,-1,-SQRT(2)/2.,-1,-2,-3,-4],Y=[-4,-3,-2,-1,-SQRT(2)/2.,0,0,0,0],FILLCOL=NEAREST_LOW(mylevs,freq_phase(1))+2,THICK=0
GPLOT,X=[-1,-2,-3,-4,-4,-3,-2,-1,-SQRT(2)/2.],Y=[0,0,0,0,4,3,2,1,SQRT(2)/2.],FILLCOL=NEAREST_LOW(mylevs,freq_phase(8))+2,THICK=0
GPLOT,X=[-4,-3,-2,-1,-SQRT(2)/2.,0,0,0,0],Y=[4,3,2,1,SQRT(2)/2.,1,2,3,4],FILLCOL=NEAREST_LOW(mylevs,freq_phase(7))+2,THICK=0
GPLOT,X=[0,0,0,0,4,3,2,1,SQRT(2)/2.],Y=[1,2,3,4,4,3,2,1,SQRT(2)/2.],FILLCOL=NEAREST_LOW(mylevs,freq_phase(6))+2,THICK=0
GPLOT,X=[4,3,2,1,SQRT(2)/2.,1,2,3,4],Y=[4,3,2,1,SQRT(2)/2.,0,0,0,0],FILLCOL=NEAREST_LOW(mylevs,freq_phase(5))+2,THICK=0

white=FSC_COLOR("white",28)
GPLOT,X=x,Y=y,FILLCOL=28

title_xpos=[-1.25,-0.30, 0.30, 1.20,1.50,0.30,-0.30,-1.25]
title_ypos=[-0.60,-1.50,-1.50,-0.60,0.50,1.50, 1.50, 0.50]
title_align=[1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0]
column_xpos=[-2.4,-0.3,1.7,3.8,3.8,1.7,-0.3,-2.4]
column_start_ypos=[-0.7,-2.0,-2.0,-0.7,2.0,3.5,3.5,2.0]
row_spacing=0.3

u_vector=[-3./SQRT(2),-1,1,3./SQRT(2),3./SQRT(2),1,-1,-3./SQRT(2)]*SIGN(lag_offset)
v_vector=[-1,-3./SQRT(2),-3./SQRT(2),-1,1,3./SQRT(2),3./SQRT(2),1]*SIGN(lag_offset)
x_vector=[-0.95,-0.4,0.4,0.95,0.95,0.4,-0.4,-0.95]
y_vector=[-0.4,-0.95,-0.95,-0.4,0.4,0.95,0.95,0.4]

gentext_xpos=[-0.75,-0.30,0.30,0.75,0.67,0.28,-0.28,-0.67]
gentext_ypos=[-0.30,-0.75,-0.75,-0.30,0.28,0.67,0.67,0.28]
gentext_orientation=[-60,-30,30,60,-60,-30,30,60]

GPLOT,X=0,Y=0.1,TEXT='Freq. = '+STRMID(STRTRIM(STRING(freq_phase(0)),1),0,5),CHARSIZE=85
GPLOT,X=0,Y=-0.1,TEXT='Persist = '+STRMID(STRTRIM(STRING(freq_phase_persist(0)),1),0,5),CHARSIZE=85

FOR i=0,7 DO BEGIN
   phase=i+1
   IF freq_phase(phase) gt 0.030 and freq_phase(phase) lt 0.050 THEN BEGIN
      text_color=28
   ENDIF ELSE $
      text_color=1
   GPLOT,X=title_xpos(i),Y=title_ypos(i),TEXT='Phase '+STRTRIM(STRING(phase),1),ALIGN=title_align(i),COL=text_color
   GPLOT,X=column_xpos(i),Y=column_start_ypos(i),TEXT='Freq. = '+STRMID(STRTRIM(STRING(freq_phase(phase)),1),0,5),ALIGN=1.0,COL=text_color
   GPLOT,X=column_xpos(i),Y=column_start_ypos(i)-row_spacing,TEXT='Persist = '+STRMID(STRTRIM(STRING(freq_phase_persist(phase)),1),0,5),ALIGN=1.0,COL=text_color   
   GPLOT,X=column_xpos(i),Y=column_start_ypos(i)-2*row_spacing,TEXT='Prev = '+STRMID(STRTRIM(STRING(freq_phase_prev(phase)),1),0,5),ALIGN=1.0,COL=text_color
   GPLOT,X=column_xpos(i),Y=column_start_ypos(i)-3*row_spacing,TEXT='Next = '+STRMID(STRTRIM(STRING(freq_phase_next(phase)),1),0,5),ALIGN=1.0,COL=text_color
   IF KEYWORD_SET(prev_day) THEN BEGIN
      GPLOT,X=column_xpos(i),Y=column_start_ypos(i)-4*row_spacing,TEXT='Birth = '+STRMID(STRTRIM(STRING(freq_phase_death(phase)),1),0,5),ALIGN=1.0,COL=text_color
      GPLOT,X=column_xpos(i),Y=column_start_ypos(i)-5*row_spacing,TEXT='Well = '+STRMID(STRTRIM(STRING(freq_phase_ill(phase)),1),0,5),ALIGN=1.0,COL=text_color
   ENDIF ELSE BEGIN
      GPLOT,X=column_xpos(i),Y=column_start_ypos(i)-4*row_spacing,TEXT='Death = '+STRMID(STRTRIM(STRING(freq_phase_death(phase)),1),0,5),ALIGN=1.0,COL=text_color
      GPLOT,X=column_xpos(i),Y=column_start_ypos(i)-5*row_spacing,TEXT='Ill = '+STRMID(STRTRIM(STRING(freq_phase_ill(phase)),1),0,5),ALIGN=1.0,COL=text_color
   ENDELSE
   u(0,0)=u_vector(i)
   v(0,0)=v_vector(i)
   IF KEYWORD_SET(prev_day) THEN BEGIN
      VECT,U=REVERSE(u),V=REVERSE(v),X=x_vector(i),Y=y_vector(i),TYPE=3,/NOLEGEND,LENGTH=FLOOR(freq_phase_gen(phase)*1000.)+100
   ENDIF ELSE $
      VECT,U=u,V=v,X=x_vector(i),Y=y_vector(i),TYPE=3,/NOLEGEND,LENGTH=FLOOR(freq_phase_gen(phase)*1000.)+100
   GPLOT,X=gentext_xpos(i),Y=gentext_ypos(i),TEXT=STRMID(STRTRIM(STRING(freq_phase_gen(phase)),1),0,4),ORIENTATION=gentext_orientation(i),CHARSIZE=80
ENDFOR

GPLOT,X=0,Y=-4.5,TEXT='Indian Ocean',ALIGN=0.5,CHARSIZE=100
GPLOT,X=4.5,Y=0,TEXT='Maritime Continent',ALIGN=0.5,CHARSIZE=100,ORIENTATION=90
GPLOT,X=0,Y=4.25,TEXT='Western Pacific',ALIGN=0.5,CHARSIZE=100
GPLOT,X=-4.5,Y=0,TEXT='Western Hemisphere and Africa',ALIGN=0.5,CHARSIZE=100,ORIENTATION=90
;COLBAR,COORDS=[22000,3000,22500,18000],LEVS=STRMID(STRTRIM(STRING(mylevs),1),0,5),TITLE='Fractional frequency',/LOWER,/UPPER
AXES,XSTEP=1,YSTEP=1,XMINOR=0.5,YMINOR=0.5,XTITLE='RMM1',YTITLE='RMM2'

PSCLOSE,/NOVIEW

psfile='/home/ss901165/idl/mjo_indices/mjo_indices_wh04_probabilities.lag_composite.'+runid+'.ps'
PSOPEN,file=psfile,FONT=2,CHARSIZE=120,MARGIN=2000,SPACE3=1500,XOFFSET=1000,YOFFSET=500,TFONT=2,TCHARSIZE=100,CB_WIDTH=110,$
       XSIZE=16000,YSIZE=16000
GSET,XMIN=-2,XMAX=2,YMIN=-2,YMAX=2,TITLE='Lag composite evolution of strong MJO activity in each phase - '+runid
GPLOT,X=REPLICATE(0,2),Y=[-2,-1],STYLE=0,THICK=80
GPLOT,X=REPLICATE(0,2),Y=[1,2],STYLE=0,THICK=80
GPLOT,X=[1,2],Y=REPLICATE(0,2),STYLE=0,THICK=80
GPLOT,X=[-1,-2],Y=REPLICATE(0,2),STYLE=0,THICK=80
GPLOT,X=[SQRT(2)/2.,1,2],Y=[SQRT(2)/2.,1,2],STYLE=0,THICK=80
GPLOT,X=[-SQRT(2)/2.,-1,-2],Y=[SQRT(2)/2.,1,2],STYLE=0,THICK=80
GPLOT,X=[-SQRT(2)/2.,-1,-2],Y=[-SQRT(2)/2.,-1,-2],STYLE=0,THICK=80
GPLOT,X=[SQRT(2)/2.,1,2],Y=[-SQRT(2)/2.,-1,-2],STYLE=0,THICK=80
x=COS(points)
y=SIN(points)
white=FSC_COLOR("white",28)
GPLOT,X=x,Y=y,FILLCOL=28

CS,SCALE=26,NCOLS=9
;FOR i=0,n_lag_days-1 DO BEGIN
;lag_composites(*,i,0)=lag_composites(*,i,0)*1.02^(i*1.3)*1.319
;lag_composites(*,i,1)=lag_composites(*,i,1)*1.02^(i*1.3)*1.328
;ENDFOR
FOR i=1,8 DO BEGIN
   GPLOT,X=lag_composites(i,0,0),Y=lag_composites(i,0,1),THICK=200,SYM=1,/NOLINES 
   GPLOT,X=lag_composites(i,*,0),Y=lag_composites(i,*,1),THICK=200,STYLE=0,COL=2+i
   ;GPLOT,X=lag_composites_genesis(i,*,0),Y=lag_composites_genesis(i,*,1),THICK=200,STYLE=1,COL=2+i
   FOR j=0,n_lag_days-1,5 DO $
      GPLOT,X=lag_composites(i,j,0),Y=lag_composites(i,j,1),THICK=150,SYM=2,/NOLINES,SIZE=80   
ENDFOR

GPLOT,X=0,Y=-2.35,TEXT='Indian Ocean',ALIGN=0.5,CHARSIZE=100
GPLOT,X=2.25,Y=0,TEXT='Maritime Continent',ALIGN=0.5,CHARSIZE=100,ORIENTATION=90
GPLOT,X=0,Y=2.15,TEXT='Western Pacific',ALIGN=0.5,CHARSIZE=100
GPLOT,X=-2.45,Y=0,TEXT='Western Hemisphere and Africa',ALIGN=0.5,CHARSIZE=100,ORIENTATION=90
AXES,XSTEP=0.5,YSTEP=0.5,XMINOR=0.25,YMINOR=0.25,NDECS=2
GPLOT,X=0,Y=-2.5,TEXT='RMM1'
GPLOT,X=-2.6,Y=0,TEXT='RMM2',ORIENTATION=90
PSCLOSE,/NOVIEW
;STOP

IF KEYWORD_SET(obs_input_file) THEN BEGIN

   obs_ncid_in=NCDF_OPEN(obs_input_file)
   obs_varid=NCDF_VARID(obs_ncid_in,'rmm1')
   obs_varstruct=NCDF_VARINQ(obs_ncid_in,obs_varid)
   
   ; Read missing value
   NCDF_ATTGET,obs_ncid_in,obs_varid,'missing_value',missing

   NCDF_DIMINQ,obs_ncid_in,obs_varstruct.dim(1),obs_time_name,obs_ndays_per_year
   IF KEYWORD_SET(obs_nyear) THEN BEGIN
      our_obs_nyear=obs_nyear
   ENDIF ELSE $
      NCDF_DIMINQ,obs_ncid_in,obs_varstruct.dim(0),obs_year_name,our_obs_nyear
   IF KEYWORD_SET(obs_start_year) THEN BEGIN
      our_obs_start_year=obs_start_year
   ENDIF ELSE $
      our_obs_start_year=0

   freq_obs_phase=fltarr(9)
   freq_obs_phase(*)=0.
   freq_obs_phase_persist=fltarr(9)
   freq_obs_phase_persist(*)=0.
   freq_obs_phase_next=fltarr(9)
   freq_obs_phase_next(*)=0.
   freq_obs_phase_prev=fltarr(9)
   freq_obs_phase_prev(*)=0.
   freq_obs_phase_death=fltarr(9)
   freq_obs_phase_death(*)=0.
   freq_obs_phase_ill=fltarr(9)
   freq_obs_phase_ill(*)=0.
   freq_obs_phase_gen=fltarr(9)
   freq_obs_phase_gen(*)=0.
   obs_ndays=obs_nyear*obs_ndays_per_year
   phase=REFORM(OPEN_AND_EXTRACT(obs_input_file,'phase',$
                                       offset=[our_obs_start_year,0],count=[our_obs_nyear,obs_ndays_per_year]))
   amplitude=REFORM(OPEN_AND_EXTRACT(obs_input_file,'amplitude',$
                                     offset=[our_obs_start_year,0],count=[our_obs_nyear,obs_ndays_per_year]))
   obs_phase_ts=intarr(obs_ndays)
   obs_amplitude_ts=fltarr(obs_ndays)
   FOR i=0,our_obs_nyear-1 DO BEGIN
      start=(i*obs_ndays_per_year)
      stop=(i+1)*obs_ndays_per_year-1
      obs_phase_ts(start:stop)=phase(i,*)
      obs_amplitude_ts(start:stop)=amplitude(i,*)
   ENDFOR
   
   FOR j=1,obs_ndays-2 DO BEGIN
      IF obs_amplitude_ts(j) eq missing THEN BEGIN
         obs_ndays=obs_ndays-1
      ENDIF ELSE BEGIN
         IF obs_amplitude_ts(j) lt 1 THEN BEGIN
            freq_obs_phase(0)=freq_obs_phase(0)+1
            IF obs_amplitude_ts(j+lag_offset) lt 1 THEN $
               freq_obs_phase_persist(0)=freq_obs_phase_persist(0)+1
         ENDIF ELSE BEGIN
            freq_obs_phase(obs_phase_ts(j))=freq_obs_phase(obs_phase_ts(j))+1
            IF obs_amplitude_ts(j+lag_offset) gt 1 and obs_phase_ts(j+lag_offset) eq obs_phase_ts(j) THEN $
               freq_obs_phase_persist(obs_phase_ts(j))=freq_obs_phase_persist(obs_phase_ts(j))+1            
            IF obs_phase_ts(j) ne 8 THEN BEGIN
               IF obs_amplitude_ts(j+lag_offset) gt 1 and obs_phase_ts(j+lag_offset)-obs_phase_ts(j) eq 1 THEN $
                  freq_obs_phase_next(obs_phase_ts(j))=freq_obs_phase_next(obs_phase_ts(j))+1
            ENDIF ELSE $
               IF obs_amplitude_ts(j+lag_offset) gt 1 and obs_phase_ts(j+lag_offset) eq 1 THEN $
                  freq_obs_phase_next(8)=freq_obs_phase_next(8)+1
            IF obs_phase_ts(j) ne 1 THEN BEGIN
               IF obs_amplitude_ts(j+lag_offset) gt 1 and obs_phase_ts(j+lag_offset)-obs_phase_ts(j) eq -1 THEN $
                  freq_obs_phase_prev(obs_phase_ts(j))=freq_obs_phase_prev(obs_phase_ts(j))+1
            ENDIF ELSE $
               IF obs_amplitude_ts(j+lag_offset) gt 1 and obs_phase_ts(j+lag_offset) eq 8 THEN $
                  freq_obs_phase_prev(1)=freq_obs_phase_prev(1)+1        
            IF obs_amplitude_ts(j+lag_offset) lt 1 and j le obs_ndays-12 THEN BEGIN
               count=0
               FOR k=j+lag_offset,j+10*SIGN(lag_offset),SIGN(lag_offset) DO BEGIN
                  IF obs_amplitude_ts(k) lt 1 THEN $
                     count=count+1
               ENDFOR
               IF count eq 10 THEN BEGIN
                  freq_obs_phase_death(obs_phase_ts(j))=freq_obs_phase_death(obs_phase_ts(j))+1
                  this_amplitude=0.
                  k=j+lag_offset
                  WHILE this_amplitude le 1 DO BEGIN
                     k=k+lag_offset
                     this_amplitude=obs_amplitude_ts(k)                  
                     code=0
                     IF k eq obs_ndays-2 or k eq 1 THEN BEGIN
                        code=-99
                        this_amplitude=99
                     ENDIF
                  ENDWHILE
                  IF code eq 0 THEN $
                     freq_obs_phase_gen(obs_phase_ts(k))=freq_obs_phase_gen(obs_phase_ts(k))+1
               ENDIF ELSE IF count lt 10 THEN $
                  freq_obs_phase_ill(obs_phase_ts(j))=freq_obs_phase_ill(obs_phase_ts(j))+1               
            ENDIF
         ENDELSE
      ENDELSE
   ENDFOR
   freq_obs_phase_persist=freq_obs_phase_persist/FLOAT(freq_obs_phase)
   freq_obs_phase_ill=freq_obs_phase_ill/FLOAT(freq_obs_phase)
   freq_obs_phase_death=freq_obs_phase_death/FLOAT(freq_obs_phase)
   freq_obs_phase_prev=freq_obs_phase_prev/FLOAT(freq_obs_phase)
   freq_obs_phase_next=freq_obs_phase_next/FLOAT(freq_obs_phase)
   freq_obs_phase=freq_obs_phase/FLOAT(obs_ndays-1)
   freq_obs_phase_gen=freq_obs_phase_gen/FLOAT(TOTAL(freq_obs_phase_gen))

   IF KEYWORD_SET(prev_day) THEN BEGIN
      psfile='/home/ss901165/idl/mjo_indices/mjo_indices_wh04_probabilities.freq_phase.'+runid+'_lines_compare_obs.prev_day.ps'
   ENDIF ELSE $
      psfile='/home/ss901165/idl/mjo_indices/mjo_indices_wh04_probabilities.freq_phase.'+runid+'_lines_compare_obs.ps'
   PSOPEN,file=psfile,FONT=2,CHARSIZE=110,MARGIN=1500,SPACE2=100,XOFFSET=1200,YOFFSET=6000,TFONT=2,TCHARSIZE=100,CB_WIDTH=110,$
          XSIZE=23000,YSIZE=12000,SPACE3=100   
   GSET,XMIN=0.5,XMAX=8.5,YMIN=0,YMAX=0.18,TITLE='Comparison of MJO phase frequency statistics for model run '+runid+' against observations'
   phase_numbers=indgen(8)+1
   black=FSC_COLOR("black",2)
   red=FSC_COLOR("red",3)
   blue=FSC_COLOR("blue",4)
   orange=FSC_COLOR("orange",5)
   brown=FSC_COLOR("brown",6)
   cyan=FSC_COLOR("cyan",7)
   purple=FSC_COLOR("purple",8)
   GPLOT,X=phase_numbers,Y=freq_obs_phase(1:8),STYLE=0,SYM=6,COL=2
   GPLOT,X=phase_numbers,Y=freq_phase(1:8),STYLE=2,SYM=4,COL=2
   GPLOT,X=phase_numbers,Y=freq_obs_phase_prev(1:8),STYLE=0,SYM=6,COL=6
   GPLOT,X=phase_numbers,Y=freq_phase_prev(1:8),STYLE=2,SYM=4,COL=6
   GPLOT,X=phase_numbers,Y=freq_obs_phase_next(1:8),STYLE=0,SYM=6,COL=3
   GPLOT,X=phase_numbers,Y=freq_phase_next(1:8),STYLE=2,SYM=4,COL=3
   GPLOT,X=phase_numbers,Y=freq_obs_phase_ill(1:8),STYLE=0,SYM=6,COL=7
   GPLOT,X=phase_numbers,Y=freq_phase_ill(1:8),STYLE=2,SYM=4,COL=7   
   GPLOT,X=phase_numbers,Y=freq_obs_phase_death(1:8),STYLE=0,SYM=6,COL=4
   GPLOT,X=phase_numbers,Y=freq_phase_death(1:8),STYLE=2,SYM=4,COL=4
   AXES,XVALS=[1,2,3,4,5,6,7,8],XLABELS=['1','2','3','4','5','6','7','8'],YSTEP=0.02,$
        YTITLE='Probability',XTITLE='MJO phase (Wheeler and Hendon 2004)',NDECS=2,/NORIGHT
   
   GSET,XMIN=0.5,XMAX=8.5,YMIN=0.3,YMAX=0.901
   GPLOT,X=phase_numbers,Y=freq_obs_phase_persist(1:8),STYLE=0,SYM=6,COL=5
   GPLOT,X=phase_numbers,Y=freq_phase_persist(1:8),STYLE=2,SYM=4,COL=5
   GPLOT,X=phase_numbers,Y=freq_obs_phase_gen(1:8)+0.5,STYLE=0,SYM=6,COL=8
   GPLOT,X=phase_numbers,Y=freq_phase_gen(1:8)+0.5,STYLE=2,SYM=4,COL=8
   AXES,XVALS=[1,2,3,4,5,6,7,8],XLABELS=['1','2','3','4','5','6','7','8'],YSTEP=0.05,YTITLE='Probability',NDECS=2,/ONLYRIGHT
   
   IF KEYWORD_SET(prev_day) THEN BEGIN
      items=['Normalized prob. (+0.5) of MJO death (right axis; normalized by number of death events)',$
             'Probability of MJO birth (amplitude < 1 for previous >= 10 days)',$
             'Probability of MJO recovery (amplitude < 1 for previous < 10 days)',$
             'Probability of transition from next phase (with amplitude > 1) on previous day',$
             'Probabilitiy of transition from previous phase (with amplitude > 1) on previous day',$
             'Probability of same phase with amplitude > 1 on previous day (lag 1 persist, right axis)',$
             'Probability of occurence of phase with amplitude > 1',$
             'Model '+runid+' ('+STRTRIM(STRING(our_nyear),1)+' years)','Observations ('+STRTRIM(STRING(obs_nyear),1)+' years)']
   ENDIF ELSE BEGIN
      items=['Normalized prob. (+0.5) of MJO genesis, after death (right axis; normalized by number of genesis events)',$
             'Probability of MJO death (amplitude < 1 for next >= 10 days)',$
             'Probability of MJO illness (amplitude < 1 for next < 10 days)',$
             'Probability of transition to next phase (with amplitude > 1) on next day',$
             'Probabilitiy of transition to previous phase (with amplitude > 1) on next day',$
             'Probability of same phase with amplitude > 1 on next day (lag 1 persist, right axis)',$
             'Probability of occurence of phase with amplitude > 1',$
             'Model '+runid+' ('+STRTRIM(STRING(our_nyear),1)+' years)','Observations ('+STRTRIM(STRING(obs_nyear),1)+' years)']
   ENDELSE
   colors=[8,4,7,3,6,5,2,2,2]
   style=[0,0,0,0,0,0,0,2,0]
   GLEGEND,labels=items,COL=colors,STYLE=style,SYM=sym,LEGXOFFSET=1000,LEGYOFFSET=-1000
   PSCLOSE,/NOVIEW
ENDIF

;STOP

END
