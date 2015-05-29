PRO mjo_indices_spacetime_power_spectra,infile,runid,varname,box

; Procedure to plot space-time power spectra (pre-calculated) from a
; netCDF file.  User needs to provide a "box" of the form
;
; box=[min_wavenumber,min_frequency,max_wavenumber,max_frequency]
;
; for plotting.

frequency=OPEN_AND_EXTRACT(infile,'frequency')
wavenumber=OPEN_AND_EXTRACT(infile,'wavenumber')
powerspec=OPEN_AND_EXTRACT(infile,'powerspec')

n_frequency=N_ELEMENTS(frequency)
n_wavenumber=N_ELEMENTS(wavenumber)

CASE varname OF
   'olr' : BEGIN
      mylevs=['0.004','0.008','0.016','0.032','0.064','0.128','0.256','0.512','1.024']
;      mylevs=['0.0008','0.0016','0.0032','0.0064','0.0128','0.0256','0.0512','0.1024','0.2048']
   END
   'u200' : BEGIN
      mylevs=['0.001','0.002','0.004','0.008','0.016','0.032','0.064','0.128','0.256']
   END
   'u850' : BEGIN
      mylevs=['0.0004','0.0008','0.0016','0.0032','0.0064','0.0128','0.0256','0.0512','0.1024']
   END
   'precip' : BEGIN
      obs_file='/home/ss901165/datasets/MJO_INDICES/noaa_eraint_trmm_obs.jan-dec_dmeans_anom-3harm_space-time_powerspec.1999-2008.precip.2.5x2.5.nc'
      mylevs=['0.0004','0.0008','0.0016','0.0032','0.0064','0.0128','0.0256','0.0512','0.1024']
   END
ENDCASE

min_wave=NEAREST(wavenumber,box(0))
max_wave=NEAREST(wavenumber,box(2))
min_freq=NEAREST(frequency,box(1))
max_freq=NEAREST(frequency,box(3))

period=1./frequency
min_period=period(min_freq)
max_period=period(max_freq)

yvals=[2,3,4,5,6,8,10,15,20,25,30,40,50,60,80,100,120,150,200,250,300,350]
yvals_plot=[0]
FOR i=0,N_ELEMENTS(yvals)-1 DO BEGIN
   IF yvals(i) le min_period and yvals(i) ge max_period THEN $
      yvals_plot=[yvals_plot,yvals(i)]
ENDFOR

psfile='/home/ss901165/idl/mjo_indices/plots/mjo_indices_spacetime_power_spectra.'+runid+'_'+varname+'.ps'
PSOPEN,file=psfile,FONT=2,CHARSIZE=130,MARGIN=2000,SPACE2=800,XOFFSET=1200,YOFFSET=500,TFONT=2,TCHARSIZE=100,SPACE3=50
GSET,XMIN=box(0),YMIN=min_period,XMAX=box(2),YMAX=max_period,TITLE='Space-time power spectrum for latitude-averaged '+varname+' in '+runid,/YLOG
LEVS,MANUAL=mylevs
CS,SCALE=26,NCOLS=N_ELEMENTS(mylevs)+1,white=[2]
CON,X=wavenumber(min_wave:max_wave),Y=period(min_freq:max_freq),FIELD=powerspec(min_wave:max_wave,min_freq:max_freq),/NOLINELABELS,$
    CB_TITLE='Spectral power in '+varname,/NOLINES
GPLOT,X=REPLICATE(0,max_freq-min_freq+1),Y=period(min_freq:max_freq),THICK=250
AXES,XSTEP=1,yvals=yvals_plot(1:N_ELEMENTS(yvals_plot)-1),NDECS=2,$
     XTITLE='<----    Westward          Wavenumber          Eastward    ---->',$ 
     YTITLE='Period (days)'
PSCLOSE,/NOVIEW

END
