  PROGRAM wk99_power_symmetric
    USE nicks_subroutines
    USE netcdf
    IMPLICIT NONE

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! NOTE :
! Daehyun Kim modified Dr. Matthew C. Wheeler's code for the MJO
! Simulation Metrics calculation system
! Note that any bugs or errors might be caused by modifications.
! Contact : Daehyun Kim (kim@climate.snu.ac.kr)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Further modifications by Nick Klingaman to read/write netCDF.
!
!--------------------------------------------------------------------

      integer i,j,t,n,pt,pn,NT,NL,NM,itr,m,nlat,hnlat,nskip,nvar
      integer iy, nyr,ntime
      parameter (itr=1,NT=num_s,NL=num_x,NM=sel_y,nlat=sel_y)
      parameter (nskip=0,nvar=1)
      
      CHARACTER*(*) directory,input,output,varname
      CHARACTER(4) time_name
      parameter(directory=directory)
      parameter(input=directory//'/'//input_file)
      parameter(output=directory//'/'//output_file)      
      parameter(varname=variable_name)
      INTEGER ncid,ncid_out,in_varids(3),out_varids(3),out_dimids(2),status,&
           time_dimid

! Make NT and NL even numbers (it's just easier that way!)
! If they are not even, then code will need to be changed.

      real olr(NL,NM,NT)
      complex EE(NL,NT),CEE1(NL),CEE2(NT)
      real WSAVE1(4*NL+15),WSAVE2(4*NT+15)
      real PEE(NL+1,NT/2+1),totvar,globmean,totvar2,globmean2
      real sumPEE(NL+1,NT/2+1)
      real yrPEE(NL+1,NT/2+1)
      real totvarhalf
      real Rd,eif,s,tt,xx,xxm,ll,Pi,g,ps,k,ss(NL+1),ff(NT/2+1)
      real rand

! EE(n,t) = the initial (real) data set that later becomes the 
!        (complex) space-time spectrum
! PEE(n,t) = the (real) power spectrum
! ss() and ff() are arrays of wavenumbers and frequencies (cycles per day)
!          corresponding to PEE(,)
! totvar = the total variance of the dataset (about the global mean)
! globmean = the global time mean (average in space and time)
! NL = Number of longitudes in 'test' dataset
! NT = Number of times in 'test' dataset

! dimensional quantities are in MKS units.
! itr = No. of discrete samples per day
! s = planetary zonal wavenumber
! k = zonal wavenumber = 2 Pi s / L  (m-1)  
! ps = dimensional phase speed w.r.t. the earth 
! ll = distance around latitude circle (m)
! xx = longitude in degrees
! xxm = longitudinal distance from Greenwich in metres.
! id = time in days  (integer)
! tt = time in seconds (real)
! eif = dimensional frequency
!--------------------------------------------------------------------
! Variables for plotting routine
      integer il,iasf(13),li
      real sleft,sright,fftop,ffbot
      integer xmi,xma,mm,ymi,yma,nn
      parameter (sleft=-15.,sright=15.,fftop=.20,ffbot=0.)
! sleft/sright/fftop/ffbot are the zonal wavenumbers and frequencies
! (in cycles per day) to plot between.
      parameter (xmi=nint(sleft)+NL/2+1)
      parameter (xma=nint(sright)+NL/2+1)
      parameter (ymi=nint(ffbot*NT/(itr))+1)
      parameter (yma=nint(fftop*NT/(itr))+1)
      parameter (mm=xma-xmi+1,nn=yma-ymi+1)         
      real vpl,vpr,vpb,vpt,ul,ur,ub,ut,ffl,cfux,cfuy
      character*3 dd3
      character*2 dd2
      character*1 dd1
      character*80 label1
      data iasf/13*1/
!--------------------------------------------------------------------
!------------------------MAIN PROGRAM--------------------------------

!      print*,'xmi',xmi,' xma=',xma,' ymi=',ymi,' yma=',yma

      Pi = 3.1415926
      g = 9.81
      ll = 2.*Pi*6.37E06           ! at the equator
      totvar = 0.
      globmean = 0.

! Input 
!      open(1,file='homedir/level_2/wk99/variable
!     */data/seg96_over60_sym.gdat', access='direct',
!     *   form='unformatted',recl=NL*NM*linux_recl,status='old')

!      open(11,file='homedir/level_2/wk99/variable
!     */power/power.sym.gdat',access='direct',
!     *   form='unformatted',recl=(NL+1)*(NT/2+1)*linux_recl,
!     *status='unknown')

      status=NF90_OPEN(input,NF90_NOWRITE,ncid)
      IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_OPEN on '//input)
      status=NF90_INQ_VARID(ncid,'longitude',in_varids(1))
      IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_INQ_VARID on longitude')    
      status=NF90_INQ_VARID(ncid,'latitude',in_varids(2))
      IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_INQ_VARID on latitude')
      status=NF90_INQ_VARID(ncid,varname,in_varids(3))
      IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_INQ_VARID on '//varname)
      status=NF90_INQ_DIMID(ncid,'time',time_dimid)
      status=NF90_INQUIRE_DIMENSION(ncid,time_dimid,time_name,ntime)
      nyr=ntime/NT
      WRITE(6,*) 'Number of segments = ',nyr

      status=NF90_CREATE(output,0,ncid_out)
      IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_CREATE')
      status=NF90_DEF_DIM(ncid_out,'wavenumber',nl+1,out_dimids(1))
      IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
      status=NF90_DEF_DIM(ncid_out,'frequency',nt/2+1,out_dimids(2))
      IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
      status=NF90_DEF_VAR(ncid_out,'wavenumber',NF90_FLOAT,(/out_dimids(1)/),out_varids(1))
      IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
      status=NF90_DEF_VAR(ncid_out,'frequency',NF90_FLOAT,(/out_dimids(2)/),out_varids(2))
      IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
      status=NF90_DEF_VAR(ncid_out,'power',NF90_FLOAT,&
           (/out_dimids(1),out_dimids(2)/),out_varids(3))
      IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
      status=NF90_ENDDEF(ncid_out)
!
      do t = 1, NT/2+1
         do n = 1, NL+1
            yrPEE(n,t)=0.
         enddo
      enddo
!    
      do 2003 iy=1,nyr
         
         status=NF90_GET_VAR(ncid,in_varids(3),olr,start=(/1,1,(iy-1)*NT+1/),count=(/NL,nlat,NT/))
         !do t=1,NT
  	 !read(1,rec=nskip+(iy-1)*NT+t) ((olr(n,m,t),n=1,NL),m=1,nlat)
         !enddo
!
         do t = 1, NT/2+1
            do n = 1, NL+1
               sumPEE(n,t)=0.
            enddo
         enddo
         !
!! summed over lat (15S-15N)
!
         do 1000 m=1,nlat

            totvar=0.
            globmean=0.
            
            do 100 t=1,NT
               do 100 n=1,NL
                  EE(n,t) = olr(n,m,t)     ! symmetric / asymmetric
                  globmean = globmean + REAL(EE(n,t))/float(NT*NL)
100               continue
                  !      print*,'eif=',eif,' k=',k
                  
                  do t=1,NT
                     do n=1,NL
                        totvar = totvar + ((REAL(EE(n,t))-globmean)**2)/float(NT*NL)
                     enddo
                  enddo
                  !--------------------------------------------------------------------
                  !------------------COMPUTING SPACE-TIME SPECTRUM---------------------
                  
                  call cffti(NL,WSAVE1)            ! Initialize the (complex) FFT
                  do 200 t=1,NT
                     do 150 n=1,NL
                        CEE1(n) = EE(n,t)
                        ! CEE1(n) contains the grid values around a latitude circle
150                     continue
                        call cfftf(NL,CEE1,WSAVE1)
                        do 170 n=1,NL
                           EE(n,t) = CEE1(n)/float(NL)
170                        continue        
200                        continue
                           
                           ! Now the array EE(n,t) contains the Fourier coefficients (in planetary
                           ! wavenumber space) for each time.
                           
                           call cffti(NT,WSAVE2)      ! Initialize FFT for a different length. 
                           do 300 n=1,NL
                              do 250 t=1,NT
                                 CEE2(t) = EE(n,t)
                                 ! CEE2(t) contains a time-series of the coefficients for a single
                                 ! planetary zonal wavenumber
250                              continue
                                 call cfftf(NT,CEE2,WSAVE2)
                                 do 270 t=1,NT
                                    EE(n,t) = CEE2(t)/float(NT)
270                                 continue
300                                 continue
                                    
                                    ! Now the array EE(n,t) contains the (complex) space-time spectrum.
                                    
                                    ! Check some quantities such as the total variance.
                                    globmean2 = EE(1,1)
                                    totvar2 = 0.
                                    totvarhalf = 0.
                                    do t=1,NT
                                       do n=1,NL
                                          if(t.eq.1.and.n.eq.1) then
                                             continue
                                          else
                                             totvar2 = totvar2+(ABS(EE(n,t)))**2
                                             if(t.gt.1.and.t.lt.NT/2+1) then
                                                totvarhalf = totvarhalf+(ABS(SQRT(2.)*EE(n,t)))**2
                                             elseif(t.eq.1.or.t.eq.NT/2+1) then
                                                totvarhalf = totvarhalf+(ABS(EE(n,t)))**2
                                             endif
                                             !      N.b. summing for totvarhalf only works if NT is even.
                                          endif
                                       enddo
                                    enddo
                                    
                                    !     print*,'globmean =',globmean,'  globmean2=',globmean2
                                    !     print*,'totvar =',totvar,'  totvar2=',totvar2,' totvarhalf=',
                                    !    *     totvarhalf
                                    
                                    ! Create array PEE(NL+1,NT/2+1) which contains the (real) power spectrum.
                                    ! In this array, the westward moving waves will be from pn=1 to NL/2;
                                    ! The eastward waves will be for pn=NL/2+2 to NL+1.
                                    ! Positive frequencies will be from pt=1,NT/2+1.
                                    ! Information about zonal mean (wavenumber 0) will be for pn=NL/2+1.
                                    ! Information about time mean will be for pt=1.
                                    ! Information about the Nyquist Frequency is at pt=NT/2+1
                                    
                                    do pt=1,NT/2+1
                                       ff(pt)=float(pt-1)/float(NT)*float(itr) 
                                       do pn=1,NL+1
                                          ss(pn)=float(pn-1-NL/2)
                                          if(pn.le.NL/2) then
                                             n=NL/2+2-pn
                                             t=pt
                                          elseif(pn.eq.NL/2+1) then
                                             n=1
                                             t=pt
                                          elseif(pn.ge.NL/2+2) then
                                             n=pn-NL/2
                                             if (pt.eq.1) then
                                                t=pt
                                             else
                                                t=NT+2-pt
                                             endif
                                          endif
                                          PEE(pn,pt) = (CABS(EE(n,t)))**2
                                       enddo
                                    enddo
                                    
                                    do t = 1, NT/2+1
                                       do n = 1, NL+1
                                          if (mod(nlat,2).ne.0) then 
                                             if(m.eq.nlat/2+1) then
                                             else
                                                sumPEE(n,t)=sumPEE(n,t)+PEE(n,t)   ! summed over lat.
                                             endif
                                          else
                                             sumPEE(n,t)=sumPEE(n,t)+PEE(n,t)   ! summed over lat.
                                          endif
                                       enddo
                                    enddo
1000                                continue
                                    
                                    do t = 1, NT/2+1
                                       do n = 1, NL+1
                                          yrPEE(n,t)=yrPEE(n,t)+sumPEE(n,t)/float(nyr)
                                       enddo
                                    enddo
2003                                continue
                                    
                                    status=NF90_PUT_VAR(ncid_out,out_varids(1),ss)
                                    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_PUT_VAR on wavenumber')
                                    status=NF90_PUT_VAR(ncid_out,out_varids(2),ff)
                                    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_PUT_VAR on frequency')
                                    status=NF90_PUT_VAR(ncid_out,out_varids(3),yrPEE)
                                    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_PUT_VAR on power')

                                    status=NF90_CLOSE(ncid)
                                    status=NF90_CLOSE(ncid_out)
                                    !write (11,rec=1) yrPEE
                                    
                                  END PROGRAM
    
