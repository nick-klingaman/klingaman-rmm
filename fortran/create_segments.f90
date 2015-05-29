  PROGRAM create_segments
    USE netcdf
    USE nicks_subroutines
    IMPLICIT NONE
    
    integer nx,ny,ndays_per_year,nyears,ntot,nt,ii,jj,jump_y,overlap  
    parameter (nx=num_x, ny=num_y, ndays_per_year=ndpy, nyears=num_yr, &
         ntot=ndays_per_year*nyears, nt=num_s, jump_y=yoff, overlap=num_o)
    parameter (ii=num_x, jj=sel_y)
    real olr(nx,ny,nt), c_olr(ii,jj,nt), sym(ii,jj,nt), asy(ii,jj,nt),&
         longitude(nx),latitude(ny)
    real,allocatable :: time(:)
    CHARACTER*(*) directory,input,output_sym,output_asym,varname
    parameter(directory=directory)
    parameter(input=directory//'/'//input_file)
    parameter(output_sym=directory//'/'//output_file_sym)
    parameter(output_asym=directory//'/'//output_file_asym)
    parameter(varname=variable_name)
    INTEGER ncid,ncid_sym,in_varids(3),status,sym_dimids(3),sym_varids(4),&
         ncid_asym,asym_dimids(3),asym_varids(4),i,irec,iseg,it1,it2,j,j1,&
         j2,k,nskip,this_year,this_day,left

    status=NF90_OPEN(input,NF90_NOWRITE,ncid)
    IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_OPEN on '//input)
    status=NF90_INQ_VARID(ncid,'longitude',in_varids(1))
    IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_INQ_VARID on longitude')    
    status=NF90_INQ_VARID(ncid,'latitude',in_varids(2))
    IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_INQ_VARID on latitude')
    status=NF90_INQ_VARID(ncid,varname,in_varids(3))
    IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_INQ_VARID on '//varname)
    status=NF90_GET_VAR(ncid,in_varids(1),longitude)
    IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_GET_VAR on longitude')
    status=NF90_GET_VAR(ncid,in_varids(2),latitude)
    IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_GET_VAR on latitude')

    status=NF90_CREATE(output_sym,0,ncid_sym)
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_CREATE')
    status=NF90_DEF_DIM(ncid_sym,'longitude',nx,sym_dimids(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_DIM(ncid_sym,'latitude',jj,sym_dimids(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_DIM(ncid_sym,'time',NF90_UNLIMITED,sym_dimids(3))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_VAR(ncid_sym,'longitude',NF90_FLOAT,(/sym_dimids(1)/),sym_varids(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(ncid_sym,'latitude',NF90_FLOAT,(/sym_dimids(2)/),sym_varids(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(ncid_sym,'time',NF90_FLOAT,(/sym_dimids(3)/),sym_varids(3))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')   
    status=NF90_DEF_VAR(ncid_sym,'olr_sym',NF90_FLOAT,&
         (/sym_dimids(1),sym_dimids(2),sym_dimids(3)/),sym_varids(4))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_ENDDEF(ncid_sym)
    status=NF90_PUT_VAR(ncid_sym,sym_varids(1),longitude)
    status=NF90_PUT_VAR(ncid_sym,sym_varids(2),latitude(jump_y+1:jump_y+jj))
    
    status=NF90_CREATE(output_asym,0,ncid_asym)
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_CREATE')
    status=NF90_DEF_DIM(ncid_asym,'longitude',nx,asym_dimids(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_DIM(ncid_asym,'latitude',jj,asym_dimids(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_DIM(ncid_asym,'time',NF90_UNLIMITED,asym_dimids(3))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_VAR(ncid_asym,'longitude',NF90_FLOAT,(/asym_dimids(1)/),asym_varids(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(ncid_asym,'latitude',NF90_FLOAT,(/asym_dimids(2)/),asym_varids(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(ncid_asym,'time',NF90_FLOAT,(/asym_dimids(3)/),asym_varids(3))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')   
    status=NF90_DEF_VAR(ncid_asym,'olr_asym',NF90_FLOAT,&
         (/asym_dimids(1),asym_dimids(2),asym_dimids(3)/),asym_varids(4))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_ENDDEF(ncid_asym)
    status=NF90_PUT_VAR(ncid_asym,asym_varids(1),longitude)
    status=NF90_PUT_VAR(ncid_asym,asym_varids(2),latitude(jump_y+1:jump_y+jj))

!    open (1,file='homedir/level_1/variable/data/
!    *daily.anom.period.gdat',
!    *  form='unformatted',access='direct',recl=nx*ny*linux_recl,
!    *  status='old')
    
!    open (11,file='homedir/level_2/wk99/
!    *variable/data/seg96_over60_sym.gdat',
!    *  form='unformatted',access='direct',recl=ii*jj*linux_recl,
!    *  status='unknown')
!    open (12,file='homedir/level_2/wk99/
!    *variable/data/seg96_over60_asy.gdat',
!    *  form='unformatted',access='direct',recl=ii*jj*linux_recl,
!    *  status='unknown')
    
    iseg=0
    irec=0
    nskip=0    
    
    it1= nskip+1
    it2= it1+nt-1
    
10  continue	
!    jrec=0
!    do it=it1,it2
    !       jrec=jrec+1
!       read (1,rec=it) ((olr(i,j,jrec),i=1,nx),j=1,ny)
!    enddo
    
    !missing check for TRMM
    !       do jr = 1, jrec
    !        do j = 1, jj
    !        do i = 1, nx 
    !          j1 = j+jump_y
    !         if (olr(i,j1,jr).eq.undef) then
    !          print*,'ok'
    !           olr(i,j1,jr) = 0.5*(olr(i,j1,jr-1)+olr(i,j1,jr+1))
    !         endif
    !        enddo
    !        enddo
    !       enddo

    this_year=(it1/360)+1
    this_day=MOD(it1,360)
    IF (this_day+(it2-it1+1) .gt. 360) THEN 
       status=NF90_GET_VAR(ncid,in_varids(3),olr(:,:,1:360-this_day),start=(/1,1,this_year,this_day/),&
            count=(/nx,ny,1,360-this_day/))
       IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_GET_VAR on '//varname//' wrap 1')      
       left=(it2-it1+1)-(360-this_day)
       status=NF90_GET_VAR(ncid,in_varids(3),olr(:,:,360-this_day+1:(it2-it1+1)),&
            start=(/1,1,this_year+1,1/),count=(/nx,ny,1,left/))
       IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_GET_VAR on '//varname//' wrap 2')    
    ELSE
       status=NF90_GET_VAR(ncid,in_varids(3),olr,start=(/1,1,this_year,this_day/),count=(/nx,ny,1,it2-it1+1/))
       IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_GET_VAR on '//varname//' no wrap')
    ENDIF

    iseg=iseg+1
    !print*,iseg,'=iseg'

! sym & asy
    do k=1,nt
       do j=1,jj
          j1 = j+jump_y
          j2 = ny+1-j1
          
          do i=1,nx
             sym(i,j,k) = 0.5*(olr(i,j1,k)+olr(i,j2,k))
             asy(i,j,k) = 0.5*(olr(i,j1,k)-olr(i,j2,k))
          enddo ! i
          
       enddo ! j
    enddo ! k

! detrend	
    call detrend (sym,nx,jj,nt)
    call detrend (asy,nx,jj,nt)
    ! tapering both ends to zero
    !      apply a tapering to the first and last 5 days by multiplication
    !      by a segment of the cosine curve so that the ends of the series
    !      taper towards zero (mean has been removed anyway)!
    call tapertozero (sym,nx,jj,nt,1,nt,5)
    call tapertozero (asy,nx,jj,nt,1,nt,5)
    ! write_seg
    status=NF90_PUT_VAR(ncid_sym,sym_varids(4),sym,start=(/1,1,irec+1/),count=(/nx,jj,nt/))
    IF (status.NE.0) CALL ERROR_HANDLER(status,'NF90_PUT_VAR on symmetric')
    
    status=NF90_PUT_VAR(ncid_asym,asym_varids(4),asy,start=(/1,1,irec+1/),count=(/nx,jj,nt/))
    IF (status.NE.0) CALL ERROR_HANDLER(status,'NF90_PUT_VAR on anti-symmetric')
!	do it=1,nt
    irec=irec+nt
!	 write (11,rec=irec) ((sym(i,j,it),i=1,ii),j=1,jj)
!	 write (12,rec=irec) ((asy(i,j,it),i=1,ii),j=1,jj)
!	enddo
	
    it1=it2-overlap+1
    it2=it1+nt-1
    !print*,'it1=',it1,'it2=',it2
    if(it2.gt.ntot) goto 20
    goto 10
    
20  continue
    
    allocate(time(irec))
    CALL RRAPID_FILL(time,1.0,0.5,.TRUE.)
    status=NF90_PUT_VAR(ncid_asym,asym_varids(3),time)
    status=NF90_PUT_VAR(ncid_sym,sym_varids(3),time)
    status=NF90_CLOSE(ncid)
    IF (status .ne. 0) CALL ERROR_HANDLER(status,'NF90_CLOSE on '//input)
    status=NF90_CLOSE(ncid_sym)
    status=NF90_CLOSE(ncid_asym)
    print*,iseg,'= total 96-day segments'
   
    stop
  end PROGRAM create_segments

!-----------------------------------------------------------------------
      subroutine detrend (tmp2,nl,nlat,nt)
      implicit none
      integer NL,nlat,NT,n,m,t, tmean,sum
      real tmp2(NL,nlat,NT),y(NT)
      real a, b
        do 10 n=1,nl
        do 10 m=1,nlat
         do 20 t=1,nt
           y(t)=tmp2(n,m,t)
  20     continue
        tmean=sum(y,nt)
         call reg (y,nt,a,b)
         do 30 t=1,nt
           tmp2(n,m,t)=tmp2(n,m,t)-(a+b*float(t-1))
!           tmp2(n,m,t)=tmp2(n,m,t)-tmean/float(nt)
  30     continue
  10  continue
      return
    end subroutine detrend
!-----------------------------------------------------------------------
      subroutine reg (y,n,a,b)
      implicit none
      integer i,n
      real a,b,xbar,ybar,sum
      real y(n)
      real xy(n),xx(n)
!
      xbar=n*(n+1)/2.
      ybar=sum (y,n)/float(n)
      do 100 i=1,n
         xy(i)=(i-xbar)*(y(i)-ybar)
         xx(i)=(i-xbar)**2
100      continue
         b=sum(xy,n)/sum(xx,n)
         a=ybar-b*xbar         
         return
      end subroutine reg
!-----------------------------------------------------------------------
      function sum (x,n)
      implicit none
      integer n,i
      real x(n),sumx, sum
      sumx=0.
      do 100 i=1,n
       sumx=sumx+x(i)
100    continue
       sum=sumx
       return
    end function sum

!-----------------------------------------------------------------
    subroutine tapertozero(ts,imax,jmax,N,nmi,nn,tp)
      implicit none
      
! "taper" the first and last 'tp' members of 'ts' by multiplication by
! a segment of the cosine curve so that the ends of the series
! taper toward zero. This satisfies the
! periodic requirement of the FFT.
! Only the data from nmi to nmi+nn-1 is deemed useful, and
! the rest is set to zero. (There are nn points of useful data)
       integer i,j,N,tp,nmi,nn, imax, jmax, ii, jj
       real Pi,ts(imax,jmax,N)
       parameter (Pi=3.1415926)   !tp is number to taper on each end.

       if (N.lt.25.or.nn.lt.25) then
        print*,'No use doing the tapering if less than 25 values!'
        STOP
       endif
       if (nmi+nn-1.gt.N) then
        print*,'ERROR as nmi+nn-1 > N'
        STOP
       endif

       do jj = 1, jmax
       do ii = 1, imax

       do i=1,N
        j=i-nmi+1
        if (j.le.0.or.j.gt.nn) then
         ts(ii,jj,i)=0.
        elseif (j.le.tp) then
         ts(ii,jj,i)= ts(ii,jj,i)*.5*(1.-COS((j-1)*Pi/float(tp)))
        elseif (j.gt.(nn-tp).and.j.le.nn) then
         ts(ii,jj,i)= ts(ii,jj,i)*.5*(1.-COS((nn-j)*Pi/float(tp)))
        else
         ts(ii,jj,i)= ts(ii,jj,i)
        endif
       enddo

       enddo
       enddo

       return
     end subroutine tapertozero
