  PROGRAM wk99_smoothing
    USE nicks_subroutines
    USE netcdf
    IMPLICIT NONE

    ! Modified by Nick Klingaman to use netCDF input and output.
    ! 16/4/14
    
    integer nl,nt
    parameter (NL=num_x,NT=num_s)
    integer pt, pn, i, j, n, t

    CHARACTER*(*) directory,input_sym,input_asym,&
         output_sym,output_asym,output_bkgrnd,varname_sym,varname_asym         
    CHARACTER(4) time_name
    parameter(directory=directory)
    parameter(input_sym=directory//'/'//input_file_sym)
    parameter(input_asym=directory//'/'//input_file_asym)
    parameter(output_sym=directory//'/'//output_file_sym)
    parameter(output_asym=directory//'/'//output_file_asym)
    parameter(output_bkgrnd=directory//'/'//output_file_bkgrnd)
    parameter(varname_sym=variable_name_sym)
    parameter(varname_asym=variable_name_asym)    
    INTEGER ncid_sym,ncid_asym,varid_sym,varid_asym,&
         out_varids_sym(3),out_varids_asym(3),out_varids_bkgrnd(3),&
         out_dimids_sym(2),out_dimids_asym(2),out_dimids_bkgrnd(2),&
         out_ncid_sym,out_ncid_asym,out_ncid_bkgrnd,status
    
    real yrPEE1(NL+1,NT/2+1), yrPEE2(NL+1,NT/2+1)
    real sumPEE(NL+1,NT/2+1), log_sum(NL+1,NT/2+1)
    real log_sym(NL+1,NT/2+1), log_asy(NL+1,NT/2+1)
    real norm_sym(NL+1,NT/2+1), norm_asy(NL+1,NT/2+1)
    real ff(NT/2+1), ss(NL+1), mm(NL+1), frq(NT/2+1)
    real dmiss
    data dmiss /-9999./
    logical olr
    data olr /.false./
    
    do pt=1,NT/2+1
       ff(pt)=float(pt-1)/float(NT)
       do pn=1,NL+1
          ss(pn)=float(pn-1-NL/2)
       enddo
    enddo
    
    status=NF90_OPEN(input_sym,NF90_NOWRITE,ncid_sym)
    IF (status.NE.0) CALL ERROR_HANDLER(status,'NF90_OPEN on symmetric input file')
    status=NF90_OPEN(input_asym,NF90_NOWRITE,ncid_asym)
    IF (status.NE.0) CALL ERROR_HANDLER(status,'NF90_OPEN on asymmetric input file')    
    status=NF90_INQ_VARID(ncid_sym,varname_sym,varid_sym)
    IF (status.NE.0) CALL ERROR_HANDLER(status,'NF90_INQ_VARID on symmetric variable')
    status=NF90_INQ_VARID(ncid_asym,varname_asym,varid_asym)
    IF (status.NE.0) CALL ERROR_HANDLER(status,'NF90_INQ_VARID on asymmetric variable')

    status=NF90_GET_VAR(ncid_sym,varid_sym,yrPEE1)
    IF (status.NE.0) CALL ERROR_HANDLER(status,'NF90_GET_VAR on symmetric power')
    status=NF90_GET_VAR(ncid_asym,varid_asym,yrPEE2)
    IF (status.NE.0) CALL ERROR_HANDLER(status,'NF90_GET_VAR on asymmetric power')
    
    status=NF90_CREATE(output_sym,0,out_ncid_sym)
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_CREATE')
    status=NF90_DEF_DIM(out_ncid_sym,'wavenumber',nl+1,out_dimids_sym(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_DIM(out_ncid_sym,'frequency',nt/2+1,out_dimids_sym(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_VAR(out_ncid_sym,'wavenumber',NF90_FLOAT,(/out_dimids_sym(1)/),out_varids_sym(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(out_ncid_sym,'frequency',NF90_FLOAT,(/out_dimids_sym(2)/),out_varids_sym(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(out_ncid_sym,'power',NF90_FLOAT,&
         (/out_dimids_sym(1),out_dimids_sym(2)/),out_varids_sym(3))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_ENDDEF(out_ncid_sym)
        
    status=NF90_CREATE(output_asym,0,out_ncid_asym)
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_CREATE')
    status=NF90_DEF_DIM(out_ncid_asym,'wavenumber',nl+1,out_dimids_asym(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_DIM(out_ncid_asym,'frequency',nt/2+1,out_dimids_asym(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_VAR(out_ncid_asym,'wavenumber',NF90_FLOAT,(/out_dimids_asym(1)/),out_varids_asym(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(out_ncid_asym,'frequency',NF90_FLOAT,(/out_dimids_asym(2)/),out_varids_asym(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(out_ncid_asym,'power',NF90_FLOAT,&
         (/out_dimids_asym(1),out_dimids_asym(2)/),out_varids_asym(3))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_ENDDEF(out_ncid_asym)

    status=NF90_CREATE(output_bkgrnd,0,out_ncid_bkgrnd)
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_CREATE')
    status=NF90_DEF_DIM(out_ncid_bkgrnd,'wavenumber',nl+1,out_dimids_bkgrnd(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_DIM(out_ncid_bkgrnd,'frequency',nt/2+1,out_dimids_bkgrnd(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_DIM')
    status=NF90_DEF_VAR(out_ncid_bkgrnd,'wavenumber',NF90_FLOAT,(/out_dimids_bkgrnd(1)/),out_varids_bkgrnd(1))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(out_ncid_bkgrnd,'frequency',NF90_FLOAT,(/out_dimids_bkgrnd(2)/),out_varids_bkgrnd(2))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_DEF_VAR(out_ncid_bkgrnd,'power',NF90_FLOAT,&
         (/out_dimids_bkgrnd(1),out_dimids_bkgrnd(2)/),out_varids_bkgrnd(3))
    IF (status.ne.0) CALL ERROR_HANDLER(status,'NF90_DEF_VAR')
    status=NF90_ENDDEF(out_ncid_bkgrnd)

    status=NF90_PUT_VAR(out_ncid_sym,out_varids_sym(1),ss)
    status=NF90_PUT_VAR(out_ncid_sym,out_varids_sym(2),ff)
    status=NF90_PUT_VAR(out_ncid_asym,out_varids_asym(1),ss)
    status=NF90_PUT_VAR(out_ncid_asym,out_varids_asym(2),ff)
    status=NF90_PUT_VAR(out_ncid_bkgrnd,out_varids_bkgrnd(1),ss)
    status=NF90_PUT_VAR(out_ncid_bkgrnd,out_varids_bkgrnd(2),ff)
    
    !open (1,file='homedir/level_2/wk99/variable/
    !*power/power.sym.gdat', access='direct',
    !*   form='unformatted',recl=(NL+1)*(NT/2+1)*linux_recl,
    !*status='unknown')
    
    !open (2,file='homedir/level_2/wk99/variable/
    !*power/power.asy.gdat', access='direct',
    !*   form='unformatted',recl=(NL+1)*(NT/2+1)*linux_recl,
    !*status='unknown')

    !open (11,file='homedir/level_2/wk99/variable/
    !*power/norm.sym.gdat',access='direct',
    !*   form='unformatted',recl=(NL+1)*(NT/2+1)*linux_recl,
    !*status='unknown')
    
    !open (12,file='homedir/level_2/wk99/variable/
    !*power/norm.asy.gdat',access='direct',
    !*   form='unformatted',recl=(NL+1)*(NT/2+1)*linux_recl,
    !*status='unknown')
    
    !open (13,file='homedir/level_2/wk99/variable/
    !*power/back.gdat',access='direct',
    !*   form='unformatted',recl=(NL+1)*(NT/2+1)*linux_recl,
    !*status='unknown')

    !read (1,rec=1) yrPEE1
    !read (2,rec=1) yrPEE2

    do j=1,nt/2+1
       do i=1,nl+1
	  sumPEE(i,j)=0.5*(yrPEE1(i,j)+yrPEE2(i,j))
   log_sym(i,j) = log10(yrPEE1(i,j))
   log_asy(i,j) = log10(yrPEE2(i,j))
enddo
enddo

! Remove the satellite aliases for olr
if (olr) then
   
   do j = 1, nt/2+1
      do i = 1, nl+1
         if((ss(i).ge.13..and.ss(i).le.15.).and.&
              ff(j).gt.0.09.and.ff(j).lt.0.14) then
            sumPEE(i,j)=-9999.
         elseif((ss(i).ge.13..and.ss(i).le.15.).and.&
              ff(j).ge.0.20.and.ff(j).lt.0.23) then
            sumPEE(i,j)=-9999.
         endif
      enddo
   enddo
endif

! smoothing
! This smoothing DOES include wavenumber zero
do t=1, nt/2+1
   do n=1, nl+1
      mm(n)=sumPEE(n,t)
   enddo
   
   if (ff(t).lt.0.1) then
      do i=1,5
         call smooth121(mm,nl+1,nl+1)
      enddo
   elseif (ff(t).lt.0.2) then
      do i=1,10
         call smooth121(mm,nl+1,nl+1)
      enddo
   elseif (ff(t).lt.0.3) then
      do i=1,20
         call smooth121(mm,nl+1,nl+1)
      enddo
   else
      do i=1,40
         call smooth121(mm,nl+1,nl+1)
      enddo
   endif
   
   do n=1, nl+1
      sumPEE(n,t)=mm(n)
   enddo
enddo

do n=1, nl+1
   do t=1,NT/2-1
      frq(t)=sumPEE(n,t+1)
   enddo
   do i=1,10
      call smooth121(frq,NT/2+1,NT/2-1)
   enddo
   do t=1,NT/2-1
      sumPEE(n,t+1)=frq(t)
   enddo
enddo


!write (13,rec=1) ((sumPEE(i,j),i=1,nl+1),j=1,nt/2+1)
status=NF90_PUT_VAR(out_ncid_bkgrnd,out_varids_bkgrnd(3),sumPEE)
IF (status .NE. 0) CALL ERROR_HANDLER(status,'NF90_PUT_VAR on background power')

do i=1,nl+1
   do j=1,nt/2+1
      if(sumPEE(i,j).ne.dmiss) then
         log_sum(i,j)=log10(sumPEE(i,j))
      else
         log_sum(i,j)=dmiss
      endif
      if (ff(j).eq.0) then
         log_sum(i,j) = dmiss
      endif
      
      if (log_sum(i,j).ne.dmiss) then
         norm_sym(i,j) = log_sym(i,j)-log_sum(i,j)
         norm_asy(i,j) = log_asy(i,j)-log_sum(i,j)
      else
         norm_sym(i,j) = dmiss
         norm_asy(i,j) = dmiss
      endif
   enddo
enddo

!write (11,rec=1) ((norm_sym(i,j),i=1,nl+1),j=1,nt/2+1)
!write (12,rec=1) ((norm_asy(i,j),i=1,nl+1),j=1,nt/2+1)
status=NF90_PUT_VAR(out_ncid_sym,out_varids_sym(3),norm_sym)
IF (status .NE. 0) CALL ERROR_HANDLER(status,'NF90_PUT_VAR on symmetric power')
status=NF90_PUT_VAR(out_ncid_asym,out_varids_asym(3),norm_asy)
IF (status .NE. 0) CALL ERROR_HANDLER(status,'NF90_PUT_VAR on asymmetric power')

status=NF90_CLOSE(ncid_sym)
status=NF90_CLOSE(ncid_asym)
status=NF90_CLOSE(out_ncid_bkgrnd)
status=NF90_CLOSE(out_ncid_sym)
status=NF90_CLOSE(out_ncid_asym)

999 stop
end PROGRAM wk99_smoothing

!-----------------------------------------------------------------
      subroutine smooth121(vv,vn,nn)
      implicit none
! Smooths vv by passing it through a 1-2-1 filter.
! The first and last points are given 3-1 (1st) or 1-3 (last)
! weightings (Note that this conserves the total sum).
! The routine also skips-over missing data (assigned to be
! a value of 1.E36).
! There are 'nn' pieces of useful information, which may be less
! than or equal to 'vn'.

      integer nn,i,vn
      real spv,vv(vn),dum(5000)

      if (nn.gt.5000) then
       print*,'need to increase 5000 in smooth121.f'
       STOP
      endif

      spv=-9999.
      i=0
 10   continue
      i=i+1
      if (vv(i).eq.spv) then
       dum(i)=spv
      elseif(i.eq.1.or.vv(i-1).eq.spv) then
       dum(i)=(3.*vv(i)+vv(i+1))/4.
      elseif(i.eq.nn.or.vv(i+1).eq.spv) then
       dum(i)=(vv(i-1)+3.*vv(i))/4.
      else
       dum(i)=(1.*vv(i-1)+2.*vv(i)+1.*vv(i+1))/4.
      endif
      if (i.ne.nn) goto 10

      do i=1,nn
       vv(i)=dum(i)
      enddo
      RETURN
    END subroutine smooth121
  
    
