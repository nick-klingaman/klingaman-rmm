 PROGRAM projectRMM
  USE netcdf
  USE nicks_power_toys, ONLY : error_handler,irapid_fill,rrapid_fill
  IMPLICIT NONE
 
! Originally written by Matthew Wheeler in ~2003
! Modified in February 2008 to work on ASCII input data (for the EOF sturctures)
! Modified by Nicholas Klingaman 12/2/2010 for NetCDF input data

! This program computes the RMM1 and RMM2 indices for any day (or series
! of days) that we have *all* sources of data for.
  !
! The required input data is daily-averaged anomalies of u850, u200, and olr.
! These fields must have already had the seasonal cycle removed as well as 
! the mean of the previous 120 days (as an estimate of interannual variability).
! They should be averages from 15S-15N on a 2.5 degree grid (i.e. 144 points).
! The input data should be centred on the longitudes 0deg, 2.5deg, ....., 357.5deg.
! 
! We must also read-in the file that contains the EOFs so that we can make
! the projection onto those structures.
!
! No missing data is allowed.

  !  Variables for input EOF file
  integer NX,TNX
  parameter (NX=144,TNX=3*NX)
  real eigval(2),eigvec(TNX,2),norm(TNX)
  integer imd,isp,it
  
  !  Variables for input data
  INTEGER :: olr_id,u200_id,u850_id,olr_varid,u200_varid,u850_varid,out_id,&
       out_dimids(2),out_varids(12)
    
  !  Other variables
!  integer fstyr,fstmo,fstda      !day to start producing RMM1 and RMM2
  INTEGER :: nyr,nmo,nda1,nda2,nda3
  INTEGER :: i,j,k,error,start,count
  INTEGER, parameter :: n_years=-1,n_days_per_year=-1,n_days_removed=-1
  REAL :: time_out(n_days_per_year)
  REAL :: olr(NX,n_days_per_year),u850(NX,n_days_per_year),u200(NX,n_days_per_year)
  REAL :: datmat(TNX),pc(2,n_days_per_year),amplitude(n_days_per_year),&
       contrib_olr(2,n_days_per_year),contrib_u200(2,n_days_per_year),contrib_u850(2,n_days_per_year)
  INTEGER :: phase(n_days_per_year)
  character*(*) directory,inpute,input1,input2,input3,output,olr_name,u200_name,u850_name

  ! INPUT EOF FILE    (for the spatial structure to project onto)
  parameter(inpute='/home/ss901165/src/mjo/WH04_EOFstruc.txt')
  
  ! INPUT DATA FILES  (daily anomalies with interannual removed)
  parameter(directory=directory)
  parameter(input1=directory//'/'//olr_file)
  parameter(input2=directory//'/'//u850_file)
  parameter(input3=directory//'/'//u200_file)
  parameter(output=directory//'/'//output_file)
  parameter(olr_name=olr_name)
  parameter(u200_name=u200_name)
  parameter(u850_name=u850_name)

  open(unit=11,file=inpute,status='old')
  
  error=NF90_OPEN(input1,NF90_NOWRITE,olr_id)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_OPEN on OLR file')
  error=NF90_OPEN(input2,NF90_NOWRITE,u850_id)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_OPEN on U850 file')
  error=NF90_OPEN(input3,NF90_NOWRITE,u200_id)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_OPEN on U200 file')
  
  error=NF90_INQ_VARID(olr_id,olr_name,olr_varid)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_INQ_VARID on OLR variable')
  error=NF90_INQ_VARID(u850_id,u850_name,u850_varid)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_INQ_VARID on U850 variable')
  error=NF90_INQ_VARID(u200_id,u200_name,u200_varid)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_INQ_VARID on U200 variable')
  
  !     OUTPUT FILE
  error=NF90_CREATE(output,NF90_CLOBBER,out_id)
  error=NF90_DEF_DIM(out_id,'time',NF90_UNLIMITED,out_dimids(1))
  error=NF90_DEF_VAR(out_id,'time',NF90_FLOAT,(/out_dimids(1)/),out_varids(1))
  error=NF90_DEF_VAR(out_id,'rmm1',NF90_FLOAT,(/out_dimids(1)/),out_varids(3))
  error=NF90_DEF_VAR(out_id,'rmm2',NF90_FLOAT,(/out_dimids(1)/),out_varids(4))
  error=NF90_DEF_VAR(out_id,'contrib_rmm1_olr',NF90_FLOAT,(/out_dimids(1)/),out_varids(5))
  error=NF90_DEF_VAR(out_id,'contrib_rmm1_u850',NF90_FLOAT,(/out_dimids(1)/),out_varids(6))
  error=NF90_DEF_VAR(out_id,'contrib_rmm1_u200',NF90_FLOAT,(/out_dimids(1)/),out_varids(7))
  error=NF90_DEF_VAR(out_id,'contrib_rmm2_olr',NF90_FLOAT,(/out_dimids(1)/),out_varids(8))
  error=NF90_DEF_VAR(out_id,'contrib_rmm2_u850',NF90_FLOAT,(/out_dimids(1)/),out_varids(9))
  error=NF90_DEF_VAR(out_id,'contrib_rmm2_u200',NF90_FLOAT,(/out_dimids(1)/),out_varids(10))
  error=NF90_DEF_VAR(out_id,'amplitude',NF90_FLOAT,(/out_dimids(1)/),out_varids(11))
  error=NF90_DEF_VAR(out_id,'phase',NF90_SHORT,(/out_dimids(1)/),out_varids(12))

  CALL RRAPID_FILL(time_out,1.0,0.,.FALSE.)

  error=NF90_PUT_ATT(out_id,out_varids(3),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(4),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(5),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(6),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(7),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(8),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(9),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(10),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(11),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(12),'missing_value',-9999)
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on missing_value')
  error=NF90_PUT_ATT(out_id,out_varids(3),'long_name','Real-time Multivariate Index 1 (Wheeler and Hendon 2004)')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
  error=NF90_PUT_ATT(out_id,out_varids(4),'long_name','Real-time Multivariate Index 2 (Wheeler and Hendon 2004)')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
  error=NF90_PUT_ATT(out_id,out_varids(5),'long_name','Contribution to RMM1 from OLR field')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
  error=NF90_PUT_ATT(out_id,out_varids(6),'long_name','Contribution to RMM1 from U850 field')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
  error=NF90_PUT_ATT(out_id,out_varids(7),'long_name','Contribution to RMM1 from U200 field')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
    error=NF90_PUT_ATT(out_id,out_varids(8),'long_name','Contribution to RMM2 from OLR field')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
  error=NF90_PUT_ATT(out_id,out_varids(9),'long_name','Contribution to RMM2 from U850 field')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
  error=NF90_PUT_ATT(out_id,out_varids(10),'long_name','Contribution to RMM2 from U200 field')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
  error=NF90_PUT_ATT(out_id,out_varids(11),'long_name','Amplitude of RMM indices (RMM1**2 + RMM2**2)**(1/2)')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
  error=NF90_PUT_ATT(out_id,out_varids(12),'long_name','Phase of the MJO from RMM indices (Wheeler and Hendon 2004)')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on long_name')
  error=NF90_PUT_ATT(out_id,NF90_GLOBAL,'calculation','Model data projected onto observed EOFs')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on calculation')
  error=NF90_PUT_ATT(out_id,out_varids(1),'units','days since forecast initialisation')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on calculation')
  error=NF90_PUT_ATT(out_id,out_varids(1),'cal','365-day')
  IF (error.NE.0) CALL ERROR_HANDLER(error,'NF90_PUT_ATT on calculation')

  error=NF90_ENDDEF(out_id)

  error=NF90_PUT_VAR(out_id,out_varids(1),time_out)

  !     Read in the EOF spatial structures
  !     ----------------------------------
  
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*) (eigval(imd),imd=1,2)  ! eigenvalues
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  do isp=1,TNX
     read(11,*) (eigvec(isp,imd),imd=1,2)  ! eigenvectors
  enddo
  read(11,*)
  read(11,*)
  do isp=1,TNX
     read(11,*) norm(isp)   ! normalization factor 
  enddo

  start=n_days_removed+1
  count=n_days_per_year-n_days_removed
  error=NF90_GET_VAR(olr_id,olr_varid,olr(:,n_days_removed+1:n_days_per_year),&
       start=(/1,start/),count=(/NX,count/))
  IF (error .NE. 0) CALL ERROR_HANDLER(error,'NF90_GET_VAR on OLR')
  error=NF90_GET_VAR(u850_id,u850_varid,u850(:,n_days_removed+1:n_days_per_year),&
       start=(/1,start/),count=(/NX,count/))
  IF (error .NE. 0) CALL ERROR_HANDLER(error,'NF90_GET_VAR on U850')
  error=NF90_GET_VAR(u200_id,u200_varid,u200(:,n_days_removed+1:n_days_per_year),&
       start=(/1,start/),count=(/NX,count/))
  IF (error .NE. 0) CALL ERROR_HANDLER(error,'NF90_GET_VAR on U200')
  pc(1,1:n_days_removed)=-9999
  pc(2,1:n_days_removed)=-9999
  contrib_olr(:,1:n_days_removed)=-9999
  contrib_u850(:,1:n_days_removed)=-9999
  contrib_u200(:,1:n_days_removed)=-9999
  amplitude(1:n_days_removed)=-9999
  phase(1:n_days_removed)=-9999        
  
  DO j=start,n_days_per_year
     
     !     -------
     !     Create extended vector of OLR, u850, and u200 data.
     DO i=1,TNX
        if(i.le.NX) then
           datmat(i)=olr(i,j)
        elseif(i.le.(NX*2)) then
           datmat(i)=u850(i-NX,j)
        elseif(i.le.(NX*3)) then
           datmat(i)=u200(i-(NX*2),j)
        endif
        !     Divide datmat by the "norm" as calculated in the EOF program.
        datmat(i)=datmat(i)/norm(i) ! I normalize by the "global" s.d. of each field.
     ENDDO
     !     -------
     
     !     Compute RMM1 and RMM2  (the first two normalized PCs)
     !     -----------------------------------------------------
     do imd=1,2
        pc(imd,j)=0.
        contrib_olr(imd,j)=0.
        contrib_u850(imd,j)=0.
        contrib_u200(imd,j)=0.
        do isp=1,TNX
           pc(imd,j)=pc(imd,j)+(eigvec(isp,imd)*datmat(isp))           
        enddo
        do isp=1,NX
           contrib_olr(imd,j)=contrib_olr(imd,j)+(eigvec(isp,imd)*datmat(isp))
           contrib_u850(imd,j)=contrib_u850(imd,j)+(eigvec(NX+isp,imd)*datmat(NX+isp))
           contrib_u200(imd,j)=contrib_u200(imd,j)+(eigvec(NX*2+isp,imd)*datmat(NX*2+isp))
        enddo
     ENDDO
     !     Now normalize (by EOF-calcualted s.d.) the newly calculated PCs
     do imd=1,2
        pc(imd,j)=pc(imd,j)/sqrt(eigval(imd))
        contrib_olr(imd,j)=contrib_olr(imd,j)/sqrt(eigval(imd))
        contrib_u850(imd,j)=contrib_u850(imd,j)/sqrt(eigval(imd))
        contrib_u200(imd,j)=contrib_u200(imd,j)/sqrt(eigval(imd))
     enddo
     
     amplitude(j)=(pc(1,j)**2+pc(2,j)**2)**(1./2.)
     IF (pc(1,j) .gt. 0) THEN
        IF (pc(2,j) .gt. 0) THEN 
           IF (ABS(pc(1,j)) .gt. ABS(pc(2,j))) THEN
              phase(j)=5
              !                 WRITE(6,*) 'j= ',j,' setting phase = ',phase(j),' pc1=',pc(1,j),'pc2=',pc(2,j)
           ELSE
              phase(j)=6
              !                 WRITE(6,*) 'j= ',j,' setting phase = ',phase(j),' pc1=',pc(1,j),'pc2=',pc(2,j)
           ENDIF
        ELSE
           IF (ABS(pc(1,j)) .gt. ABS(pc(2,j))) THEN
              phase(j)=4
              !                 WRITE(6,*) 'j= ',j,' setting phase = ',phase(j),' pc1=',pc(1,j),'pc2=',pc(2,j)
           ELSE
              phase(j)=3
              !                 WRITE(6,*) 'j= ',j,' setting phase = ',phase(j),' pc1=',pc(1,j),'pc2=',pc(2,j)
           ENDIF
        ENDIF
     ELSE
        IF (pc(2,j) .gt. 0) THEN
           IF (ABS(pc(1,j)) .gt. ABS(pc(2,j))) THEN
              phase(j)=8
              !                 WRITE(6,*) 'j= ',j,' setting phase = ',phase(j),' pc1=',pc(1,j),'pc2=',pc(2,j)
           ELSE
              phase(j)=7
              !                 WRITE(6,*) 'j= ',j,' setting phase = ',phase(j),' pc1=',pc(1,j),'pc2=',pc(2,j)
           ENDIF
        ELSE
           IF (ABS(pc(1,j)) .gt. ABS(pc(2,j))) THEN
              phase(j)=1
              !                 WRITE(6,*) 'j= ',j,' setting phase = ',phase(j),' pc1=',pc(1,j),'pc2=',pc(2,j)
           ELSE
              phase(j)=2
              !                 WRITE(6,*) 'j= ',j,' setting phase = ',phase(j),' pc1=',pc(1,j),'pc2=',pc(2,j)
           ENDIF
        ENDIF
     ENDIF
  ENDDO
  
  error=NF90_PUT_VAR(out_id,out_varids(3),pc(1,1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))
  error=NF90_PUT_VAR(out_id,out_varids(4),pc(2,1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))
  error=NF90_PUT_VAR(out_id,out_varids(5),contrib_olr(1,1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))
  error=NF90_PUT_VAR(out_id,out_varids(6),contrib_u850(1,1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))
  error=NF90_PUT_VAR(out_id,out_varids(7),contrib_u200(1,1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))
  error=NF90_PUT_VAR(out_id,out_varids(8),contrib_olr(2,1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))
  error=NF90_PUT_VAR(out_id,out_varids(9),contrib_u850(2,1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))
  error=NF90_PUT_VAR(out_id,out_varids(10),contrib_u200(2,1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))     
  error=NF90_PUT_VAR(out_id,out_varids(11),amplitude(1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))
  error=NF90_PUT_VAR(out_id,out_varids(12),phase(1:n_days_per_year),start=(/1/),count=(/n_days_per_year/))

  error=NF90_CLOSE(out_id)

!888 print*,'last input/output on ',nyr,nmo,' and day',nda1,' or ',nda2,' or ',nda3
END PROGRAM projectRMM
