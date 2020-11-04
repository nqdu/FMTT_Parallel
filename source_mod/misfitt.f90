!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to calculate an estimate of model 
! roughness and variance by dicing up a cubic B-spline grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM misfit
IMPLICIT NONE
INTEGER :: i,j,k,l,m,n,i1,j1,k1
INTEGER :: nvr,nvt,nvp,conr,cont,conp
INTEGER :: ndr,ndt,ndp,nnr,nnt,nnp,str,stt,stp
INTEGER :: checkstat
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: gor,got,gop,gsr,gst,gsp,u,mvar,mrough
REAL(KIND=i10) :: rgsr,rgst,rgsp,rdm,sumi,sumj,sumk,earth
REAL(KIND=i10) :: rdm1,sumi1,sumj1,sumk1,dp,dt,dr,ri,risti
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi,wi
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: velsm,velrm,velr
REAL(KIND=i10), PARAMETER :: pi=3.14159265359
CHARACTER (LEN=25) :: rmfile,smfile
!
! nvr = number of vertices in radial direction
! nvt = number of vertices in theta (N-S) direction
! nvp = number of vertices in phi (E-W) direction
! gor = grid origin in radius
! got = grid origin in theta (N-S)
! gop = grid origin in phi (E-W)
! gsr = grid separation in radius
! gst = grid separation in theta (N-S)
! gsp = grid separation in phi (E-W)
! smfile = solution model file 
! rmfile = reference model file
! velr = velocity at refined grid point
! velsm = velocity at spline vertex for solution model
! velrm = velocity at spline vertex for reference model
! ndr,ndt,ndp = node dicing level in r,theta,phi
! nnr,nnt,nnp = number of diced nodes in r,theta,phi
! u = Cubic spline independent variable
! ui,vi,wi = Cubic spline basis functions
! rgsr,rgst,rgsp = Refined node spacing in r,theta,phi
! sumi,sumj,sumk = Summation variables for constructing spline
! conr,cont,conp = Counters for refining grid
! str,stt,stp = Refined grid location in r,theta,phi
! mvar = Model variance
! mrough = Model roughness
! sumi1,sumj1,sumk1= Summation variables for reference spline
! dp,dt,dr = Contributions to roughness Laplacian
! ri,risti = Denomenators for difference operator
! earth = Earth radius
!
OPEN(UNIT=10,FILE='misfitt.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a25)')smfile
READ(10,'(a25)')rmfile
READ(10,*)ndr,ndt,ndp
READ(10,*)earth
CLOSE(10)
!
! Read in B-spline grid of solution model
!
OPEN(UNIT=20,FILE=smfile,STATUS='old')
READ(20,*)nvr,nvt,nvp
READ(20,*)gor,got,gop
READ(20,*)gsr,gst,gsp
ALLOCATE(velsm(0:nvr+1,0:nvt+1,0:nvp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL velv'
ENDIF
DO i=0,nvp+1
   DO j=0,nvt+1
      DO k=0,nvr+1
         READ(20,*)velsm(k,j,i)
      ENDDO
   ENDDO
ENDDO
CLOSE(20)
!
! Read in B-spline grid of reference model
!
OPEN(UNIT=20,FILE=rmfile,STATUS='old')
READ(20,*)nvr,nvt,nvp
READ(20,*)gor,got,gop
READ(20,*)gsr,gst,gsp
got=(90.0-got)*pi/180.0
gop=gop*pi/180.0
gst=gst*pi/180.0
gsp=gsp*pi/180.0
ALLOCATE(velrm(0:nvr+1,0:nvt+1,0:nvp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL velv'
ENDIF
DO i=0,nvp+1
   DO j=0,nvt+1
      DO k=0,nvr+1
         READ(20,*)velrm(k,j,i)
      ENDDO
   ENDDO
ENDDO
CLOSE(20)
! Calculate total numer of refined nodes in r,theta,phi
! and the refined grid spacing.
!
nnr=(nvr-1)*ndr+1
nnt=(nvt-1)*ndt+1
nnp=(nvp-1)*ndp+1
rgsr=gsr/ndr
rgst=gst/ndt
rgsp=gsp/ndp
ALLOCATE(velr(nnr,nnt,nnp), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL velr'
ENDIF
!
! Calculate the values of the basis functions
!
ALLOCATE(ui(nvr+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL ui'
ENDIF
DO i=1,ndr+1
   u=ndr
   u=(i-1)/u
   ui(i,1)=(1.0-u)**3/6.0
   ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   ui(i,4)=u**3/6.0
ENDDO
ALLOCATE(vi(nvt+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL vi'
ENDIF
DO i=1,ndt+1
   u=ndt
   u=(i-1)/u
   vi(i,1)=(1.0-u)**3/6.0
   vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   vi(i,4)=u**3/6.0
ENDDO
ALLOCATE(wi(nvp+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM cbsvel: REAL wi'
ENDIF
DO i=1,ndp+1
   u=ndp
   u=(i-1)/u
   wi(i,1)=(1.0-u)**3/6.0
   wi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   wi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   wi(i,4)=u**3/6.0
ENDDO
!
! Calculate velocity values on refined grid
!
mvar=0.0
DO i=1,nvp-1
   conp=ndp
   IF(i==nvp-1)conp=ndp+1
   DO j=1,nvt-1
      cont=ndt
      IF(j==nvt-1)cont=ndt+1
      DO k=1,nvr-1
         conr=ndr
         IF(k==nvr-1)conr=ndr+1
         DO l=1,conp
            stp=ndp*(i-1)+l
            DO m=1,cont
               stt=ndt*(j-1)+m
               DO n=1,conr
                  str=ndr*(k-1)+n
                  sumi=0.0
                  sumi1=0.0
                  DO i1=1,4
                     sumj=0.0
                     sumj1=0.0
                     DO j1=1,4
                        sumk=0.0
                        sumk1=0.0
                        DO k1=1,4
                           rdm=ui(n,k1)*velsm(k-2+k1,j-2+j1,i-2+i1)
                           rdm1=ui(n,k1)*velrm(k-2+k1,j-2+j1,i-2+i1)
                           sumk=sumk+rdm
                           sumk1=sumk1+rdm1
                        ENDDO
                        sumj=sumj+vi(m,j1)*sumk
                        sumj1=sumj1+vi(m,j1)*sumk1
                     ENDDO
                     sumi=sumi+wi(l,i1)*sumj
                     sumi1=sumi1+wi(l,i1)*sumj1
                  ENDDO
                  velr(str,stt,stp)=sumi-sumi1
                  mvar=mvar+(sumi-sumi1)**2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
mvar=mvar/(nnp*nnt*nnr-1.0)
!
! Now estimate model roughness
!
mrough=0.0
DO i=2,nnp-1
   DO j=2,nnt-1
      DO k=2,nnr-1
         ri=gor-(k-1)*rgsr+earth
         risti=ri*sin(got+(j-1)*rgst)
         dp=velr(k,j,i+1)-2.0*velr(k,j,i)+velr(k,j,i-1)
         dp=dp/((risti*rgsp)**2)
         dt=velr(k,j+1,i)-2.0*velr(k,j,i)+velr(k,j-1,i)
         dt=dt/((ri*rgst)**2)
         dr=velr(k+1,j,i)-2.0*velr(k,j,i)+velr(k-1,j,i)
         dr=dr/(rgsr**2)
         mrough=mrough+ABS(dp)+ABS(dt)+ABS(dr)
      ENDDO
   ENDDO
ENDDO
mrough=mrough/((nnp-2)*(nnt-2)*(nnr-2))
WRITE(6,*)'Model variance in (km/s)**2 is ',mvar
WRITE(6,*)'Model roughness in (kms)**(-1) is ',mrough
DEALLOCATE(ui,vi,wi, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM cbsvel: REAL ui,vi,wi'
ENDIF
DEALLOCATE(velsm,velrm,velr, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM cbsvel: REAL velsm,velrm,velr'
ENDIF
STOP
END PROGRAM misfit
