!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program calculates the RMS traveltime residual and
! variance of the current model.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module valid_rays
!nqdu
! self-defined struct
type,public :: FmttPair
   integer :: ns
   integer,ALLOCATABLE  :: valid_rays(:)
end type
type(FmttPair),ALLOCATABLE    :: pairs(:)

contains

subroutine FmttPair_dealloc(pair)
   implicit none
   type(FmttPair) :: pair

   DEALLOCATE(pair%valid_rays)

end subroutine FmttPair_dealloc
!===============================================
end module valid_rays

PROGRAM residual
use valid_rays
IMPLICIT NONE
INTEGER :: i,ii,j,idm,rmtr,isum,itern,stnm,invst
INTEGER :: nr,ns,nt,iscn,checkstat,idum,nra
INTEGER, DIMENSION (:), ALLOCATABLE :: ntsrc,ots
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: mtr,rsum,var,rdum
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: stt,srcm
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: mtim,rtim,otim
CHARACTER (LEN=27) :: mtimes,rtimes,otimes,sources,receivers,subinv
CHARACTER (LEN=27) :: stfile,itfile

!
! rtimes = input file of reference times
! otimes = output file of traveltime residuals
! sources = file containing source information
! receivers = file containing receiver information
! mtimes = file containing model traveltimes
! subinv = file containing subspace parameters
! ns= number of sources
! nr = number of receivers
! nt = number of traveltimes
! mtr = set of traveltime residuals
! var = variance
! rmtr = remove model residual traveltimes (0=no,1=yes)
! stfile = File containing station terms
! itfile = File containing inversion iteration number
! invst = Invert for station terms (0=no, 1=yes)
! stt = Station traveltime term
! itern = Inversion iteration number
! stnm = Station number
! srcm = model traveltime residual mean for each source
! ntsrc = Number of traveltimes for each source
! iscn = Source number
! nra = Number of receiver arrays
! mtim = Model traveltime array
! rtim = Reference traveltime array
! otim = Observed traveltime array
! ots = Observed traveltime status (0 or 1)
!
! Read in the input parameters
!
OPEN(UNIT=10,FILE='residualst.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a27)')otimes
READ(10,'(a27)')rtimes
READ(10,'(a27)')subinv
READ(10,'(a27)')stfile
READ(10,'(a27)')itfile
READ(10,'(a27)')sources
READ(10,'(a27)')receivers
READ(10,'(a27)')mtimes
CLOSE(10)

! read valid_rays
OPEN(UNIT=50,FILE=sources,STATUS='old')
READ(50,*)nra
open(99,file='valid_rays.dat')
ALLOCATE(pairs(nra))
do ii=1,nra
   read(50,*)idum
   ALLOCATE(pairs(ii)%valid_rays(idum))
   pairs(ii)%ns = idum
   do i=1,idum
      read(50,*)
      read(50,*)
      read(99,*)pairs(ii)%valid_rays(i)
   enddo
enddo
close(50)
close(99)

!
! Determine the number of receiver arrays
!
OPEN(UNIT=50,FILE=sources,STATUS='old')
READ(50,*)nra
OPEN(UNIT=60,FILE=receivers,STATUS='old')
READ(60,*)idum
IF(idum.NE.nra)THEN
   WRITE(6,*)'ERROR!!!'
   WRITE(6,*)'Source and receiver files are'
   WRITE(6,*)'inconsistent!!!!'
   WRITE(6,*)'First line of each file should'
   WRITE(6,*)'be identical!!!'
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
!
! Determine whether to remove residual mean and
! whether station terms are to be included
!
OPEN(UNIT=10,FILE=subinv,STATUS='old')
DO i=1,8
   READ(10,*)
ENDDO
READ(10,*)rmtr
DO i=1,8
   READ(10,*)
ENDDO
READ(10,*)invst
CLOSE(10)
OPEN(UNIT=10,FILE=mtimes,STATUS='old')
OPEN(UNIT=20,FILE=rtimes,STATUS='old')
OPEN(UNIT=30,FILE=otimes,STATUS='old')
IF(invst.EQ.1)THEN
   OPEN(UNIT=40,FILE=itfile,STATUS='old')
   READ(40,*)itern
   CLOSE(40)
   IF(itern.GT.1)THEN
      OPEN(UNIT=40,FILE=stfile,STATUS='old')
   ENDIF
ENDIF
rsum=0.0
isum=0
DO ii=1,nra
   READ(50,*)ns
   READ(60,*)nr
   DO i=1,ns
      READ(50,*)
      READ(50,*)
   ENDDO
   DO i=1,nr
      READ(60,*)
   ENDDO
   nt=ns*nr
   !
   ! If station terms are applied, then read in
   ! residual.
   !
   IF(invst.EQ.1)THEN
      ALLOCATE(stt(nr))
      stt=0.0
      IF(itern.GT.1)THEN
         DO i=1,nr
            READ(40,*)stt(i)
         ENDDO
      ENDIF
   ENDIF
   ALLOCATE(mtim(nt))
   ALLOCATE(rtim(nt))
   ALLOCATE(otim(nt))
   ALLOCATE(ots(nt))
   DO i=1,ns
      do j=1,nr
         READ(10,*)mtim((i-1)*nr+j)
         READ(20,*)rtim((i-1)*nr + j)
         READ(30,*)ots((i-1) * nr + j),otim((i-1) * nr + j)
         ots((i-1)*nr + j) =  ots((i-1)*nr + j) * pairs(ii)%valid_rays(i)
      enddo
   enddo
!
!  Compute average residual if required
!
   IF(rmtr.EQ.1)THEN
      ALLOCATE(srcm(ns),ntsrc(ns))
      srcm=0.0
      ntsrc=0
      DO i=1,nt
         IF(invst.EQ.1)THEN
            stnm=MOD(i,nr)
            IF(stnm.EQ.0)stnm=nr
         ENDIF
         rdum=mtim(i)
         IF(invst.EQ.1)rdum=rdum+stt(stnm) 
         IF(ots(i).EQ.1)THEN
            iscn=INT((i-1)/nr)+1
            srcm(iscn)=srcm(iscn)+rdum-rtim(i)
            ntsrc(iscn)=ntsrc(iscn)+1
         ENDIF
      ENDDO
      DO i=1,ns
         srcm(i)=srcm(i)/REAL(ntsrc(i))
      ENDDO
   ENDIF
!
!  Take difference between model and reference traveltimes
!  and then subtract the observed residuals.
!
   DO i=1,nt
      IF(invst.EQ.1)THEN
         stnm=MOD(i,nr)
         IF(stnm.EQ.0)stnm=nr
      ENDIF
      rdum=mtim(i)
      IF(invst.EQ.1)rdum=rdum+stt(stnm)
      IF(ots(i).EQ.1)THEN
         IF(rmtr.EQ.1)THEN
            iscn=INT((i-1)/nr)+1
            mtr=rdum-rtim(i)-otim(i)-srcm(iscn)
         ELSE
            mtr=rdum-rtim(i)-otim(i)
         ENDIF
         rsum=rsum+mtr**2
         isum=isum+1
      ENDIF
   ENDDO
   DEALLOCATE(rtim,mtim,otim,ots)
   IF(rmtr.EQ.1)DEALLOCATE(srcm,ntsrc)
   IF(invst.EQ.1)DEALLOCATE(stt)
ENDDO
nt=isum
CLOSE(10)
CLOSE(20)
CLOSE(30)
IF(invst.EQ.1)THEN
   IF(itern.GT.1)CLOSE(40)
ENDIF
CLOSE(50)
CLOSE(60)
var=rsum/REAL(nt-1)
rsum=SQRT(rsum/REAL(nt))
rsum=1000.0*rsum
WRITE(6,'(f9.2,f9.5)')rsum,var

do i=1,nra
   call FmttPair_dealloc(pairs(i))
enddo
DEALLOCATE(pairs)
END PROGRAM residual
