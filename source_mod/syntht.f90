!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program constructs a synthetic set of traveltime
! residuals, with the option of adding noise.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM syntht
IMPLICIT NONE
INTEGER :: i,ii
INTEGER :: agn,rmfr,nr,ns,nt,checkstat,rseed,dray
INTEGER :: tswitch,nra,iscn,idum
INTEGER, DIMENSION (:), ALLOCATABLE :: ntsrc,ots
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: sdgn,dsd,rsum
REAL(KIND=i10), DIMENSION(:), ALLOCATABLE :: tres,srcm
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: mtim,rtim,otim
!REAL(KIND=i10), EXTERNAL :: gasdev
REAL, EXTERNAL :: gasdev
real(i10) tmp
CHARACTER (LEN=27) :: observed,rtimes,otimes,sources,receivers,obsray
!
! observed = input file of synthetic observed times
! rtimes = input file of reference times
! otimes = output file of traveltime residuals
! sources = file containing source information
! receivers = file containing receiver information
! agn = add Gaussian noise? (0=no,1=yes)
! rmfr = remove mean from residuals? (0=no,1=yes)
! sdgn = standard deviation of Gaussian noise
! dsd = default standard deviation for data covariance
! ns= number of sources
! nr = number of receivers
! nt = number of traveltimes
! tres = set of traveltime residuals
! tswitch = If no pick, then set to zero (not for synthetics)
! rseed = Random seed for Gaussian noise generation
! gasdev = Gaussian noise value
! dray = Delete rays? (0=no,1=yes)
! obsray = File containing rays for deletion
! nra = Number of receiver arrays
! srcm = model traveltime residual mean for each source
! ntsrc = Number of traveltimes for each source
! iscn = Source number
! nra = Number of receiver arrays
! rtim = Reference traveltime array
! otim = Observed traveltime array
! ots = Observed traveltime status (0 or 1)
!
! Read in the input parameters
!
tswitch=1
OPEN(UNIT=10,FILE='syntht.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a27)')observed
READ(10,'(a27)')rtimes
READ(10,'(a27)')sources
READ(10,'(a27)')receivers
READ(10,'(a27)')otimes
READ(10,*)agn
READ(10,*)sdgn
READ(10,*)rseed
READ(10,*)dsd
READ(10,*)rmfr
READ(10,*)dray
READ(10,'(a27)')obsray
CLOSE(10)
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
OPEN(UNIT=10,FILE=observed,STATUS='old')
OPEN(UNIT=20,FILE=rtimes,STATUS='old')
IF(dray.EQ.1)OPEN(UNIT=30,FILE=obsray,STATUS='old')
OPEN(UNIT=40,FILE=otimes,STATUS='unknown')
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
   ALLOCATE(tres(nt), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM syntht: REAL tres'
   ENDIF
   ALLOCATE(otim(nt))
   ALLOCATE(rtim(nt))
   ALLOCATE(ots(nt))
   ots=1
   DO i=1,nt
      READ(10,*)otim(i)
      READ(20,*)rtim(i)
      IF(dray.EQ.1)READ(30,*)ots(i)
   ENDDO
!
!  Compute average residuals if required
!
   IF(rmfr.EQ.1)THEN
      ALLOCATE(srcm(ns),ntsrc(ns))
      srcm=0.0
      ntsrc=0
      DO i=1,nt
         IF(ots(i).EQ.1)THEN
            iscn=INT((i-1)/nr)+1
            srcm(iscn)=srcm(iscn)+otim(i)-rtim(i)
            ntsrc(iscn)=ntsrc(iscn)+1
         ENDIF
      ENDDO
      DO i=1,ns
         srcm(i)=srcm(i)/REAL(ntsrc(i))
      ENDDO
   ENDIF
!
!  Compute residual and remove mean if required.
!
   DO i=1,nt
      tres(i)=otim(i)-rtim(i)
      IF(rmfr.EQ.1)THEN
         iscn=INT((i-1)/nr)+1
         tres(i)=tres(i)-srcm(iscn)
      ENDIF
   ENDDO
!
!  Add Gaussian  noise if required
!
   IF(agn.EQ.1)THEN
      dsd=sdgn
      DO i=1,nt
         tmp = gasdev(rseed)
         !tres(i)=tres(i)+gasdev(rseed)*sdgn
         tres(i)=tres(i)+tmp*sdgn
      ENDDO
   ENDIF
!
!  Write output to file
!
   DO i=1,nt
      WRITE(40,1)ots(i),tres(i),dsd
   ENDDO
   1 FORMAT(i4,f11.5,f8.4)
   DEALLOCATE(tres, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: PROGRAM syntht: REAL tres'
   ENDIF
   DEALLOCATE(otim,rtim,ots)
   IF(rmfr.EQ.1)DEALLOCATE(srcm,ntsrc) 
ENDDO
CLOSE(10)
CLOSE(20)
IF(dray.EQ.1)CLOSE(30)
CLOSE(40)
CLOSE(50)
CLOSE(60)
END PROGRAM syntht

REAL FUNCTION gasdev(idum)
IMPLICIT NONE
INTEGER :: iset,i
INTEGER, PARAMETER :: imax=100000
INTEGER :: idum
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: fac,rsq,v1,v2
REAL(KIND=i10), SAVE :: gset
!REAL(KIND=i10), EXTERNAL :: ran1
REAL, EXTERNAL :: ran1
real(i10) tmp
iset=0
IF(iset.EQ.0)then
   DO i=1,imax
      tmp = ran1(idum)
      !v1=2.0*ran1(idum)-1.
      !v2=2.0*ran1(idum)-1.
      v1 = 2.0 * tmp -1
      tmp = ran1(idum)
      v2 = 2.0 * tmp -1
      rsq=v1**2+v2**2
      if(rsq.LT.1.AND.rsq.NE.0.)EXIT
   ENDDO
   fac=sqrt(-2.0*LOG(rsq)/rsq)
   gset=v1*fac
   gasdev=v2*fac
   iset=1
ELSE
   gasdev=gset
   iset=0
ENDIF
END FUNCTION gasdev

REAL FUNCTION ran1(idum)
IMPLICIT NONE
INTEGER :: idum
INTEGER, PARAMETER :: ia=16807,im=2147483647,iq=127773
INTEGER, PARAMETER :: ir=2836,ntab=32,ndiv=1+(im-1)/ntab
INTEGER :: j,k
INTEGER, SAVE :: iy
INTEGER, DIMENSION (:), SAVE ::iv(ntab)
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10), PARAMETER :: eps=1.2e-7,rnmx=1.0-eps,am=1./im
iv=ntab*0
iy=0
IF(idum.LE.0.OR.iy.EQ.0)THEN
   DO j=ntab+8,1,-1
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      IF(idum.LT.0)idum=idum+im
      IF(j.LE.ntab)iv(j)=idum
   ENDDO
   iy=iv(1)
ENDIF
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
IF(idum.LT.0)idum=idum+im
j=1+iy/ndiv
iy=iv(j)
iv(j)=idum
ran1=min(am*iy,rnmx)
END FUNCTION ran1
