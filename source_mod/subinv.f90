!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! USE in any subroutine or function or other module.
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
INTEGER :: checkstat
INTEGER, SAVE :: subdim,ntr,nmp,asds,nvr,nvt,nvp
INTEGER, SAVE :: invst,nmpv,nmps
INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: cnfe,fcoln,tcnfe,tfcoln
REAL(KIND=i5), SAVE :: epsilon,eta,earth,stdf
REAL(KIND=i5), SAVE :: gor,got,gop,dnr,dnt,dnp
REAL(KIND=i5), DIMENSION (:), ALLOCATABLE, SAVE :: mo,mc,cm,dm
REAL(KIND=i5), DIMENSION (:), ALLOCATABLE, SAVE :: dobs,dmod,cd
REAL(KIND=i5), DIMENSION (:), ALLOCATABLE, SAVE :: frech,tfrech
REAL(KIND=i5), PARAMETER :: pi=3.1415926535898
!
! epsilon = The damping factor used in the inversion
! eta = The smoothing factor used in the inversion
! subdim = The subspace dimension used in the inversion
! ntr = Number of traveltimes
! nmp = Number of model parameters
! nmpv = Number of model parameters for velocity
! nmps = Number of model parameters for station terms
! asds = Apply second derivative smoothing (0=no,1=yes)
! mo = Reference model parameters
! mc = Current model parameters
! cm = Diagonal elements of a priori covariance matrix
! dobs = Observed traveltime
! dmod = Model traveltimes
! cd = Diagonal elements of data covariance matrix
! frech = Frechet derivatives
! cnfe = Cumulative number of nonzero elements of G for row i
! fcoln = Pointer to column number of G
! tfrech,tcnfe,tfcoln = Same as frech,cnfe,fcoln but for transpose
! dm = Model perturbation
! gor,got,gop = Grid origin in r,theta and phi
! dnr,dnt,dnp = Grid spacing in r, theta and phi
! nvr,nvt,nvp = Number of vertices in r, theta, phi
! earth = Earth radius in km.
! invst = Invert for station terms (0=no, 1=yes)
! stdf = Station term damping factor
!
END MODULE globalp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program will perform a single iteration of an 
! n-dimensional subspace inversion method using teleseismic
! data and forward computations from the 3-D fast marching
! program fmm3d. LU decomposition is used to perform the matrix
! inversion and SVD is used to construct an orthonormal basis
! for the projection matrix A. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM subinv
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,maxf,istep,jstep,jup,nrow,cstp,nra,idum
INTEGER :: invstep,nsrc,nrc,tstat,rmtr,stnm,iscn,isum
INTEGER, DIMENSION (:), ALLOCATABLE :: stpv,ntsrc,idst,idrt
INTEGER, PARAMETER :: maxss=50
integer,ALLOCATABLE     ::pick(:)
REAL(KIND=i5) :: dref,fracg
REAL(KIND=i5), DIMENSION (:), ALLOCATABLE :: srcm
CHARACTER (LEN=25) :: gridi,gridc,otimes,mtimes,rtimes
CHARACTER (LEN=25) :: frechet,sources,receive
CHARACTER (LEN=25) :: subiter,stfile
!
! gridi = File containing starting velocity field
! gridc = File containing current velocity field
! otimes = File containing observed traveltime residuals
! rtimes = File containing reference traveltimes
! mtimes = File containing model traveltimes
! frechet = File containing Frechet derivatives
! sources = File containing source information
! receive = File containing receiver information
! subiter = File indicating current inversion step
! invstep = Inversion step
! nsrc = Number of sources
! nrc = Number of receivers
! maxf = maximum size of frechet matrix
! fracg = Fraction of maximum size of Frechet matrix
! dref = Reference traveltimes
! istep,jstep = counting parameters
! jup = jstep+cnfe(istep)-1
! nrow = number of nonzero Frechet derivatives in row i
! stpv,cstp = variables for determining transpose of Frechet matrix
! tstat = status of observed time pick (0=none, 1=picked)
! maxss = Maximum recommended dimension of subspace
! rmtr = remove model traveltime residual (0=no,1=yes)
! stnm = station number for application of station term
! stfile = station term file
! srcm = model traveltime residual mean for each source
! ntsrc = Number of traveltimes for each source
! iscn = Source number
! nra = Number of receiver arrays
! idst = ID of source for each traveltime value
! idrt = ID of receiver for each traveltime value
!
! Read in the input parameters
!
OPEN(UNIT=10,FILE="subinv.in",STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)gridi
READ(10,1)gridc
READ(10,1)otimes
READ(10,1)rtimes
READ(10,1)mtimes
READ(10,*)rmtr
READ(10,1)frechet
READ(10,1)sources
READ(10,1)receive
READ(10,1)subiter
READ(10,*)epsilon
READ(10,*)subdim
READ(10,*)asds
READ(10,*)eta
READ(10,*)invst
READ(10,1)stfile
READ(10,*)stdf
READ(10,*)fracg
READ(10,*)earth
CLOSE(10)
1 FORMAT(a25)
IF(subdim.GT.maxss)THEN
   WRITE(6,*)'Subspace dimension is larger than ',maxss
   WRITE(6,*)'If you make the subspace dimension very'
   WRITE(6,*)'large, then the small improvements in convergence'
   WRITE(6,*)'will not justify the additional compute time.'
   WRITE(6,*)'Recommend using a subspace dimension'
   WRITE(6,*)'no larger than ',maxss
ELSE IF(subdim.LT.1)THEN
   WRITE(6,*)'Subspace dimension cannot be less than one!!!'
   WRITE(6,*)'Terminating program!!!'
   STOP
ENDIF
IF(subdim.EQ.1.AND.invst.EQ.1)THEN
   WRITE(6,*)'ERROR!!!!!!'
   WRITE(6,*)'Cannot have a subspace dimension of 1'
   WRITE(6,*)'and invert for velocity AND station terms.'
   WRITE(6,*)'Set subspace dimension to value greater than 1.'
   WRITE(6,*)'TERMINATING PROGRAM!!!!!'
ENDIF
!
! Determine the iteration number
!
OPEN(UNIT=10,FILE=subiter,STATUS='old')
READ(10,*)invstep
CLOSE(10)
!
! The maximum number of traveltime residuals is equal to the
! sum of the product of the number of receivers and sources for each 
! receiver array.
!
OPEN(UNIT=10,FILE=sources,STATUS='old')
OPEN(UNIT=20,FILE=receive,STATUS='old')
READ(10,*)nra
READ(20,*)idum
IF(idum.NE.nra)THEN
   WRITE(6,*)'ERROR!!!'
   WRITE(6,*)'Source and receiver files are'
   WRITE(6,*)'inconsistent!!!!'
   WRITE(6,*)'First line of each file should'
   WRITE(6,*)'be identical!!!'
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
ntr=0
nmps=0
nsrc=0
DO i=1,nra
   READ(10,*)idum
   READ(20,*)nrc
   ntr=ntr+nrc*idum
   nmps=nmps+nrc
   nsrc=nsrc+idum
   DO j=1,idum
      READ(10,*)
      READ(10,*)
   ENDDO
   DO j=1,nrc
      READ(20,*)
   ENDDO
ENDDO
CLOSE(10)
CLOSE(20)

! read valid_rays
ALLOCATE(pick(ntr))
istep = 1
OPEN(UNIT=10,FILE=sources,STATUS='old')
OPEN(UNIT=20,FILE=receive,STATUS='old')
open(99,file='valid_rays.dat')
READ(10,*)nra
READ(20,*)idum
do i=1,nra
   READ(10,*)idum
   READ(20,*)nrc
   do j=1,idum
      read(99,*)isum
      pick(istep:istep + nrc-1) = isum
      istep = istep + nrc
      READ(10,*)
      READ(10,*)
   enddo
   DO j=1,nrc
      READ(20,*)
   enddo
enddo
close(10)
close(20)
close(99)


!
! If station terms are included, or the model
! traveltime residual is to be removed, then
! we need to associate each traveltime with a
! source
!
IF(rmtr.EQ.1.OR.invst.EQ.1)THEN
   ALLOCATE(idst(ntr))
   ALLOCATE(idrt(ntr))
   OPEN(UNIT=10,FILE=sources,STATUS='old')
   OPEN(UNIT=20,FILE=receive,STATUS='old')
   READ(10,*)
   READ(20,*)
   isum=1
   nsrc=0
   nmps=0
   DO i=1,nra
      READ(10,*)idum
      READ(20,*)nrc
      DO j=1,idum
         READ(10,*)
         READ(10,*)
      ENDDO
      DO j=1,nrc
         READ(20,*)
      ENDDO
      DO j=1,idum
         DO k=1,nrc
            idst(isum)=j+nsrc
            idrt(isum)=k+nmps
            isum=isum+1
         ENDDO
      ENDDO
      nsrc=nsrc+idum
      nmps=nmps+nrc
   ENDDO
   CLOSE(10)
   CLOSE(20) 
ENDIF
!
! Now read in reference grid and current grid (if required)
!
OPEN(UNIT=10,FILE=gridi,status='old')
READ(10,*)nvr,nvt,nvp
nmp=(nvr+2)*(nvt+2)*(nvp+2)
READ(10,*)gor,got,gop
READ(10,*)dnr,dnt,dnp
nmpv=nmp
IF(invst.EQ.1)nmp=nmpv+nmps
!
! Allocate memory to model arrays
!
ALLOCATE(mo(nmp), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL mo'
ENDIF
ALLOCATE(mc(nmp), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL mc'
ENDIF
ALLOCATE(cm(nmp), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL cm'
ENDIF
DO i=1,nmpv
   READ(10,*)mo(i),cm(i)
   cm(i)=cm(i)**2
   IF(invstep.LE.1)mc(i)=mo(i)
ENDDO
CLOSE(10)
IF(invstep.GT.1)THEN
   OPEN(UNIT=10,FILE=gridc,status='old')
   READ(10,*)
   READ(10,*)
   READ(10,*)
   DO i=1,nmpv
      READ(10,*)mc(i)
   ENDDO
   CLOSE(10)
ENDIF
!
! If required, include station terms.
!
IF(invst.EQ.1)THEN
   DO i=nmpv+1,nmp
      mo(i)=0.0
      cm(i)=1.0
      IF(invstep.LE.1)mc(i)=mo(i)
   ENDDO
   IF(invstep.GT.1)THEN
      OPEN(UNIT=10,FILE=stfile,STATUS='old')
      DO i=nmpv+1,nmp
         READ(10,*)mc(i)
      ENDDO
      CLOSE(10)
   ENDIF
ENDIF
!
! Now read in the observed and model traveltime residuals, plus
! the Frechet matrix.
!
maxf=int(fracg*ntr*nmp)
ALLOCATE(dobs(ntr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL dobs'
ENDIF
ALLOCATE(cd(ntr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL cd'
ENDIF
ALLOCATE(dmod(ntr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL dmod'
ENDIF
ALLOCATE(frech(maxf), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL frech'
ENDIF
ALLOCATE(fcoln(maxf), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL fcoln'
ENDIF
ALLOCATE(cnfe(0:ntr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL cnfe'
ENDIF
OPEN(UNIT=10,FILE=otimes,STATUS='old')
OPEN(UNIT=20,FILE=mtimes,STATUS='old')
OPEN(UNIT=30,FILE=rtimes,STATUS='old')
OPEN(UNIT=40,FILE=frechet,FORM='unformatted',STATUS='old')
istep=1
jstep=0
cnfe(0)=0
IF(rmtr.EQ.1)THEN
   ALLOCATE(srcm(nsrc),ntsrc(nsrc))
   srcm=0.0
   ntsrc=0
ENDIF
DO i=1,ntr
   READ(10,*)tstat,dobs(istep),cd(istep)
   tstat = tstat * pick(i)
   READ(20,*)dmod(istep)
   READ(30,*)dref
   READ(40)nrow
!
!  If station terms are included, identify receiver
!  and apply correction to dmod
!
   IF(invst.EQ.1)THEN
      stnm=idrt(i)
      dmod(istep)=dmod(istep)+mc(nmpv+stnm)
      nrow=nrow+1
   ENDIF
   cd(istep)=cd(istep)**2
   cnfe(istep)=jstep+nrow
   IF(nrow.GT.0)THEN
      jstep=jstep+1
      jup=jstep+nrow-1
      IF(jup.GT.maxf)THEN
         WRITE(6,*)'Frechet matrix too large for last input'
         WRITE(6,*)'parameter in file subinv.in'
         WRITE(6,*)'STOPPING PROGRAM!!!'
         STOP
      ENDIF
      IF(invst.EQ.1)THEN
         IF(nrow.GT.1)THEN
            DO j=jstep,jup-1
               READ(40)fcoln(j),frech(j)
            ENDDO
         ENDIF
         fcoln(jup)=nmpv+stnm
         frech(jup)=1.0
      ELSE
         DO j=jstep,jup
            READ(40)fcoln(j),frech(j)
         ENDDO
      ENDIF
      IF(tstat.EQ.1)THEN
         jstep=jup
      ELSE
         jstep=jstep-1
      ENDIF
   ENDIF
   IF(tstat.EQ.1)THEN
      dmod(istep)=dmod(istep)-dref
      IF(rmtr.EQ.1)THEN
         iscn=idst(i)
         srcm(iscn)=srcm(iscn)+dmod(istep)
         ntsrc(iscn)=ntsrc(iscn)+1
      ENDIF
      istep=istep+1
   ENDIF
ENDDO
ntr=istep-1
IF(rmtr.EQ.1)THEN
!
!  Remove mean from model residuals if required.
!
   istep=1
   DO i=1,nsrc
      srcm(i)=srcm(i)/REAL(ntsrc(i))
      DO j=1,ntsrc(i)
         dmod(istep)=dmod(istep)-srcm(i)
         istep=istep+1
      ENDDO
   ENDDO
   DEALLOCATE(srcm,ntsrc)
ENDIF      
maxf=jstep
CLOSE(10)
CLOSE(20)
CLOSE(30)
CLOSE(40)
!
! Now construct the transpose of the Frechet matrix
!
ALLOCATE(tfrech(maxf), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL tfrech'
ENDIF
ALLOCATE(tfcoln(maxf), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL tfcoln'
ENDIF
ALLOCATE(tcnfe(0:nmp), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL tcnfe'
ENDIF
ALLOCATE(stpv(ntr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL tcnfe'
ENDIF
stpv=1
jstep=0
tcnfe(0)=0
DO i=1,nmp
   DO j=1,ntr
      IF(stpv(j).LE.cnfe(j)-cnfe(j-1))THEN
         cstp=fcoln(cnfe(j-1)+stpv(j))
         IF(cstp.EQ.i)THEN
            jstep=jstep+1
            tfrech(jstep)=frech(cnfe(j-1)+stpv(j))
            tfcoln(jstep)=j
            stpv(j)=stpv(j)+1
         ENDIF
      ENDIF
   ENDDO
   tcnfe(i)=jstep
ENDDO
!
! Allocate memory to model perturbation
!
ALLOCATE(dm(nmp), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subinv: REAL dm'
ENDIF
!
! Now we have set up all the vectors and matrices we require for
! the subspace inversion scheme. Call a subroutine for
! performing the inversion
!
nvr=nvr+2
nvt=nvt+2
nvp=nvp+2
CALL subspace
nvr=nvr-2
nvt=nvt-2
nvp=nvp-2
!
! Write new model to file
!
OPEN(UNIT=10,FILE=gridc,STATUS='unknown')
print*,'min and max perturbation',minval(dm),maxval(dm)
WRITE(10,*)nvr,nvt,nvp
WRITE(10,*)gor,got,gop
WRITE(10,*)dnr,dnt,dnp
DO i=1,nmpv
   WRITE(10,*)mc(i)+dm(i)
ENDDO
CLOSE(10)
!
! If required, write station terms to file.
!
IF(invst.EQ.1)THEN
   OPEN(UNIT=10,FILE=stfile,STATUS='unknown')
   DO i=nmpv+1,nmp
      WRITE(10,*)mc(i)+dm(i)
   ENDDO
   CLOSE(10)
ENDIF
DEALLOCATE(mo,mc,cm,dm, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM subinv: mo,mc,cm,dm'
ENDIF
DEALLOCATE(dobs,cd,dmod, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM subinv: dobs,cd,dmod'
ENDIF
DEALLOCATE(cnfe,fcoln,frech, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM subinv: cnfe,fcoln,frech'
ENDIF
DEALLOCATE(tcnfe,tfcoln,tfrech,stpv, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM subinv: tcnfe,tfcoln,tfrech,stpv'
ENDIF
IF(rmtr.EQ.1.OR.invst.EQ.1)THEN
   DEALLOCATE(idst,idrt)
ENDIF
DEALLOCATE(pick)
END PROGRAM subinv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates the model perturbation using
! the subspace inversion method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE subspace
USE globalp
IMPLICIT NONE
INTEGER :: ii,jj,i,j,k,l,is1,is2,is3,npc,nits
INTEGER, DIMENSION (2) :: inds, inde
REAL(KIND=i5) :: rsum1,rdnt,rdnp,rgot,rgop,ddem,thet,rad,epst
REAL(KIND=i5), DIMENSION (:,:), ALLOCATABLE :: a,mati
REAL(KIND=i5), DIMENSION (:), ALLOCATABLE :: gamma,dtrav,ga,dem
REAL(KIND=i5), DIMENSION (subdim) :: r1mat,r2mat
!
! a = Projection matrix
! gamma = gradient vector
! dtrav = weighted differential traveltime residuals
! ga = Ga, the product of Frechet matrix and vector of a
! dem = Finite difference estimate of spatial derivative
! ddem = dem premultipled by transform of D
! is1,is2,is3 = counters for smoothness operators
! mati = matrix for inversion
! r1mat,r2mat = RHS vectors of matrix equation for model perturbation
! rdnt,rndp,rgot,rgop = dnt,dnp,got,gop in radians
! thet,rad = value of theta and radius at node point
! epst = substitution variable for damping factor epsilon
! npc = Number of parameter classes
! inds(i), inde(i) = Index start and end for each parameter class.
! nits = Number of iterations for subspace vector generation
!
  rdnt=dnt*pi/180.0
  rdnp=dnp*pi/180.0
  rgot=got*pi/180.0
  rgop=gop*pi/180.0
!
! The first step is to compute the projection matrix a. Allocate
! memory to this array.
!
ALLOCATE(a(nmp,subdim), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subspace: REAL a'
ENDIF
a=0.0
!
! The first subspace dimension is given by the gradient vector in
! model space. Allocate memory to this vector.
!
ALLOCATE(gamma(nmp), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subspace: REAL gamma'
ENDIF
!
! Allocate memory to weighted differential traveltime residuals
!
ALLOCATE(dtrav(ntr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subspace: REAL dtrav'
ENDIF
DO i=1,ntr
   dtrav(i)=(dmod(i)-dobs(i))/cd(i)
ENDDO
!
! Now compute gamma and the first subspace vector. There are
! two options here, depending on whether smoothing is
! applied or not.
!
epst=epsilon
DO i=1,nmp
   rsum1=0.0
   IF(tcnfe(i).GT.tcnfe(i-1))THEN
      DO j=tcnfe(i-1)+1,tcnfe(i)
         rsum1=rsum1+tfrech(j)*dtrav(tfcoln(j))
      ENDDO
   ENDIF
   IF(i.GT.nmpv)epst=stdf
   gamma(i)=rsum1+epst*(mc(i)-mo(i))/cm(i)
ENDDO
!
! Apply second derivative smoothing if required
!
IF(asds.EQ.1)THEN
   ALLOCATE(dem(nmpv), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM subspace: REAL dem'
   ENDIF
   dem=0.0
   DO i=1,nvp
      DO j=1,nvt
         thet=sin(rgot+(j-1)*rdnt)
         DO k=1,nvr
            rad=gor-(k-1)*dnr+earth
            is2=(i-1)*nvt*nvr+(j-1)*nvr+k 
            IF(k.NE.1.AND.k.NE.nvr)THEN
               is1=(i-1)*nvt*nvr+(j-1)*nvr+k-1
               is3=(i-1)*nvt*nvr+(j-1)*nvr+k+1
               rsum1=mc(is1)-2.0*mc(is2)+mc(is3)
               rsum1=rsum1-mo(is1)+2.0*mo(is2)-mo(is3)
               dem(is2)=dem(is2)+rsum1
            ENDIF
            IF(j.NE.1.AND.j.NE.nvt)THEN
               is1=(i-1)*nvt*nvr+(j-2)*nvr+k
               is3=(i-1)*nvt*nvr+(j)*nvr+k
               rsum1=mc(is1)-2.0*mc(is2)+mc(is3)
               rsum1=rsum1-mo(is1)+2.0*mo(is2)-mo(is3)
               dem(is2)=dem(is2)+rsum1*dnr**2/(rad*rdnt)**2
            ENDIF
            IF(i.NE.1.AND.i.NE.nvp)THEN
               is1=(i-2)*nvt*nvr+(j-1)*nvr+k
               is3=(i)*nvt*nvr+(j-1)*nvr+k
               rsum1=mc(is1)-2.0*mc(is2)+mc(is3)
               rsum1=rsum1-mo(is1)+2.0*mo(is2)-mo(is3)
               dem(is2)=dem(is2)+rsum1*dnr**2/(rad*thet*rdnp)**2
            ENDIF
         ENDDO
      ENDDO
   ENDDO
!
!  Now repeat the loop to premultiply by the transform of D
!
   DO i=1,nvp
      DO j=1,nvt
         thet=sin(rgot+(j-1)*rdnt)
         DO k=1,nvr
            ddem=0.0
            rad=gor-(k-1)*dnr+earth
            is2=(i-1)*nvt*nvr+(j-1)*nvr+k
            IF(k.NE.1.AND.k.NE.nvr)THEN
               ddem=ddem-2.0*dem(is2)
            ENDIF
            IF(k-2.GE.1)THEN
               is1=(i-1)*nvt*nvr+(j-1)*nvr+k-1
               ddem=ddem+dem(is1)
            ENDIF
            IF(k+2.LE.nvr)THEN
               is3=(i-1)*nvt*nvr+(j-1)*nvr+k+1
               ddem=ddem+dem(is3)
            ENDIF
            IF(j.NE.1.AND.j.NE.nvt)THEN
               ddem=ddem-2.0*dem(is2)*dnr**2/(rad*rdnt)**2
            ENDIF
            IF(j-2.GE.1)THEN
               is1=(i-1)*nvt*nvr+(j-2)*nvr+k
               ddem=ddem+dem(is1)*dnr**2/(rad*rdnt)**2
            ENDIF
            IF(j+2.LE.nvt)THEN
               is3=(i-1)*nvt*nvr+(j)*nvr+k
               ddem=ddem+dem(is3)*dnr**2/(rad*rdnt)**2
            ENDIF
            IF(i.NE.1.AND.i.NE.nvp)THEN
               ddem=ddem-2.0*dem(is2)*dnr**2/(rad*thet*rdnp)**2
            ENDIF
            IF(i-2.GE.1)THEN
               is1=(i-2)*nvt*nvr+(j-1)*nvr+k
               ddem=ddem+dem(is1)*dnr**2/(rad*thet*rdnp)**2
            ENDIF
            IF(i+2.LE.nvt)THEN
               is3=(i)*nvt*nvr+(j-1)*nvr+k
               ddem=ddem+dem(is3)*dnr**2/(rad*thet*rdnp)**2
            ENDIF
            gamma(is2)=gamma(is2)+eta*ddem
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
! Determine the number of parameter classes and the index
! of the start and end of each class
!
npc=1
inds(1)=1
inde(1)=nmpv
IF(invst.EQ.1)THEN
   npc=2
   inds(2)=nmpv+1
   inde(2)=nmp
ENDIF
!
! Now compute the first npc basis vectors
!
DO i=1,npc
   DO j=inds(i),inde(i)
      a(j,i)=gamma(j)*cm(j)
   ENDDO
ENDDO
!
! Normalize the subspace vector
!
DO i=1,npc
   rsum1=0.0
   DO j=inds(i),inde(i)
      rsum1=rsum1+a(j,i)**2
   ENDDO
   DO j=inds(i),inde(i)
      a(j,i)=a(j,i)/SQRT(rsum1)
   ENDDO
ENDDO
!
! Now compute the remaining subspace vectors if
! required
!
ALLOCATE(ga(ntr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subspace: REAL nmp'
ENDIF
IF(subdim.GT.npc)THEN
   nits=npc
   is1=npc
   is2=npc
   DO ii=1,subdim
      DO i=1,nits
!
!        First, calculate Ga weighted by Cd
!
         DO j=1,ntr
            rsum1=0.0
            IF(cnfe(j).GT.cnfe(j-1))THEN
               DO k=cnfe(j-1)+1,cnfe(j)
                  rsum1=rsum1+frech(k)*a(fcoln(k),is2-nits+i)
               ENDDO
            ENDIF
            ga(j)=rsum1/cd(j)
         ENDDO
!
!        Now calculate the new subspace vector
!
         DO jj=1,npc
            is1=is1+1
            IF(is1.GT.subdim)EXIT
            DO j=inds(jj),inde(jj)
               rsum1=0.0
               IF(tcnfe(j).GT.tcnfe(j-1))THEN
                  DO k=tcnfe(j-1)+1,tcnfe(j)
                     rsum1=rsum1+tfrech(k)*ga(tfcoln(k))
                  ENDDO
               ENDIF
               IF(j.LE.nmpv)THEN
                  a(j,is1)=rsum1*cm(j)+epsilon*a(j,is2-nits+i)
               ELSE
                  a(j,is1)=rsum1*cm(j)+stdf*a(j,is2-nits+i)
               ENDIF
            ENDDO
            !
            ! Normalize the subspace vector
            !
            rsum1=0.0
            DO j=inds(jj),inde(jj)
               rsum1=rsum1+a(j,is1)**2
            ENDDO
            DO j=inds(jj),inde(jj)
               a(j,is1)=a(j,is1)/SQRT(rsum1)
            ENDDO
         ENDDO
         IF(is1.GT.subdim)EXIT
      ENDDO
      IF(is1.GT.subdim)EXIT
      nits=nits*npc
      is2=is2+nits
   ENDDO
ENDIF
!
! If subdim is greater than 1, then we need to orthogonalize
! the vectors that make up the projection matrix A to
! avoid interdependence. This can be done using SVD.
!
CALL svdcmp(a)
!
! Now we have our projection matrix A. Construct the matrix that
! is to be inverted.
!
ALLOCATE(mati(subdim,subdim), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM subspace: REAL mati'
ENDIF
mati=0.0
IF(asds.NE.1)THEN
   DO i=1,subdim
      DO j=1,ntr
         rsum1=0.0
         IF(cnfe(j).GT.cnfe(j-1))THEN
            DO k=cnfe(j-1)+1,cnfe(j)
               rsum1=rsum1+frech(k)*a(fcoln(k),i)
            ENDDO
         ENDIF
         ga(j)=rsum1/cd(j)
      ENDDO
      epst=epsilon
      DO j=1,nmp
         rsum1=0.0
         IF(tcnfe(j).GT.tcnfe(j-1))THEN
            DO k=tcnfe(j-1)+1,tcnfe(j)
               rsum1=rsum1+tfrech(k)*ga(tfcoln(k))
            ENDDO
         ENDIF
         IF(j.GT.nmpv)epst=stdf
         rsum1=rsum1+epst*a(j,i)/cm(j)
         DO k=1,subdim
            mati(k,i)=mati(k,i)+rsum1*a(j,k)
         ENDDO
      ENDDO
   ENDDO
ELSE 
   DO ii=1,subdim
      DO j=1,ntr
         rsum1=0.0
         IF(cnfe(j).GT.cnfe(j-1))THEN
            DO k=cnfe(j-1)+1,cnfe(j)
               rsum1=rsum1+frech(k)*a(fcoln(k),ii)
            ENDDO
         ENDIF
         ga(j)=rsum1/cd(j)
      ENDDO
!
!     Now premultiply by smoothing operator
!
      dem=0.0
      DO i=1,nvp
         DO j=1,nvt
            thet=sin(rgot+(j-1)*rdnt)
            DO k=1,nvr
               rad=gor-(k-1)*dnr+earth
               is2=(i-1)*nvt*nvr+(j-1)*nvr+k
               IF(k.NE.1.AND.k.NE.nvr)THEN
                  is1=(i-1)*nvt*nvr+(j-1)*nvr+k-1
                  is3=(i-1)*nvt*nvr+(j-1)*nvr+k+1
                  rsum1=a(is1,ii)-2.0*a(is2,ii)+a(is3,ii)
                  dem(is2)=dem(is2)+rsum1
               ENDIF
               IF(j.NE.1.AND.j.NE.nvt)THEN
                  is1=(i-1)*nvt*nvr+(j-2)*nvr+k
                  is3=(i-1)*nvt*nvr+(j)*nvr+k
                  rsum1=a(is1,ii)-2.0*a(is2,ii)+a(is3,ii)
                  dem(is2)=dem(is2)+rsum1*dnr**2/(rad*rdnt)**2
               ENDIF
               IF(i.NE.1.AND.i.NE.nvp)THEN
                  is1=(i-2)*nvt*nvr+(j-1)*nvr+k
                  is3=(i)*nvt*nvr+(j-1)*nvr+k
                  rsum1=a(is1,ii)-2.0*a(is2,ii)+a(is3,ii)
                  dem(is2)=dem(is2)+rsum1*dnr**2/(rad*thet*rdnp)**2
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      DO i=1,nvp
         DO j=1,nvt
            thet=sin(rgot+(j-1)*rdnt)
            DO k=1,nvr
               rad=gor-(k-1)*dnr+earth
               is2=(i-1)*nvt*nvr+(j-1)*nvr+k
               rsum1=0.0
               IF(tcnfe(is2).GT.tcnfe(is2-1))THEN
                  DO l=tcnfe(is2-1)+1,tcnfe(is2)
                     rsum1=rsum1+tfrech(l)*ga(tfcoln(l))
                  ENDDO
               ENDIF
!
!              Now calculate element due to smoothing
!
               ddem=0.0
               is2=(i-1)*nvt*nvr+(j-1)*nvr+k
               IF(k.NE.1.AND.k.NE.nvr)THEN
                  ddem=ddem-2.0*dem(is2)
               ENDIF
               IF(k-2.GE.1)THEN
                  is1=(i-1)*nvt*nvr+(j-1)*nvr+k-1
                  ddem=ddem+dem(is1)
               ENDIF
               IF(k+2.LE.nvr)THEN
                  is3=(i-1)*nvt*nvr+(j-1)*nvr+k+1
                  ddem=ddem+dem(is3)
               ENDIF
               IF(j.NE.1.AND.j.NE.nvt)THEN
                  ddem=ddem-2.0*dem(is2)*dnr**2/(rad*rdnt)**2
               ENDIF
               IF(j-2.GE.1)THEN
                  is1=(i-1)*nvt*nvr+(j-2)*nvr+k
                  ddem=ddem+dem(is1)*dnr**2/(rad*rdnt)**2
               ENDIF
               IF(j+2.LE.nvt)THEN
                  is3=(i-1)*nvt*nvr+(j)*nvr+k
                  ddem=ddem+dem(is3)*dnr**2/(rad*rdnt)**2
               ENDIF
               IF(i.NE.1.AND.i.NE.nvp)THEN
                  ddem=ddem-2.0*dem(is2)*dnr**2/(rad*thet*rdnp)**2
               ENDIF
               IF(i-2.GE.1)THEN
                  is1=(i-2)*nvt*nvr+(j-1)*nvr+k
                  ddem=ddem+dem(is1)*dnr**2/(rad*thet*rdnp)**2
               ENDIF
               IF(i+2.LE.nvt)THEN
                  is3=(i)*nvt*nvr+(j-1)*nvr+k
                  ddem=ddem+dem(is3)*dnr**2/(rad*thet*rdnp)**2
               ENDIF
               rsum1=rsum1+epsilon*a(is2,ii)/cm(is2)+eta*ddem
               DO l=1,subdim
                  mati(l,ii)=mati(l,ii)+rsum1*a(is2,l)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     If station terms are included, then we need to add their
!     contribution to the matrix to be inverted
!
      IF(invst.EQ.1)THEN
         DO j=nmpv+1,nmp
            rsum1=0.0
            IF(tcnfe(j).GT.tcnfe(j-1))THEN
               DO k=tcnfe(j-1)+1,tcnfe(j)
                  rsum1=rsum1+tfrech(k)*ga(tfcoln(k))
               ENDDO
            ENDIF
            rsum1=rsum1+stdf*a(j,ii)/cm(j)
            DO k=1,subdim
               mati(k,ii)=mati(k,ii)+rsum1*a(j,k)
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDIF
!
! Now perform the matrix inversion using LU decomposition
!
CALL luinv(subdim,mati)
!
! Calculate the model perturbation
!
r1mat=0.0
DO i=1,subdim
   DO j=1,nmp
      r1mat(i)=r1mat(i)+a(j,i)*gamma(j)
   ENDDO
ENDDO
r2mat=0.0
DO i=1,subdim
   DO j=1,subdim
      r2mat(i)=r2mat(i)+mati(i,j)*r1mat(j)
   ENDDO
ENDDO
dm=0.0
DO i=1,nmp
   DO j=1,subdim
      dm(i)=dm(i)+a(i,j)*r2mat(j)
   ENDDO
   dm(i)=-dm(i)
ENDDO
DEALLOCATE(a,gamma,dtrav, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM subspace: REAL a,gamma,dtrav'
ENDIF
DEALLOCATE(ga,mati, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM subspace: REAL ga,mati'
ENDIF
IF(asds.EQ.1)THEN
   DEALLOCATE(dem, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: PROGRAM subspace: REAL dm'
   ENDIF
ENDIF
END SUBROUTINE subspace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine (adapted from Numerical Recipes) computes
! the singular value decomposition of the projection
! matrix a to produce an orthonormal basis. The dimension of
! the subspace is automatically reduced if the subspace does
! not span subdim.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE svdcmp(a)
USE globalp
IMPLICIT NONE
INTEGER :: m,n
INTEGER :: i,its,j,jj,k,l,nm,isw,isum
INTEGER, DIMENSION (subdim) :: wnz
INTEGER, PARAMETER :: nmax=500
REAL(KIND=i5), DIMENSION (nmp,subdim) :: a
REAL(KIND=i5), DIMENSION (subdim,subdim) :: v
REAL(KIND=i5), DIMENSION (subdim) :: w
REAL(KIND=i5), DIMENSION (nmax) :: rv1
REAL(KIND=i5) :: anorm,c,f,g,h,s,scale,x,y,z
REAL(KIND=i5), EXTERNAL :: pythag
REAL(KIND=i5), PARAMETER :: wtol=1.0e-7
m=nmp
n=subdim
g=0.0
scale=0.0
anorm=0.0
DO i=1,n
   l=i+1
   rv1(i)=scale*g
   g=0.0
   s=0.0
   scale=0.0
   If(i.LE.m)THEN
      DO k=i,m
         scale=scale+ABS(a(k,i))
      ENDDO
      If(scale.NE.0.0)THEN
         DO k=i,m
            a(k,i)=a(k,i)/scale
            s=s+a(k,i)*a(k,i)
         ENDDO
         f=a(i,i)
         g=-SIGN(SQRT(s),f)
         h=f*g-s
         a(i,i)=f-g
         DO j=l,n
            s=0.0
            DO k=i,m
               s=s+a(k,i)*a(k,j)
            ENDDO
            f=s/h
            DO k=i,m
               a(k,j)=a(k,j)+f*a(k,i)
            ENDDO
         ENDDO
         DO k=i,m
            a(k,i)=scale*a(k,i)
         ENDDO
      ENDIF
   ENDIF
   w(i)=scale*g
   g=0.0
   s=0.0
   scale=0.0
   IF((i.LE.m).AND.(i.NE.n))THEN
      DO k=l,n
         scale=scale+ABS(a(i,k))
      ENDDO
      IF(scale.NE.0.0)THEN
         DO k=l,n
            a(i,k)=a(i,k)/scale
            s=s+a(i,k)*a(i,k)
         ENDDO
         f=a(i,l)
         g=-SIGN(SQRT(s),f)
         h=f*g-s
         a(i,l)=f-g
         DO k=l,n
            rv1(k)=a(i,k)/h
         ENDDO
         DO j=l,m
            s=0.0
            DO k=l,n
               s=s+a(j,k)*a(i,k)
            ENDDO
            DO k=l,n
               a(j,k)=a(j,k)+s*rv1(k)
            ENDDO
         ENDDO
         DO k=l,n
            a(i,k)=scale*a(i,k)
         ENDDO
      ENDIF
   ENDIF
   anorm=MAX(anorm,(ABS(w(i))+abs(rv1(i))))
ENDDO
DO i=n,1,-1
   IF(i.LT.n)THEN
      IF(g.NE.0.0)THEN
         DO j=l,n
            v(j,i)=(a(i,j)/a(i,l))/g
         ENDDO
         DO j=l,n
            s=0.0
            DO k=l,n
               s=s+a(i,k)*v(k,j)
            ENDDO
            DO k=l,n
               v(k,j)=v(k,j)+s*v(k,i)
            ENDDO
         ENDDO
      ENDIF
      DO j=l,n
         v(i,j)=0.0
         v(j,i)=0.0
      ENDDO
   ENDIF
   v(i,i)=1.0
   g=rv1(i)
   l=i
ENDDO
DO i=MIN(m,n),1,-1
   l=i+1
   g=w(i)
   DO j=l,n
      a(i,j)=0.0
   ENDDO
   IF(g.NE.0.0)THEN
      g=1.0/g
      DO j=l,n
         s=0.0
         DO k=l,m
            s=s+a(k,i)*a(k,j)
         ENDDO
         f=(s/a(i,i))*g
         DO k=i,m
            a(k,j)=a(k,j)+f*a(k,i)
         ENDDO
      ENDDO
      DO j=i,m
         a(j,i)=a(j,i)*g
      ENDDO
   ELSE
      DO j=i,m
         a(j,i)=0.0
      ENDDO
   ENDIF
   a(i,i)=a(i,i)+1.0
ENDDO
DO k=n,1,-1
   DO its=1,30
      isw=0
      DO l=k,1,-1
         nm=l-1
         IF((ABS(rv1(l))+anorm).EQ.anorm)THEN
            isw=1
            EXIT
         ENDIF
         IF((ABS(w(nm))+anorm).EQ.anorm)EXIT
      ENDDO
      IF(isw.EQ.0)THEN
         c=0.0
         s=1.0
         DO i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            IF((ABS(f)+anorm).EQ.anorm)EXIT
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c=g*h
            s=-f*h
            DO j=1,m
               y=a(j,nm)
               z=a(j,i)
               a(j,nm)=y*c+z*s
               a(j,i)=-y*s+z*c
            ENDDO
         ENDDO
      ENDIF
      z=w(k)
      IF(l.EQ.k)THEN
         IF(z.LT.0.0)THEN
            w(k)=-z
            DO j=1,n
               v(j,k)=-v(j,k)
            ENDDO
         ENDIF
         EXIT
      ENDIF
      IF(its.EQ.30)THEN
         WRITE(6,*)'No convergence in svdcmp!'
         WRITE(6,*)'Check the integrity of your'
         WRITE(6,*)'input files!!!!!!!!'
         WRITE(6,*)'Terminating program subinv!'
         STOP
      ENDIF
      x=w(l)
      nm=k-1
      y=w(nm)
      g=rv1(nm)
      h=rv1(k) 
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
      g=pythag(f,1.0)
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x
      c=1.0
      s=1.0
      DO j=l,nm
         i=j+1
         g=rv1(i)
         y=w(i)
         h=s*g
         g=c*g
         z=pythag(f,h)
         rv1(j)=z
         c=f/z
         s=h/z
         f=x*c+g*s
         g=-x*s+g*c
         h=y*s
         y=y*c
         DO jj=1,n
            x=v(jj,j)
            z=v(jj,i)
            v(jj,j)=x*c+z*s
            v(jj,i)=-x*s+z*c
         ENDDO
         z=pythag(f,h)
         w(j)=z
         IF(z.NE.0.0)THEN
            z=1.0/z
            c=f*z
            s=h*z
         ENDIF
         f=c*g+s*y
         x=-s*g+c*y
         DO jj=1,m
            y=a(jj,j)
            z=a(jj,i)
            a(jj,j)=y*c+z*s
            a(jj,i)=-y*s+z*c
         ENDDO
      ENDDO
      rv1(l)=0.0
      rv1(k)=f
      w(k)=x
   ENDDO
ENDDO
!
! Finally, eliminate any vectors with very small w's.
!
isum=0
DO i=1,n
   IF(ABS(w(i)).GE.wtol)THEN
      isum=isum+1
      wnz(isum)=i
   ENDIF
ENDDO
IF(isum.LT.subdim)THEN
   WRITE(6,*)'Due to redundancy, the subspace dimension'
   WRITE(6,*)'is being reduced from ',subdim,' to ',isum
   WRITE(6,*)'Message from SVD orthogonalization algorithm'
   DO i=1,isum
      DO j=1,nmp
         a(j,i)=a(j,wnz(i))
      ENDDO
   ENDDO
   subdim=isum
ENDIF
END SUBROUTINE svdcmp

REAL FUNCTION pythag(a,b)
USE globalp
IMPLICIT NONE
REAL(KIND=i5), INTENT(IN) :: a,b
REAL(KIND=i5) :: absa,absb
absa=ABS(a)
absb=ABS(b)
IF(absa.GT.absb)THEN
   pythag=absa*SQRT(1.0+(absb/absa)**2)
ELSE
   IF(absb.EQ.0.0)THEN
      pythag=0.0
   ELSE
      pythag=absb*SQRT(1.0+(absa/absb)**2)
   ENDIF
ENDIF
END FUNCTION pythag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine (adapted from Numerical Recipes) inverts
! the given matrix a using LU decomposition.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE luinv(n,a)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,n
INTEGER, DIMENSION(n) :: indx
REAL(KIND=i5), DIMENSION (n,n) :: a,y
DO i=1,n
   DO j=1,n
      y(i,j)=0.0
   ENDDO
   y(i,i)=1.0
ENDDO
CALL ludcmp(a,n,indx)
DO j=1,n
   CALL lubksb(a,n,indx,y(1,j))
ENDDO
DO i=1,n
   DO j=1,n
      a(i,j)=y(i,j)
   ENDDO
ENDDO
END SUBROUTINE luinv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine (adapted from Numerical Recipes) replaces
! a given matrix by its LU decomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ludcmp(a,n,indx)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,n,imax
INTEGER, DIMENSION(n) :: indx
REAL(KIND=i5), DIMENSION (n) :: vv
REAL(KIND=i5), DIMENSION (n,n) :: a
REAL(KIND=i5) :: aamax,sum,dum
REAL(KIND=i5), PARAMETER :: tiny=1.0e-20
DO i=1,n
   aamax=0.0
   DO j=1,n
      IF(ABS(a(i,j)).GT.aamax)aamax=ABS(a(i,j))
   enddo
   IF(aamax.EQ.0.0)THEN
      WRITE(6,*)'Singular matrix. Try C.G. method'
      STOP
   ENDIF
   vv(i)=1.0/aamax
ENDDO
DO j=1,n
   DO i=1,j-1
      sum=a(i,j)
      DO k=1,i-1
         sum=sum-a(i,k)*a(k,j)
      ENDDO
      a(i,j)=sum
   ENDDO
   aamax=0.0
   DO i=j,n
      sum=a(i,j)
      DO k=1,j-1
         sum=sum-a(i,k)*a(k,j)
      ENDDO
      a(i,j)=sum
      dum=vv(i)*ABS(sum)
      IF(dum.GE.aamax)THEN
         imax=i
         aamax=dum
      ENDIF
   ENDDO
   IF(j.NE.imax)THEN
      DO k=1,n
         dum=a(imax,k)
         a(imax,k)=a(j,k)
         a(j,k)=dum
      ENDDO
      vv(imax)=vv(j)
   ENDIF
   indx(j)=imax
   IF(a(j,j).EQ.0.0)a(j,j)=tiny
   IF(j.NE.n)THEN
      dum=1.0/a(j,j)
      DO i=j+1,n
         a(i,j)=a(i,j)*dum
      ENDDO
   ENDIF
ENDDO
END SUBROUTINE ludcmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine (adapted from Numerical Recipes) solves the
! set of n linear equations AX=B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lubksb(a,n,indx,b)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,n,ii,ll
INTEGER, DIMENSION (n) :: indx
REAL(KIND=i5), DIMENSION (n) :: b
REAL(KIND=i5), DIMENSION (n,n) :: a
REAL(KIND=i5) :: sum
ii=0
DO i=1,n
   ll=indx(i)
   sum=b(ll)
   b(ll)=b(i)
   IF(ii.NE.0)THEN
      DO j=ii,i-1
         sum=sum-a(i,j)*b(j)
      ENDDO
   ELSE IF(sum.NE.0.0)THEN
      ii=1
   ENDIF
   b(i)=sum
ENDDO
DO i=n,1,-1
   sum=b(i)
   DO j=i+1,n
      sum=sum-a(i,j)*b(j)
   ENDDO
   b(i)=sum/a(i,i)
ENDDO
END SUBROUTINE lubksb
