!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! use in any subroutine or function or other module.
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
IMPLICIT NONE

! precision defined
INTEGER, PARAMETER :: i10= SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: i6 = SELECTED_REAL_KIND(6,20)

!nqdu
! self-defined struct
type,public :: FmttPair
  integer :: ns
  integer,ALLOCATABLE  :: valid_rays(:)
end type

INTEGER, SAVE :: nnr,nnth,nnph,nstth,nstph
INTEGER, SAVE :: ndl,pors
REAL(KIND=i10), SAVE :: gor,goth,goph,dnr,dnth,dnph
REAL(KIND=i10), SAVE :: earth,ddep
REAL(KIND=i6), SAVE :: sgoth,sgoph,sdnth,sdnph
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE, SAVE :: veld
REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE, SAVE :: tgbs,thgbs,phgbs
REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE, SAVE :: ttgb
REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE, SAVE :: ttgs,azgs,iags
REAL(KIND=i10), PARAMETER :: pi=3.1415926535898
type(FmttPair),ALLOCATABLE    :: pairs(:)
!
! nnr,nnth,nnph = Number of nodes in radius, theta, phi
! gor,goth,goph = Origin of grid (radius, theta, phi)
! dnr,dnth,dnph = Node separation in radius, theta, phi
! earth = Radius of Earth (in km)
! ndl = Number of depth layers for discretized 1-D model
! veld = Velocity in each layer of 1-D model
! ddep = Depth discretization of 1-D model (in km)
! pors = P (0) or S (1) wavespeed reference model
! tgbs = Traveltime to grid base for ray i,j
! thgbs = Theta value at grid base for ray i,j
! phgbs = Phi value at grid base for ray i,j
! ttgs = Traveltime of ray i,j to grid surface
! ttgb = True traveltime to node i,j, at grid base
! azgs = Azimuth from grid surface for ray i,j
! iags = Inclination angle to vertical for ray i,j
! nstth,nstph = Number of surface traveltimes in lat,lon
! sgoth,sgoph = Origin of surface ray grid (lat, lon)
! sdnth,sdnph = Surface ray grid spacing (lat,lon)
!

contains
subroutine FmttPair_dealloc(pair)
   implicit none
   type(FmttPair) :: pair

   DEALLOCATE(pair%valid_rays)

end subroutine FmttPair_dealloc

END MODULE globalp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed the reference 1-D velocity model
! and discretizes it as a series of constant velocity
! spherical shells of thickness ddep. Linear interpolation
! is used to determine velocities.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE vellay(velfile)
USE globalp
IMPLICIT NONE
INTEGER :: i
REAL(KIND=i10) :: rdep,d1,d2,v1,v2,rdum
CHARACTER (LEN=20) :: velfile
!
! velfile = File containing 1-D reference model
! rdep = Required depth at which to calculate velocity
! d1,d2 = Bounds of a given depth interval
! v1,v2 = Corresponding velocities at d1,d2
! rdum = dummy variable
!
!
! Determine the number of model layers
!
ndl=INT(dnr*(nnr-1)/ddep-gor)
!
! Allocate memory to layered 1-D velocity field
!
ALLOCATE(veld(ndl))
!
! Open file containing 1-D reference model and read in first
! two values.
!
OPEN(UNIT=10,FILE=velfile,STATUS='old')
IF(pors.EQ.0)THEN
   READ(10,*)d1,v1
   READ(10,*)d2,v2
ELSE
   READ(10,*)d1,rdum,v1
   READ(10,*)d2,rdum,v2
ENDIF
!
! Make sure that the grid surface does not lie
! above the 1-D model
!
rdep=-ddep/2.0 
IF(rdep.LT.d1)THEN
   WRITE(6,*)'Top of 3-D grid lies above radial extent of 1-D model!!'
   WRITE(6,*)'Terminating program'
   STOP
ENDIF
!
! Now calculate the velocity of each layer
!
DO i=1,ndl
   rdep=ddep/2.0+(i-1)*ddep
!
!  Ensure that rdep lies within the bounds of d1 and d2
!
   DO
      IF(rdep.GT.d2)THEN
         d1=d2
         v1=v2
         IF(pors.EQ.0)THEN
           READ(10,*)d2,v2
         ELSE
            READ(10,*)d2,rdum,v2
         ENDIF
      ELSE
         EXIT
      ENDIF
   ENDDO
!
!  Now calculate the velocity at the specified depth
!
   veld(i)=v1+(v2-v1)*(rdep-d1)/(d2-d1)
ENDDO
CLOSE(10)
END SUBROUTINE vellay

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine traces rays from surface grid points to the
! base of the local 3-D model by using Snell's law in
! spherical geometry. Traveltime is also calculated.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE rtrace(srcid,arrid)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,srcid,arrid
REAL(KIND=i10) :: r1,r2,v1,v2,omeg1d,omeg2,omeg2d
REAL(KIND=i10) :: delta,tt,delta1,delta2,rgb,rdum1,rdum2
!
! r1,r2 = Radius of two bounding interfaces
! v1,v2 = Velocity below r1 and r2
! omega1d = Ray inclination beneath r1
! omeg2, omeg2d = Ray inclination above and below r2
! delta1, delta2 = Angular distance from source at r1 and r2
! tt = traveltime of ray from surface
! rgb = Radial distance to grid base
! srcid = id number of current source
! arrid = id number of receiver array
!
! Begin a loop through all the grid points to determine ray
! paths and travel times.
!
DO i=1,nstph ! for each  lon
   DO j=1,nstth ! for each lat
      omeg2d=iags(j,i)
      r2=earth
      delta2=0.0
      tt=0.0
      rgb=earth+gor-(nnr-1)*dnr
!
!     Begin a loop over all the layers
!
      DO k=1,ndl-1
         r1=r2
         r2=earth-k*ddep
         v1=veld(k)
         v2=veld(k+1)
         delta1=delta2
         omeg1d=omeg2d
!
!        Calculate omeg2
!
         omeg2=r1*SIN(omeg1d)/r2
         IF(ABS(omeg2).GT.1.0)THEN
            WRITE(6,*)'Ray path bottoms out above grid base'
            WRITE(6,*)'Check proximity of source to receivers'
            WRITE(6,*)'Problem source is number ',srcid
            WRITE(6,*)'Problem receiver arrays is ',arrid
            ! time(nstth,nstph,nsrc)
            pairs(arrid)%valid_rays(srcid) = 0 ! 
            return;
            WRITE(6,*)'Terminating program itimes'
            STOP
         ENDIF
         omeg2=ASIN(omeg2)
!
!        Calculate delta2
!
         delta2=delta1+omeg2-omeg1d
!
!        Calculate omeg2d
!
         omeg2d=ASIN(v2*r1*sin(omeg1d)/(v1*r2))
!
!        Calculate traveltime from surface
! 
         tt=tt+SQRT(r1**2+r2**2-2.0*r1*r2*COS(delta2-delta1))/v1
      ENDDO
!
!     Now calculate angular distance and traveltime to base 
!     of grid
!
      r1=r2
      r2=rgb
      v1=v2
      delta1=delta2
      omeg1d=omeg2d
      omeg2=ASIN(r1*SIN(omeg1d)/r2)
      delta2=delta1+omeg2-omeg1d
      tt=tt+SQRT(r1**2+r2**2-2.0*r1*r2*COS(delta2-delta1))/v1
      delta=delta2
      tgbs(j,i)=ttgs(j,i)-tt
!
!     Now calculate (theta,phi) location of ray end point
!
      rdum1=(sgoth-(j-1)*sdnth)*pi/180.0
      thgbs(j,i)=SIN(rdum1)*COS(delta1)+COS(rdum1)*COS(azgs(j,i))*SIN(delta1)
      IF(thgbs(j,i).LT.-1.0)THEN
         thgbs(j,i)=-1.0
      ELSE IF(thgbs(j,i).GT.1.0)THEN
         thgbs(j,i)=1.0
      ENDIF
      thgbs(j,i)=ASIN(thgbs(j,i))
      rdum2=(sgoph+(i-1)*sdnph)*pi/180.0
      phgbs(j,i)=COS(delta1)-SIN(thgbs(j,i))*SIN(rdum1)
      phgbs(j,i)=phgbs(j,i)/(COS(thgbs(j,i))*COS(rdum1))
      IF(phgbs(j,i).LT.-1.0)THEN
         phgbs(j,i)=-1.0
      ELSE IF(phgbs(j,i).GT.1.0)THEN
         phgbs(j,i)=1.0
      ENDIF
      phgbs(j,i)=ACOS(phgbs(j,i))
      IF(azgs(j,i).GE.0.0.AND.azgs(j,i).LE.pi)THEN
         phgbs(j,i)=ABS(phgbs(j,i))
      ELSE
         phgbs(j,i)=-ABS(phgbs(j,i))
      ENDIF
      phgbs(j,i)=rdum2+phgbs(j,i)
   ENDDO
ENDDO
END SUBROUTINE rtrace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine interpolates the traveltimes of the ray
! end points onto the base of the 3-D velocity grid. Linear
! interpolation is used via triangular cells of constant
! traveltime gradient (defined parametrically).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE interp
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,l,m
INTEGER :: i1,i2,i3,j1,j2,j3,imin,imax,jmin,jmax
REAL(KIND=i10) :: tha,thb,thc,pha,phb,phc
REAL(KIND=i10) :: u,v,ud,un,vd,vn
REAL(KIND=i10) :: rgoth,rgoph,rdnth,rdnph
!
! tha,thb,thc = Differential theta values for triangle nodes
! pha,phb,phc = Differential phi values for triangle nodes
! i1,i2,i3 = Nearest grid points in phi to triangle nodes
! j1,j2,j3 = Nearest grid points in theta to triangle nodes
! imin,imax = Bounding grid points in phi for triangle
! jmin,jmax = bounding grid points in theta for triangle
! u,v = independent variables of parametric triangle
! ud,un = denominator and numerator of formula for u
! vd,vn = denominator and numerator of formula for v
! rgoth,rgoph,rdnth,rdnph = goth,goph,dnth,dnph in radians
!
rgoth=goth*pi/180.0
rgoph=goph*pi/180.0
rdnth=dnth*pi/180.0
rdnph=dnph*pi/180.0
!
! Set all traveltimes to grid points at base of 3-D grid to
! an arbitrary negative value.
!
  ttgb=-100.0
!
! Start a loop through all the triangular cells
!
DO i=1,nstph-1
   DO j=1,nstth-1
      DO k=0,1
         tha=thgbs(j,i+k)-rgoth
         thb=thgbs(j+1,i)-thgbs(j,i+k)
         thc=thgbs(j+k,i+1)-thgbs(j,i+k)
         pha=phgbs(j,i+k)-rgoph
         phb=phgbs(j+1,i)-phgbs(j,i+k)
         phc=phgbs(j+k,i+1)-phgbs(j,i+k)
!
!        Find the the coordinates of the three points
!        closest to the vertices of the triangle
!
         i1=NINT(pha/rdnph)+1
         j1=-NINT(tha/rdnth)+1
         i2=NINT((pha+phc)/rdnph)+1
         j2=-NINT((tha+thc)/rdnth)+1
         i3=NINT((pha+phb)/rdnph)+1
         j3=-NINT((tha+thb)/rdnth)+1
!
!        Now identify the bounding points of
!        the local rectangular grid
!        
         imin=i1
         IF(i2.LT.imin)imin=i2
         IF(i3.LT.imin)imin=i3
         imax=i1
         IF(i2.GT.imax)imax=i2
         IF(i3.GT.imax)imax=i3
         jmin=j1
         IF(j2.LT.jmin)jmin=j2
         IF(j3.LT.jmin)jmin=j3
         jmax=j1
         IF(j2.GT.jmax)jmax=j2
         IF(j3.GT.jmax)jmax=j3
!
!        If local grid does not overlap global
!        grid then move to the next triangle
!
         IF(imax.LT.1.OR.jmax.LT.1)THEN
            CYCLE
         ELSE IF(imin.GT.nnph.OR.jmin.GT.nnth)THEN
            CYCLE
         ENDIF
!
!        At this point, at least part of the local grid
!        overlaps the global grid. If necessary,
!        adjust the local grid so that it fits completely
!        within the global grid.
!
         IF(imin.LT.1)imin=1
         IF(jmin.LT.1)jmin=1
         IF(imax.GT.nnph)imax=nnph
         IF(jmax.GT.nnth)jmax=nnth
!
!        Now loop through the points of the local grid and
!        calculate their u and v values. If 0.le.u+v.le.1,
!        then the grid point lies within the triangle
!        and interpolation can occur.
!
         ud=phc*thb-thc*phb
         vd=thc*phb-phc*thb
         un=pha*thc-tha*phc
         vn=pha*thb-tha*phb
         DO l=imin,imax
            DO m=jmin,jmax
               u=(un-(m-1)*rdnth*phc-(l-1)*rdnph*thc)/ud
               v=(vn-(m-1)*rdnth*phb-(l-1)*rdnph*thb)/vd
               IF(u.LT.0.0.OR.v.LT.0.0)THEN
                  CYCLE
               ELSE IF(u+v.GT.1.0)THEN
                  CYCLE
               ENDIF
!
!              At this point, the grid node lies within
!              the triangle, so we can calculate the
!              interpolated traveltime value.
!
               ttgb(m,l)=tgbs(j,i+k)+(tgbs(j+1,i)-tgbs(j,i+k))*u
               ttgb(m,l)=ttgb(m,l)+(tgbs(j+k,i+1)-tgbs(j,i+k))*v
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to backtrack traveltimes calculated
! at the surface of the Earth to the base of a local 3-D grid
! of nodes. Programs such as ttimes can be easily used to
! calculate traveltimes from a distant source at depth to
! a receiver located on the surface; however, they normally
! do not track rays to points located below the surface.
! In teleseismic tomography, this is what is required. Any
! global 1-D reference model can be provided, but ak135 is
! used by default
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM itimes
USE globalp
IMPLICIT NONE
INTEGER :: ii,i,j,k,nsrc,nrfr,nrfth,nrfph,nra
integer :: idum
REAL(KIND=i6) :: spttgs,spazgs,spiags
CHARACTER (LEN=25) :: velfile,srcfile,grdfile
CHARACTER (LEN=25) :: stmfile,otpfile,pfile
!
! velfile = File containing 1-D reference model
! srcfile = File containing list of source info
! grdfile = File containing 3-D velocity grid info
! stmfile = File containing 1-D times at surface
! otpfile = Output file containing times to base of model
! pfile = FMM parameter file
! ddep = Depth discretization of 1-D model (in km)
! nsrc = Number of sources
! spttgs,spazgs,spiags = Single precision forms of ttgs,azgs,iags
! nrfr,nrfth,nrfph = B-spline dicing in r, theta, phi
! nra = Number of receiver arrays
!
! Start by reading in parameter file
!
OPEN(UNIT=10,FILE='itimes.in',STATUS='old')
READ(10,1)velfile
READ(10,1)srcfile
READ(10,1)grdfile
READ(10,1)pfile !fm3dt.in
READ(10,1)stmfile !aktsurf.out
READ(10,1)otpfile !itimes.out
READ(10,*)pors
READ(10,*)ddep
READ(10,*)earth
1   FORMAT(a25)
CLOSE(10)
!
! Read in the number of event types (i.e. the number of receiver
! arrays).
!
OPEN(UNIT=30,FILE=srcfile,STATUS='old')
READ(30,*)nra
ALLOCATE(pairs(nra))
!
! Read in the 3-D velocity grid information
!
OPEN(UNIT=10,FILE=grdfile,STATUS='old')
READ(10,*)nnr,nnth,nnph
READ(10,*)gor,goth,goph
READ(10,*)dnr,dnth,dnph
CLOSE(10)
!
! Read in parameter file information
!
OPEN(UNIT=10,FILE=pfile,STATUS='old')
DO i=1,7
   READ(10,*)
ENDDO
READ(10,*)nrfr,nrfth,nrfph ! dicing factor
CLOSE(10)
nnr=(nnr-1)*nrfr+1
nnth=(nnth-1)*nrfth+1
nnph=(nnph-1)*nrfph+1  ! refined grid
dnr=dnr/nrfr
dnth=dnth/nrfth
dnph=dnph/nrfph
!
! Call a subroutine that formulates the discretized
! 1-D model
!
CALL vellay(velfile)
!
! Allocate memory to arrays describing traveltime
! to underside of grid, location of ray end point,
! and traveltime to actual grid points at the base
! of the grid.
!
ALLOCATE(ttgb(nnth,nnph))
!
! Open stmfile and otpfile and begin a loop through
! each station.
!
OPEN(UNIT=10,FILE=stmfile,FORM='unformatted',STATUS='old') !aktsurf.out
OPEN(UNIT=20,FILE=otpfile,FORM='unformatted',STATUS='unknown')
DO ii=1,nra ! loop around each array
   READ(10)sgoth,sgoph
   READ(10)nstth,nstph !nlat,nlon
   READ(10)sdnth,sdnph
   ALLOCATE(ttgs(nstth,nstph))
   ALLOCATE(azgs(nstth,nstph))
   ALLOCATE(iags(nstth,nstph))
   ALLOCATE(tgbs(nstth,nstph))
   ALLOCATE(thgbs(nstth,nstph))
   ALLOCATE(phgbs(nstth,nstph))
   READ(30,*)nsrc
   pairs(ii)%ns = nsrc
   ALLOCATE(pairs(ii)%valid_rays(nsrc))
   pairs(ii)%valid_rays(:) = 1
   DO i=1,nsrc
      READ(30,*)
      READ(30,*)
   ENDDO
   DO i=1,nsrc ! loop around sources
!
!     Read in traveltimes, azimuths and inclination angles
!     to all surface grid nodes for source i
!
      DO j=1,nstph
         DO k=1,nstth
            READ(10)spttgs,spazgs,spiags !time, ray azimuth and incidence angle
            ttgs(k,j)=spttgs
            azgs(k,j)=spazgs
            iags(k,j)=spiags
            azgs(k,j)=azgs(k,j)*pi/180.0
            iags(k,j)=iags(k,j)*pi/180.0
         ENDDO
      ENDDO
!
!     Call a subroutine which calculates traveltimes and ray
!     end points
!
      CALL rtrace(i,ii)
      if(pairs(ii)%valid_rays(i) == 1) then
         !  
         !Finally, call a subroutine which interpolates the
         !traveltimes of the ray end points onto the nodes
         !of the grid.
         !
         CALL interp
         !
         !     Write output to output file. Note that negative
         !     traveltime values denote nodes which do not
         !     have traveltime values associated with them.
         !
         idum = 0
         DO j=1,nnph
            DO k=1,nnth
               if (ttgb(k,j) > 0.) idum = idum + 1
            ENDDO
         ENDDO
         if(idum > 0) then 
            DO j=1,nnph
               DO k=1,nnth
                  WRITE(20)ttgb(k,j)
               ENDDO
            ENDDO
         else
            pairs(ii)%valid_rays(i) = 0
         endif
      endif

   ENDDO
   DEALLOCATE(ttgs,azgs,iags)
   DEALLOCATE(tgbs,thgbs,phgbs)
ENDDO
CLOSE(10)
CLOSE(20)
CLOSE(30)

! save valid rays
open(10,file='valid_rays.dat')
!write(10,*)nra
do ii=1,nra
   !WRITE(10,*)pairs(ii)%ns
   do i=1,pairs(ii)%ns
      WRITE(10,*) pairs(ii)%valid_rays(i)
   enddo
enddo
!
! Deallocate array memory
!
DEALLOCATE(veld,ttgb)
do i=1,nra
   call FmttPair_dealloc(pairs(i))
enddo
DEALLOCATE(pairs)
END PROGRAM itimes
