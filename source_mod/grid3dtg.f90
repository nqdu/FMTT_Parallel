!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to generate a 3-D grid of regularly
! spaced velocity vertices in spherical coordinates. Velocity
! is allowed to vary with depth only. Velocity vertices are
! to be interpolated by cubic B-splines. The will also add, if
! required, random structure and/or a checkerboard.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM gridder
IMPLICIT NONE
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER :: i,j,k,l,ars,umocvg,checkstat
INTEGER :: nvr,nvt,nvp,pors,aapmc,adch,nch,usp
INTEGER :: vusp,vusp1,vusp2,vuspo,vusp1o,vusp2o
INTEGER :: asp,nsp
INTEGER, DIMENSION (100) :: ispr, ispt,ispp
REAL(KIND=i10) :: gor,got,gop,gsr,gst,gsp,rdum
REAL(KIND=i10) :: vtop,vbot,vgr,vel,rssf,velp
REAL(KIND=i10) :: d1,d2,v1,v2,rdep,decm,mpc
REAL(KIND=i10) :: chvp,chvp1,chvp2
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: velm
REAL(KIND=i10), DIMENSION (100) :: spa,spr,spt,spp
CHARACTER (LEN=20) :: ofile,mfile
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
! ofile = Name of output grid file
! mfile = Name of velocity model file
! vtop,vbot = Velocity at top and bottom of grid
! vgr = Vertical velocity gradient
! vel = velocity of grid node
! ars = Add random structure (0=no,1=yes)
! rssf = Random structure scaling factor
! velp = Velocity perturbation obtained randomly
! umocvg = Use velocity model (0) or constant gradient (1)
! velm = Velocity discretized at depth from a 1-D model
! pors = P (0) or S (1) velocity model
! r1,r2 = radius at which v1,v2 occur for 1-D model
! v1,v2 = model velocities at r1,r2
! rdep = Required depth at which node requires model velocity
! aapmc = Add a priori model covariance (0=no, 1=yes)
! decm = Diagonal elements of covariance matrix
! adch = Add checkerboard (0=no,1=yes)
! mpc = Maximum perturbation of checkerboard
! nch = size of checkerboard cell
! chvp = Checkerboard velocity perturbation
! chvp1,chvp2 = Dummy checkerboard variables
! usp = use spacing for checkerboard (0=no,1=yes)
! vusp = Checkerboard spacing variable
! vusp1,vusp2 = Dummy spacing variables
! vuspo,vusp1o,vusp2o = Previous values of above
! asp = Apply spike (0=no,1=yes)
! nsp = Number of spikes
! spa = Amplitude of spikes
! spr,spr,spp = Coordinates of spikes
! ispr,ispt,ispp = Grid locations of spikes
! 
OPEN(UNIT=10,FILE='grid3dtg.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)nvr ! no. of r direction
READ(10,*)nvt ! no. of latitude
READ(10,*)nvp ! no. of longitude
READ(10,*)gor ! origin of r
READ(10,*)got ! origin of theta
READ(10,*)gop ! origin of phi
READ(10,*)gsr ! dr
READ(10,*)gst ! dtheta
READ(10,*)gsp ! dphi
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)umocvg ! model(0) or gradient(1)
READ(10,*)pors   ! P or S
READ(10,1)mfile  ! model file
READ(10,*)vtop
READ(10,*)vbot
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)ars
READ(10,*)rssf
READ(10,*)aapmc
READ(10,*)decm
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)adch  ! ad check
READ(10,*)mpc   ! max perturbation velocity
READ(10,*)nch   ! cb size
READ(10,*)usp   ! use spacing
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)asp
READ(10,*)nsp
DO i=1,nsp
   READ(10,*)spa(i)
   READ(10,*)spr(i),spt(i),spp(i)
ENDDO
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)ofile
CLOSE(10)
!
! If spikes are required, compute grid locations
!
IF(asp.EQ.1)THEN
   DO i=1,nsp
      ispr(i)=(gor-spr(i))/gsr+1
      ispt(i)=(got-spt(i))/gst+1
      ispp(i)=(spp(i)-gop)/gsp+1
   ENDDO
ENDIF
1 FORMAT(a20)
!
! Determine velocity gradient (umocvg=1)
! or velocity values from model (umocvg=0)
!
ALLOCATE(velm(nvr+2), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM gridder: REAL velm'
ENDIF
IF(umocvg.EQ.1)THEN
   vgr=(vbot-vtop)/((nvr-1)*gsr)
ELSE
!
!  Open file containing 1-D reference model and read in first
!  two values.
!
   OPEN(UNIT=10,FILE=mfile,STATUS='old')
   IF(pors.EQ.0)THEN
      READ(10,*)d1,v1
      READ(10,*)d2,v2
   ELSE
      READ(10,*)d1,rdum,v1
      READ(10,*)d2,rdum,v2
   ENDIF
!
!  Make sure that the grid surface does not lie
!  above the 1-D model
!
   rdep=-(gor+gsr)
   IF(rdep.LT.d1)THEN
      WRITE(6,*)'Top of 3-D grid lies above radial extent of 1-D model!!'
      WRITE(6,*)'Terminating program'
      STOP
   ENDIF
!
!  Now calculate the velocity at each depth node
!
   DO i=1,nvr+2
      rdep=-(gor+gsr-(i-1)*gsr)
!
!     Ensure that rdep lies within the bounds of d1 and d2
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
!     Now calculate the velocity at the specified depth
!
      velm(i)=v1+(v2-v1)*(rdep-d1)/(d2-d1)
   ENDDO
   CLOSE(10)
ENDIF
!
! Now write output to file, noting that
! we automatically add a cushion of
! boundary nodes to the grid.
!
OPEN(UNIT=20,FILE=ofile,STATUS='unknown')
WRITE(20,*)nvr,nvt,nvp
WRITE(20,'(3f12.6)')gor,got,gop
WRITE(20,'(3f12.6)')gsr,gst,gsp
WRITE(20,'(1X)')
IF(adch.EQ.1)THEN
   chvp1=mpc
   chvp2=mpc
   chvp=mpc
   vusp1=-1
   vusp2=-1
   vusp=-1
ENDIF
DO i=0,nvp+1
   IF(adch.EQ.1)THEN
      IF(MOD(i,nch).EQ.0)THEN
         chvp1=-chvp1
         IF(usp.EQ.1)THEN
            IF(vusp1.EQ.0)THEN
               IF(vusp1o.EQ.-1)THEN
                  vusp1=1
               ELSE
                  vusp1=-1
               ENDIF
            ELSE
               vusp1o=vusp1
               vusp1=0
            ENDIF
         ENDIF
      ENDIF
      chvp2=chvp1
      vusp2=1
      vusp2o=1
   ENDIF
   DO j=0,nvt+1
      IF(adch.EQ.1)THEN
         IF(MOD(j,nch).EQ.0)THEN
            chvp2=-chvp2
            IF(usp.EQ.1)THEN
               IF(vusp2.EQ.0)THEN
                  IF(vusp2o.EQ.-1)THEN
                     vusp2=1
                  ELSE
                     vusp2=-1
                  ENDIF
               ELSE
                  vusp2o=vusp2
                  vusp2=0
               ENDIF
            ENDIF
         ENDIF
         chvp=chvp2
         vusp=1
         vuspo=1
      ENDIF
      DO k=0,nvr+1
         IF(umocvg.EQ.1)THEN
            vel=vtop+gsr*(k-1)*vgr
         ELSE
            vel=velm(k+1)
         ENDIF
!
!        Add random structure if required.
!
         IF(ars.Eq.1)THEN
            CALL RANDOM_NUMBER(velp)
            vel=vel+(velp-0.5)*rssf
         ENDIF
!
!        Add checkerboard if required
!
         IF(adch.EQ.1)THEN
            IF(MOD(k,nch).EQ.0)THEN
               chvp=-chvp
               IF(usp.EQ.1)THEN
                  IF(vusp.EQ.0)THEN
                     IF(vuspo.EQ.-1)THEN
                        vusp=1
                     ELSE
                        vusp=-1
                     ENDIF
                  ELSE
                     vuspo=vusp
                     vusp=0
                  ENDIF
               ENDIF
            ENDIF
            vel=vel+vusp1*vusp2*vusp*chvp
         ENDIF
!
!        Apply spikes if required
!
         IF(asp.EQ.1)THEN
            DO l=1,nsp
               IF(i.EQ.ispp(l).OR.i.EQ.ispp(l)+1)THEN
                  IF(j.EQ.ispt(l).OR.j.EQ.ispt(l)+1)THEN
                     IF(k.EQ.ispr(l).OR.k.EQ.ispr(l)+1)THEN
                        vel=vel+spa(l)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
!
!        Write out specified covariance value
!        if required
!
         IF(aapmc.EQ.1)THEN
            WRITE(20,'(2f12.8)')vel,decm
         ELSE
            WRITE(20,'(f12.8)')vel
         ENDIF
      ENDDO
      WRITE(20,'(1X)')
   ENDDO
   WRITE(20,'(1X)')
ENDDO
CLOSE(20)
DEALLOCATE(velm, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM gridder: velm'
ENDIF
STOP
END PROGRAM gridder
