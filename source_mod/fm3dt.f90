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

!nqdu
! self-defined struct
type,public :: FmttPair
  integer :: ns
  integer,ALLOCATABLE  :: valid_rays(:)
end type

INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER :: checkstat
INTEGER, SAVE :: nnr,nnx,nnz,nrc,fom,ltfr
INTEGER, SAVE :: nvr,nvx,nvz,nrfr,nrfx,nrfz
INTEGER, SAVE :: nsnn,nsns,nsne,nsnw
REAL(KIND=i10), SAVE :: gor,gox,goz,dnr,dnx,dnz,snb,earth
REAL(KIND=i10), SAVE :: dvr,dvx,dvz
REAL(KIND=i10), SAVE :: goxd,gozd,dnxd,dnzd
REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: veln
REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: ttn 
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE, SAVE :: rcr,rcx,rcz
REAL(KIND=i10), PARAMETER :: pi=3.1415926535898
type(FmttPair),ALLOCATABLE :: pairs(:)
!
! nnr,nnx,nnz = Number of nodes of refined grid in r,x and z
! nvr,nvx,nvz = Number of B-spline vertices in r,x and z
! gor,gox,goz = Origin of grid (radius,theta,phi)
! dnr,dnx,dnz = Node separation of refined grid in r, x and z
! dvr,dvx,dvz = Node separation of B-spline grid in r, x and z
! nrfr,nrfx,nrfz = B-spline dicing level in r, x and z
! veln(i,j,k) = velocity values on a refined grid of nodes
! ttn(i,j,k) = traveltime field on the refined grid of nodes
! checkstat = check status of memory allocation
! fom = use first-order(0) or mixed-order(1) scheme
! snb = Maximum size of narrow band as fraction of nnx*nnz
! nrc = number of receivers
! rcr(i),rcx(i),rcz(i) = (r-earth,x,z) coordinates of receivers
! earth = radius of Earth (in km)
! goxd,gozd = gox,goz in radians
! dnzd,dnzd = dnx,dnz in radians
! dfr,dfx,dfz = B-spline dicing factor in r, theta and phi
! ltfr = Limit traveltime field to receiver array (0=no, 1=yes)
! nsnn,nsns = Latitude bounds (N-S) of surface receiver grid
! nsne,nsnw = Longitude bounds (E-W) of surface receiver grid
!
contains
subroutine FmttPair_dealloc(pair)
   implicit none
   type(FmttPair) :: pair

   DEALLOCATE(pair%valid_rays)

end subroutine FmttPair_dealloc

END MODULE globalp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module contains all the subroutines used to calculate
! the first-arrival traveltime field through the grid.
! Subroutines are:
! (1) travel
! (2) fouds1
! (3) fouds2
! (4) addtree
! (5) downtree
! (6) updtree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE traveltime
USE globalp
IMPLICIT NONE
INTEGER ntr
INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: nsts
TYPE backpointer
   INTEGER(KIND=2) :: pr,px,pz
END TYPE backpointer
TYPE(backpointer), DIMENSION (:), ALLOCATABLE :: btg
!
! nsts(i,j,k) = node status (-1=far,0=alive,>0=close)
! btg = backpointer to relate gird nodes to binary tree entries
! pr = grid point in r
! px = grid-point in x
! pz = grid-point in z
! ntr = number of entries in binary tree
!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE travel(srcid)
IMPLICIT NONE
INTEGER :: i,j,ir,ix,iz,maxbt,sw,srcid,ittmin,jttmin
INTEGER :: isum,tnrn
REAL(KIND=i10) :: rd1,ttmin
! ir,ix,iz = i,j,k position of "close" point with minimum traveltime
! maxbt = maximum size of narrow band binary tree
! rd1 = substitution variable
! srcid = id number of source
! sw = switch (0 or 1)
! ittmin,jttmin = grid point of minimum traveltime at grid base
! ttmin = minimum traveltime along grid base
! tnrn = total number of receiver nodes
!
! Allocate nsts and set all elements of array nsts equal to 
! (-1) since to begin with, all points are "far" points.
!
ALLOCATE(nsts(0:nnz+1,0:nnx+1,0:nnr+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE travel: INTEGER nsts'
ENDIF
nsts=-1
IF(ltfr.EQ.1)THEN
   tnrn=(nsns+1-nsnn)*(nsne+1-nsnw)
   ttn=0.0
   isum=0
ENDIF
!
! Allocate memory for binary tree
!
maxbt=nint(snb*nnr*nnx*nnz)
ALLOCATE(btg(maxbt), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE travel: TYPE btg'
ENDIF
!
! set initial size of tree to zero
!
ntr=0
!
! Read in traveltimes calculated a priori to base of grid,
! find point with minimum traveltime,
!
sw=0
DO i=1,nnz
   DO j=1,nnx
      READ(20)rd1
      IF(rd1.GT.0.0)THEN
         
         IF(sw.eq.0)THEN
            ttmin=rd1
            ittmin=i
            jttmin=j
            sw=1
         ELSE IF(rd1.LT.ttmin)THEN
            ttmin=rd1
            ittmin=i
            jttmin=j
         ENDIF
         ttn(i,j,nnr)=rd1
         nsts(i,j,nnr)=0
      ENDIF
   ENDDO
ENDDO
IF(sw.EQ.0)THEN
   WRITE(6,*)'No initial traveltimes for source ',srcid
   WRITE(6,*)'Terminating program'
   STOP
ENDIF
!
! Initiate starting point of the narrow band
!
nsts(ittmin,jttmin,nnr)=-1
CALL addtree(ittmin,jttmin,nnr)
!
! Now calculate the first-arrival traveltimes at the
! remaining grid points. This is done via a loop which
! repeats the procedure of finding the first-arrival
! of all "close" points, adding it to the set of "alive"
! points and updating the points surrounding the new "alive"
! point. The process ceases when the binary tree is empty,
! in which case all grid points are "alive".
!
DO WHILE(ntr.gt.0)
!
! Set the "close" point with minimum traveltime
! to "alive"
!
   ir=btg(1)%pr
   ix=btg(1)%px
   iz=btg(1)%pz
   nsts(iz,ix,ir)=0
!
!  If FMM is restricted to targeting the receiver grid,
!  then check if the receiver grid is full
!
   IF(ltfr.EQ.1)THEN
      IF(ir.EQ.1)THEN
         IF(ix.GE.nsnn.AND.ix.LE.nsns)THEN
            IF(iz.GE.nsnw.AND.iz.LE.nsne)THEN
               isum=isum+1
               IF(isum.EQ.tnrn)EXIT
            ENDIF
         ENDIF
      ENDIF
   ENDIF
!
! Update the binary tree by removing the root and
! sweeping down the tree.
!
   CALL downtree
!
! Now update or find values of up to six grid points
! that surround the new "alive" point.
!
! Test points that vary in r
!
   DO i=ir-1,ir+1,2
      IF(i.ge.1.and.i.le.nnr)THEN
         IF(nsts(iz,ix,i).eq.-1)THEN
!
! This option occurs when a far point is added to the list
! of "close" points
!
            IF(fom.eq.0)THEN
               CALL fouds1(iz,ix,i)
            ELSE
               CALL fouds2(iz,ix,i)
            ENDIF
            CALL addtree(iz,ix,i)
         ELSE IF(nsts(iz,ix,i).gt.0)THEN
!
! This happens when a "close" point is updated
!
            IF(fom.eq.0)THEN
               CALL fouds1(iz,ix,i)
            ELSE
               CALL fouds2(iz,ix,i)
            ENDIF
            CALL updtree(iz,ix,i)
         ENDIF
      ENDIF
   ENDDO
!
! Test points that vary in x
!
   DO i=ix-1,ix+1,2
      IF(i.ge.1.and.i.le.nnx)THEN
         IF(nsts(iz,i,ir).eq.-1)THEN
!
! This option occurs when a far point is added to the list
! of "close" points
!
            IF(fom.eq.0)THEN
               CALL fouds1(iz,i,ir)
            ELSE
               CALL fouds2(iz,i,ir)
            ENDIF
            CALL addtree(iz,i,ir)
         ELSE IF(nsts(iz,i,ir).gt.0)THEN
!
! This happens when a "close" point is updated
!
            IF(fom.eq.0)THEN
               CALL fouds1(iz,i,ir)
            ELSE
               CALL fouds2(iz,i,ir)
            ENDIF
            CALL updtree(iz,i,ir)
         ENDIF
      ENDIF
   ENDDO
!
! Test points that vary in z
!
   DO i=iz-1,iz+1,2
      IF(i.ge.1.and.i.le.nnz)THEN
         IF(nsts(i,ix,ir).eq.-1)THEN
!
! This option occurs when a far point is added to the list
! of "close" points
!
            IF(fom.eq.0)THEN
               CALL fouds1(i,ix,ir)
            ELSE
               CALL fouds2(i,ix,ir)
            ENDIF
            CALL addtree(i,ix,ir)
         ELSE IF(nsts(i,ix,ir).gt.0)THEN
!
! This happens when a "close" point is updated
!
            IF(fom.eq.0)THEN
               CALL fouds1(i,ix,ir)
            ELSE
               CALL fouds2(i,ix,ir)
            ENDIF
            CALL updtree(i,ix,ir)
         ENDIF
      ENDIF
   ENDDO
ENDDO
DEALLOCATE(nsts,btg, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE travel: nsts,btg'
ENDIF
END SUBROUTINE travel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates a trial first-arrival traveltime
! at a given node from surrounding nodes using the
! First-Order Upwind Difference Scheme (FOUDS) of
! Sethian and Popovici (1999).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE fouds1(iz,ix,ir)
IMPLICIT NONE
INTEGER :: i,j,k,ir,ix,iz,tsw1,swsol
REAL(KIND=i10) :: trav,travm,slown,tdsh,tref
REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti
REAL(KIND=i10) :: rd1,rd2,rd3
!
! ir = Radial position of node coordinate for determination
! ix = NS position of node coordinate for determination
! iz = EW vertical position of node coordinate for determination
! trav = traveltime calculated for trial node
! travm = minimum traveltime calculated for trial node
! slown = slowness at (iz,ix,ir)
! tsw1 = traveltime switch (0=first time,1=previously)
! a,b,c,u,v,w,em,en = Convenience variables for solving quadratic
! tdsh = local traveltime from neighbouring node
! tref = reference traveltime at neighbouring node
! ri = Radial distance to ir
! risti = ri*sin(theta) at point (iz,ix,ir)
! rd1,rd2,rd3 = dummy variables
! swsol = switch for solution (0=no solution, 1=solution)
!
! Inspect each of the eight quadrants for the minimum time
! solution.
!
tsw1=0
slown=1.0/veln(iz,ix,ir)
ri=gor-(ir-1)*dnr+earth
risti=ri*sin(gox+(ix-1)*dnx)
DO i=ir-1,ir+1,2
   DO j=ix-1,ix+1,2
      DO k=iz-1,iz+1,2 
         IF(i.GE.1.AND.i.LE.nnr)THEN
            IF(j.GE.1.AND.j.LE.nnx)THEN
               IF(k.GE.1.AND.k.LE.nnz)THEN
!
!                 There are seven solution options in
!                 each quadrant.
!
                  swsol=0
                  IF(nsts(iz,ix,i).EQ.0)THEN
                     swsol=1
                     IF(nsts(iz,j,ir).EQ.0.AND.nsts(k,ix,ir).EQ.0)THEN
                        u=dnr
                        v=ri*dnx
                        w=risti*dnz
                        em=ttn(iz,j,ir)-ttn(iz,ix,i)
                        en=ttn(k,ix,ir)-ttn(iz,ix,i)
                        rd1=v**2*w**2
                        rd2=u**2*w**2
                        rd3=u**2*v**2
                        a=rd1+rd2+rd3
                        b=-2.0*(em*rd2+en*rd3)
                        c=em**2*rd2+en**2*rd3-slown**2*u**2*rd1
                        tref=ttn(iz,ix,i)
                     ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                        u=dnr
                        v=ri*dnx
                        em=ttn(iz,j,ir)-ttn(iz,ix,i)
                        a=u**2+v**2
                        b=-2.0*u**2*em
                        c=u**2*(em**2-v**2*slown**2)
                        tref=ttn(iz,ix,i)
                     ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                        u=dnr
                        v=risti*dnz
                        em=ttn(k,ix,ir)-ttn(iz,ix,i)
                        a=u**2+v**2
                        b=-2.0*u**2*em
                        c=u**2*(em**2-v**2*slown**2)
                        tref=ttn(iz,ix,i)
                     ELSE
                        a=1.0
                        b=0.0
                        c=-slown**2*dnr**2
                        tref=ttn(iz,ix,i)
                     ENDIF
                  ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                     swsol=1
                     IF(nsts(k,ix,ir).EQ.0)THEN
                        u=ri*dnx
                        v=risti*dnz
                        em=ttn(k,ix,ir)-ttn(iz,j,ir)
                        a=u**2+v**2
                        b=-2.0*u**2*em
                        c=u**2*(em**2-v**2*slown**2)
                        tref=ttn(iz,j,ir)
                     ELSE
                        a=1.0
                        b=0.0
                        c=-slown**2*ri**2*dnx**2
                        tref=ttn(iz,j,ir)
                     ENDIF
                  ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                     swsol=1
                     a=1.0
                     b=0.0
                     c=-(slown*risti*dnz)**2
                     tref=ttn(k,ix,ir)
                  ENDIF
!
!                 Now find the solution of the quadratic equation
!
                  IF(swsol.EQ.1)THEN
                     rd1=b**2-4.0*a*c
                     IF(rd1.LT.0.0)rd1=0.0
                     tdsh=(-b+sqrt(rd1))/(2.0*a)
                     trav=tref+tdsh
                     IF(tsw1.EQ.1)THEN
                        travm=MIN(trav,travm)
                     ELSE
                        travm=trav
                        tsw1=1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO
ttn(iz,ix,ir)=travm
END SUBROUTINE fouds1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates a trial first-arrival traveltime
! at a given node from surrounding nodes using the
! Mixed-Order (2nd) Upwind Difference Scheme (FOUDS) of
! Popovici and Sethian (2002).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE fouds2(iz,ix,ir)
IMPLICIT NONE
INTEGER :: i,j,k,i2,j2,k2,ir,ix,iz,tsw1
INTEGER :: swi,swj,swk,swsol
REAL(KIND=i10) :: trav,travm,slown,tdsh,tref,tdiv
REAL(KIND=i10) :: a,b,c,u,v,w,em,en,ri,risti,rd1
!
! ir = Radial position of node coordinate for determination
! ix = NS position of node coordinate for determination
! iz = EW vertical position of node coordinate for determination
! trav = traveltime calculated for trial node
! travm = minimum traveltime calculated for trial node
! slown = slowness at (iz,ix,ir)
! tsw1 = traveltime switch (0=first time,1=previously)
! a,b,c,u,v,w,em,en = Convenience variables for solving quadratic
! tdsh = local traveltime from neighbouring node
! tref = reference traveltime at neighbouring node
! ri = Radial distance to ir
! risti = ri*sin(theta) at point (iz,ix,ir)
! swi,swj,swk = switches for second order operators
! tdiv = term to divide tref by depending on operator order
! swsol = switch for solution (0=no solution, 1=solution)
!
! Inspect each of the eight quadrants for the minimum time
! solution.
!
tsw1=0
slown=1.0/veln(iz,ix,ir)
ri=gor-(ir-1)*dnr+earth
risti=ri*sin(gox+(ix-1)*dnx)
DO i=ir-1,ir+1,2
   swi=-1
   IF(i.eq.ir-1)THEN
      i2=i-1
      IF(i2.GE.1)THEN
         IF(nsts(iz,ix,i2).EQ.0)swi=0
      ENDIF
   ELSE
      i2=i+1
      IF(i2.LE.nnr)THEN
         IF(nsts(iz,ix,i2).EQ.0)swi=0
      ENDIF
   ENDIF
   IF(nsts(iz,ix,i).EQ.0.AND.swi.EQ.0)THEN
      swi=-1
      IF(ttn(iz,ix,i).GT.ttn(iz,ix,i2))THEN
         swi=0
      ENDIF
   ELSE
      swi=-1
   ENDIF
   DO j=ix-1,ix+1,2
      swj=-1
      IF(j.eq.ix-1)THEN
         j2=j-1
         IF(j2.GE.1)THEN
            IF(nsts(iz,j2,ir).EQ.0)swj=0
         ENDIF
      ELSE
         j2=j+1
         IF(j2.LE.nnx)THEN
            IF(nsts(iz,j2,ir).EQ.0)swj=0
         ENDIF
      ENDIF
      IF(nsts(iz,j,ir).EQ.0.AND.swj.EQ.0)THEN
         swj=-1
         IF(ttn(iz,j,ir).GT.ttn(iz,j2,ir))THEN
            swj=0
         ENDIF
      ELSE
         swj=-1
      ENDIF
      DO k=iz-1,iz+1,2
         swk=-1
         IF(k.eq.iz-1)THEN
            k2=k-1
            IF(k2.GE.1)THEN
               IF(nsts(k2,ix,ir).EQ.0)swk=0
            ENDIF
         ELSE
            k2=k+1
            IF(k2.LE.nnz)THEN
               IF(nsts(k2,ix,ir).EQ.0)swk=0
            ENDIF
         ENDIF
         IF(nsts(k,ix,ir).EQ.0.AND.swk.EQ.0)THEN
            swk=-1
            IF(ttn(k,ix,ir).GT.ttn(k2,ix,ir))THEN
               swk=0
            ENDIF
         ELSE
            swk=-1
         ENDIF
         IF(i.GE.1.AND.i.LE.nnr)THEN
            IF(j.GE.1.AND.j.LE.nnx)THEN
               IF(k.GE.1.AND.k.LE.nnz)THEN
!
!                 There are 26 solution options in
!                 each quadrant.
!
                  swsol=0
                  IF(swi.EQ.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(iz,j,ir)
                           em=em+ttn(iz,j2,ir)
                           en=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(k,ix,ir)
                           en=en+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           u=risti*dnz
                           v=2.0*ri*dnx
                           w=2.0*dnr
                           em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn(k,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*dnr
                           v=2.0*ri*dnx
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(iz,j,ir)
                           em=em+ttn(iz,j2,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2) 
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ENDIF
                     ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=3.0*ttn(iz,j,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           u=ri*dnx
                           v=2.0*dnr
                           w=risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           en=ttn(iz,j,ir)-ttn(k,ix,ir)
                           a=w**2*v**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE
                           u=ri*dnx
                           v=2.0*dnr
                           em=3.0*ttn(iz,j,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=2.0*dnr
                           v=2.0*risti*dnz
                           em=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)-4.0*ttn(k,ix,ir)
                           em=em+ttn(k2,ix,ir)
                           a=u**2+v**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           u=risti*dnz
                           v=2.0*dnr
                           em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,ix,i)+ttn(iz,ix,i2)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*dnr
                           a=1.0
                           b=0.0
                           c=-u**2*slown**2
                           tref=4.0*ttn(iz,ix,i)-ttn(iz,ix,i2)
                           tdiv=3.0
                        ENDIF
                     ENDIF
                  ELSE IF(nsts(iz,ix,i).EQ.0)THEN
                     swsol=1
                     IF(swj.EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=2.0*risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+9.0*u**2*v**2
                           b=6.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           u=dnr
                           v=2.0*ri*dnx
                           w=risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           en=ttn(iz,ix,i)-ttn(k,ix,ir)
                           a=v**2*w**2+9.0*u**2*w**2+u**2*v**2
                           b=u**2*(6.0*em*w**2+2.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=2.0*ri*dnx
                           em=3.0*ttn(iz,ix,i)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=2.0*risti*dnz
                           em=ttn(iz,ix,i)-ttn(iz,j,ir)
                           en=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2*w**2+u**2*w**2+9.0*u**2*v**2
                           b=u**2*(2.0*em*w**2+6.0*en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-v**2*w**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           u=dnr
                           v=ri*dnx
                           w=risti*dnz
                           em=ttn(iz,j,ir)-ttn(iz,ix,i)
                           en=ttn(k,ix,ir)-ttn(iz,ix,i)
                           a=v**2*w**2+u**2*w**2+u**2*v**2
                           b=-2.0*u**2*(em*w**2+en*v**2)
                           c=u**2*(em**2*w**2+en**2*v**2-slown**2*v**2*w**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           u=dnr
                           v=ri*dnx
                           em=ttn(iz,j,ir)-ttn(iz,ix,i)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           u=dnr
                           v=2.0*risti*dnz
                           em=3.0*ttn(iz,ix,i)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           u=dnr
                           v=risti*dnz
                           em=ttn(k,ix,ir)-ttn(iz,ix,i)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ELSE
                           a=1.0
                           b=0.0
                           c=-slown**2*dnr**2
                           tref=ttn(iz,ix,i)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ELSE
                     IF(swj.EQ.0)THEN
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=2.0*ri*dnx
                           v=2.0*risti*dnz
                           em=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)-4.0*ttn(k,ix,ir)
                           em=em+ttn(k2,ix,ir)
                           a=v**2+u**2
                           b=2.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)
                           tdiv=3.0
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           u=risti*dnz
                           v=2.0*ri*dnx
                           em=3.0*ttn(k,ix,ir)-4.0*ttn(iz,j,ir)+ttn(iz,j2,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-slown**2*v**2)
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ELSE
                           u=2.0*ri*dnx
                           a=1.0
                           b=0.0
                           c=-u**2*slown**2
                           tref=4.0*ttn(iz,j,ir)-ttn(iz,j2,ir)
                           tdiv=3.0
                        ENDIF
                     ELSE IF(nsts(iz,j,ir).EQ.0)THEN
                        swsol=1
                        IF(swk.EQ.0)THEN
                           u=ri*dnx
                           v=2.0*risti*dnz
                           em=3.0*ttn(iz,j,ir)-4.0*ttn(k,ix,ir)+ttn(k2,ix,ir)
                           a=v**2+9.0*u**2
                           b=6.0*em*u**2
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           u=ri*dnx
                           v=risti*dnz
                           em=ttn(k,ix,ir)-ttn(iz,j,ir)
                           a=u**2+v**2
                           b=-2.0*u**2*em
                           c=u**2*(em**2-v**2*slown**2)
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ELSE
                           a=1.0
                           b=0.0
                           c=-slown**2*ri**2*dnx**2
                           tref=ttn(iz,j,ir)
                           tdiv=1.0
                        ENDIF
                     ELSE
                        IF(swk.EQ.0)THEN
                           swsol=1
                           u=2.0*risti*dnz
                           a=1.0
                           b=0.0
                           c=-u**2*slown**2
                           tref=4.0*ttn(k,ix,ir)-ttn(k2,ix,ir)
                           tdiv=3.0
                        ELSE IF(nsts(k,ix,ir).EQ.0)THEN
                           swsol=1
                           a=1.0
                           b=0.0
                           c=-slown**2*risti**2*dnz**2
                           tref=ttn(k,ix,ir)
                           tdiv=1.0
                        ENDIF
                     ENDIF
                  ENDIF
!
!                 Now find the solution of the quadratic equation
!
                  IF(swsol.EQ.1)THEN
                     rd1=b**2-4.0*a*c
                     IF(rd1.LT.0.0)rd1=0.0
                     tdsh=(-b+sqrt(rd1))/(2.0*a)
                     trav=(tref+tdsh)/tdiv
                     IF(tsw1.EQ.1)THEN
                        travm=MIN(trav,travm)
                     ELSE
                        travm=trav
                        tsw1=1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO
ttn(iz,ix,ir)=travm
END SUBROUTINE fouds2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine adds a value to the binary tree by
! placing a value at the bottom and pushing it up
! to its correct position.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE addtree(iz,ix,ir)
IMPLICIT NONE
INTEGER :: ix,iz,ir,tpp,tpc
TYPE(backpointer) :: exch
!
! ir,ix,iz = grid position of new addition to tree
! tpp = tree position of parent
! tpc = tree position of child
! exch = dummy to exchange btg values
!
! First, increase the size of the tree by one.
!
ntr=ntr+1
!
! Put new value at base of tree
!
nsts(iz,ix,ir)=ntr
btg(ntr)%pr=ir
btg(ntr)%px=ix
btg(ntr)%pz=iz
!
! Now filter the new value up to its correct position
!
tpc=ntr
tpp=tpc/2
DO WHILE(tpp.gt.0)
   IF(ttn(iz,ix,ir).lt.ttn(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr))THEN
      nsts(iz,ix,ir)=tpp
      nsts(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)=tpc
      exch=btg(tpc)
      btg(tpc)=btg(tpp)
      btg(tpp)=exch
      tpc=tpp
      tpp=tpc/2
   ELSE
      tpp=0
   ENDIF
ENDDO
END SUBROUTINE addtree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine updates the binary tree after the root
! value has been used. The root is replaced by the value
! at the bottom of the tree, which is then filtered down
! to its correct position. This ensures that the tree remains
! balanced.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE downtree
IMPLICIT NONE
INTEGER :: tpp,tpc
REAL(KIND=i10) :: rd1,rd2
TYPE(backpointer) :: exch
!
! tpp = tree position of parent
! tpc = tree position of child
! exch = dummy to exchange btg values
! rd1,rd2 = substitution variables
!
! Replace root of tree with its last value
!
IF(ntr.EQ.1)THEN
   ntr=ntr-1
   RETURN
ENDIF
nsts(btg(ntr)%pz,btg(ntr)%px,btg(ntr)%pr)=1
btg(1)=btg(ntr)
!
! Reduce size of tree by one
!
ntr=ntr-1
!
! Now filter new root down to its correct position
!
tpp=1
tpc=2*tpp
DO WHILE(tpc.lt.ntr)
!
! Check which of the two children is smallest - use the smallest
!
   rd1=ttn(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)
   rd2=ttn(btg(tpc+1)%pz,btg(tpc+1)%px,btg(tpc+1)%pr)
   IF(rd1.gt.rd2)THEN
      tpc=tpc+1
   ENDIF
!
!  Check whether the child is smaller than the parent; if so, then swap,
!  if not, then we are done
!
   rd1=ttn(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)
   rd2=ttn(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)
   IF(rd1.lt.rd2)THEN
      nsts(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)=tpc
      nsts(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)=tpp
      exch=btg(tpc)
      btg(tpc)=btg(tpp)
      btg(tpp)=exch
      tpp=tpc
      tpc=2*tpp
   ELSE
      tpc=ntr+1
   ENDIF
ENDDO
!
! If ntr is an even number, then we still have one more test to do
!
IF(tpc.eq.ntr)THEN
   rd1=ttn(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)
   rd2=ttn(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)
   IF(rd1.lt.rd2)THEN
      nsts(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)=tpc
      nsts(btg(tpc)%pz,btg(tpc)%px,btg(tpc)%pr)=tpp
      exch=btg(tpc)
      btg(tpc)=btg(tpp)
      btg(tpp)=exch
   ENDIF
ENDIF
END SUBROUTINE downtree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine updates a value on the binary tree. The FMM
! should only produce updated values that are less than their
! prior values.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE updtree(iz,ix,ir)
IMPLICIT NONE
INTEGER :: ir,ix,iz,tpp,tpc
TYPE(backpointer) :: exch
!
! ir,ix,iz = grid position of new addition to tree
! tpp = tree position of parent
! tpc = tree position of child
! exch = dummy to exchange btg values
!
! Filter the updated value to its correct position
!
tpc=nsts(iz,ix,ir)
tpp=tpc/2
DO WHILE(tpp.gt.0)
   IF(ttn(iz,ix,ir).lt.ttn(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr))THEN
      nsts(iz,ix,ir)=tpp
      nsts(btg(tpp)%pz,btg(tpp)%px,btg(tpp)%pr)=tpc
      exch=btg(tpc)
      btg(tpc)=btg(tpp)
      btg(tpp)=exch
      tpc=tpp
      tpp=tpc/2
   ELSE
      tpp=0
   ENDIF
ENDDO
END SUBROUTINE updtree

END MODULE traveltime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM   
! CODE: FORTRAN 90
! This program is designed to implement the Fast Marching
! Method (FMM) for calculating first-arrival traveltimes
! through a 3-D continuous velocity medium. It is written in
! Fortran 90, although it is probably more accurately
! described as Fortran 77 with some of the Fortran 90
! extensions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM fmmin3d
USE globalp
USE traveltime
IMPLICIT NONE
CHARACTER (LEN=20) :: sources,itimes,receivers,gridv
CHARACTER (LEN=20) :: travelt,rtravel,wrays,cdum
CHARACTER (LEN=20) :: frechet,pht
INTEGER :: ii,i,j,k,l,nsrc,wttf,fsrt,wrgf,cfd,tnr,nra,idum
INTEGER :: awttf,awrgf
REAL(KIND=i10) :: cslat,cslong,abw,abe,abn,abs
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: scr,scx,scz
!
! sources = File containing source locations
! itimes = File containing initial traveltimes
! receivers = File containing receiver locations
! gridv = File containing grid of velocity vertices
! travelt = File name for storage of traveltime field
! wttf = Write traveltimes to file? (0=no,>0=source id)
! awttf = Receiver array ID for traveltimes
! fom = Use first-order(0) or mixed-order(1) scheme
! nsrc = number of sources
! scr,scx,scz = source location in r,x,z
! fsrt = find source-receiver traveltimes? (0=no,1=yes)
! rtravel = output file for source-receiver traveltimes
! cdum = dummy character variable
! wrgf = write ray geometries to file? (<0=all,0=no,>0=source id)
! awrgf = Receiver array ID for ray paths
! wrays = file containing raypath geometries
! cfd = calculate Frechet derivatives? (0=no,1=yes)
! frechet = output file containing matrix of frechet derivatives
! pht = phase type (e.g. P, S, pP etc.)
! tnr = total number of rays
! nra = number of receiver arrays
! cslat,cslong = Latidude and longitude width of receiver cushion
! abw,abe,abn,abs = Limits of receiver grid (W, E, N, S)
!
OPEN(UNIT=10,FILE='fm3dt.in',STATUS='old')
READ(10,1)cdum
READ(10,1)cdum
READ(10,1)cdum
READ(10,1)sources
READ(10,1)itimes
READ(10,1)receivers
READ(10,1)gridv
READ(10,*)nrfr,nrfx,nrfz
READ(10,*)earth
READ(10,*)fom
READ(10,*)snb
READ(10,*)ltfr
READ(10,*)cslat,cslong
READ(10,1)cdum
READ(10,1)cdum
READ(10,1)cdum
READ(10,*)fsrt
READ(10,1)rtravel
READ(10,*)cfd
READ(10,1)frechet
READ(10,*)wttf,awttf
READ(10,1)travelt
READ(10,*)wrgf,awrgf
READ(10,1)wrays
1   FORMAT(a20)
CLOSE(10)
!
! Call a subroutine which reads in the velocity grid
!
CALL gridder(gridv)
!
! Compute the total number of ray paths if
! necessary
!
IF(wrgf.NE.0)THEN
   OPEN(UNIT=10,FILE=sources,STATUS='old')
   OPEN(UNIT=20,FILE=receivers,STATUS='old')
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
   tnr=0
   DO i=1,nra
      READ(10,*)nsrc
      READ(20,*)nrc
      IF(wrgf.GT.0)THEN
         IF(awrgf.EQ.i)tnr=nrc
      ELSE
         tnr=tnr+nsrc*nrc
      ENDIF
      DO j=1,nsrc
         READ(10,*)
         READ(10,*)
      ENDDO
      DO j=1,nrc
         READ(20,*)
      ENDDO
   ENDDO
   CLOSE(10)
   CLOSE(20)
ENDIF

! read valid_rays file
open(77,file='valid_rays.dat')
open(60,file=sources)
read(60,*)nra
ALLOCATE(pairs(nra))
do ii=1,nra
   read(60,*)nsrc
   ! read valid rays
   ALLOCATE(pairs(ii)%valid_rays(nsrc))
   pairs(ii)%ns = nsrc
   do i=1,nsrc
      read(77,*)pairs(ii)%valid_rays(i)
   enddo
   do i=1,nsrc
      read(60,*)
      read(60,*)
   enddo
enddo
close(77)
close(60)

!
! Read in the number of source types (i.e. the number
! of receiver arrays).
!
Open(UNIT=60,FILE=sources,STATUS='old')
READ(60,*)nra
!
! Open receiver file if required.
!
IF(fsrt.eq.1)THEN
   OPEN(UNIT=70,FILE=receivers,status='old')
   READ(70,*)idum
   IF(idum.NE.nra)THEN
      WRITE(6,*)'ERROR!!!'
      WRITE(6,*)'Source and receiver files are'
      WRITE(6,*)'inconsistent!!!!'
      WRITE(6,*)'First line of each file should'
      WRITE(6,*)'be identical!!!'
      WRITE(6,*)'TERMINATING PROGRAM!!!'
      STOP
   ENDIF
ENDIF
!
! Now work out, source by source, the first-arrival traveltime
! field plus traveltimes and ray paths from the base of the 3-D grid 
! if required. First, allocate memory to the
! traveltime field array
!
ALLOCATE(ttn(nnz,nnx,nnr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL ttn'
ENDIF
!
! Open file for source-receiver traveltime output if required.
!
IF(fsrt.eq.1)THEN
   OPEN(UNIT=10,FILE=rtravel,STATUS='unknown')
ENDIF
!
! Open file for ray path output if required
!
IF(wrgf.NE.0)THEN
   OPEN(UNIT=40,FILE=wrays,FORM='unformatted',STATUS='unknown')
   WRITE(40)tnr
ENDIF
!
! Open file for Frechet derivative output if required.
!
IF(cfd.EQ.1)THEN
   OPEN(UNIT=50,FILE=frechet,FORM='unformatted',STATUS='unknown')
ENDIF
!
! Open file containing traveltimes at base of model
!
OPEN(UNIT=20,FILE=itimes,FORM='unformatted',STATUS='old')

!
! Loop over each receiver array
!
DO ii=1,nra
   READ(60,*)nsrc

   ALLOCATE(scr(nsrc),scx(nsrc),scz(nsrc), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL scr,scx,scz'
   ENDIF
   DO i=1,nsrc
      READ(60,*)scr(i),scx(i),scz(i)
      READ(60,'(a8)')pht
!
!     Convert source coordinates in degrees to radians
!
      scx(i)=(90.0-scx(i))*pi/180.0
      scz(i)=scz(i)*pi/180.0
   ENDDO
   IF(fsrt.EQ.1.or.ltfr.EQ.1)THEN
      READ(70,*)nrc
      ALLOCATE(rcr(nrc),rcx(nrc),rcz(nrc), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL rcr,rcx,rcz'
      ENDIF
      DO i=1,nrc
         READ(70,*)rcr(i),rcx(i),rcz(i)
         IF(ltfr.EQ.1)THEN
            IF(i.EQ.1)THEN
               abe=rcz(i)
               abw=abe
               abn=rcx(i)
               abs=abn
            ELSE
               IF(rcz(i).LT.abw)THEN
                  abw=rcz(i)
               ELSE IF(rcz(i).GT.abe)THEN
                  abe=rcz(i)
               ENDIF
               IF(rcx(i).LT.abs)THEN
                  abs=rcx(i)
               ELSE IF(rcx(i).GT.abn)THEN
                  abn=rcx(i)
               ENDIF
            ENDIF
         ENDIF
!
!        Convert receiver coordinates in degrees to radians
!
         rcx(i)=(90.0-rcx(i))*pi/180.0
         rcz(i)=rcz(i)*pi/180.0
      ENDDO
      IF(ltfr.EQ.1)THEN
!
!        Add cushion values to bounds of receiver
!        array, and then locate surface grid which
!        is used to terminate the fast marching process.
!
         abw=(abw-cslong)*pi/180.0
         abe=(abe+cslong)*pi/180.0
         abs=abs-cslat
         abn=abn+cslat
         abs=(90.0-abs)*pi/180.0
         abn=(90.0-abn)*pi/180.0
         nsnn=INT((abn-gox)/dnx)
         IF(nsnn.LT.1)nsnn=1
         nsns=INT((abs-gox)/dnx)+2
         IF(nsns.GT.nnx)nsns=nnx
         nsne=INT((abe-goz)/dnz)+2
         IF(nsne.GT.nnz)nsne=nnz
         nsnw=INT((abw-goz)/dnz)
         IF(nsnw.LT.1)nsnw=1
      ENDIF
   ENDIF
   DO i=1,nsrc
      ! only tackle valid rays
      if(pairs(ii)%valid_rays(i) == 1) then
         !
         !     Call a subroutine that works out the first-arrival traveltime
         !     field.
         !
         CALL travel(i)
         !
         !     Find traveltimes from base of grid if required
         !
         IF(fsrt.eq.1)THEN
            CALL srtimes
         ENDIF
         !
         !     Calculate raypath geometries and write to file if required.
         !     Calculate Frechet derivatives with the same subroutine 
         !     if required.
         !
         IF((wrgf.EQ.i.AND.awrgf.EQ.ii).OR.wrgf.LT.0.OR.cfd.EQ.1)THEN
            CALL rpaths(wrgf,i,awrgf,ii,cfd)
         ENDIF
         !
         !     If required, write traveltime field to file
         !
         IF(wttf.EQ.i.AND.awttf.EQ.ii)THEN
            OPEN(UNIT=30,FILE=travelt,FORM='unformatted',STATUS='unknown')
            WRITE(30)gor,goxd,gozd
            WRITE(30)nnr,nnx,nnz
            WRITE(30)dnr,dnxd,dnzd
            DO j=1,nnz
               DO k=1,nnx
                  DO l=1,nnr
                     WRITE(30)ttn(j,k,l)
                  ENDDO
               ENDDO
            ENDDO
            CLOSE(30)
         ENDIF
      else
         do j=1,nrc
            write(10,*)0.0
         enddo
         k = 0
         do j=1,nrc
            WRITE(50)k
         enddo
      endif
   ENDDO
   DEALLOCATE(scr,scx,scz)
   IF(fsrt.EQ.1.or.ltfr.EQ.1)THEN
      DEALLOCATE(rcr,rcx,rcz)
   ENDIF
ENDDO
CLOSE(60)
!
! Close rtravel if required
!
IF(fsrt.eq.1)THEN
   CLOSE(70)
ENDIF
CLOSE(20)
IF(wrgf.NE.0)THEN
   CLOSE(40)
ENDIF
IF(cfd.EQ.1)THEN
   CLOSE(50)
ENDIF
DEALLOCATE (veln,ttn, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin3d: final deallocate'
ENDIF
WRITE(6,*)'Program fm3dt has finished successfully!'
do i=1,nra
   call FmttPair_dealloc(pairs(i))
enddo
DEALLOCATE(pairs)
END PROGRAM fmmin3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed the name of the velocity
! grid file (grid) and reads in the velocity values.
! The gridded values are globally shared via
! a MODULE statement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gridder(gridv)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,l,m,n,i1,j1,k1
INTEGER :: conr,conx,conz,str,stx,stz
REAL(KIND=i10) :: u,sumi,sumj,sumk
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: velv
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi,wi
CHARACTER (LEN=20) :: gridv
!
! velv = B-spline velocity grid values
! u = Cubic spline independent variable
! ui,vi,wi = Cubic spline basis functions
! sumi,sumj,sumk = Summation variables for constructing spline
! conr,conx,conz = Counters for refining grid in r,theta,phi
! str,stx,stz = Refined grid location in r,theta,phi
! gridv = Input B-spline vertex file
!
! Read in B-spline grid
!
OPEN(UNIT=10,FILE=gridv,STATUS='old')
READ(10,*)nvr,nvx,nvz
READ(10,*)gor,goxd,gozd
READ(10,*)dvr,dvx,dvz
ALLOCATE(velv(0:nvr+1,0:nvx+1,0:nvz+1))
DO i=0,nvz+1
   DO j=0,nvx+1
      DO k=0,nvr+1
         READ(10,*)velv(k,j,i)
      ENDDO
   ENDDO
ENDDO
CLOSE(10)
!
! Calculate total numer of refined nodes in r,theta,phi
! and the refined grid spacing.
!
nnr=(nvr-1)*nrfr+1
nnx=(nvx-1)*nrfx+1
nnz=(nvz-1)*nrfz+1
dnr=dvr/nrfr
dnxd=dvx/nrfx
dnzd=dvz/nrfz
ALLOCATE(veln(nnz,nnx,nnr))
!
! Calculate the values of the basis functions
!
ALLOCATE(ui(nvr+1,4))
DO i=1,nrfr+1
   u=nrfr
   u=(i-1)/u
   ui(i,1)=(1.0-u)**3/6.0
   ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   ui(i,4)=u**3/6.0
ENDDO
ALLOCATE(vi(nvx+1,4))
DO i=1,nrfx+1
   u=nrfx
   u=(i-1)/u
   vi(i,1)=(1.0-u)**3/6.0
   vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   vi(i,4)=u**3/6.0
ENDDO
ALLOCATE(wi(nvz+1,4))
DO i=1,nrfz+1
   u=nrfz
   u=(i-1)/u
   wi(i,1)=(1.0-u)**3/6.0
   wi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   wi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   wi(i,4)=u**3/6.0
ENDDO
!
! Calculate velocity values on refined grid
!
DO i=1,nvz-1
   conz=nrfz
   IF(i==nvz-1)conz=nrfz+1
   DO j=1,nvx-1
      conx=nrfx
      IF(j==nvx-1)conx=nrfx+1
      DO k=1,nvr-1
         conr=nrfr
         IF(k==nvr-1)conr=nrfr+1
         DO l=1,conz
            stz=nrfz*(i-1)+l
            DO m=1,conx
               stx=nrfx*(j-1)+m
               DO n=1,conr
                  str=nrfr*(k-1)+n
                  sumi=0.0
                  DO i1=1,4
                     sumj=0.0
                     DO j1=1,4
                        sumk=0.0
                        DO k1=1,4
                           sumk=sumk+ui(n,k1)*velv(k-2+k1,j-2+j1,i-2+i1)
                        ENDDO
                        sumj=sumj+vi(m,j1)*sumk
                     ENDDO
                     sumi=sumi+wi(l,i1)*sumj
                  ENDDO
                  veln(stz,stx,str)=sumi
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
! Finally, convert all parameters in degrees to radians
!
dnx=dnxd*pi/180.0
dnz=dnzd*pi/180.0
gox=(90.0-goxd)*pi/180.0
goz=gozd*pi/180.0
dvx=dvx*pi/180.0
dvz=dvz*pi/180.0
DEALLOCATE(velv)
DEALLOCATE(ui,vi,wi)
END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates all receiver traveltimes for
! a given source and writes the results to file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE srtimes
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,l,irr,irx,irz,sw
REAL :: trr
REAL(KIND=i10) :: drr,drx,drz,produ
!
! irr,irx,irz = Coordinates of cell containing receiver
! trr = traveltime value at receiver
! produ = dummy multiplier
! drr,drx,drz = receiver distance from (i,j,k) grid node
!
! Determine source-receiver traveltimes one at a time.
!
DO i=1,nrc
! 
!  The first step is to locate the receiver in the grid.
!
   irr=INT((gor-rcr(i))/dnr)+1
   irx=INT((rcx(i)-gox)/dnx)+1
   irz=INT((rcz(i)-goz)/dnz)+1
   sw=0
   IF(irr.lt.1.or.irr.gt.nnr)sw=1
   IF(irx.lt.1.or.irx.gt.nnx)sw=1
   IF(irz.lt.1.or.irz.gt.nnz)sw=1
   IF(sw.eq.1)then
      WRITE(6,*)"Receiver lies outside model (rr,xr,zr)= ",irr,irx,irz
   ENDIF
   IF(irr.eq.nnr)irr=irr-1
   IF(irx.eq.nnx)irx=irx-1
   IF(irz.eq.nnz)irz=irz-1
!
!  Location of receiver successfully found within the grid. Now approximate
!  traveltime at receiver using trilinear interpolation from eight
!  surrounding grid points
!
   drx=(rcx(i)-gox)-(irx-1)*dnx
   drz=(rcz(i)-goz)-(irz-1)*dnz
   drr=(gor-rcr(i))-(irr-1)*dnr
   trr=0.0
   DO j=1,2
      DO k=1,2
         DO l=1,2
            produ=(1.0-ABS(((l-1)*dnz-drz)/dnz))*(1.0-ABS(((k-1)*dnx-drx)/dnx))
            produ=produ*(1.0-ABS(((j-1)*dnr-drr)/dnr))
            trr=trr+ttn(irz-1+l,irx-1+k,irr-1+j)*produ
         ENDDO
      ENDDO
   ENDDO
   WRITE(10,*)trr
ENDDO
END SUBROUTINE srtimes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates ray path geometries for each
! receiver from the base of the model. It will also compute
! Frechet derivatives using these ray paths if required.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE rpaths(wrgf,csid,awrgf,caid,cfd)
USE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
INTEGER :: i,j,k,l,m,n,ipr,ipx,ipz,nrp,sw
INTEGER :: wrgf,cfd,csid,ipro,ipxo,ipzo,caid,awrgf
INTEGER :: ivr,ivx,ivz,ivro,ivxo,ivzo,nhp
INTEGER :: ivrt,ivxt,ivzt,iprt,ipxt,ipzt,isum
INTEGER, DIMENSION (4) :: chp
INTEGER, PARAMETER :: maxrp=100000
REAL(KIND=i5), PARAMETER :: ftol=1.0e-6
REAL(KIND=i5) :: rayr,rayx,rayz
REAL(KIND=i10) :: dpl,crad,rd1,rd2,atio,ri,xi,zi,vel,velo
REAL(KIND=i10) :: u,v,w,rigz,rigx,rigr,dinc
REAL(KIND=i10) :: dtr,dtx,dtz,drr,drx,drz,produ
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: rgr,rgx,rgz
REAL(KIND=i5), DIMENSION (:,:,:), ALLOCATABLE :: fdm
REAL(KIND=i10), DIMENSION (4) :: vrat,ui,vi,wi,uio,vio,wio
!
! ipr,ipx,ipz = Coordinates of cell containing current point
! ipro,ipxo,ipzo = Coordinates of previous point
! rgr,rgx,rgz = (r,x,z) coordinates of ray geometry
! ivr,ivx,ivz = Coordinates of B-spline vertex containing current point
! ivro,ivxo,ivzo = Coordinates of previous point
! maxrp = maximum number of ray points
! nrp = number of points to describe ray
! dpl = incremental path length of ray
! crad = current radius of ray point
! atio = ratio for locating ray endpoint
! ri,xi,zi = edge of model coordinates
! dtr,dtx,dtz = components of gradT
! wrgf = Write out raypaths? (<0=all,0=no,>0=souce id)
! awrgf = Array ID if wrgf>0
! cfd = calculate Frechet derivatives? (0=no,1=yes)
! csid = current source id
! fdm = Frechet derivative matrix
! nhp = Number of ray segment-B-spline cell hit points
! vrat = length ratio of ray sub-segment
! chp = pointer to incremental change in r,x or z cell
! drr,drx,drz = distance from reference node of cell
! produ = variable for trilinear interpolation
! vel = velocity at current point
! velo = velocity at previous point
! u,v,w = local variables of r,x,z
! ui,vi,wi = B-spline basis functions at current point
! uio,vio,wio = ui,vi,wi for previous point
! ivrt,ivxt,ivzt = temporary ivr,ivx,ivz values
! rigr,rigx,rigz = end point of sub-segment of ray path
! iprt,ipxt,ipzt = temporary ipr,ipx,ipz values
! dinc = path length of ray sub-segment
! rayr,rayx,rayz = ray path coordinates in single precision
! caid = current receiver array id.
!
! Allocate memory to arrays for storing ray path geometry
!
ALLOCATE(rgr(nnx*nnr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgr'
ENDIF
ALLOCATE(rgx(nnx*nnr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgx'
ENDIF
ALLOCATE(rgz(nnx*nnr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgz'
ENDIF
!
! Allocate memory to partial derivative array
!
IF(cfd.EQ.1)THEN
   ALLOCATE(fdm(0:nvz+1,0:nvx+1,0:nvr+1), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE : REAL veln'
   ENDIF
ENDIF
!
! Set ray incremental path length equal to grid separation
! in depth.
!
  dpl=dnr
!
! Loop through all the receivers
!
DO i=1,nrc
   IF(cfd.EQ.1)THEN
      fdm=0.0
   ENDIF
!
!  The first step is to locate the receiver in the grid.
!
   ipr=INT((gor-rcr(i))/dnr)+1
   ipx=INT((rcx(i)-gox)/dnx)+1
   ipz=INT((rcz(i)-goz)/dnz)+1
   sw=0
   IF(ipr.lt.1.or.ipr.ge.nnr)sw=1
   IF(ipx.lt.1.or.ipx.ge.nnx)sw=1
   IF(ipz.lt.1.or.ipz.ge.nnz)sw=1
   IF(sw.eq.1)then
      WRITE(6,*)"Receiver lies outside model (rr,xr,zr)= ",ipr,ipx,ipz
   ENDIF
   IF(ipr.eq.nnr)ipr=ipr-1
   IF(ipx.eq.nnx)ipx=ipx-1
   IF(ipz.eq.nnz)ipz=ipz-1
!  
!  First point of the ray path is the receiver
!
   rgr(1)=rcr(i)
   rgx(1)=rcx(i)
   rgz(1)=rcz(i)
!
!  Now trace ray from receiver to "source"
!
   DO j=1,maxrp
!
!     Calculate traveltime gradient vector for current cell
!
      dtr=ttn(ipz,ipx,ipr+1)-ttn(ipz,ipx,ipr)
      dtr=dtr+ttn(ipz+1,ipx,ipr+1)-ttn(ipz+1,ipx,ipr)
      dtr=dtr+ttn(ipz+1,ipx+1,ipr+1)-ttn(ipz+1,ipx+1,ipr)
      dtr=dtr+ttn(ipz,ipx+1,ipr+1)-ttn(ipz,ipx+1,ipr)
      dtr=dtr/(4.0*dnr)
      dtx=ttn(ipz,ipx+1,ipr)-ttn(ipz,ipx,ipr)
      dtx=dtx+ttn(ipz+1,ipx+1,ipr)-ttn(ipz+1,ipx,ipr)
      dtx=dtx+ttn(ipz+1,ipx+1,ipr+1)-ttn(ipz+1,ipx,ipr+1)
      dtx=dtx+ttn(ipz,ipx+1,ipr+1)-ttn(ipz,ipx,ipr+1)
      dtx=dtx/(4.0*(earth+gor-(ipr-1)*dnr)*dnx)
      dtz=ttn(ipz+1,ipx,ipr)-ttn(ipz,ipx,ipr)
      dtz=dtz+ttn(ipz+1,ipx+1,ipr+1)-ttn(ipz,ipx+1,ipr+1)
      dtz=dtz+ttn(ipz+1,ipx+1,ipr)-ttn(ipz,ipx+1,ipr)
      dtz=dtz+ttn(ipz+1,ipx,ipr+1)-ttn(ipz,ipx,ipr+1)
      dtz=dtz/(4.0*(earth+gor-(ipr-1)*dnr)*SIN(rgx(j))*dnz)
!
!     Calculate the next ray path point
!
      crad=earth+rgr(j)
      rd1=SQRT(dtr**2+dtx**2+dtz**2)
      rgr(j+1)=rgr(j)+dpl*dtr/rd1
      rgx(j+1)=rgx(j)-dpl*dtx/(crad*rd1)
      rgz(j+1)=rgz(j)-dpl*dtz/(crad*SIN(rgx(j))*rd1)
!
!     Determine which cell the new ray endpoint
!     lies in.
!
      ipro=ipr
      ipxo=ipx
      ipzo=ipz
      ipr=INT((gor-rgr(j+1))/dnr)+1
      ipx=INT((rgx(j+1)-gox)/dnx)+1
      ipz=INT((rgz(j+1)-goz)/dnz)+1
      sw=0
      IF(ipr.LT.1.OR.ipr.GE.nnr)sw=1
      IF(ipx.LT.1.OR.ipx.GE.nnx)sw=1
      IF(ipz.LT.1.OR.ipz.GE.nnz)sw=1
!
!     If sw.ne.0 then we are done; find
!     the intersection point of the ray segment
!     with the boundary of the model.
!
      IF(sw.NE.0)THEN
         sw=0
         IF(ipr.LT.1.OR.ipr.GE.nnr)THEN
            IF(ipr.LT.1)THEN
               ri=gor
               ipr=1
            ELSE
               ri=gor-(nnr-1)*dnr
               ipr=nnr-1
            ENDIF
            atio=(ri-rgr(j))/(rgr(j+1)-rgr(j))
            sw=1
         ENDIF
         IF(ipx.LT.1.OR.ipx.GE.nnx)THEN
            IF(ipx.LT.1)THEN
               xi=gox
               ipx=1
            ELSE
               xi=(gox+(nnx-1)*dnx)
               ipx=nnx-1
            ENDIF
            rd1=(xi-rgx(j))/(rgx(j+1)-rgx(j))
            IF(sw.eq.1)then
               IF(rd1.LT.atio)THEN
                  atio=rd1
                  sw=2
               ENDIF
            ELSE
               atio=rd1
               sw=2
            ENDIF
         ENDIF
         IF(ipz.LT.1.OR.ipz.GE.nnz)THEN
            IF(ipz.LT.1)THEN
               zi=goz
               ipz=1
            ELSE
               zi=goz+(nnz-1)*dnz
               ipz=nnz-1
            ENDIF
            rd1=(zi-rgz(j))/(rgz(j+1)-rgz(j))
            IF(sw.NE.0)then
               IF(rd1.LT.atio)THEN
                  atio=rd1
                  sw=3
               ENDIF
            ELSE
               atio=rd1
               sw=3
            ENDIF
         ENDIF
         IF(sw.EQ.1)THEN
            ipx=ipxo
            ipz=ipzo
         ELSE IF(sw.EQ.2)THEN
            ipr=ipro
            ipz=ipzo
         ELSE IF(sw.EQ.3)THEN
            ipr=ipro
            ipx=ipxo
         ENDIF
         rgz(j+1)=rgz(j)+atio*(rgz(j+1)-rgz(j))
         rgx(j+1)=rgx(j)+atio*(rgx(j+1)-rgx(j))
         rgr(j+1)=rgr(j)+atio*(rgr(j+1)-rgr(j))
         nrp=j+1
         IF(cfd.EQ.1)THEN
            sw=1
         ELSE
            EXIT
         ENDIF
      ENDIF
!
!     Calculate the Frechet derivatives if required.
!
      IF(cfd.EQ.1)THEN
!
!        First determine which B-spline cell the refined cells
!        containing the ray path segment lies in. If they lie
!        in more than one, then we need to divide the problem
!        into separate parts (up to four).
!
         ivr=INT((ipr-1)/nrfr)+1
         ivx=INT((ipx-1)/nrfx)+1
         ivz=INT((ipz-1)/nrfz)+1
         ivro=INT((ipro-1)/nrfr)+1
         ivxo=INT((ipxo-1)/nrfx)+1
         ivzo=INT((ipzo-1)/nrfz)+1
!
!        Calculate up to three hit points between straight
!        ray segment and cell faces.
!
         nhp=0
         IF(ivr.NE.ivro)THEN
            nhp=nhp+1
            IF(ivr.GT.ivro)THEN
               ri=gor-(ivr-1)*dvr
            ELSE
               ri=gor-ivr*dvr
            ENDIF
            vrat(nhp)=(ri-rgr(j))/(rgr(j+1)-rgr(j))
            chp(nhp)=1
         ENDIF
         IF(ivx.NE.ivxo)THEN
            nhp=nhp+1
            IF(ivx.GT.ivxo)THEN
               xi=gox+(ivx-1)*dvx
            ELSE
               xi=gox+ivx*dvx
            ENDIF
            rd1=(xi-rgx(j))/(rgx(j+1)-rgx(j))
            IF(nhp.EQ.1)THEN
               vrat(nhp)=rd1
               chp(nhp)=2
            ELSE
               IF(rd1.GE.vrat(nhp-1))THEN
                  vrat(nhp)=rd1
                  chp(nhp)=2
               ELSE
                  vrat(nhp)=vrat(nhp-1)
                  chp(nhp)=chp(nhp-1)
                  vrat(nhp-1)=rd1
                  chp(nhp-1)=2
               ENDIF
            ENDIF
         ENDIF
         IF(ivz.NE.ivzo)THEN
            nhp=nhp+1
            IF(ivz.GT.ivzo)THEN
               zi=goz+(ivz-1)*dvz 
            ELSE
               zi=goz+ivz*dvz
            ENDIF
            rd1=(zi-rgz(j))/(rgz(j+1)-rgz(j))
            IF(nhp.EQ.1)THEN
               vrat(nhp)=rd1
               chp(nhp)=3
            ELSE IF(nhp.EQ.2)THEN
               IF(rd1.GE.vrat(nhp-1))THEN
                  vrat(nhp)=rd1
                  chp(nhp)=3
               ELSE
                  vrat(nhp)=vrat(nhp-1)
                  chp(nhp)=chp(nhp-1)
                  vrat(nhp-1)=rd1
                  chp(nhp-1)=3
               ENDIF
            ELSE
               IF(rd1.GE.vrat(nhp-1))THEN
                  vrat(nhp)=rd1
                  chp(nhp)=3
               ELSE IF(rd1.GE.vrat(nhp-2))THEN
                  vrat(nhp)=vrat(nhp-1)
                  chp(nhp)=chp(nhp-1)
                  vrat(nhp-1)=rd1
                  chp(nhp-1)=3
               ELSE
                  vrat(nhp)=vrat(nhp-1)
                  chp(nhp)=chp(nhp-1)
                  vrat(nhp-1)=vrat(nhp-2)
                  chp(nhp-1)=chp(nhp-2)
                  vrat(nhp-2)=rd1
                  chp(nhp-2)=3
               ENDIF
            ENDIF
         ENDIF
         nhp=nhp+1
         vrat(nhp)=1.0
         chp(nhp)=0
!
!        Calculate the velocity, u,v and w values of the
!        first point
!
         drr=(gor-rgr(j))-(ipro-1)*dnr
         drx=(rgx(j)-gox)-(ipxo-1)*dnx
         drz=(rgz(j)-goz)-(ipzo-1)*dnz
         vel=0.0
         DO k=1,2
            DO l=1,2
               DO m=1,2
                  produ=(1.0-ABS(((m-1)*dnz-drz)/dnz))
                  produ=produ*(1.0-ABS(((l-1)*dnx-drx)/dnx))
                  produ=produ*(1.0-ABS(((k-1)*dnr-drr)/dnr))
         IF(ipzo-1+m.LE.nnz.AND.ipxo-1+l.LE.nnx.AND.ipro-1+k.LE.nnr)THEN
                  vel=vel+veln(ipzo-1+m,ipxo-1+l,ipro-1+k)*produ
         ENDIF
               ENDDO
            ENDDO
         ENDDO
         drr=(gor-rgr(j))-(ivro-1)*dvr
         drx=(rgx(j)-gox)-(ivxo-1)*dvx
         drz=(rgz(j)-goz)-(ivzo-1)*dvz
         u=drr/dvr
         v=drx/dvx
         w=drz/dvz
!
!        Calculate the 12 basis values at the point
!
         ui(1)=(1.0-u)**3/6.0
         ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
         ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
         ui(4)=u**3/6.0
         vi(1)=(1.0-v)**3/6.0
         vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
         vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
         vi(4)=v**3/6.0
         wi(1)=(1.0-w)**3/6.0
         wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
         wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
         wi(4)=w**3/6.0
         ivrt=ivro
         ivxt=ivxo
         ivzt=ivzo
!
!        Now loop through the one or more sub-segments of the
!        ray path segment and calculate partial derivatives
!
         DO k=1,nhp
            velo=vel
            uio=ui
            vio=vi
            wio=wi
            IF(k.GT.1)THEN
               IF(chp(k-1).EQ.1)THEN
                  ivrt=ivr
               ELSE IF(chp(k-1).EQ.2)THEN
                  ivxt=ivx
               ELSE IF(chp(k-1).EQ.3)THEN
                  ivzt=ivz
               ENDIF
            ENDIF
!
!           Calculate the velocity, u,v and w values of the
!           new point
!
            rigz=rgz(j)+vrat(k)*(rgz(j+1)-rgz(j))
            rigx=rgx(j)+vrat(k)*(rgx(j+1)-rgx(j))
            rigr=rgr(j)+vrat(k)*(rgr(j+1)-rgr(j))
            iprt=INT((gor-rigr)/dnr)+1
            ipxt=INT((rigx-gox)/dnx)+1
            ipzt=INT((rigz-goz)/dnz)+1
            drr=(gor-rigr)-(iprt-1)*dnr
            drx=(rigx-gox)-(ipxt-1)*dnx
            drz=(rigz-goz)-(ipzt-1)*dnz
            vel=0.0
            DO l=1,2
               DO m=1,2
                  DO n=1,2
                     produ=(1.0-ABS(((n-1)*dnz-drz)/dnz))
                     produ=produ*(1.0-ABS(((m-1)*dnx-drx)/dnx))
                     produ=produ*(1.0-ABS(((l-1)*dnr-drr)/dnr))
            IF(ipzt-1+n.LE.nnz.AND.ipxt-1+m.LE.nnx.AND.iprt-1+l.LE.nnr)THEN
                     vel=vel+veln(ipzt-1+n,ipxt-1+m,iprt-1+l)*produ
            ENDIF
                  ENDDO
               ENDDO
            ENDDO
            drr=(gor-rigr)-(ivrt-1)*dvr
            drx=(rigx-gox)-(ivxt-1)*dvx
            drz=(rigz-goz)-(ivzt-1)*dvz
            u=drr/dvr
            v=drx/dvx
            w=drz/dvz
!
!           Calculate the 12 basis values at the new point
!
            ui(1)=(1.0-u)**3/6.0
            ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
            ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
            ui(4)=u**3/6.0
            vi(1)=(1.0-v)**3/6.0
            vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
            vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
            vi(4)=v**3/6.0
            wi(1)=(1.0-w)**3/6.0
            wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
            wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
            wi(4)=w**3/6.0
!
!           Calculate the incremental path length
!     
            IF(k.EQ.1)THEN
               dinc=vrat(k)*dpl
            ELSE 
               dinc=(vrat(k)-vrat(k-1))*dpl
            ENDIF
!
!           Now compute the 64 contributions to the partial
!           derivatives.
!
            DO l=1,4
               DO m=1,4
                  DO n=1,4
                     rd1=ui(n)*vi(m)*wi(l)/vel**2
                     rd2=uio(n)*vio(m)*wio(l)/velo**2
                     rd1=-(rd1+rd2)*dinc/2.0
                     rd2=fdm(ivzt-2+l,ivxt-2+m,ivrt-2+n)
                     fdm(ivzt-2+l,ivxt-2+m,ivrt-2+n)=rd1+rd2
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF(sw.EQ.1)EXIT
      ENDIF
   ENDDO
!
!  Write ray paths to output file
!
   IF((wrgf.EQ.csid.AND.awrgf.EQ.caid).OR.wrgf.LT.0)THEN
      WRITE(40)nrp
      DO j=1,nrp
         rayr=rgr(j)
         rayx=(pi/2-rgx(j))*180.0/pi
         rayz=rgz(j)*180.0/pi
         WRITE(40)rayr,rayx,rayz
      ENDDO
   ENDIF
!
!  Write partial derivatives to output file
!
   IF(cfd.EQ.1)THEN
!
!     Determine the number of non-zero elements.
!
      isum=0
      DO j=0,nvz+1
         DO k=0,nvx+1
            DO l=0,nvr+1
               IF(ABS(fdm(j,k,l)).GE.ftol)isum=isum+1
            ENDDO
         ENDDO
      ENDDO
      WRITE(50)isum
      isum=0
      DO j=0,nvz+1
         DO k=0,nvx+1
            DO l=0,nvr+1
               isum=isum+1
               IF(ABS(fdm(j,k,l)).GE.ftol)WRITE(50)isum,fdm(j,k,l)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDDO
IF(cfd.EQ.1)THEN
   DEALLOCATE(fdm, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE rpaths: fdm'
   ENDIF
ENDIF
DEALLOCATE(rgr,rgx,rgz, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE rpaths: rgr,rgx,rgz'
ENDIF
END SUBROUTINE rpaths
