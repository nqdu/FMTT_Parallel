!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to take depth slices through
! the FMM computed traveltime field and put them in a
! form suitable for contour plotting by GMT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM slice
IMPLICIT NONE
INTEGER :: i,j,k,l,m,n,i1,j1,k1,extrds,extrnss,extrews,idp,iew,ins
INTEGER :: checkstat,isw,conp,cont,conr,id1,id2
INTEGER :: nnr,nnt,nnp,pvgd,pvgns,pvgew,prpd,prpns,prpew,nre,nr
INTEGER :: idm1,idm2,idm3,nnx,nnz,stp,stt,str
INTEGER :: ptfd,ptfns,ptfew,rnode
INTEGER :: ddt,ddp,dnsr,dnst,dewr,dewp
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
REAL(KIND=i5) :: rdep,rlat,rlon
REAL(KIND=i10) :: gor,got,gop,rgsr,rgst,rgsp,sldep,slns,slew
REAL(KIND=i10) :: tt,ttt,ttb,sumi,sumj,sumk
REAL(KIND=i10) :: lft,rgt,btm,top
REAL(KIND=i10) :: rdm,rd1,rd2,rd3,u
REAL(KIND=i10), DIMENSION(:) :: wi(4)
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi,vela
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: ttn,veln
CHARACTER (LEN=30) :: ifilet,ifilev,irfile,ifilevr,sep
CHARACTER (LEN=30) :: ofiledb,ofilensb,ofileewb
CHARACTER (LEN=30) :: ofiledv,ofilensv,ofileewv
CHARACTER (LEN=30) :: ofiledt,ofilenst,ofileewt
CHARACTER (LEN=30) :: ofiledr,ofilensr,ofileewr
!
! sldep = slice depth
! slns = NS slice
! slew = EW slice
! ifilet = input 3-D traveltime grid file
! ifilev = input velocity grid file
! ifilevr = input reference velocity grid file
! irfile = input ray path file
! ofiledt = output depth slice file for traveltimes
! ofilenst = output N-S slice file for traveltimes
! ofileewt = output E-W slice file for traveltimes
! ofiledb = bounds for output depth slice file
! ofilensb = bounds for output N-S slice file
! ofileewb = bounds for output E-W slice file
! ofiledv = output velocity file for depth slice
! ofilensv = output velocity file for N-S slice
! ofileewv = output velocity file for E-W slice
! ofiledr = output ray path file in depth
! ofilensr = output ray path file in N-S
! ofileewr = output ray path file in E-W
! nnr,nnt,nnp = number of diced nodes in r,theta,phi
! gor = grid origin in radius
! got = grid origin in theta (N-S)
! gop = grid origin in phi (E-W)
! rgsr,rgst,rgsp = Refined node spacing in r,theta,phi
! idp = radial index of slice depth
! ins = NS index for latitude slice
! iew = EW index for longitude slice
! tt = traveltime at or near a node
! ttn = traveltime grid values
! ttt = traveltime of node above slice 
! ttb =  traveltime of node below slice
! lft,rgt,btm,top = plotting bounds
! pvgd = plot velocity depth slice (0=no,1=yes)
! pvgns = plot velocity N-S slice (0=no, 1=yes)
! pvgew = plot velocity E-W slice (0=no, 1=yes)
! ptfd = plot traveltime field depth slice? (0=no,1=yes)
! ptfns = plot traveltime field N-S slice? (0=no, 1=yes)
! ptfew = plot traveltime field E-W slice? (0=no,1=yes)
! veln = velocity grid values
! prpd = plot ray path in depth? (0=no,1=yes)
! prpns = plot ray path in N-S? (0=no,1=yes)
! prpew = plot ray path in E-W? (0=no,1=yes)
! nre = number of ray elements
! rdep,rlat,rlon = ray point depth, latitude, longitude
! nr = number of receivers
! sep = character marker for separating rays
! extrds = extract depth slice? (0=no,1=yes)
! extrnss = extract N-S slice? (0=no,1=yes)
! extrews = extract E-W slice? (0=no,1=yes)
! ddt,ddp = dicing of velocity grid for depth slice
! dnsr,dnst = dicing of velocity grid for N-S slice
! dewr,dewp = dicing of velocity grid for E-W slice
! u = bspline independent coordinate
! ui,vi = bspline basis functions
! vela = diced velocity values
! nnx,nnz = dimensions of vela
! conr,conp,cont = variables for edge of bspline grid
! str,stp,stt = counters for vela grid points
! rnode = reference node for slice
!
OPEN(UNIT=10,FILE='gmtslicet.in',STATUS='old')
!
! Read in input file names
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a26)')ifilev
READ(10,'(a26)')ifilevr
READ(10,'(a26)')ifilet
READ(10,'(a26)')irfile
!
! Read in slice parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)extrds
READ(10,*)sldep
READ(10,'(a22)')ofiledb
READ(10,*)extrnss
READ(10,*)slns
READ(10,'(a22)')ofilensb
READ(10,*)extrews
READ(10,*)slew
READ(10,'(a22)')ofileewb
!
! Calculate GMT bounds files. Start off by reading in
! velocity grid.
!
OPEN(UNIT=20,FILE=ifilev,status='old')
READ(20,*)nnr,nnt,nnp
READ(20,*)gor,got,gop
READ(20,*)rgsr,rgst,rgsp
ALLOCATE(veln(0:nnr+1,0:nnt+1,0:nnp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL veln'
ENDIF
DO i=0,nnp+1
   DO j=0,nnt+1
      DO k=0,nnr+1
         READ(20,*)veln(k,j,i)
      ENDDO
   ENDDO
ENDDO
CLOSE(20)
!
! Now read in velocity grid parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)pvgd
READ(10,*)ddt,ddp
READ(10,'(a22)')ofiledv
READ(10,*)pvgns
READ(10,*)dnsr,dnst
READ(10,'(a22)')ofilensv
READ(10,*)pvgew
READ(10,*)dewr,dewp
READ(10,'(a22)')ofileewv
!
! Calculate GMT bounds file for depth slice if required.
!
IF(extrds.EQ.1)THEN
   lft=gop
   rgt=gop+(nnp-1)*rgsp
   btm=got-(nnt-1)*rgst
   top=got
   OPEN(UNIT=50,FILE=ofiledb,STATUS='unknown')
   WRITE(50,'(f16.11)')lft
   WRITE(50,'(f16.11)')rgt
   WRITE(50,'(f16.11)')btm
   WRITE(50,'(f16.11)')top
   id1=(nnp-1)*ddp+1
   id2=(nnt-1)*ddt+1
   WRITE(50,*)id1
   WRITE(50,*)id2
ENDIF
!
! Calculate GMT bounds file for N-S slice if required.
!
IF(extrnss.EQ.1)THEN
   rgt=got-(nnt-1)*rgst
   lft=got
   btm=gor-(nnr-1)*rgsr
   top=gor
   OPEN(UNIT=60,FILE=ofilensb,STATUS='unknown')
   WRITE(60,'(f16.11)')rgt
   WRITE(60,'(f16.11)')lft
   WRITE(60,'(f16.11)')btm
   WRITE(60,'(f16.11)')top
   id1=(nnt-1)*dnst+1
   id2=(nnr-1)*dnsr+1
   WRITE(60,*)id1
   WRITE(60,*)id2
ENDIF
!
! Calculate GMT bounds file for E-W slice if required.
!
IF(extrews.EQ.1)THEN
   lft=gop
   rgt=gop+(nnp-1)*rgsp
   btm=gor-(nnr-1)*rgsr
   top=gor
   OPEN(UNIT=70,FILE=ofileewb,STATUS='unknown')
   WRITE(70,'(f16.11)')lft
   WRITE(70,'(f16.11)')rgt
   WRITE(70,'(f16.11)')btm
   WRITE(70,'(f16.11)')top
   id1=(nnp-1)*dewp+1
   id2=(nnr-1)*dewr+1
   WRITE(70,*)id1
   WRITE(70,*)id2
ENDIF
!
! Read in reference velocity grid file if required
!
IF(pvgd.EQ.1.OR.pvgns.EQ.1.OR.pvgew.EQ.1)THEN
   isw=0
   OPEN(UNIT=20,FILE=ifilevr,status='old')
   READ(20,*)idm1,idm2,idm3
   IF(idm1.NE.nnr.OR.idm2.NE.nnt.OR.idm3.NE.nnp)isw=1
   READ(20,*)rd1,rd2,rd3
   IF(rd1.NE.gor.OR.rd2.NE.got.OR.rd3.NE.gop)isw=1
   READ(20,*)rd1,rd2,rd3
   !if((abs(rd1 - rgsr) > 1.0e-4) .or.(abs(rd2 - rgst) > 1.0e-4) &
   !   .or. (abs(rd3 - rgsp) > 1.0e-4) ) then
   !      isw = 1
   !endif
   !IF(rd1.NE.rgsr.OR.rd2.NE.rgst.OR.rd3.NE.rgsp)isw=1
   IF(isw.EQ.1)THEN
      WRITE(6,*)'ERROR! Actual velocity grid and reference'
      WRITE(6,*)'velocity grid have different dimensions or'
      WRITE(6,*)'different numbers of grid points!'
      WRITE(6,*)'TERMINATING PROGRAM!!!'
      STOP
   ENDIF
   DO i=0,nnp+1
      DO j=0,nnt+1
         DO k=0,nnr+1
            READ(20,*)rd1
!
!           This gives the velocity anomaly.
!
            veln(k,j,i)=100 * (veln(k,j,i)-rd1) / rd1
         ENDDO
      ENDDO
   ENDDO
   CLOSE(20)
ENDIF
!
! Extract depth slice if required
!
IF(pvgd.EQ.1)THEN
!
! Allocate memory to velocity grid array
!
  nnx=(nnt-1)*ddt+1
  nnz=(nnp-1)*ddp+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Compute the values of the basis functions
!
  ALLOCATE(ui(ddt+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,ddt+1
     u=ddt
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(ddp+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,ddp+1
     u=ddp
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO
!
! Work out the value of u for the given depth
!
  rnode=INT((gor-sldep)/rgsr)+1
  u=ABS(gor-sldep-(rnode-1)*rgsr)/rgsr
  wi(1)=(1.0-u)**3/6.0
  wi(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
  wi(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
  wi(4)=u**3/6.0
  DO i=1,nnp-1
     conp=ddp
     IF(i==nnp-1)conp=ddp+1
     DO j=1,nnt-1
        cont=ddt
        IF(j==nnt-1)cont=ddt+1
        DO l=1,conp
           stp=ddp*(i-1)+l
           DO m=1,cont
              stt=ddt*(j-1)+m
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumk=0.0
                    DO k1=1,4
                       rdm=wi(k1)*veln(rnode-2+k1,j-2+j1,i-2+i1)
                       sumk=sumk+rdm
                    ENDDO
                    sumj=sumj+ui(m,j1)*sumk
                 ENDDO
                 sumi=sumi+vi(l,i1)*sumj
              ENDDO
              vela(stt,stp)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofiledv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF
!
! Extract N-S slice if required
!
IF(pvgns.EQ.1)THEN
!
! Allocate memory to velocity grid array
!
  nnx=(nnr-1)*dnsr+1
  nnz=(nnt-1)*dnst+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Compute the values of the basis functions
!
  ALLOCATE(ui(dnsr+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,dnsr+1
     u=dnsr
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(dnst+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,dnst+1
     u=dnst
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO
!
! Work out the value of u for the given longitude
!
  rnode=INT((slns-gop)/rgsp)+1
  u=ABS(slns-gop-(rnode-1)*rgsp)/rgsp
  wi(1)=(1.0-u)**3/6.0
  wi(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
  wi(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
  wi(4)=u**3/6.0
  DO j=1,nnt-1
     cont=dnst
     IF(j==nnt-1)cont=dnst+1
     DO k=1,nnr-1
        conr=dnsr
        IF(k==nnr-1)conr=dnsr+1
        DO m=1,cont
           stt=dnst*(j-1)+m
           DO n=1,conr
              str=dnsr*(k-1)+n
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumk=0.0
                    DO k1=1,4
                       rdm=ui(n,k1)*veln(k-2+k1,j-2+j1,rnode-2+i1)
                       sumk=sumk+rdm
                    ENDDO
                    sumj=sumj+vi(m,j1)*sumk
                 ENDDO
                 sumi=sumi+wi(i1)*sumj
              ENDDO
              vela(str,stt)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofilensv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF
!
! Extract E-W slice if required
!
IF(pvgew.EQ.1)THEN
!
! Allocate memory to velocity grid array
!
  nnx=(nnr-1)*dewr+1
  nnz=(nnp-1)*dewp+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Compute the values of the basis functions
!
  ALLOCATE(ui(dewr+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,dewr+1
     u=dewr
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(dewp+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,dewp+1
     u=dewp
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO
!
! Work out the value of u for the given latitude
!
  rnode=INT((got-slew)/rgst)+1
  u=ABS(got-slew-(rnode-1)*rgst)/rgst
  wi(1)=(1.0-u)**3/6.0
  wi(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
  wi(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
  wi(4)=u**3/6.0
  DO i=1,nnp-1
     conp=dewp
     IF(i==nnp-1)conp=dewp+1
     DO k=1,nnr-1
        conr=dewr
        IF(k==nnr-1)conr=dewr+1
        DO l=1,conp
           stp=dewp*(i-1)+l
           DO n=1,conr
              str=dewr*(k-1)+n
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumk=0.0
                    DO k1=1,4
                       rdm=ui(n,k1)*veln(k-2+k1,rnode-2+j1,i-2+i1)
                       sumk=sumk+rdm
                    ENDDO
                    sumj=sumj+wi(j1)*sumk
                 ENDDO
                 sumi=sumi+vi(l,i1)*sumj
              ENDDO
              vela(str,stp)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofileewv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF
DEALLOCATE(veln, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL veln'
ENDIF
!
! Read in input parameters for traveltime grid
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)ptfd
READ(10,'(a22)')ofiledt
READ(10,*)ptfns
READ(10,'(a22)')ofilenst
READ(10,*)ptfew
READ(10,'(a22)')ofileewt
!
! Open traveltime field file if required
!
IF(ptfd.EQ.1.OR.ptfns.EQ.1.OR.ptfew.EQ.1)THEN
   OPEN(UNIT=20,FILE=ifilet,FORM='unformatted',STATUS='old')
   READ(20)gor,got,gop
   READ(20)nnr,nnt,nnp
   READ(20)rgsr,rgst,rgsp
   ALLOCATE(ttn(nnp,nnt,nnr), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL ttn'
   ENDIF
   DO i=1,nnp
      DO j=1,nnt
         DO k=1,nnr
            READ(20)ttn(i,j,k)
         ENDDO
      ENDDO
   ENDDO
   CLOSE(20)
ENDIF
!
! Finish off GMT boundary files
!
WRITE(50,*)nnp
WRITE(50,*)nnt
WRITE(60,*)nnt
WRITE(60,*)nnr
WRITE(70,*)nnp
WRITE(70,*)nnr
CLOSE(50)
CLOSE(60)
CLOSE(70)
!
! Extract depth slice if required.
!
IF(ptfd.EQ.1)THEN
   idp=INT((gor-sldep)/rgsr)+1
   IF(idp.LT.1.OR.idp.GT.nnr)THEN
      WRITE(6,*)'Depth requested lies outside model'
      STOP
   ENDIF
!
!  Write out slice to GMT xyz file
!
   OPEN(UNIT=30,FILE=ofiledt,STATUS='unknown')
   DO i=1,nnp
      DO j=nnt,1,-1
         ttt=ttn(i,j,idp)
         IF(idp.EQ.nnr)THEN
            ttb=ttt
         ELSE
            ttb=ttn(i,j,idp+1)
         ENDIF
!
!        Apply linear interpolation to get travltime
!
         tt=ttt+(ttb-ttt)*((gor-(idp-1)*rgsr)-sldep)/rgsr
         WRITE(30,*)tt
      ENDDO
   ENDDO
   CLOSE(30)
ENDIF
!
! Extract N-S slice if required.
!
IF(ptfns.EQ.1)THEN
   ins=INT((slns-gop)/rgsp)+1
   IF(ins.LT.1.OR.ins.GT.nnp)THEN
      WRITE(6,*)'Longitude requested lies outside model'
      STOP
   ENDIF
   OPEN(UNIT=30,FILE=ofilenst,STATUS='unknown')
   DO i=1,nnt
      DO j=nnr,1,-1
         ttt=ttn(ins,i,j)
         IF(ins.EQ.nnp)THEN
            ttb=ttt
         ELSE
            ttb=ttn(ins+1,i,j)
         ENDIF
!
!        Apply linear interpolation to get velocity
!
         tt=ttt+(ttb-ttt)*(slns-(gop+(ins-1)*rgsp))/rgsp
         WRITE(30,*)tt
      ENDDO
   ENDDO
   CLOSE(30)
ENDIF
!
! Extract E-W slice if required.
!
IF(ptfew.EQ.1)THEN
   iew=INT((got-slew)/rgst)+1
   IF(iew.LT.1.OR.iew.GT.nnp)THEN
      WRITE(6,*)'Latitude requested lies outside model'
      STOP
   ENDIF
   OPEN(UNIT=30,FILE=ofileewt,STATUS='unknown')
   DO i=1,nnp
      DO j=nnr,1,-1
         ttt=ttn(i,iew,j)
         IF(iew.EQ.nnt)THEN
            ttb=ttt
         ELSE
            ttb=ttn(i,iew+1,j)
         ENDIF
!
!        Apply linear interpolation to get velocity
!
         tt=ttt+(ttb-ttt)*((got-(iew-1)*rgst)-slew)/rgst
         WRITE(30,*)ttt
      ENDDO
   ENDDO
   CLOSE(30)
ENDIF
IF(ptfd.EQ.1.OR.ptfns.EQ.1.OR.ptfew.EQ.1)THEN
   DEALLOCATE(ttn, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL ttn'
   ENDIF
ENDIF
!
! Read in ray path parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)prpd
READ(10,'(a22)')ofiledr
READ(10,*)prpns
READ(10,'(a22)')ofilensr
READ(10,*)prpew
READ(10,'(a22)')ofileewr
CLOSE(10)
!
! Plot raypaths in depth.
!
IF(prpd.EQ.1)THEN
   sep='>'
   OPEN(UNIT=20,FILE=irfile,FORM='unformatted',STATUS='old')
   READ(20)nr
   OPEN(UNIT=30,FILE=ofiledr,STATUS='unknown')
   DO i=1,nr
      READ(20)nre
      DO j=1,nre
         READ(20)rdep,rlat,rlon
         WRITE(30,'(2f10.4)')rlon,rlat
      ENDDO
      WRITE(30,'(a1)')sep
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF
!
! Plot raypaths in N-S.
!
IF(prpns.EQ.1)THEN
   sep='>'
   OPEN(UNIT=20,FILE=irfile,FORM='unformatted',STATUS='old')
   READ(20)nr
   OPEN(UNIT=30,FILE=ofilensr,STATUS='unknown')
   DO i=1,nr
      READ(20)nre
      DO j=1,nre
         READ(20)rdep,rlat,rlon
         WRITE(30,'(2f10.4)')rlat,rdep
      ENDDO
      WRITE(30,'(a1)')sep
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF
!
! Plot raypaths in E-W.
!
IF(prpew.EQ.1)THEN
   sep='>'
   OPEN(UNIT=20,FILE=irfile,FORM='unformatted',STATUS='old')
   READ(20)nr
   OPEN(UNIT=30,FILE=ofileewr,STATUS='unknown')
   DO i=1,nr
      READ(20)nre
      DO j=1,nre
         READ(20)rdep,rlat,rlon
         WRITE(30,'(2f10.4)')rlon,rdep
      ENDDO
      WRITE(30,'(a1)')sep
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF
STOP
END PROGRAM slice
