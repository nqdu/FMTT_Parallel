module FMTTArray
implicit none

type,public :: FmttPair
  integer :: ns,nr,counter,num_data
  real,ALLOCATABLE    :: sourceloc(:,:)
  CHARACTER(8),ALLOCATABLE :: phase(:)
  real,ALLOCATABLE    :: stationloc(:,:)
end type

contains
subroutine pair_dealloc(p)
  implicit none
  type(FmttPair)  :: p
  DEALLOCATE(p%sourceloc,p%stationloc,p%phase)
end

end module
program inversion
  use lsmrModule
  use FMTTArray,only  : FmttPair,pair_dealloc
  implicit none

  real(4),PARAMETER   :: pi = 3.1415926535898
  CHARACTER (LEN=25)  :: gridi,gridc,otimes,mtimes,rtimes
  CHARACTER (LEN=25)  :: frechet,sources,receive
  CHARACTER (LEN=25)  :: subiter
  integer             :: i,j,k,rmtr,localsize,invstep,istep,istep1
  real(4)             :: epsilon, smooth, frac, mean,weight
  !integer             :: ridx,xidx,zidx

  ! velocity model defined
  real(4),ALLOCATABLE :: v(:),dv(:)
  integer             :: nx,nz,nr
  real(4)             :: dx,dz,dr,x0,z0,r0

  ! matrix defined
  real(4),ALLOCATABLE :: val(:) ! value of kernel matrix
  integer,ALLOCATABLE :: rw(:),col(:) ! row and col of kernel matrix
  real(4),ALLOCATABLE :: dref(:),dobs(:),dtrav(:),d(:),cd(:)
  integer,ALLOCATABLE :: pick(:)
  integer             :: nar,m,n,num,nrow,cidx,rwc,narray
  !integer,ALLOCATABLE :: nrc(:)
  real(4)             :: fdm

  ! temporary variables
  real(4)             :: temp,sintht,rgot,rgop,rdnt,rdnp,rad
  real(4),PARAMETER   :: earth=6371.0
  integer             :: nt,ns,nsta,nn,idm

  ! station-source pair 
  type(FmttPair),ALLOCATABLE :: pairs(:)

  !--------------------------------
  ! read control parameters
  !-------------------------------
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
  READ(10,*)epsilon ! damping parameter
  READ(10,*)localsize ! localsize
  READ(10,*)
  READ(10,*)smooth
  READ(10,*)
  READ(10,*)
  READ(10,*)
  READ(10,*)frac
  CLOSE(10)
  1 FORMAT(a)

  !-----------------------------
  ! read current inversion step
  !-----------------------------
  OPEN(UNIT=10,FILE=subiter,STATUS='old')
  READ(10,*)invstep
  CLOSE(10)

  !--------------------------
  ! read current model
  !--------------------------
  if(invstep == 1) then
    open(10,file=gridi)
  else
    open(10,file=gridc)
  endif
  READ(10,*)nr,nx,nz
  n = (nx+2) * (nz+2) * (nr+2)
  allocate(v(n),dv(n))
  read(10,*) r0,x0,z0
  read(10,*) dr,dx,dz

  ! compute corner points and intervals in rad
  rgop = z0 * pi / 180.
  rgot = (90.0 - x0) * pi / 180.
  rdnt = dx * pi / 180.
  rdnp = dz * pi / 180.

  ! read current model
  istep = 1
  do k=1,nz+2
    do j=1,nx+2
      do i=1,nr+2
        if(invstep == 1) then
          read(10,*) v(istep),temp
        else
          read(10,*) v(istep)
        endif
        istep = istep + 1
      enddo
    enddo
  enddo
  close(10)

  !------------------------------
  ! read source and receiver
  !-----------------------------
  m = 1
  open(11,file=sources)
  open(12,file=receive)
  read(11,*)narray
  read(12,*)
  ALLOCATE(pairs(narray))
  do i=1,narray
    read(11,*) ns
    pairs(i)%counter = m
    pairs(i)%ns= ns
    do j=1,ns
      read(11,*)
      read(11,*)
    enddo 
    read(12,*)nsta
    pairs(i)%nr = nsta
    m = m +  nsta * ns
    pairs(i)%num_data = nsta * ns
    do j=1,nsta
      read(12,*)
    enddo
  enddo
  close(11)
  close(12)
  m = m -1 ! no. of data in total!

  !-----------------------
  ! read no. of data
  !----------------------
  allocate(dref(m),dobs(m),cd(m),pick(m),dtrav(m))
  ! read reftimes
  open(10,file=rtimes)
  do i=1,m
    read(10,*)dref(i)
  enddo
  close(10)

  ! read residuals and weights
  open(10,file=otimes)
  do i=1,m
    read(10,*) pick(i),dobs(i),cd(i)
    !write(*,*)dres(i,1),dres(i,2),dres(i,3)
  enddo
  close(10)

  open(10,file=mtimes)
  do i=1,m
    read(10,*)dtrav(i)
  enddo
  close(10)
  
  ! check all valid rays
  open(10,file='valid_rays.dat')
  istep = 1
  do i=1,narray
    do j=1,pairs(i)%ns
      read(10,*)idm
      pick(istep:istep + pairs(i)%nr - 1) = pick(istep:istep + pairs(i)%nr -1 ) * idm
      istep = istep  + pairs(i)%nr
    enddo
  enddo
  close(10)

  
  !-----------------------------
  ! remove mean if required
  !----------------------------
  dtrav(1:m) = dtrav(1:m) - dref(1:m)
  dtrav(1:m) = dtrav(1:m) * real(pick(1:m))
  if(rmtr == 1) then
    istep = 1
    do i=1,narray
      istep = pairs(i)%counter
      istep1 = pairs(i)%nr
      do j=1,pairs(i)%ns
        nn = sum(pick(istep:istep+istep1 - 1))
        mean = sum(dtrav(istep:istep+istep1-1))/ real(nn)
        dtrav(istep:istep+istep1-1) =  dtrav(istep:istep+istep1-1)  - mean
        istep = istep + istep1
      enddo
    enddo
  endif

  ! add weights to residuals data
  cd(1:m) = cd(1:m)**2
  nt = sum(pick) ! no. of data used 
  ALLOCATE(d(nt+n))
  istep = 0
  do i=1,m
    if(pick(i) == 1) then
      istep = istep + 1
      d(istep) = (dobs(i) - dtrav(i)) / cd(i)
      if(cd(i) == 0.0) then
        print*, 'hhahahhaha',pick(i),cd(i)
        stop
      endif
    endif
  enddo

  ! -------------------
  ! read frechet kernel
  !--------------------
  num = int(frac * nt * n)
  allocate(rw(num),col(num),val(num))
  rw(1:num) = 0
  col(1:num) = 0
  val(1:num) = 0.0
  OPEN(UNIT=50,FILE=frechet,FORM='unformatted',STATUS='old')

  istep = 0 ! matrix index
  istep1 = 0 ! data index
  idm = 0 ! used-data index
  do i=1,narray
    do k=1,pairs(i)%num_data
      istep1 = istep1 + 1
      if(pick(istep1) == 0) then
        read(50) nrow
        do j=1,nrow
          read(50)cidx,fdm 
        enddo
      else
        idm = idm + 1
        weight = cd(istep1)
        read(50) nrow
        do j=1,nrow
          ! read one column of current frechet kernel
          read(50)cidx,fdm
          istep = istep + 1
          rw(istep) = idm
          col(istep) = cidx
          val(istep) = fdm / weight
        enddo
      endif
    enddo
  enddo
  close(50)
  nar = istep

  !open(40,file='invfre.dat')
  !do i=1,nar
  !  write(40,*)rw(i),col(i),val(i)
  !enddo
  !close(40)

  ! regularization terms
  istep = nt
  do i=1,nz+2
    do j=1,nx+2
      sintht = sin(rgot+(j-1)*rdnt)
      do k=1,nr+2
        rad =r0-(k-1)*dr+ earth
        if(i==1 .or. i==nz+2 .or. j==1.or. j==nx+2 .or. &
          k==1 .or. k==nr+2) then
          ! absorb variations at boundary points
          nar = nar + 1
          istep = istep + 1
          cidx = (i-1) * (nx+2) *(nr+2) + (j-1) * (nr+2) + k
          col(nar) = cidx
          val(nar) = 5.0 * smooth
          rw(nar) = istep
        else
          ! regularization in spherical coordinates
          istep = istep + 1
          cidx = (i-1) * (nx+2) *(nr+2) + (j-1) * (nr+2) + k
          rwc = istep

          ! r direction
          temp = 2.0*smooth *(1.0 +dr**2/(rad*rdnt)**2 &
                + dr**2/(rad*sintht*rdnp))
          val(nar + 1) = temp
          rw(nar + 1:nar+7) = rwc
          col(nar + 1) = cidx
          col(nar +2) = cidx-1
          val(nar +2) = -smooth
          col(nar +3) = cidx+1
          val(nar +3) = -smooth

          ! x direction
          col(nar + 4) = cidx-(nr+2)
          val(nar + 4) = -smooth*dr**2/(rad*rdnt)**2
          col(nar + 5) = cidx+(nr+2)
          val(nar + 5) = -smooth*dr**2/(rad*rdnt)**2

          ! z direction
          col(nar+6) = cidx-(nr+2)*(nx+2)
          val(nar + 6) = -smooth*dr**2/(rad*sintht*rdnp)**2
          col(nar+7) = cidx+(nr+2)*(nx+2)
          val(nar+7) = -smooth*dr**2/(rad*sintht*rdnp)**2

          ! renew nar
          nar = nar + 7
        endif
      enddo
    enddo
  enddo

  !--------------------------
  ! solve by lsmr
  !--------------------------
  tmp : block
  integer  :: itnlim,istop,itn
  real(4) :: atol,btol,conlim,anorm,acond,arnorm,xnorm,rnorm
  !real(4) :: wt(n)
  atol = 1.0e-6
  btol = 1.0e-6
  conlim = 1.0e6
  istop = 0
  anorm = 0.0
  acond = 0.0
  arnorm = 0.0
  xnorm = 0.0
  rnorm = 0.0;
  dv = 0.0
  itnlim = n * 2
  localsize = 50

  !wt(1:n) = 0.0

  m= nt
  call LSMR(nt+n,n,nar,rw,col,val,d,epsilon,atol,btol,conlim,itnlim,localsize,&
            dv,istop,itn,anorm,acond,rnorm,arnorm,xnorm)

  end block tmp

  ! renew model
  write(*,*)'max and min change ',maxval(dv),' ',minval(dv)
  where(dv < -0.5)
    dv = -0.5
  elsewhere(dv > 0.5)
    dv= 0.5
  endwhere
  v(1:n) = v(1:n) + dv(1:n)

  istep = 1
  ! write model
  open(40,file=gridc)
  write(40,*)nr,nx,nz
  write(40,*) r0,x0,z0
  write(40,*)dr,dx,dz
  do k=1,nz+2
    do j=1,nx+2
      do i=1,nr+2
        write(40,'(f10.5)') v(istep)
        istep = istep + 1
      enddo
    enddo
  enddo
  close(40)


  ! deallocate
  deallocate(v,dv)
  deallocate(dref,d,dobs,pick,cd,dtrav,pairs)
end
