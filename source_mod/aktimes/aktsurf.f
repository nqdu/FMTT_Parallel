ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Program akts
      implicit none
c
c     Determines ak135 traveltimes for a specified phase to
c     all surface grid nodes specified in 
c
c     By Nick Rawlinson (2/4/2004)
c     Australian National University
c
c     NOTE: Makes use of subroutines from the freeware
c           program ttimes.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer ii,i,j,k,l,n,phid,maxnp,gsn,nsrc,nrc
      integer ndep,nlat,nlon,nrfdep,nrflat,nrflon
      integer nra,idum,pors,ltfr
      parameter(maxnp=60)
      real tt(maxnp),dtdd(maxnp),dtdh(maxnp),dddp(maxnp)
      real evlat,evlon,evdep,ndlat,ndlon,rslat,usrc(2)
      real deltas,cazim,bazim,edist,bazr,etcor,tph,azima
      real odep,olat,olon,ddep,dlat,dlon,angi,pi,er,kmpd
      real*8 dodep,dolat,dolon,dddep,ddlat,ddlon,vels,velsp,velss
      real*8 cslat,cslong,abe,abw,abn,abs,rdum1,rdum2,rdum3
      parameter(pi=3.1415926535,er=6371.0)
      character*8 phcd(maxnp),phlst(maxnp),speph,sname
      character*25 modnam,sfile,pfile,efile,rfile,ofile
      logical prnt(2)
c
c     tt = traveltime of phase
c     dtdd, dtdh,dddp = partial derivatives
c     phcd = Phase id
c     phlst = Phase list
c     prnt = debugging print flag
c     modnam = Name of velocity model
c     evlat = event latitude
c     evlon = event longitude
c     evdep = event depth
c     ndlat = node latitude
c     ndlon = node longitude
c     rslat = station co-latitude
c     deltas = event-station angular distance
c     cazim,bazim,azima = azimuth information
c     n = number of phases found
c     edist = event-station distance
c     bazr = adjusted bazim
c     etcor = elliptical correction for travel time
c     phid = Phase id of specific phase
c     sfile = File containing grid vertex coordinates
c     pfile = File containing dicing factor information
c     efile = File containing event coordinates
c     rfile = File containing receiver coordinates
c     ofile = Output file containing ak135 traveltimes
c     sname = Station name
c     tph = traveltime of specific phase
c     maxnp = maximum number of ak135 phases
c     gsn = number of stations
c     speph = specific phase name for given event
c     nsrc = number of sources
c     ndep,nlat,nlon = number of nodes in depth,lat,lon
c     odep,olat,olon = grid origin in depth,lat,lon
c     ddep,dlat,dlon = grid spacing in depth,lat,lon
c     angi = angle of incidence of ray at surface
c     pi = pi
c     er = Earth radius
c     kmpd = Number of km per great circle degree
c     vels = velocity at surface
c     velsp = velocity at surface for P-wave
c     velss = velocity at surface for S-wave
c     nrfdep,nrflat,nrflon = B-spline dicing in r, theta, phi
c     nra = number of receiver arrays
c     nrc = number of receivers
c     abe,abw,abn,abs = horizontal array bounds (E, W, N, S)
c     pors = use P (0) or S (1) wavespeeds
c     ltfr = limit traveltime field to receiver array (0=no, 1=yes)
c
c     Set velocity at surface. NOTE that this assumes an
c     ak135 input model!!!!!
c
      velsp=5.8
      velss=3.46
      kmpd=2.0*pi*er/360.0
c
c     Open and read in the input parameter file
c
      open(unit=10,file='aktsurf.in',status='old')
      read(10,'(a25)')sfile
      read(10,'(a25)')pfile
      read(10,'(a25)')efile
      read(10,'(a25)')rfile
      read(10,*)pors
      read(10,'(a25)')ofile
      close(10)
c
c     Open the event file and read in the number of event types
c     (i.e. the number of receiver arrays)
c
      open(unit=20,file=efile,status='old')
      read(20,*)nra
c
c     Read in grid information
c
      open(unit=10,file=sfile,status='old')
      read(10,*)ndep,nlat,nlon
      read(10,*)dodep,dolat,dolon
      read(10,*)dddep,ddlat,ddlon
      close(10)
c
c     Open parameter file and read in dicing factors and
c     receiver cushion paprameters
c
      open(unit=10,file=pfile,status='old')
      do i=1,7
         read(10,*)
      enddo
      read(10,*)nrfdep,nrflat,nrflon
      do i=1,3
         read(10,*)
      enddo
      read(10,*)ltfr
      read(10,*)cslat,cslong
      close(10)
      ndep=(ndep-1)*nrfdep+1
      dddep=dddep/nrfdep
      ddlat=ddlat/nrflat
      ddlon=ddlon/nrflon
      odep=dodep
      ddep=dddep
c
c     Specify model type
c
      modnam='ak135'
c
c     Open travel time tables
c
      prnt(1)=.false.
      prnt(2)=.false.
      phlst(2)='  '
      call tabin(10,modnam)
c
c     Begin a loop over all events (sources) to generate
c     nsrc surface grids of ak135 traveltimes
c
      open(unit=30,file=ofile,form='unformatted',status='unknown')
      open(unit=40,file=rfile,status='old')
      read(40,*)idum
      if(idum.ne.nra)then
         write(6,*)'Receiver file is incompatible with'
         write(6,*)'the source file!!!!'
         write(6,*)'First line of both files should be'
         write(6,*)'identical!!!!'
         write(6,*)'TERMINATING PROGRAM!!!'
         stop
      endif
      do ii=1,nra
c
c        Read in the list of receiver coordinates
c        and identify the horizontal limits of the
c        array
c
         if(ltfr.eq.1)then
            read(40,*)nrc
            do i=1,nrc
               read(40,*)rdum1,rdum2,rdum3
               if(i.eq.1)then
                  abe=rdum3
                  abw=abe
                  abn=rdum2
                  abs=abn
               else
                  if(rdum3.lt.abw)then
                     abw=rdum3
                  else if(rdum3.gt.abe)then
                     abe=rdum3
                  endif
                  if(rdum2.lt.abs)then
                     abs=rdum2
                  else if(rdum2.gt.abn)then
                     abn=rdum2
                  endif
              endif
            enddo
c
c           Add cushion value
c
            abw=abw-cslong
            abe=abe+cslong
            abs=abs-cslat
            abn=abn+cslat
c
c           Now set origin of the grid
c 
            olat=abn
            olon=abw
c
c           Determine the number of points in latitude
c           and longitude which span the array
c
            nlat=INT((abn-abs)/ddlat)+2
            nlon=INT((abe-abw)/ddlon)+2
            dlat=(abn-abs)/REAL(nlat-2)
            dlon=(abe-abw)/REAL(nlon-2)
         else
            olat=dolat
            olon=dolon
            nlat=nlat*nrflat-1
            nlon=nlon*nrflon-1
            dlat=ddlat
            dlon=ddlon
         endif
         read(20,*)nsrc
         write(30)olat,olon
         write(30)nlat,nlon
         write(30)dlat,dlon
         If(pors.eq.0)then
            vels=velsp
         else
            vels=velss
         endif
         do i=1,nsrc
c
c           Read in source location and phase type
c
            read(20,*)evlat,evlon,evdep
            read(20,*)speph
            phlst(1)=speph
            call brnset(1,phlst,prnt)
            rslat = (90.-evlat)*0.01745329252
            call ellref(rslat)
            call depset(evdep,usrc)
c
c           Now loop through all the grid points
c
            do j=1,nlon
               ndlon=olon+(j-1)*dlon
c
c              This is a bit of a hack fix
c              for the case when event and grid node
c              longitude are equal. For some reason
c              the ellipticity corrections are wrong
c              for this case
c
               if(ndlon.eq.evlon)then
                  ndlon=ndlon+dlon/50.0
               endif
               do k=1,nlat
                  ndlat=olat-(k-1)*dlat
c
c                 Determine phase traveltimes.
c
          call ydist(evlat,evlon,ndlat,ndlon,deltas,cazim,bazim,azima)
                  call trtm(deltas,maxnp,n,tt,dtdd,dtdh,dddp,phcd)
c
c                 Extract specified phase.
c
                  l=1
                  phid=0
                  do while(l.le.n.and.phid.eq.0)
                     if(phcd(l).eq.speph)then
                        phid=l
                     else
                        l=l+1
                     endif
                  enddo
                  if(phid.eq.0)then
                     write(6,*)'Requested phase ',speph, 'does not'
                     write(6,*)'exist for delta = ',deltas
                     write(6,*)'Station number = ',i
                     write(6,*)'Array number = ',ii
                     write(6,*)'Terminating program'
                     write(6,*)'source coordinates',evlat,evlon
                     write(6,*)'receiver coordinates',rdum2,rdum3
                     stop
                  endif
                  tph=tt(phid)
c
c                 Apply elliptical correction
c
                  edist=deltas*0.017453292
                  bazr=bazim*0.017453292
                  call ellcor(edist, bazr, evdep, speph, etcor)
                  tph=tph+etcor
c
c                 Output traveltime, ray azimuth and incidence 
c                 angle to file
c
                  angi=asin(vels*dtdd(phid)/kmpd)*180.0/pi
                  write(30)tph,azima,angi
               enddo
            enddo
         enddo
      enddo
      close(20)
      close(30)
      close(40)
      stop
      end
