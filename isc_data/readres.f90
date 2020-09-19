!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This fortran90 code is designed to extract .res file in ISC-EHB Bulletin:
!               http://www.isc.ac.uk/isc-ehb/
! Written by Nanqiao Du, refer to .res format
module iscdata
implicit none

character(1000)       :: isc_path = '/home/nqdu/Documents/indonesia/iscdata/'

contains
subroutine same_earthquake(elat1,elon1,edep1,elat2,elon2,edep2,flag)
  implicit none
  real elat1,elon1,edep1,elat2,elon2,edep2

  real ela,elo,eld
  integer flag
  ela = abs(elat1 - elat2)
  elo = abs(elon1 - elon2)
  eld = abs(edep1-edep2)

  if((ela.le.1.0e-4).and.(elo.le.1.0e-4).and.(eld.le.1.0e-4)) then
      flag = 1
  else
      flag = 0
  endif
end

subroutine research_area(slat,slon,latmin,latmax,lonmin,lonmax,flag)
  implicit none
  real slat,slon,lonmin,latmin,lonmax,latmax
  integer flag

  if ((slat.le.latmax).and.(slat.ge.latmin).and. &
      (slon.le.lonmax).and.(slon.ge.lonmin)) then
      ! in research area 
      flag = 1
  else
      flag = 0
  endif
end

function geocen2geogra(colat) result(x)
  !! convert geocentric colat to geographic lat
  !! current coordinate system is WGS84
  implicit none
  real(4)               :: colat
  real(4)               :: x
  real(8)               :: temp
  real(8),PARAMETER     :: pi = atan(1.0d0) * 4.0
  real(8),PARAMETER     :: deg2rad = pi / 180.0
  real(8),PARAMETER     :: f = 1.0/298.257223563 !WGS84 flattening
  real(8)               :: latrad

  latrad = (90.0 - dble(colat)) * deg2rad
  if(-abs(latrad) + pi * 0.5 >= 0.05 ) then ! far away from poles
    temp = atan(tan(latrad) / (1-f)**2)
    x = real(temp/deg2rad)
  elseif(latrad > 0) then
    temp = pi * 0.5 - (1-f)**2 * (pi * 0.5 - latrad)
    x = real(temp / deg2rad)
  else
    temp = -pi * 0.5 + (1-f)**2 * (pi * 0.5 + latrad)
    x = real(temp / deg2rad)
  endif
  return;
end

subroutine read_res(startyear,endyear,evrg,strg,istele,degmin,degmax,phase)
  implicit none

  real     :: degmin,degmax
  real     :: evrg(4),strg(4) ! 1-4 (latmin,latmax,lonmin,lonmax)
  integer  :: startyear,endyear,istele

  integer  :: year
  integer  :: rnev,riyr,rimon,riday,rihold,rihr,rimin,&
              rntot,rntel,iphj,iphi,ipho,iprec,iflg,&
              eof,ievt
  real     :: ropenaz2,rropenaz2,topenaz2,rsec,elat,elon,rdepth,&
              fmb,fms,slat,slon,elev,delta,azim,rdtdd,rdelta,&
              razim,dbot,gblat,gblon,stadel,bdep,tbath,twater,&
              obstt,prett,rawres,ecor,scor,elcor,resid,wgt,&
              tdelta,ttime,delisc,resisc
  character*100 outfile,name1
  character*4  cTemp 
  character*3  isol
  character*6  sta
  character*3  risol
  character*2  riseq,comp
  character*1  onset,w
  character*8  phasej
  character*1000 fmt

  CHARACTER(8) phase

  integer flag,flag1
  real elat1,elon1,edep1,prec
  integer nout1,nout2
  integer success

  fmt = '(i7,1x,a3,a2,3f6.1,i5,2i3,i2,2i3,f6.2,2f8.3,f6.1,2f4.1,2i5,5x,&
          1x,a6,2f8.3,f7.3,2f8.3,5x,1x,a2,1x,a1,1x,a8,3i4,5x,&
          f8.4,2f8.3,f7.1,5x,3f8.3,f7.3,2f7.2,5x,f10.2,i3,f10.2,f7.2,5x,&
          4f7.2,i2,f5.2,5x,f8.3,f10.2,5x,f8.3,f7.2,1x,a1)'

  do year = startyear,endyear
    nout1 = 10
    nout2 = 30

    write(*,'(a,i5)')'begin converting year ',year
    write(cTemp,'(i4)') year
    name1=trim(adjustl(isc_path)) // trim(adjustl(cTemp))//'.res'
    outfile=trim(adjustl(cTemp))//'.out'
    open(nout1,file=outfile)
    open(nout2,file=name1)
    elat1 = 123;elon1 =-12312421; edep1= 3.1415926

    do
      read(nout2,fmt,iostat=eof) rnev,risol,riseq,ropenaz2,rropenaz2,topenaz2,&
              riyr,rimon,riday,rihold,rihr,rimin,rsec,&
              elat,elon,rdepth,fmb,fms,rntot,rntel,&
              sta,slat,slon,elev,delta,azim,&
              comp,onset,phasej,iphj,iphi,ipho,&
              rdtdd,rdelta,razim,dbot,&
              gblat,gblon,stadel,bdep,tbath,twater,&
              obstt,iprec,prett,rawres,&
              ecor,scor,elcor,resid,iflg,wgt,tdelta,ttime,&
              delisc,resisc,w,ievt
      if(eof/=0) exit

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! some restrictions for phase, station and earthquake
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! change geocentric latitude to geographic latitude
      elat = geocen2geogra(elat)
      slat = geocen2geogra(slat)

      ! check phase
      if( (phasej/=phase)) cycle
      if(iflg /= 0 ) cycle
      !write(*,*)phasej

      ! check location and source mechanism quality
      call research_area(elat,elon,evrg(1),evrg(2),evrg(3),evrg(4),flag)
      call research_area(slat,slon,strg(1),strg(2),strg(3),strg(4),flag1)
      if(istele.eq.1) then
          ! for teleseis case
          if((flag1.eq.0).or.(flag.eq.1).or.(isol=='XEQ'))  cycle
          if((delta<degmin).or.(delta>degmax))  cycle
      else
          ! for local and regional case
          if((flag1.eq.0).or.(flag.eq.0).or.(isol=='XEQ'))  cycle
          if((delta<degmin).or.(delta>degmax))  cycle
      endif

      ! check precision
      if(iprec.eq.-1) then
        prec = 0.1
      elseif(iprec.eq.-2) then
        prec = 0.01
      else if(iprec.eq.0) then
        prec = 1.0
      else
        cycle
      endif
      
      ! save earthquake information
      call same_earthquake(elat,elon,rdepth,elat1,elon1,edep1,flag)
      if(flag.eq.0) then
        write(nout1,'(a2,2f8.3,f6.1)') '#',elat,elon,rdepth
      endif

      ! save station information
      write(nout1,'(a6,a8,2f8.3,f7.3,f8.3,f10.2,4f7.2)')&
          sta,phasej,slat,slon,elev,delta,&
          obstt,rawres,prec,ecor,elcor
      ! renew
      elat1 = elat; elon1 = elon; edep1 = rdepth;
    enddo 

  ! close files
  close(nout1)
  close(nout2)
  enddo
  success = 1
end subroutine read_res

end module

program  main
  use iscdata
  implicit none
  integer         :: startyear,endyear,istele
  real            :: stlomin,stlamin,stlomax,stlamax
  real            :: evlomin,evlamin,evlomax,evlamax
  real            :: degmin,degmax

  real            :: evrg(4),strg(4)
  INTEGER         :: i
  CHARACTER(8)    :: phase

  ! read parameters
  open(12,file='restable.in')
  do i=1,3
    read(12,*)
  enddo
  read(12,*)startyear,endyear
  read(12,*)stlomin,stlomax,stlamin,stlamax
  read(12,*)evlomin,evlomax,evlamin,evlamax
  read(12,*)degmin,degmax
  read(12,*)istele
  read(12,*)phase
  close(12)

  ! construct event and station range
  evrg = (/evlamin,evlamax,evlomin,evlomax/)
  strg = (/stlamin,stlamax,stlomin,stlomax/)

  call read_res(startyear,endyear,evrg,strg,istele,degmin,degmax,phase)

  
end program  main
