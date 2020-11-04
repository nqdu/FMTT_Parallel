!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program outputs traveltime residual data in a form
! suitable for input to GMT.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM resplot
IMPLICIT NONE
INTEGER :: ii,i,j,idm,iorf,isum,invst,itern,rmtr,stnm,iscn
INTEGER :: nr,ns,nt,nra,idum,checkstat
INTEGER, DIMENSION (:), ALLOCATABLE :: ntsrc,ots
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: tres,prs,rd1,rd2,rd3,rdum
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: stt,srcm
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: mtim,rtim,otim
CHARACTER (LEN=40) :: mtimes,rtimes,otimes,sources,receivers,ofile
CHARACTER (LEN=40) :: stfile,itfile,subinv
CHARACTER (LEN=10), DIMENSION(:), ALLOCATABLE :: recn
!
! rtimes = input file of reference times
! otimes = output file of traveltime residuals
! sources = file containing source information
! receivers = file containing receiver information
! mtimes = file containing model traveltimes
! ofile = output file for plotting
! ns= number of sources
! nr = number of receivers
! nt = number of traveltimes
! iorf = initial (0) or final (1) residuals plotted
! prs = print residuals larger than this size
! recn = recorder name
! stfile = File containing station terms
! itfile = FIle containing inversion iteration number
! subinv = file containing subspace parameters
! invst = Invert for station terms (0=no, 1=yes)
! stt = Station traveltime term
! itern = Inversion iteration number
! rmtr = Remove traveltime residual (0=no, 1=yes)
! srcm = model traveltime residual mean for each source
! ntsrc = Number of traveltimes for each source
! stnm = Station number
! iscn = Source number
! nra = Number of receiver arrays
! mtim = Model traveltime array
! rtim = Reference traveltime array
! otim = Observed traveltime array
! ots = Observed traveltime status (0 or 1)
!
! Read in the input parameters
!
OPEN(UNIT=10,FILE='resplott.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)iorf
READ(10,'(a32)')ofile
READ(10,*)prs
READ(10,'(a32)')otimes
READ(10,'(a32)')rtimes
READ(10,'(a32)')subinv
READ(10,'(a32)')stfile
READ(10,'(a32)')itfile
READ(10,'(a32)')sources
READ(10,'(a32)')receivers
READ(10,'(a32)')mtimes
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
IF(invst.EQ.1)THEN
   OPEN(UNIT=10,FILE=itfile,STATUS='old')
   READ(10,*)itern
   CLOSE(10)
   IF(itern.GT.1)THEN
      OPEN(UNIT=70,FILE=stfile,STATUS='old')
   ENDIF
ENDIF
OPEN(UNIT=10,FILE=mtimes,STATUS='old')
OPEN(UNIT=20,FILE=rtimes,STATUS='old')
OPEN(UNIT=30,FILE=otimes,STATUS='old')
OPEN(UNIT=40,FILE=ofile,STATUS='unknown')
DO ii=1,nra
   READ(50,*)ns
   READ(60,*)nr
   ALLOCATE(recn(nr), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM resplot: CHAR recn'
   ENDIF
   DO i=1,nr
      READ(60,*)rd1,rd2,rd3,recn(i)
   ENDDO
   nt=ns*nr
!
!  If station terms are applied, then read in
!  residual.
!
   IF(invst.EQ.1)THEN
      ALLOCATE(stt(nr))
      stt=0.0
      IF(itern.GT.1)THEN
         DO i=1,nr
            READ(70,*)stt(i)
         ENDDO
      ENDIF
   ENDIF
!
!  Read in the model and reference traveltimes and take
!  the difference and then subtract the observed residuals.
!
   ALLOCATE(mtim(nt))
   ALLOCATE(rtim(nt))
   ALLOCATE(otim(nt))
   ALLOCATE(ots(nt))
   DO i=1,nt
      READ(10,*)mtim(i)
      READ(20,*)rtim(i)
      READ(30,*)ots(i),otim(i)
   ENDDO
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
!  Read in the model and reference traveltimes and take
!  the difference and then subtract the observed residuals.
!
   isum=1
   DO i=1,ns
      DO j=1,nr
         rdum=mtim(isum)
         IF(invst.EQ.1)rdum=rdum+stt(j)
         IF(ots(isum).EQ.1)THEN
            IF(iorf.EQ.0)THEN
               tres=otim(isum)
            ELSE
               IF(rmtr.EQ.1)THEN
                  tres=otim(isum)-(rdum-rtim(isum))+srcm(i)
               ELSE
                  tres=otim(isum)-(rdum-rtim(isum))
               ENDIF
            ENDIF
            WRITE(40,*)tres
            IF(ABS(tres).GT.prs)THEN
               WRITE(6,1)recn(j),ii,i,tres
1              FORMAT('Recorder ',a7,' in array ',i3, &
              &' for source ',i5,' has a &
              &residual size of ',f7.4,' s.')
            ENDIF
         ENDIF
         isum=isum+1
      ENDDO
   ENDDO
   DEALLOCATE(recn, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: PROGRAM resplot: CHAR recn'
   ENDIF
   DEALLOCATE(mtim,rtim,otim,ots)
   IF(invst.EQ.1)DEALLOCATE(stt)
   IF(rmtr.EQ.1)DEALLOCATE(srcm,ntsrc)
ENDDO
CLOSE(10)
CLOSE(20)
CLOSE(30)
CLOSE(40)
CLOSE(50)
CLOSE(60)
IF(invst.EQ.1)THEN
   IF(itern.GT.1)CLOSE(70)
ENDIF
END PROGRAM resplot
