SUBROUTINE output
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-04  Time: 21:21:36

!     *************************************************************************

!     OUTPUT
!     ======

!     AUTHOR
!     ------
!     R.S.CANT - CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     28-DEC-2003:  RSC MODIFIED FOR SENGA2
!     11-JUL-2009:  RSC ADD A DUMP FORMAT SWITCH; REPORT THE DUMP
!     29-AUG-2009:  RSC UPDATE NUMBER OF PROCESSORS

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     PROCESSES THE RESULTS

!     FILES: REPORT FILE       - UNIT NCREPT - NAME FNREPT  - FORMATTED
!            DUMP OUTPUT FILES - UNIT NCDMPO - NAMES FNDMPO - UNFORMATTED
!            STATISTICS FILE   - UNIT NCSTAT - NAME FNSTAT  - FORMATTED

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
#ifdef HDF5
use hdf5io
#endif
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL PARAMETERS
!     ================
!     DIAGNOSTICS
CHARACTER (LEN=4) :: pnxres
PARAMETER(pnxres = '.res')
CHARACTER (LEN=2) :: pnxryf
PARAMETER(pnxryf = 'yf')
INTEGER :: ncdiag
PARAMETER(ncdiag = 11)
INTEGER :: iddump

!     LOCAL DATA
!     ==========
!     DIAGNOSTICS
real(kind=8) :: deltag,fornow
real(kind=8) :: ttemp(nxsize,nysize,nzsize)
real(kind=8) :: ptemp(nxsize,nysize,nzsize)
real(kind=8) :: ytemp(nxsize,nysize,nzsize,nspec)
real(kind=8) :: vort1(nxsize,nysize,nzsize)
real(kind=8) :: vort2(nxsize,nysize,nzsize)
real(kind=8) :: vort3(nxsize,nysize,nzsize)
real(kind=8) :: enstro(nxsize,nysize,nzsize)
real(kind=8) :: tkel(nxsize,nysize,nzsize)

real(kind=8) :: tkes,tkeg,enstrs,enstrg

INTEGER :: ispec
INTEGER :: ic,jc,kc
INTEGER :: igofst,ix
CHARACTER (LEN=21) :: fndiag
!     RSC UPDATE NUMBER OF PROCESSORS
CHARACTER (LEN=6) :: pnproc
CHARACTER (LEN=11) :: strqty
CHARACTER (LEN=2) :: strspc
LOGICAL :: bcflag

CHARACTER (LEN=60) :: fname
CHARACTER (LEN=4) :: proc
CHARACTER (LEN=5) :: ipdump

!     BEGIN
!     =====

!     =========================================================================

!     UMOD END
!     =====
!     FULL DUMP OUTPUT
!     ================
  IF(MOD(itime,ntdump) == 0)THEN

!   OPS writing to HDF5
    call print_output()
  
!       CARRY OUT A FULL DUMP
!       ---------------------
!       USE THE DUMP FILE INDICATED BY IDFLAG
!       RSC 11-JUL-2009 ADD A DUMP FORMAT SWITCH
#ifndef HDF5
  IF(ndofmt == 0)THEN
   
    write(*,*) "Timestep: ", itime, "  Writing output to: ",fndmpo(idflag+1)
 
!         UNFORMATTED DUMP OUTPUT
    OPEN(UNIT=ncdmpo,FILE=fndmpo(idflag+1),STATUS='OLD', FORM='UNFORMATTED')
    REWIND(ncdmpo)
    WRITE(ncdmpo)nxnode,nynode,nznode,nspec, drun,urun,vrun,wrun,erun,yrun,  &
        etime,tstep,errold,errldr
    CLOSE(ncdmpo)
    
  ELSE
    
!         FORMATTED DUMP OUTPUT
    OPEN(UNIT=ncdmpo,FILE=fndmpo(idflag+1),STATUS='OLD', FORM='FORMATTED')
    REWIND(ncdmpo)
    WRITE(ncdmpo,*)nxnode,nynode,nznode,nspec
    DO kc = 1,nznode
      DO jc = 1,nynode
        DO ic = 1,nxnode
          WRITE(ncdmpo,*)drun(ic,jc,kc),  &
              urun(ic,jc,kc),vrun(ic,jc,kc),wrun(ic,jc,kc), erun(ic,jc,kc),  &
              (yrun(ic,jc,kc,ispec),ispec=1,nspec)
        END DO
      END DO
    END DO
    WRITE(ncdmpo,*)etime,tstep,errold,errldr
    CLOSE(ncdmpo)
    
  END IF
  
#else
  CALL write_h5_dumpfile
#endif

!       REPORT THE DUMP
!       RSC 11-JUL-2009
IF(iproc == 0)THEN
  
  OPEN(UNIT=ncrept,FILE=fnrept,STATUS='OLD',FORM='FORMATTED')
 3000      CONTINUE
 READ(ncrept,9000,END=3010)
 GO TO 3000
 3010      BACKSPACE(ncrept)
 WRITE(ncrept,9120)fndmpo(idflag+1)
CLOSE(ncrept)
  
END IF

!       RESET THE DUMP FLAG
idflag = MOD(idflag+1,2)

END IF

!     TIME STEP HISTORY
IF(iproc == 0)THEN
  WRITE(*,'(I7,1PE12.4,I5)')itime,tstep,inderr
END IF

 9000  FORMAT(a)
 9100  FORMAT('Time step number: ',i7)
 9110  FORMAT('Elapsed time: ',1PE12.4,';',2X,'next time step:',1PE12.4)
 9120  FORMAT('Dump completed: ',a)
 9200  FORMAT(i5,/ 5X,6(1PE12.4)/  &
    5X,6(1PE12.4)/ 5X,6(1PE12.4)/  &
   5X,4(1PE12.4)/ 5X,3(1PE12.4)/  &
    5X,3(1PE12.4)/ 5X,3(1PE12.4)/  &
    5X,3(1PE12.4)/ 5X,3(1PE12.4)/  &
    5X,2(1PE12.4))
 9300  FORMAT(2(1PE15.7))

END SUBROUTINE output
