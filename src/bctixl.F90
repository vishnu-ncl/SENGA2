SUBROUTINE bctixl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:58:49

!     *************************************************************************

!     BCTIXL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     30-MAR-2006:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     INITIALISES TURBULENT INFLOW

!     X-DIRECTION LEFT-HAND END

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     PARAMETERS
!     ==========
CHARACTER (LEN=4) :: pntixl,pntcxl
CHARACTER (LEN=4) :: pnxdat
PARAMETER(pntixl='tixl', pntcxl='tcxl', pnxdat='.dat')


!     LOCAL DATA
!     ==========
INTEGER :: icproc
INTEGER :: ic,jc,kc
INTEGER :: kxmodd,ixmodd
INTEGER :: itotxl
CHARACTER (LEN=6) :: pnproc
LOGICAL :: fxdump


!     BEGIN
!     =====

!     =========================================================================

!     BUILD THE FILENAMES FOR THE INLET TURBULENT VELOCITY FIELD
WRITE(pnproc,'(I6.6)')iproc
!     RESTART FILE
fntixl = pntixl//pnproc//pnxdat
!     COLD START FILE
fntcxl = pntcxl//pnproc//pnxdat

nctixl = 3

!     =========================================================================

!     INLET COLD START SWITCH
!     PARAMETER I2=0
IF(nxlprm(2) == 0)THEN
  
!       =======================================================================
  
!       INLET COLD START
!       ----------------
  
!       CHECK AND INITIALISE RESTART FILE
  INQUIRE(FILE=fntixl,EXIST=fxdump)
  IF(.NOT.fxdump)THEN
    OPEN(UNIT=nctixl,FILE=fntixl,STATUS='REPLACE',FORM='UNFORMATTED')
    CLOSE(nctixl)
  END IF
  
!       READ THE COLD START INLET TURBULENT VELOCITY FIELD
  OPEN(UNIT=nctixl,FILE=fntcxl,STATUS='OLD',FORM='UNFORMATTED')
  READ(nctixl)store1,store4,store5,store6
  CLOSE(nctixl)
  
!       SET THE REAL PARTS
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        urun(ic,jc,kc) = store4(ic,jc,kc)
        vrun(ic,jc,kc) = store5(ic,jc,kc)
        wrun(ic,jc,kc) = store6(ic,jc,kc)
        
      END DO
    END DO
  END DO
  
!       ZERO THE IMAGINARY PARTS
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        utmp(ic,jc,kc) = zero
        vtmp(ic,jc,kc) = zero
        wtmp(ic,jc,kc) = zero
        
      END DO
    END DO
  END DO
  
!       PARTIAL (X-WISE) FOURIER TRANSFORM
  CALL buftxl
  
!       MEAN INLET VELOCITY
!       SCANNING PLANE LOCATION AND VELOCITY
!       PARAMETER R1=MEAN VEL, R2=EXTRA SCAN VEL, R3=INITIAL LOCATION
  bvelxl = rxlprm(1)
  svelxl = bvelxl + rxlprm(2)
  elocxl = xgdlen - rxlprm(3)
  
!       =======================================================================
  
ELSE
  
!       =======================================================================
  
!       INLET RESTART
!       -------------
!       READ THE RESTART INLET TURBULENT VELOCITY FIELD
  OPEN(UNIT=nctixl,FILE=fntixl,STATUS='OLD',FORM='UNFORMATTED')
  READ(nctixl)ufxl,vfxl,wfxl,elocxl,svelxl,bvelxl
  CLOSE(nctixl)
  
!       =======================================================================
  
END IF

!     =========================================================================

!     INITIALISE RUNNING LOCATION
slocxl = elocxl

!     INITIALISE SCALE FACTORS
tpovxg = two*pi/xgdlen
scauxl = two/REAL(nxglbl,kind=8)
scduxl = -two*pi*scauxl*REAL(nxglbl-1,kind=8)/REAL(nxglbl,kind=8)/xgdlen

!     INITIALISE FLAGS AND INDICES FOR INLET PLANE DFT
fllixl = .false.
fltrxl = .false.
istaxl = 1
istoxl = npmapx(ixproc)

kminxl = 0
DO icproc = 0, ixproc-1
  kminxl = kminxl + npmapx(icproc)
END DO
kxmodd = MOD(kminxl,2)
itotxl = kminxl + npmapx(ixproc)

ixmodd = MOD(istoxl,2)

IF(kxmodd == 0)THEN
  
!       EVEN NUMBER OF POINTS TO LH SIDE OF THIS PROCESSOR
!       LH PROCESSOR HAS TRAILING REAL VALUE
!       LOCAL PROCESSOR MUST HAVE LEADING IMAGINARY VALUE
!       EXCEPT FOR FIRST POINT
  IF(kminxl > 0)fllixl = .true.
  kminxl = kminxl/2
  istaxl = 2
  IF(ixmodd == 0)THEN
!         EVEN NUMBER OF POINTS ON LOCAL PROCESSOR
!         LOCAL PROCESSOR MUST HAVE TRAILING REAL VALUE
    istoxl = istoxl - 1
    IF(itotxl == nxglbl)THEN
!           END OF PENCIL: IGNORE TRAILING REAL VALUE
      fltrxl = .false.
    ELSE
      fltrxl = .true.
    END IF
  ELSE
!         ODD NUMBER OF POINTS ON LOCAL PROCESSOR
!         LOCAL PROCESSOR HAS NO TRAILING REAL VALUE
    fltrxl = .false.
  END IF
  
ELSE
  
!       ODD NUMBER OF POINTS TO LH SIDE OF THIS PROCESSOR
!       LH PROCESSOR HAS NO TRAILING REAL VALUE
!       LOCAL PROCESSOR MUST HAVE NO LEADING IMAGINARY VALUE
  fllixl = .false.
  kminxl = kminxl/2 + 1
  istaxl = 1
  IF(ixmodd == 0)THEN
!         EVEN NUMBER OF POINTS ON LOCAL PROCESSOR
!         LOCAL PROCESSOR HAS NO TRAILING REAL VALUE
    fltrxl = .false.
  ELSE
!         ODD NUMBER OF POINTS ON LOCAL PROCESSOR
!         LOCAL PROCESSOR MUST HAVE TRAILING REAL VALUE
    istoxl = istoxl - 1
    IF(itotxl == nxglbl)THEN
!           END OF PENCIL: IGNORE TRAILING REAL VALUE
      fltrxl = .false.
    ELSE
      fltrxl = .true.
    END IF
  END IF
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bctixl
