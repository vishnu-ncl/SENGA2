SUBROUTINE bcinit
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:24:34

!     *************************************************************************

!     BCINIT
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     26-MAR-2003:  CREATED
!     26-OCT-2008:  RSC/TDD BUG FIX BCTIXL
!     23-JAN-2015:  RSC BUG FIX WALL MASS FLUX BC
!     19-MAY-2015:  RSC UPDATE WALL DIFFUSIVE BC

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     CARRIES OUT INITIALISATION OF BOUNDARY CONDITIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     BEGIN
!     =====

!     =========================================================================

!     APPLY GLOBAL BOUNDARY CONDITIONS TO LOCAL PROCESSOR
!     ===================================================

!     X-LEFT
IF(ixproc == 0)THEN
  nsbcxl = ngbcxl
  nendxl = nbound
  IF(nsbcxl == nsperi)THEN
    nendxl = nobc
    IF(nxproc == 1)nendxl = nperi
  END IF
ELSE
  nsbcxl = nsnull
  nendxl = nobc
END IF

!     X-RIGHT
IF(ixproc == nxprm1)THEN
  nsbcxr = ngbcxr
  nendxr = nbound
  IF(nsbcxr == nsperi)THEN
    nendxr = nobc
    IF(nxproc == 1)nendxr = nperi
  END IF
ELSE
  nsbcxr = nsnull
  nendxr = nobc
END IF

!     Y-LEFT
IF(iyproc == 0)THEN
  nsbcyl = ngbcyl
  nendyl = nbound
  IF(nsbcyl == nsperi)THEN
    nendyl = nobc
    IF(nyproc == 1)nendyl = nperi
  END IF
ELSE
  nsbcyl = nsnull
  nendyl = nobc
END IF

!     Y-RIGHT
IF(iyproc == nyprm1)THEN
  nsbcyr = ngbcyr
  nendyr = nbound
  IF(nsbcyr == nsperi)THEN
    nendyr = nobc
    IF(nyproc == 1)nendyr = nperi
  END IF
ELSE
  nsbcyr = nsnull
  nendyr = nobc
END IF

!     Z-LEFT
IF(izproc == 0)THEN
  nsbczl = ngbczl
  nendzl = nbound
  IF(nsbczl == nsperi)THEN
    nendzl = nobc
    IF(nzproc == 1)nendzl = nperi
  END IF
ELSE
  nsbczl = nsnull
  nendzl = nobc
END IF

!     Z-RIGHT
IF(izproc == nzprm1)THEN
  nsbczr = ngbczr
  nendzr = nbound
  IF(nsbczr == nsperi)THEN
    nendzr = nobc
    IF(nzproc == 1)nendzr = nperi
  END IF
ELSE
  nsbczr = nsnull
  nendzr = nobc
END IF

!     =========================================================================

!     INITIALISE PARALLEL TRANSFER BC FLAGS
!     =====================================
prgoxl = .false.
IF(nendxl == nobc)prgoxl = .true.
prgoxr = .false.
IF(nendxr == nobc)prgoxr = .true.
prgoyl = .false.
IF(nendyl == nobc)prgoyl = .true.
prgoyr = .false.
IF(nendyr == nobc)prgoyr = .true.
prgozl = .false.
IF(nendzl == nobc)prgozl = .true.
prgozr = .false.
IF(nendzr == nobc)prgozr = .true.

!     =========================================================================

!     INITIALISE BC FLAGS
!     ===================
!     FLAG = .TRUE. INDICATES THAT TERMS ARE TO BE SET TO ZERO ON THE BOUNDARY
!     F--CNV CONVECTIVE TERMS
!     F--VSN NORMAL VISCOUS TERMS
!     F--VST TANGENTIAL VISCOUS TERMS
!     F--CON THERMAL CONDUCTION TERMS
!     F--ADB OTHER HEAT FLUX TERMS
!     F--DIF MASS DIFFUSION TERMS

fxlcnv = .false.
fxlvsn = .false.
fxlvst = .false.
fxlcon = .false.
fxladb = .false.
fxldif = .false.
!     RSC 19-MAY-2015
fxldfw = .false.
!     RSC 28-JUN-2015
fxlcnw = .false.
fxladw = .false.

fxrcnv = .false.
fxrvsn = .false.
fxrvst = .false.
fxrcon = .false.
fxradb = .false.
fxrdif = .false.
!     RSC 19-MAY-2015
fxrdfw = .false.
!     RSC 28-JUN-2015
fxrcnw = .false.
fxradw = .false.

fylcnv = .false.
fylvsn = .false.
fylvst = .false.
fylcon = .false.
fyladb = .false.
fyldif = .false.
!     RSC 19-MAY-2015
fyldfw = .false.
!     RSC 28-JUN-2015
fylcnw = .false.
fyladw = .false.

fyrcnv = .false.
fyrvsn = .false.
fyrvst = .false.
fyrcon = .false.
fyradb = .false.
fyrdif = .false.
!     RSC 19-MAY-2015
fyrdfw = .false.
!     RSC 28-JUN-2015
fyrcnw = .false.
fyradw = .false.

fzlcnv = .false.
fzlvsn = .false.
fzlvst = .false.
fzlcon = .false.
fzladb = .false.
fzldif = .false.
!     RSC 19-MAY-2015
fzldfw = .false.
!     RSC 28-JUN-2015
fzlcnw = .false.
fzladw = .false.

fzrcnv = .false.
fzrvsn = .false.
fzrvst = .false.
fzrcon = .false.
fzradb = .false.
fzrdif = .false.
!     RSC 19-MAY-2015
fzrdfw = .false.
!     RSC 28-JUN-2015
fzrcnw = .false.
fzradw = .false.

!     =========================================================================

!     FLAGS FOR TURBULENT INFLOW VELOCITY FIELD
!     =========================================
fxltrb = .false.
fxrtrb = .false.
fyltrb = .false.
fyrtrb = .false.
fzltrb = .false.
fzrtrb = .false.

!     =========================================================================

!     INITIALISE RK COUNTERS
!     ======================
istald = istal
istold = istol
jstald = jstal
jstold = jstol
kstald = kstal
kstold = kstol

istalu = istal
istolu = istol
jstalu = jstal
jstolu = jstol
kstalu = kstal
kstolu = kstol

istalv = istal
istolv = istol
jstalv = jstal
jstolv = jstol
kstalv = kstal
kstolv = kstol

istalw = istal
istolw = istol
jstalw = jstal
jstolw = jstol
kstalw = kstal
kstolw = kstol

istale = istal
istole = istol
jstale = jstal
jstole = jstol
kstale = kstal
kstole = kstol

istaly = istal
istoly = istol
jstaly = jstal
jstoly = jstol
kstaly = kstal
kstoly = kstol

istalt = istal-nhalox
istolt = istol+nhalox
jstalt = jstal-nhaloy
jstolt = jstol+nhaloy
kstalt = kstal-nhaloz
kstolt = kstol+nhaloz

!     =========================================================================

!     SET BC FLAGS, RK COUNTERS AND PHYSICAL PARAMETERS
!     =================================================

!     =========================================================================

!     OUTFLOW BC No 1
!     ---------------
!     SUBSONIC NON-REFLECTING OUTFLOW
!     WITH OPTION TO SET PRESSURE AT INFINITY
!     EQUIVALENT TO NSCBC OUTFLOW 1
!     PARAMETERS: R1=P AT INFINITY, R2=COEFF, R3=MACH NO
!     X-LEFT
IF(nsbcxl == nsbco1)THEN
  fxlcnv = .true.
  fxlvst = .true.
  fxlcon = .true.
  fxladb = .true.
  fxldif = .true.
!       RSC 19-MAY-2015
!       RSC 28-JUN-2015
!        FXLDFW = .TRUE.
  istalt = istal
  cobcxl = rxlprm(2)*(one-rxlprm(3))/(two*xgdlen)
  pinfxl = rxlprm(1)
END IF

!     X-RIGHT
IF(nsbcxr == nsbco1)THEN
  fxrcnv = .true.
  fxrvst = .true.
  fxrcon = .true.
  fxradb = .true.
  fxrdif = .true.
!       RSC 19-MAY-2015
!       RSC 28-JUN-2015
!        FXRDFW = .TRUE.
  istolt = istol
  cobcxr = rxrprm(2)*(one-rxrprm(3))/(two*xgdlen)
  pinfxr = rxrprm(1)
END IF

!     Y-LEFT
IF(nsbcyl == nsbco1)THEN
  fylcnv = .true.
  fylvst = .true.
  fylcon = .true.
  fyladb = .true.
  fyldif = .true.
!       RSC 19-MAY-2015
!       RSC 28-JUN-2015
!        FYLDFW = .TRUE.
  jstalt = jstal
  cobcyl = rylprm(2)*(one-rylprm(3))/(two*ygdlen)
  pinfyl = rylprm(1)
END IF

!     Y-RIGHT
IF(nsbcyr == nsbco1)THEN
  fyrcnv = .true.
  fyrvst = .true.
  fyrcon = .true.
  fyradb = .true.
  fyrdif = .true.
!       RSC 19-MAY-2015
!       RSC 28-JUN-2015
!        FYRDFW = .TRUE.
  jstolt = jstol
  cobcyr = ryrprm(2)*(one-ryrprm(3))/(two*ygdlen)
  pinfyr = ryrprm(1)
END IF

!     Z-LEFT
IF(nsbczl == nsbco1)THEN
  fzlcnv = .true.
  fzlvst = .true.
  fzlcon = .true.
  fzladb = .true.
  fzldif = .true.
!       RSC 19-MAY-2015
!       RSC 28-JUN-2015
!        FZLDFW = .TRUE.
  kstalt = kstal
  cobczl = rzlprm(2)*(one-rzlprm(3))/(two*zgdlen)
  pinfzl = rzlprm(1)
END IF

!     Z-RIGHT
IF(nsbczr == nsbco1)THEN
  fzrcnv = .true.
  fzrvst = .true.
  fzrcon = .true.
  fzradb = .true.
  fzrdif = .true.
!       RSC 19-MAY-2015
!       RSC 28-JUN-2015
!        FZRDFW = .TRUE.
  kstolt = kstol
  cobczr = rzrprm(2)*(one-rzrprm(3))/(two*zgdlen)
  pinfzr = rzrprm(1)
END IF

!     =========================================================================

!     INFLOW BC no 1
!     --------------
!     SUBSONIC NON-REFLECTING LAMINAR INFLOW
!     EQUIVALENT TO NSCBC INFLOW 4

!     X-LEFT
IF(nsbcxl == nsbci1)THEN
  fxlcnv = .true.
  fxlvsn = .true.
  istalt = istal
END IF

!     X-RIGHT
IF(nsbcxr == nsbci1)THEN
  fxrcnv = .true.
  fxrvsn = .true.
  istolt = istol
END IF

!     Y-LEFT
IF(nsbcyl == nsbci1)THEN
  fylcnv = .true.
  fylvsn = .true.
  jstalt = jstal
END IF

!     Y-RIGHT
IF(nsbcyr == nsbci1)THEN
  fyrcnv = .true.
  fyrvsn = .true.
  jstolt = jstol
END IF

!     Z-LEFT
IF(nsbczl == nsbci1)THEN
  fzlcnv = .true.
  fzlvsn = .true.
  kstalt = kstal
END IF

!     Z-RIGHT
IF(nsbczr == nsbci1)THEN
  fzrcnv = .true.
  fzrvsn = .true.
  kstolt = kstol
END IF

!     =========================================================================

!     INFLOW BC No 2
!     --------------
!     SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
!     WITH OPTION FOR INFLOW TURBULENCE
!     EQUIVALENT TO NSCBC INFLOW 1
!     AS IMPLEMENTED BY SUTHERLAND+KENNEDY

!     X-LEFT

!     GLOBAL BC SUPPORT
!     TURBULENT INFLOW VELOCITY FIELD
IF(ngbcxl == nsbci2)THEN
  IF(nxlprm(1) == 3)THEN
    fxltrb = .true.
    CALL bctixl
  END IF
END IF

!     LOCAL BC SUPPORT
IF(nsbcxl == nsbci2)THEN
  
  fxlcnv = .true.
  istalu = istap1
  istalv = istap1
  istalw = istap1
  istale = istap1
  istaly = istap1
  istalt = istal
  fxltrb = .false.
  
END IF

!     X-RIGHT
IF(nsbcxr == nsbci2)THEN
  
  fxrcnv = .true.
  istolu = istom1
  istolv = istom1
  istolw = istom1
  istole = istom1
  istoly = istom1
  istolt = istol
  
END IF

!     Y-LEFT
IF(nsbcyl == nsbci2)THEN
  
  fylcnv = .true.
  jstalu = jstap1
  jstalv = jstap1
  jstalw = jstap1
  jstale = jstap1
  jstaly = jstap1
  jstalt = jstal
  
END IF

!     Y-RIGHT
IF(nsbcyr == nsbci2)THEN
  
  fyrcnv = .true.
  jstolu = jstom1
  jstolv = jstom1
  jstolw = jstom1
  jstole = jstom1
  jstoly = jstom1
  jstolt = jstol
  
END IF

!     Z-LEFT
IF(nsbczl == nsbci2)THEN
  
  fzlcnv = .true.
  kstalu = kstap1
  kstalv = kstap1
  kstalw = kstap1
  kstale = kstap1
  kstaly = kstap1
  kstalt = kstal
  
END IF

!     Z-RIGHT
IF(nsbczr == nsbci2)THEN
  
  fzrcnv = .true.
  kstolu = kstom1
  kstolv = kstom1
  kstolw = kstom1
  kstole = kstom1
  kstoly = kstom1
  kstolt = kstol
  
END IF

!     RSC/TDD BUG FIX BCTIXL
!C     INITIALISE INLET TURBULENCE GENERATOR
!      IF(FXLTRB)CALL BCTIXL

!     =========================================================================

!     INFLOW BC No 3
!     --------------
!     SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
!     WITH OPTION FOR INFLOW TURBULENCE
!     EQUIVALENT TO NSCBC INFLOW 2

!     X-LEFT

!     GLOBAL BC SUPPORT
!     X-LEFT TURBULENT INFLOW VELOCITY FIELD
IF(ngbcxl == nsbci3)THEN
  IF(nxlprm(1) == 3)THEN
    fxltrb = .true.
    CALL bctixl
  END IF
END IF

!     LOCAL BC SUPPORT
IF(nsbcxl == nsbci3)THEN
  
  fxlcnv = .true.
  fxlvsn = .true.
  istald = istap1
  istalu = istap1
  istalv = istap1
  istalw = istap1
  istaly = istap1
  istalt = istal
  fxltrb = .false.
  
END IF

!     X-RIGHT
IF(nsbcxr == nsbci3)THEN
  
  fxrcnv = .true.
  fxrvsn = .true.
  istold = istom1
  istolu = istom1
  istolv = istom1
  istolw = istom1
  istoly = istom1
  istolt = istol
  
END IF

!     Y-LEFT
IF(nsbcyl == nsbci3)THEN
  
  fylcnv = .true.
  fylvsn = .true.
  jstald = jstap1
  jstalu = jstap1
  jstalv = jstap1
  jstalw = jstap1
  jstaly = jstap1
  jstalt = jstal
  
END IF

!     Y-RIGHT
IF(nsbcyr == nsbci3)THEN
  
  fyrcnv = .true.
  fyrvsn = .true.
  jstold = jstom1
  jstolu = jstom1
  jstolv = jstom1
  jstolw = jstom1
  jstoly = jstom1
  jstolt = jstol
  
END IF

!     Z-LEFT
IF(nsbczl == nsbci3)THEN
  
  fzlcnv = .true.
  fzlvsn = .true.
  kstald = kstap1
  kstalu = kstap1
  kstalv = kstap1
  kstalw = kstap1
  kstaly = kstap1
  kstalt = kstal
  
END IF

!     Z-RIGHT
IF(nsbczr == nsbci3)THEN
  
  fzrcnv = .true.
  fzrvsn = .true.
  kstold = kstom1
  kstolu = kstom1
  kstolv = kstom1
  kstolw = kstom1
  kstoly = kstom1
  kstolt = kstol
  
END IF

!     =========================================================================

!     WALL BC No 1
!     ------------
!     NO-SLIP WALL - ADIABATIC

!     X-LEFT
IF(nsbcxl == nsbcw1)THEN
  
  fxlcnv = .true.
  fxlcon = .true.
  fxladb = .true.
!       RSC 23-JAN-2015
  fxldif = .true.
!       RSC 19-MAY-2015
  fxlcnw = .true.
  fxldfw = .true.
  istow  = istap4
  istalu = istap1
  istalv = istap1
  istalw = istap1
  istalt = istal
  
END IF

!     X-RIGHT
IF(nsbcxr == nsbcw1)THEN
  
  fxrcnv = .true.
  fxrcon = .true.
  fxradb = .true.
!       RSC 23-JAN-2015
  fxrdif = .true.
!       RSC 19-MAY-2015
  fxrcnw = .true.
  fxrdfw = .true.
  istaw  = istom4
  istolu = istom1
  istolv = istom1
  istolw = istom1
  istolt = istol
  
END IF

!     Y-LEFT
IF(nsbcyl == nsbcw1)THEN
  
  fylcnv = .true.
  fylcon = .true.
  fyladb = .true.
!       RSC 23-JAN-2015
  fyldif = .true.
!       RSC 19-MAY-2015
  fylcnw = .true.
  fyldfw = .true.
  jstow  = jstap4
  jstalu = jstap1
  jstalv = jstap1
  jstalw = jstap1
  jstalt = jstal
  
END IF

!     Y-RIGHT
IF(nsbcyr == nsbcw1)THEN
  
  fyrcnv = .true.
  fyrcon = .true.
  fyradb = .true.
!       RSC 23-JAN-2015
  fyrdif = .true.
!       RSC 19-MAY-2015
  fyrcnw = .true.
  fyrdfw = .true.
  jstaw  = jstom4
  jstolu = jstom1
  jstolv = jstom1
  jstolw = jstom1
  jstolt = jstol
  
END IF

!     Z-LEFT
IF(nsbczl == nsbcw1)THEN
  
  fzlcnv = .true.
  fzlcon = .true.
  fzladb = .true.
!       RSC 23-JAN-2015
  fzldif = .true.
!       RSC 19-MAY-2015
  fzlcnw = .true.
  fzldfw = .true.
  kstow  = kstap4
  kstalu = kstap1
  kstalv = kstap1
  kstalw = kstap1
  kstalt = kstal
  
END IF

!     Z-RIGHT
IF(nsbczr == nsbcw1)THEN
  
  fzrcnv = .true.
  fzrcon = .true.
  fzradb = .true.
!       RSC 23-JAN-2015
  fzrdif = .true.
!       RSC 19-MAY-2015
  fzrcnw = .true.
  fzrdfw = .true.
  kstaw  = kstom4
  kstolu = kstom1
  kstolv = kstom1
  kstolw = kstom1
  kstolt = kstol
  
END IF

!     =========================================================================

!     WALL BC No 2
!     ------------
!     NO-SLIP WALL - ISOTHERMAL

!     X-LEFT
IF(nsbcxl == nsbcw2)THEN
  
  fxlcnv = .true.
!       RSC 23-JAN-2015
  fxldif = .true.
!       RSC 19-MAY-2015
  fxladw = .true.
  fxldfw = .true.
  istow  = istap4
  istalu = istap1
  istalv = istap1
  istalw = istap1
  istale = istap1
  istalt = istal
  
END IF

!     X-RIGHT
IF(nsbcxr == nsbcw2)THEN
  
  fxrcnv = .true.
!       RSC 23-JAN-2015
  fxrdif = .true.
!       RSC 19-MAY-2015
  fxradw = .true.
  fxrdfw = .true.
  istaw  = istom4
  istolu = istom1
  istolv = istom1
  istolw = istom1
  istole = istom1
  istolt = istol
  
END IF

!     Y-LEFT
IF(nsbcyl == nsbcw2)THEN
  
  fylcnv = .true.
!       RSC 23-JAN-2015
  fyrdif = .true.
!       RSC 19-MAY-2015
  fyladw = .true.
  fyldfw = .true.
  jstow  = jstap4
  jstalu = jstap1
  jstalv = jstap1
  jstalw = jstap1
  jstale = jstap1
  jstalt = jstal
  
END IF

!     Y-RIGHT
IF(nsbcyr == nsbcw2)THEN
  
  fyrcnv = .true.
!       RSC 23-JAN-2015
  fyrdif = .true.
!       RSC 19-MAY-2015
  fyradw = .true.
  fyrdfw = .true.
  jstaw  = jstom4
  jstolu = jstom1
  jstolv = jstom1
  jstolw = jstom1
  jstole = jstom1
  jstolt = jstol
  
END IF

!     Z-LEFT
IF(nsbczl == nsbcw2)THEN
  
  fzlcnv = .true.
!       RSC 23-JAN-2015
  fzldif = .true.
!       RSC 19-MAY-2015
  fzladw = .true.
  fzldfw = .true.
  kstow  = kstap4
  kstalu = kstap1
  kstalv = kstap1
  kstalw = kstap1
  kstale = kstap1
  kstalt = kstal
  
END IF

!     Z-RIGHT
IF(nsbczr == nsbcw2)THEN
  
  fzrcnv = .true.
!       RSC 23-JAN-2015
  fzrdif = .true.
!       RSC 19-MAY-2015
  fzradw = .true.
  fzrdfw = .true.
  kstaw  = kstom4
  kstolu = kstom1
  kstolv = kstom1
  kstolw = kstom1
  kstole = kstom1
  kstolt = kstol
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcinit
