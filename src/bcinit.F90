SUBROUTINE bcinit

!   *************************************************************************

!   BCINIT
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   26-MAR-2003:  CREATED
!   26-OCT-2008:  RSC/TDD BUG FIX BCTIXL
!   23-JAN-2015:  RSC BUG FIX WALL MASS FLUX BC
!   19-MAY-2015:  RSC UPDATE WALL DIFFUSIVE BC

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   CARRIES OUT INITIALISATION OF BOUNDARY CONDITIONS

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
    use com_senga
!   -------------------------------------------------------------------------

!   BEGIN
!   =====

!   =========================================================================

!   APPLY GLOBAL BOUNDARY CONDITIONS TO LOCAL PROCESSOR
!   ===================================================

!   X-LEFT
    nsbcxl = ngbcxl
    nendxl = nbound
    IF(nsbcxl == nsperi)THEN
        nendxl = nperi
    END IF

!   X-RIGHT
    nsbcxr = ngbcxr
    nendxr = nbound
    IF(nsbcxr == nsperi)THEN
        nendxr = nperi
    END IF

!   Y-LEFT
    nsbcyl = ngbcyl
    nendyl = nbound
    IF(nsbcyl == nsperi)THEN
        nendyl = nperi
    END IF

!   Y-RIGHT
    nsbcyr = ngbcyr
    nendyr = nbound
    IF(nsbcyr == nsperi)THEN
        nendyr = nperi
    END IF

!   Z-LEFT
    nsbczl = ngbczl
    nendzl = nbound
    IF(nsbczl == nsperi)THEN
        nendzl = nperi
    END IF

!   Z-RIGHT
    nsbczr = ngbczr
    nendzr = nbound
    IF(nsbczr == nsperi)THEN
        nendzr = nperi
    END IF

    IF (iproc == 0) THEN
        print *, "nendxl:",nendxl, " nendxr:",nendxr
        print *, "nendyl:",nendyl, " nendyr:",nendyr
        print *, "nendzl:",nendzl, " nendzr:",nendzr
    END IF

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
END IF

!     X-RIGHT
IF(nsbcxr == nsbci1)THEN
  fxrcnv = .true.
  fxrvsn = .true.
END IF

!     Y-LEFT
IF(nsbcyl == nsbci1)THEN
  fylcnv = .true.
  fylvsn = .true.
END IF

!     Y-RIGHT
IF(nsbcyr == nsbci1)THEN
  fyrcnv = .true.
  fyrvsn = .true.
END IF

!     Z-LEFT
IF(nsbczl == nsbci1)THEN
  fzlcnv = .true.
  fzlvsn = .true.
END IF

!     Z-RIGHT
IF(nsbczr == nsbci1)THEN
  fzrcnv = .true.
  fzrvsn = .true.
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
  fxltrb = .false.

END IF

!     X-RIGHT
IF(nsbcxr == nsbci2)THEN

  fxrcnv = .true.

END IF

!     Y-LEFT
IF(nsbcyl == nsbci2)THEN

  fylcnv = .true.

END IF

!     Y-RIGHT
IF(nsbcyr == nsbci2)THEN

  fyrcnv = .true.

END IF

!     Z-LEFT
IF(nsbczl == nsbci2)THEN

  fzlcnv = .true.

END IF

!     Z-RIGHT
IF(nsbczr == nsbci2)THEN

  fzrcnv = .true.

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
  fxltrb = .false.

END IF

!     X-RIGHT
IF(nsbcxr == nsbci3)THEN

  fxrcnv = .true.
  fxrvsn = .true.

END IF

!     Y-LEFT
IF(nsbcyl == nsbci3)THEN

  fylcnv = .true.
  fylvsn = .true.

END IF

!     Y-RIGHT
IF(nsbcyr == nsbci3)THEN

  fyrcnv = .true.
  fyrvsn = .true.

END IF

!     Z-LEFT
IF(nsbczl == nsbci3)THEN

  fzlcnv = .true.
  fzlvsn = .true.

END IF

!     Z-RIGHT
IF(nsbczr == nsbci3)THEN

  fzrcnv = .true.
  fzrvsn = .true.

END IF

!     =========================================================================

!     INFLOW BC no 4
!     --------------
!     SUBSONIC NON-REFLECTING INFLOW
!     FY - IMPLEMENTATION AS PER LODATO
!     SPECIFY TEMPERATURE, VELOCITY AND MASS FRACTIONS

!     X-LEFT
IF(nsbcxl == nsbci4)THEN
  fxlcnv = .true.
  fxlvsn = .true.
  fxltrb = .false.
END IF

!     X-RIGHT
IF(nsbcxr == nsbci4)THEN
  fxrcnv = .true.
  fxrvsn = .true.
END IF

!     Y-LEFT
IF(nsbcyl == nsbci4)THEN
  fylcnv = .true.
  fylvsn = .true.
END IF

!     Y-RIGHT
IF(nsbcyr == nsbci4)THEN
  fyrcnv = .true.
  fyrvsn = .true.
END IF

!     Z-LEFT
IF(nsbczl == nsbci4)THEN
  fzlcnv = .true.
  fzlvsn = .true.
END IF

!     Z-RIGHT
IF(nsbczr == nsbci4)THEN
  fzrcnv = .true.
  fzrvsn = .true.
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

END IF

!     X-RIGHT
IF(nsbcxr == nsbcw2)THEN

  fxrcnv = .true.
!       RSC 23-JAN-2015
  fxrdif = .true.
!       RSC 19-MAY-2015
  fxradw = .true.
  fxrdfw = .true.

END IF

!     Y-LEFT
IF(nsbcyl == nsbcw2)THEN

  fylcnv = .true.
!       RSC 23-JAN-2015
  fyldif = .true.
!       RSC 19-MAY-2015
  fyladw = .true.
  fyldfw = .true.

END IF

!     Y-RIGHT
IF(nsbcyr == nsbcw2)THEN

  fyrcnv = .true.
!       RSC 23-JAN-2015
  fyrdif = .true.
!       RSC 19-MAY-2015
  fyradw = .true.
  fyrdfw = .true.

END IF

!     Z-LEFT
IF(nsbczl == nsbcw2)THEN

  fzlcnv = .true.
!       RSC 23-JAN-2015
  fzldif = .true.
!       RSC 19-MAY-2015
  fzladw = .true.
  fzldfw = .true.

END IF

!     Z-RIGHT
IF(nsbczr == nsbcw2)THEN

  fzrcnv = .true.
!       RSC 23-JAN-2015
  fzrdif = .true.
!       RSC 19-MAY-2015
  fzradw = .true.
  fzrdfw = .true.

END IF

!     =========================================================================


RETURN
END SUBROUTINE bcinit
