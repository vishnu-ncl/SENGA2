SUBROUTINE espini

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_espect
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   ESPINI
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  DEPARTMENT OF MECHANICAL ENGINEERING, UMIST

!   CHANGE RECORD
!   -------------
!   11-APR-2004:  CREATED
!   09-JUN-2008:  RSC BUG FIX TURBULENCE INTEGRAL SCALING
!   09-JUN-2008:  RSC ENHANCEMENTS TO LENGTH SCALE AND DISSIPATION DATA

!   DESCRIPTION
!   -----------
!   INITIALISES CONSTANTS FOR THE TURBULENT ENERGY SPECTRUM FUNCTION ESPECT

!   REFERENCES
!   ----------
!   1) BATCHELOR, G.K., TOWNSEND, A.A.: JFM 88(4), 685-709, 1948.

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   PARAMETERS
!   ==========
!   RSC 09-JUN-2008 ENHANCEMENTS
    real(kind=8), parameter :: ten = 10.0_8

!   FUNCTIONS
!   =========
    real(kind=8) :: espect
    EXTERNAL espect,espovk,espksq

!   LOCAL DATA
!   ==========
    real(kind=8) :: tenrgy,tinscl,tdisip
    real(kind=8) :: speclo,spechi
    real(kind=8) :: argmnt
!   RSC 09-JUN-2008 ENHANCEMENTS
    real(kind=8) :: viscos,viskin,tkolmo,taylor,rtin,dvrin
    integer(kind=4) :: imodes,nmodes
!   RSC 09-JUN-2008 ENHANCEMENTS
    integer(kind=4) :: ispec

!   BEGIN
!   =====
!   =========================================================================

!   SET THE SPECTRUM PARAMETERS
!   ---------------------------
    const0 = sparam(1)
    ck0 = sparam(2)
    ovk0 = one/ck0
    covk0 = const0/ck0

!   =========================================================================

!   CHECK THE ENERGY SPECTRUM
!   -------------------------
    IF(iproc == 0)THEN

        WRITE(ncrept,*)
        WRITE(ncrept,*)'Energy spectrum function data'
        WRITE(ncrept,*)'-----------------------------'
        WRITE(ncrept,*)

!       SET MAX NO OF MODES FOR SPECTRUM EVALUATION
        nmodes = ngzmax/2

!       EVALUATE THE SPECTRUM
        WRITE(ncrept,*)'spectrum function:'
        WRITE(ncrept,*)'mode       k          E(k)'
        DO imodes = 0, nmodes
!           RSC 09-JUN-2008 BUG FIX
!           ARGMNT = TWO*PI*REAL(IMODES)
            argmnt = REAL(imodes,kind=8)
            WRITE(ncrept,'(I5,2(1PE12.4))')imodes,argmnt,espect(argmnt)
        END DO
        WRITE(ncrept,*)

!       INTEGRATE THE SPECTRUM
!       LIMITS
        speclo = zero
!       RSC 09-JUN-2008 BUG FIX
!       PECHI = TWO*PI*REAL(NMODES)
        spechi = REAL(nmodes,kind=8)

!       INTEGRATE FOR THE ENERGY
        call integf(espect,speclo,spechi,tenrgy)
!       INTEGRATE FOR THE LENGTH SCALE
        call integf(espovk,speclo,spechi,tinscl)
!       INTEGRATE FOR THE DISSIPATION
        call integf(espksq,speclo,spechi,tdisip)

!       WRITE OUT THE RESULTS
!       WRITE OUT THE RAW INTEGRALS
!       APPLIES TO SPHERICAL UNIT VOLUME IN K-SPACE
        WRITE(ncrept,*)'raw integrals:'
        WRITE(ncrept,*)'energy, length scale, dissipation'
        WRITE(ncrept,'(3(1PE12.4))')tenrgy,tinscl,tdisip
        WRITE(ncrept,*)

!       CONVERT TO LENGTH SCALE AND DISSIPATION
!       RSC 09-JUN-2008 BUG FIX
!       TINSCL = TINSCL*THREE*PI/TENRGY/TWO/TWO
!       TDISIP = TWO*TDISIP
        tinscl = tinscl*three/tenrgy/two/two/two
        tdisip = two*two*two*pi*pi*tdisip
        WRITE(ncrept,*)'scaled values:'
        WRITE(ncrept,*)'energy, length scale, dissipation/viscosity:'
        WRITE(ncrept,'(3(1PE12.4))')tenrgy,tinscl,tdisip
        WRITE(ncrept,*)

!       RSC 09-JUN-2008 ENHANCEMENTS
!       CONVERT TO DIMENSIONAL LENGTH SCALE AND DISSIPATION/VISCOSITY
        tinscl = tinscl*xgdlen
        tdisip = tdisip/(xgdlen*xgdlen)
        WRITE(ncrept,*)'dimensional values:'
        WRITE(ncrept,*)'energy, length scale, dissipation/viscosity:'
        WRITE(ncrept,'(3(1PE12.4))')tenrgy,tinscl,tdisip
        WRITE(ncrept,*)

!       RSC 09-JUN-2008 ENHANCEMENTS
!       COMPUTE VISCOSITIES
        viscos = alamda*EXP(rlamda*LOG(trin))
        viscos = viscos*prantl
!       REFERENCE DENSITY
        rtin = zero
        DO ispec = 1,nspec
            rtin = rtin + rgspec(ispec)*yrin(ispec)
        END DO
        rtin = rtin*trin
        dvrin = prin/rtin
        viskin = viscos/dvrin

        WRITE(ncrept,*)'viscosities: dynamic, kinematic:'
        WRITE(ncrept,'(2(1PE12.4))')viscos,viskin
        WRITE(ncrept,*)
!       RSC 09-JUN-2008 ENHANCEMENTS
!       CONVERT TO DIMENSIONAL DISSIPATION
        tdisip = tdisip*viskin
        WRITE(ncrept,*)'energy, length scale, dissipation:'
        WRITE(ncrept,'(3(1PE12.4))')tenrgy,tinscl,tdisip
        WRITE(ncrept,*)

!       RSC 09-JUN-2008 ENHANCEMENTS
!       COMPUTE TAYLOR LENGTH SCALE
        taylor = ten*viskin*tenrgy/tdisip
        taylor = SQRT(taylor)
!       COMPUTE KOLMOGOROV LENGTH SCALE
        tkolmo = viskin*viskin*viskin/tdisip
        tkolmo = SQRT(SQRT(tkolmo))
        WRITE(ncrept,*)'length scales: int, Taylor, Kolmo:'
        WRITE(ncrept,'(3(1PE12.4))')tinscl,taylor,tkolmo
        WRITE(ncrept,*)

        WRITE(ncrept,*)'domain/length scale ratios: int, Taylor, Kolmo:'
        WRITE(ncrept,'(3(1PE12.4))')xgdlen/tinscl,xgdlen/taylor, xgdlen/tkolmo
        WRITE(ncrept,*)

!       RSC 09-JUN-2008 ENHANCEMENTS
        WRITE(ncrept,*)'------------------------------------'
        WRITE(ncrept,*)

    END IF

!   =========================================================================

END SUBROUTINE espini
