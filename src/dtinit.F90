SUBROUTINE dtinit

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   DTINIT
!   ======

!   AUTHOR
!   ------
!   R.S.CANT -- CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   18-MAY-2003:  CREATED
!   07-JUL-2009:  RSC BUG FIX ERROR NORMS

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   INITIALISES TIME-STEPPING FOR ERK SCHEME

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    real(kind=8) :: acofrk(nrkmax),bcofrk(nrkmax),bhatrk(nrkmax)
    integer(kind=4) :: ic,jc,kc
    integer(kind=4) :: ispec
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   VALUES DEFINED IN THIS BLOCK FOR

!   ERK COEFFICIENTS
!   ERK ERROR ESTIMATION NORMALISING VALUES
!   ERK TIME STEP CONTROLLER COEFFICIENTS
!   ERK TOLERANCES AND LIMITS

!   -------------------------------------------------------------------------
!
!   ERK SCHEME RK3(2)4[2R+]C
!   ------------------------
!   NRKSTP = 4

!   ACOFRK(1) = 11847461282814D0/36547543011857D0
!   ACOFRK(2) = 3943225443063D0/7078155732230D0
!   ACOFRK(3) = -346793006927D0/4029903576067D0
!   ACOFRK(4) = ZERO

!   BCOFRK(1) = 1017324711453D0/9774461848756D0
!   BCOFRK(2) = 8237718856693D0/13685301971492D0
!   BCOFRK(3) = 57731312506979D0/19404895981398D0
!   BCOFRK(4) = -101169746363290D0/37734290219643D0

!   BHATRK(1) = 15763415370699D0/46270243929542D0
!   BHATRK(2) = 514528521746D0/5659431552419D0
!   BHATRK(3) = 27030193851939D0/9429696342944D0
!   BHATRK(4) = -69544964788955D0/30262026368149D0

!   -------------------------------------------------------------------------

!   ERK SCHEME RK4(3)5[2R+]C
!   ------------------------
    nrkstp = 5

    acofrk(1) = 970286171893.0_8/4311952581923.0_8
    acofrk(2) = 6584761158862.0_8/12103376702013.0_8
    acofrk(3) = 2251764453980.0_8/15575788980749.0_8
    acofrk(4) = 26877169314380.0_8/34165994151039.0_8
    acofrk(5) = zero

    bcofrk(1) = 1153189308089.0_8/22510343858157.0_8
    bcofrk(2) = 1772645290293.0_8/4653164025191.0_8
    bcofrk(3) = -1672844663538.0_8/4480602732383.0_8
    bcofrk(4) = 2114624349019.0_8/3568978502595.0_8
    bcofrk(5) = 5198255086312.0_8/14908931495163.0_8

    bhatrk(1) = 1016888040809.0_8/7410784769900.0_8
    bhatrk(2) = 11231460423587.0_8/58533540763752.0_8
    bhatrk(3) = -1563879915014.0_8/6823010717585.0_8
    bhatrk(4) = 606302364029.0_8/971179775848.0_8
    bhatrk(5) = 1097981568119.0_8/3980877426909.0_8

!   -------------------------------------------------------------------------

!   SET TOLERANCES AND LIMITS FOR THE RK TIME STEP CONTROLLER
!   ---------------------------------------------------------
!   RSC 07-JUL-2009 BUG FIX ERROR NORMS
!   ERRTOL = 1.0D-3
    errtol = 0.00010_8
    errlow = 0.0000000000000000000000000000010_8
    IF(ncdmpi == 0) THEN
        errold = errtol
        errldr = errtol
    END IF
    trmax = 1.01_8
    trmin = 0.010_8
    tsmax = 1.0_8
    tsmin = 0.0000000000000010_8

!   -------------------------------------------------------------------------

!   INITIALISE THE NORMALISING VALUES FOR RK ERROR ESTIMATION
!   ---------------------------------------------------------
    erdnrm = drin
    erunrm = one
    IF(ABS(urin) > errtol) erunrm = drin*ABS(urin)
    ervnrm = one
    IF(ABS(vrin) > errtol) ervnrm = drin*ABS(vrin)
    erwnrm = one
    IF(ABS(wrin) > errtol) erwnrm = drin*ABS(wrin)
    erenrm = one
    IF(ABS(erin) > errtol) erenrm = drin*ABS(erin)
    DO ispec = 1,nspec
        erynrm(ispec) = one
    END DO

!   RSC 07-JUL-2009 BUG FIX ERROR NORMS
    erdnrm = one/erdnrm
    erunrm = one/erunrm
    ervnrm = one/ervnrm
    erwnrm = one/erwnrm
    erenrm = one/erenrm

!   RSC 23-AUG-2009 REFORMULATE ERROR NORMS
    erdnrm = zero
    erunrm = 0.0000010_8
    ervnrm = 0.0000010_8
    erwnrm = 0.0000010_8
    erenrm = 0.010_8
    DO ispec = 1,nspec
        erynrm(ispec) = 0.00000000010_8
    END DO

!   -------------------------------------------------------------------------

!   SET THE COEFFICIENTS FOR THE RK TIME STEP CONTROLLER
!   ----------------------------------------------------
!   I CONTROLLER
!   CTALPH = 1.0/ORDER OF EMBEDDED SCHEME+1
    ctmult = 0.90_8
    ctalph = one/four

!   PI CONTROLLER
!   CTALPH = 0.7/ORDER OF EMBEDDED SCHEME
!   CTBETA = 0.4/ORDER OF EMBEDDED SCHEME
    ctmult = 0.90_8
    ctalph = 0.70_8/three
    ctbeta = 0.40_8/three

!   PID CONTROLLER
!   CTALPH = 0.49/ORDER OF EMBEDDED SCHEME
!   CTBETA = 0.34/ORDER OF EMBEDDED SCHEME
!   CTGAMA = 0.10/ORDER OF EMBEDDED SCHEME
    ctmult = 0.90_8
    ctalph = 0.49_8/three
    ctbeta = 0.34_8/three
    ctgama = 0.10_8/three

!   -------------------------------------------------------------------------

!   SET TIME STEP ADAPTION FLAG
!   ---------------------------
    fladpt = .true.
    IF(nstpsw == 0) fladpt = .false.

!   -------------------------------------------------------------------------

!   END OF BLOCK FOR

!   ERK COEFFICIENTS
!   ERK ERROR ESTIMATION NORMALISING VALUES
!   ERK TIME STEP CONTROLLER COEFFICIENTS
!   ERK TOLERANCES AND LIMITS

!   =========================================================================

!   NO OF SUBSTEPS MINUS ONE
!   ------------------------
    nrksm1 = nrkstp - 1

!   =========================================================================

!   RK COEFFICIENTS
!   ---------------
    DO irkstp = 1,nrkstp
        rklhs(irkstp) = bcofrk(irkstp)
        rkrhs(irkstp) = acofrk(irkstp)
        rkerr(irkstp) = bcofrk(irkstp) - bhatrk(irkstp)
    END DO

!   =========================================================================

!   TIME-ADVANCEMENT COEFFICIENTS
!   -----------------------------
    DO irkstp = 1,nrkstp
        rklhs(irkstp) = rklhs(irkstp)*tstep
        rkrhs(irkstp) = rkrhs(irkstp)*tstep
        rkerr(irkstp) = rkerr(irkstp)*tstep
    END DO

!   =========================================================================

!   TIME LEVELS FOR EACH RK SUBSTEP
!   -------------------------------
    rktim(1) = zero
    DO irkstp = 2,nrkstp
        rktim(irkstp) = acofrk(irkstp-1)
        DO jc = 1,irkstp-2
            rktim(irkstp) = rktim(irkstp) + bcofrk(jc)
        END DO
        rktim(irkstp) = rktim(irkstp)*tstep
    END DO

!   =========================================================================

!   INITIALISE ERK ERROR ARRAYS
!   ---------------------------
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    DO ispec = 1,nspec
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))

    END DO

!   =========================================================================

!   INITIALISE ERK SUBSTEP ERROR NORMS
!   ----------------------------------
    DO irkstp = 1,nrkstp

        erdrhs(irkstp) = zero
        erurhs(irkstp) = zero
        ervrhs(irkstp) = zero
        erwrhs(irkstp) = zero
        ererhs(irkstp) = zero
        DO ispec = 1,nspec
            eryrhs(ispec,irkstp) = zero
        END DO

    END DO

!   =========================================================================

END SUBROUTINE dtinit
