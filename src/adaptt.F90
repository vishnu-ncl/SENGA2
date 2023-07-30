SUBROUTINE adaptt

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   ADAPTT
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   19-JAN-2003:  CREATED
!   23-AUG-2009:  RSC REVISE ERROR NORM EVALUATION
!   08-AUG-2012:  RSC EVALUATE ALL SPECIES
!   09-AUG-2012   RSC/RACG USE GLOBAL ERROR

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   COMPUTES NEW TIMESTEP FOR ERK SCHEME

!   *************************************************************************


!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    real(kind=8) :: erytot(nspcmx)
    real(kind=8) :: erdtot,erutot,ervtot,erwtot,eretot
    real(kind=8) :: errmax,tratio,tstold
!   RSC/RACG 09-AUG-2012 USE GLOBAL ERROR
!   real(kind=8) TSTLOC
    integer(kind=4) :: ispec
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   CHECK ADAPTION FLAG
!   -------------------
    IF(fladpt) THEN

!   =======================================================================

!   INITIALISE THE ERROR NORM TOTALS
!   --------------------------------
        erdtot = zero
        erutot = zero
        ervtot = zero
        erwtot = zero
        eretot = zero
!       RSC 08-AUG-2012 EVALUATE ALL SPECIES
!       DO ISPEC = 1,NSPM1
        DO ispec = 1,nspec
            erytot(ispec) = zero
        END DO

!       =======================================================================

!       ERK ERROR EVALUATION
!       --------------------
!       USING ERK ERROR ARRAYS
!       RSC 23 AUG-2009 REVISE ERROR NORM EVALUATION

!       EVALUATE ERROR NORMS
        rangexyz(1) = 1
        IF (nsbcxl == nsbci3) rangexyz(1) = 2
        rangexyz(2) = nxglbl
        IF (nsbcxr == nsbci3) rangexyz(2) = nxglbl-1
        rangexyz(3) = 1
        IF (nsbcyl == nsbci3) rangexyz(3) = 2
        rangexyz(4) = nyglbl
        IF (nsbcyr == nsbci3) rangexyz(4) = nyglbl-1
        rangexyz(5) = 1
        IF (nsbczl == nsbci3) rangexyz(5) = 2
        rangexyz(6) = nzglbl
        IF (nsbczr == nsbci3) rangexyz(6) = nzglbl-1

        call ops_par_loop(adaptt_kernel_err_eval, "EVALUATE ERROR NORMS", senga_grid, 3, rangexyz,  &
                        &  ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_gbl(erdnrm, 1, "real(kind=8)", OPS_READ), &
                        &  ops_arg_reduce(h_erdtot, 1, "real(kind=8)", OPS_MAX))
        call ops_reduction_result(h_erdtot, erdtot)

        rangexyz(1) = 1
        IF (nsbcxl == nsbci2 .or. nsbcxl == nsbci3 .or. nsbcxl == nsbcw1 .or. nsbcxl == nsbcw2) rangexyz(1) = 2
        rangexyz(2) = nxglbl
        IF (nsbcxr == nsbci2 .or. nsbcxr == nsbci3 .or. nsbcxr == nsbcw1 .or. nsbcxr == nsbcw2) rangexyz(2) = nxglbl-1
        rangexyz(3) = 1
        IF (nsbcyl == nsbci2 .or. nsbcyl == nsbci3 .or. nsbcyl == nsbcw1 .or. nsbcyl == nsbcw2) rangexyz(3) = 2
        rangexyz(4) = nyglbl
        IF (nsbcyr == nsbci2 .or. nsbcyr == nsbci3 .or. nsbcyr == nsbcw1 .or. nsbcyr == nsbcw2) rangexyz(4) = nyglbl-1
        rangexyz(5) = 1
        IF (nsbczl == nsbci2 .or. nsbczl == nsbci3 .or. nsbczl == nsbcw1 .or. nsbczl == nsbcw2) rangexyz(5) = 2
        rangexyz(6) = nzglbl
        IF (nsbczr == nsbci2 .or. nsbczr == nsbci3 .or. nsbczr == nsbcw1 .or. nsbczr == nsbcw2) rangexyz(6) = nzglbl-1

        call ops_par_loop(adaptt_kernel_err_eval, "EVALUATE ERROR NORMS", senga_grid, 3, rangexyz,  &
                        &  ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_gbl(erunrm, 1, "real(kind=8)", OPS_READ), &
                        &  ops_arg_reduce(h_erutot, 1, "real(kind=8)", OPS_MAX))
        call ops_reduction_result(h_erutot, erutot)

        call ops_par_loop(adaptt_kernel_err_eval, "EVALUATE ERROR NORMS", senga_grid, 3, rangexyz,  &
                        &  ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_gbl(ervnrm, 1, "real(kind=8)", OPS_READ), &
                        &  ops_arg_reduce(h_ervtot, 1, "real(kind=8)", OPS_MAX))
        call ops_reduction_result(h_ervtot, ervtot)

        call ops_par_loop(adaptt_kernel_err_eval, "EVALUATE ERROR NORMS", senga_grid, 3, rangexyz,  &
                        &  ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_gbl(erwnrm, 1, "real(kind=8)", OPS_READ), &
                        &  ops_arg_reduce(h_erwtot, 1, "real(kind=8)", OPS_MAX))
        call ops_reduction_result(h_erwtot, erwtot)

        rangexyz(1) = 1
        IF (nsbcxl == nsbci2 .or. nsbcxl == nsbcw2) rangexyz(1) = 2
        rangexyz(2) = nxglbl
        IF (nsbcxr == nsbci2 .or. nsbcxr == nsbcw2) rangexyz(2) = nxglbl-1
        rangexyz(3) = 1
        IF (nsbcyl == nsbci2 .or. nsbcyl == nsbcw2) rangexyz(3) = 2
        rangexyz(4) = nyglbl
        IF (nsbcyr == nsbci2 .or. nsbcyr == nsbcw2) rangexyz(4) = nyglbl-1
        rangexyz(5) = 1
        IF (nsbczl == nsbci2 .or. nsbczl == nsbcw2) rangexyz(5) = 2
        rangexyz(6) = nzglbl
        IF (nsbczr == nsbci2 .or. nsbczr == nsbcw2) rangexyz(6) = nzglbl-1

        call ops_par_loop(adaptt_kernel_err_eval, "EVALUATE ERROR NORMS", senga_grid, 3, rangexyz,  &
                        &  ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_gbl(erenrm, 1, "real(kind=8)", OPS_READ), &
                        &  ops_arg_reduce(h_eretot, 1, "real(kind=8)", OPS_MAX))
        call ops_reduction_result(h_eretot, eretot)

!       RSC 08-AUG-2012 EVALUATE ALL SPECIES
!       DO ISPEC = 1,NSPM1
        rangexyz(1) = 1
        IF (nsbcxl == nsbci2 .or. nsbcxl == nsbci3) rangexyz(1) = 2
        rangexyz(2) = nxglbl
        IF (nsbcxr == nsbci2 .or. nsbcxr == nsbci3) rangexyz(2) = nxglbl-1
        rangexyz(3) = 1
        IF (nsbcyl == nsbci2 .or. nsbcyl == nsbci3) rangexyz(3) = 2
        rangexyz(4) = nyglbl
        IF (nsbcyr == nsbci2 .or. nsbcyr == nsbci3) rangexyz(4) = nyglbl-1
        rangexyz(5) = 1
        IF (nsbczl == nsbci2 .or. nsbczl == nsbci3) rangexyz(5) = 2
        rangexyz(6) = nzglbl
        IF (nsbczr == nsbci2 .or. nsbczr == nsbci3) rangexyz(6) = nzglbl-1

        DO ispec = 1,nspec

            call ops_par_loop(adaptt_kernel_err_eval_MD, "EVALUATE ERROR NORMS - MULTIDIM", senga_grid, 3, rangexyz,  &
                        &  ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                        &  ops_arg_gbl(erynrm, nspcmx, "real(kind=8)", OPS_READ), &
                        &  ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                        &  ops_arg_reduce(h_erytot, 1, "real(kind=8)", OPS_MAX))
            call ops_reduction_result(h_erytot, erytot(ispec))

        END DO

!       =======================================================================

!       FIND THE MAXIMUM
!       ----------------
        errmax = zero
        IF(erdtot > errmax)THEN
            errmax = erdtot
            inderr = -4
        END IF

        IF(erutot > errmax)THEN
            errmax = erutot
            inderr = -1
            END IF

        IF(ervtot > errmax)THEN
            errmax = ervtot
            inderr = -2
        END IF

        IF(erwtot > errmax)THEN
            errmax = erwtot
            inderr = -3
        END IF

        IF(eretot > errmax)THEN
            errmax = eretot
            inderr = 0
        END IF

!       RSC 08-AUG-2012 EVALUATE ALL SPECIES
!       DO ISPEC = 1,NSPM1
        DO ispec = 1,nspec
            IF(erytot(ispec) > errmax)THEN
                errmax = erytot(ispec)
                inderr = ispec
            END IF
        END DO

!       =======================================================================

        IF(iproc == 0) THEN
            write(*,'(a,1PE12.4)') "MPI_MAX error: ",errmax
        END IF

!       =======================================================================

!       EVALUATE THE NEW TIME STEP
!       --------------------------
!       ZERO CHECK
        IF(errmax < errlow)errmax = errlow

!       -----------------------------------------------------------------------

!       -----------------------------------------------------------------------

!       PID-CONTROLLER
        tratio = ctmult*EXP(ctalph*LOG(errtol/errmax) +ctbeta*LOG(errold/errtol)  &
                +ctgama*LOG(errtol/errldr))
        errldr = errold
        errold = errmax

!       -----------------------------------------------------------------------

!       LIMIT CHANGES TO TIME STEP
        IF(tratio > trmax) tratio = trmax
        IF(tratio < trmin) tratio = trmin

!       -----------------------------------------------------------------------

!       SAVE THE OLD TIME STEP
        tstold = tstep

!       SET THE NEW TIME STEP
        tstep = tstep*tratio

!       LIMIT THE TIME STEP
        IF(tstep > tsmax) tstep = tsmax
        IF(tstep < tsmin) tstep = tsmin

!       =======================================================================

!       PARALLEL TRANSFER TO SET NEW GLOBAL TIME STEP
!       ---------------------------------------------
!       NEW TIME STEP IS THE GLOBAL MINIMUM OVER ALL PROCESSORS
!       RSC/RACG 09-AUG-2012 USE GLOBAL ERROR
!       TSTLOC = TSTEP
!       CALL P_GMIN(TSTLOC,TSTEP)

!       =======================================================================

!       UPDATE THE TIME ADVANCEMENT COEFFICIENTS
!       AND THE RK SUBSTEP TIME LEVELS
        tratio = tstep/tstold
        DO irkstp = 1, nrkstp
            rklhs(irkstp) = rklhs(irkstp)*tratio
            rkrhs(irkstp) = rkrhs(irkstp)*tratio
            rkerr(irkstp) = rkerr(irkstp)*tratio
            rktim(irkstp) = rktim(irkstp)*tratio
        END DO

!       =======================================================================

!       (RE)INITIALISE ERK ERROR ARRAYS
!       -------------------------------
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!       RSC 08-AUG-2012 EVALUATE ALL SPECIES
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        DO ispec = 1,nspec
            call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))

        END DO

!       =======================================================================

!       (RE)INITIALISE ERK SUBSTEP ERROR NORMS
!       --------------------------------------
        DO irkstp = 1,nrkstp

            erdrhs(irkstp) = zero
            erurhs(irkstp) = zero
            ervrhs(irkstp) = zero
            erwrhs(irkstp) = zero
            ererhs(irkstp) = zero
!           RSC 08-AUG-2012 EVALUATE ALL SPECIES
!           DO ISPEC = 1,NSPM1
            DO ispec = 1,nspec
                eryrhs(ispec,irkstp) = zero
            END DO

        END DO

!   =======================================================================

    END IF
!   ADAPTION FLAG

!   =========================================================================

END SUBROUTINE adaptt
