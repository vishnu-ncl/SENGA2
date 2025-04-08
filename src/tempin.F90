SUBROUTINE tempin

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   TEMPIN
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   15-MAY-2004:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   INITIALISES TEMPERATURE AND PRESSURE
!   USES BISECTION METHOD FOR ROBUSTNESS

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   PARAMETERS
!   ==========
    real(kind=8), parameter :: toltmp=0.00010_8
    real(kind=8), parameter :: tininc=50.0_8
    real(kind=8), parameter :: tlimlo=200.0_8
    real(kind=8), parameter :: tlimhi=3000.0_8

!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)
    real(kind=8) :: tempor,tupper,tlower
    integer(kind=4) :: iindex,ipower,ispec
    logical :: fnconv

!   BEGIN
!   =====

!   =========================================================================

!   TEMPERATURE AND PRESSURE
!   ------------------------

!   TEMPERATURE AND PRESSURE ARE PARALLEL

    rangexyz(1) = 1-nhalox
    IF (nsbcxl == nsbco1 .or. nsbcxl == nsbci1 .or. nsbcxl == nsbci2 .or. &
        nsbcxl == nsbci3 .or. nsbcxl == nsbcw1 .or. nsbcxl == nsbcw2) rangexyz(1) = 1

    rangexyz(2) = nxglbl+nhalox
    IF (nsbcxr == nsbco1 .or. nsbcxr == nsbci1 .or. nsbcxr == nsbci2 .or. &
        nsbcxr == nsbci3 .or. nsbcxr == nsbcw1 .or. nsbcxr == nsbcw2) rangexyz(2) = nxglbl

    rangexyz(3) = 1-nhaloy
    IF (nsbcyl == nsbco1 .or. nsbcyl == nsbci1 .or. nsbcyl == nsbci2 .or. &
        nsbcyl == nsbci3 .or. nsbcyl == nsbcw1 .or. nsbcyl == nsbcw2) rangexyz(3) = 1

    rangexyz(4) = nyglbl+nhaloy
    IF (nsbcyr == nsbco1 .or. nsbcyr == nsbci1 .or. nsbcyr == nsbci2 .or. &
        nsbcyr == nsbci3 .or. nsbcyr == nsbcw1 .or. nsbcyr == nsbcw2) rangexyz(4) = nyglbl

    rangexyz(5) = 1-nhaloz
    IF (nsbczl == nsbco1 .or. nsbczl == nsbci1 .or. nsbczl == nsbci2 .or. &
        nsbczl == nsbci3 .or. nsbczl == nsbcw1 .or. nsbczl == nsbcw2) rangexyz(5) = 1

    rangexyz(6) = nzglbl+nhaloz
    IF (nsbczr == nsbco1 .or. nsbczr == nsbci1 .or. nsbczr == nsbci2 .or. &
        nsbczr == nsbci3 .or. nsbczr == nsbcw1 .or. nsbczr == nsbcw2) rangexyz(6) = nzglbl

    DO ispec = 1,nspcmx
        call ops_par_loop(copy_kernel_sdim_to_mdim, "A_multidim(ispec) = B", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrhs_mdim, 9, s3d_000, "real(kind=8)", OPS_RW), &
                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
    END DO

    call ops_par_loop(tempin_kernel_main, "tempin kernel", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_yrhs_mdim, 9, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(amascp, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(amasct, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                    ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                    ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                    ops_arg_gbl(trin, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(nspec, 1, "integer(kind=4)", OPS_READ), &
                    ops_arg_gbl(iproc, 1, "integer(kind=4)", OPS_READ), &
                    ops_arg_idx())

    DO ispec = 1,nspcmx
        call ops_par_loop(copy_kernel_mdim_to_sdim, "A = B_multidim(ispec)", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs_mdim, 9, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
    END DO

!   CONSTRUCT THE TEMPERATURE INTERVAL INDEX
!   EVALUATE PRESSURE
!   EVALUATE MIXTURE SPECIFIC HEAT CP
    DO iindex = 1,nintmx
        call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
    END DO

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    DO ispec = 1, nspec
!       SET THE TEMPERATURE INTERVAL INDEX
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1

        call ops_par_loop(temper_kernel_eqE, "temper eq E", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC), &
                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(amascp, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ))

!       EVALUATE (DENSITY TIMES) MIXTURE GAS CONSTANT FOR PRESSURE
        call ops_par_loop(temper_kernel_eqC, "temper eq C", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END DO

!   =========================================================================

    call ops_par_loop(temper_kernel_eqF, "temper eq F", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

!   =========================================================================

!    call ops_free_dat(d_yrhs_mdim)

END SUBROUTINE tempin
