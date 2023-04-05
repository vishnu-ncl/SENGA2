SUBROUTINE temper
 
    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   TEMPER
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   16-NOV-2002:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   COMPUTES TEMPERATURE AND PRESSURE

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

    integer :: rangexyz(6)

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

    call ops_par_loop(temper_kernel_main, "temper kernel", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_dat(d_yrhs, 2, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_gbl(amascp, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                    ops_arg_gbl(amasct, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                    ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                    ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                    ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                    ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                    ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                    ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ), &
                    ops_arg_gbl(nspec, 1, "integer", OPS_READ), &
                    ops_arg_idx())

!   =========================================================================

END SUBROUTINE temper
