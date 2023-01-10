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

    rangexyz = (/istalt,istolt,jstalt,jstolt,kstalt,kstolt/)
    call ops_par_loop(temper_kernel_main, "temper kernel", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_WRITE), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_gbl(amascp, 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(amasct, 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(ncenth, 1, "integer", OPS_READ), &
                    ops_arg_gbl(ncpom1, 1, "integer", OPS_READ), &
                    ops_arg_gbl(ncpoly, 1, "integer", OPS_READ), &
                    ops_arg_gbl(tinthi, 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rgspec, 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(ntint, 1, "integer", OPS_READ), &
                    ops_arg_gbl(nctmax, 1, "integer", OPS_READ), &
                    ops_arg_gbl(ncofmx, 1, "integer", OPS_READ), &
                    ops_arg_gbl(ntinmx, 1, "integer", OPS_READ), &
                    ops_arg_gbl(nspcmx, 1, "integer", OPS_READ), &
                    ops_arg_gbl(nintmx, 1, "integer", OPS_READ), &
                    ops_arg_gbl(nspec, 1, "integer", OPS_READ), &
                    ops_arg_gbl(nspimx, 1, "integer", OPS_READ), &
                    ops_arg_gbl(ntbase, 1, "integer", OPS_READ), &
                    ops_arg_gbl(nctmm1, 1, "integer", OPS_READ), &
                    ops_arg_idx())

!   =========================================================================

END SUBROUTINE temper
