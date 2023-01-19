SUBROUTINE bcytxl

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga
 
!   *************************************************************************

!   BCYTXL
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   30-DEC-2003:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
!   AND THEIR TIME DERIVATIVES

!   X-DIRECTION LEFT-HAND END

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    integer :: ispec
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!   =========================================================================

!   EVALUATE AND RETURN STRYXL,DYDTXL
    DO ispec = 1,nspec
        rangexyz = (/1,1,1,nysize,1,nzsize/)
        call ops_par_loop(bcyt_kernel_xdir, "bcyt_kernel_xdir", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stryxl, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                            ops_arg_dat(d_dydtxl, 9, s3d_000_strid3d_yz, "real(8)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                            ops_arg_gbl(yrin(ispec), 1, "real(8)", OPS_READ))

    END DO

!   =========================================================================

END SUBROUTINE bcytxl
