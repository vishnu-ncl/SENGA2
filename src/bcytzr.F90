SUBROUTINE bcytzr
 
    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCYTZR
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   26-OCT-2013:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
!   AND THEIR TIME DERIVATIVES

!   Z-DIRECTION RIGHT-HAND END

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

!   EVALUATE AND RETURN STRYZR,DYDTZR
    DO ispec = 1,nspec
        rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
        call ops_par_loop(bcyt_kernel_zdir, "bcyt_kernel_zdir", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stryzr, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                            ops_arg_dat(d_dydtzr, 9, s3d_000_strid3d_xy, "real(8)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                            ops_arg_gbl(yrin(ispec), 1, "real(8)", OPS_READ))

    END DO

!   =========================================================================

END SUBROUTINE bcytzr
