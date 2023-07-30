SUBROUTINE bcytzl

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCYTZL
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

!   Z-DIRECTION LEFT-HAND END

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    integer(kind=4) :: ispec
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!   =========================================================================

!   EVALUATE AND RETURN STRYZL,DYDTZL
    DO ispec = 1,nspec
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(bcyt_kernel_zdir, "bcyt_kernel_zdir", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dydtzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(yrin, nspcmx, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END DO

!   =========================================================================

END SUBROUTINE bcytzl
