SUBROUTINE bcytyr

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCYTYR
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

!   Y-DIRECTION RIGHT-HAND END

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

!   EVALUATE AND RETURN STRYYR,DYDTYR
    DO ispec = 1,nspec
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(bcyt_kernel_ydir, "bcyt_kernel_ydir", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dydtyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(yrin, nspcmx, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
    END DO

!   =========================================================================

END SUBROUTINE bcytyr
