SUBROUTINE bcytxr

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCYTXR
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

!   X-DIRECTION RIGHT-HAND END

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

!   EVALUATE AND RETURN STRYXR,DYDTXR
    DO ispec = 1,nspec
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(bcyt_kernel_xdir_eqA, "bcyt_kernel_xdir_eqA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dydtxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(yrin, nspcmx, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END DO

!   =========================================================================

END SUBROUTINE bcytxr
