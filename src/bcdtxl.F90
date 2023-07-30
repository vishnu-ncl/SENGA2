SUBROUTINE bcdtxl

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCDTXL
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
!   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR DENSITY
!   AND ITS TIME DERIVATIVE

!   X-DIRECTION LEFT-HAND END

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!   =========================================================================

!   EVALUATE AND RETURN STRDXL,DDDTXL
    rangexyz = [1,1,1,nyglbl,1,nzglbl]
    call ops_par_loop(bcdt_kernel_xdir, "bcdt_kernel_xdir", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                    ops_arg_dat(d_dddtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(drin, 1, "real(kind=8)", OPS_READ))

!   =========================================================================

END SUBROUTINE bcdtxl
