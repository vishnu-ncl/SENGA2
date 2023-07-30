SUBROUTINE bcdtyr

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCDTYR
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
!   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR DENSITY
!   AND ITS TIME DERIVATIVE

!   Y-DIRECTION RIGHT-HAND END

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

!   EVALUATE AND RETURN STRDYR,DDDTYR
    rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
    call ops_par_loop(bcdt_kernel_ydir, "bcdt_kernel_ydir", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                    ops_arg_dat(d_dddtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(drin, 1, "real(kind=8)", OPS_READ))

!   =========================================================================

END SUBROUTINE bcdtyr
