SUBROUTINE bcdtyl

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga
 
!   *************************************************************************

!   BCDTYL
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

!   Y-DIRECTION LEFT-HAND END

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!   =========================================================================

!   EVALUATE AND RETURN STRDYL,DDDTYL
    rangexyz = (/istal,istol,1,1,kstal,kstol/)
    call ops_par_loop(bcdt_kernel_ydir, "bcdt_kernel_ydir", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                    ops_arg_dat(d_dddtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(drin, 1, "real(dp)", OPS_READ))

!   =========================================================================

END SUBROUTINE bcdtyl
