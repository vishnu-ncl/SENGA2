SUBROUTINE bcttyr

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCTTYR
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   26-OCT-2013:  CREATED
!   09-MAY-2015:  RSC MODIFIED FOR ISOTHERMAL WALL

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR TEMPERATURE
!   AND ITS TIME DERIVATIVE

!   Y-DIRECTION RIGHT-HAND END

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

!   EVALUATE AND RETURN STRTYR,DTDTYR
    rangexyz = (/istal,istol,1,1,kstal,kstol/)
    call ops_par_loop(bcdt_kernel_ydir, "bcdt_kernel_ydir", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                    ops_arg_dat(d_dtdtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(trin, 1, "real(dp)", OPS_READ))

!   =========================================================================

!   ISOTHERMAL WALL
    IF(nsbcyr == nsbcw2) THEN
        rangexyz = (/istal,istol,1,1,kstal,kstol/)
        call ops_par_loop(bcdt_kernel_ydir, "bcdt_kernel_ydir", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                    ops_arg_dat(d_dtdtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(ryrprm(1), 1, "real(dp)", OPS_READ))

    END IF

!   =========================================================================

END SUBROUTINE bcttyr
