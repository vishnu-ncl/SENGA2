SUBROUTINE bcutzr

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCUTZR
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
!   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
!   AND THEIR TIME DERIVATIVES

!   Z-DIRECTION RIGHT-HAND END

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

!   CONSTANT W-VELOCITY
!   PARAMETER I1=1, R1=W-VELOCITY
    IF(nzrprm(1) == 1) THEN
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(bcut_kernel_zdir, "bcut_kernel_zdir", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_dudtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_dvdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_dwdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_gbl(rzrprm, nbcprr, "real(kind=8)", OPS_READ))

    END IF

!   =========================================================================

END SUBROUTINE bcutzr
