SUBROUTINE bcutyl

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCUTYL
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

!   Y-DIRECTION LEFT-HAND END

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!     -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!   =========================================================================

!   CONSTANT V-VELOCITY
!   PARAMETER I1=1, R1=V-VELOCITY
    IF(nylprm(1) == 1) THEN
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(bcut_kernel_ydir, "bcut_kernel_ydir", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_dudtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_dwdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_gbl(rylprm, nbcprr, "real(kind=8)", OPS_READ))

    END IF

!   =========================================================================

END SUBROUTINE bcutyl
