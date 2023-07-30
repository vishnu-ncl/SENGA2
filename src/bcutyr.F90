SUBROUTINE bcutyr

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCUTYR
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

!   CONSTANT V-VELOCITY
!   PARAMETER I1=1, R1=V-VELOCITY
    IF(nyrprm(1) == 1) THEN
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(bcut_kernel_ydir, "bcut_kernel_ydir", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_dudtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_dvdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_dwdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_gbl(ryrprm, nbcprr, "real(kind=8)", OPS_READ))

    END IF

!   =========================================================================

END SUBROUTINE bcutyr
