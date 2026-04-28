
! Auto-generated at 2026-04-28 18:43:09.284497 by ops-translator
SUBROUTINE bcutzl

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE bcut_kernel_zdir_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   BCUTZL
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

  !   Z-DIRECTION LEFT-HAND END

  !   *************************************************************************

  !   GLOBAL DATA
  !   ===========
  !   -------------------------------------------------------------------------
  !   -------------------------------------------------------------------------

  !   LOCAL DATA
  !   ==========
  INTEGER(KIND = 4) :: rangexyz(6)

  !   BEGIN
  !   =====

  !   =========================================================================

  !   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

  !   =========================================================================

  !   CONSTANT W-VELOCITY
  !   PARAMETER I1=1, R1=W-VELOCITY
  IF (nzlprm(1) == 1) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
    CALL bcut_kernel_zdir_host("bcut_kernel_zdir", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dudtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dvdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(rzlprm, nbcprr, "real(kind=8)", OPS_READ))

  END IF

  !   =========================================================================

END SUBROUTINE bcutzl