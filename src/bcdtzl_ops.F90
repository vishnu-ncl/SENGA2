
! Auto-generated at 2026-04-28 18:43:09.243556 by ops-translator
SUBROUTINE bcdtzl

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE bcdt_kernel_zdir_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   BCDTZL
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

  !   EVALUATE AND RETURN STRDZL,DDDTZL
  rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
  CALL bcdt_kernel_zdir_host("bcdt_kernel_zdir", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dddtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(drin, 1, "real(kind=8)", OPS_READ))

  !   =========================================================================

END SUBROUTINE bcdtzl