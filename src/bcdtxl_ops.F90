
! Auto-generated at 2026-04-28 18:43:09.228393 by ops-translator
SUBROUTINE bcdtxl

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE bcdt_kernel_xdir_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

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
  INTEGER(KIND = 4) :: rangexyz(6)

  !   BEGIN
  !   =====

  !   =========================================================================

  !   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

  !   =========================================================================

  !   EVALUATE AND RETURN STRDXL,DDDTXL
  rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
  CALL bcdt_kernel_xdir_host("bcdt_kernel_xdir", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dddtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(drin, 1, "real(kind=8)", OPS_READ))

  !   =========================================================================

END SUBROUTINE bcdtxl