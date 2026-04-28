
! Auto-generated at 2026-04-28 18:43:09.240536 by ops-translator
SUBROUTINE bcdtyr

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE bcdt_kernel_ydir_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

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
  INTEGER(KIND = 4) :: rangexyz(6)

  !   BEGIN
  !   =====

  !   =========================================================================

  !   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

  !   =========================================================================

  !   EVALUATE AND RETURN STRDYR,DDDTYR
  rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
  CALL bcdt_kernel_ydir_host("bcdt_kernel_ydir", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dddtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(drin, 1, "real(kind=8)", OPS_READ))

  !   =========================================================================

END SUBROUTINE bcdtyr