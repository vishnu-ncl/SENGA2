
! Auto-generated at 2026-04-28 18:43:09.310222 by ops-translator
SUBROUTINE bcytzr

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE bcyt_kernel_zdir_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   BCYTZR
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
  !   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
  !   AND THEIR TIME DERIVATIVES

  !   Z-DIRECTION RIGHT-HAND END

  !   *************************************************************************

  !   GLOBAL DATA
  !   ===========
  !   -------------------------------------------------------------------------
  !   -------------------------------------------------------------------------

  !   LOCAL DATA
  !   ==========
  INTEGER(KIND = 4) :: ispec
  INTEGER(KIND = 4) :: rangexyz(6)

  !   BEGIN
  !   =====

  !   =========================================================================

  !   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

  !   =========================================================================

  !   EVALUATE AND RETURN STRYZR,DYDTZR
  DO ispec = 1, nspec
    rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
    CALL bcyt_kernel_zdir_host("bcyt_kernel_zdir", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dydtzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(yrin, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

  END DO

  !   =========================================================================

END SUBROUTINE bcytzr