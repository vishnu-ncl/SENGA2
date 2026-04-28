
! Auto-generated at 2026-04-28 18:43:09.042740 by ops-translator
SUBROUTINE d2fdx2(functn, fderiv)

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE d2fdx2_kernel_main_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   D2FDX2
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT

  !   CHANGE RECORD
  !   -------------
  !   01-AUG-1996:  CREATED
  !   06-APR-2003:  RSC MODIFIED FOR SENGA2

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   EVALUATES SECOND X-DERIVATIVE OF SPECIFIED FUNCTION
  !   EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
  !   EXPLICIT 8TH,6TH,4TH,4TH ORDER END CONDITIONS

  !   *************************************************************************

  !   ARGUMENTS
  !   =========
  TYPE(ops_dat) :: functn, fderiv

  !   LOCAL DATA
  !   ==========
  INTEGER(KIND = 4) :: rangexyz(6)

  !   BEGIN
  !   =====

  !   =========================================================================

  !   INTERIOR SCHEME
  !   ===============

  !   TENTH ORDER EXPLICIT DIFFERENCES

  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL d2fdx2_kernel_main_host("d2fdx2_main_scheme", senga_grid, 3, rangexyz, &
ops_arg_dat(functn, 1, s3d_p500_to_m500_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(nxglbl_ops, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(nendxl, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(nendxr, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(nbound_ops, 1, "integer(kind=4)", OPS_READ), &
ops_arg_idx())

  !   =========================================================================

END SUBROUTINE d2fdx2