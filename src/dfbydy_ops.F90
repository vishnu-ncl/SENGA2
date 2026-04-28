
! Auto-generated at 2026-04-28 18:43:09.046275 by ops-translator
SUBROUTINE dfbydy(functn, fderiv)

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE dfbydy_kernel_null_module
  USE dfbydy_kernel_main_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   DFBYDY
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT

  !   CHANGE RECORD
  !   -------------
  !   01-AUG-1996:  CREATED
  !   11-APR-2003:  RSC MODIFIED FOR SENGA2
  !   10-OCT-2004:  RSC NULL VERSION

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   EVALUATES FIRST Y-DERIVATIVE OF SPECIFIED FUNCTION

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

  IF (nyglbl == 1) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL dfbydy_kernel_null_host("dfbydy_null", senga_grid, 3, rangexyz, &
ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  ELSE
    !       INTERIOR SCHEME
    !       ===============

    !       TENTH ORDER EXPLICIT DIFFERENCES
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL dfbydy_kernel_main_host("dfbydy_main_scheme", senga_grid, 3, rangexyz, &
ops_arg_dat(functn, 1, s3d_p050_to_m050_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(nyglbl_ops, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(nendyl, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(nendyr, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(nbound_ops, 1, "integer(kind=4)", OPS_READ), &
ops_arg_idx())

  END IF

  !   =========================================================================

END SUBROUTINE dfbydy