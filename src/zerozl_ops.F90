
! Auto-generated at 2026-04-28 18:43:09.067941 by ops-translator
SUBROUTINE zerozl(farray)

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE set_zero_kernel_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   ARGUMENTS
  !   =========
  TYPE(ops_dat) :: farray

  !   LOCAL DATA
  !   ==========
  INTEGER(KIND = 4) :: rangexyz(6)

  !   BEGIN
  !   =====

  !   =========================================================================

  rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(farray, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  !   =========================================================================

END SUBROUTINE zerozl