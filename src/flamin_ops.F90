
! Auto-generated at 2026-04-28 18:43:09.528712 by ops-translator
SUBROUTINE flamin

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE copy_kernel_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   FLAMIN
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

  !   CHANGE RECORD
  !   -------------
  !   28-DEC-2003:  CREATED
  !   08-JAN-2005:  RSC INITIAL 1D LAMINAR FLAME PROFILE

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   SETS INITIAL THERMOCHEMICAL FIELD
  !   1D LAMINAR FLAME PROFILE (LEFT OR RIGHT FACING)
  !   SPECIAL FOR 21 STEP HYDROGEN MECHAMISM

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

  !   POPULATE DATA FROM THE FLAMIN INDATA GENERATED FROM STANDALONE ROUTINE
  !   ======================================================================
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_drun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

  CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_urun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

  CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_vrun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

  CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_wrun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

  CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_trun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

    DO ispec = 1, nspcmx
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrun_dump(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
  END DO

  !   =========================================================================

END SUBROUTINE flamin