
! Auto-generated at 2026-04-28 18:43:09.255345 by ops-translator
SUBROUTINE bcttyl

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE bcdt_kernel_ydir_module
  USE bcdt_kernel_ydir_eqa_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   ************************************************************************

  !   BCTTYL
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

  !   CHANGE RECORD
  !   -------------
  !   26-OCT-2013:  CREATED
  !   09-MAY-2015:  RSC MODIFIED FOR ISOTHERMAL WALL

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR TEMPERATURE
  !   AND ITS TIME DERIVATIVE

  !   Y-DIRECTION LEFT-HAND END

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

  !   EVALUATE AND RETURN STRTYL,DTDTYL
  rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
  CALL bcdt_kernel_ydir_host("bcdt_kernel_ydir", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dtdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(trin, 1, "real(kind=8)", OPS_READ))

    !   =========================================================================

    !   ISOTHERMAL WALL
    IF (nsbcyl == nsbcw2) THEN
    rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
    CALL bcdt_kernel_ydir_eqa_host("bcdt_kernel_ydir", senga_grid, 3, rangexyz, &
ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dtdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(rylprm, nbcprr, "real(kind=8)", OPS_READ))

  END IF

  !   =========================================================================

END SUBROUTINE bcttyl