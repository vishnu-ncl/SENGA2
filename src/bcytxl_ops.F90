
! Auto-generated at 2026-04-28 18:43:09.291838 by ops-translator
SUBROUTINE bcytxl

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE bcyt_kernel_xdir_eqa_module
  USE bcyt_kernel_xdir_eqb_module
  USE set_zero_kernel_xdir_module
  USE bcyt_kernel_xdir_eqc_module
  USE bcyt_kernel_xdir_eqd_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   BCYTXL
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
  !   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
  !   AND THEIR TIME DERIVATIVES

  !   X-DIRECTION LEFT-HAND END

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

  !   EVALUATE AND RETURN STRYXL,DYDTXL
  DO ispec = 1, nspec
    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL bcyt_kernel_xdir_eqa_host("bcyt_kernel_xdir_eqA", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dydtxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(yrin, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

  END DO

    !   VM: SYNTHETIC SCALAR INFLOW
    !   VM: NXLPRM(2)=1 IMPLIES THAT THE SCALAR SYTHETIC DIGITAL FILTERING IS ON
    IF ((nxlprm(2) == 1) .AND. (nxlprm(1) == 4) .AND. (ngbcxl == 12)) THEN
    DO ispec = 1, nspec
      rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
      CALL bcyt_kernel_xdir_eqb_host("bcyt_kernel_xdir_eqB", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dydtxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yinf2(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_yinf1(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_gbl(yrin, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(tstep, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
    END DO

    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL set_zero_kernel_xdir_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_totyxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
    DO ispec = 1, nspec - 1
      rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
      CALL bcyt_kernel_xdir_eqc_host("bcyt_kernel_xdir_eqC", senga_grid, 3, rangexyz, &
ops_arg_dat(d_totyxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
    END DO

    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL bcyt_kernel_xdir_eqd_host("bcyt_kernel_xdir_eqD", senga_grid, 3, rangexyz, &
ops_arg_dat(d_stryxl(nspec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dydtxl(nspec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yinf2(nspec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_yinf1(nspec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_totyxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_gbl(tstep, 1, "real(kind=8)", OPS_READ))

  END IF

  !   =========================================================================

END SUBROUTINE bcytxl