
! Auto-generated at 2026-04-28 18:43:09.275434 by ops-translator
SUBROUTINE bcutxr

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE bcut_kernel_xdir_eqf_module
  USE bcut_kernel_xdir_eqg_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   BCUTXR
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
  !   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
  !   AND THEIR TIME DERIVATIVES

  !   X-DIRECTION RIGHT-HAND END

  !   *************************************************************************

  !   GLOBAL DATA
  !   ===========
  !   -------------------------------------------------------------------------
  !   -------------------------------------------------------------------------

  !   LOCAL DATA
  !   ==========
  !KA   FIX INFLOW BC
  !KA      real(kind=8) BTIME
  REAL(KIND = 8) :: fornow, argmnt
  INTEGER(KIND = 4) :: jc, kc
  INTEGER(KIND = 4) :: rangexyz(6)
  REAL(KIND = 8) :: init_val1, init_val2
  REAL(KIND = 8) :: rxrprm_1

  !   BEGIN
  !   =====

  !   =========================================================================

  !   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
  !KA   FIX INFLOW BC
  !KA      BTIME = ETIME + RKTIM(IRKSTP)

  !   =========================================================================

  !   CONSTANT U-VELOCITY
  !   PARAMETER I1=1, R1=U-VELOCITY
  IF (nxrprm(1) == 1) THEN
    rxrprm_1 = rxrprm(1)
    rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL bcut_kernel_xdir_eqf_host("bcut_kernel_xdir_eqF", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dudtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dvdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dwdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(rxrprm_1, 1, "real(kind=8)", OPS_READ))

  END IF

    !   =========================================================================

    !   SINUSOIDAL U-VELOCITY
    !   PARAMETER I1=2, R1=AMPLITUDE, R2=PERIOD
    IF (nxrprm(1) == 2) THEN
    fornow = two * pi / rxrprm(2)
    argmnt = fornow * btime
    init_val1 = rxrprm(1) * SIN(argmnt)
    init_val2 = fornow * rxrprm(1) * COS(argmnt)

    rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL bcut_kernel_xdir_eqg_host("bcut_kernel_xdir_eqG", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dudtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dvdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dwdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(init_val1, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(init_val2, 1, "real(kind=8)", OPS_READ))

  END IF

  !   =========================================================================

END SUBROUTINE bcutxr