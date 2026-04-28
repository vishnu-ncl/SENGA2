
! Auto-generated at 2026-04-28 18:43:09.637838 by ops-translator
SUBROUTINE tempin

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE copy_kernel_sdim_to_mdim_module
  USE tempin_kernel_main_module
  USE copy_kernel_mdim_to_sdim_module
  USE set_zero_kernel_int_module
  USE set_zero_kernel_module
  USE temper_kernel_eqe_module
  USE temper_kernel_eqc_module
  USE temper_kernel_eqf_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   TEMPIN
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

  !   CHANGE RECORD
  !   -------------
  !   15-MAY-2004:  CREATED

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   INITIALISES TEMPERATURE AND PRESSURE
  !   USES BISECTION METHOD FOR ROBUSTNESS

  !   *************************************************************************

  !   GLOBAL DATA
  !   ===========
  !   -------------------------------------------------------------------------
  !   -------------------------------------------------------------------------

  !   PARAMETERS
  !   ==========
  REAL(KIND = 8), PARAMETER :: toltmp = 0.00010_8
  REAL(KIND = 8), PARAMETER :: tininc = 50.0_8
  REAL(KIND = 8), PARAMETER :: tlimlo = 200.0_8
  REAL(KIND = 8), PARAMETER :: tlimhi = 3000.0_8

  !   LOCAL DATA
  !   ==========
  INTEGER(KIND = 4) :: rangexyz(6)
  REAL(KIND = 8) :: tempor, tupper, tlower
  INTEGER(KIND = 4) :: iindex, ipower, ispec
  LOGICAL :: fnconv

  !   BEGIN
  !   =====

  !   =========================================================================

  !   TEMPERATURE AND PRESSURE
  !   ------------------------

  !   TEMPERATURE AND PRESSURE ARE PARALLEL

  rangexyz(1) = 1 - nhalox
  IF (nsbcxl == nsbco1 .OR. nsbcxl == nsbci1 .OR. nsbcxl == nsbci2 .OR. nsbcxl == nsbci3 .OR. nsbcxl == nsbcw1 .OR. nsbcxl == nsbcw2) rangexyz(1) = 1

  rangexyz(2) = nxglbl + nhalox
  IF (nsbcxr == nsbco1 .OR. nsbcxr == nsbci1 .OR. nsbcxr == nsbci2 .OR. nsbcxr == nsbci3 .OR. nsbcxr == nsbcw1 .OR. nsbcxr == nsbcw2) rangexyz(2) = nxglbl

  rangexyz(3) = 1 - nhaloy
  IF (nsbcyl == nsbco1 .OR. nsbcyl == nsbci1 .OR. nsbcyl == nsbci2 .OR. nsbcyl == nsbci3 .OR. nsbcyl == nsbcw1 .OR. nsbcyl == nsbcw2) rangexyz(3) = 1

  rangexyz(4) = nyglbl + nhaloy
  IF (nsbcyr == nsbco1 .OR. nsbcyr == nsbci1 .OR. nsbcyr == nsbci2 .OR. nsbcyr == nsbci3 .OR. nsbcyr == nsbcw1 .OR. nsbcyr == nsbcw2) rangexyz(4) = nyglbl

  rangexyz(5) = 1 - nhaloz
  IF (nsbczl == nsbco1 .OR. nsbczl == nsbci1 .OR. nsbczl == nsbci2 .OR. nsbczl == nsbci3 .OR. nsbczl == nsbcw1 .OR. nsbczl == nsbcw2) rangexyz(5) = 1

  rangexyz(6) = nzglbl + nhaloz
  IF (nsbczr == nsbco1 .OR. nsbczr == nsbci1 .OR. nsbczr == nsbci2 .OR. nsbczr == nsbci3 .OR. nsbczr == nsbcw1 .OR. nsbczr == nsbcw2) rangexyz(6) = nzglbl

    DO ispec = 1, nspcmx
    CALL copy_kernel_sdim_to_mdim_host("A_multidim(ispec) = B", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs_mdim, 22, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
  END DO

  CALL tempin_kernel_main_host("tempin kernel", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_yrhs_mdim, 22, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(amascp, ncofmx * ntinmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(amasct, ncofmx * ntinmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ncpoly, ntinmx * nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ncenth, ntinmx * nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(tinthi, ntinmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(trin, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(nspec, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(iproc, 1, "integer(kind=4)", OPS_READ), &
ops_arg_idx())

    DO ispec = 1, nspcmx
    CALL copy_kernel_mdim_to_sdim_host("A = B_multidim(ispec)", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs_mdim, 22, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
  END DO

    !   CONSTRUCT THE TEMPERATURE INTERVAL INDEX
    !   EVALUATE PRESSURE
    !   EVALUATE MIXTURE SPECIFIC HEAT CP
    DO iindex = 1, nintmx
    CALL set_zero_kernel_int_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
  END DO

  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE))
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    DO ispec = 1, nspec
    !       SET THE TEMPERATURE INTERVAL INDEX
    iindex = 1 + (ispec - 1) / nspimx
    ipower = ispec - (iindex - 1) * nspimx - 1

    CALL temper_kernel_eqe_host("temper eq E", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(amascp, ncofmx * ntinmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ncpoly, ntinmx * nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ncpom1, ntinmx * nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(tinthi, ntinmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ))

    !       EVALUATE (DENSITY TIMES) MIXTURE GAS CONSTANT FOR PRESSURE
    CALL temper_kernel_eqc_host("temper eq C", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

  END DO

  !   =========================================================================

  CALL temper_kernel_eqf_host("temper eq F", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

  !   =========================================================================

  !    call ops_free_dat(d_yrhs_mdim)

END SUBROUTINE tempin