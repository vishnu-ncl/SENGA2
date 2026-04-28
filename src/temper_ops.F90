
! Auto-generated at 2026-04-28 18:43:09.624700 by ops-translator
SUBROUTINE temper

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE temper_kernel_eqa_module
  USE set_zero_kernel_md5_module
  USE set_zero_kernel_module
  USE temper_kernel_eqb_module
  USE temper_kernel_eqc_module
  USE temper_kernel_eqd_module
  USE set_zero_kernel_int_module
  USE temper_kernel_eqe_module
  USE temper_kernel_eqf_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   TEMPER
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

  !   CHANGE RECORD
  !   -------------
  !   16-NOV-2002:  CREATED

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   COMPUTES TEMPERATURE AND PRESSURE

  !   *************************************************************************

  !   GLOBAL DATA
  !   ===========
  !   -------------------------------------------------------------------------
  !   -------------------------------------------------------------------------

  INTEGER(KIND = 4) :: icp, ispec, rangexyz(6)
  INTEGER(KIND = 4) :: iindex, ipower, icoef1, icoef2

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

  !    IF (ops_is_root() == 1) THEN
  !        print *, "x: ",rangexyz(1),"-",rangexyz(2), "y: ",rangexyz(3),"-",rangexyz(4), "z: ",rangexyz(5),"-",rangexyz(6)
  !    END IF

  !   INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
  !   AND ITS DERIVATIVE
  CALL temper_kernel_eqa_host("temper eq A", senga_grid, 3, rangexyz, &
ops_arg_dat(d_tcoeff, 6, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

  !    DO icp = 1, nctmax
  !        ! OPS_RW is used instead of OPS_WRITE to get correct tiling result
  !        call ops_par_loop(set_zero_kernel_MD6, "set_zero tcoeff", senga_grid, 3, rangexyz,  &
  !                        ops_arg_dat(d_tcoeff, 6, s3d_000, "real(kind=8)", OPS_RW), &
  !                        ops_arg_gbl(icp+1, 1, "integer(kind=4)", OPS_READ))
  !    END DO

  CALL set_zero_kernel_md5_host("set_zero tderiv", senga_grid, 3, rangexyz, &
ops_arg_dat(d_tderiv, 5, s3d_000, "real(kind=8)", OPS_WRITE))

  !   USE STORE7 TO ACCUMULATE MIXTURE SPECIFIC GAS CONSTANT
  !   INITIALISE STORE7
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    !   ===================================================================

    !   RUN THROUGH ALL SPECIES
    DO ispec = 1, nspec
    !       LOCATE TEMPERATURE IN AN INTERVAL
    iindex = 1 + (ispec - 1) / nspimx
    ipower = ispec - (iindex - 1) * nspimx - 1
    icoef2 = ntbase ** ipower
    icoef1 = icoef2 * ntbase

    !       =================================================================

    !       CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
    CALL temper_kernel_eqb_host("temper eq B", senga_grid, 3, rangexyz, &
ops_arg_dat(d_tcoeff, 6, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_tderiv, 5, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(amascp, ncofmx * ntinmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(amasct, ncofmx * ntinmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ncpoly, ntinmx * nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ncenth, ntinmx * nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    !       =================================================================

    !       USE STORE7
    !       TO ACCUMULATE (DENSITY TIMES) MIXTURE SPECIFIC GAS CONSTANT
    CALL temper_kernel_eqc_host("temper eq C", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    !       =================================================================

  END DO

  !   END OF RUN THROUGH ALL SPECIES

  !   ===================================================================

  !   SOLVE FOR TEMPERATURE
  CALL temper_kernel_eqd_host("temper eq D", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_tcoeff, 6, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tderiv, 5, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_idx())

    !   FOR ALL SPECIES RELOCATE TEMPERATURE IN AN INTERVAL
    !   EVALUATE MIXTURE SPECIFIC HEAT CP
    DO iindex = 1, nintmx
    CALL set_zero_kernel_int_host("set_zero itndex", senga_grid, 3, rangexyz, &
ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
  END DO

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

  END DO

  !   =========================================================================

  CALL temper_kernel_eqf_host("temper eq F", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

  !   =========================================================================

END SUBROUTINE temper