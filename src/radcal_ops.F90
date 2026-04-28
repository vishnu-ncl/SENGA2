
! Auto-generated at 2026-04-28 18:43:09.513869 by ops-translator
SUBROUTINE radcal

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE set_zero_kernel_module
  USE radcal_kernel_meancoef_module
  USE radcal_kernel_addspecies_module
  USE radcal_kernel_addradiation_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   RADCAL
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

  !   CHANGE RECORD
  !   -------------
  !   14-JUL-2013:  CREATED

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   RADIATION TREATMENT
  !   USING OPTICALLY THIN ASSUMPTION: Ju et al: JFM 342, 315-334, 1997.
  !   AFTER TOM DUNSTAN 2012

  !   *************************************************************************

  !   GLOBAL DATA
  !   ===========
  !   -------------------------------------------------------------------------
  !   -------------------------------------------------------------------------

  !   LOCAL DATA
  !   ==========
  REAL(KIND = 8) :: plspec, fornow
  INTEGER(KIND = 4) :: ispec, jspec, icp
  INTEGER(KIND = 4) :: rangexyz(6)

  !   BEGIN
  !   =====

  !   =========================================================================

  !   BUILD THE PLANCK MEAN ABSORPTION COEFFICIENT OF THE MIXTURE
  !   -----------------------------------------------------------

  !   INITIALISE THE ACCUMULATOR
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    !   -------------------------------------------------------------------------

    !   RUN THROUGH ALL RADIATING SPECIES
    DO jspec = 1, nsprad

    !       PLANCK MEAN ABSORPTION COEFFICIENT OF EACH SPECIES
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL radcal_kernel_meancoef_host("PLANCK MEAN ABSORPTION COEF", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(akprad, ncfrmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(nkprad, nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(nkprm1, nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(jspec, 1, "integer(kind=4)", OPS_READ))

    !       SPECIES ID
    ispec = nsprid(jspec)

    !       ADD THE SPECIES CONTRIBUTION
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL radcal_kernel_addspecies_host("ADD THE SPECIES CONTRIBUTION", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

  END DO

  !   =========================================================================

  !   INCLUDE THE RADIATION TERM IN THE ENERGY EQUATION
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL radcal_kernel_addradiation_host("INCLUDE THE RADIATION TERM IN THE ENERGY EQUATION", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

  !   =========================================================================

END SUBROUTINE radcal