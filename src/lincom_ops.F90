
! Auto-generated at 2026-04-28 18:43:09.208985 by ops-translator
SUBROUTINE lincom

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE lincom_kernel_main_module
  USE lincom_kernel_eqa_module
  USE copy_kernel_module
  USE lincom_kernel_eqb_module
  USE lincom_kernel_eqc_module
  USE lincom_kernel_eqd_module
  USE lincom_kernel_eqe_module
  USE lincom_kernel_eqf_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   LINCOM
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

  !   CHANGE RECORD
  !   -------------
  !   15-JAN-2003:  CREATED
  !   08-AUG-2012:  RSC EVALUATE ALL SPECIES

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   COMPUTES INTERMEDIATE SOLUTION VALUES IN ERK SCHEME
  !   BY DOING LINEAR COMBINATIONS OF LEFT- AND RIGHT-HAND SIDES

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

  !   ERK SUBSTEP
  !   ===========

  !   -------------------------------------------------------------------------
  !   NOTE: ALL ERK ERROR ARRAYS ARE INITIALISED TO ZERO IN SUBROUTINE ADAPTT
  !   -------------------------------------------------------------------------

  !   DENSITY
  !   -------
  rangexyz(1) = 1
  IF (nsbcxl == nsbci3) rangexyz(1) = 2
  rangexyz(2) = nxglbl
  IF (nsbcxr == nsbci3) rangexyz(2) = nxglbl - 1
  rangexyz(3) = 1
  IF (nsbcyl == nsbci3) rangexyz(3) = 2
  rangexyz(4) = nyglbl
  IF (nsbcyr == nsbci3) rangexyz(4) = nyglbl - 1
  rangexyz(5) = 1
  IF (nsbczl == nsbci3) rangexyz(5) = 2
  rangexyz(6) = nzglbl
  IF (nsbczr == nsbci3) rangexyz(6) = nzglbl - 1

  CALL lincom_kernel_main_host("lincom_main", senga_grid, 3, rangexyz, &
ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_gbl(rkerr, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rklhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rkrhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(irkstp, 1, "integer(kind=4)", OPS_READ))

  !   -------------------------------------------------------------------------
  !   U-VELOCITY
  !   ----------
  rangexyz(1) = 1
  IF (nsbcxl == nsbci2 .OR. nsbcxl == nsbci3 .OR. nsbcxl == nsbcw1 .OR. nsbcxl == nsbcw2) rangexyz(1) = 2
  rangexyz(2) = nxglbl
  IF (nsbcxr == nsbci2 .OR. nsbcxr == nsbci3 .OR. nsbcxr == nsbcw1 .OR. nsbcxr == nsbcw2) rangexyz(2) = nxglbl - 1
  rangexyz(3) = 1
  IF (nsbcyl == nsbci2 .OR. nsbcyl == nsbci3 .OR. nsbcyl == nsbcw1 .OR. nsbcyl == nsbcw2) rangexyz(3) = 2
  rangexyz(4) = nyglbl
  IF (nsbcyr == nsbci2 .OR. nsbcyr == nsbci3 .OR. nsbcyr == nsbcw1 .OR. nsbcyr == nsbcw2) rangexyz(4) = nyglbl - 1
  rangexyz(5) = 1
  IF (nsbczl == nsbci2 .OR. nsbczl == nsbci3 .OR. nsbczl == nsbcw1 .OR. nsbczl == nsbcw2) rangexyz(5) = 2
  rangexyz(6) = nzglbl
  IF (nsbczr == nsbci2 .OR. nsbczr == nsbci3 .OR. nsbczr == nsbcw1 .OR. nsbczr == nsbcw2) rangexyz(6) = nzglbl - 1

  CALL lincom_kernel_main_host("lincom_main", senga_grid, 3, rangexyz, &
ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_gbl(rkerr, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rklhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rkrhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(irkstp, 1, "integer(kind=4)", OPS_READ))

  !   -------------------------------------------------------------------------
  !   V-VELOCITY
  !   ----------
  CALL lincom_kernel_main_host("lincom_main", senga_grid, 3, rangexyz, &
ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_gbl(rkerr, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rklhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rkrhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(irkstp, 1, "integer(kind=4)", OPS_READ))

  !   -------------------------------------------------------------------------
  !   W-VELOCITY
  !   ----------
  CALL lincom_kernel_main_host("lincom_main", senga_grid, 3, rangexyz, &
ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_gbl(rkerr, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rklhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rkrhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(irkstp, 1, "integer(kind=4)", OPS_READ))

  !   -------------------------------------------------------------------------
  !   STAGNATION INTERNAL ENERGY
  !   --------------------------
  rangexyz(1) = 1
  IF (nsbcxl == nsbci2 .OR. nsbcxl == nsbcw2) rangexyz(1) = 2
  rangexyz(2) = nxglbl
  IF (nsbcxr == nsbci2 .OR. nsbcxr == nsbcw2) rangexyz(2) = nxglbl - 1
  rangexyz(3) = 1
  IF (nsbcyl == nsbci2 .OR. nsbcyl == nsbcw2) rangexyz(3) = 2
  rangexyz(4) = nyglbl
  IF (nsbcyr == nsbci2 .OR. nsbcyr == nsbcw2) rangexyz(4) = nyglbl - 1
  rangexyz(5) = 1
  IF (nsbczl == nsbci2 .OR. nsbczl == nsbcw2) rangexyz(5) = 2
  rangexyz(6) = nzglbl
  IF (nsbczr == nsbci2 .OR. nsbczr == nsbcw2) rangexyz(6) = nzglbl - 1

  CALL lincom_kernel_main_host("lincom_main", senga_grid, 3, rangexyz, &
ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_gbl(rkerr, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rklhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rkrhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(irkstp, 1, "integer(kind=4)", OPS_READ))

  !   -------------------------------------------------------------------------
  !   SPECIES MASS FRACTIONS
  !   ----------------------
  !   RSC 08-AUG-2012 EVALUATE ALL SPECIES
  !   DO ISPEC = 1,NSPM1
  rangexyz(1) = 1
  IF (nsbcxl == nsbci2 .OR. nsbcxl == nsbci3) rangexyz(1) = 2
  rangexyz(2) = nxglbl
  IF (nsbcxr == nsbci2 .OR. nsbcxr == nsbci3) rangexyz(2) = nxglbl - 1
  rangexyz(3) = 1
  IF (nsbcyl == nsbci2 .OR. nsbcyl == nsbci3) rangexyz(3) = 2
  rangexyz(4) = nyglbl
  IF (nsbcyr == nsbci2 .OR. nsbcyr == nsbci3) rangexyz(4) = nyglbl - 1
  rangexyz(5) = 1
  IF (nsbczl == nsbci2 .OR. nsbczl == nsbci3) rangexyz(5) = 2
  rangexyz(6) = nzglbl
  IF (nsbczr == nsbci2 .OR. nsbczr == nsbci3) rangexyz(6) = nzglbl - 1

    DO ispec = 1, nspec

    CALL lincom_kernel_main_host("lincom_main", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_gbl(rkerr, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rklhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rkrhs, nrkmax, "real(kind=8)", OPS_READ), &
ops_arg_gbl(irkstp, 1, "integer(kind=4)", OPS_READ))

  END DO

    !   -------------------------------------------------------------------------

    !   VM & NC: GRADIENT OF SPECIES AT WALL EQUAL TO ZERO
    IF ((nsbcxl == nsbcw2) .OR. (nsbcxl == nsbcw1)) THEN
    DO ispec = 1, nspec

      rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
      CALL lincom_kernel_eqa_host("lincom_kernel_eqA", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ))
      CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF

    IF ((nsbcxr == nsbcw2) .OR. (nsbcxr == nsbcw1)) THEN
    DO ispec = 1, nspec

      rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL lincom_kernel_eqb_host("lincom_kernel_eqB", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ))
      CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
    END DO
  END IF

    IF ((nsbcyl == nsbcw2) .OR. (nsbcyl == nsbcw1)) THEN
    DO ispec = 1, nspec

      rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
      CALL lincom_kernel_eqc_host("lincom_kernel_eqC", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ))
      CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
    END DO
  END IF

    IF ((nsbcyr == nsbcw2) .OR. (nsbcyr == nsbcw1)) THEN
    DO ispec = 1, nspec

      rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
      CALL lincom_kernel_eqd_host("lincom_kernel_eqD", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ))
      CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF

    IF ((nsbczl == nsbcw2) .OR. (nsbczl == nsbcw1)) THEN
    DO ispec = 1, nspec

      rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
      CALL lincom_kernel_eqe_host("lincom_kernel_eqE", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ))
      CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF

    IF ((nsbczr == nsbcw2) .OR. (nsbczr == nsbcw1)) THEN
    DO ispec = 1, nspec

      rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
      CALL lincom_kernel_eqf_host("lincom_kernel_eqF", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ))
      CALL copy_kernel_host("copy", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF

  !   -------------------------------------------------------------------------

END SUBROUTINE lincom