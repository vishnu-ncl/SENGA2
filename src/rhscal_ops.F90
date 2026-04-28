
! Auto-generated at 2026-04-28 18:43:09.164210 by ops-translator
SUBROUTINE rhscal

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE boundary_kernel_cpandgas_xdir_module
  USE boundary_kernel_cpandgas_ydir_module
  USE boundary_kernel_cpandgas_zdir_module
  USE maths_kernel_eqx_module
  USE maths_kernel_eqbr_xdir_module
  USE maths_kernel_eqbr_ydir_module
  USE maths_kernel_eqbr_zdir_module
  USE maths_kernel_eqt_module
  USE boundary_kernel_internalenergy_xdir_module
  USE boundary_kernel_internalenergy_ydir_module
  USE boundary_kernel_internalenergy_zdir_module
  USE maths_kernel_eqw_module
  USE maths_kernel_eqz_module
  USE maths_kernel_eqas_module
  USE boundary_kernel_temperature_xdir_module
  USE boundary_kernel_temperature_ydir_module
  USE boundary_kernel_temperature_zdir_module
  USE maths_kernel_eqav_module
  USE maths_kernel_eqbc_module
  USE set_zero_kernel_module
  USE maths_kernel_eqbd_module
  USE maths_kernel_eqbe_module
  USE maths_kernel_eqah_module
  USE hf_kernel_eqa_module
  USE hf_kernel_eqb_module
  USE hf_kernel_eqc_module
  USE hf_kernel_eqd_module
  USE hf_kernel_eqe_module
  USE hf_kernel_eqf_module
  USE maths_kernel_eqaf_module
  USE copy_kernel_module
  USE copy_kernel_xdir_module
  USE copy_kernel_ydir_module
  USE copy_kernel_zdir_module
  USE maths_kernel_eqa_module
  USE maths_kernel_eqap_module
  USE maths_kernel_eqaq_module
  USE boundary_kernel_mass_xdir_module
  USE boundary_kernel_mass_ydir_module
  USE boundary_kernel_mass_zdir_module
  USE maths_kernel_eqat_module
  USE maths_kernel_eqh_module
  USE copy_kernel_sdim_to_mdim_module
  USE maths_kernel_eqbfg_module
  USE maths_kernel_eqbh_module
  USE maths_kernel_eqbs_module
  USE maths_kernel_eqaa_fused_module
  USE maths_kernel_eqbl_module
  USE boundary_kernel_speciesh_xdir_module
  USE boundary_kernel_speciesh_ydir_module
  USE boundary_kernel_speciesh_zdir_module
  USE maths_kernel_eql_module
  USE maths_kernel_eqba_module
  USE maths_kernel_eqai_module
  USE hf_kernel_eqs_module
  USE hf_kernel_eqt_module
  USE hf_kernel_equ_module
  USE hf_kernel_eqv_module
  USE hf_kernel_eqw_module
  USE hf_kernel_eqx_module
  USE maths_kernel_eqbb_module
  USE maths_kernel_eqag_module
  USE maths_kernel_eqo_module
  USE maths_kernel_eqy_module
  USE maths_kernel_eqm_module
  USE hf_kernel_eqm_module
  USE hf_kernel_eqn_module
  USE hf_kernel_eqo_module
  USE hf_kernel_eqp_module
  USE hf_kernel_eqq_module
  USE hf_kernel_eqr_module
  USE hf_kernel_eqg_module
  USE hf_kernel_eqh_module
  USE hf_kernel_eqi_module
  USE hf_kernel_eqj_module
  USE hf_kernel_eqk_module
  USE hf_kernel_eql_module
  USE maths_kernel_eqam_module
  USE maths_kernel_eqbijk_module
  USE maths_kernel_eqan_module
  USE boundary_kernel_density_xdir_module
  USE boundary_kernel_density_ydir_module
  USE boundary_kernel_density_zdir_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   RHSCAL
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

  !   CHANGE RECORD
  !   -------------
  !   12-NOV-2002:  CREATED
  !   26-OCT-2008:  RSC/TDD BUG FIX FZLCON
  !   08-AUG-2012:  RSC EVALUATE ALL SPECIES
  !   17-APR-2013:  RSC MIXTURE AVERAGED TRANSPORT
  !   14-JUL-2013:  RSC RADIATION HEAT LOSS
  !   08-JUN-2015:  RSC REMOVE Nth SPECIES TREATMENT
  !   08-JUN-2015:  RSC UPDATED WALL BCS
  !   01-DEC-2022:  VM CONSERVATIVE-PRIMITIVE FIX/SORET FIX

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   COMPUTES RIGHT-HAND-SIDES FOR TIME INTEGRATION OF SCALAR PDEs
  !   INCLUDES MULTIPLE SCALARS AND MULTI-STEP CHEMISTRY
  !   ENERGY EQUATION REQUIRES PRESSURE-WORK AND VISCOUS WORK TERMS
  !   COMPUTED IN RHSVEL

  !   *************************************************************************

  !   GLOBAL DATA
  !   ===========
  !   -------------------------------------------------------------------------

  !   -------------------------------------------------------------------------

  !   LOCAL DATA
  !   ==========
  INTEGER(KIND = 4) :: ispec, jspec
  INTEGER(KIND = 4) :: iindex, ipower, icoef1, icoef2
  LOGICAL :: flmtds
  INTEGER(KIND = 4) :: rangexyz(6)
  REAL(KIND = 8) :: rgspec_ispec

  !    BEGIN
  !    =====

  !   =========================================================================
  !   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !   =========================================================================

  !   EVALUATE THE TEMPERATURE
  !   ------------------------
  !   ALSO PRESSURE, MIXTURE CP AND MIXTURE GAS CONSTANT
  CALL temper
  !    call temper_fused      - need to update when nintmx=3, currently working when it is 2

    !                                                             PRUN,TRUN = P,T
    !                                                         STORE7 = RHO*MIX RG
    !   =========================================================================

    !   COLLECT MIXTURE CP AND GAS CONSTANT FOR BCs
    !   -------------------------------------------

    !   X-DIRECTION
    IF (fxlcnv) THEN
    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL boundary_kernel_cpandgas_xdir_host("COLLECT CP AND GAS FOR BCs - RHSCAL 75", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

  END IF
  IF (fxrcnv) THEN
    rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL boundary_kernel_cpandgas_xdir_host("COLLECT CP AND GAS FOR BCs - RHSCAL 85", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

  END IF

    !   Y-DIRECTION
    IF (fylcnv) THEN
    rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
    CALL boundary_kernel_cpandgas_ydir_host("COLLECT CP AND GAS FOR BCs - RHSCAL 97", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

  END IF
  IF (fyrcnv) THEN
    rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
    CALL boundary_kernel_cpandgas_ydir_host("COLLECT CP AND GAS FOR BCs - RHSCAL 107", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

  END IF

    !   Z-DIRECTION
    IF (fzlcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
    CALL boundary_kernel_cpandgas_zdir_host("COLLECT CP AND GAS FOR BCs - RHSCAL 119", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

  END IF
  IF (fzrcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
    CALL boundary_kernel_cpandgas_zdir_host("COLLECT CP AND GAS FOR BCs - RHSCAL 129", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

  END IF

  !                                                            ALL STORES CLEAR
  !   =========================================================================

  !   MASS FLUX DIVERGENCE
  !   --------------------
  !   URHS,VRHS,WRHS CONTAIN RHO U, RHO V, RHO W

  CALL dfbydx(d_urhs, d_store1)
  CALL dfbydy(d_vrhs, d_store2)
  CALL dfbydz(d_wrhs, d_store3)

  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL maths_kernel_eqx_host("A = B+C+D - RHSCAL 150", senga_grid, 3, rangexyz, &
ops_arg_dat(d_divm, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    !   TRANSVERSE TERMS FOR BOUNDARIES (IMPLEMENTED BY NC)
    !   X-DIRECTION
    IF (fxlcnv) THEN
    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqbr_xdir_host("RHSCAL 161", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t1bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))
  END IF

    IF (fxrcnv) THEN
    rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqbr_xdir_host("RHSCAL 169", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t1bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))
  END IF

    !   Y-DIRECTION
    IF (fylcnv) THEN
    rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
    CALL maths_kernel_eqbr_ydir_host("RHSCAL 178", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t1byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))
  END IF

    IF (fyrcnv) THEN
    rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqbr_ydir_host("RHSCAL 186", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t1byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))
  END IF

    !   Z-DIRECTION
    IF (fzlcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
    CALL maths_kernel_eqbr_zdir_host("RHSCAL 195", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t1bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
  END IF

    IF (fzrcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
    CALL maths_kernel_eqbr_zdir_host("RHSCAL 195", senga_grid, 3, rangexyz, &
ops_arg_dat(d_t1bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))
  END IF

  !=========================REFER TO EQ. 3.74 OF LODATO'S THESIS=================
  !                                                            ALL STORES CLEAR
  !   =========================================================================
  !   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !   =========================================================================

  !   INTERNAL ENERGY EQUATION
  !   ========================

  !   CONVERT INTERNAL ENERGY
  !   -----------------------

  !   ERHS CONTAINS RHO E: CONVERT TO E
  !   E IS PARALLEL
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
  CALL maths_kernel_eqt_host("A = A/B - RHSCAL 192", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

    !                                                            ALL STORES CLEAR
    !   =========================================================================

    !   COLLECT INTERNAL ENERGY FOR BCs
    !   -------------------------------

    !   X-DIRECTION
    IF (fxlcnv) THEN
    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL boundary_kernel_internalenergy_xdir_host("COLLECT INTERNAL ENERGY FOR BCs - RHSCAL 205", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

  END IF
  IF (fxrcnv) THEN
    rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL boundary_kernel_internalenergy_xdir_host("COLLECT INTERNAL ENERGY FOR BCs - RHSCAL 212", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

  END IF

    !   Y-DIRECTION
    IF (fylcnv) THEN
    rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
    CALL boundary_kernel_internalenergy_ydir_host("COLLECT INTERNAL ENERGY FOR BCs - RHSCAL 221", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

  END IF
  IF (fyrcnv) THEN
    rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
    CALL boundary_kernel_internalenergy_ydir_host("COLLECT INTERNAL ENERGY FOR BCs - RHSCAL 228", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

  END IF

    !   Z-DIRECTION
    IF (fzlcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
    CALL boundary_kernel_internalenergy_zdir_host("COLLECT INTERNAL ENERGY FOR BCs - RHSCAL 237", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

  END IF
  IF (fzrcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
    CALL boundary_kernel_internalenergy_zdir_host("COLLECT INTERNAL ENERGY FOR BCs - RHSCAL 244", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

  END IF

  !                                                            ALL STORES CLEAR
  !   =========================================================================

  !   E EQUATION: CONVECTIVE TERMS
  !   ----------------------------
  !   HALF E DIV RHO U

  !   COLLECT E DIV RHO U IN STORE4 FOR NOW
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL maths_kernel_eqw_host("A = B*C - RHSCAL 259", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_divm, 1, s3d_000, "real(kind=8)", OPS_READ))

  !                                                        STORE4 = E DIV RHO U
  !   =========================================================================

  !   E EQUATION: CONVECTIVE TERMS
  !   ----------------------------
  !   HALF DIV RHO U E

  !   D/DX RHO U E
  !   RHO U E IS PARALLEL
  rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
  CALL maths_kernel_eqw_host("A = B*C - RHSCAL 274", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ))

  CALL dfbydx(d_store7, d_store1)

  !   D/DY RHO V E
  !   RHO V E IS PARALLEL
  rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
  CALL maths_kernel_eqw_host("A = B*C - RHSCAL 284", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ))

  CALL dfbydy(d_store7, d_store2)

  !   D/DZ RHO W E
  !   RHO W E IS PARALLEL
  rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
  CALL maths_kernel_eqw_host("A = B*C - RHSCAL 294", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ))

  CALL dfbydz(d_store7, d_store3)

  !   COLLECT DIV RHO U E IN STORE4 FOR NOW
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL maths_kernel_eqz_host("A = A+B+C+D - RHSCAL 303", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

  !                                          STORE4 = E DIV RHO U + DIV RHO U E
  !   =========================================================================

  !   E EQUATION: CONVECTIVE TERMS
  !   ----------------------------
  !   HALF RHO U.DEL E

  CALL dfbydx(d_erhs, d_store1)
  CALL dfbydy(d_erhs, d_store2)
  CALL dfbydz(d_erhs, d_store3)

  !   COLLECT ALL CONVECTIVE TERMS IN ERHS
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL maths_kernel_eqas_host("A = -half*(B+C*D+E*F+G*H) - RHSCAL 322", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ))

  !   -------------------------------------
  !   E EQUATION: CONVECTIVE TERMS COMPLETE
  !   -------------------------------------
  !                                                            ALL STORES CLEAR
  !   =========================================================================

  !   E-EQUATION: HEAT FLUX TERMS
  !   ---------------------------

  !   TEMPERATURE GRADIENTS
  CALL dfbydx(d_trun, d_store1)
  CALL dfbydy(d_trun, d_store2)
  CALL dfbydz(d_trun, d_store3)

    !                                                       STORE1,2,3 = DTDX,Y,Z
    !   =========================================================================

    !   COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs
    !   ---------------------------------------------

    !   X-DIRECTION
    IF (fxlcnv) THEN
    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL boundary_kernel_temperature_xdir_host("COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs - RHSCAL 355", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcltxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

  END IF
  IF (fxrcnv) THEN
    rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL boundary_kernel_temperature_xdir_host("COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs - RHSCAL 364", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcltxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

  END IF

    !   Y-DIRECTION
    IF (fylcnv) THEN
    rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
    CALL boundary_kernel_temperature_ydir_host("COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs - RHSCAL 375", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcltyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

  END IF
  IF (fyrcnv) THEN
    rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
    CALL boundary_kernel_temperature_ydir_host("COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs - RHSCAL 384", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcltyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

  END IF

    !   Z-DIRECTION
    IF (fzlcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
    CALL boundary_kernel_temperature_zdir_host("COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs - RHSCAL 395", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcltzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

  END IF
  IF (fzrcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
    CALL boundary_kernel_temperature_zdir_host("COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs - RHSCAL 404", senga_grid, 3, rangexyz, &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcltzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

  END IF

  !                                                       STORE1,2,3 = DTDX,Y,Z
  !   =========================================================================

  !   E-EQUATION: HEAT FLUX TERMS
  !   ---------------------------

  !   THERMAL CONDUCTIVITY
  !   ANALYTICAL FUNCTION OF TEMPERATURE
  !   TRANSP CONTAINS MIXTURE CP
  !   STORE CONDUCTIVITY/CP IN TRANSP FOR USE IN DIFFUSIVITY AND VISCOSITY
  !   STORE CONDUCTIVITY IN STORE7 FOR NOW

  !   THERMAL CONDUCTIVITY IS PARALLEL
  rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
  CALL maths_kernel_eqav_host("THERMAL CONDUCTIVITY - RHSCAL 426", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ))

    !   MIXTURE AVERAGED TRANSPORT
    !   RSC 17-APR-2013
    !   THERMAL CONDUCTIVITY

    IF (flmavt) THEN

    !#if defined(OPS_LAZY) || defined(OPS_WITH_CUDAFOR) || defined(OPS_WITH_OMPOFFLOADFOR)

    rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
    CALL maths_kernel_eqbc_host("THERMAL CONDUCTIVITY - part 1 - RHSCAL 437", senga_grid, 3, rangexyz, &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(tdifgb, 1, "real(kind=8)", OPS_READ))

    CALL set_zero_kernel_host("set_zero - RHSCAL 442", senga_grid, 3, rangexyz, &
ops_arg_dat(d_combo1, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_host("set_zero - RHSCAL 444", senga_grid, 3, rangexyz, &
ops_arg_dat(d_combo2, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    CALL set_zero_kernel_host("set_zero - RHSCAL 446", senga_grid, 3, rangexyz, &
ops_arg_dat(d_combo3, 1, s3d_000, "real(kind=8)", OPS_WRITE))

      !       CONDUCTIVITY FOR EACH SPECIES
      DO ispec = 1, nspec
      CALL maths_kernel_eqbd_host("THERMAL CONDUCTIVITY - part 2 - RHSCAL 451", senga_grid, 3, rangexyz, &
ops_arg_dat(d_combo1, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_combo2, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_combo3, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(condco, nccfmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ncocon, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ncocm1, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
    END DO

    CALL maths_kernel_eqbe_host("THERMAL CONDUCTIVITY - part 3 - RHSCAL 464", senga_grid, 3, rangexyz, &
ops_arg_dat(d_combo1, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_combo2, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_combo3, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

    !#else
    !   Ashutosh - 22 July 24, observed precision difference when using fused equation, need to debug

    !        rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    !        DO jspec = 1, nspec
    !            call ops_par_loop(copy_kernel_sdim_to_mdim, "A_multidim(ispec) = B", senga_grid, 3, rangexyz,  &
    !                            ops_arg_dat(d_yrhs_mdim, 31, s3d_000, "real(kind=8)", OPS_RW), &
    !                            ops_arg_dat(d_yrhs(jspec), 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                            ops_arg_gbl(jspec, 1, "integer(kind=4)", OPS_READ))
    !        END DO

    !        call ops_par_loop(maths_kernel_eqBCDE, "THERMAL CONDUCTIVITY - RHSCAL 527", senga_grid, 3, rangexyz, &
    !                        ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
    !                        ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
    !                        ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
    !                        ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                        ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                        ops_arg_dat(d_yrhs_mdim, 31, s3d_000, "real(kind=8)", OPS_READ), &
    !                        ops_arg_gbl(condco, nccfmx*nspcmx, "real(kind=8)", OPS_READ), &
    !                        ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
    !                        ops_arg_gbl(tdifgb, 1, "real(kind=8)", OPS_READ), &
    !                        ops_arg_gbl(ncocon, 1, "integer(kind=4)", OPS_READ), &
    !                        ops_arg_gbl(ncocm1, 1, "integer(kind=4)", OPS_READ))

    !#endif

  END IF

  !   CONDUCTIVITY GRADIENTS
  CALL dfbydx(d_store7, d_store4)
  CALL dfbydy(d_store7, d_store5)
  CALL dfbydz(d_store7, d_store6)

  !   BOUNDARY CONDITIONS
  !   BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fxlcon) CALL zeroxl(d_store4)
  IF (fxrcon) CALL zeroxr(d_store4)

  !   BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fylcon) CALL zeroyl(d_store5)
  IF (fyrcon) CALL zeroyr(d_store5)

  !   BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  !   RSC/TDD BUG FIX FZLCON
  IF (fzlcon) CALL zerozl(d_store6)
  IF (fzrcon) CALL zerozr(d_store6)

  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL maths_kernel_eqah_host("A = A+(B*C+D*E+F*G) - RHSCAL 494", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    !                                                       STORE1,2,3 = DTDX,Y,Z
    !                                                       STORE7 = CONDUCTIVITY
    !   =========================================================================

    !   E-EQUATION: HEAT FLUX TERMS
    !   ---------------------------
    !   WALL BC: THERMAL CONDUCTION TERMS
    IF (fxlcnw) THEN
    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL hf_kernel_eqa_host("hf_kernel_eqA - RHSCAL 512", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

  END IF
  IF (fxrcnw) THEN
    rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL hf_kernel_eqb_host("hf_kernel_eqB - RHSCAL 521", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

  END IF
  IF (fylcnw) THEN
    rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
    CALL hf_kernel_eqc_host("hf_kernel_eqC - RHSCAL 530", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store2, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

  END IF
  IF (fyrcnw) THEN
    rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
    CALL hf_kernel_eqd_host("hf_kernel_eqD - RHSCAL 539", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store2, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

  END IF
  IF (fzlcnw) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
    CALL hf_kernel_eqe_host("hf_kernel_eqE - RHSCAL 548", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store3, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

  END IF
  IF (fzrcnw) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
    CALL hf_kernel_eqf_host("hf_kernel_eqF - RHSCAL 557", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store3, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))
  END IF

  !   =========================================================================

  !   E-EQUATION: HEAT FLUX TERMS
  !   ---------------------------
  !   SECOND DERIVATIVE TERMS

  !   TEMPERATURE SECOND DERIVATIVES
  CALL d2fdx2(d_trun, d_store1)
  CALL d2fdy2(d_trun, d_store2)
  CALL d2fdz2(d_trun, d_store3)

  !   BOUNDARY CONDITIONS
  !   BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fxlcon) CALL zeroxl(d_store1)
  IF (fxrcon) CALL zeroxr(d_store1)

  !   BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fylcon) CALL zeroyl(d_store2)
  IF (fyrcon) CALL zeroyr(d_store2)

  !   BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  !   RSC 28-JUN-2015 BUG FIX FZLCON
  IF (fzlcon) CALL zerozl(d_store3)
  IF (fzrcon) CALL zerozr(d_store3)

  !   COLLECT CONDUCTIVITY TERMS
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL maths_kernel_eqaf_host("A = A+(B+C+D)*E - RHSCAL 591", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

  !   ---------------------------------------------------
  !   E-EQUATION: FURTHER HEAT FLUX TERMS EVALUATED BELOW
  !   ---------------------------------------------------
  !   E-EQUATION: PRESSURE-WORK AND VISCOUS WORK TERMS
  !               EVALUATED IN SUBROUTINE RHSVEL
  !   ---------------------------------------------------
  !                                                            ALL STORES CLEAR
  !   =========================================================================

  !   E-EQUATION: RADIATION HEAT LOSS
  !   -------------------------------
  IF (flradn) CALL radcal

  !   =========================================================================

  !   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !   =========================================================================

  !   SPECIES MASS FRACTION EQUATIONS
  !   ===============================

  !   REACTION RATE FOR ALL SPECIES
  !   -----------------------------
  CALL chrate
  !---UA
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  DO ispec = 1, nspec
    CALL copy_kernel_host("copy - RHSCAL 625", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rrte(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

  END DO
  !---end-UA
    !                                                        RATE = REACTION RATE
    !   =========================================================================

    !   COLLECT REACTION RATE FOR BCs
    !   -----------------------------

    !   X-DIRECTION
    IF (fxlcnv) THEN
    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    DO ispec = 1, nspec
      CALL copy_kernel_xdir_host("copy - RHSCAL 641", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF
  IF (fxrcnv) THEN
    rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
    DO ispec = 1, nspec
      CALL copy_kernel_xdir_host("copy - RHSCAL 650", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF

    !   Y-DIRECTION
    IF (fylcnv) THEN
    rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
    DO ispec = 1, nspec
      CALL copy_kernel_ydir_host("copy - RHSCAL 661", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF
  IF (fyrcnv) THEN
    rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
    DO ispec = 1, nspec
      CALL copy_kernel_ydir_host("copy - RHSCAL 670", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF

    !   Z-DIRECTION
    IF (fzlcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
    DO ispec = 1, nspec
      CALL copy_kernel_zdir_host("copy - RHSCAL 681", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF
  IF (fzrcnv) THEN
    rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
    DO ispec = 1, nspec
      CALL copy_kernel_zdir_host("copy - RHSCAL 690", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

    END DO
  END IF

  !                                                        RATE = REACTION RATE
  !   =========================================================================

  !   ZERO THE ACCUMULATORS FOR THE DIFFUSION CORRECTION VELOCITY
  !   AND ITS DIVERGENCE

  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL set_zero_kernel_host("set_zero - RHSCAL 704", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero - RHSCAL 707", senga_grid, 3, rangexyz, &
ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero - RHSCAL 710", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_WRITE))

  CALL set_zero_kernel_host("set_zero - RHSCAL 713", senga_grid, 3, rangexyz, &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))


  !   ZERO THE ACCUMULATOR FOR THE MIXTURE ENTHALPY
  !   MIXTURE H IS PARALLEL

  rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
  CALL set_zero_kernel_host("set_zero - RHSCAL 721", senga_grid, 3, rangexyz, &
ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    !                                                        RATE = REACTION RATE
    !                                                         VTMP = DIV CORR VEL
    !                                                     WTMP = MIXTURE H
    !   =========================================================================

    !   MIXTURE AVERAGED TRANSPORT
    !   RSC 17-APR-2013
    !   EVALUATE FIRST AND SECOND DERIVATIVES
    !   OF LN(MIXTURE MOLAR MASS), LN(PRESSURE) AND LN(TEMPERATURE)

    !   MIXTURE MOLAR MASS
    IF (flmixw) THEN

    rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
    CALL maths_kernel_eqa_host("A = log(B) - RHSCAL 738", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_READ))

    CALL dfbydx(d_store7, d_wd1x)
    CALL dfbydy(d_store7, d_wd1y)
    CALL dfbydz(d_store7, d_wd1z)
    CALL d2fdx2(d_store7, d_wd2x)
    CALL d2fdy2(d_store7, d_wd2y)
    CALL d2fdz2(d_store7, d_wd2z)

  END IF

    !   PRESSURE
    IF (flmixp) THEN

    rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
    CALL maths_kernel_eqa_host("A = log(B) - RHSCAL 755", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ))

    CALL dfbydx(d_store7, d_pd1x)
    CALL dfbydy(d_store7, d_pd1y)
    CALL dfbydz(d_store7, d_pd1z)
    CALL d2fdx2(d_store7, d_pd2x)
    CALL d2fdy2(d_store7, d_pd2y)
    CALL d2fdz2(d_store7, d_pd2z)

  END IF

    !   TEMPERATURE
    IF (flmixt) THEN
    !       TRANSP CONTAINS LN(T/TDIFGB)
    CALL dfbydx(d_transp, d_td1x)
    CALL dfbydy(d_transp, d_td1y)
    CALL dfbydz(d_transp, d_td1z)
    CALL d2fdx2(d_transp, d_td2x)
    CALL d2fdy2(d_transp, d_td2y)
    CALL d2fdz2(d_transp, d_td2z)

  END IF

    !   =========================================================================

    !   RUN THROUGH ALL SPECIES
    !   -----------------------
    !   RSC 08-AUG-2012 EVALUATE ALL SPECIES
    !   RSC 08-JUN-2015 REMOVE Nth SPECIES TREATMENT
    !   VM: CONSERVATIVE-PRIMITIVE FIX
    DO ispec = 1, nspec

    !   =======================================================================

    !       YRHS CONTAINS RHO Y: CONVERT TO Y
    !       Y IS PARALLEL
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

    CALL maths_kernel_eqt_host("A = A/B - RHSCAL 816", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

  END DO

    DO ispec = 1, nspec

    !       =======================================================================

    !       Y EQUATION: CONVECTIVE TERMS
    !       ----------------------------
    !       HALF Y DIV RHO U

    !       COLLECT Y SOURCE TERMS IN RATE FOR NOW
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqap_host("A = A - half*B*C - RHSCAL 828", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_divm, 1, s3d_000, "real(kind=8)", OPS_READ))

    !                                                         RATE = Y SOURCE TERMS
    !                                                           VTMP = DIV CORR VEL
    !                                                       WTMP = MIXTURE H
    !       =======================================================================

    !       Y EQUATION: CONVECTIVE TERMS
    !       ----------------------------
    !       HALF DIV RHO U Y

    !       D/DX RHO U Y
    !       RHO U Y IS PARALLEL
    rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
    CALL maths_kernel_eqw_host("A = B*C", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ))

    CALL dfbydx(d_store7, d_store1)

    !       D/DY RHO V Y
    !       RHO V Y IS PARALLEL
    rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
    CALL maths_kernel_eqw_host("A = B*C - RHSCAL 855", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ))

    CALL dfbydy(d_store7, d_store2)

    !       D/DZ RHO W Y
    !       RHO W Y IS PARALLEL
    rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
    CALL maths_kernel_eqw_host("A = B*C - RHSCAL 865", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ))

    CALL dfbydz(d_store7, d_store3)

    !       COLLECT DIV RHO U Y IN RATE FOR NOW
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqaq_host("A = A - half*(B+C+D) - RHSCAL 874", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    !                                                  RATE = Y SOURCE TERMS
    !                                                           VTMP = DIV CORR VEL
    !                                                       WTMP = MIXTURE H
    !       =======================================================================

    !       SPECIES MASS FRACTION GRADIENT TERMS
    !       ------------------------------------

    !       SPECIES MASS FRACTION GRADIENTS
    !        rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    !        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz, &
    !                        ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
    !                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
    !
    !        call dfbydx(d_store7,d_store1)
    !        call dfbydy(d_store7,d_store2)
    !        call dfbydz(d_store7,d_store3)
    CALL dfbydx(d_yrhs(ispec), d_store1)
    CALL dfbydy(d_yrhs(ispec), d_store2)
    CALL dfbydz(d_yrhs(ispec), d_store3)



      !                                                         STORE1,2,3 = DYDX,Y,Z
      !                                                         RATE = Y SOURCE TERMS
      !                                                           VTMP = DIV CORR VEL
      !                                                       WTMP = MIXTURE H
      !       =======================================================================

      !       COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs
      !       -------------------------------------------------------

      !       X-DIRECTION: DYDX
      IF (fxlcnv) THEN
      rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
      CALL boundary_kernel_mass_xdir_host("COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs - RHSCAL 915", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_t6bxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF
    IF (fxrcnv) THEN
      rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL boundary_kernel_mass_xdir_host("COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs - RHSCAL 924", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bclyxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_t6bxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF

      !       Y-DIRECTION: DYDY
      IF (fylcnv) THEN
      rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
      CALL boundary_kernel_mass_ydir_host("COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs - RHSCAL 935", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_stryyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bclyyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_t6byl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF
    IF (fyrcnv) THEN
      rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
      CALL boundary_kernel_mass_ydir_host("COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs - RHSCAL 944", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bclyyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_t6byr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF

      !       Z-DIRECTION: DYDZ
      IF (fzlcnv) THEN
      rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
      CALL boundary_kernel_mass_zdir_host("COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs - RHSCAL 955", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bclyzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_t6bzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF
    IF (fzrcnv) THEN
      rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
      CALL boundary_kernel_mass_zdir_host("COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs - RHSCAL 964", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bclyzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_t6bzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF
    !                                                         STORE1,2,3 = DYDX,Y,Z
    !                                                         RATE = Y SOURCE TERMS
    !                                                           VTMP = DIV CORR VEL
    !                                                       WTMP = MIXTURE H
    !       =======================================================================

    !       Y EQUATION: CONVECTIVE TERMS
    !       ----------------------------
    !       HALF RHO U.DEL Y

    !       COLLECT HALF RHO U.DEL Y IN RATE FOR NOW
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqat_host("A = A - half*(B*C+D*E+F*G) - RHSCAL 983", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_READ))

    !       -------------------------------------
    !       Y-EQUATION: CONVECTIVE TERMS COMPLETE
    !       -------------------------------------
    !                                                         STORE1,2,3 = DYDX,Y,Z
    !                                                         RATE = Y SOURCE TERMS
    !                                                           VTMP = DIV CORR VEL
    !                                                       WTMP = MIXTURE H
    !       =======================================================================

    !       Y-EQUATION: DIFFUSIVE TERMS
    !       ---------------------------
    !       E-EQUATION: FURTHER HEAT FLUX TERMS


    !       MASS DIFFUSIVITY FOR Y
    !       ----------------------
    !       ANALYTICAL FUNCTION OF TEMPERATURE
    !       TRANSP CONTAINS CONDUCTIVITY/CP
    !       STORE DIFFUSIVITY IN STORE7 FOR NOW
    !       Y DIFFUSIVITY IS PARALLEL
    rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
    CALL maths_kernel_eqh_host("A = B*var(indx) - RHSCAL 1013", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(olewis, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

      !                                                         STORE1,2,3 = DYDX,Y,Z
      !                                                          STORE7 = DIFFUSIVITY
      !                                                         RATE = Y SOURCE TERMS
      !                                                           VTMP = DIV CORR VEL
      !                                                       WTMP = MIXTURE H
      !       -----------------------------------------------------------------------

      !       MIXTURE AVERAGED TRANSPORT
      !       RSC 17-APR-2013
      !       TRANSP CONTAINS LN(T/TDIFGB)
      IF (flmavt) THEN

      !           MASS DIFFUSIVITY FOR EACH SPECIES
      !           RELATIVE TO CURRENT SPECIES
      !#if defined(OPS_LAZY) || defined(OPS_WITH_CUDAFOR) || defined(OPS_WITH_OMPOFFLOADFOR)
      !
      !            rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
      !            call ops_par_loop(set_zero_kernel, "set_zero - RHSCAL 1034", senga_grid, 3, rangexyz, &
      !                            ops_arg_dat(d_combo1, 1, s3d_000, "real(kind=8)", OPS_WRITE))
      !            call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
      !                            ops_arg_dat(d_combo2, 1, s3d_000, "real(kind=8)", OPS_WRITE))
      !
      !            DO jspec = 1, nspec
      !!               COMBINATION RULE FOR MASS DIFFUSIVITY
      !                call ops_par_loop(maths_kernel_eqBF, "MASS DIFFUSIVITY FOR EACH SPECIES - part 1 - RHSCAL 1041", senga_grid, 3, rangexyz, &
      !                                ops_arg_dat(d_ctrans(jspec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
      !                                ops_arg_dat(d_combo1, 1, s3d_000, "real(kind=8)", OPS_RW), &
      !                                ops_arg_dat(d_combo2, 1, s3d_000, "real(kind=8)", OPS_RW), &
      !                                ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                                ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                                ops_arg_dat(d_yrhs(jspec), 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                                ops_arg_gbl(diffco, ndcfmx*nspcmx*nspcmx, "real(kind=8)", OPS_READ), &
      !                                ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
      !                                ops_arg_gbl(pdifgb, 1, "real(kind=8)", OPS_READ), &
      !                                ops_arg_gbl(dfctol_ops, 1, "real(kind=8)", OPS_READ), &
      !                                ops_arg_gbl(ncodif, 1, "integer(kind=4)", OPS_READ), &
      !                                ops_arg_gbl(ncodm1, 1, "integer(kind=4)", OPS_READ), &
      !                                ops_arg_gbl(jspec, 1, "integer(kind=4)", OPS_READ), &
      !                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
      !            END DO
      !
      !            call ops_par_loop(maths_kernel_eqBG, "MASS DIFFUSIVITY FOR EACH SPECIES - part 2 - RHSCAL 1058", senga_grid, 3, rangexyz, &
      !                            ops_arg_dat(d_combo1, 1, s3d_000, "real(kind=8)", OPS_RW), &
      !                            ops_arg_dat(d_combo2, 1, s3d_000, "real(kind=8)", OPS_RW), &
      !                            ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
      !                            ops_arg_dat(d_ctrans(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
      !                            ops_arg_gbl(dfctol_ops, 1, "real(kind=8)", OPS_READ), &
      !                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
      !
      !#else

      rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
      DO jspec = 1, nspec
        CALL copy_kernel_sdim_to_mdim_host("A_multidim(ispec) = B", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs_mdim, 22, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_yrhs(jspec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(jspec, 1, "integer(kind=4)", OPS_READ))
      END DO

      CALL maths_kernel_eqbfg_host("MASS DIFFUSIVITY FOR EACH SPECIES - RHSCAL 1175", senga_grid, 3, rangexyz, &
ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_yrhs_mdim, 22, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(diffco, ndcfmx * nspcmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(pdifgb, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(dfctol_ops, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ncodif, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ncodm1, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

      !#endif

    END IF

      !       -----------------------------------------------------------------------

      !       MIXTURE AVERAGED TRANSPORT
      !       RSC 17-APR-2013
      !       TRANSP CONTAINS LN(T/TDIFGB)
      IF (flmtdr(ispec)) THEN

      !           THERMAL DIFFUSION RATIO FOR EACH SPECIES
      !           RELATIVE TO CURRENT SPECIES

      rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
      CALL set_zero_kernel_host("set_zero", senga_grid, 3, rangexyz, &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        DO jspec = 1, nspec

        flmtds = flmtdr(jspec) .AND. (ispec /= jspec)
        IF (flmtds) THEN

          rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
          CALL maths_kernel_eqbh_host("THERMAL DIFFUSION RATIO FOR EACH SPECIES - RHSCAL 1093", senga_grid, 3, rangexyz, &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_yrhs(jspec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(tdrcco, ndcfmx * nspcmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(tdifgb, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ncotdr, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(jspec, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ncotm1, 1, "integer(kind=4)", OPS_READ))

        END IF

      END DO

      !           VM WITH HELP FROM RSC: CORRECTION IN THERMAL DIFFUSION TO ADD MOLE-FRACTION OF
      !           ORIGINAL SPECIES
      !           TDRCCO IS THERMAL DIFFUSION FACTOR (MASON ET AL. 1968)
      !           NEED TO MULTIPLY BY X_ALPHA AS WELL
      rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
      CALL maths_kernel_eqbs_host("maths_kernel_eqBS", senga_grid, 3, rangexyz, &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END IF

    !       =======================================================================

    !       DIFFUSION CORRECTION VELOCITY
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    !        call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 1223", senga_grid, 3, rangexyz,  &
    !                        ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_INC), &
    !                        ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

    !        call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 1228", senga_grid, 3, rangexyz,  &
    !                        ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
    !                        ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ))

    !        call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 1233", senga_grid, 3, rangexyz,  &
    !                        ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
    !                        ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                        ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    CALL maths_kernel_eqaa_fused_host("A = A+B*C - RHSCAL 1238", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    !                                                         STORE1,2,3 = DYDX,Y,Z
    !                                                          STORE7 = DIFFUSIVITY
    !                                                         RATE = Y SOURCE TERMS
    !                                                           VTMP = DIV CORR VEL
    !                                                       WTMP = MIXTURE H
    !       =======================================================================

    !       SPECIES ENTHALPY
    !       ----------------

    !       TEMPERATURE INTERVAL INDEXING
    iindex = 1 + (ispec - 1) / nspimx
    ipower = ispec - (iindex - 1) * nspimx - 1
    icoef2 = ntbase ** ipower
    icoef1 = icoef2 * ntbase

    !       SPECIES H IS PARALLEL
    !       STORE SPECIES H IN UTMP FOR NOW
    !       STORE MIXTURE H IN WTMP FOR NOW
    rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
    CALL maths_kernel_eqbl_host("SPECIES ENTHALPY - RHSCAL 1151", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(amasch, ncofmx * ntinmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ncpoly, ntinmx * nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ncpom1, ntinmx * nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ncenth, ntinmx * nspcmx, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

      !                                                         STORE1,2,3 = DYDX,Y,Z
      !                                                          STORE7 = DIFFUSIVITY
      !                                                  RATE = Y SOURCE TERMS
      !                                                       UTMP = SPECIES H
      !                                                           VTMP = DIV CORR VEL
      !                                                       WTMP = MIXTURE H
      !       =======================================================================

      !       COLLECT SPECIES H FOR BCs
      !       -------------------------

      !       X-DIRECTION
      IF (fxlcnv) THEN
      rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
      CALL boundary_kernel_speciesh_xdir_host("COLLECT SPECIES H FOR BCs - RHSCAL 1179", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF
    IF (fxrcnv) THEN
      rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL boundary_kernel_speciesh_xdir_host("COLLECT SPECIES H FOR BCs - RHSCAL 1186", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF

      !       Y-DIRECTION
      IF (fylcnv) THEN
      rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
      CALL boundary_kernel_speciesh_ydir_host("COLLECT SPECIES H FOR BCs - RHSCAL 1195", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF
    IF (fyrcnv) THEN
      rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
      CALL boundary_kernel_speciesh_ydir_host("COLLECT SPECIES H FOR BCs - RHSCAL 1202", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF

      !       Z-DIRECTION
      IF (fzlcnv) THEN
      rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
      CALL boundary_kernel_speciesh_zdir_host("COLLECT SPECIES H FOR BCs - RHSCAL 1211", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF
    IF (fzrcnv) THEN
      rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
      CALL boundary_kernel_speciesh_zdir_host("COLLECT SPECIES H FOR BCs - RHSCAL 1218", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF

      !                                                         STORE1,2,3 = DYDX,Y,Z
      !                                                          STORE7 = DIFFUSIVITY
      !                                                  RATE = Y SOURCE TERMS
      !                                                       UTMP = SPECIES H
      !                                                           VTMP = DIV CORR VEL
      !                                                       WTMP = MIXTURE H
      !       =======================================================================

      !       MIXTURE AVERAGED TRANSPORT
      !       RSC 23-APR-2013
      !       ADD DUFOUR EFFECT TERMS TO SPECIES ENTHALPY
      IF (flmduf(ispec)) THEN
      rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
      CALL maths_kernel_eql_host("A = A+var(indx)*B*C - RHSCAL 1237", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END IF

    !       =======================================================================

    !       Y EQUATION: DIFFUSIVE TERMS
    !       ---------------------------
    !       E EQUATION: FURTHER HEAT FLUX TERMS
    !       DIFFUSIVITY GRADIENT TERMS

    !       DIFFUSIVITY GRADIENTS
    CALL dfbydx(d_store7, d_store4)
    CALL dfbydy(d_store7, d_store5)
    CALL dfbydz(d_store7, d_store6)

    !       BOUNDARY CONDITIONS
    !       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fxldif) CALL zeroxl(d_store4)
    IF (fxrdif) CALL zeroxr(d_store4)
    !       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fyldif) CALL zeroyl(d_store5)
    IF (fyrdif) CALL zeroyr(d_store5)
    !       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fzldif) CALL zerozl(d_store6)
    IF (fzrdif) CALL zerozr(d_store6)

    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqba_host("maths_kernel_eqBA - RHSCAL 1270", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ))

    !       BOUNDARY CONDITIONS
    !       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF (fxladb) CALL zeroxl(d_store4)
    IF (fxradb) CALL zeroxr(d_store4)

    !       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF (fyladb) CALL zeroyl(d_store5)
    IF (fyradb) CALL zeroyr(d_store5)

    !       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF (fzladb) CALL zerozl(d_store6)
    IF (fzradb) CALL zerozr(d_store6)

    !       E EQUATION
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqai_host("A = A+(B*C+D*E+F*G)*H - RHSCAL 1295", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ))

    !                                                         STORE1,2,3 = DYDX,Y,Z
    !                                                          STORE7 = DIFFUSIVITY
    !                                                         RATE = Y SOURCE TERMS
    !                                                              UTMP = SPECIES H
    !                                                           VTMP = DIV CORR VEL
    !                                                              WTMP = MIXTURE H
    !       =======================================================================

    !       E-EQUATION: FURTHER HEAT FLUX TERMS
    !       -----------------------------------
    !       SPECIES ENTHALPY GRADIENT TERMS

    !       SPECIES ENTHALPY GRADIENTS
    CALL dfbydx(d_utmp, d_store4)
    CALL dfbydy(d_utmp, d_store5)
    CALL dfbydz(d_utmp, d_store6)

    !       BOUNDARY CONDITIONS
    !       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fxldif) CALL zeroxl(d_store4)
    IF (fxrdif) CALL zeroxr(d_store4)

    !       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fyldif) CALL zeroyl(d_store5)
    IF (fyrdif) CALL zeroyr(d_store5)

    !       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fzldif) CALL zerozl(d_store6)
    IF (fzrdif) CALL zerozr(d_store6)

    !       BOUNDARY CONDITIONS
    !       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF (fxladb) CALL zeroxl(d_store4)
    IF (fxradb) CALL zeroxr(d_store4)

    !       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF (fyladb) CALL zeroyl(d_store5)
    IF (fyradb) CALL zeroyr(d_store5)

    !       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF (fzladb) CALL zerozl(d_store6)
    IF (fzradb) CALL zerozr(d_store6)

    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqai_host("A = A+(B*C+D*E+F*G)*H - RHSCAL 1349", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

      !                                                         STORE1,2,3 = DYDX,Y,Z
      !                                                          STORE7 = DIFFUSIVITY
      !                                                         RATE = Y SOURCE TERMS
      !                                                              UTMP = SPECIES H
      !                                                           VTMP = DIV CORR VEL
      !                                                              WTMP = MIXTURE H
      !       =======================================================================

      !       Y-EQUATION: DIFFUSIVE TERMS
      !       ---------------------------
      !       WALL BC: MASS DIFFUSION TERMS
      !       E-EQUATION: HEAT FLUX TERMS
      !       WALL BC: ENTHALPY DIFFUSION TERMS
      IF (fxldfw) THEN
      rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
      CALL hf_kernel_eqs_host("hf_kernel_eqS - RHSCAL 1374", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

    END IF
    IF (fxrdfw) THEN
      rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL hf_kernel_eqt_host("hf_kernel_eqT - RHSCAL 1386", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

    END IF
    IF (fyldfw) THEN
      rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
      CALL hf_kernel_equ_host("hf_kernel_eqU - RHSCAL 1398", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

    END IF
    IF (fyrdfw) THEN
      rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
      CALL hf_kernel_eqv_host("hf_kernel_eqV - RHSCAL 1410", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

    END IF
    IF (fzldfw) THEN
      rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
      CALL hf_kernel_eqw_host("hf_kernel_eqW - RHSCAL 1422", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

    END IF
    IF (fzrdfw) THEN
      rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
      CALL hf_kernel_eqx_host("hf_kernel_eqX - RHSCAL 1434", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

    END IF

    !       =======================================================================

    !       Y-EQUATION: DIFFUSIVE TERMS
    !       ---------------------------
    !       E-EQUATION: FURTHER HEAT FLUX TERMS
    !       SECOND DERIVATIVE TERMS

    !       SPECIES MASS FRACTION SECOND DERIVATIVES
    !       MOVE DIFFUSIVITY TO STORE4
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL copy_kernel_host("copy - RHSCAL 1455", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

    !       MOVE MASS FRACTION TO STORE7
    !        rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    !        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
    !                        ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
    !                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
    !
    !        call d2fdx2(d_store7,d_store1)
    !        call d2fdy2(d_store7,d_store2)
    !        call d2fdz2(d_store7,d_store3)
    CALL d2fdx2(d_yrhs(ispec), d_store1)
    CALL d2fdy2(d_yrhs(ispec), d_store2)
    CALL d2fdz2(d_yrhs(ispec), d_store3)


    !       BOUNDARY CONDITIONS
    !       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fxldif) CALL zeroxl(d_store1)
    IF (fxrdif) CALL zeroxr(d_store1)

    !       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fyldif) CALL zeroyl(d_store2)
    IF (fyrdif) CALL zeroyr(d_store2)

    !       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fzldif) CALL zerozl(d_store3)
    IF (fzrdif) CALL zerozr(d_store3)

    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqbb_host("maths_kernel_eqBB - RHSCAL 1487", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ))

    !       BOUNDARY CONDITIONS
    !       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF (fxladb) CALL zeroxl(d_store1)
    IF (fxradb) CALL zeroxr(d_store1)

    !       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF (fyladb) CALL zeroyl(d_store2)
    IF (fyradb) CALL zeroyr(d_store2)

    !       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF (fzladb) CALL zerozl(d_store3)
    IF (fzradb) CALL zerozr(d_store3)

    !       E EQUATION
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqag_host("A = A+(B+C+D)*E*F - RHSCAL 1510", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ))

      !                                                         RATE = Y SOURCE TERMS
      !                                                           VTMP = DIV CORR VEL
      !                                                              WTMP = MIXTURE H
      !       =======================================================================

      !       MIXTURE AVERAGED TRANSPORT
      !       RSC 23-APR-2013
      !       MOLAR MASS TERMS, PRESSURE TERMS, SORET EFFECT

      !       MIXTURE MOLAR MASS TERMS
      IF (flmixw) THEN
      !           FIRST AND SECOND DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED
      rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
      CALL maths_kernel_eqw_host("A = B*C - RHSCAL 1531", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))

      !           DIFFUSION CORRECTION VELOCITY
      !           FIRST DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED
      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      !            call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 1646", senga_grid, 3, rangexyz,  &
      !                            ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_INC), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_wd1x, 1, s3d_000, "real(kind=8)", OPS_READ))

      !            call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 1651", senga_grid, 3, rangexyz,  &
      !                            ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_wd1y, 1, s3d_000, "real(kind=8)", OPS_READ))

      !            call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 1656", senga_grid, 3, rangexyz,  &
      !                            ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_wd1z, 1, s3d_000, "real(kind=8)", OPS_READ))

      CALL maths_kernel_eqaa_fused_host("A = A+B*C - RHSCAL 1661", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1z, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           Y EQUATION: DIFFUSIVE TERMS
      !           E EQUATION: FURTHER HEAT FLUX TERMS

      !           DIFFUSIVITY GRADIENT TERMS

      !           DIFFUSIVITY GRADIENTS
      CALL dfbydx(d_store7, d_store1)
      CALL dfbydy(d_store7, d_store2)
      CALL dfbydz(d_store7, d_store3)

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fxldif) CALL zeroxl(d_store1)
      IF (fxrdif) CALL zeroxr(d_store1)

      !           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fyldif) CALL zeroyl(d_store2)
      IF (fyrdif) CALL zeroyr(d_store2)

      !           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fzldif) CALL zerozl(d_store3)
      IF (fzrdif) CALL zerozr(d_store3)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqba_host("maths_kernel_eqBA - RHSCAL 1578", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1z, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fxladb) CALL zeroxl(d_store1)
      IF (fxradb) CALL zeroxr(d_store1)

      !           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fyladb) CALL zeroyl(d_store2)
      IF (fyradb) CALL zeroyr(d_store2)

      !           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fzladb) CALL zerozl(d_store3)
      IF (fzradb) CALL zerozr(d_store3)

      !           E EQUATION
      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqai_host("A = A+(B*C+D*E+F*G)*H - RHSCAL 1603", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           E-EQUATION: FURTHER HEAT FLUX TERMS
      !           SPECIES ENTHALPY GRADIENT TERMS

      !           SPECIES ENTHALPY GRADIENTS ALREADY IN STORE4,5,6

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fxldif) CALL zeroxl(d_store4)
      IF (fxrdif) CALL zeroxr(d_store4)

      !           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fyldif) CALL zeroyl(d_store5)
      IF (fyrdif) CALL zeroyr(d_store5)

      !           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fzldif) CALL zerozl(d_store6)
      IF (fzrdif) CALL zerozr(d_store6)

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fxladb) CALL zeroxl(d_store4)
      IF (fxradb) CALL zeroxr(d_store4)

      !           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fyladb) CALL zeroyl(d_store5)
      IF (fyradb) CALL zeroyr(d_store5)

      !           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fzladb) CALL zerozl(d_store6)
      IF (fzradb) CALL zerozr(d_store6)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqai_host("A = A+(B*C+D*E+F*G)*H - RHSCAL 1654", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

        !           ====================================================================

        !           Y-EQUATION: DIFFUSIVE TERMS
        !           ---------------------------
        !           WALL BC: MOLAR MASS TERMS
        !           E-EQUATION: HEAT FLUX TERMS
        !           WALL BC: ENTHALPY DIFFUSION TERMS
        IF (fxldfw) THEN
        rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
        CALL hf_kernel_eqs_host("hf_kernel_eqS - RHSCAL 1664", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1x, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fxrdfw) THEN
        rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
        CALL hf_kernel_eqt_host("hf_kernel_eqT - RHSCAL 1676", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1x, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fyldfw) THEN
        rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
        CALL hf_kernel_equ_host("hf_kernel_eqU - RHSCAL 1688", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1y, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fyrdfw) THEN
        rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
        CALL hf_kernel_eqv_host("hf_kernel_eqV - RHSCAL 1700", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1y, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fzldfw) THEN
        rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
        CALL hf_kernel_eqw_host("hf_kernel_eqW - RHSCAL 1712", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1z, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fzrdfw) THEN
        rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
        CALL hf_kernel_eqx_host("hf_kernel_eqX - RHSCAL 1724", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd1z, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

      END IF

      !           ====================================================================

      !           Y-EQUATION: DIFFUSIVE TERMS
      !           E-EQUATION: FURTHER HEAT FLUX TERMS
      !           SECOND DERIVATIVE TERMS
      !           SECOND DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fxldif) CALL zeroxl(d_wd2x)
      IF (fxrdif) CALL zeroxr(d_wd2x)

      !           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fyldif) CALL zeroyl(d_wd2y)
      IF (fyrdif) CALL zeroyr(d_wd2y)

      !           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fzldif) CALL zerozl(d_wd2z)
      IF (fzrdif) CALL zerozr(d_wd2z)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqbb_host("maths_kernel_eqBB - RHSCAL 1756", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_wd2x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd2y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd2z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fxladb) CALL zeroxl(d_wd2x)
      IF (fxradb) CALL zeroxr(d_wd2x)

      !           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fyladb) CALL zeroyl(d_wd2y)
      IF (fyradb) CALL zeroyr(d_wd2y)

      !           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fzladb) CALL zerozl(d_wd2z)
      IF (fzradb) CALL zerozr(d_wd2z)

      !           E EQUATION
      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqag_host("A = A+(B+C+D)*E*F - RHSCAL 1779", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_wd2x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd2y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wd2z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF
    !flmixw
      !       MIXTURE MOLAR MASS TERMS

      !       =======================================================================

      !       PRESSURE DIFFUSION TERMS
      IF (flmixp) THEN
      !           FIRST AND SECOND DERIVATIVES OF LN(PRESSURE) ALREADY STORED

      rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
      CALL maths_kernel_eqo_host("A = B*C*(one-var(ispec)/D) - RHSCAL 1797", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wmomix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(wmolar, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

      !           DIFFUSION CORRECTION VELOCITY
      !           FIRST DERIVATIVES OF LN(PRESSURE) ALREADY STORED
      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      !            call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 1924", senga_grid, 3, rangexyz,  &
      !                            ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_INC), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_pd1x, 1, s3d_000, "real(kind=8)", OPS_READ))

      !            call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 1929", senga_grid, 3, rangexyz,  &
      !                            ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_pd1y, 1, s3d_000, "real(kind=8)", OPS_READ))

      !            call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 1934", senga_grid, 3, rangexyz,  &
      !                            ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_pd1z, 1, s3d_000, "real(kind=8)", OPS_READ))

      CALL maths_kernel_eqaa_fused_host("A = A+B*C - RHSCAL 1939", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1z, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           Y EQUATION: DIFFUSIVE TERMS
      !           E EQUATION: FURTHER HEAT FLUX TERMS

      !           DIFFUSIVITY GRADIENT TERMS

      !           DIFFUSIVITY GRADIENTS
      CALL dfbydx(d_store7, d_store1)
      CALL dfbydy(d_store7, d_store2)
      CALL dfbydz(d_store7, d_store3)

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fxldif) CALL zeroxl(d_store1)
      IF (fxrdif) CALL zeroxr(d_store1)

      !           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fyldif) CALL zeroyl(d_store2)
      IF (fyrdif) CALL zeroyr(d_store2)

      !           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fzldif) CALL zerozl(d_store3)
      IF (fzrdif) CALL zerozr(d_store3)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqba_host("maths_kernel_eqBA - RHSCAL 1874", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1z, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fxladb) CALL zeroxl(d_store1)
      IF (fxradb) CALL zeroxr(d_store1)

      !           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fyladb) CALL zeroyl(d_store2)
      IF (fyradb) CALL zeroyr(d_store2)

      !           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fzladb) CALL zerozl(d_store3)
      IF (fzradb) CALL zerozr(d_store3)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqai_host("A = A+(B*C+D*E+F*G)*H - RHSCAL 1871", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           E-EQUATION: FURTHER HEAT FLUX TERMS
      !           SPECIES ENTHALPY GRADIENT TERMS

      !           SPECIES ENTHALPY GRADIENTS ALREADY IN STORE4,5,6

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fxldif) CALL zeroxl(d_store4)
      IF (fxrdif) CALL zeroxr(d_store4)

      !           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fyldif) CALL zeroyl(d_store5)
      IF (fyrdif) CALL zeroyr(d_store5)

      !           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fzldif) CALL zerozl(d_store6)
      IF (fzrdif) CALL zerozr(d_store6)

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fxladb) CALL zeroxl(d_store4)
      IF (fxradb) CALL zeroxr(d_store4)

      !           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fyladb) CALL zeroyl(d_store5)
      IF (fyradb) CALL zeroyr(d_store5)

      !           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fzladb) CALL zerozl(d_store6)
      IF (fzradb) CALL zerozr(d_store6)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqai_host("A = A+(B*C+D*E+F*G)*H - RHSCAL 1913", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

        !           ====================================================================

        !           Y-EQUATION: DIFFUSIVE TERMS
        !           ---------------------------
        !           WALL BC: PRESSURE TERMS
        !           E-EQUATION: HEAT FLUX TERMS
        !           WALL BC: ENTHALPY DIFFUSION TERMS
        IF (fxldfw) THEN
        rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
        CALL hf_kernel_eqs_host("hf_kernel_eqS - RHSCAL 1932", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1x, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fxrdfw) THEN
        rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
        CALL hf_kernel_eqt_host("hf_kernel_eqT - RHSCAL 1944", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1x, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fyldfw) THEN
        rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
        CALL hf_kernel_equ_host("hf_kernel_eqU - RHSCAL 1956", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1y, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fyrdfw) THEN
        rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
        CALL hf_kernel_eqv_host("hf_kernel_eqV - RHSCAL 1968", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1y, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fzldfw) THEN
        rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
        CALL hf_kernel_eqw_host("hf_kernel_eqW - RHSCAL 1980", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1z, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fzrdfw) THEN
        rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
        CALL hf_kernel_eqx_host(" hf_kernel_eqX - RHSCAL 1992", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd1z, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

      END IF

      !           ====================================================================

      !           Y-EQUATION: DIFFUSIVE TERMS
      !           E-EQUATION: FURTHER HEAT FLUX TERMS
      !           SECOND DERIVATIVE TERMS
      !           SECOND DERIVATIVES OF LN(PRESSURE) ALREADY STORED

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fxldif) CALL zeroxl(d_pd2x)
      IF (fxrdif) CALL zeroxr(d_pd2x)

      !           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fyldif) CALL zeroyl(d_pd2y)
      IF (fyrdif) CALL zeroyr(d_pd2y)

      !           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fzldif) CALL zerozl(d_pd2z)
      IF (fzrdif) CALL zerozr(d_pd2z)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqbb_host("maths_kernel_eqBB - RHSCAL 2024", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_pd2x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd2y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd2z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fxladb) CALL zeroxl(d_pd2x)
      IF (fxradb) CALL zeroxr(d_pd2x)

      !           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fyladb) CALL zeroyl(d_pd2y)
      IF (fyradb) CALL zeroyr(d_pd2y)

      !           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fzladb) CALL zerozl(d_pd2z)
      IF (fzradb) CALL zerozr(d_pd2z)

      !           E EQUATION
      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqag_host("A = A+(B+C+D)*E*F - RHSCAL 2047", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_pd2x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd2y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_pd2z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF
    !       PRESSURE DIFFUSION TERMS

      !       =======================================================================

      !       SORET EFFECT (THERMAL DIFFUSION) TERMS
      IF (flmsor(ispec)) THEN
      !           FIRST AND SECOND DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED

      rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
      CALL maths_kernel_eqy_host("A = B*C*D - RHSCAL 2065", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           DIFFUSION CORRECTION VELOCITY
      !           FIRST DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED
      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      !            call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 2199", senga_grid, 3, rangexyz,  &
      !                            ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_INC), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_td1x, 1, s3d_000, "real(kind=8)", OPS_READ))

      !            call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 2204", senga_grid, 3, rangexyz,  &
      !                            ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_td1y, 1, s3d_000, "real(kind=8)", OPS_READ))

      !            call ops_par_loop(maths_kernel_eqAA, "A = A+B*C - RHSCAL 2209", senga_grid, 3, rangexyz,  &
      !                            ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
      !                            ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
      !                            ops_arg_dat(d_td1z, 1, s3d_000, "real(kind=8)", OPS_READ))

      CALL maths_kernel_eqaa_fused_host("A = A+B*C - RHSCAL 2214", senga_grid, 3, rangexyz, &
ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           Y EQUATION: DIFFUSIVE TERMS
      !           E EQUATION: FURTHER HEAT FLUX TERMS

      !           DIFFUSIVITY GRADIENT TERMS

      !           DIFFUSIVITY GRADIENTS
      CALL dfbydx(d_store7, d_store1)
      CALL dfbydy(d_store7, d_store2)
      CALL dfbydz(d_store7, d_store3)

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fxldif) CALL zeroxl(d_store1)
      IF (fxrdif) CALL zeroxr(d_store1)

      !           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fyldif) CALL zeroyl(d_store2)
      IF (fyrdif) CALL zeroyr(d_store2)

      !           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fzldif) CALL zerozl(d_store3)
      IF (fzrdif) CALL zerozr(d_store3)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqba_host("maths_kernel_eqBA - RHSCAL 2113", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000, "real(kind=8)", OPS_READ))

        !           SUBTRACT DUFOUR EFFECT TERMS TO RESTORE SPECIES ENTHALPY
        !           RSC 08-JUN-2015 BUG FIX
        IF (flmduf(ispec)) THEN
        rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
        CALL maths_kernel_eqm_host("A = A-var(indx)*B*C - RHSCAL 2127", senga_grid, 3, rangexyz, &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

      END IF

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fxladb) CALL zeroxl(d_store1)
      IF (fxradb) CALL zeroxr(d_store1)

      !           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fyladb) CALL zeroyl(d_store2)
      IF (fyradb) CALL zeroyr(d_store2)

      !           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fzladb) CALL zerozl(d_store3)
      IF (fzradb) CALL zerozr(d_store3)

      !           E EQUATION
      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqai_host("A = A+(B*C+D*E+F*G)*H - RHSCAL 2151", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           E-EQUATION: FURTHER HEAT FLUX TERMS
      !           SPECIES ENTHALPY GRADIENT TERMS

      !           EVALUATE SPECIES ENTHALPY GRADIENTS USING STORE4,5,6
      CALL dfbydx(d_utmp, d_store4)
      CALL dfbydy(d_utmp, d_store5)
      CALL dfbydz(d_utmp, d_store6)

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fxldif) CALL zeroxl(d_store4)
      IF (fxrdif) CALL zeroxr(d_store4)

      !           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fyldif) CALL zeroyl(d_store5)
      IF (fyrdif) CALL zeroyr(d_store5)

      !           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fzldif) CALL zerozl(d_store6)
      IF (fzrdif) CALL zerozr(d_store6)

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fxladb) CALL zeroxl(d_store4)
      IF (fxradb) CALL zeroxr(d_store4)

      !           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fyladb) CALL zeroyl(d_store5)
      IF (fyradb) CALL zeroyr(d_store5)

      !           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fzladb) CALL zerozl(d_store6)
      IF (fzradb) CALL zerozr(d_store6)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqai_host("A = A+(B*C+D*E+F*G)*H - RHSCAL 2196", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store5, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store6, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

        !           ====================================================================

        !           Y-EQUATION: DIFFUSIVE TERMS
        !           ---------------------------
        !           WALL BC: SORET EFFECT TERMS
        !           E-EQUATION: HEAT FLUX TERMS
        !           WALL BC: ENTHALPY DIFFUSION TERMS
        IF (fxldfw) THEN
        rgspec_ispec = rgspec(ispec)
        rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
        CALL hf_kernel_eqm_host("hf_kernel_eqM - RHSCAL 2215", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fxrdfw) THEN
        rgspec_ispec = rgspec(ispec)
        rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
        CALL hf_kernel_eqn_host("hf_kernel_eqN - RHSCAL 2231", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fyldfw) THEN
        rgspec_ispec = rgspec(ispec)
        rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
        CALL hf_kernel_eqo_host("hf_kernel_eqO - RHSCAL 2247", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fyrdfw) THEN
        rgspec_ispec = rgspec(ispec)
        rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
        CALL hf_kernel_eqp_host("hf_kernel_eqP - RHSCAL 2263", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fzldfw) THEN
        rgspec_ispec = rgspec(ispec)
        rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
        CALL hf_kernel_eqq_host("hf_kernel_eqQ - RHSCAL 2279", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ))

      END IF
      IF (fzrdfw) THEN
        rgspec_ispec = rgspec(ispec)
        rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
        CALL hf_kernel_eqr_host("hf_kernel_eqR - RHSCAL 2295", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ))

      END IF

        !           E-EQUATION: HEAT FLUX TERMS
        !           WALL BC: SORET AND DUFOUR TERMS
        IF (flmduf(ispec)) THEN
        !               E-EQUATION: HEAT FLUX TERMS
          !               WALL BC: ADIABATIC WALL
          IF (fxlcnw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
          CALL hf_kernel_eqg_host("hf_kernel_eqG - RHSCAL 2521", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))
        END IF
        IF (fxrcnw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
          CALL hf_kernel_eqh_host("hf_kernel_eqH - RHSCAL 2533", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))
        END IF
        IF (fylcnw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
          CALL hf_kernel_eqi_host("hf_kernel_eqI - RHSCAL 2545", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))
        END IF
        IF (fyrcnw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
          CALL hf_kernel_eqj_host("hf_kernel_eqJ - RHSCAL 2557", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

        END IF
        IF (fzlcnw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
          CALL hf_kernel_eqk_host("hf_kernel_eqK - RHSCAL 2570", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

        END IF
        IF (fzrcnw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
          CALL hf_kernel_eql_host("hf_kernel_eqL - RHSCAL 2583", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

        END IF

          !               E-EQUATION: HEAT FLUX TERMS
          !               WALL BC: ISOTHERMAL WALL
          IF (fxladw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
          CALL hf_kernel_eqg_host("hf_kernel_eqG - RHSCAL 2521", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

        END IF
        IF (fxradw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
          CALL hf_kernel_eqh_host("hf_kernel_eqH - RHSCAL 2533", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1x, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltax, 1, "real(kind=8)", OPS_READ))

        END IF
        IF (fyladw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
          CALL hf_kernel_eqi_host("hf_kernel_eqI - RHSCAL 2545", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

        END IF
        IF (fyradw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
          CALL hf_kernel_eqj_host("hf_kernel_eqJ - RHSCAL 2557", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1y, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltay, 1, "real(kind=8)", OPS_READ))

        END IF
        IF (fzladw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
          CALL hf_kernel_eqk_host("hf_kernel_eqK - RHSCAL 2570", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

        END IF
        IF (fzradw) THEN
          rgspec_ispec = rgspec(ispec)
          rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
          CALL hf_kernel_eql_host("hf_kernel_eqL - RHSCAL 2583", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_trun, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_tdrmix, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td1z, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltaz, 1, "real(kind=8)", OPS_READ))

        END IF

      END IF

      !           ====================================================================

      !           Y-EQUATION: DIFFUSIVE TERMS
      !           E-EQUATION: FURTHER HEAT FLUX TERMS
      !           SECOND DERIVATIVE TERMS
      !           SECOND DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fxldif) CALL zeroxl(d_td2x)
      IF (fxrdif) CALL zeroxr(d_td2x)

      !           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fyldif) CALL zeroyl(d_td2y)
      IF (fyrdif) CALL zeroyr(d_td2y)

      !           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
      IF (fzldif) CALL zerozl(d_td2z)
      IF (fzrdif) CALL zerozr(d_td2z)

      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqbb_host("maths_kernel_eqBB - RHSCAL 2498", senga_grid, 3, rangexyz, &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_td2x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td2y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td2z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ))

      !           BOUNDARY CONDITIONS
      !           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fxladb) CALL zeroxl(d_td2x)
      IF (fxradb) CALL zeroxr(d_td2x)

      !           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fyladb) CALL zeroyl(d_td2y)
      IF (fyradb) CALL zeroyr(d_td2y)

      !           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
      IF (fzladb) CALL zerozl(d_td2z)
      IF (fzradb) CALL zerozr(d_td2z)

      !           E EQUATION
      rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL maths_kernel_eqag_host("A = A+(B+C+D)*E*F - RHSCAL 2521", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_td2x, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td2y, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_td2z, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_utmp, 1, s3d_000, "real(kind=8)", OPS_READ))

    END IF

    !                                                      RATE = Y SOURCE TERMS
    !                                                       VTMP = DIV CORR VEL
    !                                                          WTMP = MIXTURE H
    !   =======================================================================

    !   ----------------------------------------------------------------
    !   E-EQUATION: DIFFUSION CORRECTION VELOCITY TERMS EVALUATED BELOW
    !   Y-EQUATION: DIFFUSION CORRECTION VELOCITY TERMS EVALUATED BELOW
    !   ----------------------------------------------------------------

  END DO
  ! END of ispec loop
  !   RSC 08-AUG-2012 EVALUATE ALL SPECIES
  !   END OF RUN THROUGH ALL SPECIES

  !   =========================================================================

  !   EVALUATE DIFFUSION CORRECTION VELOCITY TERMS
  !   --------------------------------------------

  !   E-EQUATION: FURTHER HEAT FLUX TERMS
  !   -----------------------------------
  !   MIXTURE ENTHALPY GRADIENTS
  CALL dfbydx(d_wtmp, d_store1)
  CALL dfbydy(d_wtmp, d_store2)
  CALL dfbydz(d_wtmp, d_store3)

  !   BOUNDARY CONDITIONS
  !   BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF (fxldif) CALL zeroxl(d_store1)
  IF (fxrdif) CALL zeroxr(d_store1)

  !   BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF (fyldif) CALL zeroyl(d_store2)
  IF (fyrdif) CALL zeroyr(d_store2)

  !   BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF (fzldif) CALL zerozl(d_store3)
  IF (fzrdif) CALL zerozr(d_store3)

  !   BOUNDARY CONDITIONS
  !   BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fxladb) CALL zeroxl(d_store1)
  IF (fxradb) CALL zeroxr(d_store1)

  !   BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fyladb) CALL zeroyl(d_store2)
  IF (fyradb) CALL zeroyr(d_store2)

  !   BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fzladb) CALL zerozl(d_store3)
  IF (fzradb) CALL zerozr(d_store3)

  !   TRANSFER DIV CORR VEL TO TEMPORARY STORE
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL copy_kernel_host("copy - RHSCAL 2585", senga_grid, 3, rangexyz, &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ))

  !   BOUNDARY CONDITIONS
  !   BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fxladb) CALL zeroxl(d_store4)
  IF (fxradb) CALL zeroxr(d_store4)

  !   BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fyladb) CALL zeroyl(d_store4)
  IF (fyradb) CALL zeroyr(d_store4)

  !   BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF (fzladb) CALL zerozl(d_store4)
  IF (fzradb) CALL zerozr(d_store4)

  !   DIV RHO VCORR HMIX
  rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
  CALL maths_kernel_eqam_host("A = A-B*C-D*E-F*G-H*I - RHSCAL 2604", senga_grid, 3, rangexyz, &
ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
ops_arg_dat(d_wtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store4, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_READ))

    !                                                       RATE = Y SOURCE TERMS
    !                                                         VTMP = DIV CORR VEL
    !   =========================================================================

    !   MIXTURE AVERAGED TRANSPORT
    !   EVALUATE THE VISCOSITY

    !   RSC 17-APR-2013
    !   TRANSP CONTAINS LN(T)
    !   STORE VISCOSITY IN DIFMIX FOR NOW
    IF (flmavt) THEN

    !#if defined(OPS_LAZY) || defined(OPS_WITH_CUDAFOR) || defined(OPS_WITH_OMPOFFLOADFOR)
    !
    !        rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    !
    !        DO ispec = 1, nspec
    !            call ops_par_loop(maths_kernel_eqBI, "STORE VISCOSITY IN DIFMIX - part 1 - RHSCAL 2630", senga_grid, 3, rangexyz, &
    !                            ops_arg_dat(d_ctrans(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
    !                            ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                            ops_arg_gbl(viscco, nvcfmx*nspcmx, "real(kind=8)", OPS_READ), &
    !                            ops_arg_gbl(ncovis, 1, "integer(kind=4)", OPS_READ), &
    !                            ops_arg_gbl(ncovm1, 1, "integer(kind=4)", OPS_READ), &
    !                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
    !       END DO
    !
    !        call ops_par_loop(set_zero_kernel, "set_zero - RHSCAL 2639", senga_grid, 3, rangexyz, &
    !                        ops_arg_dat(d_combo1, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    !
    !        DO ispec = 1, nspec
    !            call ops_par_loop(set_zero_kernel, "set_zero - RHSCAL 2643", senga_grid, 3, rangexyz, &
    !                            ops_arg_dat(d_combo2, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    !
    !            DO jspec = 1, nspec
    !                call ops_par_loop(maths_kernel_eqBJ, "STORE VISCOSITY IN DIFMIX - part 2 - RHSCAL 2647", senga_grid, 3, rangexyz, &
    !                                ops_arg_dat(d_combo2, 1, s3d_000, "real(kind=8)", OPS_RW), &
    !                                ops_arg_dat(d_ctrans(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                                ops_arg_dat(d_ctrans(jspec), 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                                ops_arg_dat(d_yrhs(jspec), 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                                ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
    !                                ops_arg_gbl(wilko1, nspcmx*nspcmx, "real(kind=8)", OPS_READ), &
    !                                ops_arg_gbl(wilko2, nspcmx*nspcmx, "real(kind=8)", OPS_READ), &
    !                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
    !                                ops_arg_gbl(jspec, 1, "integer(kind=4)", OPS_READ))
    !            END DO
    !
    !            call ops_par_loop(maths_kernel_eqBK, "STORE VISCOSITY IN DIFMIX - part 2 - RHSCAL 2659", senga_grid, 3, rangexyz, &
    !                            ops_arg_dat(d_combo1, 1, s3d_000, "real(kind=8)", OPS_RW), &
    !                            ops_arg_dat(d_ctrans(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                            ops_arg_dat(d_combo2, 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                            ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
    !                            ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
    !                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
    !
    !        END DO
    !
    !        call ops_par_loop(copy_kernel, "copy - RHSCAL 2669", senga_grid, 3, rangexyz, &
    !                        ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
    !                        ops_arg_dat(d_combo1, 1, s3d_000, "real(kind=8)", OPS_READ))
    !
    !#else
    !
    rangexyz = [1 - nhalox, nxglbl + nhalox, 1 - nhaloy, nyglbl + nhaloy, 1 - nhaloz, nzglbl + nhaloz]
    DO jspec = 1, nspec
      CALL copy_kernel_sdim_to_mdim_host("A_multidim(ispec) = B", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs_mdim, 22, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_yrhs(jspec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(jspec, 1, "integer(kind=4)", OPS_READ))
    END DO

    CALL maths_kernel_eqbijk_host("STORE VISCOSITY IN DIFMIX - RHSCAL 2849", senga_grid, 3, rangexyz, &
ops_arg_dat(d_difmix, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_transp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_yrhs_mdim, 22, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_gbl(viscco, nvcfmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(wilko1, nspcmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(wilko2, nspcmx * nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ovwmol, nspcmx, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ncovis, 1, "integer(kind=4)", OPS_READ), &
ops_arg_gbl(ncovm1, 1, "integer(kind=4)", OPS_READ))

    !#endif

  END IF

    !   =========================================================================

    !   RUN THROUGH ALL SPECIES
    !   -----------------------
    !   RSC 08-AUG-2012 EVALUATE ALL SPECIES
    !   RSC 08-JUN-2015 REMOVE Nth SPECIES TREATMENT
    DO ispec = 1, nspec
    !       Y-EQUATION: DIFFUSIVE TERMS
    !       ---------------------------
    !       RECOMPUTE SPECIES MASS FRACTION GRADIENTS
    !        rangexyz = [1-nhalox,nxglbl+nhalox,1-nhaloy,nyglbl+nhaloy,1-nhaloz,nzglbl+nhaloz]
    !        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz, &
    !                        ops_arg_dat(d_store7, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
    !                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
    !
    !        call dfbydx(d_store7,d_store1)
    !        call dfbydy(d_store7,d_store2)
    !        call dfbydz(d_store7,d_store3)
    CALL dfbydx(d_yrhs(ispec), d_store1)
    CALL dfbydy(d_yrhs(ispec), d_store2)
    CALL dfbydz(d_yrhs(ispec), d_store3)

    !       BOUNDARY CONDITIONS
    !       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fxldif) CALL zeroxl(d_store1)
    IF (fxrdif) CALL zeroxr(d_store1)

    !       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fyldif) CALL zeroyl(d_store2)
    IF (fyrdif) CALL zeroyr(d_store2)

    !       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF (fzldif) CALL zerozl(d_store3)
    IF (fzrdif) CALL zerozr(d_store3)

    !       DIV RHO VCORR Y
    !       STORE Y SOURCE TERMS IN YRHS
    rangexyz = [1, nxglbl, 1, nyglbl, 1, nzglbl]
    CALL maths_kernel_eqan_host("A = B-A*C-D*E-F*G-H*I - RHSCAL 2713", senga_grid, 3, rangexyz, &
ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_RW), &
ops_arg_dat(d_rate(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vtmp, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_ucor, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vcor, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_wcor, 1, s3d_000, "real(kind=8)", OPS_READ))

  END DO

    !   RSC 08-AUG-2012 EVALUATE ALL SPECIES
    !   END OF RUN THROUGH ALL SPECIES
    !                                                            ALL STORES CLEAR
    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !   ------------------------------------------------
    !   Y-EQUATION: SOURCE TERMS COMPLETE
    !   ------------------------------------------------
    !   E-EQUATION: PRESSURE-WORK AND VISCOUS WORK TERMS
    !               EVALUATED IN SUBROUTINE RHSVEL
    !   ------------------------------------------------

    !   =========================================================================
    !   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !   =========================================================================

    !   COLLECT DENSITY AND ITS GRADIENTS FOR BCs
    !   -----------------------------------------

    !   X-DIRECTION: DRHODX
    IF (fxlcnv .OR. fxrcnv) THEN

    CALL dfbydx(d_drhs, d_store1)

      IF (fxlcnv) THEN
      rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
      CALL boundary_kernel_density_xdir_host("COLLECT DENSITY AND ITS GRADIENTS FOR BCs - RHSCAL 2752", senga_grid, 3, rangexyz, &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF
    IF (fxrcnv) THEN
      rangexyz = [nxglbl, nxglbl, 1, nyglbl, 1, nzglbl]
      CALL boundary_kernel_density_xdir_host("COLLECT DENSITY AND ITS GRADIENTS FOR BCs - RHSCAL 2761", senga_grid, 3, rangexyz, &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

    END IF

  END IF

    !   Y-DIRECTION: DRHODY
    IF (fylcnv .OR. fyrcnv) THEN

    CALL dfbydy(d_drhs, d_store2)

      IF (fylcnv) THEN
      rangexyz = [1, nxglbl, 1, 1, 1, nzglbl]
      CALL boundary_kernel_density_ydir_host("COLLECT DENSITY AND ITS GRADIENTS FOR BCs - RHSCAL 2778", senga_grid, 3, rangexyz, &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF
    IF (fyrcnv) THEN
      rangexyz = [1, nxglbl, nyglbl, nyglbl, 1, nzglbl]
      CALL boundary_kernel_density_ydir_host("COLLECT DENSITY AND ITS GRADIENTS FOR BCs - RHSCAL 2787", senga_grid, 3, rangexyz, &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

    END IF

  END IF

    !   Z-DIRECTION: DRHODZ
    IF (fzlcnv .OR. fzrcnv) THEN

    CALL dfbydz(d_drhs, d_store3)

      IF (fzlcnv) THEN
      rangexyz = [1, nxglbl, 1, nyglbl, 1, 1]
      CALL boundary_kernel_density_zdir_host("COLLECT DENSITY AND ITS GRADIENTS FOR BCs - RHSCAL 2804", senga_grid, 3, rangexyz, &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF
    IF (fzrcnv) THEN
      rangexyz = [1, nxglbl, 1, nyglbl, nzglbl, nzglbl]
      CALL boundary_kernel_density_zdir_host("COLLECT DENSITY AND ITS GRADIENTS FOR BCs - RHSCAL 2813", senga_grid, 3, rangexyz, &
ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_store3, 1, s3d_000, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

    END IF

  END IF

  !   =========================================================================

END SUBROUTINE rhscal