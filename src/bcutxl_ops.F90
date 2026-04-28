
! Auto-generated at 2026-04-28 18:43:09.271854 by ops-translator
SUBROUTINE bcutxl

  USE ops_fortran_declarations
  USE ops_fortran_rt_support
  USE bcut_kernel_xdir_eqf_module
  USE bcut_kernel_xdir_eqg_module
  USE bcut_kernel_xdir_eqh_module
  USE bcut_kernel_xdir_eqi_module

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   BCUTXL
  !   ======

  !   AUTHOR
  !   ------
  !   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

  !   CHANGE RECORD
  !   -------------
  !   30-DEC-2003:  CREATED
  !   04-JAN-2007:  RSC REVISE PARALLEL RECEIVES

  !   DESCRIPTION
  !   -----------
  !   DNS CODE SENGA2
  !   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
  !   AND THEIR TIME DERIVATIVES

  !   X-DIRECTION LEFT-HAND END

  !   *************************************************************************

  !   GLOBAL DATA
  !   ===========
  !   -------------------------------------------------------------------------
  !   -------------------------------------------------------------------------

  !   LOCAL DATA
  !   ==========
  !KA   FIX INFLOW BUG, BTIME IS DEFINED IN COM_SENGA2.H
  !KA      real(kind=8)) BTIME
  REAL(KIND = 8) :: fornow, argmnt, argval, realkx
  REAL(KIND = 8) :: cosval, sinval, costht, sintht
  REAL(KIND = 8) :: pcount
  INTEGER(KIND = 4) :: ic, jc, kc
  INTEGER(KIND = 4) :: iic, iim, kx, kxbase
  INTEGER(KIND = 4) :: icproc, ncount, irproc, irtag
  INTEGER(KIND = 4) :: rangexyz(6)
  REAL(KIND = 8) :: init_val1, init_val2
  REAL(KIND = 8) :: rxlprm_1, rxlprc_7

  !   VM: SYNTHETIC DIGITAL FILTERING METHOD
  REAL(KIND = 8) :: vfac, dvfdt, coflow

  !   FY - FOR NON-REFLECTING INLOW (INFLOW OPTION 4)
  REAL(KIND = 8) :: deltagy, ycoord
  INTEGER(KIND = 4) :: igofsty, iy
  REAL(KIND = 8) :: lambda
  ! SOUND WAVELENGTH
  REAL(KIND = 8) :: pulrat
  ! RATIO OF PULSE WIDTH TO SOUND WAVELENGTH
  REAL(KIND = 8) :: ptly
  ! BASE WIDTH OF PULSE
  REAL(KIND = 8) :: widthp
  ! ADDITIONAL WIDTH PARAMETER - DEFAULT 0.5
  REAL(KIND = 8) :: slope
  ! LARGER VALUE RESULTS IN SHARPER SLOPE

  !   BEGIN
  !   =====

  !   lambda = 3.74_8 * 10.0_8**-3
  lambda = 3.74_8 * (1.0_8 / (10.0_8 * 10.0_8 * 10.0_8))
  pulrat = 0.1_8
  ptly = pulrat * lambda
  widthp = 0.5_8
  slope = 2.0_8 * 10.0_8 ** 4

    !   =========================================================================

    !KA   THIS WAS MOVED TO BOUNDT & BOUNTT TO FIX INFLOW SCANNING LOCATION
    !     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
    !KA      BTIME = ETIME + RKTIM(IRKSTP)

    !   =========================================================================

    !   CONSTANT U-VELOCITY
    !   PARAMETER I1=1, R1=U-VELOCITY
    IF (nxlprm(1) == 1) THEN
    rxlprm_1 = rxlprm(1)
    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL bcut_kernel_xdir_eqf_host("bcut_kernel_xdir_eqF", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(rxlprm_1, 1, "real(kind=8)", OPS_READ))

  END IF

    !   =========================================================================

    !   SINUSOIDAL U-VELOCITY
    !   PARAMETER I1=2, R1=AMPLITUDE, R2=PERIOD
    IF (nxlprm(1) == 2) THEN

    fornow = two * pi / rxlprc(8)
    argmnt = fornow * (etime + rxlprc(9) * rxlprc(8))

      !       NXLPRM(4) = 0 - OPTION FOR STANDARD FLAT SINUSOIDAL INFLOW VELOCITY IN TIME ON WHOLE XL FACE
      IF (nxlprm(4) == 0) THEN

      init_val1 = rxlprm(1) + rxlprc(7) * SIN(argmnt)
      init_val2 = fornow * rxlprc(7) * COS(argmnt)

      rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
      CALL bcut_kernel_xdir_eqg_host("bcut_kernel_xdir_eqG", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(init_val1, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(init_val2, 1, "real(kind=8)", OPS_READ))

      !       NXLPRM(4) = 4 - OPTION FOR SINUSOIDAL VELOCITY IN TIME FOR PART OF XL FACE - ACTS AS A POINT SOURCE
    ELSE IF (nxlprm(4) == 4) THEN

      deltagy = ygdlen / REAL(nyglbl - 1, kind = 8)
      rxlprm_1 = rxlprm(1)
      rxlprc_7 = rxlprc(7)

      rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
      CALL bcut_kernel_xdir_eqh_host("bcut_kernel_xdir_eqH", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_gbl(slope, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(widthp, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ptly, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rxlprm_1, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rxlprc_7, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(fornow, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(argmnt, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(deltagy, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(ygdlen, 1, "real(kind=8)", OPS_READ), &
ops_arg_idx())

    END IF
  END IF

    !   =========================================================================

    !   GENERATING TURBULENT FIELD USING SYNTHETIC DIGITAL FILTERING
    !   METHOD
    !   VM: NXLPRM(1)=4 IMPLIES THAT THE VELOCITY SYTHETIC DIGITAL FILTERING
    !   IS ON

    IF (nxlprm(1) == 4) THEN

    vfac = one
    coflow = rxlprm(4)

      IF (vfac < 1.0_8) THEN
      dvfdt = 10000.0_8
    ELSE
      dvfdt = 0.0_8
    END IF

    rxlprm_1 = rxlprm(1)

    rangexyz = [1, 1, 1, nyglbl, 1, nzglbl]
    CALL bcut_kernel_xdir_eqi_host("bcut_kernel_xdir_eqI", senga_grid, 3, rangexyz, &
ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
ops_arg_dat(d_ustead, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_uinf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vinf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_winf1, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_uinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_vinf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_dat(d_winf2, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
ops_arg_gbl(rxlprm_1, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(vfac, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(coflow, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(tstep, 1, "real(kind=8)", OPS_READ), &
ops_arg_gbl(dvfdt, 1, "real(kind=8)", OPS_READ))

  END IF

  !   =========================================================================

END SUBROUTINE bcutxl