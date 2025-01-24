SUBROUTINE bounds

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BOUNDS
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   13-JUL-2003:  RSC MODIFIED FOR SENGA2
!   08-AUG-2012:  RSC EVALUATE ALL SPECIES
!   26-OCT-2013:  RSC ACTIVATE ALL BCS ON ALL SIDES

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   COMPUTES CHARACTERISTIC BOUNDARY CONDITIONS FOR ALL PDES

!   *************************************************************************

!   VISHNU MOHAN (VM), KHALIL ABO AMSHA (KAA) AND NILANJAN CHAKRABORTY (NC)
!   IMPLEMENTED THE TRANSVERSE DERIVATIVE TERMS FOR BOUNDARY CONDITION AS PER METHOD
!   OF YOO AND IM


!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!     -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    integer(kind=4) :: ispec
    integer(kind=4) :: rangexyz(6)
    real(kind=8), parameter :: bet=0.050_8
    real(kind=8), target :: bet_ops=bet
!   VM: SWITCH TO SET BET=LOCAL_MACH
!   FLAG_BET=0 POINTS TO GLOBAL VALUE OF BET
!   FLAG_BET=1 POINTS TO LOCAL VALUE OF BET
    integer(kind=4) :: flag_bet_xl,flag_bet_xr,flag_bet_yl,flag_bet_yr
    integer(kind=4) :: flag_bet_zl,flag_bet_zr
!   VM: SWITCH TO TURN POINTWISE INFLOW-OUTFLOW ON/OFF
!   0 is OFF, 1 is ON
    integer(kind=4) :: flag_pio_xl,flag_pio_xr,flag_pio_yl,flag_pio_yr
    integer(kind=4) :: flag_pio_zl,flag_pio_zr

!   FY - FOR NR INFLOW
    real(kind=8), parameter :: m2max=0.001_8
    real(kind=8), target :: m2max_ops=m2max
    real(kind=8) :: nrieta1,nrieta2,nrieta3,nrieta4,nrieta5,nrieta6

!   BEGIN
!   =====

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   X-DIRECTION LEFT-HAND END
!   -------------------------
    IF(fxlcnv) THEN

!       =======================================================================

!       STR ARRAYS CONTAIN STORED VALUES
!       STRUXL = PRIMITIVE U-VELOCITY COMPONENT
!       STRVXL = PRIMITIVE V-VELOCITY COMPONENT
!       STRWXL = PRIMITIVE W-VELOCITY COMPONENT
!       STRPXL = PRESSURE
!       STRDXL = DENSITY
!       STRTXL = TEMPERATURE
!       STREXL = INTERNAL ENERGY
!       STRGXL = MIXTURE CP
!       STRRXL = MIXTURE SPECIFIC GAS CONSTANT
!       STRYXL(ISPEC) = SPECIES MASS FRACTION
!       RATEXL(ISPEC) = SPECIES REACTION RATE
!       STRHXL(ISPEC) = SPECIES ENTHALPY

!       BCL ARRAYS CONTAIN FIRST DERIVATIVES
!       BCL1XL = DUDX
!       BCL2XL = DRHODX
!       BCL3XL = DVDX
!       BCL4XL = DWDX
!       BCL5XL = DPDX
!       BCLYXL(ISPEC) = DYDX

!       =======================================================================

!       REDUCED SPECIES ENTHALPY
!       ------------------------
        DO ispec = 1,nspec
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqA_xdir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC),  &
                    ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

        END DO

!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqB_xdir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_ovgmxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!       SPEED OF SOUND
!       --------------
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqC_xdir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqAA_xdir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_tt1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_t1bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t2bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t51bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t52bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t4bxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        DO ispec = 1,nspec
            call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_tt6xl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_t6bxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
        END DO

!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbcxl == nsbco1) THEN

            flag_pio_xl = nxlprm(1) ! 0
            flag_bet_xl = nxlprm(2) ! 1

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L5X AS REQUIRED
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqH_xl, "OUTFLOW BC1 XL", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobcxl, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfxl, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(xgdlen, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(bet_ops, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_bet_xl, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_gbl(flag_pio_xl, 1, "integer(kind=4)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqP_xl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_pio_xl, 1, "integer(kind=4)", OPS_READ))

!           RSC 08-AUG-2012 EVALUATE ALL SPECIES
!           DO ISPEC = 1,NSPM1
            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqQ_xl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(flag_pio_xl, 1, "integer(kind=4)", OPS_READ))
            END DO

        END IF

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbcxl == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L2X-L5X
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqI_xl, "L2X to L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobcxl, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfxl, 1, "real(kind=8)", OPS_READ))

!           LYX
            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqJ_xl, "LYX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqR_xl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqS_xl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcxl == nsbci2) THEN

!           INFLOW BOUNDARY CONDITION No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqE_xdir, "A_yz = A_yz + B_mulditim_yz*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_dydtxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
                call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqG_xdir, "A_yz = A_yz/B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqF_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X,L2X,L5X
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqK_xl, "L1X L2X L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dtdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqT_xl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcxl == nsbci3)THEN

!           INFLOW BOUNDARY CONDITION No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           SPECIFY L's AS REQUIRED
!           L1X-L5X
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqL_xl, "L1X to L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dddtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqU_xl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqV_xl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_dydtxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================
        IF(nsbcxl == nsbci4)THEN

!           INFLOW BC No 4 - IMPLEMENTED BY FY
!           SUBSONIC NON-REFLECTING INFLOW AS PER LODATO WITH TRANSVERSE TERMS INCLUDED

!           ASSIGN NON-REFLECTING INFLOW COEFFICIENTS READ FROM CONTROL FILE
            nrieta1 = rxlprc(1)
            nrieta2 = rxlprc(2)
            nrieta3 = rxlprc(3)
            nrieta4 = rxlprc(4)
            nrieta5 = rxlprc(5)
            nrieta6 = rxlprc(6)

!           STORE VELOCITY, TEMPERATURE MASS FRACTION VALUES AT THE BOUNDARY
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(copy_kernel_xxdir, "copy - BOUNDS 511", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stluxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            call ops_par_loop(copy_kernel_xxdir, "copy - BOUNDS 514", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stlvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            call ops_par_loop(copy_kernel_xxdir, "copy - BOUNDS 517", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stlwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            call ops_par_loop(copy_kernel_xxdir, "copy - BOUNDS 520", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stltxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            DO ispec = 1,nspec
                call ops_par_loop(copy_kernel_xxdir, "copy - BOUNDS 524", senga_grid, 3, rangexyz, &
                                ops_arg_dat(d_stlyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                                ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            END DO

!           SET VELOCITY, TEMPERATURE, MASS FRACTION TARGET VALUES
            call bcutxl
            call bcttxl
            call bcytxl

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAB_xdir, "bounds_kernel_eqAB_xdir", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), & 
                            ops_arg_dat(d_tt2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_stluxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_stlvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_stlwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_stltxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(xgdlen, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(nrieta2, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(nrieta3, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(nrieta4, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(nrieta5, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(m2max, 1, "real(kind=8)", OPS_READ))

!           LYX
            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqAC_xdir, "bounds_kernel_eqAC_xdir", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_stlyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt6xl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_stluxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(xgdlen, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(nrieta6, 1, "real(kind=8)", OPS_READ))
            END DO

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAD_xdir, "bounds_kernel_eqAD_xdir", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_stluxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_stlvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_stlwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqAE_xdir, "bounds_kernel_eqAE_xdir", senga_grid, 3, rangexyz, &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stlyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            END DO
        END IF

!       =======================================================================

!       WALL BOUNDARY CONDITIONS
!       ------------------------

        IF(nsbcxl == nsbcw1) THEN

!           WALL BOUNDARY CONDITION No 1
!           NO-SLIP WALL - ADIABATIC

!           VELOCITY AND TEMPERATURE IMPOSED AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_xl, "bounds_kernel_eqAF_xl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_RW))

!           PRECOMPUTE CHEMISTRY TERMS 
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(bounds_kernel_eqAG_xdir, "bounds_kernel_eqAG_xdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            END DO

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAH_xdir, "bounds_kernel_eqAH_xdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X-L5X
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqM_xl, "bounds_kernel_eqM_xl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           LYX: OLD VALUE OF LYX
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqW_xl, "bounds_kernel_eqW_xl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqX_xl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcxl == nsbcw2) THEN

!           WALL BOUNDARY CONDITION No 2
!           NO-SLIP WALL - ISOTHERMAL

!           VELOCITY AND TEMPERATURE IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratexl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X-L5X
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqN_xl, "bounds_kernel_eqN_xl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           LYX
            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqO_xl, "bounds_kernel_eqO_xl", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqY_xl, "bounds_kernel_eqY_xl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqZ_xl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

    END IF
!   X-DIRECTION LEFT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   X-DIRECTION RIGHT-HAND END
!   --------------------------
    IF(fxrcnv) THEN

!       =======================================================================

!       STR ARRAYS CONTAIN STORED VALUES
!       STRUXR = PRIMITIVE U-VELOCITY COMPONENT
!       STRVXR = PRIMITIVE V-VELOCITY COMPONENT
!       STRWXR = PRIMITIVE W-VELOCITY COMPONENT
!       STRPXR = PRESSURE
!       STRDXR = DENSITY
!       STRTXR = TEMPERATURE
!       STREXR = INTERNAL ENERGY
!       STRGXR = MIXTURE CP
!       STRRXR = MIXTURE SPECIFIC GAS CONSTANT
!       STRYXR(ISPEC) = SPECIES MASS FRACTION
!       RATEXR(ISPEC) = SPECIES REACTION RATE
!       STRHXR(ISPEC) = SPECIES ENTHALPY

!       BCL ARRAYS CONTAIN FIRST DERIVATIVES
!       BCL1XR = DUDX
!       BCL2XR = DRHODX
!       BCL3XR = DVDX
!       BCL4XR = DWDX
!       BCL5XR = DPDX
!       BCLYXR(ISPEC) = DYDX

!       =======================================================================

!       REDUCED SPECIES ENTHALPY
!       ------------------------
        DO ispec = 1,nspec
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqA_xdir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

        END DO

!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqB_xdir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_ovgmxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!       SPEED OF SOUND
!       --------------
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqC_xdir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqAA_xdir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_tt1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_t1bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t2bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t51bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t52bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t4bxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        DO ispec = 1,nspec
            call ops_par_loop(copy_kernel_xxdir, "copy", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_tt6xr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_t6bxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
        END DO


!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbcxr == nsbco1) THEN

            flag_pio_xr = nxrprm(1) ! 0
            flag_bet_xr = nxrprm(2) ! 1

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L1X AS REQUIRED
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqH_xr, "OUTFLOW BC1 XR", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobcxr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfxr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(xgdlen, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(bet_ops, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_bet_xr, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_gbl(flag_pio_xr, 1, "integer(kind=4)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqP_xr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_pio_xr, 1, "integer(kind=4)", OPS_READ))

!           RSC 08-AUG-2012 EVALUATE ALL SPECIES
!           DO ISPEC = 1,NSPM1
            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqQ_xr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bclyxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(flag_pio_xr, 1, "integer(kind=4)", OPS_READ))
            END DO

        END IF

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbcxr == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X-L4X
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqI_xr, "L1X-L4X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobcxr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfxr, 1, "real(kind=8)", OPS_READ))

!           LYX
            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqJ_xr, "LYX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqR_xr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqS_xr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcxr == nsbci2) THEN

!           INFLOW BOUNDARY CONDITION No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqE_xdir, "A_yz = A_yz + B_mulditim_yz*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_dydtxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
                call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqG_xdir, "A_yz = A_yz/B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqF_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X,L2X,L5X
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqK_xr, "L1X L2X L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dudtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dtdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqT_xr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcxr == nsbci3) THEN

!           INFLOW BOUNDARY CONDITION No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           SPECIFY L's AS REQUIRED
!           L1X-L5X
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqL_xr, "L1X to L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dudtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dddtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dvdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dwdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqU_xr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqV_xr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_dydtxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

!       WALL BOUNDARY CONDITIONS
!       ------------------------

        IF(nsbcxr == nsbcw1) THEN

!           WALL BOUNDARY CONDITION No 1
!           NO-SLIP WALL - ADIABATIC

!           VELOCITY AND TEMPERATURE IMPOSED AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_xr, "bounds_kernel_eqAF_xr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_RW))

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(bounds_kernel_eqAG_xdir, "bounds_kernel_eqAG_xdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
            END DO

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAH_xdir, "bounds_kernel_eqAH_xdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X,L3X-L5X
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqM_xr, "bounds_kernel_eqM_xr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           LYX: OLD VALUE OF LYX
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqW_xr, "bounds_kernel_eqW_xr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqX_xr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcxr == nsbcw2) THEN

!           WALL BOUNDARY CONDITION No 2
!           NO-SLIP WALL - ISOTHERMAL

!           VELOCITY AND TEMPERATURE IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratexr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X-L5X
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqN_xr, "bounds_kernel_eqN_xr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           LYX
            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqO_xr, "bounds_kernel_eqO_xr", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqY_xr, "bounds_kernel_eqY_xr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqZ_xr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2xr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

    END IF
!   X-DIRECTION RIGHT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Y-DIRECTION LEFT-HAND END
!   -------------------------
    IF(fylcnv) THEN

!       =======================================================================

!       STR ARRAYS CONTAIN STORED VALUES
!       STRUYL = PRIMITIVE U-VELOCITY COMPONENT
!       STRVYL = PRIMITIVE V-VELOCITY COMPONENT
!       STRWYL = PRIMITIVE W-VELOCITY COMPONENT
!       STRPYL = PRESSURE
!       STRDYL = DENSITY
!       STRTYL = TEMPERATURE
!       STREYL = INTERNAL ENERGY
!       STRGYL = MIXTURE CP
!       STRRYL = MIXTURE SPECIFIC GAS CONSTANT
!       STRYYL(ISPEC) = SPECIES MASS FRACTION
!       RATEYL(ISPEC) = SPECIES REACTION RATE
!       STRHYL(ISPEC) = SPECIES ENTHALPY

!       BCL ARRAYS CONTAIN FIRST DERIVATIVES
!       BCL1YL = DVDY
!       BCL2YL = DRHODY
!       BCL3YL = DUDY
!       BCL4YL = DWDY
!       BCL5YL = DPDY
!       BCLYYL(ISPEC) = DYDY

!       =======================================================================

!       REDUCED SPECIES ENTHALPY
!       ------------------------
        DO ispec = 1,nspec
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqA_ydir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

        END DO

!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqB_ydir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_ovgmyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!       SPEED OF SOUND
!       --------------
        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqC_ydir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqAA_ydir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_tt1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_t1byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t2byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t51byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t52byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(copy_kernel_yydir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        call ops_par_loop(copy_kernel_yydir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t4byl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,1,1,nzglbl]
        DO ispec = 1,nspec
            call ops_par_loop(copy_kernel_yydir, "copy", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_tt6yl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_t6byl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))
        END DO


!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbcyl == nsbco1) THEN

            flag_pio_yl = nylprm(1) ! 0
            flag_bet_yl = nylprm(2) ! 1

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1, nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L5Y AS REQUIRED
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqH_yl, "OUTFLOW BC1 YL", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobcyl, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfyl, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ygdlen, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(bet_ops, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_bet_yl, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_gbl(flag_pio_yl, 1, "integer(kind=4)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqP_yl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_pio_yl, 1, "integer(kind=4)", OPS_READ))

!           RSC 08-AUG-2012 EVALUATE ALL SPECIES
!           DO ISPEC = 1,NSPM1
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqQ_yl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bclyyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(flag_pio_yl, 1, "integer(kind=4)", OPS_READ))
            END DO

        END IF

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbcyl == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L2Y-L5Y
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqI_yl, "L2Y to L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobcyl, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfyl, 1, "real(kind=8)", OPS_READ))

!           LYY
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqJ_yl, "LYY", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqR_yl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqS_yl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcyl == nsbci2) THEN

!           INFLOW BOUNDARY CONDITION No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqE_ydir, "A_xz = A_xz + B_mulditim_xz*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_dydtyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
                call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqG_ydir, "A_xz = A_xz/B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqF_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y,L2Y,L5Y
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqK_yl, "L1Y L2Y L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dtdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqT_yl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcyl == nsbci3)THEN

!           INFLOW BOUNDARY CONDITION No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           SPECIFY L's AS REQUIRED
!           L1Y-L5Y
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqL_yl, "L1Y to L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dudtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dddtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dwdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqU_yl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqV_yl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_dydtyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

!       WALL BOUNDARY CONDITIONS
!       ------------------------

        IF(nsbcyl == nsbcw1) THEN

!           WALL BOUNDARY CONDITION No 1
!           NO-SLIP WALL - ADIABATIC

!           VELOCITY AND TEMPERATURE IMPOSED AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_yl, "bounds_kernel_eqAF_yl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_RW))

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(bounds_kernel_eqAG_ydir, "bounds_kernel_eqAG_ydir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))
            END DO

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAH_ydir, "bounds_kernel_eqAH_ydir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y,L3Y-L5Y
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqM_yl, "bounds_kernel_eqM_yl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           LYY: OLD VALUE OF LYY
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqW_yl, "bounds_kernel_eqW_yl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqX_yl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcyl == nsbcw2) THEN

!           WALL BOUNDARY CONDITION No 2
!           NO-SLIP WALL - ISOTHERMAL

!           VELOCITY AND TEMPERATURE IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_rateyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y-L5Y
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqN_yl, "bounds_kernel_eqN_yl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           LYY
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqO_yl, "bounds_kernel_eqO_yl", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqY_yl, "bounds_kernel_eqY_yl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqZ_yl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

    END IF
!   Y-DIRECTION LEFT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Y-DIRECTION RIGHT-HAND END
!   --------------------------
    IF(fyrcnv) THEN

!       =======================================================================

!       STR ARRAYS CONTAIN STORED VALUES
!       STRUYR = PRIMITIVE U-VELOCITY COMPONENT
!       STRVYR = PRIMITIVE V-VELOCITY COMPONENT
!       STRWYR = PRIMITIVE W-VELOCITY COMPONENT
!       STRPYR = PRESSURE
!       STRDYR = DENSITY
!       STRTYR = TEMPERATURE
!       STREYR = INTERNAL ENERGY
!       STRGYR = MIXTURE CP
!       STRRYR = MIXTURE SPECIFIC GAS CONSTANT
!       STRYYR(ISPEC) = SPECIES MASS FRACTION
!       RATEYR(ISPEC) = SPECIES REACTION RATE
!       STRHYR(ISPEC) = SPECIES ENTHALPY

!       BCL ARRAYS CONTAIN FIRST DERIVATIVES
!       BCL1YR = DVDY
!       BCL2YR = DRHODY
!       BCL3YR = DUDY
!       BCL4YR = DWDY
!       BCL5YR = DPDY
!       BCLYYR(ISPEC) = DYDY

!       =======================================================================

!       REDUCED SPECIES ENTHALPY
!       ------------------------
        DO ispec = 1,nspec
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqA_ydir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

        END DO

!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqB_ydir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_ovgmyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!       SPEED OF SOUND
!       --------------
        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqC_ydir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(bounds_kernel_eqAA_ydir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_tt1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_t1byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t2byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t51byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t52byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel_yydir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel_yydir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t4byr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
        DO ispec = 1,nspec
            call ops_par_loop(copy_kernel_yydir, "copy", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_tt6yr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_t6byr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))
        END DO


!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbcyr == nsbco1) THEN

            flag_pio_yr = nyrprm(1) ! 0
            flag_bet_yr = nyrprm(2) ! 1

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L1Y AS REQUIRED
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqH_yr, "OUTFLOW BC1 YR", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobcyr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfyr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ygdlen, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(bet_ops, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_bet_yr, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_gbl(flag_pio_yr, 1, "integer(kind=4)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqP_yr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_pio_yr, 1, "integer(kind=4)", OPS_READ))

!           RSC 08-AUG-2012 EVALUATE ALL SPECIES
!           DO ISPEC = 1,NSPM1
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqQ_yr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bclyyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(flag_pio_yr, 1, "integer(kind=4)", OPS_READ))
            END DO

        END IF

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbcyr == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y-L4Y
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqI_yr, "L1Y-L4Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobcyr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfyr, 1, "real(kind=8)", OPS_READ))

!           LYY
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqJ_yr, "LYY", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqR_yr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqS_yr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcyr == nsbci2) THEN

!           INFLOW BOUNDARY CONDITION No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqE_ydir, "A_xz = A_xz + B_mulditim_xz*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_dydtyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
                call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqG_ydir, "A_xz = A_xz/B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqF_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y,L2Y,L5Y
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqK_yr, "L1Y L2Y L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dvdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dtdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqT_yr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcyr == nsbci3) THEN

!           INFLOW BOUNDARY CONDITION No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           SPECIFY L's AS REQUIRED
!           L1Y-L5Y
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqL_yr, "L1Y to L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dudtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dddtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dvdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dwdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqU_yr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqV_yr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_dydtyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

!       WALL BOUNDARY CONDITIONS
!       ------------------------

        IF(nsbcyr == nsbcw1)THEN

!           WALL BOUNDARY CONDITION No 1
!           NO-SLIP WALL - ADIABATIC

!           VELOCITY AND TEMPERATURE IMPOSED AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_yr, "bounds_kernel_eqAF_yr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_RW))

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(bounds_kernel_eqAG_ydir, "bounds_kernel_eqAG_ydir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))
            END DO

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAH_ydir, "bounds_kernel_eqAH_ydir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y,L3Y-L5Y
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqM_yr, "bounds_kernel_eqM_yr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           LYY: OLD VALUE OF LYY
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))
            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqW_yr, "bounds_kernel_eqW_yr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqX_yr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcyr == nsbcw2) THEN

!           WALL BOUNDARY CONDITION No 2
!           NO-SLIP WALL - ISOTHERMAL

!           VELOCITY AND TEMPERATURE IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_rateyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y-L5Y
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqN_yr, "bounds_kernel_eqN_yr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           LYY
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqO_yr, "bounds_kernel_eqO_yr", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqY_yr, "bounds_kernel_eqY_yr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bounds_kernel_eqZ_yr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2yr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

    END IF
!   Y-DIRECTION RIGHT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Z-DIRECTION LEFT-HAND END
!   -------------------------
    IF(fzlcnv) THEN

!       =======================================================================

!       STR ARRAYS CONTAIN STORED VALUES
!       STRUZL = PRIMITIVE U-VELOCITY COMPONENT
!       STRVZL = PRIMITIVE V-VELOCITY COMPONENT
!       STRWZL = PRIMITIVE W-VELOCITY COMPONENT
!       STRPZL = PRESSURE
!       STRDZL = DENSITY
!       STRTZL = TEMPERATURE
!       STREZL = INTERNAL ENERGY
!       STRGZL = MIXTURE CP
!       STRRZL = MIXTURE SPECIFIC GAS CONSTANT
!       STRYZL(ISPEC) = SPECIES MASS FRACTION
!       RATEZL(ISPEC) = SPECIES REACTION RATE
!       STRHZL(ISPEC) = SPECIES ENTHALPY

!       BCL ARRAYS CONTAIN FIRST DERIVATIVES
!       BCL1ZL = DWDZ
!       BCL2ZL = DRHODZ
!       BCL3ZL = DUDZ
!       BCL4ZL = DVDZ
!       BCL5ZL = DPDZ
!       BCLYZL(ISPEC) = DYDZ

!       =======================================================================

!       REDUCED SPECIES ENTHALPY
!       ------------------------
        DO ispec = 1,nspec
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqA_zdir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

        END DO

!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(bounds_kernel_eqB_zdir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_ovgmzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!       SPEED OF SOUND
!       --------------
        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(bounds_kernel_eqC_zdir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(bounds_kernel_eqAA_zdir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_tt1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_t1bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t2bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t51bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t52bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(copy_kernel_zzdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        call ops_par_loop(copy_kernel_zzdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t4bzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,nyglbl,1,1]
        DO ispec = 1,nspec
            call ops_par_loop(copy_kernel_zzdir, "copy", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_tt6zl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_t6bzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))
        END DO

!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbczl == nsbco1) THEN

            flag_pio_zl = nzlprm(1) ! 0
            flag_bet_zl = nzlprm(2) ! 1

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqF_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L5Z AS REQUIRED
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqH_zl, "OUTFLOW BC1 ZL", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobczl, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfzl, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(zgdlen, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(bet_ops, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_bet_zl, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_gbl(flag_pio_zl, 1, "integer(kind=4)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqP_zl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_pio_zl, 1, "integer(kind=4)", OPS_READ))

!           RSC 08-AUG-2012 EVALUATE ALL SPECIES
!           DO ISPEC = 1,NSPM1
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqQ_zl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bclyzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(flag_pio_zl, 1, "integer(kind=4)", OPS_READ))
            END DO

        END IF

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbczl == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqF_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L2Z-L5Z
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqI_zl, "L2Z to L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobczl, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfzl, 1, "real(kind=8)", OPS_READ))

!           LYZ
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqJ_zl, "LYZ", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqR_zl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqS_zl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbczl == nsbci2) THEN

!           INFLOW BOUNDARY CONDITION No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqE_zdir, "A_xy = A_xy + B_mulditim_xy*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_dydtzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
                call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqG_zdir, "A_xy = A_xy/B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqF_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z,L2Z,L5Z
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqK_zl, "L1Z L2Z L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dtdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqT_zl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbczl == nsbci3) THEN

!           INFLOW BOUNDARY CONDITION No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           SPECIFY L's AS REQUIRED
!           L1Z-L5Z
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqL_zl, "L1Z to L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dudtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dddtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dvdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqU_zl, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqV_zl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_dydtzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

!       WALL BOUNDARY CONDITIONS
!       ------------------------

        IF(nsbczl == nsbcw1) THEN

!           WALL BOUNDARY CONDITION No 1
!           NO-SLIP WALL - ADIABATIC

!           VELOCITY AND TEMPERATURE IMPOSED AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqAF_zl, "bounds_kernel_eqAF_zl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_RW))

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            DO ispec = 1,nspec
                call ops_par_loop(bounds_kernel_eqAG_zdir, "bounds_kernel_eqAG_zdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))
            END DO

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqAH_zdir, "bounds_kernel_eqAH_zdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z,L3Z-L5Z
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqM_zl, "bounds_kernel_eqM_zl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           LYZ: OLD VALUE OF LYZ
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            DO ispec = 1,nspec
                call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqW_zl, "bounds_kernel_eqW_zl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), & 
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqX_zl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbczl == nsbcw2) THEN

!           WALL BOUNDARY CONDITION No 2
!           NO-SLIP WALL - ISOTHERMAL

!           VELOCITY AND TEMPERATURE IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratezl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqF_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z-L5Z
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqN_zl, "bounds_kernel_eqN_zl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           LYZ
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqO_zl, "bounds_kernel_eqO_zl", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqY_zl, "bounds_kernel_eqY_zl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bounds_kernel_eqZ_zl, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

    END IF
!   Z-DIRECTION LEFT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Z-DIRECTION RIGHT-HAND END
!   --------------------------
    IF(fzrcnv) THEN

!       =======================================================================

!       STR ARRAYS CONTAIN STORED VALUES
!       STRUZR = PRIMITIVE U-VELOCITY COMPONENT
!       STRVZR = PRIMITIVE V-VELOCITY COMPONENT
!       STRWZR = PRIMITIVE W-VELOCITY COMPONENT
!       STRPZR = PRESSURE
!       STRDZR = DENSITY
!       STRTZR = TEMPERATURE
!       STREZR = INTERNAL ENERGY
!       STRGZR = MIXTURE CP
!       STRRZR = MIXTURE SPECIFIC GAS CONSTANT
!       STRYZR(ISPEC) = SPECIES MASS FRACTION
!       RATEZR(ISPEC) = SPECIES REACTION RATE
!       STRHZR(ISPEC) = SPECIES ENTHALPY

!       BCL ARRAYS CONTAIN FIRST DERIVATIVES
!       BCL1ZR = DWDR
!       BCL2ZR = DRHODR
!       BCL3ZR = DUDR
!       BCL4ZR = DVDR
!       BCL5ZR = DPDR
!       BCLYZR(ISPEC) = DYDR

!       =======================================================================

!       REDUCED SPECIES ENTHALPY
!       ------------------------
        DO ispec = 1,nspec
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqA_zdir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

        END DO

!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(bounds_kernel_eqB_zdir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC),  &
                        ops_arg_dat(d_ovgmzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!       SPEED OF SOUND
!       --------------
        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(bounds_kernel_eqC_zdir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(bounds_kernel_eqAA_zdir, "SPEED OF SOUND", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_tt1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_tt5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_t1bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t2bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t51bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_t52bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(copy_kernel_zzdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t3bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        call ops_par_loop(copy_kernel_zzdir, "copy", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_tt4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_t4bzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
        DO ispec = 1,nspec
            call ops_par_loop(copy_kernel_zzdir, "copy", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_tt6zr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_t6bzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))
        END DO

!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbczr == nsbco1) THEN

            flag_pio_zr = nzrprm(1) ! 0
            flag_bet_zr = nzrprm(2) ! 1

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L1Z AS REQUIRED
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqH_zr, "OUTFLOW BC1 ZR", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_tt1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobczr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfzr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(zgdlen, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(bet_ops, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_bet_zr, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_gbl(flag_pio_zr, 1, "integer(kind=4)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqP_zr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(flag_pio_zr, 1, "integer(kind=4)", OPS_READ))

!           RSC 08-AUG-2012 EVALUATE ALL SPECIES
!           DO ISPEC = 1,NSPM1
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqQ_zr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bclyzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(flag_pio_zr, 1, "integer(kind=4)", OPS_READ))
            END DO

        END IF

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbczr == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z-L4Z
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqI_zr, "L1Z-L4Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(cobczr, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(pinfzr, 1, "real(kind=8)", OPS_READ))

!           LYZ
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqJ_zr, "LYZ", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqR_zr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqS_zr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbczr == nsbci2) THEN

!           INFLOW BOUNDARY CONDITION No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqE_zdir, "A_xy = A_xy + B_mulditim_xy*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_dydtzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
                call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqG_zdir, "A_xy = A_xy/B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqF_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z,L2Z,L5Z
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqK_zr, "L1Z L2Z L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dwdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dtdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqT_zr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbczr == nsbci3) THEN

!           INFLOW BOUNDARY CONDITION No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           SPECIFY L's AS REQUIRED
!           L1Z-L5Z
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqL_zr, "L1Z to L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dudtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dddtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dvdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_dwdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqU_zr, "ADD TO CONSERVATIVE SOURCE TERMS", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ovgmzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqV_zr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_bclyzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                                ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_dydtzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

!       WALL BOUNDARY CONDITIONS
!       ------------------------

        IF(nsbczr == nsbcw1) THEN

!           WALL BOUNDARY CONDITION No 1
!           NO-SLIP WALL - ADIABATIC

!           VELOCITY AND TEMPERATURE IMPOSED AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_zr, "bounds_kernel_eqAF_zr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_RW))

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(bounds_kernel_eqAG_zdir, "bounds_kernel_eqAG_zdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))
            END DO

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqAH_zdir, "bounds_kernel_eqAH_zdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z,L3Z-L5Z
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqM_zr, "bounds_kernel_eqM_zr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           LYZ:  OLD VALUE OF LYZ
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            DO ispec = 1,nspec
                call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))
            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqW_zr, "bounds_kernel_eqW_zr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqX_zr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbczr == nsbcw2) THEN

!           WALL BOUNDARY CONDITION No 2
!           NO-SLIP WALL - ISOTHERMAL

!           VELOCITY AND TEMPERATURE IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_strhzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ratezr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqF_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z-L5Z
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqN_zr, "bounds_kernel_eqN_zr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl3zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl4zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_bcl5zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_acouzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           LYZ
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqO_zr, "bounds_kernel_eqO_zr", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE))

            END DO

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqY_zr, "bounds_kernel_eqY_zr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bounds_kernel_eqZ_zr, "EVALUATE ALL SPECIES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_INC), &
                                ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl1zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_ova2zr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

    END IF
!   Z-DIRECTION RIGHT-HAND END

!   =========================================================================

END SUBROUTINE bounds
