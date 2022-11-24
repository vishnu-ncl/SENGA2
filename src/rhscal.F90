SUBROUTINE rhscal

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

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
    real(kind=dp) :: ctrans(nspcmx)
    real(kind=dp) :: fornow,combo1,combo2,combo3
    integer :: ic,jc,kc,ispec
    integer :: itint,icp,iindex,ipower,icoef1,icoef2
    logical :: flmtds
    integer :: rangexyz(6)

!    BEGIN
!    =====

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   EVALUATE THE TEMPERATURE
!   ------------------------
!   ALSO PRESSURE, MIXTURE CP AND MIXTURE GAS CONSTANT
    call temper

!                                                             PRUN,TRUN = P,T
!                                                         STORE7 = RHO*MIX RG
!   =========================================================================

!   COLLECT MIXTURE CP AND GAS CONSTANT FOR BCs
!   -------------------------------------------

!   X-DIRECTION
    IF(fxlcnv)THEN
        rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
        call ops_par_loop(boundary_kernel_CPandGAS_xdir, "COLLECT CP AND GAS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

    END IF
    IF(fxrcnv)THEN
        rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(boundary_kernel_CPandGAS_xdir, "COLLECT CP AND GAS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

    END IF

!   Y-DIRECTION
    IF(fylcnv)THEN
        rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
        call ops_par_loop(boundary_kernel_CPandGAS_ydir, "COLLECT CP AND GAS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

    END IF
    IF(fyrcnv)THEN
        rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
        call ops_par_loop(boundary_kernel_CPandGAS_ydir, "COLLECT CP AND GAS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

    END IF

!   Z-DIRECTION
    IF(fzlcnv)THEN
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
        call ops_par_loop(boundary_kernel_CPandGAS_zdir, "COLLECT CP AND GAS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

    END IF
    IF(fzrcnv)THEN
        rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
        call ops_par_loop(boundary_kernel_CPandGAS_zdir, "COLLECT CP AND GAS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

    END IF

!                                                            ALL STORES CLEAR
!   =========================================================================

!   MASS FLUX DIVERGENCE
!   --------------------
!   URHS,VRHS,WRHS CONTAIN RHO U, RHO V, RHO W

    call dfbydx(d_urhs,d_store1)
    call dfbydy(d_vrhs,d_store2)
    call dfbydz(d_wrhs,d_store3)

    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    
    call ops_par_loop(math_kernel_eqL, "A=B+C+D", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_divm, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ)) 

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
    rangexyz = (/istalt,istolt,jstalt,jstolt,kstalt,kstolt/)
    call ops_par_loop(math_kernel_eqS, "A=A/B", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ))

!                                                            ALL STORES CLEAR
!   =========================================================================

!   COLLECT INTERNAL ENERGY FOR BCs
!   -------------------------------

!   X-DIRECTION
    IF(fxlcnv)THEN
        rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
        call ops_par_loop(boundary_kernel_internalenergy_xdir, "COLLECT INTERNAL ENERGY FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

    END IF
    IF(fxrcnv)THEN
        rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(boundary_kernel_internalenergy_xdir, "COLLECT INTERNAL ENERGY FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

    END IF

!   Y-DIRECTION
    IF(fylcnv)THEN
        rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
        call ops_par_loop(boundary_kernel_internalenergy_ydir, "COLLECT INTERNAL ENERGY FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

    END IF
    IF(fyrcnv)THEN
        rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
        call ops_par_loop(boundary_kernel_internalenergy_ydir, "COLLECT INTERNAL ENERGY FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

    END IF

!   Z-DIRECTION
    IF(fzlcnv)THEN
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
        call ops_par_loop(boundary_kernel_internalenergy_zdir, "COLLECT INTERNAL ENERGY FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

    END IF
    IF(fzrcnv)THEN
        rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
        call ops_par_loop(boundary_kernel_internalenergy_zdir, "COLLECT INTERNAL ENERGY FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

    END IF

!                                                            ALL STORES CLEAR
!   =========================================================================

!   E EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   HALF E DIV RHO U

!   COLLECT E DIV RHO U IN STORE4 FOR NOW
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(math_kernel_eqV, "A=B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_divm, 1, s3d_000, "real(dp)", OPS_READ))

!                                                        STORE4 = E DIV RHO U
!   =========================================================================

!   E EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   HALF DIV RHO U E

!   D/DX RHO U E
!   RHO U E IS PARALLEL
    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(math_kernel_eqV, "A=B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_READ))

    call dfbydx(d_store7,d_store1)

!   D/DY RHO V E
!   RHO V E IS PARALLEL
    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(math_kernel_eqV, "A=B*C", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_READ))
    
    call dfbydy(d_store7,d_store2)

!   D/DZ RHO W E
!   RHO W E IS PARALLEL
    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(math_kernel_eqV, "A=B*C", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_READ))

    call dfbydz(d_store7,d_store3)

!   COLLECT DIV RHO U E IN STORE4 FOR NOW
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(math_kernel_eqM, "A=A+B+C+D", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ))

!                                          STORE4 = E DIV RHO U + DIV RHO U E
!   =========================================================================

!   E EQUATION: CONVECTIVE TERMS
!   ----------------------------
!   HALF RHO U.DEL E

    call dfbydx(d_erhs,d_store1)
    call dfbydy(d_erhs,d_store2)
    call dfbydz(d_erhs,d_store3)

!   COLLECT ALL CONVECTIVE TERMS IN ERHS
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(math_kernel_eqAD, "A = -half*(B+C*D+E*F+G*H)", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_READ))

!   -------------------------------------
!   E EQUATION: CONVECTIVE TERMS COMPLETE
!   -------------------------------------
!                                                            ALL STORES CLEAR
!   =========================================================================

!   E-EQUATION: HEAT FLUX TERMS
!   ---------------------------

!   TEMPERATURE GRADIENTS
    call dfbydx(d_trun,d_store1)
    call dfbydy(d_trun,d_store2)
    call dfbydz(d_trun,d_store3)

!                                                       STORE1,2,3 = DTDX,Y,Z
!   =========================================================================

!   COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs
!   ---------------------------------------------

!   X-DIRECTION
    IF(fxlcnv) THEN
        rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
        call ops_par_loop(boundary_kernel_temperature_xdir, "COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_bcltxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

    END IF
    IF(fxrcnv) THEN
        rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(boundary_kernel_temperature_xdir, "COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_bcltxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

    END IF

!   Y-DIRECTION
    IF(fylcnv)THEN
        rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
        call ops_par_loop(boundary_kernel_temperature_ydir, "COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_bcltyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

    END IF
    IF(fyrcnv)THEN
        rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
        call ops_par_loop(boundary_kernel_temperature_ydir, "COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_bcltyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

    END IF

!   Z-DIRECTION
    IF(fzlcnv)THEN
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
        call ops_par_loop(boundary_kernel_temperature_zdir, "COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_bcltzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

    END IF
    IF(fzrcnv)THEN
        rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
        call ops_par_loop(boundary_kernel_temperature_zdir, "COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_bcltzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

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
    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(math_kernel_eqAP, "THERMAL CONDUCTIVITY", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_transp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ))

!   MIXTURE AVERAGED TRANSPORT
!   RSC 17-APR-2013
!   THERMAL CONDUCTIVITY

    IF(flmavt)THEN
    DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
!             CONDUCTIVITY FOR EACH SPECIES
        transp(ic,jc,kc) = LOG(trun(ic,jc,kc)/tdifgb)
        DO ispec = 1, nspec
          fornow = condco(ncocon,ispec)
          DO icp = ncocm1,1,-1
            fornow = fornow*transp(ic,jc,kc) + condco(icp,ispec)
          END DO
          ctrans(ispec) = EXP(fornow)
        END DO
        
!             COMBINATION RULE FOR CONDUCTIVITY
        combo1 = zero
        combo2 = zero
        combo3 = zero
        DO ispec = 1, nspec
          fornow = yrhs(ispec,ic,jc,kc)*ovwmol(ispec)
          combo1 = combo1 + fornow*ctrans(ispec)
          combo2 = combo2 + fornow/ctrans(ispec)
          combo3 = combo3 + fornow
        END DO
!             RSC/GVN 08-MAR-2014 BUG FIX
!              COMBO3 = DRHS(IC,JC,KC)/COMBO3
!              COMBO1 = COMBO1*COMBO3
!              COMBO2 = COMBO2*COMBO3
!              STORE7(IC,JC,KC) = HALF*(COMBO1 + ONE/COMBO2)
!              WMOMIX(IC,JC,KC) = COMBO3
        combo3 = one/combo3
        combo1 = combo1*combo3
        combo2 = combo2*combo3
        store7(ic,jc,kc) = half*(combo1 + one/combo2)
        wmomix(ic,jc,kc) = drhs(ic,jc,kc)*combo3
        
      END DO
    END DO
    END DO
  
    END IF

!   CONDUCTIVITY GRADIENTS
    call dfbydx(d_store7,d_store4)
    call dfbydy(d_store7,d_store5)
    call dfbydz(d_store7,d_store6)

!   BOUNDARY CONDITIONS
!   BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxlcon) call zeroxl(d_store4)
    IF(fxrcon) call zeroxr(d_store4)

!   BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fylcon) call zeroyl(d_store5)
    IF(fyrcon) call zeroyr(d_store5)

!   BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
!   RSC/TDD BUG FIX FZLCON
    IF(fzlcon) call zerozl(d_store6)
    IF(fzrcon) call zerozr(d_store6)

    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(math_kernel_eqAA, "A = A+(B*C+D*E+F*G)", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store5, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store6, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ))

!                                                       STORE1,2,3 = DTDX,Y,Z
!                                                       STORE7 = CONDUCTIVITY
!   =========================================================================

!   E-EQUATION: HEAT FLUX TERMS
!   ---------------------------
!   WALL BC: THERMAL CONDUCTION TERMS
    IF(fxlcnw) THEN
        rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
        call ops_par_loop(heat_flux_kernel_thermal_fxlcnw, "HEAT FLUX: Thermal fxlcnw", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                        ops_arg_gbl(acbcxl, 1, "real(dp)", OPS_READ))

    END IF
    IF(fxrcnw) THEN
        rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(heat_flux_kernel_thermal_fxrcnw, "HEAT FLUX: Thermal fxrcnw", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                    ops_arg_gbl(acbcxr, 1, "real(dp)", OPS_READ))

    END IF
    IF(fylcnw) THEN
        rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
        call ops_par_loop(heat_flux_kernel_thermal_fylcnw, "HEAT FLUX: Thermal fylcnw", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store2, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                    ops_arg_gbl(acbcyl, 1, "real(dp)", OPS_READ))

    END IF
    IF(fyrcnw) THEN
        rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
        call ops_par_loop(heat_flux_kernel_thermal_fyrcnw, "HEAT FLUX: Thermal fyrcnw", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store2, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                    ops_arg_gbl(acbcyr, 1, "real(dp)", OPS_READ))

    END IF
    IF(fzlcnw) THEN
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
        call ops_par_loop(heat_flux_kernel_thermal_fzlcnw, "HEAT FLUX: Thermal fzlcnw", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store3, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                    ops_arg_gbl(acbczl, 1, "real(dp)", OPS_READ))    

    END IF
    IF(fzrcnw)THEN
        rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
        call ops_par_loop(heat_flux_kernel_thermal_fzrcnw, "HEAT FLUX: Thermal fzrcnw", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store3, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                    ops_arg_gbl(acbczr, 1, "real(dp)", OPS_READ))
    END IF

!   =========================================================================

!   E-EQUATION: HEAT FLUX TERMS
!   ---------------------------
!   SECOND DERIVATIVE TERMS

!   TEMPERATURE SECOND DERIVATIVES
    call d2fdx2(d_trun,d_store1)
    call d2fdy2(d_trun,d_store2)
    call d2fdz2(d_trun,d_store3)

!   BOUNDARY CONDITIONS
!   BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxlcon) call zeroxl(d_store1)
    IF(fxrcon) call zeroxr(d_store1)

!   BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fylcon) call zeroyl(d_store2)
    IF(fyrcon) call zeroyr(d_store2)

!   BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
!   RSC 28-JUN-2015 BUG FIX FZLCON
    IF(fzlcon) call zerozl(d_store3)
    IF(fzrcon) call zerozr(d_store3)

!   COLLECT CONDUCTIVITY TERMS
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(math_kernel_eqY, "A = A+(B+C+D)*E", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ))

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
    IF(flradn) call radcal

!   =========================================================================

!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   SPECIES MASS FRACTION EQUATIONS
!   ===============================

!   REACTION RATE FOR ALL SPECIES
!   -----------------------------
    call chrate
!---UA
    DO ispec = 1,nspec
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_MD_kernel_eqB, "A_multidim = B_multidim", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_rrte, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

    END DO
!---end-UA
!                                                        RATE = REACTION RATE
!   =========================================================================

!   COLLECT REACTION RATE FOR BCs
!   -----------------------------

!   X-DIRECTION
    IF(fxlcnv)THEN
        DO ispec = 1,nspec
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_reaction_xdir, "COLLECT REACTION RATE FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ratexl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END DO
    END IF
    IF(fxrcnv)THEN
        DO ispec = 1,nspec
            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_reaction_xdir, "COLLECT REACTION RATE FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ratexr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END DO
    END IF

!   Y-DIRECTION
    IF(fylcnv)THEN
        DO ispec = 1,nspec
            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(boundary_kernel_reaction_ydir, "COLLECT REACTION RATE FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_rateyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END DO
    END IF
    IF(fyrcnv)THEN
        DO ispec = 1,nspec
            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_reaction_ydir, "COLLECT REACTION RATE FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_rateyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END DO
    END IF

!   Z-DIRECTION
    IF(fzlcnv)THEN
        DO ispec = 1,nspec
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(boundary_kernel_reaction_zdir, "COLLECT REACTION RATE FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ratezl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END DO
    END IF
    IF(fzrcnv)THEN
        DO ispec = 1,nspec
            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(boundary_kernel_reaction_zdir, "COLLECT REACTION RATE FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ratezr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END DO
    END IF

!                                                        RATE = REACTION RATE
!   =========================================================================

!   ZERO THE ACCUMULATORS FOR THE DIFFUSION CORRECTION VELOCITY
!   AND ITS DIVERGENCE

    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE))
    

!   ZERO THE ACCUMULATOR FOR THE MIXTURE ENTHALPY
!   MIXTURE H IS PARALLEL

    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(dp)", OPS_WRITE))

!                                                        RATE = REACTION RATE
!                                                         VTMP = DIV CORR VEL
!                                                     WTMP = MIXTURE H
!   =========================================================================

!   MIXTURE AVERAGED TRANSPORT
!   RSC 17-APR-2013
!   EVALUATE FIRST AND SECOND DERIVATIVES
!   OF LN(MIXTURE MOLAR MASS), LN(PRESSURE) AND LN(TEMPERATURE)

!   MIXTURE MOLAR MASS
    IF(flmixw) THEN

        rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
        call ops_par_loop(math_kernel_eqA, "A=log(B)", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_wmomix, 1, s3d_000, "real(dp)", OPS_READ))

        call dfbydx(d_store7,d_wd1x)
        call dfbydy(d_store7,d_wd1y)
        call dfbydz(d_store7,d_wd1z)
        call d2fdx2(d_store7,d_wd2x)
        call d2fdy2(d_store7,d_wd2y)
        call d2fdz2(d_store7,d_wd2z)

    END IF

!   PRESSURE
    IF(flmixp) THEN
    
        rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
        call ops_par_loop(math_kernel_eqA, "A=log(B)", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_prun, 1, s3d_000, "real(dp)", OPS_READ))

        call dfbydx(d_store7,d_pd1x)
        call dfbydy(d_store7,d_pd1y)
        call dfbydz(d_store7,d_pd1z)
        call d2fdx2(d_store7,d_pd2x)
        call d2fdy2(d_store7,d_pd2y)
        call d2fdz2(d_store7,d_pd2z)

    END IF

!   TEMPERATURE
    IF(flmixt)THEN
!       TRANSP CONTAINS LN(T/TDIFGB)
        call dfbydx(d_transp,d_td1x)
        call dfbydy(d_transp,d_td1y)
        call dfbydz(d_transp,d_td1z)
        call d2fdx2(d_transp,d_td2x)
        call d2fdy2(d_transp,d_td2y)
        call d2fdz2(d_transp,d_td2z)

    END IF

!   =========================================================================

!   RUN THROUGH ALL SPECIES
!   -----------------------
!   RSC 08-AUG-2012 EVALUATE ALL SPECIES
!   RSC 08-JUN-2015 REMOVE Nth SPECIES TREATMENT
    DO ispec = 1,nspec
  
!   =======================================================================
  
!       YRHS CONTAINS RHO Y: CONVERT TO Y
!       Y IS PARALLEL
        rangexyz = (/istalt,istolt,jstalt,jstolt,kstalt,kstolt/)
        call ops_par_loop(math_MD_kernel_eqE, "A_multidim = A_multidim/B", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!       =======================================================================
  
!       Y EQUATION: CONVECTIVE TERMS
!       ----------------------------
!       HALF Y DIV RHO U
  
!       COLLECT Y SOURCE TERMS IN RATE FOR NOW
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_MD_kernel_eqI, "A_multidim = A_multidim - half*B_multidim*C", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_divm, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!       =======================================================================
  
!       Y EQUATION: CONVECTIVE TERMS
!       ----------------------------
!       HALF DIV RHO U Y
  
!       D/DX RHO U Y
!       RHO U Y IS PARALLEL
        rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
        call ops_par_loop(math_MD_kernel_eqD, "A = B_multidim*C", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        call dfbydx(d_store7,d_store1)
  
!       D/DY RHO V Y
!       RHO V Y IS PARALLEL
        rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
        call ops_par_loop(math_MD_kernel_eqD, "A = B_multidim*C", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        call dfbydy(d_store7,d_store2)
  
!       D/DZ RHO W Y
!       RHO W Y IS PARALLEL
        rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
        call ops_par_loop(math_MD_kernel_eqD, "A = B_multidim*C", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        call dfbydz(d_store7,d_store3)
  
!       COLLECT DIV RHO U Y IN RATE FOR NOW
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_MD_kernel_eqH, "A_multidim = A_multidim - half*(B+C+D)", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!                                                  RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!       =======================================================================
  
!       SPECIES MASS FRACTION GRADIENT TERMS
!       ------------------------------------
  
!       SPECIES MASS FRACTION GRADIENTS
        rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
        call ops_par_loop(math_MD_kernel_eqA, "A = B_multidim", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        call dfbydx(d_store7,d_store1)
        call dfbydy(d_store7,d_store2)
        call dfbydz(d_store7,d_store3)

!                                                         STORE1,2,3 = DYDX,Y,Z
!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!       =======================================================================
  
!       COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs
!       -------------------------------------------------------
  
!       X-DIRECTION: DYDX
        IF(fxlcnv) THEN
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_mass_xdir, "COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_stryxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bclyxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fxrcnv) THEN
            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_mass_xdir, "COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_stryxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bclyxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF

!       Y-DIRECTION: DYDY
        IF(fylcnv) THEN
            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(boundary_kernel_mass_ydir, "COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_stryyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bclyyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fyrcnv) THEN
            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_mass_ydir, "COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_stryyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bclyyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF

!       Z-DIRECTION: DYDZ
        IF(fzlcnv) THEN
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(boundary_kernel_mass_zdir, "COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_stryzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bclyzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fzrcnv) THEN
            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(boundary_kernel_mass_zdir, "COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_stryzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bclyzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

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
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_MD_kernel_eqJ, "A_multidim = A_multidim - half*(B*C+D*E+F*G)", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

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
        rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
        call ops_par_loop(math_kernel_eqD, "A=B*val", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_transp, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(olewis(ispec), 1, "real(dp)", OPS_READ))

!                                                         STORE1,2,3 = DYDX,Y,Z
!                                                          STORE7 = DIFFUSIVITY
!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!       -----------------------------------------------------------------------
  
!       MIXTURE AVERAGED TRANSPORT
!       RSC 17-APR-2013
!       TRANSP CONTAINS LN(T/TDIFGB)
        IF(flmavt) THEN
    
!       MASS DIFFUSIVITY FOR EACH SPECIES
!       RELATIVE TO CURRENT SPECIES
        DO kc = kstab,kstob
        DO jc = jstab,jstob
        DO ic = istab,istob
          
          DO jspec = 1, nspec
            fornow = diffco(ncodif,jspec,ispec)
            DO icp = ncodm1,1,-1
              fornow = fornow*transp(ic,jc,kc) + diffco(icp,jspec,ispec)
            END DO
            ctrans(jspec) = EXP(fornow)*pdifgb/prun(ic,jc,kc)
          END DO
          
!               COMBINATION RULE FOR MASS DIFFUSIVITY
          combo1 = zero
          combo2 = zero
          DO jspec = 1, nspec
            fornow = yrhs(jspec,ic,jc,kc) + dfctol
            combo1 = combo1 + fornow
            combo2 = combo2 + fornow*ovwmol(jspec)/ctrans(jspec)
          END DO
          fornow = yrhs(ispec,ic,jc,kc) + dfctol
          combo1 = combo1 - fornow
          combo2 = combo2 - fornow*ovwmol(ispec)/ctrans(ispec)
          combo2 = combo2*wmomix(ic,jc,kc)
          difmix(ic,jc,kc) = drhs(ic,jc,kc)*combo1/combo2
          store7(ic,jc,kc) = difmix(ic,jc,kc)
          
        END DO
        END DO
            END DO
    
        END IF
  
!       -----------------------------------------------------------------------
  
!       MIXTURE AVERAGED TRANSPORT
!       RSC 17-APR-2013
!       TRANSP CONTAINS LN(T/TDIFGB)
        IF(flmtdr(ispec))THEN
    
!           THERMAL DIFFUSION RATIO FOR EACH SPECIES
!           RELATIVE TO CURRENT SPECIES
  
            rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
            call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_WRITE))
    
            DO jspec = 1, nspec
      
                flmtds = flmtdr(jspec).AND.(ispec /= jspec)
                IF(flmtds)THEN
        
                DO kc = kstab,kstob
                DO jc = jstab,jstob
                DO ic = istab,istob
              
!               THERMAL DIFFUSION RATIO FOR THIS SPECIES PAIR
                combo2 = trun(ic,jc,kc)/tdifgb
                fornow = tdrcco(ncotdr,jspec,ispec)
                DO icp = ncotm1,1,-1
                fornow = fornow*combo2 + tdrcco(icp,jspec,ispec)
                END DO
                ctrans(jspec) = fornow
              
!               COMBINATION RULE FOR THERMAL DIFFUSIION RATIO
                fornow = yrhs(jspec,ic,jc,kc)*ovwmol(jspec)
                tdrmix(ic,jc,kc) = tdrmix(ic,jc,kc)  &
                  + fornow*wmomix(ic,jc,kc)*ctrans(jspec)
              
                END DO
                END DO
                END DO
        
                END IF
      
            END DO
    
        END IF
  
!       =======================================================================
  
!       DIFFUSION CORRECTION VELOCITY
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &    
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ))

        call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ))

        call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ))

!                                                         STORE1,2,3 = DYDX,Y,Z
!                                                          STORE7 = DIFFUSIVITY
!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!       =======================================================================
  
!       SPECIES ENTHALPY
!       ----------------
  
!       TEMPERATURE INTERVAL INDEXING
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1
        icoef2 = ntbase**ipower
        icoef1 = icoef2*ntbase
  
!       SPECIES H IS PARALLEL
!       STORE SPECIES H IN UTMP FOR NOW
!       STORE MIXTURE H IN WTMP FOR NOW
        DO kc = kstab,kstob
        DO jc = jstab,jstob
        DO ic = istab,istob
        
        itint = 1 + MOD(itndex(ic,jc,kc,iindex),icoef1)/icoef2
        fornow = amasch(ncpoly(itint,ispec),itint,ispec)
        DO icp = ncpom1(itint,ispec),1,-1
          fornow = fornow*trun(ic,jc,kc) + amasch(icp,itint,ispec)
        END DO
        utmp(ic,jc,kc) = amasch(ncenth(itint,ispec),itint,ispec)  &
            + fornow*trun(ic,jc,kc)
        
!             MIXTURE H
        wtmp(ic,jc,kc) = wtmp(ic,jc,kc) + utmp(ic,jc,kc)*yrhs(ispec,ic,jc,kc)
        
        END DO
        END DO
        END DO
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
        IF(fxlcnv) THEN
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_speciesH_xdir, "COLLECT SPECIES H FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strhxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fxrcnv) THEN
            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_speciesH_xdir, "COLLECT SPECIES H FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strhxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
  
!       Y-DIRECTION
        IF(fylcnv) THEN
            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(boundary_kernel_speciesH_ydir, "COLLECT SPECIES H FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strhyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fyrcnv) THEN
            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_speciesH_ydir, "COLLECT SPECIES H FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strhyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
  
!       Z-DIRECTION
        IF(fzlcnv) THEN
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(boundary_kernel_speciesH_zdir, "COLLECT SPECIES H FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strhzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fzrcnv) THEN
            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(boundary_kernel_speciesH_zdir, "COLLECT SPECIES H FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strhzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

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
        IF(flmduf(ispec)) THEN
            rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
            call ops_par_loop(math_kernel_eqF, "A = A+var*B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
!       Y EQUATION: DIFFUSIVE TERMS
!       ---------------------------
!       E EQUATION: FURTHER HEAT FLUX TERMS
!       DIFFUSIVITY GRADIENT TERMS
  
!       DIFFUSIVITY GRADIENTS
        call dfbydx(d_store7,d_store4)
        call dfbydy(d_store7,d_store5)
        call dfbydz(d_store7,d_store6)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fxldif) call zeroxl(d_store4)
        IF(fxrdif) call zeroxr(d_store4)
!       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fyldif) call zeroyl(d_store5)
        IF(fyrdif) call zeroyr(d_store5)
!       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fzldif) call zerozl(d_store6)
        IF(fzrdif) call zerozr(d_store6)
 
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_MD_kernel_eqL, "multiple math equations", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(fxladb) call zeroxl(d_store4)
        IF(fxradb) call zeroxr(d_store4)

!       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(fyladb) call zeroyl(d_store5)
        IF(fyradb) call zeroyr(d_store5)

!       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(fzladb) call zerozl(d_store6)
        IF(fzradb) call zerozr(d_store6)

!       E EQUATION
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_kernel_eqAB, "A = A+(B*C+D*E+F*G)*H", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ))

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
        call dfbydx(d_utmp,d_store4)
        call dfbydy(d_utmp,d_store5)
        call dfbydz(d_utmp,d_store6)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fxldif) call zeroxl(d_store4)
        IF(fxrdif) call zeroxr(d_store4)

!       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fyldif) call zeroyl(d_store5)
        IF(fyrdif) call zeroyr(d_store5)

!       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fzldif) call zerozl(d_store6)
        IF(fzrdif) call zerozr(d_store6)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(fxladb) call zeroxl(d_store4)
        IF(fxradb) call zeroxr(d_store4)

!       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(fyladb) call zeroyl(d_store5)
        IF(fyradb) call zeroyr(d_store5)

!       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(fzladb) call zerozl(d_store6)
        IF(fzradb) call zerozr(d_store6)

        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_kernel_eqAB, "A = A+(B*C+D*E+F*G)*H", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store5, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store6, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ))

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
        IF(fxldfw) THEN
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fxldfw, "HEAT FLUX: Enthalpy fxldfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store1, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcxl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fxrdfw) THEN
            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fxrdfw, "HEAT FLUX: Enthalpy fxrdfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store1, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcxr, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fyldfw) THEN
            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fyldfw, "HEAT FLUX: Enthalpy fyldfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store2, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcyl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fyrdfw) THEN
            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fyrdfw, "HEAT FLUX: Enthalpy fyrdfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store2, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcyr, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fzldfw) THEN
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fzldfw, "HEAT FLUX: Enthalpy fzldfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbczl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
        IF(fzrdfw) THEN
            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fzrdfw, "HEAT FLUX: Enthalpy fzrdfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbczr, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END IF
  
!       =======================================================================
  
!       Y-EQUATION: DIFFUSIVE TERMS
!       ---------------------------
!       E-EQUATION: FURTHER HEAT FLUX TERMS
!       SECOND DERIVATIVE TERMS
  
!       SPECIES MASS FRACTION SECOND DERIVATIVES
!       MOVE DIFFUSIVITY TO STORE4
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ))

!       MOVE MASS FRACTION TO STORE7
        rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
        call ops_par_loop(math_MD_kernel_eqA, "A = B_multidim", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        call d2fdx2(d_store7,d_store1)
        call d2fdy2(d_store7,d_store2)
        call d2fdz2(d_store7,d_store3)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fxldif) call zeroxl(d_store1)
        IF(fxrdif) call zeroxr(d_store1)

!       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fyldif) call zeroyl(d_store2)
        IF(fyrdif) call zeroyr(d_store2)

!       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fzldif) call zerozl(d_store3)
        IF(fzrdif) call zerozr(d_store3)
 
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_MD_kernel_eqM, "multiple math equations", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(fxladb) call zeroxl(d_store1)
        IF(fxradb) call zeroxr(d_store1)

!       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(fyladb) call zeroyl(d_store2)
        IF(fyradb) call zeroyr(d_store2)

!       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
        IF(fzladb) call zerozl(d_store3)
        IF(fzradb) call zerozr(d_store3)
  
!       E EQUATION
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_kernel_eqZ, "A=A+(B+C+D)*E*F", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ))

!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                              WTMP = MIXTURE H
!       =======================================================================
  
!       MIXTURE AVERAGED TRANSPORT
!       RSC 23-APR-2013
!       MOLAR MASS TERMS, PRESSURE TERMS, SORET EFFECT
  
!       MIXTURE MOLAR MASS TERMS
        IF(flmixw) THEN
!           FIRST AND SECOND DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED
            rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
            call ops_par_loop(math_MD_kernel_eqC, "A = B*C_multidim", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_difmix, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!           DIFFUSION CORRECTION VELOCITY
!           FIRST DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1x, 1, s3d_000, "real(dp)", OPS_READ))

            call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1y, 1, s3d_000, "real(dp)", OPS_READ))

            call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1z, 1, s3d_000, "real(dp)", OPS_READ)) 
    
!           Y EQUATION: DIFFUSIVE TERMS
!           E EQUATION: FURTHER HEAT FLUX TERMS
    
!           DIFFUSIVITY GRADIENT TERMS
    
!           DIFFUSIVITY GRADIENTS
            call dfbydx(d_store7,d_store1)
            call dfbydy(d_store7,d_store2)
            call dfbydz(d_store7,d_store3)
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fxldif) call zeroxl(d_store1)
            IF(fxrdif) call zeroxr(d_store1)

!           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fyldif) call zeroyl(d_store2)
            IF(fyrdif) call zeroyr(d_store2)

!           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fzldif) call zerozl(d_store3)
            IF(fzrdif) call zerozr(d_store3)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_MD_kernel_eqL, "multiple math equations", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wd1x, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wd1y, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wd1z, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))            
 
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fxladb) call zeroxl(d_store1)
            IF(fxradb) call zeroxr(d_store1)

!           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fyladb) call zeroyl(d_store2)
            IF(fyradb) call zeroyr(d_store2)

!           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fzladb) call zerozl(d_store3)
            IF(fzradb) call zerozr(d_store3)

!           E EQUATION
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqAB, "A = A+(B*C+D*E+F*G)*H", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1x, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1y, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1z, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ))

!           E-EQUATION: FURTHER HEAT FLUX TERMS
!           SPECIES ENTHALPY GRADIENT TERMS
    
!           SPECIES ENTHALPY GRADIENTS ALREADY IN STORE4,5,6
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fxldif) call zeroxl(d_store4)
            IF(fxrdif) call zeroxr(d_store4)

!           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fyldif) call zeroyl(d_store5)
            IF(fyrdif) call zeroyr(d_store5)

!           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fzldif) call zerozl(d_store6)
            IF(fzrdif) call zerozr(d_store6)
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fxladb) call zeroxl(d_store4)
            IF(fxradb) call zeroxr(d_store4)

!           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fyladb) call zeroyl(d_store5)
            IF(fyradb) call zeroyr(d_store5)

!           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fzladb) call zerozl(d_store6)
            IF(fzradb) call zerozr(d_store6)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqAB, "A = A+(B*C+D*E+F*G)*H", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1x, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store5, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1y, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store6, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1z, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ))

!           ====================================================================
    
!           Y-EQUATION: DIFFUSIVE TERMS
!           ---------------------------
!           WALL BC: MOLAR MASS TERMS
!           E-EQUATION: HEAT FLUX TERMS
!           WALL BC: ENTHALPY DIFFUSION TERMS
            IF(fxldfw) THEN
                rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fxldfw, "HEAT FLUX: Enthalpy fxldfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1x, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcxl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fxrdfw) THEN
                rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fxrdfw, "HEAT FLUX: Enthalpy fxrdfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1x, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcxr, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fyldfw) THEN
                rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fyldfw, "HEAT FLUX: Enthalpy fyldfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1y, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcyl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fyrdfw) THEN
                rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fyrdfw, "HEAT FLUX: Enthalpy fyrdfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1y, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcyr, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fzldfw) THEN
                rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fzldfw, "HEAT FLUX: Enthalpy fzldfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1z, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbczl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fzrdfw) THEN
                rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fzrdfw, "HEAT FLUX: Enthalpy fzrdfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd1z, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbczr, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
    
!           ====================================================================
    
!           Y-EQUATION: DIFFUSIVE TERMS
!           E-EQUATION: FURTHER HEAT FLUX TERMS
!           SECOND DERIVATIVE TERMS
!           SECOND DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fxldif) call zeroxl(d_wd2x)
            IF(fxrdif) call zeroxr(d_wd2x)

!           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fyldif) call zeroyl(d_wd2y)
            IF(fyrdif) call zeroyr(d_wd2y)

!           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fzldif) call zerozl(d_wd2z)
            IF(fzrdif) call zerozr(d_wd2z)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_MD_kernel_eqM, "multiple math equations", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_wd2x, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wd2y, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wd2z, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))            

!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fxladb) call zeroxl(d_wd2x)
            IF(fxradb) call zeroxr(d_wd2x)

!           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fyladb) call zeroyl(d_wd2y)
            IF(fyradb) call zeroyr(d_wd2y)

!           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fzladb) call zerozl(d_wd2z)
            IF(fzradb) call zerozr(d_wd2z)
    
!           E EQUATION
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqZ, "A=A+(B+C+D)*E*F", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_wd2x, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd2y, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wd2z, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ))

        END IF  !flmixw
!       MIXTURE MOLAR MASS TERMS
  
!       =======================================================================
  
!       PRESSURE DIFFUSION TERMS
        IF(flmixp) THEN
!           FIRST AND SECOND DERIVATIVES OF LN(PRESSURE) ALREADY STORED

            rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
            call ops_par_loop(math_MD_kernel_eqG, "A = B*C_multidim*(one-const_val/D)", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_difmix, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_wmomix, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_gbl(wmolar(ispec), 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!           DIFFUSION CORRECTION VELOCITY
!           FIRST DERIVATIVES OF LN(PRESSURE) ALREADY STORED
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1x, 1, s3d_000, "real(dp)", OPS_READ))

            call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1y, 1, s3d_000, "real(dp)", OPS_READ))

            call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1z, 1, s3d_000, "real(dp)", OPS_READ))

!           Y EQUATION: DIFFUSIVE TERMS
!           E EQUATION: FURTHER HEAT FLUX TERMS
    
!           DIFFUSIVITY GRADIENT TERMS
    
!           DIFFUSIVITY GRADIENTS
            call dfbydx(d_store7,d_store1)
            call dfbydy(d_store7,d_store2)
            call dfbydz(d_store7,d_store3)
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fxldif) call zeroxl(d_store1)
            IF(fxrdif) call zeroxr(d_store1)

!           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fyldif) call zeroyl(d_store2)
            IF(fyrdif) call zeroyr(d_store2)

!           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fzldif) call zerozl(d_store3)
            IF(fzrdif) call zerozr(d_store3)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_MD_kernel_eqL, "multiple math equations", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_pd1x, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_pd1y, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_pd1z, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fxladb) call zeroxl(d_store1)
            IF(fxradb) call zeroxr(d_store1)

!           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fyladb) call zeroyl(d_store2)
            IF(fyradb) call zeroyr(d_store2)

!           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fzladb) call zerozl(d_store3)
            IF(fzradb) call zerozr(d_store3)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqAB, "A = A+(B*C+D*E+F*G)*H", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1x, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1y, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1z, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ))

!           E-EQUATION: FURTHER HEAT FLUX TERMS
!           SPECIES ENTHALPY GRADIENT TERMS
    
!           SPECIES ENTHALPY GRADIENTS ALREADY IN STORE4,5,6
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fxldif) call zeroxl(d_store4)
            IF(fxrdif) call zeroxr(d_store4)

!           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fyldif) call zeroyl(d_store5)
            IF(fyrdif) call zeroyr(d_store5)

!           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fzldif) call zerozl(d_store6)
            IF(fzrdif) call zerozr(d_store6)
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fxladb) call zeroxl(d_store4)
            IF(fxradb) call zeroxr(d_store4)

!           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fyladb) call zeroyl(d_store5)
            IF(fyradb) call zeroyr(d_store5)

!           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fzladb) call zerozl(d_store6)
            IF(fzradb) call zerozr(d_store6)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqAB, "A = A+(B*C+D*E+F*G)*H", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1x, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store5, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1y, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store6, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1z, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ))

!           ====================================================================
    
!           Y-EQUATION: DIFFUSIVE TERMS
!           ---------------------------
!           WALL BC: PRESSURE TERMS
!           E-EQUATION: HEAT FLUX TERMS
!           WALL BC: ENTHALPY DIFFUSION TERMS
            IF(fxldfw) THEN
                rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fxldfw, "HEAT FLUX: Enthalpy fxldfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1x, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcxl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fxrdfw) THEN
                rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fxrdfw, "HEAT FLUX: Enthalpy fxrdfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1x, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcxr, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fyldfw) THEN
                rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fyldfw, "HEAT FLUX: Enthalpy fyldfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1y, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcyl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fyrdfw) THEN
                rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fyrdfw, "HEAT FLUX: Enthalpy fyrdfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1y, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbcyr, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fzldfw) THEN
                rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fzldfw, "HEAT FLUX: Enthalpy fzldfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1z, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbczl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fzrdfw) THEN
                rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(heat_flux_kernel_enthalpy2_fzrdfw, "HEAT FLUX: Enthalpy fzrdfw", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd1z, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                            ops_arg_gbl(acbczr, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
    
!           ====================================================================
    
!           Y-EQUATION: DIFFUSIVE TERMS
!           E-EQUATION: FURTHER HEAT FLUX TERMS
!           SECOND DERIVATIVE TERMS
!           SECOND DERIVATIVES OF LN(PRESSURE) ALREADY STORED
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fxldif) call zeroxl(d_pd2x)
            IF(fxrdif) call zeroxr(d_pd2x)

!           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fyldif) call zeroyl(d_pd2y)
            IF(fyrdif) call zeroyr(d_pd2y)

!           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fzldif) call zerozl(d_pd2z)
            IF(fzrdif) call zerozr(d_pd2z)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_MD_kernel_eqM, "multiple math equations", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_pd2x, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_pd2y, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_pd2z, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))
 
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fxladb) call zeroxl(d_pd2x)
            IF(fxradb) call zeroxr(d_pd2x)

!           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fyladb) call zeroyl(d_pd2y)
            IF(fyradb) call zeroyr(d_pd2y)

!           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fzladb) call zerozl(d_pd2z)
            IF(fzradb) call zerozr(d_pd2z)

!           E EQUATION
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqZ, "A=A+(B+C+D)*E*F", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_pd2x, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd2y, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_pd2z, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
!       PRESSURE DIFFUSION TERMS

!       =======================================================================
  
!       SORET EFFECT (THERMAL DIFFUSION) TERMS
        IF(flmsor(ispec))THEN
!           FIRST AND SECOND DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED

            rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
            call ops_par_loop(math_MD_kernel_eqF, "A = B*C_multidim*D", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_difmix, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!           DIFFUSION CORRECTION VELOCITY
!           FIRST DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td1x, 1, s3d_000, "real(dp)", OPS_READ))

            call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td1y, 1, s3d_000, "real(dp)", OPS_READ))

            call ops_par_loop(math_kernel_eqN, "A=A+B*C", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td1z, 1, s3d_000, "real(dp)", OPS_READ))

!           Y EQUATION: DIFFUSIVE TERMS
!           E EQUATION: FURTHER HEAT FLUX TERMS
    
!           DIFFUSIVITY GRADIENT TERMS
    
!           DIFFUSIVITY GRADIENTS
            call dfbydx(d_store7,d_store1)
            call dfbydy(d_store7,d_store2)
            call dfbydz(d_store7,d_store3)
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fxldif) call zeroxl(d_store1)
            IF(fxrdif) call zeroxr(d_store1)

!           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fyldif) call zeroyl(d_store2)
            IF(fyrdif) call zeroyr(d_store2)

!           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fzldif) call zerozl(d_store3)
            IF(fzrdif) call zerozr(d_store3)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_MD_kernel_eqL, "multiple math equations", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_td1x, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_td1y, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_td1z, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

!           SUBTRACT DUFOUR EFFECT TERMS TO RESTORE SPECIES ENTHALPY
!           RSC 08-JUN-2015 BUG FIX
            IF(flmduf(ispec))THEN
                rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
                call ops_par_loop(math_kernel_eqH, "A = A-var*B*C", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

            END IF
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fxladb) call zeroxl(d_store1)
            IF(fxradb) call zeroxr(d_store1)

!           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fyladb) call zeroyl(d_store2)
            IF(fyradb) call zeroyr(d_store2)

!           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fzladb) call zerozl(d_store3)
            IF(fzradb) call zerozr(d_store3)

!           E EQUATION
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqAB, "A = A+(B*C+D*E+F*G)*H", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td1x, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td1y, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td1z, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ))

!           E-EQUATION: FURTHER HEAT FLUX TERMS
!           SPECIES ENTHALPY GRADIENT TERMS
    
!           EVALUATE SPECIES ENTHALPY GRADIENTS USING STORE4,5,6
            call dfbydx(d_utmp,d_store4)
            call dfbydy(d_utmp,d_store5)
            call dfbydz(d_utmp,d_store6)
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fxldif) call zeroxl(d_store4)
            IF(fxrdif) call zeroxr(d_store4)

!           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fyldif) call zeroyl(d_store5)
            IF(fyrdif) call zeroyr(d_store5)

!           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fzldif) call zerozl(d_store6)
            IF(fzrdif) call zerozr(d_store6)
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fxladb) call zeroxl(d_store4)
            IF(fxradb) call zeroxr(d_store4)

!           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fyladb) call zeroyl(d_store5)
            IF(fyradb) call zeroyr(d_store5)

!           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fzladb) call zerozl(d_store6)
            IF(fzradb) call zerozr(d_store6)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqAB, "A = A+(B*C+D*E+F*G)*H", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td1x, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store5, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td1y, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store6, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td1z, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ))

!           ====================================================================
    
!           Y-EQUATION: DIFFUSIVE TERMS
!           ---------------------------
!           WALL BC: SORET EFFECT TERMS
!           E-EQUATION: HEAT FLUX TERMS
!           WALL BC: ENTHALPY DIFFUSION TERMS
            IF(fxldfw) THEN
                rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
                call ops_par_loop(heat_flux_kernel_enthalpy_fxldfw, "HEAT FLUX: Enthalpy fxldfw", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_store7, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_td1x, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                                ops_arg_gbl(acbcxl, 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fxrdfw) THEN
                rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
                call ops_par_loop(heat_flux_kernel_enthalpy_fxrdfw, "HEAT FLUX: Enthalpy fxrdfw", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_store7, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_td1x, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                                ops_arg_gbl(acbcxr, 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fyldfw) THEN
                rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
                call ops_par_loop(heat_flux_kernel_enthalpy_fyldfw, "HEAT FLUX: Enthalpy fyldfw", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_store7, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_td1y, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                                ops_arg_gbl(acbcyl, 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fyrdfw) THEN
                rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
                call ops_par_loop(heat_flux_kernel_enthalpy_fyrdfw, "HEAT FLUX: Enthalpy fyrdfw", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_store7, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_td1y, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                                ops_arg_gbl(acbcyr, 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fzldfw) THEN
                rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
                call ops_par_loop(heat_flux_kernel_enthalpy_fzldfw, "HEAT FLUX: Enthalpy fzldfw", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_store7, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_td1z, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                                ops_arg_gbl(acbczl, 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
            IF(fzrdfw) THEN
                rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
                call ops_par_loop(heat_flux_kernel_enthalpy_fzrdfw, "HEAT FLUX: Enthalpy fzrdfw", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_store7, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_td1z, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                                ops_arg_gbl(acbczr, 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END IF
    
!           E-EQUATION: HEAT FLUX TERMS
!           WALL BC: SORET AND DUFOUR TERMS
            IF(flmduf(ispec)) THEN
!               E-EQUATION: HEAT FLUX TERMS
!               WALL BC: ADIABATIC WALL
                IF(fxlcnw) THEN
                    rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
                    call ops_par_loop(heat_flux_kernel_soret_fxlcnw, "HEAT FLUX: Soret and DUFOUR fxlcnw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1x,   1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbcxl,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fxrcnw) THEN
                    rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
                    call ops_par_loop(heat_flux_kernel_soret_fxrcnw, "HEAT FLUX: Soret and DUFOUR fxrcnw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1x,   1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbcxr,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fylcnw) THEN
                    rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
                    call ops_par_loop(heat_flux_kernel_soret_fylcnw, "HEAT FLUX: Soret and DUFOUR fylcnw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1y,   1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbcyl,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fyrcnw) THEN
                    rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
                    call ops_par_loop(heat_flux_kernel_soret_fyrcnw, "HEAT FLUX: Soret and DUFOUR fyrcnw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1y,   1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbcyr,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fzlcnw) THEN
                    rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
                    call ops_par_loop(heat_flux_kernel_soret_fzlcnw, "HEAT FLUX: Soret and DUFOUR fzlcnw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1z,   1, s3d_p001_to_p004_z, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbczl,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fzrcnw) THEN
                    rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
                    call ops_par_loop(heat_flux_kernel_soret_fzrcnw, "HEAT FLUX: Soret and DUFOUR fzrcnw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1z,   1, s3d_m001_to_m004_z, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbczr,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
      
!               E-EQUATION: HEAT FLUX TERMS
!               WALL BC: ISOTHERMAL WALL
                IF(fxladw) THEN
                    rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
                    call ops_par_loop(heat_flux_kernel_isothermal_fxladw, "HEAT FLUX: Isothermal fxladw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_000_to_p400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_000_to_p400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_000_to_p400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1x,   1, s3d_000_to_p400_x, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbcxl,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fxradw) THEN
                    rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
                    call ops_par_loop(heat_flux_kernel_isothermal_fxradw, "HEAT FLUX: Isothermal fxradw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_000_to_m400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_000_to_m400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_000_to_m400_x, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1x,   1, s3d_000_to_m400_x, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbcxr,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fyladw) THEN
                    rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
                    call ops_par_loop(heat_flux_kernel_isothermal_fyladw, "HEAT FLUX: Isothermal fyladw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_000_to_p040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_000_to_p040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_000_to_p040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1y,   1, s3d_000_to_p040_y, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbcyl,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fyradw) THEN
                    rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
                    call ops_par_loop(heat_flux_kernel_isothermal_fyradw, "HEAT FLUX: Isothermal fyradw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_000_to_m040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_000_to_m040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_000_to_m040_y, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1y,   1, s3d_000_to_m040_y, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbcyr,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fzladw) THEN
                    rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
                    call ops_par_loop(heat_flux_kernel_isothermal_fzladw, "HEAT FLUX: Isothermal fzladw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_000_to_p004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_000_to_p004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_000_to_p004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1z,   1, s3d_000_to_p004_z, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbczl,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF
                IF(fzradw)THEN
                    rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
                    call ops_par_loop(heat_flux_kernel_isothermal_fzradw, "HEAT FLUX: Isothermal fzradw", senga_grid, 3, rangexyz,  &
                                    ops_arg_dat(d_erhs,   1, s3d_000, "real(dp)", OPS_WRITE), &
                                    ops_arg_dat(d_trun,   1, s3d_000_to_m004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_tdrmix, 1, s3d_000_to_m004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_store7, 1, s3d_000_to_m004_z, "real(dp)", OPS_READ), &
                                    ops_arg_dat(d_td1z,   1, s3d_000_to_m004_z, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(acbczr,   1, "real(dp)", OPS_READ), &
                                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ))

                END IF

            END IF
    
!           ====================================================================
    
!           Y-EQUATION: DIFFUSIVE TERMS
!           E-EQUATION: FURTHER HEAT FLUX TERMS
!           SECOND DERIVATIVE TERMS
!           SECOND DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED
    
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fxldif) call zeroxl(d_td2x)
            IF(fxrdif) call zeroxr(d_td2x)

!           BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fyldif) call zeroyl(d_td2y)
            IF(fyrdif) call zeroyr(d_td2y)

!           BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
            IF(fzldif) call zerozl(d_td2z)
            IF(fzrdif) call zerozr(d_td2z)

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_MD_kernel_eqM, "multiple math equations", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_td2x, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_td2y, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_td2z, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))
 
!           BOUNDARY CONDITIONS
!           BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fxladb) call zeroxl(d_td2x)
            IF(fxradb) call zeroxr(d_td2x)

!           BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fyladb) call zeroyl(d_td2y)
            IF(fyradb) call zeroyr(d_td2y)

!           BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
            IF(fzladb) call zerozl(d_td2z)
            IF(fzradb) call zerozr(d_td2z)

!           E EQUATION
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(math_kernel_eqZ, "A=A+(B+C+D)*E*F", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_td2x, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td2y, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_td2z, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_utmp, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!                                                      RATE = Y SOURCE TERMS
!                                                       VTMP = DIV CORR VEL
!                                                          WTMP = MIXTURE H
!   =======================================================================
  
!   ----------------------------------------------------------------
!   E-EQUATION: DIFFUSION CORRECTION VELOCITY TERMS EVALUATED BELOW
!   Y-EQUATION: DIFFUSION CORRECTION VELOCITY TERMS EVALUATED BELOW
!   ----------------------------------------------------------------
  
    END DO ! END of ispec loop
!   RSC 08-AUG-2012 EVALUATE ALL SPECIES
!   END OF RUN THROUGH ALL SPECIES

!   =========================================================================

!   EVALUATE DIFFUSION CORRECTION VELOCITY TERMS
!   --------------------------------------------

!   E-EQUATION: FURTHER HEAT FLUX TERMS
!   -----------------------------------
!   MIXTURE ENTHALPY GRADIENTS
    call dfbydx(d_wtmp,d_store1)
    call dfbydy(d_wtmp,d_store2)
    call dfbydz(d_wtmp,d_store3)

!   BOUNDARY CONDITIONS
!   BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif) call zeroxl(d_store1)
    IF(fxrdif) call zeroxr(d_store1)

!   BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif) call zeroyl(d_store2)
    IF(fyrdif) call zeroyr(d_store2)

!   BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif) call zerozl(d_store3)
    IF(fzrdif) call zerozr(d_store3)

!   BOUNDARY CONDITIONS
!   BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb) call zeroxl(d_store1)
    IF(fxradb) call zeroxr(d_store1)

!   BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb) call zeroyl(d_store2)
    IF(fyradb) call zeroyr(d_store2)

!   BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb) call zerozl(d_store3)
    IF(fzradb) call zerozr(d_store3)

!   TRANSFER DIV CORR VEL TO TEMPORARY STORE
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_READ))

!   BOUNDARY CONDITIONS
!   BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb) call zeroxl(d_store4)
    IF(fxradb) call zeroxr(d_store4)

!   BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb) call zeroyl(d_store4)
    IF(fyradb) call zeroyr(d_store4)

!   BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb) call zerozl(d_store4)
    IF(fzradb) call zerozr(d_store4)

!   DIV RHO VCORR HMIX
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(math_kernel_eqAH, "A = A-B*C-D*E-F*G-H*I", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_wtmp, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_READ))

!                                                       RATE = Y SOURCE TERMS
!                                                         VTMP = DIV CORR VEL
!   =========================================================================

!   MIXTURE AVERAGED TRANSPORT
!   EVALUATE THE VISCOSITY

!   RSC 17-APR-2013
!   TRANSP CONTAINS LN(T)
!   STORE VISCOSITY IN DIFMIX FOR NOW
    IF(flmavt)THEN
  
        DO kc = kstab,kstob
        DO jc = jstab,jstob
        DO ic = istab,istob
        
!             VISCOSITY FOR EACH SPECIES
        DO ispec = 1, nspec
          fornow = viscco(ncovis,ispec)
          DO icp = ncovm1,1,-1
            fornow = fornow*transp(ic,jc,kc) + viscco(icp,ispec)
          END DO
          ctrans(ispec) = EXP(fornow)
        END DO
        
!             COMBINATION RULE FOR VISCOSITY
        combo1 = zero
        DO ispec = 1, nspec
          combo2 = zero
          DO jspec = 1, nspec
            fornow = SQRT(ctrans(ispec)/ctrans(jspec))
            fornow = one + fornow*wilko2(jspec,ispec)
            fornow = wilko1(jspec,ispec)*fornow*fornow
            combo2 = combo2 + yrhs(jspec,ic,jc,kc)*ovwmol(jspec)*fornow
          END DO
          fornow = ctrans(ispec)/combo2
          combo1 = combo1 + yrhs(ispec,ic,jc,kc)*ovwmol(ispec)*fornow
          
        END DO
        difmix(ic,jc,kc) = combo1
        
        END DO
        END DO
        END DO
  
    END IF

!   =========================================================================

!   RUN THROUGH ALL SPECIES
!   -----------------------
!   RSC 08-AUG-2012 EVALUATE ALL SPECIES
!   RSC 08-JUN-2015 REMOVE Nth SPECIES TREATMENT
    DO ispec = 1,nspec
!       Y-EQUATION: DIFFUSIVE TERMS
!       ---------------------------
!       RECOMPUTE SPECIES MASS FRACTION GRADIENTS
        rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
        call ops_par_loop(math_MD_kernel_eqA, "A = B_multidim", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        call dfbydx(d_store7,d_store1)
        call dfbydy(d_store7,d_store2)
        call dfbydz(d_store7,d_store3)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fxldif) call zeroxl(d_store1)
        IF(fxrdif) call zeroxr(d_store1)

!       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fyldif) call zeroyl(d_store2)
        IF(fyrdif) call zeroyr(d_store2)

!       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
        IF(fzldif) call zerozl(d_store3)
        IF(fzrdif) call zerozr(d_store3)
  
!       DIV RHO VCORR Y
!       STORE Y SOURCE TERMS IN YRHS
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(math_MD_kernel_eqK, "A_multidim = B_multidim - A_multidim*C - D*E - F*G - H*I", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_rate, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

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
    IF(fxlcnv.OR.fxrcnv) THEN
  
        call dfbydx(d_drhs,d_store1)
  
        IF(fxlcnv) THEN
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_density_xdir, "COLLECT DENSITY AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

        END IF
        IF(fxrcnv) THEN
            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_density_xdir, "COLLECT DENSITY AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2xr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

        END IF

    END IF

!   Y-DIRECTION: DRHODY
    IF(fylcnv.OR.fyrcnv) THEN
  
        call dfbydy(d_drhs,d_store2)
  
        IF(fylcnv) THEN
            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(boundary_kernel_density_ydir, "COLLECT DENSITY AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

        END IF
        IF(fyrcnv) THEN
            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(boundary_kernel_density_ydir, "COLLECT DENSITY AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2yr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

        END IF
  
    END IF

!   Z-DIRECTION: DRHODZ
    IF(fzlcnv.OR.fzrcnv) THEN

        call dfbydz(d_drhs,d_store3)
  
        IF(fzlcnv) THEN
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(boundary_kernel_density_zdir, "COLLECT DENSITY AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

        END IF
        IF(fzrcnv) THEN
            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(boundary_kernel_density_zdir, "COLLECT DENSITY AND ITS GRADIENTS FOR BCs", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2zr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

        END IF
  
    END IF

!   =========================================================================

END SUBROUTINE rhscal
