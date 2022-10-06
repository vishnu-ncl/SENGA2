SUBROUTINE rhscal

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga


! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:26:16

!     *************************************************************************

!     RHSCAL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     12-NOV-2002:  CREATED
!     26-OCT-2008:  RSC/TDD BUG FIX FZLCON
!     08-AUG-2012:  RSC EVALUATE ALL SPECIES
!     17-APR-2013:  RSC MIXTURE AVERAGED TRANSPORT
!     14-JUL-2013:  RSC RADIATION HEAT LOSS
!     08-JUN-2015:  RSC REMOVE Nth SPECIES TREATMENT
!     08-JUN-2015:  RSC UPDATED WALL BCS

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     COMPUTES RIGHT-HAND-SIDES FOR TIME INTEGRATION OF SCALAR PDEs
!     INCLUDES MULTIPLE SCALARS AND MULTI-STEP CHEMISTRY
!     ENERGY EQUATION REQUIRES PRESSURE-WORK AND VISCOUS WORK TERMS
!     COMPUTED IN RHSVEL

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
real(kind=dp) :: ctrans(nspcmx)
real(kind=dp) :: fornow,combo1,combo2,combo3
INTEGER :: ic,jc,kc,ispec
INTEGER :: itint,icp,iindex,ipower,icoef1,icoef2
LOGICAL :: flmtds
INTEGER :: rangexyz(6)

!     BEGIN
!     =====

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     EVALUATE THE TEMPERATURE
!     ------------------------
!     ALSO PRESSURE, MIXTURE CP AND MIXTURE GAS CONSTANT
CALL temper
!                                                               PRUN,TRUN = P,T
!                                                           STORE7 = RHO*MIX RG
!     =========================================================================

!     COLLECT MIXTURE CP AND GAS CONSTANT FOR BCs
!     -------------------------------------------

!     X-DIRECTION
IF(fxlcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      strgxl(jc,kc) = transp(istal,jc,kc)
      strrxl(jc,kc) = store7(istal,jc,kc)/drhs(istal,jc,kc)
      
    END DO
  END DO
END IF
IF(fxrcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      strgxr(jc,kc) = transp(istol,jc,kc)
      strrxr(jc,kc) = store7(istol,jc,kc)/drhs(istol,jc,kc)
      
    END DO
  END DO
END IF

!     Y-DIRECTION
IF(fylcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      strgyl(ic,kc) = transp(ic,jstal,kc)
      strryl(ic,kc) = store7(ic,jstal,kc)/drhs(ic,jstal,kc)
      
    END DO
  END DO
END IF
IF(fyrcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      strgyr(ic,kc) = transp(ic,jstol,kc)
      strryr(ic,kc) = store7(ic,jstol,kc)/drhs(ic,jstol,kc)
      
    END DO
  END DO
END IF

!     Z-DIRECTION
IF(fzlcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      strgzl(ic,jc) = transp(ic,jc,kstal)
      strrzl(ic,jc) = store7(ic,jc,kstal)/drhs(ic,jc,kstal)
      
    END DO
  END DO
END IF
IF(fzrcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      strgzr(ic,jc) = transp(ic,jc,kstol)
      strrzr(ic,jc) = store7(ic,jc,kstol)/drhs(ic,jc,kstol)
      
    END DO
  END DO
END IF
!                                                              ALL STORES CLEAR
!     =========================================================================

!     MASS FLUX DIVERGENCE
!     --------------------
!     URHS,VRHS,WRHS CONTAIN RHO U, RHO V, RHO W

CALL dfbydx(d_urhs,d_store1)
CALL dfbydy(d_vrhs,d_store2)
CALL dfbydz(d_wrhs,d_store3)

    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    
    call ops_par_loop(compute_kernel_AequalsBCDplus, "A=B+C+D", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_divm, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ)) 

!                                                              ALL STORES CLEAR
!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     INTERNAL ENERGY EQUATION
!     ========================

!     CONVERT INTERNAL ENERGY
!     -----------------------

!     ERHS CONTAINS RHO E: CONVERT TO E
!     E IS PARALLEL
    rangexyz = (/istalt,istolt,jstalt,jstolt,kstalt,kstolt/)
    call ops_par_loop(compute_kernel_AequalsAdivB, "A=A/B", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ))

!                                                              ALL STORES CLEAR
!     =========================================================================

!     COLLECT INTERNAL ENERGY FOR BCs
!     -------------------------------

!     X-DIRECTION
IF(fxlcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      strexl(jc,kc) = erhs(istal,jc,kc)
      
    END DO
  END DO
END IF
IF(fxrcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      strexr(jc,kc) = erhs(istol,jc,kc)
      
    END DO
  END DO
END IF

!     Y-DIRECTION
IF(fylcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      streyl(ic,kc) = erhs(ic,jstal,kc)
      
    END DO
  END DO
END IF
IF(fyrcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      streyr(ic,kc) = erhs(ic,jstol,kc)
      
    END DO
  END DO
END IF

!     Z-DIRECTION
IF(fzlcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      strezl(ic,jc) = erhs(ic,jc,kstal)
      
    END DO
  END DO
END IF
IF(fzrcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      strezr(ic,jc) = erhs(ic,jc,kstol)
      
    END DO
  END DO
END IF
!                                                              ALL STORES CLEAR
!     =========================================================================

!     E EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     HALF E DIV RHO U

!     COLLECT E DIV RHO U IN STORE4 FOR NOW
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(compute_kernel_AequalsBmulC, "A=B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_divm, 1, s3d_000, "real(dp)", OPS_READ))
    
!                                                          STORE4 = E DIV RHO U
!     =========================================================================

!     E EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     HALF DIV RHO U E

!     D/DX RHO U E
!     RHO U E IS PARALLEL
    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(compute_kernel_AequalsBmulC, "A=B*C", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_READ))

CALL dfbydx(d_store7,d_store1)

!     D/DY RHO V E
!     RHO V E IS PARALLEL
    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(compute_kernel_AequalsBmulC, "A=B*C", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_READ))
    
CALL dfbydy(d_store7,d_store2)

!     D/DZ RHO W E
!     RHO W E IS PARALLEL
    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(compute_kernel_AequalsBmulC, "A=B*C", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_READ))

CALL dfbydz(d_store7,d_store3)

!     COLLECT DIV RHO U E IN STORE4 FOR NOW
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(compute_kernel_AequalsABCDplus, "A=A+B+C+D", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store3, 1, s3d_000, "real(dp)", OPS_READ))

!                                            STORE4 = E DIV RHO U + DIV RHO U E
!     =========================================================================

!     E EQUATION: CONVECTIVE TERMS
!     ----------------------------
!     HALF RHO U.DEL E

CALL dfbydx(d_erhs,d_store1)
CALL dfbydy(d_erhs,d_store2)
CALL dfbydz(d_erhs,d_store3)

!     COLLECT ALL CONVECTIVE TERMS IN ERHS
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = -half*(store4(ic,jc,kc)  &
          + store1(ic,jc,kc)*urhs(ic,jc,kc) + store2(ic,jc,kc)*vrhs(ic,jc,kc)  &
          + store3(ic,jc,kc)*wrhs(ic,jc,kc))
      
    END DO
  END DO
END DO

!     -------------------------------------
!     E EQUATION: CONVECTIVE TERMS COMPLETE
!     -------------------------------------
!                                                              ALL STORES CLEAR
!     =========================================================================

!     E-EQUATION: HEAT FLUX TERMS
!     ---------------------------

!     TEMPERATURE GRADIENTS
CALL dfbydx(d_trun,d_store1)
CALL dfbydy(d_trun,d_store2)
CALL dfbydz(d_trun,d_store3)

!                                                         STORE1,2,3 = DTDX,Y,Z
!     =========================================================================

!     COLLECT TEMPERATURE AND ITS GRADIENTS FOR BCs
!     ---------------------------------------------

!     X-DIRECTION
IF(fxlcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      strtxl(jc,kc) = trun(istal,jc,kc)
      bcltxl(jc,kc) = store1(istal,jc,kc)
      
    END DO
  END DO
END IF
IF(fxrcnv)THEN
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      strtxr(jc,kc) = trun(istol,jc,kc)
      bcltxr(jc,kc) = store1(istol,jc,kc)
      
    END DO
  END DO
END IF

!     Y-DIRECTION
IF(fylcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      strtyl(ic,kc) = trun(ic,jstal,kc)
      bcltyl(ic,kc) = store2(ic,jstal,kc)
      
    END DO
  END DO
END IF
IF(fyrcnv)THEN
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      strtyr(ic,kc) = trun(ic,jstol,kc)
      bcltyr(ic,kc) = store2(ic,jstol,kc)
      
    END DO
  END DO
END IF

!     Z-DIRECTION
IF(fzlcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      strtzl(ic,jc) = trun(ic,jc,kstal)
      bcltzl(ic,jc) = store3(ic,jc,kstal)
      
    END DO
  END DO
END IF
IF(fzrcnv)THEN
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      strtzr(ic,jc) = trun(ic,jc,kstol)
      bcltzr(ic,jc) = store3(ic,jc,kstol)
      
    END DO
  END DO
END IF
!                                                         STORE1,2,3 = DTDX,Y,Z
!     =========================================================================

!     E-EQUATION: HEAT FLUX TERMS
!     ---------------------------

!     THERMAL CONDUCTIVITY
!     ANALYTICAL FUNCTION OF TEMPERATURE
!     TRANSP CONTAINS MIXTURE CP
!     STORE CONDUCTIVITY/CP IN TRANSP FOR USE IN DIFFUSIVITY AND VISCOSITY
!     STORE CONDUCTIVITY IN STORE7 FOR NOW

!     THERMAL CONDUCTIVITY IS PARALLEL
DO kc = kstab,kstob
  DO jc = jstab,jstob
    DO ic = istab,istob
      
      fornow = alamda*EXP(rlamda*LOG(trun(ic,jc,kc)))
      store7(ic,jc,kc) = fornow*transp(ic,jc,kc)
      transp(ic,jc,kc) = fornow
      
    END DO
  END DO
END DO

!     MIXTURE AVERAGED TRANSPORT
!     RSC 17-APR-2013
!     THERMAL CONDUCTIVITY

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
          fornow = yrhs(ic,jc,kc,ispec)*ovwmol(ispec)
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

!     CONDUCTIVITY GRADIENTS
CALL dfbydx(d_store7,d_store4)
CALL dfbydy(d_store7,d_store5)
CALL dfbydz(d_store7,d_store6)

!     BOUNDARY CONDITIONS
!     BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fxlcon)CALL zeroxl(d_store4)
IF(fxrcon)CALL zeroxr(d_store4)
!     BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fylcon)CALL zeroyl(d_store5)
IF(fyrcon)CALL zeroyr(d_store5)
!     BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
!     RSC/TDD BUG FIX FZLCON
IF(fzlcon)CALL zerozl(d_store6)
IF(fzrcon)CALL zerozr(d_store6)

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc) +(store4(ic,jc,kc)*store1(ic,jc,kc)  &
          + store5(ic,jc,kc)*store2(ic,jc,kc) + store6(ic,jc,kc)*store3(ic,jc,kc))
    END DO
  END DO
END DO
!                                                         STORE1,2,3 = DTDX,Y,Z
!                                                         STORE7 = CONDUCTIVITY
!     =========================================================================

!     E-EQUATION: HEAT FLUX TERMS
!     ---------------------------
!     WALL BC: THERMAL CONDUCTION TERMS
IF(fxlcnw)THEN

    rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
    call ops_par_loop(heat_flux_kernel_thermal_fxlcnw, "HEAT FLUX: Thermal fxlcnw", senga_grid, 3, rangexyz,  &
                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                ops_arg_dat(d_store1, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                ops_arg_dat(d_store7, 1, s3d_p100_to_p400_x, "real(dp)", OPS_READ), &
                ops_arg_gbl(acbcxl, 1, "real(dp)", OPS_READ))

END IF
IF(fxrcnw)THEN

    rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(heat_flux_kernel_thermal_fxrcnw, "HEAT FLUX: Thermal fxrcnw", senga_grid, 3, rangexyz,  &
                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                ops_arg_dat(d_store1, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                ops_arg_dat(d_store7, 1, s3d_m100_to_m400_x, "real(dp)", OPS_READ), &
                ops_arg_gbl(acbcxr, 1, "real(dp)", OPS_READ))

END IF
IF(fylcnw)THEN

    rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
    call ops_par_loop(heat_flux_kernel_thermal_fylcnw, "HEAT FLUX: Thermal fylcnw", senga_grid, 3, rangexyz,  &
                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                ops_arg_dat(d_store2, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                ops_arg_dat(d_store7, 1, s3d_p010_to_p040_y, "real(dp)", OPS_READ), &
                ops_arg_gbl(acbcyl, 1, "real(dp)", OPS_READ))

END IF
IF(fyrcnw)THEN

    rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
    call ops_par_loop(heat_flux_kernel_thermal_fyrcnw, "HEAT FLUX: Thermal fyrcnw", senga_grid, 3, rangexyz,  &
                ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                ops_arg_dat(d_store2, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                ops_arg_dat(d_store7, 1, s3d_m010_to_m040_y, "real(dp)", OPS_READ), &
                ops_arg_gbl(acbcyr, 1, "real(dp)", OPS_READ))

END IF
IF(fzlcnw)THEN

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

!     =========================================================================

!     E-EQUATION: HEAT FLUX TERMS
!     ---------------------------
!     SECOND DERIVATIVE TERMS

!     TEMPERATURE SECOND DERIVATIVES
CALL d2fdx2(d_trun,d_store1)
CALL d2fdy2(d_trun,d_store2)
CALL d2fdz2(d_trun,d_store3)

!     BOUNDARY CONDITIONS
!     BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fxlcon)CALL zeroxl(d_store1)
IF(fxrcon)CALL zeroxr(d_store1)
!     BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fylcon)CALL zeroyl(d_store2)
IF(fyrcon)CALL zeroyr(d_store2)
!     BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
!     RSC 28-JUN-2015 BUG FIX FZLCON
IF(fzlcon)CALL zerozl(d_store3)
IF(fzrcon)CALL zerozr(d_store3)

!     COLLECT CONDUCTIVITY TERMS
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc) +(store1(ic,jc,kc)  &
          + store2(ic,jc,kc) + store3(ic,jc,kc))  &
          *store7(ic,jc,kc)
      
    END DO
  END DO
END DO

!     ---------------------------------------------------
!     E-EQUATION: FURTHER HEAT FLUX TERMS EVALUATED BELOW
!     ---------------------------------------------------
!     E-EQUATION: PRESSURE-WORK AND VISCOUS WORK TERMS
!                 EVALUATED IN SUBROUTINE RHSVEL
!     ---------------------------------------------------
!                                                              ALL STORES CLEAR
!     =========================================================================

!     E-EQUATION: RADIATION HEAT LOSS
!     -------------------------------
IF(flradn)CALL radcal

!     =========================================================================

!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     SPECIES MASS FRACTION EQUATIONS
!     ===============================

!     REACTION RATE FOR ALL SPECIES
!     -----------------------------
CALL chrate
!---UA
DO ispec = 1,nspec
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        rrte(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec)
      END DO
    END DO
  END DO
END DO
!---end-UA
!                                                          RATE = REACTION RATE
!     =========================================================================

!     COLLECT REACTION RATE FOR BCs
!     -----------------------------

!     X-DIRECTION
IF(fxlcnv)THEN
  DO ispec = 1,nspec
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        ratexl(jc,kc,ispec) = rate(istal,jc,kc,ispec)
        
      END DO
    END DO
  END DO
END IF
IF(fxrcnv)THEN
  DO ispec = 1,nspec
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        ratexr(jc,kc,ispec) = rate(istol,jc,kc,ispec)
        
      END DO
    END DO
  END DO
END IF

!     Y-DIRECTION
IF(fylcnv)THEN
  DO ispec = 1,nspec
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        rateyl(ic,kc,ispec) = rate(ic,jstal,kc,ispec)
        
      END DO
    END DO
  END DO
END IF
IF(fyrcnv)THEN
  DO ispec = 1,nspec
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        rateyr(ic,kc,ispec) = rate(ic,jstol,kc,ispec)
        
      END DO
    END DO
  END DO
END IF

!     Z-DIRECTION
IF(fzlcnv)THEN
  DO ispec = 1,nspec
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        ratezl(ic,jc,ispec) = rate(ic,jc,kstal,ispec)
        
      END DO
    END DO
  END DO
END IF
IF(fzrcnv)THEN
  DO ispec = 1,nspec
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        ratezr(ic,jc,ispec) = rate(ic,jc,kstol,ispec)
        
      END DO
    END DO
  END DO
END IF
!                                                          RATE = REACTION RATE
!     =========================================================================

!     ZERO THE ACCUMULATORS FOR THE DIFFUSION CORRECTION VELOCITY
!     AND ITS DIVERGENCE

    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_WRITE))

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_WRITE))
    

!     ZERO THE ACCUMULATOR FOR THE MIXTURE ENTHALPY
!     MIXTURE H IS PARALLEL

    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_wtmp, 1, s3d_000, "real(dp)", OPS_WRITE))

!                                                          RATE = REACTION RATE
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!     =========================================================================

!     MIXTURE AVERAGED TRANSPORT
!     RSC 17-APR-2013
!     EVALUATE FIRST AND SECOND DERIVATIVES
!     OF LN(MIXTURE MOLAR MASS), LN(PRESSURE) AND LN(TEMPERATURE)

!     MIXTURE MOLAR MASS
IF(flmixw)THEN
  
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = LOG(wmomix(ic,jc,kc))
        
      END DO
    END DO
  END DO   
   
  CALL dfbydx(d_store7,d_wd1x)
  CALL dfbydy(d_store7,d_wd1y)
  CALL dfbydz(d_store7,d_wd1z)
  CALL d2fdx2(d_store7,d_wd2x)
  CALL d2fdy2(d_store7,d_wd2y)
  CALL d2fdz2(d_store7,d_wd2z)
  
END IF

!     PRESSURE
IF(flmixp)THEN
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = LOG(prun(ic,jc,kc))
        
      END DO
    END DO
  END DO
  
  CALL dfbydx(d_store7,d_pd1x)
  CALL dfbydy(d_store7,d_pd1y)
  CALL dfbydz(d_store7,d_pd1z)
  CALL d2fdx2(d_store7,d_pd2x)
  CALL d2fdy2(d_store7,d_pd2y)
  CALL d2fdz2(d_store7,d_pd2z)
  
END IF

!     TEMPERATURE
IF(flmixt)THEN
  
!       TRANSP CONTAINS LN(T/TDIFGB)
  CALL dfbydx(d_transp,d_td1x)
  CALL dfbydy(d_transp,d_td1y)
  CALL dfbydz(d_transp,d_td1z)
  CALL d2fdx2(d_transp,d_td2x)
  CALL d2fdy2(d_transp,d_td2y)
  CALL d2fdz2(d_transp,d_td2z)
  
END IF

!     =========================================================================

!     RUN THROUGH ALL SPECIES
!     -----------------------
!     RSC 08-AUG-2012 EVALUATE ALL SPECIES
!     RSC 08-JUN-2015 REMOVE Nth SPECIES TREATMENT
DO ispec = 1,nspec
  
!       =======================================================================
  
!       YRHS CONTAINS RHO Y: CONVERT TO Y
!       Y IS PARALLEL
  DO kc = kstalt,kstolt
    DO jc = jstalt,jstolt
      DO ic = istalt,istolt
        
        yrhs(ic,jc,kc,ispec) = yrhs(ic,jc,kc,ispec)/drhs(ic,jc,kc)
        
      END DO
    END DO
  END DO
  
!       =======================================================================
  
!       Y EQUATION: CONVECTIVE TERMS
!       ----------------------------
!       HALF Y DIV RHO U
  
!       COLLECT Y SOURCE TERMS IN RATE FOR NOW
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec)  &
            - half*yrhs(ic,jc,kc,ispec)*divm(ic,jc,kc)
        
      END DO
    END DO
  END DO
!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!       =======================================================================
  
!       Y EQUATION: CONVECTIVE TERMS
!       ----------------------------
!       HALF DIV RHO U Y
  
!       D/DX RHO U Y
!       RHO U Y IS PARALLEL
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = yrhs(ic,jc,kc,ispec)*urhs(ic,jc,kc)
        
      END DO
    END DO
  END DO
  CALL dfbydx(d_store7,d_store1)
  
!       D/DY RHO V Y
!       RHO V Y IS PARALLEL
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = yrhs(ic,jc,kc,ispec)*vrhs(ic,jc,kc)
        
      END DO
    END DO
  END DO
  CALL dfbydy(d_store7,d_store2)
  
!       D/DZ RHO W Y
!       RHO W Y IS PARALLEL
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = yrhs(ic,jc,kc,ispec)*wrhs(ic,jc,kc)
        
      END DO
    END DO
  END DO
  CALL dfbydz(d_store7,d_store3)
  
!       COLLECT DIV RHO U Y IN RATE FOR NOW
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec)  &
            - half*(store1(ic,jc,kc) + store2(ic,jc,kc)  &
            + store3(ic,jc,kc))
        
      END DO
    END DO
  END DO
!                                                  RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!       =======================================================================
  
!       SPECIES MASS FRACTION GRADIENT TERMS
!       ------------------------------------
  
!       SPECIES MASS FRACTION GRADIENTS
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = yrhs(ic,jc,kc,ispec)
        
      END DO
    END DO
  END DO
  CALL dfbydx(d_store7,d_store1)
  CALL dfbydy(d_store7,d_store2)
  CALL dfbydz(d_store7,d_store3)
!                                                         STORE1,2,3 = DYDX,Y,Z
!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!       =======================================================================
  
!       COLLECT SPECIES MASS FRACTION AND ITS GRADIENTS FOR BCs
!       -------------------------------------------------------
  
!       X-DIRECTION: DYDX
  IF(fxlcnv)THEN
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        stryxl(jc,kc,ispec) = yrhs(istal,jc,kc,ispec)
        bclyxl(jc,kc,ispec) = store1(istal,jc,kc)
        
      END DO
    END DO
  END IF
  IF(fxrcnv)THEN
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        stryxr(jc,kc,ispec) = yrhs(istol,jc,kc,ispec)
        bclyxr(jc,kc,ispec) = store1(istol,jc,kc)
        
      END DO
    END DO
  END IF
  
!       Y-DIRECTION: DYDY
  IF(fylcnv)THEN
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        stryyl(ic,kc,ispec) = yrhs(ic,jstal,kc,ispec)
        bclyyl(ic,kc,ispec) = store2(ic,jstal,kc)
        
      END DO
    END DO
  END IF
  IF(fyrcnv)THEN
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        stryyr(ic,kc,ispec) = yrhs(ic,jstol,kc,ispec)
        bclyyr(ic,kc,ispec) = store2(ic,jstol,kc)
        
      END DO
    END DO
  END IF
  
!       Z-DIRECTION: DYDZ
  IF(fzlcnv)THEN
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        stryzl(ic,jc,ispec) = yrhs(ic,jc,kstal,ispec)
        bclyzl(ic,jc,ispec) = store3(ic,jc,kstal)
        
      END DO
    END DO
  END IF
  IF(fzrcnv)THEN
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        stryzr(ic,jc,ispec) = yrhs(ic,jc,kstol,ispec)
        bclyzr(ic,jc,ispec) = store3(ic,jc,kstol)
        
      END DO
    END DO
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
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec)  &
            - half*(store1(ic,jc,kc)*urhs(ic,jc,kc)  &
            + store2(ic,jc,kc)*vrhs(ic,jc,kc) + store3(ic,jc,kc)*wrhs(ic,jc,kc))
        
      END DO
    END DO
  END DO
  
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
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = transp(ic,jc,kc)*olewis(ispec)
        
      END DO
    END DO
  END DO
!                                                         STORE1,2,3 = DYDX,Y,Z
!                                                          STORE7 = DIFFUSIVITY
!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                       WTMP = MIXTURE H
!       -----------------------------------------------------------------------
  
!       MIXTURE AVERAGED TRANSPORT
!       RSC 17-APR-2013
!       TRANSP CONTAINS LN(T/TDIFGB)
  IF(flmavt)THEN
    
!         MASS DIFFUSIVITY FOR EACH SPECIES
!         RELATIVE TO CURRENT SPECIES
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
            fornow = yrhs(ic,jc,kc,jspec) + dfctol
            combo1 = combo1 + fornow
            combo2 = combo2 + fornow*ovwmol(jspec)/ctrans(jspec)
          END DO
          fornow = yrhs(ic,jc,kc,ispec) + dfctol
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
    
!         THERMAL DIFFUSION RATIO FOR EACH SPECIES
!         RELATIVE TO CURRENT SPECIES
  
    rangexyz = (/istab,istob,jstab,jstob,kstab,kstob/)
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_tdrmix, 1, s3d_000, "real(dp)", OPS_WRITE))
    
    DO jspec = 1, nspec
      
      flmtds = flmtdr(jspec).AND.(ispec /= jspec)
      IF(flmtds)THEN
        
        DO kc = kstab,kstob
          DO jc = jstab,jstob
            DO ic = istab,istob
              
!                   THERMAL DIFFUSION RATIO FOR THIS SPECIES PAIR
              combo2 = trun(ic,jc,kc)/tdifgb
              fornow = tdrcco(ncotdr,jspec,ispec)
              DO icp = ncotm1,1,-1
                fornow = fornow*combo2 + tdrcco(icp,jspec,ispec)
              END DO
              ctrans(jspec) = fornow
              
!                   COMBINATION RULE FOR THERMAL DIFFUSIION RATIO
              fornow = yrhs(ic,jc,kc,jspec)*ovwmol(jspec)
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
    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &    
                    ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ))

    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ))

    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
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
        wtmp(ic,jc,kc) = wtmp(ic,jc,kc) + utmp(ic,jc,kc)*yrhs(ic,jc,kc,ispec)
        
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
  IF(fxlcnv)THEN
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        strhxl(jc,kc,ispec) = utmp(istal,jc,kc)
        
      END DO
    END DO
  END IF
  IF(fxrcnv)THEN
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        strhxr(jc,kc,ispec) = utmp(istol,jc,kc)
        
      END DO
    END DO
  END IF
  
!       Y-DIRECTION
  IF(fylcnv)THEN
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        strhyl(ic,kc,ispec) = utmp(ic,jstal,kc)
        
      END DO
    END DO
  END IF
  IF(fyrcnv)THEN
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        strhyr(ic,kc,ispec) = utmp(ic,jstol,kc)
        
      END DO
    END DO
  END IF
  
!       Z-DIRECTION
  IF(fzlcnv)THEN
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        strhzl(ic,jc,ispec) = utmp(ic,jc,kstal)
        
      END DO
    END DO
  END IF
  IF(fzrcnv)THEN
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        strhzr(ic,jc,ispec) = utmp(ic,jc,kstol)
        
      END DO
    END DO
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
  IF(flmduf(ispec))THEN
    
    DO kc = kstab,kstob
      DO jc = jstab,jstob
        DO ic = istab,istob
          
          utmp(ic,jc,kc) = utmp(ic,jc,kc)  &
              + rgspec(ispec)*trun(ic,jc,kc)*tdrmix(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
!       Y EQUATION: DIFFUSIVE TERMS
!       ---------------------------
!       E EQUATION: FURTHER HEAT FLUX TERMS
!       DIFFUSIVITY GRADIENT TERMS
  
!       DIFFUSIVITY GRADIENTS
  CALL dfbydx(d_store7,d_store4)
  CALL dfbydy(d_store7,d_store5)
  CALL dfbydz(d_store7,d_store6)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fxldif)CALL zeroxl(d_store4)
  IF(fxrdif)CALL zeroxr(d_store4)
!       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fyldif)CALL zeroyl(d_store5)
  IF(fyrdif)CALL zeroyr(d_store5)
!       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fzldif)CALL zerozl(d_store6)
  IF(fzrdif)CALL zerozr(d_store6)
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        fornow = store4(ic,jc,kc)*store1(ic,jc,kc)  &
            + store5(ic,jc,kc)*store2(ic,jc,kc)  &
            + store6(ic,jc,kc)*store3(ic,jc,kc)
        
!             Y EQUATION
        rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec) + fornow
        
!             DIFFUSION CORRECTION VELOCITY DIVERGENCE
        vtmp(ic,jc,kc) = vtmp(ic,jc,kc) + fornow
        
      END DO
    END DO
  END DO
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF(fxladb)CALL zeroxl(d_store4)
  IF(fxradb)CALL zeroxr(d_store4)
!       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF(fyladb)CALL zeroyl(d_store5)
  IF(fyradb)CALL zeroyr(d_store5)
!       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF(fzladb)CALL zerozl(d_store6)
  IF(fzradb)CALL zerozr(d_store6)
  
!       E EQUATION
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        fornow = store4(ic,jc,kc)*store1(ic,jc,kc)  &
            + store5(ic,jc,kc)*store2(ic,jc,kc)  &
            + store6(ic,jc,kc)*store3(ic,jc,kc)
        
        erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*utmp(ic,jc,kc)
        
      END DO
    END DO
  END DO
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
  CALL dfbydx(d_utmp,d_store4)
  CALL dfbydy(d_utmp,d_store5)
  CALL dfbydz(d_utmp,d_store6)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fxldif)CALL zeroxl(d_store4)
  IF(fxrdif)CALL zeroxr(d_store4)
!       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fyldif)CALL zeroyl(d_store5)
  IF(fyrdif)CALL zeroyr(d_store5)
!       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fzldif)CALL zerozl(d_store6)
  IF(fzrdif)CALL zerozr(d_store6)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF(fxladb)CALL zeroxl(d_store4)
  IF(fxradb)CALL zeroxr(d_store4)
!       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF(fyladb)CALL zeroyl(d_store5)
  IF(fyradb)CALL zeroyr(d_store5)
!       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF(fzladb)CALL zerozl(d_store6)
  IF(fzradb)CALL zerozr(d_store6)
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        fornow = store4(ic,jc,kc)*store1(ic,jc,kc)  &
            + store5(ic,jc,kc)*store2(ic,jc,kc)  &
            + store6(ic,jc,kc)*store3(ic,jc,kc)
        
        erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*store7(ic,jc,kc)
        
      END DO
    END DO
  END DO
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
  IF(fxldfw)THEN
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        fornow = zero
        DO ic = istap1,istow
          
          fornow = fornow + acbcxl(ic-1)*store7(ic,jc,kc)*store1(ic,jc,kc)
          
        END DO
        rate(istal,jc,kc,ispec) = rate(istal,jc,kc,ispec) + fornow
        vtmp(istal,jc,kc) = vtmp(istal,jc,kc) + fornow
        erhs(istal,jc,kc) = erhs(istal,jc,kc) + fornow*utmp(istal,jc,kc)
        
      END DO
    END DO
  END IF
  IF(fxrdfw)THEN
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        fornow = zero
        DO ic = istaw,istom1
          
          fornow = fornow + acbcxr(istol-ic)*store7(ic,jc,kc)*store1(ic,jc,kc)
          
        END DO
        rate(istol,jc,kc,ispec) = rate(istol,jc,kc,ispec) + fornow
        vtmp(istol,jc,kc) = vtmp(istol,jc,kc) + fornow
        erhs(istol,jc,kc) = erhs(istol,jc,kc) + fornow*utmp(istol,jc,kc)
        
      END DO
    END DO
  END IF
  IF(fyldfw)THEN
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        fornow = zero
        DO jc = jstap1,jstow
          
          fornow = fornow + acbcyl(jc-1)*store7(ic,jc,kc)*store2(ic,jc,kc)
          
        END DO
        rate(ic,jstal,kc,ispec) = rate(ic,jstal,kc,ispec) + fornow
        vtmp(ic,jstal,kc) = vtmp(ic,jstal,kc) + fornow
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc) + fornow*utmp(ic,jstal,kc)
        
      END DO
    END DO
  END IF
  IF(fyrdfw)THEN
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        fornow = zero
        DO jc = jstaw,jstom1
          
          fornow = fornow + acbcyr(jstol-jc)*store7(ic,jc,kc)*store2(ic,jc,kc)
          
        END DO
        rate(ic,jstol,kc,ispec) = rate(ic,jstol,kc,ispec) + fornow
        vtmp(ic,jstol,kc) = vtmp(ic,jstol,kc) + fornow
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc) + fornow*utmp(ic,jstol,kc)
        
      END DO
    END DO
  END IF
  IF(fzldfw)THEN
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        fornow = zero
        DO kc = kstap1,kstow
          
          fornow = fornow + acbczl(kc-1)*store7(ic,jc,kc)*store3(ic,jc,kc)
          
        END DO
        rate(ic,jc,kstal,ispec) = rate(ic,jc,kstal,ispec) + fornow
        vtmp(ic,jc,kstal) = vtmp(ic,jc,kstal) + fornow
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal) + fornow*utmp(ic,jc,kstal)
        
      END DO
    END DO
  END IF
  IF(fzrdfw)THEN
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        fornow = zero
        DO kc = kstaw,kstom1
          
          fornow = fornow + acbczr(kstol-kc)*store7(ic,jc,kc)*store3(ic,jc,kc)
          
        END DO
        rate(ic,jc,kstol,ispec) = rate(ic,jc,kstol,ispec) + fornow
        vtmp(ic,jc,kstol) = vtmp(ic,jc,kstol) + fornow
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol) + fornow*utmp(ic,jc,kstol)
        
      END DO
    END DO
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
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = yrhs(ic,jc,kc,ispec)
        
      END DO
    END DO
  END DO
  CALL d2fdx2(d_store7,d_store1)
  CALL d2fdy2(d_store7,d_store2)
  CALL d2fdz2(d_store7,d_store3)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fxldif)CALL zeroxl(d_store1)
  IF(fxrdif)CALL zeroxr(d_store1)
!       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fyldif)CALL zeroyl(d_store2)
  IF(fyrdif)CALL zeroyr(d_store2)
!       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fzldif)CALL zerozl(d_store3)
  IF(fzrdif)CALL zerozr(d_store3)
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        fornow = (store1(ic,jc,kc) +  store2(ic,jc,kc)  &
            +  store3(ic,jc,kc))*store4(ic,jc,kc)
        
!             Y EQUATION
        rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec) + fornow
        
!             DIFFUSION CORRECTION VELOCITY DIVERGENCE
        vtmp(ic,jc,kc) = vtmp(ic,jc,kc) + fornow
        
      END DO
    END DO
  END DO
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF(fxladb)CALL zeroxl(d_store1)
  IF(fxradb)CALL zeroxr(d_store1)
!       BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF(fyladb)CALL zeroyl(d_store2)
  IF(fyradb)CALL zeroyr(d_store2)
!       BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
  IF(fzladb)CALL zerozl(d_store3)
  IF(fzradb)CALL zerozr(d_store3)
  
!       E EQUATION
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        fornow = (store1(ic,jc,kc) +  store2(ic,jc,kc)  &
            +  store3(ic,jc,kc))*store4(ic,jc,kc)
        
        erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*utmp(ic,jc,kc)
        
      END DO
    END DO
  END DO
!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                              WTMP = MIXTURE H
!       =======================================================================
  
!       MIXTURE AVERAGED TRANSPORT
!       RSC 23-APR-2013
!       MOLAR MASS TERMS, PRESSURE TERMS, SORET EFFECT
  
!       MIXTURE MOLAR MASS TERMS
  IF(flmixw)THEN
    
!         FIRST AND SECOND DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED
    
    DO kc = kstab,kstob
      DO jc = jstab,jstob
        DO ic = istab,istob
          
          store7(ic,jc,kc) = difmix(ic,jc,kc)*yrhs(ic,jc,kc,ispec)
          
        END DO
      END DO
    END DO
    
!         DIFFUSION CORRECTION VELOCITY
!         FIRST DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED
   rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_wd1x, 1, s3d_000, "real(dp)", OPS_READ))

    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_wd1y, 1, s3d_000, "real(dp)", OPS_READ))

    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_wd1z, 1, s3d_000, "real(dp)", OPS_READ)) 
    
!         Y EQUATION: DIFFUSIVE TERMS
!         E EQUATION: FURTHER HEAT FLUX TERMS
    
!         DIFFUSIVITY GRADIENT TERMS
    
!         DIFFUSIVITY GRADIENTS
    CALL dfbydx(d_store7,d_store1)
    CALL dfbydy(d_store7,d_store2)
    CALL dfbydz(d_store7,d_store3)
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif)CALL zeroxl(d_store1)
    IF(fxrdif)CALL zeroxr(d_store1)
!         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif)CALL zeroyl(d_store2)
    IF(fyrdif)CALL zeroyr(d_store2)
!         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif)CALL zerozl(d_store3)
    IF(fzrdif)CALL zerozr(d_store3)
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = store1(ic,jc,kc)*wd1x(ic,jc,kc)  &
              + store2(ic,jc,kc)*wd1y(ic,jc,kc) + store3(ic,jc,kc)*wd1z(ic,jc,kc)
          
!               Y EQUATION
          rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec) + fornow
          
!               DIFFUSION CORRECTION VELOCITY DIVERGENCE
          vtmp(ic,jc,kc) = vtmp(ic,jc,kc) + fornow
          
        END DO
      END DO
    END DO
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb)CALL zeroxl(d_store1)
    IF(fxradb)CALL zeroxr(d_store1)
!         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb)CALL zeroyl(d_store2)
    IF(fyradb)CALL zeroyr(d_store2)
!         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb)CALL zerozl(d_store3)
    IF(fzradb)CALL zerozr(d_store3)
    
!         E EQUATION
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = store1(ic,jc,kc)*wd1x(ic,jc,kc)  &
              + store2(ic,jc,kc)*wd1y(ic,jc,kc) + store3(ic,jc,kc)*wd1z(ic,jc,kc)
          
          erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*utmp(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
    
!         E-EQUATION: FURTHER HEAT FLUX TERMS
!         SPECIES ENTHALPY GRADIENT TERMS
    
!         SPECIES ENTHALPY GRADIENTS ALREADY IN STORE4,5,6
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif)CALL zeroxl(d_store4)
    IF(fxrdif)CALL zeroxr(d_store4)
!         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif)CALL zeroyl(d_store5)
    IF(fyrdif)CALL zeroyr(d_store5)
!         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif)CALL zerozl(d_store6)
    IF(fzrdif)CALL zerozr(d_store6)
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb)CALL zeroxl(d_store4)
    IF(fxradb)CALL zeroxr(d_store4)
!         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb)CALL zeroyl(d_store5)
    IF(fyradb)CALL zeroyr(d_store5)
!         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb)CALL zerozl(d_store6)
    IF(fzradb)CALL zerozr(d_store6)
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = store4(ic,jc,kc)*wd1x(ic,jc,kc)  &
              + store5(ic,jc,kc)*wd1y(ic,jc,kc) + store6(ic,jc,kc)*wd1z(ic,jc,kc)
          
          erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*store7(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
!         ====================================================================
    
!         Y-EQUATION: DIFFUSIVE TERMS
!         ---------------------------
!         WALL BC: MOLAR MASS TERMS
!         E-EQUATION: HEAT FLUX TERMS
!         WALL BC: ENTHALPY DIFFUSION TERMS
    IF(fxldfw)THEN
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = zero
          DO ic = istap1,istow
            
            fornow = fornow + acbcxl(ic-1)*store7(ic,jc,kc)*wd1x(ic,jc,kc)
            
          END DO
          rate(istal,jc,kc,ispec) = rate(istal,jc,kc,ispec) + fornow
          vtmp(istal,jc,kc) = vtmp(istal,jc,kc) + fornow
          erhs(istal,jc,kc) = erhs(istal,jc,kc) + fornow*utmp(istal,jc,kc)
          
        END DO
      END DO
    END IF
    IF(fxrdfw)THEN
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = zero
          DO ic = istaw,istom1
            
            fornow = fornow + acbcxr(istol-ic)*store7(ic,jc,kc)*wd1x(ic,jc,kc)
            
          END DO
          rate(istol,jc,kc,ispec) = rate(istol,jc,kc,ispec) + fornow
          vtmp(istol,jc,kc) = vtmp(istol,jc,kc) + fornow
          erhs(istol,jc,kc) = erhs(istol,jc,kc) + fornow*utmp(istol,jc,kc)
          
        END DO
      END DO
    END IF
    IF(fyldfw)THEN
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = zero
          DO jc = jstap1,jstow
            
            fornow = fornow + acbcyl(jc-1)*store7(ic,jc,kc)*wd1y(ic,jc,kc)
            
          END DO
          rate(ic,jstal,kc,ispec) = rate(ic,jstal,kc,ispec) + fornow
          vtmp(ic,jstal,kc) = vtmp(ic,jstal,kc) + fornow
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc) + fornow*utmp(ic,jstal,kc)
          
        END DO
      END DO
    END IF
    IF(fyrdfw)THEN
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = zero
          DO jc = jstaw,jstom1
            
            fornow = fornow + acbcyr(jstol-jc)*store7(ic,jc,kc)*wd1y(ic,jc,kc)
            
          END DO
          rate(ic,jstol,kc,ispec) = rate(ic,jstol,kc,ispec) + fornow
          vtmp(ic,jstol,kc) = vtmp(ic,jstol,kc) + fornow
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc) + fornow*utmp(ic,jstol,kc)
          
        END DO
      END DO
    END IF
    IF(fzldfw)THEN
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = zero
          DO kc = kstap1,kstow
            
            fornow = fornow + acbczl(kc-1)*store7(ic,jc,kc)*wd1z(ic,jc,kc)
            
          END DO
          rate(ic,jc,kstal,ispec) = rate(ic,jc,kstal,ispec) + fornow
          vtmp(ic,jc,kstal) = vtmp(ic,jc,kstal) + fornow
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal) + fornow*utmp(ic,jc,kstal)
          
        END DO
      END DO
    END IF
    IF(fzrdfw)THEN
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = zero
          DO kc = kstaw,kstom1
            
            fornow = fornow + acbczr(kstol-kc)*store7(ic,jc,kc)*wd1z(ic,jc,kc)
            
          END DO
          rate(ic,jc,kstol,ispec) = rate(ic,jc,kstol,ispec) + fornow
          vtmp(ic,jc,kstol) = vtmp(ic,jc,kstol) + fornow
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol) + fornow*utmp(ic,jc,kstol)
          
        END DO
      END DO
    END IF
    
!         ====================================================================
    
!         Y-EQUATION: DIFFUSIVE TERMS
!         E-EQUATION: FURTHER HEAT FLUX TERMS
!         SECOND DERIVATIVE TERMS
!         SECOND DERIVATIVES OF LN(MIXTURE MOLAR MASS) ALREADY STORED
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif)CALL zeroxl(d_wd2x)
    IF(fxrdif)CALL zeroxr(d_wd2x)
!         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif)CALL zeroyl(d_wd2y)
    IF(fyrdif)CALL zeroyr(d_wd2y)
!         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif)CALL zerozl(d_wd2z)
    IF(fzrdif)CALL zerozr(d_wd2z)
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = (wd2x(ic,jc,kc) +  wd2y(ic,jc,kc)  &
              +  wd2z(ic,jc,kc))*store7(ic,jc,kc)
          
!               Y EQUATION
          rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec) + fornow
          
!               DIFFUSION CORRECTION VELOCITY DIVERGENCE
          vtmp(ic,jc,kc) = vtmp(ic,jc,kc) + fornow
          
        END DO
      END DO
    END DO
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb)CALL zeroxl(d_wd2x)
    IF(fxradb)CALL zeroxr(d_wd2x)
!         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb)CALL zeroyl(d_wd2y)
    IF(fyradb)CALL zeroyr(d_wd2y)
!         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb)CALL zerozl(d_wd2z)
    IF(fzradb)CALL zerozr(d_wd2z)
    
!         E EQUATION
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = (wd2x(ic,jc,kc) +  wd2y(ic,jc,kc)  &
              +  wd2z(ic,jc,kc))*store7(ic,jc,kc)
          
          erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*utmp(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
  END IF
!       MIXTURE MOLAR MASS TERMS
  
!       =======================================================================
  
!       PRESSURE DIFFUSION TERMS
  IF(flmixp)THEN
    
!         FIRST AND SECOND DERIVATIVES OF LN(PRESSURE) ALREADY STORED
    
    DO kc = kstab,kstob
      DO jc = jstab,jstob
        DO ic = istab,istob
          
          store7(ic,jc,kc) = difmix(ic,jc,kc)*yrhs(ic,jc,kc,ispec)  &
              *(one-wmolar(ispec)/wmomix(ic,jc,kc))
          
        END DO
      END DO
    END DO
    
!         DIFFUSION CORRECTION VELOCITY
!         FIRST DERIVATIVES OF LN(PRESSURE) ALREADY STORED
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_pd1x, 1, s3d_000, "real(dp)", OPS_READ))

    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_pd1y, 1, s3d_000, "real(dp)", OPS_READ))

    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_pd1z, 1, s3d_000, "real(dp)", OPS_READ))

!         Y EQUATION: DIFFUSIVE TERMS
!         E EQUATION: FURTHER HEAT FLUX TERMS
    
!         DIFFUSIVITY GRADIENT TERMS
    
!         DIFFUSIVITY GRADIENTS
    CALL dfbydx(d_store7,d_store1)
    CALL dfbydy(d_store7,d_store2)
    CALL dfbydz(d_store7,d_store3)
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif)CALL zeroxl(d_store1)
    IF(fxrdif)CALL zeroxr(d_store1)
!         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif)CALL zeroyl(d_store2)
    IF(fyrdif)CALL zeroyr(d_store2)
!         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif)CALL zerozl(d_store3)
    IF(fzrdif)CALL zerozr(d_store3)
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = store1(ic,jc,kc)*pd1x(ic,jc,kc)  &
              + store2(ic,jc,kc)*pd1y(ic,jc,kc) + store3(ic,jc,kc)*pd1z(ic,jc,kc)
          
!               Y EQUATION
          rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec) + fornow
          
!               DIFFUSION CORRECTION VELOCITY DIVERGENCE
          vtmp(ic,jc,kc) = vtmp(ic,jc,kc) + fornow
          
        END DO
      END DO
    END DO
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb)CALL zeroxl(d_store1)
    IF(fxradb)CALL zeroxr(d_store1)
!         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb)CALL zeroyl(d_store2)
    IF(fyradb)CALL zeroyr(d_store2)
!         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb)CALL zerozl(d_store3)
    IF(fzradb)CALL zerozr(d_store3)
    
!         E EQUATION
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = store1(ic,jc,kc)*pd1x(ic,jc,kc)  &
              + store2(ic,jc,kc)*pd1y(ic,jc,kc) + store3(ic,jc,kc)*pd1z(ic,jc,kc)
          
          erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*utmp(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
    
!         E-EQUATION: FURTHER HEAT FLUX TERMS
!         SPECIES ENTHALPY GRADIENT TERMS
    
!         SPECIES ENTHALPY GRADIENTS ALREADY IN STORE4,5,6
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif)CALL zeroxl(d_store4)
    IF(fxrdif)CALL zeroxr(d_store4)
!         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif)CALL zeroyl(d_store5)
    IF(fyrdif)CALL zeroyr(d_store5)
!         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif)CALL zerozl(d_store6)
    IF(fzrdif)CALL zerozr(d_store6)
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb)CALL zeroxl(d_store4)
    IF(fxradb)CALL zeroxr(d_store4)
!         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb)CALL zeroyl(d_store5)
    IF(fyradb)CALL zeroyr(d_store5)
!         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb)CALL zerozl(d_store6)
    IF(fzradb)CALL zerozr(d_store6)
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = store4(ic,jc,kc)*pd1x(ic,jc,kc)  &
              + store5(ic,jc,kc)*pd1y(ic,jc,kc) + store6(ic,jc,kc)*pd1z(ic,jc,kc)
          
          erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*store7(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
!         ====================================================================
    
!         Y-EQUATION: DIFFUSIVE TERMS
!         ---------------------------
!         WALL BC: PRESSURE TERMS
!         E-EQUATION: HEAT FLUX TERMS
!         WALL BC: ENTHALPY DIFFUSION TERMS
    IF(fxldfw)THEN
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = zero
          DO ic = istap1,istow
            
            fornow = fornow + acbcxl(ic-1)*store7(ic,jc,kc)*pd1x(ic,jc,kc)
            
          END DO
          rate(istal,jc,kc,ispec) = rate(istal,jc,kc,ispec) + fornow
          vtmp(istal,jc,kc) = vtmp(istal,jc,kc) + fornow
          erhs(istal,jc,kc) = erhs(istal,jc,kc) + fornow*utmp(istal,jc,kc)
          
        END DO
      END DO
    END IF
    IF(fxrdfw)THEN
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = zero
          DO ic = istaw,istom1
            
            fornow = fornow + acbcxr(istol-ic)*store7(ic,jc,kc)*pd1x(ic,jc,kc)
            
          END DO
          rate(istol,jc,kc,ispec) = rate(istol,jc,kc,ispec) + fornow
          vtmp(istol,jc,kc) = vtmp(istol,jc,kc) + fornow
          erhs(istol,jc,kc) = erhs(istol,jc,kc) + fornow*utmp(istol,jc,kc)
          
        END DO
      END DO
    END IF
    IF(fyldfw)THEN
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = zero
          DO jc = jstap1,jstow
            
            fornow = fornow + acbcyl(jc-1)*store7(ic,jc,kc)*pd1y(ic,jc,kc)
            
          END DO
          rate(ic,jstal,kc,ispec) = rate(ic,jstal,kc,ispec) + fornow
          vtmp(ic,jstal,kc) = vtmp(ic,jstal,kc) + fornow
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc) + fornow*utmp(ic,jstal,kc)
          
        END DO
      END DO
    END IF
    IF(fyrdfw)THEN
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = zero
          DO jc = jstaw,jstom1
            
            fornow = fornow + acbcyr(jstol-jc)*store7(ic,jc,kc)*pd1y(ic,jc,kc)
            
          END DO
          rate(ic,jstol,kc,ispec) = rate(ic,jstol,kc,ispec) + fornow
          vtmp(ic,jstol,kc) = vtmp(ic,jstol,kc) + fornow
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc) + fornow*utmp(ic,jstol,kc)
          
        END DO
      END DO
    END IF
    IF(fzldfw)THEN
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = zero
          DO kc = kstap1,kstow
            
            fornow = fornow + acbczl(kc-1)*store7(ic,jc,kc)*pd1z(ic,jc,kc)
            
          END DO
          rate(ic,jc,kstal,ispec) = rate(ic,jc,kstal,ispec) + fornow
          vtmp(ic,jc,kstal) = vtmp(ic,jc,kstal) + fornow
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal) + fornow*utmp(ic,jc,kstal)
          
        END DO
      END DO
    END IF
    IF(fzrdfw)THEN
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = zero
          DO kc = kstaw,kstom1
            
            fornow = fornow + acbczr(kstol-kc)*store7(ic,jc,kc)*pd1z(ic,jc,kc)
            
          END DO
          rate(ic,jc,kstol,ispec) = rate(ic,jc,kstol,ispec) + fornow
          vtmp(ic,jc,kstol) = vtmp(ic,jc,kstol) + fornow
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol) + fornow*utmp(ic,jc,kstol)
          
        END DO
      END DO
    END IF
    
!         ====================================================================
    
!         Y-EQUATION: DIFFUSIVE TERMS
!         E-EQUATION: FURTHER HEAT FLUX TERMS
!         SECOND DERIVATIVE TERMS
!         SECOND DERIVATIVES OF LN(PRESSURE) ALREADY STORED
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif)CALL zeroxl(d_pd2x)
    IF(fxrdif)CALL zeroxr(d_pd2x)
!         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif)CALL zeroyl(d_pd2y)
    IF(fyrdif)CALL zeroyr(d_pd2y)
!         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif)CALL zerozl(d_pd2z)
    IF(fzrdif)CALL zerozr(d_pd2z)
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = (pd2x(ic,jc,kc) +  pd2y(ic,jc,kc)  &
              +  pd2z(ic,jc,kc))*store7(ic,jc,kc)
          
!               Y EQUATION
          rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec) + fornow
          
!               DIFFUSION CORRECTION VELOCITY DIVERGENCE
          vtmp(ic,jc,kc) = vtmp(ic,jc,kc) + fornow
          
        END DO
      END DO
    END DO
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb)CALL zeroxl(d_pd2x)
    IF(fxradb)CALL zeroxr(d_pd2x)
!         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb)CALL zeroyl(d_pd2y)
    IF(fyradb)CALL zeroyr(d_pd2y)
!         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb)CALL zerozl(d_pd2z)
    IF(fzradb)CALL zerozr(d_pd2z)
    
!         E EQUATION
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = (pd2x(ic,jc,kc) +  pd2y(ic,jc,kc)  &
              +  pd2z(ic,jc,kc))*store7(ic,jc,kc)
          
          erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*utmp(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
!       SORET EFFECT (THERMAL DIFFUSION) TERMS
  IF(flmsor(ispec))THEN
    
!         FIRST AND SECOND DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED
    
    DO kc = kstab,kstob
      DO jc = jstab,jstob
        DO ic = istab,istob
          
          store7(ic,jc,kc) = difmix(ic,jc,kc)*yrhs(ic,jc,kc,ispec)  &
              *tdrmix(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
!         DIFFUSION CORRECTION VELOCITY
!         FIRST DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_ucor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_td1x, 1, s3d_000, "real(dp)", OPS_READ))

    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_td1y, 1, s3d_000, "real(dp)", OPS_READ))

    call ops_par_loop(compute_kernel_AequalsAplusBmulC, "A=A+B*C", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wcor, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_store7, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_td1z, 1, s3d_000, "real(dp)", OPS_READ))

!         Y EQUATION: DIFFUSIVE TERMS
!         E EQUATION: FURTHER HEAT FLUX TERMS
    
!         DIFFUSIVITY GRADIENT TERMS
    
!         DIFFUSIVITY GRADIENTS
    CALL dfbydx(d_store7,d_store1)
    CALL dfbydy(d_store7,d_store2)
    CALL dfbydz(d_store7,d_store3)
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif)CALL zeroxl(d_store1)
    IF(fxrdif)CALL zeroxr(d_store1)
!         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif)CALL zeroyl(d_store2)
    IF(fyrdif)CALL zeroyr(d_store2)
!         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif)CALL zerozl(d_store3)
    IF(fzrdif)CALL zerozr(d_store3)
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = store1(ic,jc,kc)*td1x(ic,jc,kc)  &
              + store2(ic,jc,kc)*td1y(ic,jc,kc) + store3(ic,jc,kc)*td1z(ic,jc,kc)
          
!               Y EQUATION
          rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec) + fornow
          
!               DIFFUSION CORRECTION VELOCITY DIVERGENCE
          vtmp(ic,jc,kc) = vtmp(ic,jc,kc) + fornow
          
        END DO
      END DO
    END DO
    
!         SUBTRACT DUFOUR EFFECT TERMS TO RESTORE SPECIES ENTHALPY
!         RSC 08-JUN-2015 BUG FIX
    IF(flmduf(ispec))THEN
      DO kc = kstab,kstob
        DO jc = jstab,jstob
          DO ic = istab,istob
            
            utmp(ic,jc,kc) = utmp(ic,jc,kc)  &
                - rgspec(ispec)*trun(ic,jc,kc)*tdrmix(ic,jc,kc)
            
          END DO
        END DO
      END DO
    END IF
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb)CALL zeroxl(d_store1)
    IF(fxradb)CALL zeroxr(d_store1)
!         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb)CALL zeroyl(d_store2)
    IF(fyradb)CALL zeroyr(d_store2)
!         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb)CALL zerozl(d_store3)
    IF(fzradb)CALL zerozr(d_store3)
    
!         E EQUATION
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = store1(ic,jc,kc)*td1x(ic,jc,kc)  &
              + store2(ic,jc,kc)*td1y(ic,jc,kc) + store3(ic,jc,kc)*td1z(ic,jc,kc)
          
          erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*utmp(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
    
!         E-EQUATION: FURTHER HEAT FLUX TERMS
!         SPECIES ENTHALPY GRADIENT TERMS
    
!         EVALUATE SPECIES ENTHALPY GRADIENTS USING STORE4,5,6
    CALL dfbydx(d_utmp,d_store4)
    CALL dfbydy(d_utmp,d_store5)
    CALL dfbydz(d_utmp,d_store6)
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif)CALL zeroxl(d_store4)
    IF(fxrdif)CALL zeroxr(d_store4)
!         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif)CALL zeroyl(d_store5)
    IF(fyrdif)CALL zeroyr(d_store5)
!         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif)CALL zerozl(d_store6)
    IF(fzrdif)CALL zerozr(d_store6)
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb)CALL zeroxl(d_store4)
    IF(fxradb)CALL zeroxr(d_store4)
!         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb)CALL zeroyl(d_store5)
    IF(fyradb)CALL zeroyr(d_store5)
!         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb)CALL zerozl(d_store6)
    IF(fzradb)CALL zerozr(d_store6)
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = store4(ic,jc,kc)*td1x(ic,jc,kc)  &
              + store5(ic,jc,kc)*td1y(ic,jc,kc) + store6(ic,jc,kc)*td1z(ic,jc,kc)
          
          erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*store7(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
!         ====================================================================
    
!         Y-EQUATION: DIFFUSIVE TERMS
!         ---------------------------
!         WALL BC: SORET EFFECT TERMS
!         E-EQUATION: HEAT FLUX TERMS
!         WALL BC: ENTHALPY DIFFUSION TERMS
    IF(fxldfw)THEN
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = zero
          DO ic = istap1,istow
            
            fornow = fornow + acbcxl(ic-1)*store7(ic,jc,kc)*td1x(ic,jc,kc)
            
          END DO
          rate(istal,jc,kc,ispec) = rate(istal,jc,kc,ispec) + fornow
          vtmp(istal,jc,kc) = vtmp(istal,jc,kc) + fornow
          erhs(istal,jc,kc) = erhs(istal,jc,kc) + fornow*(utmp(istal,jc,kc)  &
              + rgspec(ispec)*trun(istal,jc,kc)*tdrmix(istal,jc,kc))
          
        END DO
      END DO
    END IF
    IF(fxrdfw)THEN
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = zero
          DO ic = istaw,istom1
            
            fornow = fornow + acbcxr(istol-ic)*store7(ic,jc,kc)*td1x(ic,jc,kc)
            
          END DO
          rate(istol,jc,kc,ispec) = rate(istol,jc,kc,ispec) + fornow
          vtmp(istol,jc,kc) = vtmp(istol,jc,kc) + fornow
          erhs(istol,jc,kc) = erhs(istol,jc,kc) + fornow*(utmp(istol,jc,kc)  &
              + rgspec(ispec)*trun(istol,jc,kc)*tdrmix(istol,jc,kc))
          
        END DO
      END DO
    END IF
    IF(fyldfw)THEN
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = zero
          DO jc = jstap1,jstow
            
            fornow = fornow + acbcyl(jc-1)*store7(ic,jc,kc)*td1y(ic,jc,kc)
            
          END DO
          rate(ic,jstal,kc,ispec) = rate(ic,jstal,kc,ispec) + fornow
          vtmp(ic,jstal,kc) = vtmp(ic,jstal,kc) + fornow
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc) + fornow*(utmp(ic,jstal,kc)  &
              + rgspec(ispec)*trun(ic,jstal,kc)*tdrmix(ic,jstal,kc))
          
        END DO
      END DO
    END IF
    IF(fyrdfw)THEN
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = zero
          DO jc = jstaw,jstom1
            
            fornow = fornow + acbcyr(jstol-jc)*store7(ic,jc,kc)*td1y(ic,jc,kc)
            
          END DO
          rate(ic,jstol,kc,ispec) = rate(ic,jstol,kc,ispec) + fornow
          vtmp(ic,jstol,kc) = vtmp(ic,jstol,kc) + fornow
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc) + fornow*(utmp(ic,jstol,kc)  &
              + rgspec(ispec)*trun(ic,jstol,kc)*tdrmix(ic,jstol,kc))
          
        END DO
      END DO
    END IF
    IF(fzldfw)THEN
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = zero
          DO kc = kstap1,kstow
            
            fornow = fornow + acbczl(kc-1)*store7(ic,jc,kc)*td1z(ic,jc,kc)
            
          END DO
          rate(ic,jc,kstal,ispec) = rate(ic,jc,kstal,ispec) + fornow
          vtmp(ic,jc,kstal) = vtmp(ic,jc,kstal) + fornow
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal) + fornow*(utmp(ic,jc,kstal)  &
              + rgspec(ispec)*trun(ic,jc,kstal)*tdrmix(ic,jc,kstal))
          
        END DO
      END DO
    END IF
    IF(fzrdfw)THEN
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = zero
          DO kc = kstaw,kstom1
            
            fornow = fornow + acbczr(kstol-kc)*store7(ic,jc,kc)*td1z(ic,jc,kc)
            
          END DO
          rate(ic,jc,kstol,ispec) = rate(ic,jc,kstol,ispec) + fornow
          vtmp(ic,jc,kstol) = vtmp(ic,jc,kstol) + fornow
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol) + fornow*(utmp(ic,jc,kstol)  &
              + rgspec(ispec)*trun(ic,jc,kstol)*tdrmix(ic,jc,kstol))
          
        END DO
      END DO
    END IF
    
!         E-EQUATION: HEAT FLUX TERMS
!         WALL BC: SORET AND DUFOUR TERMS
    IF(flmduf(ispec))THEN
!           E-EQUATION: HEAT FLUX TERMS
!           WALL BC: ADIABATIC WALL
      IF(fxlcnw)THEN
        DO kc = kstal,kstol
          DO jc = jstal,jstol
            
            fornow = zero
            DO ic = istap1,istow
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              fornow = fornow  &
                  + acbcxl(ic-1)*combo1*store7(ic,jc,kc)*td1x(ic,jc,kc)
              
            END DO
            erhs(istal,jc,kc) = erhs(istal,jc,kc) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fxrcnw)THEN
        DO kc = kstal,kstol
          DO jc = jstal,jstol
            
            fornow = zero
            DO ic = istaw,istom1
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              fornow = fornow  &
                  + acbcxr(istol-ic)*combo1*store7(ic,jc,kc)*td1x(ic,jc,kc)
              
            END DO
            erhs(istol,jc,kc) = erhs(istol,jc,kc) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fylcnw)THEN
        DO kc = kstal,kstol
          DO ic = istal,istol
            
            fornow = zero
            DO jc = jstap1,jstow
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              fornow = fornow  &
                  + acbcyl(jc-1)*combo1*store7(ic,jc,kc)*td1y(ic,jc,kc)
              
            END DO
            erhs(ic,jstal,kc) = erhs(ic,jstal,kc) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fyrcnw)THEN
        DO kc = kstal,kstol
          DO ic = istal,istol
            
            fornow = zero
            DO jc = jstaw,jstom1
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              fornow = fornow  &
                  + acbcyr(jstol-jc)*combo1*store7(ic,jc,kc)*td1y(ic,jc,kc)
              
            END DO
            erhs(ic,jstol,kc) = erhs(ic,jstol,kc) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fzlcnw)THEN
        DO jc = jstal,jstol
          DO ic = istal,istol
            
            fornow = zero
            DO kc = kstap1,kstow
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              fornow = fornow  &
                  + acbczl(kc-1)*combo1*store7(ic,jc,kc)*td1z(ic,jc,kc)
              
            END DO
            erhs(ic,jc,kstal) = erhs(ic,jc,kstal) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fzrcnw)THEN
        DO jc = jstal,jstol
          DO ic = istal,istol
            
            fornow = zero
            DO kc = kstaw,kstom1
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              fornow = fornow  &
                  + acbczr(kstol-kc)*combo1*store7(ic,jc,kc)*td1z(ic,jc,kc)
              
            END DO
            erhs(ic,jc,kstol) = erhs(ic,jc,kstol) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      
!           E-EQUATION: HEAT FLUX TERMS
!           WALL BC: ISOTHERMAL WALL
      IF(fxladw)THEN
        DO kc = kstal,kstol
          DO jc = jstal,jstol
            
            combo2 = trun(istal,jc,kc)*tdrmix(istal,jc,kc)
            combo2 = combo2*store7(istal,jc,kc)*td1x(istal,jc,kc)
            fornow = zero
            DO ic = istap1,istow
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              combo1 = combo1*store7(ic,jc,kc)*td1x(ic,jc,kc)
              fornow = fornow + acbcxl(ic-1)*rgspec(ispec)*(combo1-combo2)
              
            END DO
            erhs(istal,jc,kc) = erhs(istal,jc,kc) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fxradw)THEN
        DO kc = kstal,kstol
          DO jc = jstal,jstol
            
            combo2 = trun(istol,jc,kc)*tdrmix(istol,jc,kc)
            combo2 = combo2*store7(istol,jc,kc)*td1x(istol,jc,kc)
            fornow = zero
            DO ic = istaw,istom1
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              combo1 = combo1*store7(ic,jc,kc)*td1x(ic,jc,kc)
              fornow = fornow + acbcxr(istol-ic)*rgspec(ispec)*(combo2-combo1)
              
            END DO
            erhs(istol,jc,kc) = erhs(istol,jc,kc) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fyladw)THEN
        DO kc = kstal,kstol
          DO ic = istal,istol
            
            combo2 = trun(ic,jstal,kc)*tdrmix(ic,jstal,kc)
            combo2 = combo2*store7(ic,jstal,kc)*td1y(ic,jstal,kc)
            fornow = zero
            DO jc = jstap1,jstow
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              combo1 = combo1*store7(ic,jc,kc)*td1y(ic,jc,kc)
              fornow = fornow + acbcyl(jc-1)*rgspec(ispec)*(combo1-combo2)
              
            END DO
            erhs(ic,jstal,kc) = erhs(ic,jstal,kc) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fyradw)THEN
        DO kc = kstal,kstol
          DO ic = istal,istol
            
            combo2 = trun(ic,jstol,kc)*tdrmix(ic,jstol,kc)
            combo2 = combo2*store7(ic,jstol,kc)*td1y(ic,jstol,kc)
            fornow = zero
            DO jc = jstaw,jstom1
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              combo1 = combo1*store7(ic,jc,kc)*td1y(ic,jc,kc)
              fornow = fornow + acbcyr(jstol-jc)*rgspec(ispec)*(combo2-combo1)
              
            END DO
            erhs(ic,jstol,kc) = erhs(ic,jstol,kc) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fzladw)THEN
        DO jc = jstal,jstol
          DO ic = istal,istol
            
            combo2 = trun(ic,jc,kstal)*tdrmix(ic,jc,kstal)
            combo2 = combo2*store7(ic,jc,kstal)*td1z(ic,jc,kstal)
            fornow = zero
            DO kc = kstap1,kstow
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              combo1 = combo1*store7(ic,jc,kc)*td1z(ic,jc,kc)
              fornow = fornow + acbczl(kc-1)*rgspec(ispec)*(combo1-combo2)
              
            END DO
            erhs(ic,jc,kstal) = erhs(ic,jc,kstal) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
      IF(fzradw)THEN
        DO jc = jstal,jstol
          DO ic = istal,istol
            
            combo2 = trun(ic,jc,kstol)*tdrmix(ic,jc,kstol)
            combo2 = combo2*store7(ic,jc,kstol)*td1z(ic,jc,kstol)
            fornow = zero
            DO kc = kstaw,kstom1
              
              combo1 = trun(ic,jc,kc)*tdrmix(ic,jc,kc)
              combo1 = combo1*store7(ic,jc,kc)*td1z(ic,jc,kc)
              fornow = fornow + acbczr(kstol-kc)*rgspec(ispec)*(combo2-combo1)
              
            END DO
            erhs(ic,jc,kstol) = erhs(ic,jc,kstol) + rgspec(ispec)*fornow
            
          END DO
        END DO
      END IF
    END IF
    
!         ====================================================================
    
!         Y-EQUATION: DIFFUSIVE TERMS
!         E-EQUATION: FURTHER HEAT FLUX TERMS
!         SECOND DERIVATIVE TERMS
!         SECOND DERIVATIVES OF LN(TEMPERATURE) ALREADY STORED
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fxldif)CALL zeroxl(d_td2x)
    IF(fxrdif)CALL zeroxr(d_td2x)
!         BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fyldif)CALL zeroyl(d_td2y)
    IF(fyrdif)CALL zeroyr(d_td2y)
!         BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
    IF(fzldif)CALL zerozl(d_td2z)
    IF(fzrdif)CALL zerozr(d_td2z)
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = (td2x(ic,jc,kc) +  td2y(ic,jc,kc)  &
              +  td2z(ic,jc,kc))*store7(ic,jc,kc)
          
!               Y EQUATION
          rate(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec) + fornow
          
!               DIFFUSION CORRECTION VELOCITY DIVERGENCE
          vtmp(ic,jc,kc) = vtmp(ic,jc,kc) + fornow
          
        END DO
      END DO
    END DO
    
!         BOUNDARY CONDITIONS
!         BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fxladb)CALL zeroxl(d_td2x)
    IF(fxradb)CALL zeroxr(d_td2x)
!         BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fyladb)CALL zeroyl(d_td2y)
    IF(fyradb)CALL zeroyr(d_td2y)
!         BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
    IF(fzladb)CALL zerozl(d_td2z)
    IF(fzradb)CALL zerozr(d_td2z)
    
!         E EQUATION
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = (td2x(ic,jc,kc) +  td2y(ic,jc,kc)  &
              +  td2z(ic,jc,kc))*store7(ic,jc,kc)
          
          erhs(ic,jc,kc) = erhs(ic,jc,kc) + fornow*utmp(ic,jc,kc)
          
        END DO
      END DO
    END DO
    
  END IF
  
!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!                                                              WTMP = MIXTURE H
!       =======================================================================
  
!       ----------------------------------------------------------------
!       E-EQUATION: DIFFUSION CORRECTION VELOCITY TERMS EVALUATED BELOW
!       Y-EQUATION: DIFFUSION CORRECTION VELOCITY TERMS EVALUATED BELOW
!       ----------------------------------------------------------------
  
END DO
!     RSC 08-AUG-2012 EVALUATE ALL SPECIES
!     END OF RUN THROUGH ALL SPECIES

!     =========================================================================

!     EVALUATE DIFFUSION CORRECTION VELOCITY TERMS
!     --------------------------------------------

!     E-EQUATION: FURTHER HEAT FLUX TERMS
!     -----------------------------------
!     MIXTURE ENTHALPY GRADIENTS
CALL dfbydx(d_wtmp,d_store1)
CALL dfbydy(d_wtmp,d_store2)
CALL dfbydz(d_wtmp,d_store3)

!     BOUNDARY CONDITIONS
!     BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
IF(fxldif)CALL zeroxl(d_store1)
IF(fxrdif)CALL zeroxr(d_store1)
!     BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
IF(fyldif)CALL zeroyl(d_store2)
IF(fyrdif)CALL zeroyr(d_store2)
!     BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
IF(fzldif)CALL zerozl(d_store3)
IF(fzrdif)CALL zerozr(d_store3)

!     BOUNDARY CONDITIONS
!     BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fxladb)CALL zeroxl(d_store1)
IF(fxradb)CALL zeroxr(d_store1)
!     BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fyladb)CALL zeroyl(d_store2)
IF(fyradb)CALL zeroyr(d_store2)
!     BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fzladb)CALL zerozl(d_store3)
IF(fzradb)CALL zerozr(d_store3)

!     TRANSFER DIV CORR VEL TO TEMPORARY STORE
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store4, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_vtmp, 1, s3d_000, "real(dp)", OPS_READ))

!     BOUNDARY CONDITIONS
!     BC IN X: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fxladb)CALL zeroxl(d_store4)
IF(fxradb)CALL zeroxr(d_store4)
!     BC IN Y: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fyladb)CALL zeroyl(d_store4)
IF(fyradb)CALL zeroyr(d_store4)
!     BC IN Z: DIFFUSIVE TERMS (HEAT FLUX) ZERO ON END POINTS
IF(fzladb)CALL zerozl(d_store4)
IF(fzradb)CALL zerozr(d_store4)

!     DIV RHO VCORR HMIX
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      erhs(ic,jc,kc) = erhs(ic,jc,kc) - wtmp(ic,jc,kc)*store4(ic,jc,kc)  &
          - store1(ic,jc,kc)*ucor(ic,jc,kc) - store2(ic,jc,kc)*vcor(ic,jc,kc)  &
          - store3(ic,jc,kc)*wcor(ic,jc,kc)
      
    END DO
  END DO
END DO
!                                                         RATE = Y SOURCE TERMS
!                                                           VTMP = DIV CORR VEL
!     =========================================================================

!     MIXTURE AVERAGED TRANSPORT
!     EVALUATE THE VISCOSITY

!     RSC 17-APR-2013
!     TRANSP CONTAINS LN(T)
!     STORE VISCOSITY IN DIFMIX FOR NOW
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
            combo2 = combo2 + yrhs(ic,jc,kc,jspec)*ovwmol(jspec)*fornow
          END DO
          fornow = ctrans(ispec)/combo2
          combo1 = combo1 + yrhs(ic,jc,kc,ispec)*ovwmol(ispec)*fornow
          
        END DO
        difmix(ic,jc,kc) = combo1
        
      END DO
    END DO
  END DO
  
END IF

!     =========================================================================

!     RUN THROUGH ALL SPECIES
!     -----------------------
!     RSC 08-AUG-2012 EVALUATE ALL SPECIES
!     RSC 08-JUN-2015 REMOVE Nth SPECIES TREATMENT
DO ispec = 1,nspec
  
!       Y-EQUATION: DIFFUSIVE TERMS
!       ---------------------------
!       RECOMPUTE SPECIES MASS FRACTION GRADIENTS
  DO kc = kstab,kstob
    DO jc = jstab,jstob
      DO ic = istab,istob
        
        store7(ic,jc,kc) = yrhs(ic,jc,kc,ispec)
        
      END DO
    END DO
  END DO
  CALL dfbydx(d_store7,d_store1)
  CALL dfbydy(d_store7,d_store2)
  CALL dfbydz(d_store7,d_store3)
  
!       BOUNDARY CONDITIONS
!       BC IN X: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fxldif)CALL zeroxl(d_store1)
  IF(fxrdif)CALL zeroxr(d_store1)
!       BC IN Y: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fyldif)CALL zeroyl(d_store2)
  IF(fyrdif)CALL zeroyr(d_store2)
!       BC IN Z: DIFFUSIVE TERMS (MASS FLUX) ZERO ON END POINTS
  IF(fzldif)CALL zerozl(d_store3)
  IF(fzrdif)CALL zerozr(d_store3)
  
!       DIV RHO VCORR Y
!       STORE Y SOURCE TERMS IN YRHS
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        yrhs(ic,jc,kc,ispec) = rate(ic,jc,kc,ispec)  &
            - yrhs(ic,jc,kc,ispec)*vtmp(ic,jc,kc)  &
            - store1(ic,jc,kc)*ucor(ic,jc,kc) - store2(ic,jc,kc)*vcor(ic,jc,kc)  &
            - store3(ic,jc,kc)*wcor(ic,jc,kc)
        
      END DO
    END DO
  END DO
  
END DO

!     RSC 08-AUG-2012 EVALUATE ALL SPECIES
!     END OF RUN THROUGH ALL SPECIES
!                                                              ALL STORES CLEAR
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ------------------------------------------------
!     Y-EQUATION: SOURCE TERMS COMPLETE
!     ------------------------------------------------
!     E-EQUATION: PRESSURE-WORK AND VISCOUS WORK TERMS
!                 EVALUATED IN SUBROUTINE RHSVEL
!     ------------------------------------------------

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     COLLECT DENSITY AND ITS GRADIENTS FOR BCs
!     -----------------------------------------

!     X-DIRECTION: DRHODX
IF(fxlcnv.OR.fxrcnv)THEN
  
  CALL dfbydx(d_drhs,d_store1)
  
  IF(fxlcnv)THEN
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        strdxl(jc,kc) = drhs(istal,jc,kc)
        bcl2xl(jc,kc) = store1(istal,jc,kc)
        
      END DO
    END DO
  END IF
  IF(fxrcnv)THEN
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        strdxr(jc,kc) = drhs(istol,jc,kc)
        bcl2xr(jc,kc) = store1(istol,jc,kc)
        
      END DO
    END DO
  END IF
  
END IF

!     Y-DIRECTION: DRHODY
IF(fylcnv.OR.fyrcnv)THEN
  
  CALL dfbydy(d_drhs,d_store2)
  
  IF(fylcnv)THEN
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        strdyl(ic,kc) = drhs(ic,jstal,kc)
        bcl2yl(ic,kc) = store2(ic,jstal,kc)
        
      END DO
    END DO
  END IF
  IF(fyrcnv)THEN
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        strdyr(ic,kc) = drhs(ic,jstol,kc)
        bcl2yr(ic,kc) = store2(ic,jstol,kc)
        
      END DO
    END DO
  END IF
  
END IF

!     Z-DIRECTION: DRHODZ
IF(fzlcnv.OR.fzrcnv)THEN
  
  CALL dfbydz(d_drhs,d_store3)
  
  IF(fzlcnv)THEN
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        strdzl(ic,jc) = drhs(ic,jc,kstal)
        bcl2zl(ic,jc) = store3(ic,jc,kstal)
        
      END DO
    END DO
  END IF
  IF(fzrcnv)THEN
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        strdzr(ic,jc) = drhs(ic,jc,kstol)
        bcl2zr(ic,jc) = store3(ic,jc,kstol)
        
      END DO
    END DO
  END IF
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE rhscal
