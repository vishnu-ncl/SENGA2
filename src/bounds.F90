SUBROUTINE bounds

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
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

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!     -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
real(kind=dp) :: fornow
INTEGER :: ic,jc,kc
    integer :: ispec
    integer :: rangexyz(6)

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
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_reduced_enthalpy_xdir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strhxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                    ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
        call ops_par_loop(bounds_kernel_reduced_energy_xdir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strexl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_ovgmxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))
                        
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      fornow = strgxl(1,jc,kc)*gam1xl(1,jc,kc)*strtxl(1,jc,kc)
      acouxl(1,jc,kc) = SQRT(fornow)
      ova2xl(1,jc,kc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
        IF(nsbcxl == nsbco1) THEN
    
!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY
    
!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))
    
            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratexl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SPECIFY L5X AS REQUIRED
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_outflowBC1_computeL_xl, "SPECIFY L5X AS REQUIRED", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_gbl(cobcxl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(pinfxl, 1, "real(dp)", OPS_READ))
    
!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)
        
        urhs(istal,jc,kc) = urhs(istal,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)*(struxl(1,jc,kc)+acouxl(1,jc,kc))
        
        vrhs(istal,jc,kc) = vrhs(istal,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)*strvxl(1,jc,kc)
        
        wrhs(istal,jc,kc) = wrhs(istal,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)*strwxl(1,jc,kc)
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
            - bcl5xl(1,jc,kc)*(ova2xl(1,jc,kc)*strexl(1,jc,kc)  &
            + struxl(1,jc,kc)/acouxl(1,jc,kc) + ovgmxl(1,jc,kc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(ispec,istal,jc,kc) = yrhs(ispec,istal,jc,kc)  &
              - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)*stryxl(ispec,1,jc,kc)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
        IF(nsbcxl == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratexl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L2X-L5X
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_inflowBC1_computeL_xl, "L2X to L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_gbl(cobcxl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(pinfxl, 1, "real(dp)", OPS_READ))
            
!           LYX
            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_inflowBC1_LYX_xl, "LYX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_ratexl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO
    
!           ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)
        
        urhs(istal,jc,kc) = urhs(istal,jc,kc) - bcl2xl(1,jc,kc)*struxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)*(struxl(1,jc,kc)+acouxl(1,jc,kc))
        
        vrhs(istal,jc,kc) = vrhs(istal,jc,kc) - bcl2xl(1,jc,kc)*strvxl(1,jc,kc)  &
            - bcl3xl(1,jc,kc)*strdxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)*strvxl(1,jc,kc)
        
        wrhs(istal,jc,kc) = wrhs(istal,jc,kc) - bcl2xl(1,jc,kc)*strwxl(1,jc,kc)  &
            - bcl4xl(1,jc,kc)*strdxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)*strwxl(1,jc,kc)
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc) - bcl2xl(1,jc,kc)*strexl(1,jc,kc)  &
            - bcl3xl(1,jc,kc)*strdxl(1,jc,kc)*strvxl(1,jc,kc)  &
            - bcl4xl(1,jc,kc)*strdxl(1,jc,kc)*strwxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*(ova2xl(1,jc,kc)*strexl(1,jc,kc)  &
            + struxl(1,jc,kc)/acouxl(1,jc,kc) + ovgmxl(1,jc,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = bclyxl(ispec,1,jc,kc)*strdxl(1,jc,kc)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc) - fornow*strhxl(ispec,1,jc,kc)
          
          yrhs(ispec,istal,jc,kc) = yrhs(ispec,istal,jc,kc)  &
              - (bcl2xl(1,jc,kc)+bcl5xl(1,jc,kc)*ova2xl(1,jc,kc))*stryxl(ispec,1,jc,kc)  &
              - fornow
          
        END DO
      END DO
      
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
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqB_xdir, "A_yz = A_yz + B_mulditim_yz*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_dydtxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))
                call ops_par_loop(bounds_kernel_eqA_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratexl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz/B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqC_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X,L2X,L5X
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_inflowBC2_computeL_xl, "L1X L2X L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dtdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sydtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)
        
      END DO
    END DO
    
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
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_inflowBC3_computeL_xl, "L1X to L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dddtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc) - bcl2xl(1,jc,kc)*strexl(1,jc,kc)  &
            - bcl3xl(1,jc,kc)*strdxl(1,jc,kc)*strvxl(1,jc,kc)  &
            - bcl4xl(1,jc,kc)*strdxl(1,jc,kc)*strwxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*(ova2xl(1,jc,kc)*strexl(1,jc,kc)  &
            + struxl(1,jc,kc)/acouxl(1,jc,kc) + ovgmxl(1,jc,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          bclyxl(ispec,1,jc,kc) = ratexl(ispec,1,jc,kc)/strdxl(1,jc,kc)  &
              - dydtxl(ispec,1,jc,kc) - struxl(1,jc,kc)*bclyxl(ispec,1,jc,kc)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
              - bclyxl(ispec,1,jc,kc)*strdxl(1,jc,kc)*strhxl(ispec,1,jc,kc)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       WALL BOUNDARY CONDITIONS
!       ------------------------
  
        IF(nsbcxl == nsbcw1)THEN
    
!           WALL BOUNDARY CONDITION No 1
!           NO-SLIP WALL - ADIABATIC
    
!           ALL VELOCITY COMPONENTS IMPOSED
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT
    
!           SPECIFY L's AS REQUIRED
!           L1X,L3X-L5X
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_wallBC1_computeL_xl, "L1X and L3X to L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), & 
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
            - bcl3xl(1,jc,kc)*strdxl(1,jc,kc)*strvxl(1,jc,kc)  &
            - bcl4xl(1,jc,kc)*strdxl(1,jc,kc)*strwxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*(ova2xl(1,jc,kc)*strexl(1,jc,kc)  &
            + struxl(1,jc,kc)/acouxl(1,jc,kc) + ovgmxl(1,jc,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(ispec,istal,jc,kc) = yrhs(ispec,istal,jc,kc)  &
              - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)*stryxl(ispec,1,jc,kc)
          
        END DO
      END DO
      
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
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratexl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X-L5X
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_wallBC2_computeL_xl, "L1X to L5X", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_gam1xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dtdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))           

!           LYX
            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_wallBC2_LYX_xl, "LYX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_bcl2xl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_ratexl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strrxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ),
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(1,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(ispec,istal,jc,kc) = yrhs(ispec,istal,jc,kc)  &
              - (bcl2xl(1,jc,kc)+bcl5xl(1,jc,kc)*ova2xl(1,jc,kc))*stryxl(ispec,1,jc,kc)
          
        END DO
      END DO
      
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
    IF(fxrcnv)THEN
  
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
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_reduced_enthalpy_xdir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                        ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))    

        END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
        call ops_par_loop(bounds_kernel_reduced_energy_xdir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strexr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_ovgmxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      fornow = strgxr(1,jc,kc)*gam1xr(1,jc,kc)*strtxr(1,jc,kc)
      acouxr(1,jc,kc) = SQRT(fornow)
      ova2xr(1,jc,kc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
        IF(nsbcxr == nsbco1) THEN

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratexr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SPECIFY L1X AS REQUIRED
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L1X
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)-strdxr(1,jc,kc)*acouxr(1,jc,kc)*bcl1xr(1,jc,kc))
        
!             SUBTRACT FROM NEW VALUE OF L1X
        bcl1xr(1,jc,kc)= half*sorpxr(1,jc,kc)  &
            + cobcxr*acouxr(1,jc,kc)*(strpxr(1,jc,kc)-pinfxr) - bcl1xr(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)
        
        urhs(istol,jc,kc) = urhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)*(struxr(1,jc,kc)-acouxr(1,jc,kc))
        
        vrhs(istol,jc,kc) = vrhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)*strvxr(1,jc,kc)
        
        wrhs(istol,jc,kc) = wrhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)*strwxr(1,jc,kc)
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*(ova2xr(1,jc,kc)*strexr(1,jc,kc)  &
            - struxr(1,jc,kc)/acouxr(1,jc,kc) + ovgmxr(1,jc,kc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(ispec,istol,jc,kc) = yrhs(ispec,istol,jc,kc)  &
              - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)*stryxr(ispec,1,jc,kc)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
        IF(nsbcxr == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratexr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!         L1X-L4X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxr(1,jc,kc)*acouxr(1,jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl2xr(1,jc,kc) = struxr(1,jc,kc)  &
            *(bcl2xr(1,jc,kc)-bcl5xr(1,jc,kc)*ova2xr(1,jc,kc))
        bcl3xr(1,jc,kc) = struxr(1,jc,kc)*bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = struxr(1,jc,kc)*bcl4xr(1,jc,kc)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
!             L5X UNCHANGED
        bcl1xr(1,jc,kc) = half*sorpxr(1,jc,kc)  &
            + cobcxr*acouxr(1,jc,kc)*(strpxr(1,jc,kc)-pinfxr) - bcl1xr(1,jc,kc)
        bcl2xr(1,jc,kc) = -bcl2xr(1,jc,kc)
        bcl3xr(1,jc,kc) = -bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = -bcl4xr(1,jc,kc)
        
      END DO
    END DO
    
!         LYX
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
!               OLD VALUE OF L's
          bclyxr(ispec,1,jc,kc) = struxr(1,jc,kc)*bclyxr(ispec,1,jc,kc)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
          bclyxr(ispec,1,jc,kc) = ratexr(ispec,1,jc,kc)/strdxr(1,jc,kc)  &
              - bclyxr(ispec,1,jc,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)  &
            - bcl2xr(1,jc,kc)
        
        urhs(istol,jc,kc) = urhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)*(struxr(1,jc,kc)-acouxr(1,jc,kc))  &
            - bcl2xr(1,jc,kc)*struxr(1,jc,kc)
        
        vrhs(istol,jc,kc) = vrhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)*strvxr(1,jc,kc)  &
            - bcl2xr(1,jc,kc)*strvxr(1,jc,kc) - bcl3xr(1,jc,kc)*strdxr(1,jc,kc)
        
        wrhs(istol,jc,kc) = wrhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)*strwxr(1,jc,kc)  &
            - bcl2xr(1,jc,kc)*strwxr(1,jc,kc) - bcl4xr(1,jc,kc)*strdxr(1,jc,kc)
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*(ova2xr(1,jc,kc)*strexr(1,jc,kc)  &
            - struxr(1,jc,kc)/acouxr(1,jc,kc) + ovgmxr(1,jc,kc))  &
            - bcl2xr(1,jc,kc)*strexr(1,jc,kc)  &
            - bcl3xr(1,jc,kc)*strdxr(1,jc,kc)*strvxr(1,jc,kc)  &
            - bcl4xr(1,jc,kc)*strdxr(1,jc,kc)*strwxr(1,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = bclyxr(ispec,1,jc,kc)*strdxr(1,jc,kc)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc) - fornow*strhxr(ispec,1,jc,kc)
          
          yrhs(ispec,istol,jc,kc) = yrhs(ispec,istol,jc,kc)  &
              - (bcl2xr(1,jc,kc)+bcl1xr(1,jc,kc)*ova2xr(1,jc,kc))*stryxr(ispec,1,jc,kc)  &
              - fornow
          
        END DO
      END DO
      
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
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))
    
            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqB_xdir, "A_yz = A_yz + B_mulditim_yz*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_dydtxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))
                call ops_par_loop(bounds_kernel_eqA_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratexr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqD_xdir, "A_yz = A_yz/B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strrxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqC_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X,L2X,L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxr(1,jc,kc)*acouxr(1,jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl2xr(1,jc,kc) = struxr(1,jc,kc)  &
            *(bcl2xr(1,jc,kc)-bcl5xr(1,jc,kc)*ova2xr(1,jc,kc))
        bcl5xr(1,jc,kc) = half*(struxr(1,jc,kc)+acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5X UNCHANGED
        bcl1xr(1,jc,kc) = bcl5xr(1,jc,kc)  &
            + strdxr(1,jc,kc)*acouxr(1,jc,kc)*dudtxr(1,jc,kc) - bcl1xr(1,jc,kc)
        bcl2xr(1,jc,kc) = gam1xr(1,jc,kc)*ova2xr(1,jc,kc)  &
            *(bcl1xr(1,jc,kc)+bcl5xr(1,jc,kc))  &
            + strdxr(1,jc,kc)*(dtdtxr(1,jc,kc)/strtxr(1,jc,kc)  &
            - sorpxr(1,jc,kc)/strpxr(1,jc,kc) + sydtxr(1,jc,kc))  &
            - bcl2xr(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl2xr(1,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxr == nsbci3)THEN
    
!         INFLOW BOUNDARY CONDITION No 3
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         SPECIFY L's AS REQUIRED
!         L1X-L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxr(1,jc,kc)*acouxr(1,jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl2xr(1,jc,kc) = struxr(1,jc,kc)  &
            *(bcl2xr(1,jc,kc)-bcl5xr(1,jc,kc)*ova2xr(1,jc,kc))
        bcl3xr(1,jc,kc) = struxr(1,jc,kc)*bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = struxr(1,jc,kc)*bcl4xr(1,jc,kc)
        bcl5xr(1,jc,kc) = half*(struxr(1,jc,kc)+acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5X UNCHANGED
        fornow = bcl5xr(1,jc,kc) + strdxr(1,jc,kc)*acouxr(1,jc,kc)*dudtxr(1,jc,kc)
        bcl1xr(1,jc,kc) = fornow - bcl1xr(1,jc,kc)
        bcl2xr(1,jc,kc) = -dddtxr(1,jc,kc)  &
            - ova2xr(1,jc,kc)*(bcl5xr(1,jc,kc)+fornow) - bcl2xr(1,jc,kc)
        bcl3xr(1,jc,kc) = -dvdtxr(1,jc,kc) - bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = -dwdtxr(1,jc,kc) - bcl4xr(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*(ova2xr(1,jc,kc)*strexr(1,jc,kc)  &
            - struxr(1,jc,kc)/acouxr(1,jc,kc) + ovgmxr(1,jc,kc))  &
            - bcl2xr(1,jc,kc)*strexr(1,jc,kc)  &
            - bcl3xr(1,jc,kc)*strdxr(1,jc,kc)*strvxr(1,jc,kc)  &
            - bcl4xr(1,jc,kc)*strdxr(1,jc,kc)*strwxr(1,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          bclyxr(ispec,1,jc,kc) = ratexr(ispec,1,jc,kc)/strdxr(1,jc,kc)  &
              - dydtxr(ispec,1,jc,kc) - struxr(1,jc,kc)*bclyxr(ispec,1,jc,kc)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
              - bclyxr(ispec,1,jc,kc)*strdxr(1,jc,kc)*strhxr(ispec,1,jc,kc)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       WALL BOUNDARY CONDITIONS
!       ------------------------
  
  IF(nsbcxr == nsbcw1)THEN
    
!         WALL BOUNDARY CONDITION No 1
!         NO-SLIP WALL - ADIABATIC
    
!         ALL VELOCITY COMPONENTS IMPOSED
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         SPECIFY L's AS REQUIRED
!         L1X,L3X-L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxr(1,jc,kc)*acouxr(1,jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl3xr(1,jc,kc) = struxr(1,jc,kc)*bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = struxr(1,jc,kc)*bcl4xr(1,jc,kc)
        bcl5xr(1,jc,kc) = half*(struxr(1,jc,kc)+acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L2X,L5X UNCHANGED
        bcl1xr(1,jc,kc) = bcl5xr(1,jc,kc)  &
            + strdxr(1,jc,kc)*acouxr(1,jc,kc)*dudtxr(1,jc,kc) - bcl1xr(1,jc,kc)
        bcl3xr(1,jc,kc) = -dvdtxr(1,jc,kc) - bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = -dwdtxr(1,jc,kc) - bcl4xr(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*(ova2xr(1,jc,kc)*strexr(1,jc,kc)  &
            + struxr(1,jc,kc)/acouxr(1,jc,kc) + ovgmxr(1,jc,kc))  &
            - bcl3xr(1,jc,kc)*strdxr(1,jc,kc)*strvxr(1,jc,kc)  &
            - bcl4xr(1,jc,kc)*strdxr(1,jc,kc)*strwxr(1,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(ispec,istol,jc,kc) = yrhs(ispec,istol,jc,kc)  &
              - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)*stryxr(ispec,1,jc,kc)
          
        END DO
      END DO
      
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
            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_xdir, "A_yz = A_yz + B_mulditim_yz*C_multidim_yz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratexr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_xdir, "A_yz = -A_yz*B_yz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1xr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1X-L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxr(1,jc,kc)*acouxr(1,jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl2xr(1,jc,kc) = struxr(1,jc,kc)  &
            *(bcl2xr(1,jc,kc)-bcl5xr(1,jc,kc)*ova2xr(1,jc,kc))
        bcl3xr(1,jc,kc) = struxr(1,jc,kc)*bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = struxr(1,jc,kc)*bcl4xr(1,jc,kc)
        bcl5xr(1,jc,kc) = half*(struxr(1,jc,kc)+acouxr(1,jc,kc))  &
            *(bcl5xr(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5X UNCHANGED
        bcl1xr(1,jc,kc) = bcl5xr(1,jc,kc)  &
            + strdxr(1,jc,kc)*acouxr(1,jc,kc)*dudtxr(1,jc,kc) - bcl1xr(1,jc,kc)
        bcl3xr(1,jc,kc) = -dvdtxr(1,jc,kc) - bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = -dwdtxr(1,jc,kc) - bcl4xr(1,jc,kc)
        bcl2xr(1,jc,kc) = gam1xr(1,jc,kc)*ova2xr(1,jc,kc)  &
            *(bcl1xr(1,jc,kc)+bcl5xr(1,jc,kc))  &
            + strdxr(1,jc,kc)*(dtdtxr(1,jc,kc)/strtxr(1,jc,kc)  &
            - sorpxr(1,jc,kc)/strpxr(1,jc,kc)) - bcl2xr(1,jc,kc)
        
      END DO
    END DO
    
!         LYX
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
!               OLD VALUE OF LYX
          bclyxr(ispec,1,jc,kc) = struxr(1,jc,kc)*bclyxr(ispec,1,jc,kc)
          
!               UPDATE L2X
          bcl2xr(1,jc,kc) = bcl2xr(1,jc,kc) + (ratexr(ispec,1,jc,kc)  &
              - strdxr(1,jc,kc)*bclyxr(ispec,1,jc,kc)) *rgspec(ispec)/strrxr(1,jc,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(1,jc,kc)*ova2xr(1,jc,kc)  &
            - bcl2xr(1,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(ispec,istol,jc,kc) = yrhs(ispec,istol,jc,kc)  &
              - (bcl2xr(1,jc,kc)+bcl1xr(1,jc,kc)*ova2xr(1,jc,kc))*stryxr(ispec,1,jc,kc)
          
        END DO
      END DO
      
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
    IF(fylcnv)THEN
  
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
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_reduced_enthalpy_ydir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))    

        END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = (/istal,istol,1,1,kstal,kstol/)
        call ops_par_loop(bounds_kernel_reduced_energy_ydir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_streyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_ovgmyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      fornow = strgyl(ic,1,kc)*gam1yl(ic,1,kc)*strtyl(ic,1,kc)
      acouyl(ic,1,kc) = SQRT(fornow)
      ova2yl(ic,1,kc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbcyl == nsbco1) THEN

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

            DO ispec = 1, nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_rateyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SPECIFY L5Y AS REQUIRED
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_outflowBC1_computeL_yl, "SPECIFY L5Y AS REQUIRED", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_gbl(cobcyl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(pinfyl, 1, "real(dp)", OPS_READ))
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)
        
        urhs(ic,jstal,kc) = urhs(ic,jstal,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)*struyl(ic,1,kc)
        
        vrhs(ic,jstal,kc) = vrhs(ic,jstal,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)*(strvyl(ic,1,kc)+acouyl(ic,1,kc))
        
        wrhs(ic,jstal,kc) = wrhs(ic,jstal,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)*strwyl(ic,1,kc)
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
            - bcl5yl(ic,1,kc)*(ova2yl(ic,1,kc)*streyl(ic,1,kc)  &
            + strvyl(ic,1,kc)/acouyl(ic,1,kc) + ovgmyl(ic,1,kc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jstal,kc) = yrhs(ispec,ic,jstal,kc)  &
              - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)*stryyl(ispec,ic,1,kc)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbcyl == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_rateyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L2Y-L5Y
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_inflowBC1_computeL_yl, "L2Y to L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_gbl(cobcyl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(pinfyl, 1, "real(dp)", OPS_READ))

!           LYY
            DO ispec = 1,nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_inflowBC1_LYY_yl, "LYY", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_rateyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)
        
        urhs(ic,jstal,kc) = urhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*struyl(ic,1,kc)  &
            - bcl3yl(ic,1,kc)*strdyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)*struyl(ic,1,kc)
        
        vrhs(ic,jstal,kc) = vrhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*strvyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)*(strvyl(ic,1,kc)+acouyl(ic,1,kc))
        
        wrhs(ic,jstal,kc) = wrhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*strwyl(ic,1,kc)  &
            - bcl4yl(ic,1,kc)*strdyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)*strwyl(ic,1,kc)
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*streyl(ic,1,kc)  &
            - bcl3yl(ic,1,kc)*strdyl(ic,1,kc)*struyl(ic,1,kc)  &
            - bcl4yl(ic,1,kc)*strdyl(ic,1,kc)*strwyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*(ova2yl(ic,1,kc)*streyl(ic,1,kc)  &
            + strvyl(ic,1,kc)/acouyl(ic,1,kc) + ovgmyl(ic,1,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = bclyyl(ispec,ic,1,kc)*strdyl(ic,1,kc)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - fornow*strhyl(ispec,ic,1,kc)
          
          yrhs(ispec,ic,jstal,kc) = yrhs(ispec,ic,jstal,kc)  &
              - (bcl2yl(ic,1,kc)+bcl5yl(ic,1,kc)*ova2yl(ic,1,kc))*stryyl(ispec,ic,1,kc)  &
              - fornow
          
        END DO
      END DO
      
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
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqB_ydir, "A_xz = A_xz + B_mulditim_xz*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_dydtyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))
                call ops_par_loop(bounds_kernel_eqA_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_rateyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz/B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqC_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y,L2Y,L5Y
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_inflowBC2_computeL_yl, "L1Y L2Y L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dtdtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sydtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)
        
      END DO
    END DO
    
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
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_inflowBC3_computeL_yl, "L1Y to L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dddtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*streyl(ic,1,kc)  &
            - bcl3yl(ic,1,kc)*strdyl(ic,1,kc)*struyl(ic,1,kc)  &
            - bcl4yl(ic,1,kc)*strdyl(ic,1,kc)*strwyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*(ova2yl(ic,1,kc)*streyl(ic,1,kc)  &
            + strvyl(ic,1,kc)/acouyl(ic,1,kc) + ovgmyl(ic,1,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          bclyyl(ispec,ic,1,kc) = rateyl(ispec,ic,1,kc)/strdyl(ic,1,kc)  &
              - dydtyl(ispec,ic,1,kc) - strvyl(ic,1,kc)*bclyyl(ispec,ic,1,kc)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
              - bclyyl(ispec,ic,1,kc)*strdyl(ic,1,kc)*strhyl(ispec,ic,1,kc)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       WALL BOUNDARY CONDITIONS
!       ------------------------
  
        IF(nsbcyl == nsbcw1)THEN
    
!           WALL BOUNDARY CONDITION No 1
!           NO-SLIP WALL - ADIABATIC
    
!           ALL VELOCITY COMPONENTS IMPOSED
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT
    
!           SPECIFY L's AS REQUIRED
!           L1Y,L3Y-L5Y
            rangexyz = (/istal,istol,1,1,kstal,kstol/)    
            call ops_par_loop(bounds_kernel_wallBC1_computeL_yl, "L1Y and L3Y to L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
            - bcl3yl(ic,1,kc)*strdyl(ic,1,kc)*struyl(ic,1,kc)  &
            - bcl4yl(ic,1,kc)*strdyl(ic,1,kc)*strwyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*(ova2yl(ic,1,kc)*streyl(ic,1,kc)  &
            + strvyl(ic,1,kc)/acouyl(ic,1,kc) + ovgmyl(ic,1,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jstal,kc) = yrhs(ispec,ic,jstal,kc)  &
              - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)*stryyl(ispec,ic,1,kc)
          
        END DO
      END DO
      
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
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_rateyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y-L5Y
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_wallBC2_computeL_yl, "L1Y to L5Y", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_gam1yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dtdtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           LYY
            DO ispec = 1,nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_wallBC2_LYY_yl, "LYY", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_bcl2yl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_rateyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strryl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,1,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jstal,kc) = yrhs(ispec,ic,jstal,kc)  &
              - (bcl2yl(ic,1,kc)+bcl5yl(ic,1,kc)*ova2yl(ic,1,kc))*stryyl(ispec,ic,1,kc)
          
        END DO
      END DO
      
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
    IF(fyrcnv)THEN
  
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
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_reduced_enthalpy_ydir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = (/istal,istol,1,1,kstal,kstol/)
        call ops_par_loop(bounds_kernel_reduced_energy_ydir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_streyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_ovgmyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      fornow = strgyr(ic,1,kc)*gam1yr(ic,1,kc)*strtyr(ic,1,kc)
      acouyr(ic,1,kc) = SQRT(fornow)
      ova2yr(ic,1,kc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbcyr == nsbco1) THEN

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_rateyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SPECIFY L1Y AS REQUIRED
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L1Y
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)-strdyr(ic,1,kc)*acouyr(ic,1,kc)*bcl1yr(ic,1,kc))
        
!             SUBTRACT FROM NEW VALUE OF L1Y
        bcl1yr(ic,1,kc)= half*sorpyr(ic,1,kc)  &
            + cobcyr*acouyr(ic,1,kc)*(strpyr(ic,1,kc)-pinfyr) - bcl1yr(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)
        
        urhs(ic,jstol,kc) = urhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)*struyr(ic,1,kc)
        
        vrhs(ic,jstol,kc) = vrhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)*(strvyr(ic,1,kc)-acouyr(ic,1,kc))
        
        wrhs(ic,jstol,kc) = wrhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)*strwyr(ic,1,kc)
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*(ova2yr(ic,1,kc)*streyr(ic,1,kc)  &
            - strvyr(ic,1,kc)/acouyr(ic,1,kc) + ovgmyr(ic,1,kc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jstol,kc) = yrhs(ispec,ic,jstol,kc)  &
              - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)*stryyr(ispec,ic,1,kc)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbcyr == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_rateyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y-L4Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyr(ic,1,kc)*acouyr(ic,1,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl2yr(ic,1,kc) = strvyr(ic,1,kc)  &
            *(bcl2yr(ic,1,kc)-bcl5yr(ic,1,kc)*ova2yr(ic,1,kc))
        bcl3yr(ic,1,kc) = strvyr(ic,1,kc)*bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = strvyr(ic,1,kc)*bcl4yr(ic,1,kc)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
!             L5Y UNCHANGED
        bcl1yr(ic,1,kc) = half*sorpyr(ic,1,kc)  &
            + cobcyr*acouyr(ic,1,kc)*(strpyr(ic,1,kc)-pinfyr) - bcl1yr(ic,1,kc)
        bcl2yr(ic,1,kc) = -bcl2yr(ic,1,kc)
        bcl3yr(ic,1,kc) = -bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = -bcl4yr(ic,1,kc)
        
      END DO
    END DO
    
!         LYY
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
!               OLD VALUE OF L's
          bclyyr(ispec,ic,1,kc) = strvyr(ic,1,kc)*bclyyr(ispec,ic,1,kc)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
          bclyyr(ispec,ic,1,kc) = rateyr(ispec,ic,1,kc)/strdyr(ic,1,kc)  &
              - bclyyr(ispec,ic,1,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)  &
            - bcl2yr(ic,1,kc)
        
        urhs(ic,jstol,kc) = urhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)*struyr(ic,1,kc)  &
            - bcl2yr(ic,1,kc)*struyr(ic,1,kc) - bcl3yr(ic,1,kc)*strdyr(ic,1,kc)
        
        vrhs(ic,jstol,kc) = vrhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)*(strvyr(ic,1,kc)-acouyr(ic,1,kc))  &
            - bcl2yr(ic,1,kc)*strvyr(ic,1,kc)
        
        wrhs(ic,jstol,kc) = wrhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)*strwyr(ic,1,kc)  &
            - bcl2yr(ic,1,kc)*strwyr(ic,1,kc) - bcl4yr(ic,1,kc)*strdyr(ic,1,kc)
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*(ova2yr(ic,1,kc)*streyr(ic,1,kc)  &
            + strvyr(ic,1,kc)/acouyr(ic,1,kc) + ovgmyr(ic,1,kc))  &
            - bcl2yr(ic,1,kc)*streyr(ic,1,kc)  &
            - bcl3yr(ic,1,kc)*strdyr(ic,1,kc)*struyr(ic,1,kc)  &
            - bcl4yr(ic,1,kc)*strdyr(ic,1,kc)*strwyr(ic,1,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = bclyyr(ispec,ic,1,kc)*strdyr(ic,1,kc)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc) - fornow*strhyr(ispec,ic,1,kc)
          
          yrhs(ispec,ic,jstol,kc) = yrhs(ispec,ic,jstol,kc)  &
              - (bcl2yr(ic,1,kc)+bcl1yr(ic,1,kc)*ova2yr(ic,1,kc))*stryyr(ispec,ic,1,kc)  &
              - fornow
          
        END DO
      END DO
      
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
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqB_ydir, "A_xz = A_xz + B_mulditim_xz*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_dydtyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))
                call ops_par_loop(bounds_kernel_eqA_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_rateyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqD_ydir, "A_xz = A_xz/B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strryr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqC_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y,L2Y,L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyr(ic,1,kc)*acouyr(ic,1,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl2yr(ic,1,kc) = strvyr(ic,1,kc)  &
            *(bcl2yr(ic,1,kc)-bcl5yr(ic,1,kc)*ova2yr(ic,1,kc))
        bcl5yr(ic,1,kc) = half*(strvyr(ic,1,kc)+acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Y UNCHANGED
        bcl1yr(ic,1,kc) = bcl5yr(ic,1,kc)  &
            + strdyr(ic,1,kc)*acouyr(ic,1,kc)*dvdtyr(ic,1,kc) - bcl1yr(ic,1,kc)
        bcl2yr(ic,1,kc) = gam1yr(ic,1,kc)*ova2yr(ic,1,kc)  &
            *(bcl1yr(ic,1,kc)+bcl5yr(ic,1,kc))  &
            + strdyr(ic,1,kc)*(dtdtyr(ic,1,kc)/strtyr(ic,1,kc)  &
            - sorpyr(ic,1,kc)/strpyr(ic,1,kc) + sydtyr(ic,1,kc))  &
            - bcl2yr(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)  &
            - bcl2yr(ic,1,kc)
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyr == nsbci3)THEN
    
!         INFLOW BOUNDARY CONDITION No 3
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         SPECIFY L's AS REQUIRED
!         L1Y-L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyr(ic,1,kc)*acouyr(ic,1,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl2yr(ic,1,kc) = strvyr(ic,1,kc)  &
            *(bcl2yr(ic,1,kc)-bcl5yr(ic,1,kc)*ova2yr(ic,1,kc))
        bcl3yr(ic,1,kc) = strvyr(ic,1,kc)*bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = strvyr(ic,1,kc)*bcl4yr(ic,1,kc)
        bcl5yr(ic,1,kc) = half*(strvyr(ic,1,kc)+acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Y UNCHANGED
        fornow = bcl5yr(ic,1,kc) + strdyr(ic,1,kc)*acouyr(ic,1,kc)*dvdtyr(ic,1,kc)
        bcl1yr(ic,1,kc) = fornow - bcl1yr(ic,1,kc)
        bcl2yr(ic,1,kc) = -dddtyr(ic,1,kc)  &
            - ova2yr(ic,1,kc)*(bcl1yr(ic,1,kc)+fornow) - bcl2yr(ic,1,kc)
        bcl3yr(ic,1,kc) = -dudtyr(ic,1,kc) - bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = -dwdtyr(ic,1,kc) - bcl4yr(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*(ova2yr(ic,1,kc)*streyr(ic,1,kc)  &
            + strvyr(ic,1,kc)/acouyr(ic,1,kc) + ovgmyr(ic,1,kc))  &
            - bcl2yr(ic,1,kc)*streyr(ic,1,kc)  &
            - bcl3yr(ic,1,kc)*strdyr(ic,1,kc)*struyr(ic,1,kc)  &
            - bcl4yr(ic,1,kc)*strdyr(ic,1,kc)*strwyr(ic,1,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          bclyyr(ispec,ic,1,kc) = rateyr(ispec,ic,1,kc)/strdyr(ic,1,kc)  &
              - dydtyr(ispec,ic,1,kc) - strvyr(ic,1,kc)*bclyyr(ispec,ic,1,kc)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
              - bclyyr(ispec,ic,1,kc)*strdyr(ic,1,kc)*strhyr(ispec,ic,1,kc)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       WALL BOUNDARY CONDITIONS
!       ------------------------
  
  IF(nsbcyr == nsbcw1)THEN
    
!         WALL BOUNDARY CONDITION No 1
!         NO-SLIP WALL - ADIABATIC
    
!         ALL VELOCITY COMPONENTS IMPOSED
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         SPECIFY L's AS REQUIRED
!         L1Y,L3Y-L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyr(ic,1,kc)*acouyr(ic,1,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl3yr(ic,1,kc) = strvyr(ic,1,kc)*bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = strvyr(ic,1,kc)*bcl4yr(ic,1,kc)
        bcl5yr(ic,1,kc) = half*(strvyr(ic,1,kc)+acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L2Y,L5Y UNCHANGED
        bcl1yr(ic,1,kc) = bcl5yr(ic,1,kc)  &
            + strdyr(ic,1,kc)*acouyr(ic,1,kc)*dvdtyr(ic,1,kc) - bcl1yr(ic,1,kc)
        bcl3yr(ic,1,kc) = -dudtyr(ic,1,kc) - bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = -dwdtyr(ic,1,kc) - bcl4yr(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*(ova2yr(ic,1,kc)*streyr(ic,1,kc)  &
            + strvyr(ic,1,kc)/acouyr(ic,1,kc) + ovgmyr(ic,1,kc))  &
            - bcl3yr(ic,1,kc)*strdyr(ic,1,kc)*struyr(ic,1,kc)  &
            - bcl4yr(ic,1,kc)*strdyr(ic,1,kc)*strwyr(ic,1,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jstol,kc) = yrhs(ispec,ic,jstol,kc)  &
              - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)*stryyr(ispec,ic,1,kc)
          
        END DO
      END DO
      
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
            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(set_zero_kernel_ydir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,1,1,kstal,kstol/)
                call ops_par_loop(bounds_kernel_eqA_ydir, "A_xz = A_xz + B_mulditim_xz*C_multidim_xz", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_rateyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,1,1,kstal,kstol/)
            call ops_par_loop(bounds_kernel_eqC_ydir, "A_xz = -A_xz*B_xz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1yr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Y-L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyr(ic,1,kc)*acouyr(ic,1,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl2yr(ic,1,kc) = strvyr(ic,1,kc)  &
            *(bcl2yr(ic,1,kc)-bcl5yr(ic,1,kc)*ova2yr(ic,1,kc))
        bcl3yr(ic,1,kc) = strvyr(ic,1,kc)*bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = strvyr(ic,1,kc)*bcl4yr(ic,1,kc)
        bcl5yr(ic,1,kc) = half*(strvyr(ic,1,kc)+acouyr(ic,1,kc))  &
            *(bcl5yr(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Y UNCHANGED
        bcl1yr(ic,1,kc) = bcl5yr(ic,1,kc)  &
            + strdyr(ic,1,kc)*acouyr(ic,1,kc)*dvdtyr(ic,1,kc) - bcl1yr(ic,1,kc)
        bcl3yr(ic,1,kc) = -dudtyr(ic,1,kc) - bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = -dwdtyr(ic,1,kc) - bcl4yr(ic,1,kc)
        bcl2yr(ic,1,kc) = gam1yr(ic,1,kc)*ova2yr(ic,1,kc)  &
            *(bcl1yr(ic,1,kc)+bcl5yr(ic,1,kc))  &
            + strdyr(ic,1,kc)*(dtdtyr(ic,1,kc)/strtyr(ic,1,kc)  &
            - sorpyr(ic,1,kc)/strpyr(ic,1,kc)) - bcl2yr(ic,1,kc)
        
      END DO
    END DO
    
!         LYY
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYY
          bclyyr(ispec,ic,1,kc) = strvyr(ic,1,kc)*bclyyr(ispec,ic,1,kc)
          
!               UPDATE L2Y
          bcl2yr(ic,1,kc) = bcl2yr(ic,1,kc) + (rateyr(ispec,ic,1,kc)  &
              - strdyr(ic,1,kc)*bclyyr(ispec,ic,1,kc)) *rgspec(ispec)/strryr(ic,1,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,1,kc)  &
            - bcl2yr(ic,1,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jstol,kc) = yrhs(ispec,ic,jstol,kc)  &
              - (bcl2yr(ic,1,kc)+bcl1yr(ic,1,kc)*ova2yr(ic,1,kc))*stryyr(ispec,ic,1,kc)
          
        END DO
      END DO
      
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
    IF(fzlcnv)THEN
  
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
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_reduced_enthalpy_zdir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

        END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = (/istal,istol,jstal,jstol,1,1/)
        call ops_par_loop(bounds_kernel_reduced_energy_zdir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strezl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_ovgmzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!       SPEED OF SOUND
!       --------------
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = strgzl(ic,jc,1)*gam1zl(ic,jc,1)*strtzl(ic,jc,1)
      acouzl(ic,jc,1) = SQRT(fornow)
      ova2zl(ic,jc,1) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbczl == nsbco1) THEN

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_eqA_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratezl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_eqC_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SPECIFY L5Z AS REQUIRED
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_outflowBC1_computeL_zl, "SPECIFY L5Z AS REQUIRED", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_gbl(cobczl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(pinfzl, 1, "real(dp)", OPS_READ))
    
!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)
        
        urhs(ic,jc,kstal) = urhs(ic,jc,kstal)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)*struzl(ic,jc,1)
        
        vrhs(ic,jc,kstal) = vrhs(ic,jc,kstal)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)*strvzl(ic,jc,1)
        
        wrhs(ic,jc,kstal) = wrhs(ic,jc,kstal)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)*(strwzl(ic,jc,1)+acouzl(ic,jc,1))
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
            - bcl5zl(ic,jc,1)*(ova2zl(ic,jc,1)*strezl(ic,jc,1)  &
            + strwzl(ic,jc,1)/acouzl(ic,jc,1) + ovgmzl(ic,jc,1))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jc,kstal) = yrhs(ispec,ic,jc,kstal)  &
              - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)*stryzl(ispec,ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbczl == nsbci1) THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_eqA_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratezl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_eqC_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L2Z-L5Z
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_inflowBC1_computeL_yl, "L2Z to L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_gbl(cobczl, 1, "real(dp)", OPS_READ), &
                            ops_arg_gbl(pinfzl, 1, "real(dp)", OPS_READ))

!           LYZ
            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_inflowBC1_LYZ_zl, "LYZ", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_ratezl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO
    
!           ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)
        
        urhs(ic,jc,kstal) = urhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*struzl(ic,jc,1)  &
            - bcl3zl(ic,jc,1)*strdzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)*struzl(ic,jc,1)
        
        vrhs(ic,jc,kstal) = vrhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*strvzl(ic,jc,1)  &
            - bcl4zl(ic,jc,1)*strdzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)*strvzl(ic,jc,1)
        
        wrhs(ic,jc,kstal) = wrhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*strwzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)*(strwzl(ic,jc,1)+acouzl(ic,jc,1))
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*strezl(ic,jc,1)  &
            - bcl3zl(ic,jc,1)*strdzl(ic,jc,1)*struzl(ic,jc,1)  &
            - bcl4zl(ic,jc,1)*strdzl(ic,jc,1)*strvzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*(ova2zl(ic,jc,1)*strezl(ic,jc,1)  &
            + strwzl(ic,jc,1)/acouzl(ic,jc,1) + ovgmzl(ic,jc,1))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = bclyzl(ispec,ic,jc,1)*strdzl(ic,jc,1)
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - fornow*strhzl(ispec,ic,jc,1)
          
          yrhs(ispec,ic,jc,kstal) = yrhs(ispec,ic,jc,kstal)  &
              - (bcl2zl(ic,jc,1)+bcl5zl(ic,jc,1)*ova2zl(ic,jc,1))*stryzl(ispec,ic,jc,1)  &
              - fornow
          
        END DO
      END DO
      
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
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_eqB_zdir, "A_xy = A_xy + B_mulditim_xy*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_dydtzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))
                call ops_par_loop(bounds_kernel_eqA_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratezl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy/B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqC_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z,L2Z,L5Z
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_inflowBC2_computeL_zl, "L1Z L2Z L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dtdtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sydtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
        IF(nsbczl == nsbci3)THEN
    
!           INFLOW BOUNDARY CONDITION No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!           VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT
    
!           SPECIFY L's AS REQUIRED
!           L1Z-L5Z
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_inflowBC3_computeL_zl, "L1Z to L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dddtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           ADD TO CONSERVATIVE SOURCE TERMS
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*strezl(ic,jc,1)  &
            - bcl3zl(ic,jc,1)*strdzl(ic,jc,1)*struzl(ic,jc,1)  &
            - bcl4zl(ic,jc,1)*strdzl(ic,jc,1)*strvzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*(ova2zl(ic,jc,1)*strezl(ic,jc,1)  &
            + strwzl(ic,jc,1)/acouzl(ic,jc,1) + ovgmzl(ic,jc,1))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          bclyzl(ispec,ic,jc,1) = ratezl(ispec,ic,jc,1)/strdzl(ic,jc,1)  &
              - dydtzl(ispec,ic,jc,1) - strwzl(ic,jc,1)*bclyzl(ispec,ic,jc,1)
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
              - bclyzl(ispec,ic,jc,1)*strdzl(ic,jc,1)*strhzl(ispec,ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       WALL BOUNDARY CONDITIONS
!       ------------------------
  
        IF(nsbczl == nsbcw1) THEN
    
!           WALL BOUNDARY CONDITION No 1
!           NO-SLIP WALL - ADIABATIC
    
!           ALL VELOCITY COMPONENTS IMPOSED
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT
    
!           SPECIFY L's AS REQUIRED
!           L1Z,L3Z-L5Z
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_wallBC1_computeL_zl, "L1Z and L3Z to L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
            - bcl3zl(ic,jc,1)*strdzl(ic,jc,1)*struzl(ic,jc,1)  &
            - bcl4zl(ic,jc,1)*strdzl(ic,jc,1)*strvzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*(ova2zl(ic,jc,1)*strezl(ic,jc,1)  &
            + strwzl(ic,jc,1)/acouzl(ic,jc,1) + ovgmzl(ic,jc,1))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jc,kstal) = yrhs(ispec,ic,jc,kstal)  &
              - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)*stryzl(ispec,ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================

        IF(nsbczl == nsbcw2)THEN

!           WALL BOUNDARY CONDITION No 2
!           NO-SLIP WALL - ISOTHERMAL

!           VELOCITY AND TEMPERATURE IMPOSED
!           AS FUNCTIONS OF TIME
!           VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!           SET IN SUBROUTINE BOUNDT

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_eqA_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratezl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_eqC_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))
 
!           SPECIFY L's AS REQUIRED
!           L1Z-L5Z
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_wallBC2_computeL_zl, "L1Z to L5Z", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_bcl1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl3zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl4zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_bcl5zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_acouzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_gam1zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_ova2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dtdtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_sorpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strpzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dudtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dvdtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_dwdtzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))            

!           LYZ
            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_wallBC2_LYZ_zl, "LYZ", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_bclyzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_bcl2zl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_ratezl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strrzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc,1)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jc,kstal) = yrhs(ispec,ic,jc,kstal)  &
              - (bcl2zl(ic,jc,1)+bcl5zl(ic,jc,1)*ova2zl(ic,jc,1))*stryzl(ispec,ic,jc,1)
          
        END DO
      END DO
      
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
    IF(fzrcnv)THEN
  
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
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_reduced_enthalpy_zdir, "REDUCED SPECIES ENTHALPY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strhzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))
    
        END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
        rangexyz = (/istal,istol,jstal,jstol,1,1/)
        call ops_par_loop(bounds_kernel_reduced_energy_zdir, "REDUCED INTERNAL ENERGY", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strezr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_ovgmzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strgzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!       SPEED OF SOUND
!       --------------
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = strgzr(ic,jc,1)*gam1zr(ic,jc,1)*strtzr(ic,jc,1)
      acouzr(ic,jc,1) = SQRT(fornow)
      ova2zr(ic,jc,1) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

        IF(nsbczr == nsbco1) THEN

!           OUTFLOW BC No 1
!           SUBSONIC NON-REFLECTING OUTFLOW
!           WITH OPTION TO SET PRESSURE AT INFINITY

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_eqA_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratezr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_eqC_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SPECIFY L1Z AS REQUIRED
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L1Z
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc,1))  &
            *(bcl5zr(ic,jc,1)-strdzr(ic,jc,1)*acouzr(ic,jc,1)*bcl1zr(ic,jc,1))
        
!             SUBTRACT FROM NEW VALUE OF L1Z
        bcl1zr(ic,jc,1)= half*sorpzr(ic,jc,1)  &
            + cobczr*acouzr(ic,jc,1)*(strpzr(ic,jc,1)-pinfzr) - bcl1zr(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)
        
        urhs(ic,jc,kstol) = urhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)*struzr(ic,jc,1)
        
        vrhs(ic,jc,kstol) = vrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)*strvzr(ic,jc,1)
        
        wrhs(ic,jc,kstol) = wrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)*(strwzr(ic,jc,1)-acouzr(ic,jc,1))
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*(ova2zr(ic,jc,1)*strezr(ic,jc,1)  &
            - strwzr(ic,jc,1)/acouzr(ic,jc,1) + ovgmzr(ic,jc,1))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jc,kstol) = yrhs(ispec,ic,jc,kstol)  &
              - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)*stryzr(ispec,ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

        IF(nsbczr == nsbci1)THEN

!           INFLOW BC No 1
!           SUBSONIC NON-REFLECTING LAMINAR INFLOW

!           PRECOMPUTE CHEMISTRY TERMS
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_eqA_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratezr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_eqC_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z-L4Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc,1)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc,1))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl2zr(ic,jc,1) = strwzr(ic,jc,1)  &
            *(bcl2zr(ic,jc,1)-bcl5zr(ic,jc,1)*ova2zr(ic,jc,1))
        bcl3zr(ic,jc,1) = strwzr(ic,jc,1)*bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = strwzr(ic,jc,1)*bcl4zr(ic,jc,1)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
!             L5Z UNCHANGED
        bcl1zr(ic,jc,1) = half*sorpzr(ic,jc,1)  &
            + cobczr*acouzr(ic,jc,1)*(strpzr(ic,jc,1)-pinfzr) - bcl1zr(ic,jc,1)
        bcl2zr(ic,jc,1) = -bcl2zr(ic,jc,1)
        bcl3zr(ic,jc,1) = -bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = -bcl4zr(ic,jc,1)
        
      END DO
    END DO
    
!         LYZ
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
!               OLD VALUE OF L's
          bclyzr(ispec,ic,jc,1) = strwzr(ic,jc,1)*bclyzr(ispec,ic,jc,1)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
          bclyzr(ispec,ic,jc,1) = ratezr(ispec,ic,jc,1)/strdzr(ic,jc,1)  &
              - bclyzr(ispec,ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)  &
            - bcl2zr(ic,jc,1)
        
        urhs(ic,jc,kstol) = urhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)*struzr(ic,jc,1)  &
            - bcl2zr(ic,jc,1)*struzr(ic,jc,1) - bcl3zr(ic,jc,1)*strdzr(ic,jc,1)
        
        vrhs(ic,jc,kstol) = vrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)*strvzr(ic,jc,1)  &
            - bcl2zr(ic,jc,1)*strvzr(ic,jc,1) - bcl4zr(ic,jc,1)*strdzr(ic,jc,1)
        
        wrhs(ic,jc,kstol) = wrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)*(strwzr(ic,jc,1)-acouzr(ic,jc,1))  &
            - bcl2zr(ic,jc,1)*strwzr(ic,jc,1)
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*(ova2zr(ic,jc,1)*strezr(ic,jc,1)  &
            + strwzr(ic,jc,1)/acouzr(ic,jc,1) + ovgmzr(ic,jc,1))  &
            - bcl2zr(ic,jc,1)*strezr(ic,jc,1)  &
            - bcl3zr(ic,jc,1)*strdzr(ic,jc,1)*struzr(ic,jc,1)  &
            - bcl4zr(ic,jc,1)*strdzr(ic,jc,1)*strvzr(ic,jc,1)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = bclyzr(ispec,ic,jc,1)*strdzr(ic,jc,1)
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol) - fornow*strhzr(ispec,ic,jc,1)
          
          yrhs(ispec,ic,jc,kstol) = yrhs(ispec,ic,jc,kstol)  &
              - (bcl2zr(ic,jc,1)+bcl5zr(ic,jc,1)*ova2zr(ic,jc,1))*stryzr(ispec,ic,jc,1)  &
              - fornow
          
        END DO
      END DO
      
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
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_eqB_zdir, "A_xy = A_xy + B_mulditim_xy*val1", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_dydtzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))
                call ops_par_loop(bounds_kernel_eqA_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratezr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_eqD_zdir, "A_xy = A_xy/B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sydtzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_strrzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))
            call ops_par_loop(bounds_kernel_eqC_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z,L2Z,L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc,1)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc,1))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl2zr(ic,jc,1) = strwzr(ic,jc,1)  &
            *(bcl2zr(ic,jc,1)-bcl5zr(ic,jc,1)*ova2zr(ic,jc,1))
        bcl5zr(ic,jc,1) = half*(strwzr(ic,jc,1)+acouzr(ic,jc,1))  &
            *(bcl5zr(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Z UNCHANGED
        bcl1zr(ic,jc,1) = bcl5zr(ic,jc,1)  &
            + strdzr(ic,jc,1)*acouzr(ic,jc,1)*dwdtzr(ic,jc,1) - bcl1zr(ic,jc,1)
        bcl2zr(ic,jc,1) = gam1zr(ic,jc,1)*ova2zr(ic,jc,1)  &
            *(bcl1zr(ic,jc,1)+bcl5zr(ic,jc,1))  &
            + strdzr(ic,jc,1)*(dtdtzr(ic,jc,1)/strtzr(ic,jc,1)  &
            - sorpzr(ic,jc,1)/strpzr(ic,jc,1) + sydtzr(ic,jc,1))  &
            - bcl2zr(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)  &
            - bcl2zr(ic,jc,1)
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczr == nsbci3)THEN
    
!         INFLOW BOUNDARY CONDITION No 3
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         SPECIFY L's AS REQUIRED
!         L1Z-L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc,1)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc,1))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl2zr(ic,jc,1) = strwzr(ic,jc,1)  &
            *(bcl2zr(ic,jc,1)-bcl5zr(ic,jc,1)*ova2zr(ic,jc,1))
        bcl3zr(ic,jc,1) = strwzr(ic,jc,1)*bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = strwzr(ic,jc,1)*bcl4zr(ic,jc,1)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Z UNCHANGED
        fornow = bcl5zr(ic,jc,1) + strdzr(ic,jc,1)*acouzr(ic,jc,1)*dwdtzr(ic,jc,1)
        bcl1zr(ic,jc,1) = fornow - bcl1zr(ic,jc,1)
        bcl2zr(ic,jc,1) = -dddtzr(ic,jc,1)  &
            - ova2zr(ic,jc,1)*(bcl1zr(ic,jc,1)+fornow) - bcl2zr(ic,jc,1)
        bcl3zr(ic,jc,1) = -dudtzr(ic,jc,1) - bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = -dvdtzr(ic,jc,1) - bcl4zr(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*(ova2zr(ic,jc,1)*strezr(ic,jc,1)  &
            + strwzr(ic,jc,1)/acouzr(ic,jc,1) + ovgmzr(ic,jc,1))  &
            - bcl2zr(ic,jc,1)*strezr(ic,jc,1)  &
            - bcl3zr(ic,jc,1)*strdzr(ic,jc,1)*struzr(ic,jc,1)  &
            - bcl4zr(ic,jc,1)*strdzr(ic,jc,1)*strvzr(ic,jc,1)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          bclyzr(ispec,ic,jc,1) = ratezr(ispec,ic,jc,1)/strdzr(ic,jc,1)  &
              - dydtzr(ispec,ic,jc,1) - strwzr(ic,jc,1)*bclyzr(ispec,ic,jc,1)
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
              - bclyzr(ispec,ic,jc,1)*strdzr(ic,jc,1)*strhzr(ispec,ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       WALL BOUNDARY CONDITIONS
!       ------------------------
  
  IF(nsbczr == nsbcw1)THEN
    
!         WALL BOUNDARY CONDITION No 1
!         NO-SLIP WALL - ADIABATIC
    
!         ALL VELOCITY COMPONENTS IMPOSED
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         SPECIFY L's AS REQUIRED
!         L1Z,L3Z-L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc,1)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc,1))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl3zr(ic,jc,1) = strwzr(ic,jc,1)*bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = strwzr(ic,jc,1)*bcl4zr(ic,jc,1)
        bcl5zr(ic,jc,1) = half*(strwzr(ic,jc,1)+acouzr(ic,jc,1))  &
            *(bcl5zr(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L2Z,L5Z UNCHANGED
        bcl1zr(ic,jc,1) = bcl5zr(ic,jc,1)  &
            + strdzr(ic,jc,1)*acouzr(ic,jc,1)*dwdtzr(ic,jc,1) - bcl1zr(ic,jc,1)
        bcl3zr(ic,jc,1) = -dudtzr(ic,jc,1) - bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = -dvdtzr(ic,jc,1) - bcl4zr(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*(ova2zr(ic,jc,1)*strezr(ic,jc,1)  &
            + strwzr(ic,jc,1)/acouzr(ic,jc,1) + ovgmzr(ic,jc,1))  &
            - bcl3zr(ic,jc,1)*strdzr(ic,jc,1)*struzr(ic,jc,1)  &
            - bcl4zr(ic,jc,1)*strdzr(ic,jc,1)*strvzr(ic,jc,1)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jc,kstol) = yrhs(ispec,ic,jc,kstol)  &
              - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)*stryzr(ispec,ic,jc,1)
          
        END DO
      END DO
      
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
            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(set_zero_kernel_zdir, "set_zero", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE))

            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,1,1/)
                call ops_par_loop(bounds_kernel_eqA_zdir, "A_xy = A_xy + B_mulditim_xy*C_multidim_xy", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                                ops_arg_dat(d_strhzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_ratezr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

            rangexyz = (/istal,istol,jstal,jstol,1,1/)
            call ops_par_loop(bounds_kernel_eqC_zdir, "A_xy = -A_xy*B_xy", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_sorpzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_WRITE), &
                            ops_arg_dat(d_gam1zr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SPECIFY L's AS REQUIRED
!           L1Z-L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc,1)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc,1))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl2zr(ic,jc,1) = strwzr(ic,jc,1)  &
            *(bcl2zr(ic,jc,1)-bcl5zr(ic,jc,1)*ova2zr(ic,jc,1))
        bcl3zr(ic,jc,1) = strwzr(ic,jc,1)*bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = strwzr(ic,jc,1)*bcl4zr(ic,jc,1)
        bcl5zr(ic,jc,1) = half*(strwzr(ic,jc,1)+acouzr(ic,jc,1))  &
            *(bcl5zr(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Z UNCHANGED
        bcl1zr(ic,jc,1) = bcl5zr(ic,jc,1)  &
            + strdzr(ic,jc,1)*acouzr(ic,jc,1)*dwdtzr(ic,jc,1) - bcl1zr(ic,jc,1)
        bcl3zr(ic,jc,1) = -dudtzr(ic,jc,1) - bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = -dvdtzr(ic,jc,1) - bcl4zr(ic,jc,1)
        bcl2zr(ic,jc,1) = gam1zr(ic,jc,1)*ova2zr(ic,jc,1)  &
            *(bcl1zr(ic,jc,1)+bcl5zr(ic,jc,1))  &
            + strdzr(ic,jc,1)*(dtdtzr(ic,jc,1)/strtzr(ic,jc,1)  &
            - sorpzr(ic,jc,1)/strpzr(ic,jc,1)) - bcl2zr(ic,jc,1)
        
      END DO
    END DO
    
!         LYZ
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYZ
          bclyzr(ispec,ic,jc,1) = strwzr(ic,jc,1)*bclyzr(ispec,ic,jc,1)
          
!               UPDATE L2Z
          bcl2zr(ic,jc,1) = bcl2zr(ic,jc,1) + (ratezr(ispec,ic,jc,1)  &
              - strdzr(ic,jc,1)*bclyzr(ispec,ic,jc,1)) *rgspec(ispec)/strrzr(ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc,1)  &
            - bcl2zr(ic,jc,1)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ispec,ic,jc,kstol) = yrhs(ispec,ic,jc,kstol)  &
              - (bcl2zr(ic,jc,1)+bcl1zr(ic,jc,1)*ova2zr(ic,jc,1))*stryzr(ispec,ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
    END IF
!   Z-DIRECTION RIGHT-HAND END

!   =========================================================================

END SUBROUTINE bounds
