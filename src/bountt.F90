SUBROUTINE bountt
 
    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BOUNTT
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   29-SEP-2003:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   SYNCHRONISES THE TIME-DEPENDENT BOUNDARY CONDITIONS

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    real(kind=dp) :: fornow
    integer :: jc,kc
    integer :: ispec
    integer :: iindex,ipower,icoef1,icoef2
    integer :: itint,icp
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   SYNCHRONISE AT CURRENT TIME STEP
!   --------------------------------
    irkstp = 1
    !KA   FIX INFLOW BC
    btime  = tstep
    fupelc = .true.
!   =========================================================================

!   X-DIRECTION LEFT-HAND END
!   -------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fxltrb) call bcutxl

!   LOCAL BC SUPPORT
    IF(fxlcnv) THEN
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
        IF(nsbcxl == nsbci2) THEN
    
!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxl
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttxl
    
!           SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        DO iindex = 1,nintmx
          itndex(iindex,istal,jc,kc) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          1000            CONTINUE
          IF(strtxl(1,jc,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 1000
            END IF
          END IF
!               END OF LOOP 1000
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,istal,jc,kc) = itndex(iindex,istal,jc,kc)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxl
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          itint = 1 +MOD(itndex(iindex,istal,jc,kc),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtxl(1,jc,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtxl(1,jc,kc)
          
          yrhs(ispec,istal,jc,kc) = drhs(istal,jc,kc)*stryxl(ispec,1,jc,kc)
          
          yrun(ispec,istal,jc,kc) = yrhs(ispec,istal,jc,kc)
          
          yerr(ispec,istal,jc,kc) = zero
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
              + (fornow-rgspec(ispec)*strtxl(1,jc,kc))*yrhs(ispec,istal,jc,kc)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))                    

        END IF
  
!       =======================================================================
  
        IF(nsbcxl == nsbci3) THEN
    
!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!           SET DENSITY AND TIME DERIVATIVE
            call bcdtxl
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxl
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqC_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxl
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bountt_kernel_eqE_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_stryxl, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF
  
!       =======================================================================
  
        IF(nsbcxl == nsbcw1)THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxl
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbcxl == nsbcw2)THEN
    
!           WALL BC No 2
!           NO-SLIP WALL - ISOTHERMAL
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxl
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttxl
    
!           SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        DO iindex = 1,nintmx
          itndex(iindex,istal,jc,kc) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          1100            CONTINUE
          IF(strtxl(1,jc,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 1100
            END IF
          END IF
!               END OF LOOP 1100
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,istal,jc,kc) = itndex(iindex,istal,jc,kc)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          itint = 1 +MOD(itndex(iindex,istal,jc,kc),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtxl(1,jc,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtxl(1,jc,kc)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
              + (fornow-rgspec(ispec)*strtxl(1,jc,kc))*yrhs(ispec,istal,jc,kc)
          
        END DO
      END DO
      
    END DO
   
            rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ)) 

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
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
        IF(nsbcxr == nsbci2) THEN
    
!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttxr
    
!           SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        DO iindex = 1,nintmx
          itndex(iindex,istol,jc,kc) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          1500            CONTINUE
          IF(strtxr(1,jc,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 1500
            END IF
          END IF
!               END OF LOOP 1500
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,istol,jc,kc) = itndex(iindex,istol,jc,kc)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))            

            call ops_par_loop(bountt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          itint = 1 +MOD(itndex(iindex,istol,jc,kc),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtxr(1,jc,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtxr(1,jc,kc)
          
          yrhs(ispec,istol,jc,kc) = drhs(istol,jc,kc)*stryxr(ispec,1,jc,kc)
          
          yrun(ispec,istol,jc,kc) = yrhs(ispec,istol,jc,kc)
          
          yerr(ispec,istol,jc,kc) = zero
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
              + (fornow-rgspec(ispec)*strtxr(1,jc,kc))*yrhs(ispec,istol,jc,kc)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbcxr == nsbci3) THEN
    
!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!           SET DENSITY AND TIME DERIVATIVE
            call bcdtxr
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxr
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqC_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
                call ops_par_loop(bountt_kernel_eqE_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_stryxr, 9, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbcxr == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxr
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbcxr == nsbcw2) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttxr
    
!           SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        DO iindex = 1,nintmx
          itndex(iindex,istol,jc,kc) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          1600            CONTINUE
          IF(strtxr(1,jc,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 1600
            END IF
          END IF
!               END OF LOOP 1600
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,istol,jc,kc) = itndex(iindex,istol,jc,kc)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(dp)", OPS_READ))

            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          itint = 1 +MOD(itndex(iindex,istol,jc,kc),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtxr(1,jc,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtxr(1,jc,kc)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
              + (fornow-rgspec(ispec)*strtxr(1,jc,kc))*yrhs(ispec,istol,jc,kc)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
    END IF
!   X-DIRECTION RIGHT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Y-DIRECTION LEFT-HAND END
!   -------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fyltrb) call bcutyl

!   LOCAL BC SUPPORT
    IF(fylcnv) THEN
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
        IF(nsbcyl == nsbci2) THEN
    
!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyl
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttyl
    
!           SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(iindex,ic,jstal,kc) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          2000            CONTINUE
          IF(strtyl(ic,1,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 2000
            END IF
          END IF
!               END OF LOOP 2000
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,ic,jstal,kc) = itndex(iindex,ic,jstal,kc)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyl
    
!           CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(iindex,ic,jstal,kc),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtyl(ic,1,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtyl(ic,1,kc)
          
          yrhs(ispec,ic,jstal,kc) = drhs(ic,jstal,kc)*stryyl(ispec,ic,1,kc)
          
          yrun(ispec,ic,jstal,kc) = yrhs(ispec,ic,jstal,kc)
          
          yerr(ispec,ic,jstal,kc) = zero
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
              + (fornow-rgspec(ispec)*strtyl(ic,1,kc))*yrhs(ispec,ic,jstal,kc)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbcyl == nsbci3) THEN
    
!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!           SET DENSITY AND TIME DERIVATIVE
            call bcdtyl
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyl
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqC_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyl
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
                call ops_par_loop(bountt_kernel_eqE_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_stryyl, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbcyl == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyl
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcyl == nsbcw2)THEN
    
!           WALL BC No 2
!           NO-SLIP WALL - ISOTHERMAL
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyl
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttyl
    
!           SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(iindex,ic,jstal,kc) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          2100            CONTINUE
          IF(strtyl(ic,1,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 2100
            END IF
          END IF
!               END OF LOOP 2100
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,ic,jstal,kc) = itndex(iindex,ic,jstal,kc)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

             call ops_par_loop(bountt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(iindex,ic,jstal,kc),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtyl(ic,1,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtyl(ic,1,kc)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
              + (fornow-rgspec(ispec)*strtyl(ic,1,kc))*yrhs(ispec,ic,jstal,kc)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
    END IF

!   Y-DIRECTION LEFT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Y-DIRECTION RIGHT-HAND END
!   --------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fyrtrb) call bcutyr

!   LOCAL BC SUPPORT
    IF(fyrcnv) THEN
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
        IF(nsbcyr == nsbci2) THEN
    
!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttyr
    
!           SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(iindex,ic,jstol,kc) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          2500            CONTINUE
          IF(strtyr(ic,1,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 2500
            END IF
          END IF
!               END OF LOOP 2500
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,ic,jstol,kc) = itndex(iindex,ic,jstol,kc)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

             call ops_par_loop(bountt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(iindex,ic,jstol,kc),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtyr(ic,1,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtyr(ic,1,kc)
          
          yrhs(ispec,ic,jstol,kc) = drhs(ic,jstol,kc)*stryyr(ispec,ic,1,kc)
          
          yrun(ispec,ic,jstol,kc) = yrhs(ispec,ic,jstol,kc)
          
          yerr(ispec,ic,jstol,kc) = zero
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
              + (fornow-rgspec(ispec)*strtyr(ic,1,kc))*yrhs(ispec,ic,jstol,kc)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbcyr == nsbci3) THEN
    
!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!           SET DENSITY AND TIME DERIVATIVE
            call bcdtyr
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyr
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqC_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
                call ops_par_loop(bountt_kernel_eqE_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_stryyr, 9, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbcyr == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyr
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

        END IF

!       =======================================================================
  
        IF(nsbcyr == nsbcw2) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttyr
    
!           SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(iindex,ic,jstol,kc) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          2600            CONTINUE
          IF(strtyr(ic,1,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 2600
            END IF
          END IF
!               END OF LOOP 2600
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,ic,jstol,kc) = itndex(iindex,ic,jstol,kc)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

             call ops_par_loop(bountt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(dp)", OPS_READ))

            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(iindex,ic,jstol,kc),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtyr(ic,1,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtyr(ic,1,kc)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
              + (fornow-rgspec(ispec)*strtyr(ic,1,kc))*yrhs(ispec,ic,jstol,kc)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
    END IF
!   Y-DIRECTION RIGHT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Z-DIRECTION LEFT-HAND END
!   -------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fzltrb) call bcutzl

!   LOCAL BC SUPPORT
    IF(fzlcnv) THEN
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
        IF(nsbczl == nsbci2) THEN
    
!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzl
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttzl
    
!           SET TEMPERATURE INTERVAL INDEX
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(iindex,ic,jc,kstal) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          3000            CONTINUE
          IF(strtzl(ic,jc,1) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 3000
            END IF
          END IF
!               END OF LOOP 3000
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,ic,jc,kstal) = itndex(iindex,ic,jc,kstal)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(bountt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

             call ops_par_loop(bountt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzl
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(iindex,ic,jc,kstal),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtzl(ic,jc,1) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtzl(ic,jc,1)
          
          yrhs(ispec,ic,jc,kstal) = drhs(ic,jc,kstal)*stryzl(ispec,ic,jc,1)
          
          yrun(ispec,ic,jc,kstal) = yrhs(ispec,ic,jc,kstal)
          
          yerr(ispec,ic,jc,kstal) = zero
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
              + (fornow-rgspec(ispec)*strtzl(ic,jc,1))*yrhs(ispec,ic,jc,kstal)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbczl == nsbci3) THEN
    
!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!           SET DENSITY AND TIME DERIVATIVE
            call bcdtzl
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzl
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(bountt_kernel_eqC_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzl
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
                call ops_par_loop(bountt_kernel_eqE_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_stryzl, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbczl == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzl
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(bountt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

        END IF

!       =======================================================================
 
        IF(nsbczl == nsbcw2)THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzl
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttzl
    
!           SET TEMPERATURE INTERVAL INDEX
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(iindex,ic,jc,kstal) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          3100            CONTINUE
          IF(strtzl(ic,jc,1) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 3100
            END IF
          END IF
!               END OF LOOP 3100
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,ic,jc,kstal) = itndex(iindex,ic,jc,kstal)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(bountt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(iindex,ic,jc,kstal),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtzl(ic,jc,1) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtzl(ic,jc,1)
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
              + (fornow-rgspec(ispec)*strtzl(ic,jc,1))*yrhs(ispec,ic,jc,kstal)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
    END IF
!   Z-DIRECTION LEFT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Z-DIRECTION RIGHT-HAND END
!   --------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fzrtrb) call bcutzr

!   LOCAL BC SUPPORT
    IF(fzrcnv) THEN
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE
  
!       =======================================================================
  
        IF(nsbczr == nsbci2) THEN
    
!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttzr
    
!           SET TEMPERATURE INTERVAL INDEX
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(iindex,ic,jc,kstol) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          3500            CONTINUE
          IF(strtzr(ic,jc,1) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 3500
            END IF
          END IF
!               END OF LOOP 3500
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,ic,jc,kstol) = itndex(iindex,ic,jc,kstol)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(bountt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(iindex,ic,jc,kstol),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtzr(ic,jc,1) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtzr(ic,jc,1)
          
          yrhs(ispec,ic,jc,kstol) = drhs(ic,jc,kstol)*stryzr(ispec,ic,jc,1)
          
          yrun(ispec,ic,jc,kstol) = yrhs(ispec,ic,jc,kstol)
          
          yerr(ispec,ic,jc,kstol) = zero
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
              + (fornow-rgspec(ispec)*strtzr(ic,jc,1))*yrhs(ispec,ic,jc,kstol)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbczr == nsbci3) THEN
    
!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!           SET DENSITY AND TIME DERIVATIVE
            call bcdtzr
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzr
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(bountt_kernel_eqC_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
                call ops_par_loop(bountt_kernel_eqE_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr, 9, s3d_000, "real(dp)", OPS_WRITE),  &
                                ops_arg_dat(d_stryzr, 9, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO
    
        END IF
  
!       =======================================================================
  
        IF(nsbczr == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzr
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(bountt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbczr == nsbcw2) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttzr
    
!           SET TEMPERATURE INTERVAL INDEX
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(iindex,ic,jc,kstol) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          3600            CONTINUE
          IF(strtzr(ic,jc,1) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 3600
            END IF
          END IF
!               END OF LOOP 3600
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(iindex,ic,jc,kstol) = itndex(iindex,ic,jc,kstol)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(bountt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_READ),  &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(dp)", OPS_READ))

            DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(iindex,ic,jc,kstol),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtzr(ic,jc,1) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtzr(ic,jc,1)
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
              + (fornow-rgspec(ispec)*strtzr(ic,jc,1))*yrhs(ispec,ic,jc,kstol)
          
        END DO
      END DO
      
    END DO

            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_READ))

        END IF
  
!       =======================================================================
  
    END IF
!   Z-DIRECTION RIGHT-HAND END

!   =========================================================================

END SUBROUTINE bountt
