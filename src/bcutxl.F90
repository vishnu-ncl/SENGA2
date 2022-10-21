SUBROUTINE bcutxl
    
    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

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
!KA      real(kind=dp)) BTIME
    real(kind=dp) :: fornow,argmnt,argval,realkx
    real(kind=dp) :: cosval,sinval,costht,sintht
    real(kind=dp) :: pcount
    integer :: ic,jc,kc
    integer :: iic,iim,kx,kxbase
    integer :: icproc,ncount,irproc,irtag
    integer :: rangexyz(6)
    real(kind=dp) :: init_val1, init_val2

!   BEGIN
!   =====

!   =========================================================================

!KA   THIS WAS MOVED TO BOUNDT & BOUNTT TO FIX INFLOW SCANNING LOCATION
!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
!KA      BTIME = ETIME + RKTIM(IRKSTP)

!   =========================================================================

!   CONSTANT U-VELOCITY
!   PARAMETER I1=1, R1=U-VELOCITY
    IF(nxlprm(1) == 1) THEN
        rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
        call ops_par_loop(bcut_kernel_xdir_const_uvel, "bcut_kernel_xdir_const_uvel", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_gbl(rxlprm(1), 1, "real(dp)", OPS_READ))
        
    END IF

!   =========================================================================

!   SINUSOIDAL U-VELOCITY
!   PARAMETER I1=2, R1=AMPLITUDE, R2=PERIOD
    IF(nxlprm(1) == 2) THEN
        fornow = two*pi/rxlprm(2)
        argmnt = fornow*btime
        init_val1 = rxlprm(1)*SIN(argmnt)
        init_val2 = fornow*rxlprm(1)*COS(argmnt)

        rangexyz = (/1,1,jstal,jstol,kstal,kstol/)
        call ops_par_loop(bcut_kernel_xdir_sinusoidal_uvel, "bcut_kernel_xdir_sinusoidal_uvel", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_dudtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_dvdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE),  &
                        ops_arg_dat(d_dwdtxl, 1, s3d_000_strid3d_yz, "real(dp)", OPS_WRITE), &
                        ops_arg_gbl(init_val1, 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(init_val2, 1, "real(dp)", OPS_READ))

    END IF

!   =========================================================================

!   TURBULENT VELOCITY FIELD
!   PARAMETER I1=3
    IF(nxlprm(1) == 3) THEN
!       INTERPOLATE STORED TURBULENT VELOCITY FIELD ONTO INLET PLANE
!       DO THE INTERPOLATION BY DFT: LOCAL-PROCESSOR CONTRIBUTION
  
!       -----------------------------------------------------------------------
  
!       UPDATE THE SCANNING PLANE LOCATION
        slocxl = elocxl - svelxl*btime
        IF(slocxl < zero) slocxl = xgdlen + slocxl
!KA         FIX INFLOW
!KA         IF(IRKSTP.EQ.NRKSTP)ELOCXL = SLOCXL
        IF(fupelc) elocxl = slocxl
  
!       INITIALISE THE PHASE ANGLE TERMS
        argmnt = tpovxg*slocxl
        costht = COS(argmnt)
        sintht = SIN(argmnt)
        kxbase = kminxl
  
!       ZERO THE LOCAL-PROCESSOR CONTRIBUTION TO THE DFT
        DO kc = kstal,kstol
            DO jc = jstal,jstol
      
                struxl(1,jc,kc) = zero
                strvxl(1,jc,kc) = zero
                strwxl(1,jc,kc) = zero
      
                dudtxl(1,jc,kc) = zero
                dvdtxl(1,jc,kc) = zero
                dwdtxl(1,jc,kc) = zero
      
            END DO
        END DO
  
!       -----------------------------------------------------------------------
  
!       SPECIAL CASE OF LEADING IMAGINARY TERM
        IF(fllixl) THEN
    
            kx = kxbase
            realkx = REAL(kx,kind=dp)
            argval = argmnt*realkx
            cosval = COS(argval)
            sinval = SIN(argval)
            iic = 1
    
            DO kc = kstal,kstol
                DO jc = jstal,jstol
        
                    struxl(1,jc,kc) = struxl(1,jc,kc) + ufxl(iic,jc,kc)*sinval
                    strvxl(1,jc,kc) = strvxl(1,jc,kc) + vfxl(iic,jc,kc)*sinval
                    strwxl(1,jc,kc) = strwxl(1,jc,kc) + wfxl(iic,jc,kc)*sinval
        
                    dudtxl(1,jc,kc) = dudtxl(1,jc,kc) - realkx*ufxl(iic,jc,kc)*cosval
                    dvdtxl(1,jc,kc) = dvdtxl(1,jc,kc) - realkx*vfxl(iic,jc,kc)*cosval
                    dwdtxl(1,jc,kc) = dwdtxl(1,jc,kc) - realkx*wfxl(iic,jc,kc)*cosval
        
                END DO
            END DO
    
            kxbase = kxbase + 1
    
        END IF
  
!       -----------------------------------------------------------------------
  
!       STANDARD LOCAL CONTRIBUTION
  
!       ZEROTH WAVENUMBER
        IF(kxbase == 0) THEN
            kx = kxbase
            realkx = REAL(kx,kind=dp)
            argval = argmnt*realkx
            cosval = COS(argval)
            sinval = SIN(argval)
            iim = 1
    
            DO kc = kstal,kstol
                DO jc = jstal,jstol
        
        
                    struxl(1,jc,kc) = struxl(1,jc,kc) + half*ufxl(iim,jc,kc)*cosval
                    strvxl(1,jc,kc) = strvxl(1,jc,kc) + half*vfxl(iim,jc,kc)*cosval
                    strwxl(1,jc,kc) = strwxl(1,jc,kc) + half*wfxl(iim,jc,kc)*cosval
        
                END DO
            END DO
    
            kxbase = kxbase + 1
    
        END IF
  
!       ALL OTHER WAVENUMBERS
        DO kc = kstal,kstol
            DO jc = jstal,jstol
      
                kx = kxbase
                realkx = REAL(kx,kind=dp)
                argval = argmnt*realkx
                cosval = COS(argval)
                sinval = SIN(argval)
      
                DO ic = istaxl,istoxl,2
        
                    iim = ic
                    iic = ic+1
        
                    struxl(1,jc,kc) = struxl(1,jc,kc) + ufxl(iim,jc,kc)*cosval  &
                                    + ufxl(iic,jc,kc)*sinval
                    strvxl(1,jc,kc) = strvxl(1,jc,kc) + vfxl(iim,jc,kc)*cosval  &
                                    + vfxl(iic,jc,kc)*sinval
                    strwxl(1,jc,kc) = strwxl(1,jc,kc) + wfxl(iim,jc,kc)*cosval  &
                                    + wfxl(iic,jc,kc)*sinval
        
                    dudtxl(1,jc,kc) = dudtxl(1,jc,kc) + realkx*(ufxl(iim,jc,kc)*sinval  &
                                    - ufxl(iic,jc,kc)*cosval)
                    dvdtxl(1,jc,kc) = dvdtxl(1,jc,kc) + realkx*(vfxl(iim,jc,kc)*sinval  &
                                    - vfxl(iic,jc,kc)*cosval)
                    dwdtxl(1,jc,kc) = dwdtxl(1,jc,kc) + realkx*(wfxl(iim,jc,kc)*sinval  &
                                    - wfxl(iic,jc,kc)*cosval)
        
                    kx = kx + 1
                    realkx = REAL(kx,kind=dp)
                    fornow = cosval
                    cosval = costht*cosval - sintht*sinval
                    sinval = sintht*fornow + costht*sinval
        
                END DO
      
            END DO
        END DO
  
!       -----------------------------------------------------------------------
  
!       SPECIAL CASE OF TRAILING REAL TERM
        IF(fltrxl) THEN
    
    kx = kxbase + istoxl/2
    realkx = REAL(kx,kind=dp)
    argval = argmnt*realkx
    cosval = COS(argval)
    sinval = SIN(argval)
    iim = istoxl + 1
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        struxl(1,jc,kc) = struxl(1,jc,kc) + ufxl(iim,jc,kc)*cosval
        strvxl(1,jc,kc) = strvxl(1,jc,kc) + vfxl(iim,jc,kc)*cosval
        strwxl(1,jc,kc) = strwxl(1,jc,kc) + wfxl(iim,jc,kc)*cosval
        
        dudtxl(1,jc,kc) = dudtxl(1,jc,kc) + realkx*ufxl(iim,jc,kc)*sinval
        dvdtxl(1,jc,kc) = dvdtxl(1,jc,kc) + realkx*vfxl(iim,jc,kc)*sinval
        dwdtxl(1,jc,kc) = dwdtxl(1,jc,kc) + realkx*wfxl(iim,jc,kc)*sinval
        
      END DO
    END DO
    
  END IF
  
!       -----------------------------------------------------------------------
  
!       PARALLEL TRANSFER
!       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
  IF(ixproc == 0)THEN
    
!         LEFTMOST PROCESSOR IN X
!         RECEIVE FROM ALL OTHER PROCESSORS IN X
    DO icproc = 1,nxprm1
      
      irproc = nprocx(icproc)
      irtag = irproc*nproc+iproc
      CALL p_recv(pcount,1,irproc,irtag)
      CALL p_recv(parray,nparay,irproc,irtag)
      
      ncount = 0
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          ncount = ncount + 1
          struxl(1,jc,kc) = struxl(1,jc,kc) + parray(ncount)
          ncount = ncount + 1
          strvxl(1,jc,kc) = strvxl(1,jc,kc) + parray(ncount)
          ncount = ncount + 1
          strwxl(1,jc,kc) = strwxl(1,jc,kc) + parray(ncount)
          ncount = ncount + 1
          dudtxl(1,jc,kc) = dudtxl(1,jc,kc) + parray(ncount)
          ncount = ncount + 1
          dvdtxl(1,jc,kc) = dvdtxl(1,jc,kc) + parray(ncount)
          ncount = ncount + 1
          dwdtxl(1,jc,kc) = dwdtxl(1,jc,kc) + parray(ncount)
          
        END DO
      END DO
      
    END DO
    
!         SCALING OF DFT
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             VELOCITIES
        struxl(1,jc,kc) = struxl(1,jc,kc)*scauxl
        strvxl(1,jc,kc) = strvxl(1,jc,kc)*scauxl
        strwxl(1,jc,kc) = strwxl(1,jc,kc)*scauxl
        
!             DERIVATIVES
        dudtxl(1,jc,kc) = dudtxl(1,jc,kc)*scduxl
        dvdtxl(1,jc,kc) = dvdtxl(1,jc,kc)*scduxl
        dwdtxl(1,jc,kc) = dwdtxl(1,jc,kc)*scduxl
        
!             ADD MEAN VELOCITY
        struxl(1,jc,kc) = struxl(1,jc,kc) + bvelxl
        
!             CONVERT SPATIAL TO TEMPORAL DERIVATIVES
        dudtxl(1,jc,kc) = dudtxl(1,jc,kc)*svelxl
        dvdtxl(1,jc,kc) = dvdtxl(1,jc,kc)*svelxl
        dwdtxl(1,jc,kc) = dwdtxl(1,jc,kc)*svelxl
        
      END DO
    END DO
    
  ELSE
    
!         NOT THE LEFTMOST PROCESSOR IN X
!         SEND TO LEFTMOST PROCESSOR IN X
    ncount = 0
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        ncount = ncount + 1
        parray(ncount) = struxl(1,jc,kc)
        ncount = ncount + 1
        parray(ncount) = strvxl(1,jc,kc)
        ncount = ncount + 1
        parray(ncount) = strwxl(1,jc,kc)
        ncount = ncount + 1
        parray(ncount) = dudtxl(1,jc,kc)
        ncount = ncount + 1
        parray(ncount) = dvdtxl(1,jc,kc)
        ncount = ncount + 1
        parray(ncount) = dwdtxl(1,jc,kc)
        
      END DO
    END DO
    
    pcount = REAL(ncount,kind=dp)
    irproc = nprocx(0)
    irtag = iproc*nproc+irproc
    CALL p_send(pcount,1,1,irproc,irtag)
    CALL p_send(parray,nparay,ncount,irproc,irtag)
    
  END IF
  
END IF

!   =========================================================================

END SUBROUTINE bcutxl
