SUBROUTINE bcutxl
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:14:46

!     *************************************************************************

!     BCUTXL
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     30-DEC-2003:  CREATED
!     04-JAN-2007:  RSC REVISE PARALLEL RECEIVES

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
!     AND THEIR TIME DERIVATIVES

!     X-DIRECTION LEFT-HAND END

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
!KA   FIX INFLOW BUG, BTIME IS DEFINED IN COM_SENGA2.H
!KA      REAL(KIND=dp) BTIME
REAL(KIND=dp) :: fornow,argmnt,argval,realkx
REAL(KIND=dp) :: cosval,sinval,costht,sintht
REAL(KIND=dp) :: pcount
INTEGER :: ic,jc,kc
INTEGER :: iic,iim,kx,kxbase
INTEGER :: icproc,ncount,irproc,irtag


!     BEGIN
!     =====

!     =========================================================================

!KA   THIS WAS MOVED TO BOUNDT & BOUNTT TO FIX INFLOW SCANNING LOCATION
!     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
!KA      BTIME = ETIME + RKTIM(IRKSTP)

!     =========================================================================

!     CONSTANT U-VELOCITY
!     PARAMETER I1=1, R1=U-VELOCITY
IF(nxlprm(1) == 1)THEN
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      struxl(1,jc,kc) = rxlprm(1)
!            FORNOW = REAL(JC)*DELTAY/(HALF*YGDLEN)
!            STRUXL(JC,KC) = RXLPRM(1)*TANH(FORNOW)
      strvxl(1,jc,kc) = zero
      strwxl(1,jc,kc) = zero
      
      dudtxl(jc,kc) = zero
      dvdtxl(jc,kc) = zero
      dwdtxl(jc,kc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     SINUSOIDAL U-VELOCITY
!     PARAMETER I1=2, R1=AMPLITUDE, R2=PERIOD
IF(nxlprm(1) == 2)THEN
  
  fornow = two*pi/rxlprm(2)
  argmnt = fornow*btime
  
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      struxl(1,jc,kc) = rxlprm(1)*SIN(argmnt)
      strvxl(1,jc,kc) = zero
      strwxl(1,jc,kc) = zero
      
      dudtxl(jc,kc) = fornow*rxlprm(1)*COS(argmnt)
      dvdtxl(jc,kc) = zero
      dwdtxl(jc,kc) = zero
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     TURBULENT VELOCITY FIELD
!     PARAMETER I1=3
IF(nxlprm(1) == 3)THEN
  
!       INTERPOLATE STORED TURBULENT VELOCITY FIELD ONTO INLET PLANE
!       DO THE INTERPOLATION BY DFT: LOCAL-PROCESSOR CONTRIBUTION
  
!       -----------------------------------------------------------------------
  
!       UPDATE THE SCANNING PLANE LOCATION
  slocxl = elocxl - svelxl*btime
  IF(slocxl < zero)slocxl = xgdlen + slocxl
!KA     FIX INFLOW
!KA        IF(IRKSTP.EQ.NRKSTP)ELOCXL = SLOCXL
  IF(fupelc)elocxl = slocxl
  
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
      
      dudtxl(jc,kc) = zero
      dvdtxl(jc,kc) = zero
      dwdtxl(jc,kc) = zero
      
    END DO
  END DO
  
!       -----------------------------------------------------------------------
  
!       SPECIAL CASE OF LEADING IMAGINARY TERM
  IF(fllixl)THEN
    
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
        
        dudtxl(jc,kc) = dudtxl(jc,kc) - realkx*ufxl(iic,jc,kc)*cosval
        dvdtxl(jc,kc) = dvdtxl(jc,kc) - realkx*vfxl(iic,jc,kc)*cosval
        dwdtxl(jc,kc) = dwdtxl(jc,kc) - realkx*wfxl(iic,jc,kc)*cosval
        
      END DO
    END DO
    
    kxbase = kxbase + 1
    
  END IF
  
!       -----------------------------------------------------------------------
  
!       STANDARD LOCAL CONTRIBUTION
  
!       ZEROTH WAVENUMBER
  IF(kxbase == 0)THEN
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        kx = kxbase
        realkx = REAL(kx,kind=dp)
        argval = argmnt*realkx
        cosval = COS(argval)
        sinval = SIN(argval)
        iim = 1
        
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
        
        dudtxl(jc,kc) = dudtxl(jc,kc) + realkx*(ufxl(iim,jc,kc)*sinval  &
            - ufxl(iic,jc,kc)*cosval)
        dvdtxl(jc,kc) = dvdtxl(jc,kc) + realkx*(vfxl(iim,jc,kc)*sinval  &
            - vfxl(iic,jc,kc)*cosval)
        dwdtxl(jc,kc) = dwdtxl(jc,kc) + realkx*(wfxl(iim,jc,kc)*sinval  &
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
  IF(fltrxl)THEN
    
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
        
        dudtxl(jc,kc) = dudtxl(jc,kc) + realkx*ufxl(iim,jc,kc)*sinval
        dvdtxl(jc,kc) = dvdtxl(jc,kc) + realkx*vfxl(iim,jc,kc)*sinval
        dwdtxl(jc,kc) = dwdtxl(jc,kc) + realkx*wfxl(iim,jc,kc)*sinval
        
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
          dudtxl(jc,kc) = dudtxl(jc,kc) + parray(ncount)
          ncount = ncount + 1
          dvdtxl(jc,kc) = dvdtxl(jc,kc) + parray(ncount)
          ncount = ncount + 1
          dwdtxl(jc,kc) = dwdtxl(jc,kc) + parray(ncount)
          
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
        dudtxl(jc,kc) = dudtxl(jc,kc)*scduxl
        dvdtxl(jc,kc) = dvdtxl(jc,kc)*scduxl
        dwdtxl(jc,kc) = dwdtxl(jc,kc)*scduxl
        
!             ADD MEAN VELOCITY
        struxl(1,jc,kc) = struxl(1,jc,kc) + bvelxl
        
!             CONVERT SPATIAL TO TEMPORAL DERIVATIVES
        dudtxl(jc,kc) = dudtxl(jc,kc)*svelxl
        dvdtxl(jc,kc) = dvdtxl(jc,kc)*svelxl
        dwdtxl(jc,kc) = dwdtxl(jc,kc)*svelxl
        
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
        parray(ncount) = dudtxl(jc,kc)
        ncount = ncount + 1
        parray(ncount) = dvdtxl(jc,kc)
        ncount = ncount + 1
        parray(ncount) = dwdtxl(jc,kc)
        
      END DO
    END DO
    
    pcount = REAL(ncount,kind=dp)
    irproc = nprocx(0)
    irtag = iproc*nproc+irproc
    CALL p_send(pcount,1,1,irproc,irtag)
    CALL p_send(parray,nparay,ncount,irproc,irtag)
    
  END IF
  
END IF

!     =========================================================================


RETURN
END SUBROUTINE bcutxl
