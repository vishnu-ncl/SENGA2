SUBROUTINE bountt
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:24:43

!     *************************************************************************

!     BOUNTT
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     29-SEP-2003:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     SYNCHRONISES THE TIME-DEPENDENT BOUNDARY CONDITIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
real(kind=dp) :: fornow
INTEGER :: jc,kc
INTEGER :: ispec
INTEGER :: iindex,ipower,icoef1,icoef2
INTEGER :: itint,icp


!     BEGIN
!     =====

!     =========================================================================

!     SYNCHRONISE AT CURRENT TIME STEP
!     --------------------------------
irkstp = 1
!KA   FIX INFLOW BC
btime  = tstep
fupelc = .true.
!     =========================================================================

!     X-DIRECTION LEFT-HAND END
!     -------------------------

!     GLOBAL BC SUPPORT
!     TURBULENT INFLOW VELOCITY FIELD
IF(fxltrb)CALL bcutxl

!     LOCAL BC SUPPORT
IF(fxlcnv)THEN
  
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
  
  IF(nsbcxl == nsbci2)THEN
    
!         INFLOW BC No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutxl
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttxl
    
!         SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        DO iindex = 1,nintmx
          itndex(istal,jc,kc,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          1000            CONTINUE
          IF(strtxl(jc,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 1000
            END IF
          END IF
!               END OF LOOP 1000
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(istal,jc,kc,iindex) = itndex(istal,jc,kc,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        urhs(istal,jc,kc) = drhs(istal,jc,kc)*struxl(jc,kc)
        vrhs(istal,jc,kc) = drhs(istal,jc,kc)*strvxl(jc,kc)
        wrhs(istal,jc,kc) = drhs(istal,jc,kc)*strwxl(jc,kc)
        
        urun(istal,jc,kc) = urhs(istal,jc,kc)
        vrun(istal,jc,kc) = vrhs(istal,jc,kc)
        wrun(istal,jc,kc) = wrhs(istal,jc,kc)
        
        uerr(istal,jc,kc) = zero
        verr(istal,jc,kc) = zero
        werr(istal,jc,kc) = zero
        
        erhs(istal,jc,kc) = half*(struxl(jc,kc)*struxl(jc,kc)  &
            + strvxl(jc,kc)*strvxl(jc,kc) + strwxl(jc,kc)*strwxl(jc,kc))
        erhs(istal,jc,kc) = drhs(istal,jc,kc)*erhs(istal,jc,kc)
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytxl
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          itint = 1 +MOD(itndex(istal,jc,kc,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtxl(jc,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtxl(jc,kc)
          
          yrhs(istal,jc,kc,ispec) = drhs(istal,jc,kc)*stryxl(jc,kc,ispec)
          
          yrun(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)
          
          yerr(istal,jc,kc,ispec) = zero
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
              + (fornow-rgspec(ispec)*strtxl(jc,kc))*yrhs(istal,jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erun(istal,jc,kc) = erhs(istal,jc,kc)
        
        eerr(istal,jc,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxl == nsbci3)THEN
    
!         INFLOW BC No 3
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!         SET DENSITY AND TIME DERIVATIVE
    CALL bcdtxl
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutxl
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = strdxl(jc,kc)
        urhs(istal,jc,kc) = strdxl(jc,kc)*struxl(jc,kc)
        vrhs(istal,jc,kc) = strdxl(jc,kc)*strvxl(jc,kc)
        wrhs(istal,jc,kc) = strdxl(jc,kc)*strwxl(jc,kc)
        
        drun(istal,jc,kc) = drhs(istal,jc,kc)
        urun(istal,jc,kc) = urhs(istal,jc,kc)
        vrun(istal,jc,kc) = vrhs(istal,jc,kc)
        wrun(istal,jc,kc) = wrhs(istal,jc,kc)
        
        derr(istal,jc,kc) = zero
        uerr(istal,jc,kc) = zero
        verr(istal,jc,kc) = zero
        werr(istal,jc,kc) = zero
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytxl
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istal,jc,kc,ispec) = strdxl(jc,kc)*stryxl(jc,kc,ispec)
          
          yrun(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)
          
          yerr(istal,jc,kc,ispec) = zero
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxl == nsbcw1)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ADIABATIC
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutxl
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        urhs(istal,jc,kc) = drhs(istal,jc,kc)*struxl(jc,kc)
        vrhs(istal,jc,kc) = drhs(istal,jc,kc)*strvxl(jc,kc)
        wrhs(istal,jc,kc) = drhs(istal,jc,kc)*strwxl(jc,kc)
        
        urun(istal,jc,kc) = urhs(istal,jc,kc)
        vrun(istal,jc,kc) = vrhs(istal,jc,kc)
        wrun(istal,jc,kc) = wrhs(istal,jc,kc)
        
        uerr(istal,jc,kc) = zero
        verr(istal,jc,kc) = zero
        werr(istal,jc,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxl == nsbcw2)THEN
    
!         WALL BC No 2
!         NO-SLIP WALL - ISOTHERMAL
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutxl
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttxl
    
!         SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        DO iindex = 1,nintmx
          itndex(istal,jc,kc,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          1100            CONTINUE
          IF(strtxl(jc,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 1100
            END IF
          END IF
!               END OF LOOP 1100
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(istal,jc,kc,iindex) = itndex(istal,jc,kc,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        urhs(istal,jc,kc) = drhs(istal,jc,kc)*struxl(jc,kc)
        vrhs(istal,jc,kc) = drhs(istal,jc,kc)*strvxl(jc,kc)
        wrhs(istal,jc,kc) = drhs(istal,jc,kc)*strwxl(jc,kc)
        
        urun(istal,jc,kc) = urhs(istal,jc,kc)
        vrun(istal,jc,kc) = vrhs(istal,jc,kc)
        wrun(istal,jc,kc) = wrhs(istal,jc,kc)
        
        uerr(istal,jc,kc) = zero
        verr(istal,jc,kc) = zero
        werr(istal,jc,kc) = zero
        
        erhs(istal,jc,kc) = half*(struxl(jc,kc)*struxl(jc,kc)  &
            + strvxl(jc,kc)*strvxl(jc,kc) + strwxl(jc,kc)*strwxl(jc,kc))
        erhs(istal,jc,kc) = drhs(istal,jc,kc)*erhs(istal,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          itint = 1 +MOD(itndex(istal,jc,kc,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtxl(jc,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtxl(jc,kc)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
              + (fornow-rgspec(ispec)*strtxl(jc,kc))*yrhs(istal,jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erun(istal,jc,kc) = erhs(istal,jc,kc)
        
        eerr(istal,jc,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
END IF
!     X-DIRECTION LEFT-HAND END

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     X-DIRECTION RIGHT-HAND END
!     --------------------------
IF(fxrcnv)THEN
  
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
  
  IF(nsbcxr == nsbci2)THEN
    
!         INFLOW BC No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutxr
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttxr
    
!         SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        DO iindex = 1,nintmx
          itndex(istol,jc,kc,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          1500            CONTINUE
          IF(strtxr(jc,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 1500
            END IF
          END IF
!               END OF LOOP 1500
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(istol,jc,kc,iindex) = itndex(istol,jc,kc,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        urhs(istol,jc,kc) = drhs(istol,jc,kc)*struxr(jc,kc)
        vrhs(istol,jc,kc) = drhs(istol,jc,kc)*strvxr(jc,kc)
        wrhs(istol,jc,kc) = drhs(istol,jc,kc)*strwxr(jc,kc)
        
        urun(istol,jc,kc) = urhs(istol,jc,kc)
        vrun(istol,jc,kc) = vrhs(istol,jc,kc)
        wrun(istol,jc,kc) = wrhs(istol,jc,kc)
        
        uerr(istol,jc,kc) = zero
        verr(istol,jc,kc) = zero
        werr(istol,jc,kc) = zero
        
        erhs(istol,jc,kc) = half*(struxr(jc,kc)*struxr(jc,kc)  &
            + strvxr(jc,kc)*strvxr(jc,kc) + strwxr(jc,kc)*strwxr(jc,kc))
        erhs(istol,jc,kc) = drhs(istol,jc,kc)*erhs(istol,jc,kc)
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytxr
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          itint = 1 +MOD(itndex(istol,jc,kc,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtxr(jc,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtxr(jc,kc)
          
          yrhs(istol,jc,kc,ispec) = drhs(istol,jc,kc)*stryxr(jc,kc,ispec)
          
          yrun(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)
          
          yerr(istol,jc,kc,ispec) = zero
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
              + (fornow-rgspec(ispec)*strtxr(jc,kc))*yrhs(istol,jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erun(istol,jc,kc) = erhs(istol,jc,kc)
        
        eerr(istol,jc,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxr == nsbci3)THEN
    
!         INFLOW BC No 3
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!         SET DENSITY AND TIME DERIVATIVE
    CALL bcdtxr
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutxr
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = strdxr(jc,kc)
        urhs(istol,jc,kc) = strdxr(jc,kc)*struxr(jc,kc)
        vrhs(istol,jc,kc) = strdxr(jc,kc)*strvxr(jc,kc)
        wrhs(istol,jc,kc) = strdxr(jc,kc)*strwxr(jc,kc)
        
        drun(istol,jc,kc) = drhs(istol,jc,kc)
        urun(istol,jc,kc) = urhs(istol,jc,kc)
        vrun(istol,jc,kc) = vrhs(istol,jc,kc)
        wrun(istol,jc,kc) = wrhs(istol,jc,kc)
        
        derr(istol,jc,kc) = zero
        uerr(istol,jc,kc) = zero
        verr(istol,jc,kc) = zero
        werr(istol,jc,kc) = zero
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytxr
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istol,jc,kc,ispec) = strdxr(jc,kc)*stryxr(jc,kc,ispec)
          
          yrun(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)
          
          yerr(istol,jc,kc,ispec) = zero
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxr == nsbcw1)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ADIABATIC
!         *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutxr
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        urhs(istol,jc,kc) = drhs(istol,jc,kc)*struxr(jc,kc)
        vrhs(istol,jc,kc) = drhs(istol,jc,kc)*strvxr(jc,kc)
        wrhs(istol,jc,kc) = drhs(istol,jc,kc)*strwxr(jc,kc)
        
        urun(istol,jc,kc) = urhs(istol,jc,kc)
        vrun(istol,jc,kc) = vrhs(istol,jc,kc)
        wrun(istol,jc,kc) = wrhs(istol,jc,kc)
        
        uerr(istol,jc,kc) = zero
        verr(istol,jc,kc) = zero
        werr(istol,jc,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxr == nsbcw2)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ISOTHERMAL
!         *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutxr
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttxr
    
!         SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        DO iindex = 1,nintmx
          itndex(istol,jc,kc,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          1600            CONTINUE
          IF(strtxr(jc,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 1600
            END IF
          END IF
!               END OF LOOP 1600
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(istol,jc,kc,iindex) = itndex(istol,jc,kc,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        urhs(istol,jc,kc) = drhs(istol,jc,kc)*struxr(jc,kc)
        vrhs(istol,jc,kc) = drhs(istol,jc,kc)*strvxr(jc,kc)
        wrhs(istol,jc,kc) = drhs(istol,jc,kc)*strwxr(jc,kc)
        
        urun(istol,jc,kc) = urhs(istol,jc,kc)
        vrun(istol,jc,kc) = vrhs(istol,jc,kc)
        wrun(istol,jc,kc) = wrhs(istol,jc,kc)
        
        uerr(istol,jc,kc) = zero
        verr(istol,jc,kc) = zero
        werr(istol,jc,kc) = zero
        
        erhs(istol,jc,kc) = half*(struxr(jc,kc)*struxr(jc,kc)  &
            + strvxr(jc,kc)*strvxr(jc,kc) + strwxr(jc,kc)*strwxr(jc,kc))
        erhs(istol,jc,kc) = drhs(istol,jc,kc)*erhs(istol,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          itint = 1 +MOD(itndex(istol,jc,kc,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtxr(jc,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtxr(jc,kc)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
              + (fornow-rgspec(ispec)*strtxr(jc,kc))*yrhs(istol,jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erun(istol,jc,kc) = erhs(istol,jc,kc)
        
        eerr(istol,jc,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
END IF
!     X-DIRECTION RIGHT-HAND END

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     Y-DIRECTION LEFT-HAND END
!     -------------------------

!     GLOBAL BC SUPPORT
!     TURBULENT INFLOW VELOCITY FIELD
IF(fyltrb)CALL bcutyl

!     LOCAL BC SUPPORT
IF(fylcnv)THEN
  
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
  
  IF(nsbcyl == nsbci2)THEN
    
!         INFLOW BC No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutyl
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttyl
    
!         SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(ic,jstal,kc,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          2000            CONTINUE
          IF(strtyl(ic,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 2000
            END IF
          END IF
!               END OF LOOP 2000
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(ic,jstal,kc,iindex) = itndex(ic,jstal,kc,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        urhs(ic,jstal,kc) = drhs(ic,jstal,kc)*struyl(ic,kc)
        vrhs(ic,jstal,kc) = drhs(ic,jstal,kc)*strvyl(ic,kc)
        wrhs(ic,jstal,kc) = drhs(ic,jstal,kc)*strwyl(ic,kc)
        
        urun(ic,jstal,kc) = urhs(ic,jstal,kc)
        vrun(ic,jstal,kc) = vrhs(ic,jstal,kc)
        wrun(ic,jstal,kc) = wrhs(ic,jstal,kc)
        
        uerr(ic,jstal,kc) = zero
        verr(ic,jstal,kc) = zero
        werr(ic,jstal,kc) = zero
        
        erhs(ic,jstal,kc) = half*(struyl(ic,kc)*struyl(ic,kc)  &
            + strvyl(ic,kc)*strvyl(ic,kc) + strwyl(ic,kc)*strwyl(ic,kc))
        erhs(ic,jstal,kc) = drhs(ic,jstal,kc)*erhs(ic,jstal,kc)
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytyl
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(ic,jstal,kc,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtyl(ic,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtyl(ic,kc)
          
          yrhs(ic,jstal,kc,ispec) = drhs(ic,jstal,kc)*stryyl(ic,kc,ispec)
          
          yrun(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)
          
          yerr(ic,jstal,kc,ispec) = zero
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
              + (fornow-rgspec(ispec)*strtyl(ic,kc))*yrhs(ic,jstal,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erun(ic,jstal,kc) = erhs(ic,jstal,kc)
        
        eerr(ic,jstal,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyl == nsbci3)THEN
    
!         INFLOW BC No 3
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!         SET DENSITY AND TIME DERIVATIVE
    CALL bcdtyl
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutyl
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = strdyl(ic,kc)
        urhs(ic,jstal,kc) = strdyl(ic,kc)*struyl(ic,kc)
        vrhs(ic,jstal,kc) = strdyl(ic,kc)*strvyl(ic,kc)
        wrhs(ic,jstal,kc) = strdyl(ic,kc)*strwyl(ic,kc)
        
        drun(ic,jstal,kc) = drhs(ic,jstal,kc)
        urun(ic,jstal,kc) = urhs(ic,jstal,kc)
        vrun(ic,jstal,kc) = vrhs(ic,jstal,kc)
        wrun(ic,jstal,kc) = wrhs(ic,jstal,kc)
        
        derr(ic,jstal,kc) = zero
        uerr(ic,jstal,kc) = zero
        verr(ic,jstal,kc) = zero
        werr(ic,jstal,kc) = zero
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytyl
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstal,kc,ispec) = strdyl(ic,kc)*stryyl(ic,kc,ispec)
          
          yrun(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)
          
          yerr(ic,jstal,kc,ispec) = zero
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyl == nsbcw1)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ADIABATIC
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutyl
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        urhs(ic,jstal,kc) = drhs(ic,jstal,kc)*struyl(ic,kc)
        vrhs(ic,jstal,kc) = drhs(ic,jstal,kc)*strvyl(ic,kc)
        wrhs(ic,jstal,kc) = drhs(ic,jstal,kc)*strwyl(ic,kc)
        
        urun(ic,jstal,kc) = urhs(ic,jstal,kc)
        vrun(ic,jstal,kc) = vrhs(ic,jstal,kc)
        wrun(ic,jstal,kc) = wrhs(ic,jstal,kc)
        
        uerr(ic,jstal,kc) = zero
        verr(ic,jstal,kc) = zero
        werr(ic,jstal,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyl == nsbcw2)THEN
    
!         WALL BC No 2
!         NO-SLIP WALL - ISOTHERMAL
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutyl
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttyl
    
!         SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(ic,jstal,kc,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          2100            CONTINUE
          IF(strtyl(ic,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 2100
            END IF
          END IF
!               END OF LOOP 2100
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(ic,jstal,kc,iindex) = itndex(ic,jstal,kc,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        urhs(ic,jstal,kc) = drhs(ic,jstal,kc)*struyl(ic,kc)
        vrhs(ic,jstal,kc) = drhs(ic,jstal,kc)*strvyl(ic,kc)
        wrhs(ic,jstal,kc) = drhs(ic,jstal,kc)*strwyl(ic,kc)
        
        urun(ic,jstal,kc) = urhs(ic,jstal,kc)
        vrun(ic,jstal,kc) = vrhs(ic,jstal,kc)
        wrun(ic,jstal,kc) = wrhs(ic,jstal,kc)
        
        uerr(ic,jstal,kc) = zero
        verr(ic,jstal,kc) = zero
        werr(ic,jstal,kc) = zero
        
        erhs(ic,jstal,kc) = half*(struyl(ic,kc)*struyl(ic,kc)  &
            + strvyl(ic,kc)*strvyl(ic,kc) + strwyl(ic,kc)*strwyl(ic,kc))
        erhs(ic,jstal,kc) = drhs(ic,jstal,kc)*erhs(ic,jstal,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(ic,jstal,kc,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtyl(ic,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtyl(ic,kc)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
              + (fornow-rgspec(ispec)*strtyl(ic,kc))*yrhs(ic,jstal,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erun(ic,jstal,kc) = erhs(ic,jstal,kc)
        
        eerr(ic,jstal,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
END IF
!     Y-DIRECTION LEFT-HAND END

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     Y-DIRECTION RIGHT-HAND END
!     --------------------------

!     GLOBAL BC SUPPORT
!     TURBULENT INFLOW VELOCITY FIELD
IF(fyrtrb)CALL bcutyr

!     LOCAL BC SUPPORT
IF(fyrcnv)THEN
  
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
  
  IF(nsbcyr == nsbci2)THEN
    
!         INFLOW BC No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutyr
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttyr
    
!         SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(ic,jstol,kc,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          2500            CONTINUE
          IF(strtyr(ic,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 2500
            END IF
          END IF
!               END OF LOOP 2500
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(ic,jstol,kc,iindex) = itndex(ic,jstol,kc,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        urhs(ic,jstol,kc) = drhs(ic,jstol,kc)*struyr(ic,kc)
        vrhs(ic,jstol,kc) = drhs(ic,jstol,kc)*strvyr(ic,kc)
        wrhs(ic,jstol,kc) = drhs(ic,jstol,kc)*strwyr(ic,kc)
        
        urun(ic,jstol,kc) = urhs(ic,jstol,kc)
        vrun(ic,jstol,kc) = vrhs(ic,jstol,kc)
        wrun(ic,jstol,kc) = wrhs(ic,jstol,kc)
        
        uerr(ic,jstol,kc) = zero
        verr(ic,jstol,kc) = zero
        werr(ic,jstol,kc) = zero
        
        erhs(ic,jstol,kc) = half*(struyr(ic,kc)*struyr(ic,kc)  &
            + strvyr(ic,kc)*strvyr(ic,kc) + strwyr(ic,kc)*strwyr(ic,kc))
        erhs(ic,jstol,kc) = drhs(ic,jstol,kc)*erhs(ic,jstol,kc)
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytyr
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(ic,jstol,kc,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtyr(ic,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtyr(ic,kc)
          
          yrhs(ic,jstol,kc,ispec) = drhs(ic,jstol,kc)*stryyr(ic,kc,ispec)
          
          yrun(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)
          
          yerr(ic,jstol,kc,ispec) = zero
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
              + (fornow-rgspec(ispec)*strtyr(ic,kc))*yrhs(ic,jstol,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erun(ic,jstol,kc) = erhs(ic,jstol,kc)
        
        eerr(ic,jstol,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyr == nsbci3)THEN
    
!         INFLOW BC No 3
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!         SET DENSITY AND TIME DERIVATIVE
    CALL bcdtyr
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutyr
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = strdyr(ic,kc)
        urhs(ic,jstol,kc) = strdyr(ic,kc)*struyr(ic,kc)
        vrhs(ic,jstol,kc) = strdyr(ic,kc)*strvyr(ic,kc)
        wrhs(ic,jstol,kc) = strdyr(ic,kc)*strwyr(ic,kc)
        
        drun(ic,jstol,kc) = drhs(ic,jstol,kc)
        urun(ic,jstol,kc) = urhs(ic,jstol,kc)
        vrun(ic,jstol,kc) = vrhs(ic,jstol,kc)
        wrun(ic,jstol,kc) = wrhs(ic,jstol,kc)
        
        derr(ic,jstol,kc) = zero
        uerr(ic,jstol,kc) = zero
        verr(ic,jstol,kc) = zero
        werr(ic,jstol,kc) = zero
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytyr
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstol,kc,ispec) = strdyr(ic,kc)*stryyr(ic,kc,ispec)
          
          yrun(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)
          
          yerr(ic,jstol,kc,ispec) = zero
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyr == nsbcw1)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ADIABATIC
!         *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutyr
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        urhs(ic,jstol,kc) = drhs(ic,jstol,kc)*struyr(ic,kc)
        vrhs(ic,jstol,kc) = drhs(ic,jstol,kc)*strvyr(ic,kc)
        wrhs(ic,jstol,kc) = drhs(ic,jstol,kc)*strwyr(ic,kc)
        
        urun(ic,jstol,kc) = urhs(ic,jstol,kc)
        vrun(ic,jstol,kc) = vrhs(ic,jstol,kc)
        wrun(ic,jstol,kc) = wrhs(ic,jstol,kc)
        
        uerr(ic,jstol,kc) = zero
        verr(ic,jstol,kc) = zero
        werr(ic,jstol,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyr == nsbcw2)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ISOTHERMAL
!         *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutyr
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttyr
    
!         SET TEMPERATURE INTERVAL INDEX
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(ic,jstol,kc,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          2600            CONTINUE
          IF(strtyr(ic,kc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 2600
            END IF
          END IF
!               END OF LOOP 2600
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(ic,jstol,kc,iindex) = itndex(ic,jstol,kc,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        urhs(ic,jstol,kc) = drhs(ic,jstol,kc)*struyr(ic,kc)
        vrhs(ic,jstol,kc) = drhs(ic,jstol,kc)*strvyr(ic,kc)
        wrhs(ic,jstol,kc) = drhs(ic,jstol,kc)*strwyr(ic,kc)
        
        urun(ic,jstol,kc) = urhs(ic,jstol,kc)
        vrun(ic,jstol,kc) = vrhs(ic,jstol,kc)
        wrun(ic,jstol,kc) = wrhs(ic,jstol,kc)
        
        uerr(ic,jstol,kc) = zero
        verr(ic,jstol,kc) = zero
        werr(ic,jstol,kc) = zero
        
        erhs(ic,jstol,kc) = half*(struyr(ic,kc)*struyr(ic,kc)  &
            + strvyr(ic,kc)*strvyr(ic,kc) + strwyr(ic,kc)*strwyr(ic,kc))
        erhs(ic,jstol,kc) = drhs(ic,jstol,kc)*erhs(ic,jstol,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(ic,jstol,kc,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtyr(ic,kc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtyr(ic,kc)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
              + (fornow-rgspec(ispec)*strtyr(ic,kc))*yrhs(ic,jstol,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erun(ic,jstol,kc) = erhs(ic,jstol,kc)
        
        eerr(ic,jstol,kc) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
END IF
!     Y-DIRECTION RIGHT-HAND END

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     Z-DIRECTION LEFT-HAND END
!     -------------------------

!     GLOBAL BC SUPPORT
!     TURBULENT INFLOW VELOCITY FIELD
IF(fzltrb)CALL bcutzl

!     LOCAL BC SUPPORT
IF(fzlcnv)THEN
  
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
  
  IF(nsbczl == nsbci2)THEN
    
!         INFLOW BC No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutzl
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttzl
    
!         SET TEMPERATURE INTERVAL INDEX
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(ic,jc,kstal,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          3000            CONTINUE
          IF(strtzl(ic,jc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 3000
            END IF
          END IF
!               END OF LOOP 3000
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(ic,jc,kstal,iindex) = itndex(ic,jc,kstal,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        urhs(ic,jc,kstal) = drhs(ic,jc,kstal)*struzl(ic,jc)
        vrhs(ic,jc,kstal) = drhs(ic,jc,kstal)*strvzl(ic,jc)
        wrhs(ic,jc,kstal) = drhs(ic,jc,kstal)*strwzl(ic,jc)
        
        urun(ic,jc,kstal) = urhs(ic,jc,kstal)
        vrun(ic,jc,kstal) = vrhs(ic,jc,kstal)
        wrun(ic,jc,kstal) = wrhs(ic,jc,kstal)
        
        uerr(ic,jc,kstal) = zero
        verr(ic,jc,kstal) = zero
        werr(ic,jc,kstal) = zero
        
        erhs(ic,jc,kstal) = half*(struzl(ic,jc)*struzl(ic,jc)  &
            + strvzl(ic,jc)*strvzl(ic,jc) + strwzl(ic,jc)*strwzl(ic,jc))
        erhs(ic,jc,kstal) = drhs(ic,jc,kstal)*erhs(ic,jc,kstal)
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytzl
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(ic,jc,kstal,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtzl(ic,jc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtzl(ic,jc)
          
          yrhs(ic,jc,kstal,ispec) = drhs(ic,jc,kstal)*stryzl(ic,jc,ispec)
          
          yrun(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)
          
          yerr(ic,jc,kstal,ispec) = zero
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
              + (fornow-rgspec(ispec)*strtzl(ic,jc))*yrhs(ic,jc,kstal,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erun(ic,jc,kstal) = erhs(ic,jc,kstal)
        
        eerr(ic,jc,kstal) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczl == nsbci3)THEN
    
!         INFLOW BC No 3
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!         SET DENSITY AND TIME DERIVATIVE
    CALL bcdtzl
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutzl
    
!         CONSERVATIVE VARIABLES
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = strdzl(ic,jc)
        urhs(ic,jc,kstal) = strdzl(ic,jc)*struzl(ic,jc)
        vrhs(ic,jc,kstal) = strdzl(ic,jc)*strvzl(ic,jc)
        wrhs(ic,jc,kstal) = strdzl(ic,jc)*strwzl(ic,jc)
        
        drun(ic,jc,kstal) = drhs(ic,jc,kstal)
        urun(ic,jc,kstal) = urhs(ic,jc,kstal)
        vrun(ic,jc,kstal) = vrhs(ic,jc,kstal)
        wrun(ic,jc,kstal) = wrhs(ic,jc,kstal)
        
        derr(ic,jc,kstal) = zero
        uerr(ic,jc,kstal) = zero
        verr(ic,jc,kstal) = zero
        werr(ic,jc,kstal) = zero
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytzl
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstal,ispec) = strdzl(ic,jc)*stryzl(ic,jc,ispec)
          
          yrun(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)
          
          yerr(ic,jc,kstal,ispec) = zero
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczl == nsbcw1)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ADIABATIC
!         *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutzl
    
!         CONSERVATIVE VARIABLES
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        urhs(ic,jc,kstal) = drhs(ic,jc,kstal)*struzl(ic,jc)
        vrhs(ic,jc,kstal) = drhs(ic,jc,kstal)*strvzl(ic,jc)
        wrhs(ic,jc,kstal) = drhs(ic,jc,kstal)*strwzl(ic,jc)
        
        urun(ic,jc,kstal) = urhs(ic,jc,kstal)
        vrun(ic,jc,kstal) = vrhs(ic,jc,kstal)
        wrun(ic,jc,kstal) = wrhs(ic,jc,kstal)
        
        uerr(ic,jc,kstal) = zero
        verr(ic,jc,kstal) = zero
        werr(ic,jc,kstal) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczl == nsbcw2)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ISOTHERMAL
!         *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutzl
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttzl
    
!         SET TEMPERATURE INTERVAL INDEX
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(ic,jc,kstal,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          3100            CONTINUE
          IF(strtzl(ic,jc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 3100
            END IF
          END IF
!               END OF LOOP 3100
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(ic,jc,kstal,iindex) = itndex(ic,jc,kstal,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        urhs(ic,jc,kstal) = drhs(ic,jc,kstal)*struzl(ic,jc)
        vrhs(ic,jc,kstal) = drhs(ic,jc,kstal)*strvzl(ic,jc)
        wrhs(ic,jc,kstal) = drhs(ic,jc,kstal)*strwzl(ic,jc)
        
        urun(ic,jc,kstal) = urhs(ic,jc,kstal)
        vrun(ic,jc,kstal) = vrhs(ic,jc,kstal)
        wrun(ic,jc,kstal) = wrhs(ic,jc,kstal)
        
        uerr(ic,jc,kstal) = zero
        verr(ic,jc,kstal) = zero
        werr(ic,jc,kstal) = zero
        
        erhs(ic,jc,kstal) = half*(struzl(ic,jc)*struzl(ic,jc)  &
            + strvzl(ic,jc)*strvzl(ic,jc) + strwzl(ic,jc)*strwzl(ic,jc))
        erhs(ic,jc,kstal) = drhs(ic,jc,kstal)*erhs(ic,jc,kstal)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(ic,jc,kstal,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtzl(ic,jc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtzl(ic,jc)
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
              + (fornow-rgspec(ispec)*strtzl(ic,jc))*yrhs(ic,jc,kstal,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erun(ic,jc,kstal) = erhs(ic,jc,kstal)
        
        eerr(ic,jc,kstal) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
END IF
!     Z-DIRECTION LEFT-HAND END

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================

!     Z-DIRECTION RIGHT-HAND END
!     --------------------------

!     GLOBAL BC SUPPORT
!     TURBULENT INFLOW VELOCITY FIELD
IF(fzrtrb)CALL bcutzr

!     LOCAL BC SUPPORT
IF(fzrcnv)THEN
  
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
  
  IF(nsbczr == nsbci2)THEN
    
!         INFLOW BC No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutzr
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttzr
    
!         SET TEMPERATURE INTERVAL INDEX
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(ic,jc,kstol,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          3500            CONTINUE
          IF(strtzr(ic,jc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 3500
            END IF
          END IF
!               END OF LOOP 3500
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(ic,jc,kstol,iindex) = itndex(ic,jc,kstol,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        urhs(ic,jc,kstol) = drhs(ic,jc,kstol)*struzr(ic,jc)
        vrhs(ic,jc,kstol) = drhs(ic,jc,kstol)*strvzr(ic,jc)
        wrhs(ic,jc,kstol) = drhs(ic,jc,kstol)*strwzr(ic,jc)
        
        urun(ic,jc,kstol) = urhs(ic,jc,kstol)
        vrun(ic,jc,kstol) = vrhs(ic,jc,kstol)
        wrun(ic,jc,kstol) = wrhs(ic,jc,kstol)
        
        uerr(ic,jc,kstol) = zero
        verr(ic,jc,kstol) = zero
        werr(ic,jc,kstol) = zero
        
        erhs(ic,jc,kstol) = half*(struzr(ic,jc)*struzr(ic,jc)  &
            + strvzr(ic,jc)*strvzr(ic,jc) + strwzr(ic,jc)*strwzr(ic,jc))
        erhs(ic,jc,kstol) = drhs(ic,jc,kstol)*erhs(ic,jc,kstol)
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytzr
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(ic,jc,kstol,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtzr(ic,jc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtzr(ic,jc)
          
          yrhs(ic,jc,kstol,ispec) = drhs(ic,jc,kstol)*stryzr(ic,jc,ispec)
          
          yrun(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)
          
          yerr(ic,jc,kstol,ispec) = zero
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
              + (fornow-rgspec(ispec)*strtzr(ic,jc))*yrhs(ic,jc,kstol,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erun(ic,jc,kstol) = erhs(ic,jc,kstol)
        
        eerr(ic,jc,kstol) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczr == nsbci3)THEN
    
!         INFLOW BC No 3
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY
    
!         SET DENSITY AND TIME DERIVATIVE
    CALL bcdtzr
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutzr
    
!         CONSERVATIVE VARIABLES
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = strdzr(ic,jc)
        urhs(ic,jc,kstol) = strdzr(ic,jc)*struzr(ic,jc)
        vrhs(ic,jc,kstol) = strdzr(ic,jc)*strvzr(ic,jc)
        wrhs(ic,jc,kstol) = strdzr(ic,jc)*strwzr(ic,jc)
        
        drun(ic,jc,kstol) = drhs(ic,jc,kstol)
        urun(ic,jc,kstol) = urhs(ic,jc,kstol)
        vrun(ic,jc,kstol) = vrhs(ic,jc,kstol)
        wrun(ic,jc,kstol) = wrhs(ic,jc,kstol)
        
        derr(ic,jc,kstol) = zero
        uerr(ic,jc,kstol) = zero
        verr(ic,jc,kstol) = zero
        werr(ic,jc,kstol) = zero
        
      END DO
    END DO
    
!         SET MASS FRACTIONS AND TIME DERIVATIVES
    CALL bcytzr
    
!         CONSERVATIVE VARIABLES
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstol,ispec) = strdzr(ic,jc)*stryzr(ic,jc,ispec)
          
          yrun(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)
          
          yerr(ic,jc,kstol,ispec) = zero
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczr == nsbcw1)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ADIABATIC
!         *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutzr
    
!         CONSERVATIVE VARIABLES
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        urhs(ic,jc,kstol) = drhs(ic,jc,kstol)*struzr(ic,jc)
        vrhs(ic,jc,kstol) = drhs(ic,jc,kstol)*strvzr(ic,jc)
        wrhs(ic,jc,kstol) = drhs(ic,jc,kstol)*strwzr(ic,jc)
        
        urun(ic,jc,kstol) = urhs(ic,jc,kstol)
        vrun(ic,jc,kstol) = vrhs(ic,jc,kstol)
        wrun(ic,jc,kstol) = wrhs(ic,jc,kstol)
        
        uerr(ic,jc,kstol) = zero
        verr(ic,jc,kstol) = zero
        werr(ic,jc,kstol) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczr == nsbcw2)THEN
    
!         WALL BC No 1
!         NO-SLIP WALL - ISOTHERMAL
!         *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***
    
!         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
    CALL bcutzr
    
!         SET TEMPERATURE AND TIME DERIVATIVE
    CALL bcttzr
    
!         SET TEMPERATURE INTERVAL INDEX
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        DO iindex = 1,nintmx
          itndex(ic,jc,kstol,iindex) = 0
        END DO
        
        DO ispec = 1,nspec
          
          itint = 1
          3600            CONTINUE
          IF(strtzr(ic,jc) > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 3600
            END IF
          END IF
!               END OF LOOP 3600
          
!               SET THE TEMPERATURE INTERVAL INDEX
          iindex = 1 + (ispec-1)/nspimx
          ipower = ispec - (iindex-1)*nspimx - 1
          itndex(ic,jc,kstol,iindex) = itndex(ic,jc,kstol,iindex)  &
              +(itint-1)*ntbase**ipower
          
        END DO
        
      END DO
    END DO
    
!         CONSERVATIVE VARIABLES
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        urhs(ic,jc,kstol) = drhs(ic,jc,kstol)*struzr(ic,jc)
        vrhs(ic,jc,kstol) = drhs(ic,jc,kstol)*strvzr(ic,jc)
        wrhs(ic,jc,kstol) = drhs(ic,jc,kstol)*strwzr(ic,jc)
        
        urun(ic,jc,kstol) = urhs(ic,jc,kstol)
        vrun(ic,jc,kstol) = vrhs(ic,jc,kstol)
        wrun(ic,jc,kstol) = wrhs(ic,jc,kstol)
        
        uerr(ic,jc,kstol) = zero
        verr(ic,jc,kstol) = zero
        werr(ic,jc,kstol) = zero
        
        erhs(ic,jc,kstol) = half*(struzr(ic,jc)*struzr(ic,jc)  &
            + strvzr(ic,jc)*strvzr(ic,jc) + strwzr(ic,jc)*strwzr(ic,jc))
        erhs(ic,jc,kstol) = drhs(ic,jc,kstol)*erhs(ic,jc,kstol)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
!           TEMPERATURE INTERVAL INDEXING
      iindex = 1 + (ispec-1)/nspimx
      ipower = ispec - (iindex-1)*nspimx - 1
      icoef2 = ntbase**ipower
      icoef1 = icoef2*ntbase
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          itint = 1 +MOD(itndex(ic,jc,kstol,iindex),icoef1)/icoef2
          fornow = amasch(ncpoly(itint,ispec),itint,ispec)
          DO icp = ncpom1(itint,ispec),1,-1
            fornow = fornow*strtzr(ic,jc) + amasch(icp,itint,ispec)
          END DO
          fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtzr(ic,jc)
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
              + (fornow-rgspec(ispec)*strtzr(ic,jc))*yrhs(ic,jc,kstol,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erun(ic,jc,kstol) = erhs(ic,jc,kstol)
        
        eerr(ic,jc,kstol) = zero
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
END IF
!     Z-DIRECTION RIGHT-HAND END

!     =========================================================================


RETURN
END SUBROUTINE bountt
