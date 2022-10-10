SUBROUTINE bounds
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:24:36

!     *************************************************************************

!     BOUNDS
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     13-JUL-2003:  RSC MODIFIED FOR SENGA2
!     08-AUG-2012:  RSC EVALUATE ALL SPECIES
!     26-OCT-2013:  RSC ACTIVATE ALL BCS ON ALL SIDES

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     COMPUTES CHARACTERISTIC BOUNDARY CONDITIONS FOR ALL PDES

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
INTEGER :: ic,jc,kc
INTEGER :: ispec


!     BEGIN
!     =====

!     =========================================================================
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!     =========================================================================


!     X-DIRECTION LEFT-HAND END
!     -------------------------
IF(fxlcnv)THEN
  
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
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        strhxl(jc,kc,ispec) = strhxl(jc,kc,ispec)  &
            - strgxl(1,jc,kc)*strtxl(1,jc,kc)*rgspec(ispec)/strrxl(1,jc,kc)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      gam1xl(jc,kc) = strgxl(1,jc,kc) - strrxl(1,jc,kc)
      strexl(1,jc,kc) = strexl(1,jc,kc) - gam1xl(jc,kc)*strtxl(1,jc,kc)
      
      gam1xl(jc,kc) = strrxl(1,jc,kc)/gam1xl(jc,kc)
      ovgmxl(jc,kc) = one/gam1xl(jc,kc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      fornow = strgxl(1,jc,kc)*gam1xl(jc,kc)*strtxl(1,jc,kc)
      acouxl(jc,kc) = SQRT(fornow)
      ova2xl(jc,kc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbcxl == nsbco1)THEN
    
!         OUTFLOW BC No 1
!         SUBSONIC NON-REFLECTING OUTFLOW
!         WITH OPTION TO SET PRESSURE AT INFINITY
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxl(jc,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          sorpxl(jc,kc) = sorpxl(jc,kc)  &
              + strhxl(jc,kc,ispec)*ratexl(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxl(jc,kc) = -sorpxl(jc,kc)*gam1xl(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L5X AS REQUIRED
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L5X
        bcl5xl(1,jc,kc) = half*(struxl(1,jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)+strdxl(1,jc,kc)*acouxl(jc,kc)*bcl1xl(1,jc,kc))
        
!             SUBTRACT FROM NEW VALUE OF L5X
        bcl5xl(1,jc,kc)= half*sorpxl(jc,kc)  &
            + cobcxl*acouxl(jc,kc)*(strpxl(1,jc,kc)-pinfxl) - bcl5xl(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl5xl(1,jc,kc)*ova2xl(jc,kc)
        
        urhs(istal,jc,kc) = urhs(istal,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(jc,kc)*(struxl(1,jc,kc)+acouxl(jc,kc))
        
        vrhs(istal,jc,kc) = vrhs(istal,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(jc,kc)*strvxl(1,jc,kc)
        
        wrhs(istal,jc,kc) = wrhs(istal,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(jc,kc)*strwxl(1,jc,kc)
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
            - bcl5xl(1,jc,kc)*(ova2xl(jc,kc)*strexl(1,jc,kc)  &
            + struxl(1,jc,kc)/acouxl(jc,kc) + ovgmxl(jc,kc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)  &
              - bcl5xl(1,jc,kc)*ova2xl(jc,kc)*stryxl(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
  IF(nsbcxl == nsbci1)THEN
    
!         INFLOW BC No 1
!         SUBSONIC NON-REFLECTING LAMINAR INFLOW
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxl(jc,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          sorpxl(jc,kc) = sorpxl(jc,kc)  &
              + strhxl(jc,kc,ispec)*ratexl(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxl(jc,kc) = -sorpxl(jc,kc)*gam1xl(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L2X-L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxl(1,jc,kc)*acouxl(jc,kc)*bcl1xl(1,jc,kc)
        bcl2xl(1,jc,kc) = struxl(1,jc,kc)  &
            *(bcl2xl(1,jc,kc)-bcl5xl(1,jc,kc)*ova2xl(jc,kc))
        bcl3xl(1,jc,kc) = struxl(1,jc,kc)*bcl3xl(1,jc,kc)
        bcl4xl(1,jc,kc) = struxl(1,jc,kc)*bcl4xl(1,jc,kc)
        bcl5xl(1,jc,kc) = half*(struxl(1,jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
!             L1X UNCHANGED
        bcl2xl(1,jc,kc) = -bcl2xl(1,jc,kc)
        bcl3xl(1,jc,kc) = -bcl3xl(1,jc,kc)
        bcl4xl(1,jc,kc) = -bcl4xl(1,jc,kc)
        bcl5xl(1,jc,kc) = half*sorpxl(jc,kc)  &
            + cobcxl*acouxl(jc,kc)*(strpxl(1,jc,kc)-pinfxl) - bcl5xl(1,jc,kc)
        
      END DO
    END DO
    
!         LYX
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
!               OLD VALUE OF L's
          bclyxl(jc,kc,ispec) = struxl(1,jc,kc)*bclyxl(jc,kc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
          bclyxl(jc,kc,ispec) = ratexl(jc,kc,ispec)/strdxl(1,jc,kc)  &
              - bclyxl(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(jc,kc)
        
        urhs(istal,jc,kc) = urhs(istal,jc,kc) - bcl2xl(1,jc,kc)*struxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(jc,kc)*(struxl(1,jc,kc)+acouxl(jc,kc))
        
        vrhs(istal,jc,kc) = vrhs(istal,jc,kc) - bcl2xl(1,jc,kc)*strvxl(1,jc,kc)  &
            - bcl3xl(1,jc,kc)*strdxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(jc,kc)*strvxl(1,jc,kc)
        
        wrhs(istal,jc,kc) = wrhs(istal,jc,kc) - bcl2xl(1,jc,kc)*strwxl(1,jc,kc)  &
            - bcl4xl(1,jc,kc)*strdxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(jc,kc)*strwxl(1,jc,kc)
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc) - bcl2xl(1,jc,kc)*strexl(1,jc,kc)  &
            - bcl3xl(1,jc,kc)*strdxl(1,jc,kc)*strvxl(1,jc,kc)  &
            - bcl4xl(1,jc,kc)*strdxl(1,jc,kc)*strwxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*(ova2xl(jc,kc)*strexl(1,jc,kc)  &
            + struxl(1,jc,kc)/acouxl(jc,kc) + ovgmxl(jc,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = bclyxl(jc,kc,ispec)*strdxl(1,jc,kc)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc) - fornow*strhxl(jc,kc,ispec)
          
          yrhs(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)  &
              - (bcl2xl(1,jc,kc)+bcl5xl(1,jc,kc)*ova2xl(jc,kc))*stryxl(jc,kc,ispec)  &
              - fornow
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxl == nsbci2)THEN
    
!         INFLOW BOUNDARY CONDITION No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sydtxl(jc,kc) = zero
        sorpxl(jc,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          sydtxl(jc,kc) = sydtxl(jc,kc) + dydtxl(jc,kc,ispec)*rgspec(ispec)
          sorpxl(jc,kc) = sorpxl(jc,kc)  &
              + strhxl(jc,kc,ispec)*ratexl(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sydtxl(jc,kc) = sydtxl(jc,kc)/strrxl(1,jc,kc)
        sorpxl(jc,kc) = -sorpxl(jc,kc)*gam1xl(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1X,L2X,L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxl(1,jc,kc)*acouxl(jc,kc)*bcl1xl(1,jc,kc)
        bcl1xl(1,jc,kc) = half*(struxl(1,jc,kc)-acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)-fornow)
        bcl2xl(1,jc,kc) = struxl(1,jc,kc)  &
            *(bcl2xl(1,jc,kc)-bcl5xl(1,jc,kc)*ova2xl(jc,kc))
        bcl5xl(1,jc,kc) = half*(struxl(1,jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1X UNCHANGED
        bcl5xl(1,jc,kc) = bcl1xl(1,jc,kc)  &
            - strdxl(1,jc,kc)*acouxl(jc,kc)*dudtxl(jc,kc) - bcl5xl(1,jc,kc)
        bcl2xl(1,jc,kc) = gam1xl(jc,kc)*ova2xl(jc,kc)  &
            *(bcl1xl(1,jc,kc)+bcl5xl(1,jc,kc))  &
            + strdxl(1,jc,kc)*(dtdtxl(jc,kc)/strtxl(1,jc,kc)  &
            - sorpxl(jc,kc)/strpxl(1,jc,kc) + sydtxl(jc,kc))  &
            - bcl2xl(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(jc,kc)
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxl == nsbci3)THEN
    
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
        fornow = strdxl(1,jc,kc)*acouxl(jc,kc)*bcl1xl(1,jc,kc)
        bcl1xl(1,jc,kc) = half*(struxl(1,jc,kc)-acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)-fornow)
        bcl2xl(1,jc,kc) = struxl(1,jc,kc)  &
            *(bcl2xl(1,jc,kc)-bcl5xl(1,jc,kc)*ova2xl(jc,kc))
        bcl3xl(1,jc,kc) = struxl(1,jc,kc)*bcl3xl(1,jc,kc)
        bcl4xl(1,jc,kc) = struxl(1,jc,kc)*bcl4xl(1,jc,kc)
        bcl5xl(1,jc,kc) = half*(struxl(1,jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1X UNCHANGED
        fornow = bcl1xl(1,jc,kc) - strdxl(1,jc,kc)*acouxl(jc,kc)*dudtxl(jc,kc)
        bcl2xl(1,jc,kc) = -dddtxl(jc,kc)  &
            - ova2xl(jc,kc)*(bcl1xl(1,jc,kc)+fornow) - bcl2xl(1,jc,kc)
        bcl3xl(1,jc,kc) = -dvdtxl(jc,kc) - bcl3xl(1,jc,kc)
        bcl4xl(1,jc,kc) = -dwdtxl(jc,kc) - bcl4xl(1,jc,kc)
        bcl5xl(1,jc,kc) = fornow - bcl5xl(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc) - bcl2xl(1,jc,kc)*strexl(1,jc,kc)  &
            - bcl3xl(1,jc,kc)*strdxl(1,jc,kc)*strvxl(1,jc,kc)  &
            - bcl4xl(1,jc,kc)*strdxl(1,jc,kc)*strwxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*(ova2xl(jc,kc)*strexl(1,jc,kc)  &
            + struxl(1,jc,kc)/acouxl(jc,kc) + ovgmxl(jc,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          bclyxl(jc,kc,ispec) = ratexl(jc,kc,ispec)/strdxl(1,jc,kc)  &
              - dydtxl(jc,kc,ispec) - struxl(1,jc,kc)*bclyxl(jc,kc,ispec)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
              - bclyxl(jc,kc,ispec)*strdxl(1,jc,kc)*strhxl(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       WALL BOUNDARY CONDITIONS
!       ------------------------
  
  IF(nsbcxl == nsbcw1)THEN
    
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
        fornow = strdxl(1,jc,kc)*acouxl(jc,kc)*bcl1xl(1,jc,kc)
        bcl1xl(1,jc,kc) = half*(struxl(1,jc,kc)-acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)-fornow)
        bcl3xl(1,jc,kc) = struxl(1,jc,kc)*bcl3xl(1,jc,kc)
        bcl4xl(1,jc,kc) = struxl(1,jc,kc)*bcl4xl(1,jc,kc)
        bcl5xl(1,jc,kc) = half*(struxl(1,jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1X,L2X UNCHANGED
        bcl3xl(1,jc,kc) = -dvdtxl(jc,kc) - bcl3xl(1,jc,kc)
        bcl4xl(1,jc,kc) = -dwdtxl(jc,kc) - bcl4xl(1,jc,kc)
        bcl5xl(1,jc,kc) = bcl1xl(1,jc,kc)  &
            - strdxl(1,jc,kc)*acouxl(jc,kc)*dudtxl(jc,kc) - bcl5xl(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl5xl(1,jc,kc)*ova2xl(jc,kc)
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
            - bcl3xl(1,jc,kc)*strdxl(1,jc,kc)*strvxl(1,jc,kc)  &
            - bcl4xl(1,jc,kc)*strdxl(1,jc,kc)*strwxl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*(ova2xl(jc,kc)*strexl(1,jc,kc)  &
            + struxl(1,jc,kc)/acouxl(jc,kc) + ovgmxl(jc,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)  &
              - bcl5xl(1,jc,kc)*ova2xl(jc,kc)*stryxl(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxl == nsbcw2)THEN
    
!         WALL BOUNDARY CONDITION No 2
!         NO-SLIP WALL - ISOTHERMAL
    
!         VELOCITY AND TEMPERATURE IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxl(jc,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          sorpxl(jc,kc) = sorpxl(jc,kc)  &
              + strhxl(jc,kc,ispec)*ratexl(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxl(jc,kc) = -sorpxl(jc,kc)*gam1xl(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1X-L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxl(1,jc,kc)*acouxl(jc,kc)*bcl1xl(1,jc,kc)
        bcl1xl(1,jc,kc) = half*(struxl(1,jc,kc)-acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)-fornow)
        bcl2xl(1,jc,kc) = struxl(1,jc,kc)  &
            *(bcl2xl(1,jc,kc)-bcl5xl(1,jc,kc)*ova2xl(jc,kc))
        bcl3xl(1,jc,kc) = struxl(1,jc,kc)*bcl3xl(1,jc,kc)
        bcl4xl(1,jc,kc) = struxl(1,jc,kc)*bcl4xl(1,jc,kc)
        bcl5xl(1,jc,kc) = half*(struxl(1,jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1X UNCHANGED
        bcl3xl(1,jc,kc) = -dvdtxl(jc,kc) - bcl3xl(1,jc,kc)
        bcl4xl(1,jc,kc) = -dwdtxl(jc,kc) - bcl4xl(1,jc,kc)
        bcl5xl(1,jc,kc) = bcl1xl(1,jc,kc)  &
            - strdxl(1,jc,kc)*acouxl(jc,kc)*dudtxl(jc,kc) - bcl5xl(1,jc,kc)
        bcl2xl(1,jc,kc) = gam1xl(jc,kc)*ova2xl(jc,kc)  &
            *(bcl1xl(1,jc,kc)+bcl5xl(1,jc,kc))  &
            + strdxl(1,jc,kc)*(dtdtxl(jc,kc)/strtxl(1,jc,kc)  &
            - sorpxl(jc,kc)/strpxl(1,jc,kc)) - bcl2xl(1,jc,kc)
        
      END DO
    END DO
    
!         LYX
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
!               OLD VALUE OF LYX
          bclyxl(jc,kc,ispec) = struxl(1,jc,kc)*bclyxl(jc,kc,ispec)
          
!               UPDATE L2X
          bcl2xl(1,jc,kc) = bcl2xl(1,jc,kc) + (ratexl(jc,kc,ispec)  &
              - strdxl(1,jc,kc)*bclyxl(jc,kc,ispec)) *rgspec(ispec)/strrxl(1,jc,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(1,jc,kc)  &
            - bcl5xl(1,jc,kc)*ova2xl(jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)  &
              - (bcl2xl(1,jc,kc)+bcl5xl(1,jc,kc)*ova2xl(jc,kc))*stryxl(jc,kc,ispec)
          
        END DO
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
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        strhxr(jc,kc,ispec) = strhxr(jc,kc,ispec)  &
            - strgxr(1,jc,kc)*strtxr(1,jc,kc)*rgspec(ispec)/strrxr(1,jc,kc)
        
      END DO
      
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      gam1xr(jc,kc) = strgxr(1,jc,kc) - strrxr(1,jc,kc)
      strexr(1,jc,kc) = strexr(1,jc,kc) - gam1xr(jc,kc)*strtxr(1,jc,kc)
      
      gam1xr(jc,kc) = strrxr(1,jc,kc)/gam1xr(jc,kc)
      ovgmxr(jc,kc) = one/gam1xr(jc,kc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      fornow = strgxr(1,jc,kc)*gam1xr(jc,kc)*strtxr(1,jc,kc)
      acouxr(jc,kc) = SQRT(fornow)
      ova2xr(jc,kc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbcxr == nsbco1)THEN
    
!         OUTFLOW BC No 1
!         SUBSONIC NON-REFLECTING OUTFLOW
!         WITH OPTION TO SET PRESSURE AT INFINITY
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxr(jc,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          sorpxr(jc,kc) = sorpxr(jc,kc)  &
              + strhxr(jc,kc,ispec)*ratexr(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxr(jc,kc) = -sorpxr(jc,kc)*gam1xr(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L1X AS REQUIRED
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L1X
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)-strdxr(1,jc,kc)*acouxr(jc,kc)*bcl1xr(1,jc,kc))
        
!             SUBTRACT FROM NEW VALUE OF L1X
        bcl1xr(1,jc,kc)= half*sorpxr(jc,kc)  &
            + cobcxr*acouxr(jc,kc)*(strpxr(1,jc,kc)-pinfxr) - bcl1xr(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(1,jc,kc)*ova2xr(jc,kc)
        
        urhs(istol,jc,kc) = urhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(jc,kc)*(struxr(1,jc,kc)-acouxr(jc,kc))
        
        vrhs(istol,jc,kc) = vrhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(jc,kc)*strvxr(1,jc,kc)
        
        wrhs(istol,jc,kc) = wrhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(jc,kc)*strwxr(1,jc,kc)
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*(ova2xr(jc,kc)*strexr(1,jc,kc)  &
            - struxr(1,jc,kc)/acouxr(jc,kc) + ovgmxr(jc,kc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)  &
              - bcl1xr(1,jc,kc)*ova2xr(jc,kc)*stryxr(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
  IF(nsbcxr == nsbci1)THEN
    
!         INFLOW BC No 1
!         SUBSONIC NON-REFLECTING LAMINAR INFLOW
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxr(jc,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          sorpxr(jc,kc) = sorpxr(jc,kc)  &
              + strhxr(jc,kc,ispec)*ratexr(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxr(jc,kc) = -sorpxr(jc,kc)*gam1xr(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1X-L4X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxr(1,jc,kc)*acouxr(jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl2xr(1,jc,kc) = struxr(1,jc,kc)  &
            *(bcl2xr(1,jc,kc)-bcl5xr(1,jc,kc)*ova2xr(jc,kc))
        bcl3xr(1,jc,kc) = struxr(1,jc,kc)*bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = struxr(1,jc,kc)*bcl4xr(1,jc,kc)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
!             L5X UNCHANGED
        bcl1xr(1,jc,kc) = half*sorpxr(jc,kc)  &
            + cobcxr*acouxr(jc,kc)*(strpxr(1,jc,kc)-pinfxr) - bcl1xr(1,jc,kc)
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
          bclyxr(jc,kc,ispec) = struxr(1,jc,kc)*bclyxr(jc,kc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
          bclyxr(jc,kc,ispec) = ratexr(jc,kc,ispec)/strdxr(1,jc,kc)  &
              - bclyxr(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(1,jc,kc)*ova2xr(jc,kc)  &
            - bcl2xr(1,jc,kc)
        
        urhs(istol,jc,kc) = urhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(jc,kc)*(struxr(1,jc,kc)-acouxr(jc,kc))  &
            - bcl2xr(1,jc,kc)*struxr(1,jc,kc)
        
        vrhs(istol,jc,kc) = vrhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(jc,kc)*strvxr(1,jc,kc)  &
            - bcl2xr(1,jc,kc)*strvxr(1,jc,kc) - bcl3xr(1,jc,kc)*strdxr(1,jc,kc)
        
        wrhs(istol,jc,kc) = wrhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(jc,kc)*strwxr(1,jc,kc)  &
            - bcl2xr(1,jc,kc)*strwxr(1,jc,kc) - bcl4xr(1,jc,kc)*strdxr(1,jc,kc)
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*(ova2xr(jc,kc)*strexr(1,jc,kc)  &
            - struxr(1,jc,kc)/acouxr(jc,kc) + ovgmxr(jc,kc))  &
            - bcl2xr(1,jc,kc)*strexr(1,jc,kc)  &
            - bcl3xr(1,jc,kc)*strdxr(1,jc,kc)*strvxr(1,jc,kc)  &
            - bcl4xr(1,jc,kc)*strdxr(1,jc,kc)*strwxr(1,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = bclyxr(jc,kc,ispec)*strdxr(1,jc,kc)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc) - fornow*strhxr(jc,kc,ispec)
          
          yrhs(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)  &
              - (bcl2xr(1,jc,kc)+bcl1xr(1,jc,kc)*ova2xr(jc,kc))*stryxr(jc,kc,ispec)  &
              - fornow
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxr == nsbci2)THEN
    
!         INFLOW BOUNDARY CONDITION No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sydtxr(jc,kc) = zero
        sorpxr(jc,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          sydtxr(jc,kc) = sydtxr(jc,kc) + dydtxr(jc,kc,ispec)*rgspec(ispec)
          sorpxr(jc,kc) = sorpxr(jc,kc)  &
              + strhxr(jc,kc,ispec)*ratexr(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sydtxr(jc,kc) = sydtxr(jc,kc)/strrxr(1,jc,kc)
        sorpxr(jc,kc) = -sorpxr(jc,kc)*gam1xr(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1X,L2X,L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxr(1,jc,kc)*acouxr(jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl2xr(1,jc,kc) = struxr(1,jc,kc)  &
            *(bcl2xr(1,jc,kc)-bcl5xr(1,jc,kc)*ova2xr(jc,kc))
        bcl5xr(1,jc,kc) = half*(struxr(1,jc,kc)+acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5X UNCHANGED
        bcl1xr(1,jc,kc) = bcl5xr(1,jc,kc)  &
            + strdxr(1,jc,kc)*acouxr(jc,kc)*dudtxr(jc,kc) - bcl1xr(1,jc,kc)
        bcl2xr(1,jc,kc) = gam1xr(jc,kc)*ova2xr(jc,kc)  &
            *(bcl1xr(1,jc,kc)+bcl5xr(1,jc,kc))  &
            + strdxr(1,jc,kc)*(dtdtxr(jc,kc)/strtxr(1,jc,kc)  &
            - sorpxr(jc,kc)/strpxr(1,jc,kc) + sydtxr(jc,kc))  &
            - bcl2xr(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl2xr(1,jc,kc)  &
            - bcl1xr(1,jc,kc)*ova2xr(jc,kc)
        
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
        fornow = strdxr(1,jc,kc)*acouxr(jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl2xr(1,jc,kc) = struxr(1,jc,kc)  &
            *(bcl2xr(1,jc,kc)-bcl5xr(1,jc,kc)*ova2xr(jc,kc))
        bcl3xr(1,jc,kc) = struxr(1,jc,kc)*bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = struxr(1,jc,kc)*bcl4xr(1,jc,kc)
        bcl5xr(1,jc,kc) = half*(struxr(1,jc,kc)+acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5X UNCHANGED
        fornow = bcl5xr(1,jc,kc) + strdxr(1,jc,kc)*acouxr(jc,kc)*dudtxr(jc,kc)
        bcl1xr(1,jc,kc) = fornow - bcl1xr(1,jc,kc)
        bcl2xr(1,jc,kc) = -dddtxr(jc,kc)  &
            - ova2xr(jc,kc)*(bcl5xr(1,jc,kc)+fornow) - bcl2xr(1,jc,kc)
        bcl3xr(1,jc,kc) = -dvdtxr(jc,kc) - bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = -dwdtxr(jc,kc) - bcl4xr(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*(ova2xr(jc,kc)*strexr(1,jc,kc)  &
            - struxr(1,jc,kc)/acouxr(jc,kc) + ovgmxr(jc,kc))  &
            - bcl2xr(1,jc,kc)*strexr(1,jc,kc)  &
            - bcl3xr(1,jc,kc)*strdxr(1,jc,kc)*strvxr(1,jc,kc)  &
            - bcl4xr(1,jc,kc)*strdxr(1,jc,kc)*strwxr(1,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          bclyxr(jc,kc,ispec) = ratexr(jc,kc,ispec)/strdxr(1,jc,kc)  &
              - dydtxr(jc,kc,ispec) - struxr(1,jc,kc)*bclyxr(jc,kc,ispec)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
              - bclyxr(jc,kc,ispec)*strdxr(1,jc,kc)*strhxr(jc,kc,ispec)
          
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
        fornow = strdxr(1,jc,kc)*acouxr(jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl3xr(1,jc,kc) = struxr(1,jc,kc)*bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = struxr(1,jc,kc)*bcl4xr(1,jc,kc)
        bcl5xr(1,jc,kc) = half*(struxr(1,jc,kc)+acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L2X,L5X UNCHANGED
        bcl1xr(1,jc,kc) = bcl5xr(1,jc,kc)  &
            + strdxr(1,jc,kc)*acouxr(jc,kc)*dudtxr(jc,kc) - bcl1xr(1,jc,kc)
        bcl3xr(1,jc,kc) = -dvdtxr(jc,kc) - bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = -dwdtxr(jc,kc) - bcl4xr(1,jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(1,jc,kc)*ova2xr(jc,kc)
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(1,jc,kc)*(ova2xr(jc,kc)*strexr(1,jc,kc)  &
            + struxr(1,jc,kc)/acouxr(jc,kc) + ovgmxr(jc,kc))  &
            - bcl3xr(1,jc,kc)*strdxr(1,jc,kc)*strvxr(1,jc,kc)  &
            - bcl4xr(1,jc,kc)*strdxr(1,jc,kc)*strwxr(1,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)  &
              - bcl1xr(1,jc,kc)*ova2xr(jc,kc)*stryxr(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcxr == nsbcw2)THEN
    
!         WALL BOUNDARY CONDITION No 2
!         NO-SLIP WALL - ISOTHERMAL
    
!         VELOCITY AND TEMPERATURE IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxr(jc,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          sorpxr(jc,kc) = sorpxr(jc,kc)  &
              + strhxr(jc,kc,ispec)*ratexr(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        sorpxr(jc,kc) = -sorpxr(jc,kc)*gam1xr(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1X-L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxr(1,jc,kc)*acouxr(jc,kc)*bcl1xr(1,jc,kc)
        bcl1xr(1,jc,kc) = half*(struxr(1,jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)-fornow)
        bcl2xr(1,jc,kc) = struxr(1,jc,kc)  &
            *(bcl2xr(1,jc,kc)-bcl5xr(1,jc,kc)*ova2xr(jc,kc))
        bcl3xr(1,jc,kc) = struxr(1,jc,kc)*bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = struxr(1,jc,kc)*bcl4xr(1,jc,kc)
        bcl5xr(1,jc,kc) = half*(struxr(1,jc,kc)+acouxr(jc,kc))  &
            *(bcl5xr(1,jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5X UNCHANGED
        bcl1xr(1,jc,kc) = bcl5xr(1,jc,kc)  &
            + strdxr(1,jc,kc)*acouxr(jc,kc)*dudtxr(jc,kc) - bcl1xr(1,jc,kc)
        bcl3xr(1,jc,kc) = -dvdtxr(jc,kc) - bcl3xr(1,jc,kc)
        bcl4xr(1,jc,kc) = -dwdtxr(jc,kc) - bcl4xr(1,jc,kc)
        bcl2xr(1,jc,kc) = gam1xr(jc,kc)*ova2xr(jc,kc)  &
            *(bcl1xr(1,jc,kc)+bcl5xr(1,jc,kc))  &
            + strdxr(1,jc,kc)*(dtdtxr(jc,kc)/strtxr(1,jc,kc)  &
            - sorpxr(jc,kc)/strpxr(1,jc,kc)) - bcl2xr(1,jc,kc)
        
      END DO
    END DO
    
!         LYX
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
!               OLD VALUE OF LYX
          bclyxr(jc,kc,ispec) = struxr(1,jc,kc)*bclyxr(jc,kc,ispec)
          
!               UPDATE L2X
          bcl2xr(1,jc,kc) = bcl2xr(1,jc,kc) + (ratexr(jc,kc,ispec)  &
              - strdxr(1,jc,kc)*bclyxr(jc,kc,ispec)) *rgspec(ispec)/strrxr(1,jc,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(1,jc,kc)*ova2xr(jc,kc)  &
            - bcl2xr(1,jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)  &
              - (bcl2xr(1,jc,kc)+bcl1xr(1,jc,kc)*ova2xr(jc,kc))*stryxr(jc,kc,ispec)
          
        END DO
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
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        strhyl(ic,kc,ispec) = strhyl(ic,kc,ispec)  &
            - strgyl(ic,1,kc)*strtyl(ic,1,kc)*rgspec(ispec)/strryl(ic,1,kc)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      gam1yl(ic,kc) = strgyl(ic,1,kc) - strryl(ic,1,kc)
      streyl(ic,1,kc) = streyl(ic,1,kc) - gam1yl(ic,kc)*strtyl(ic,1,kc)
      
      gam1yl(ic,kc) = strryl(ic,1,kc)/gam1yl(ic,kc)
      ovgmyl(ic,kc) = one/gam1yl(ic,kc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      fornow = strgyl(ic,1,kc)*gam1yl(ic,kc)*strtyl(ic,1,kc)
      acouyl(ic,kc) = SQRT(fornow)
      ova2yl(ic,kc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbcyl == nsbco1)THEN
    
!         OUTFLOW BC No 1
!         SUBSONIC NON-REFLECTING OUTFLOW
!         WITH OPTION TO SET PRESSURE AT INFINITY
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyl(ic,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1, nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          sorpyl(ic,kc) = sorpyl(ic,kc)  &
              + strhyl(ic,kc,ispec)*rateyl(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyl(ic,kc) = -sorpyl(ic,kc)*gam1yl(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L5Y AS REQUIRED
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L5Y
        bcl5yl(ic,1,kc) = half*(strvyl(ic,1,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)+strdyl(ic,1,kc)*acouyl(ic,kc)*bcl1yl(ic,1,kc))
        
!             SUBTRACT FROM NEW VALUE OF L5Y
        bcl5yl(ic,1,kc)= half*sorpyl(ic,kc)  &
            + cobcyl*acouyl(ic,kc)*(strpyl(ic,1,kc)-pinfyl) - bcl5yl(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl5yl(ic,1,kc)*ova2yl(ic,kc)
        
        urhs(ic,jstal,kc) = urhs(ic,jstal,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,kc)*struyl(ic,1,kc)
        
        vrhs(ic,jstal,kc) = vrhs(ic,jstal,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,kc)*(strvyl(ic,1,kc)+acouyl(ic,kc))
        
        wrhs(ic,jstal,kc) = wrhs(ic,jstal,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,kc)*strwyl(ic,1,kc)
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
            - bcl5yl(ic,1,kc)*(ova2yl(ic,kc)*streyl(ic,1,kc)  &
            + strvyl(ic,1,kc)/acouyl(ic,kc) + ovgmyl(ic,kc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)  &
              - bcl5yl(ic,1,kc)*ova2yl(ic,kc)*stryyl(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
  IF(nsbcyl == nsbci1)THEN
    
!         INFLOW BC No 1
!         SUBSONIC NON-REFLECTING LAMINAR INFLOW
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyl(ic,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          sorpyl(ic,kc) = sorpyl(ic,kc)  &
              + strhyl(ic,kc,ispec)*rateyl(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyl(ic,kc) = -sorpyl(ic,kc)*gam1yl(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L2Y-L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyl(ic,1,kc)*acouyl(ic,kc)*bcl1yl(ic,1,kc)
        bcl2yl(ic,1,kc) = strvyl(ic,1,kc)  &
            *(bcl2yl(ic,1,kc)-bcl5yl(ic,1,kc)*ova2yl(ic,kc))
        bcl3yl(ic,1,kc) = strvyl(ic,1,kc)*bcl3yl(ic,1,kc)
        bcl4yl(ic,1,kc) = strvyl(ic,1,kc)*bcl4yl(ic,1,kc)
        bcl5yl(ic,1,kc) = half*(strvyl(ic,1,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
!             L1Y UNCHANGED
        bcl2yl(ic,1,kc) = -bcl2yl(ic,1,kc)
        bcl3yl(ic,1,kc) = -bcl3yl(ic,1,kc)
        bcl4yl(ic,1,kc) = -bcl4yl(ic,1,kc)
        bcl5yl(ic,1,kc) = half*sorpyl(ic,kc)  &
            + cobcyl*acouyl(ic,kc)*(strpyl(ic,1,kc)-pinfyl) - bcl5yl(ic,1,kc)
        
      END DO
    END DO
    
!         LYY
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
!               OLD VALUE OF L's
          bclyyl(ic,kc,ispec) = strvyl(ic,1,kc)*bclyyl(ic,kc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
          bclyyl(ic,kc,ispec) = rateyl(ic,kc,ispec)/strdyl(ic,1,kc)  &
              - bclyyl(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,kc)
        
        urhs(ic,jstal,kc) = urhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*struyl(ic,1,kc)  &
            - bcl3yl(ic,1,kc)*strdyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,kc)*struyl(ic,1,kc)
        
        vrhs(ic,jstal,kc) = vrhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*strvyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,kc)*(strvyl(ic,1,kc)+acouyl(ic,kc))
        
        wrhs(ic,jstal,kc) = wrhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*strwyl(ic,1,kc)  &
            - bcl4yl(ic,1,kc)*strdyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,kc)*strwyl(ic,1,kc)
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*streyl(ic,1,kc)  &
            - bcl3yl(ic,1,kc)*strdyl(ic,1,kc)*struyl(ic,1,kc)  &
            - bcl4yl(ic,1,kc)*strdyl(ic,1,kc)*strwyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*(ova2yl(ic,kc)*streyl(ic,1,kc)  &
            + strvyl(ic,1,kc)/acouyl(ic,kc) + ovgmyl(ic,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = bclyyl(ic,kc,ispec)*strdyl(ic,1,kc)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - fornow*strhyl(ic,kc,ispec)
          
          yrhs(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)  &
              - (bcl2yl(ic,1,kc)+bcl5yl(ic,1,kc)*ova2yl(ic,kc))*stryyl(ic,kc,ispec)  &
              - fornow
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyl == nsbci2)THEN
    
!         INFLOW BOUNDARY CONDITION No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sydtyl(ic,kc) = zero
        sorpyl(ic,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          sydtyl(ic,kc) = sydtyl(ic,kc) + dydtyl(ic,kc,ispec)*rgspec(ispec)
          sorpyl(ic,kc) = sorpyl(ic,kc)  &
              + strhyl(ic,kc,ispec)*rateyl(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sydtyl(ic,kc) = sydtyl(ic,kc)/strryl(ic,1,kc)
        sorpyl(ic,kc) = -sorpyl(ic,kc)*gam1yl(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Y,L2Y,L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyl(ic,1,kc)*acouyl(ic,kc)*bcl1yl(ic,1,kc)
        bcl1yl(ic,1,kc) = half*(strvyl(ic,1,kc)-acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)-fornow)
        bcl2yl(ic,1,kc) = strvyl(ic,1,kc)  &
            *(bcl2yl(ic,1,kc)-bcl5yl(ic,1,kc)*ova2yl(ic,kc))
        bcl5yl(ic,1,kc) = half*(strvyl(ic,1,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y UNCHANGED
        bcl5yl(ic,1,kc) = bcl1yl(ic,1,kc)  &
            - strdyl(ic,1,kc)*acouyl(ic,kc)*dvdtyl(ic,kc) - bcl5yl(ic,1,kc)
        bcl2yl(ic,1,kc) = gam1yl(ic,kc)*ova2yl(ic,kc)  &
            *(bcl1yl(ic,1,kc)+bcl5yl(ic,1,kc))  &
            + strdyl(ic,1,kc)*(dtdtyl(ic,kc)/strtyl(ic,1,kc)  &
            - sorpyl(ic,kc)/strpyl(ic,1,kc) + sydtyl(ic,kc))  &
            - bcl2yl(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,kc)
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyl == nsbci3)THEN
    
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
        fornow = strdyl(ic,1,kc)*acouyl(ic,kc)*bcl1yl(ic,1,kc)
        bcl1yl(ic,1,kc) = half*(strvyl(ic,1,kc)-acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)-fornow)
        bcl2yl(ic,1,kc) = strvyl(ic,1,kc)  &
            *(bcl2yl(ic,1,kc)-bcl5yl(ic,1,kc)*ova2yl(ic,kc))
        bcl3yl(ic,1,kc) = strvyl(ic,1,kc)*bcl3yl(ic,1,kc)
        bcl4yl(ic,1,kc) = strvyl(ic,1,kc)*bcl4yl(ic,1,kc)
        bcl5yl(ic,1,kc) = half*(strvyl(ic,1,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y UNCHANGED
        fornow = bcl1yl(ic,1,kc) - strdyl(ic,1,kc)*acouyl(ic,kc)*dvdtyl(ic,kc)
        bcl2yl(ic,1,kc) = -dddtyl(ic,kc)  &
            - ova2yl(ic,kc)*(bcl1yl(ic,1,kc)+fornow) - bcl2yl(ic,1,kc)
        bcl3yl(ic,1,kc) = -dudtyl(ic,kc) - bcl3yl(ic,1,kc)
        bcl4yl(ic,1,kc) = -dwdtyl(ic,kc) - bcl4yl(ic,1,kc)
        bcl5yl(ic,1,kc) = fornow - bcl5yl(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - bcl2yl(ic,1,kc)*streyl(ic,1,kc)  &
            - bcl3yl(ic,1,kc)*strdyl(ic,1,kc)*struyl(ic,1,kc)  &
            - bcl4yl(ic,1,kc)*strdyl(ic,1,kc)*strwyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*(ova2yl(ic,kc)*streyl(ic,1,kc)  &
            + strvyl(ic,1,kc)/acouyl(ic,kc) + ovgmyl(ic,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          bclyyl(ic,kc,ispec) = rateyl(ic,kc,ispec)/strdyl(ic,1,kc)  &
              - dydtyl(ic,kc,ispec) - strvyl(ic,1,kc)*bclyyl(ic,kc,ispec)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
              - bclyyl(ic,kc,ispec)*strdyl(ic,1,kc)*strhyl(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       WALL BOUNDARY CONDITIONS
!       ------------------------
  
  IF(nsbcyl == nsbcw1)THEN
    
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
        fornow = strdyl(ic,1,kc)*acouyl(ic,kc)*bcl1yl(ic,1,kc)
        bcl1yl(ic,1,kc) = half*(strvyl(ic,1,kc)-acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)-fornow)
        bcl3yl(ic,1,kc) = strvyl(ic,1,kc)*bcl3yl(ic,1,kc)
        bcl4yl(ic,1,kc) = strvyl(ic,1,kc)*bcl4yl(ic,1,kc)
        bcl5yl(ic,1,kc) = half*(strvyl(ic,1,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y,L2Y UNCHANGED
        bcl3yl(ic,1,kc) = -dudtyl(ic,kc) - bcl3yl(ic,1,kc)
        bcl4yl(ic,1,kc) = -dwdtyl(ic,kc) - bcl4yl(ic,1,kc)
        bcl5yl(ic,1,kc) = bcl1yl(ic,1,kc)  &
            - strdyl(ic,1,kc)*acouyl(ic,kc)*dvdtyl(ic,kc) - bcl5yl(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl5yl(ic,1,kc)*ova2yl(ic,kc)
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
            - bcl3yl(ic,1,kc)*strdyl(ic,1,kc)*struyl(ic,1,kc)  &
            - bcl4yl(ic,1,kc)*strdyl(ic,1,kc)*strwyl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*(ova2yl(ic,kc)*streyl(ic,1,kc)  &
            + strvyl(ic,1,kc)/acouyl(ic,kc) + ovgmyl(ic,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)  &
              - bcl5yl(ic,1,kc)*ova2yl(ic,kc)*stryyl(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
  
!       =======================================================================
  
  IF(nsbcyl == nsbcw2)THEN
    
!         WALL BOUNDARY CONDITION No 2
!         NO-SLIP WALL - ISOTHERMAL
    
!         VELOCITY AND TEMPERATURE IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyl(ic,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          sorpyl(ic,kc) = sorpyl(ic,kc)  &
              + strhyl(ic,kc,ispec)*rateyl(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyl(ic,kc) = -sorpyl(ic,kc)*gam1yl(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Y-L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyl(ic,1,kc)*acouyl(ic,kc)*bcl1yl(ic,1,kc)
        bcl1yl(ic,1,kc) = half*(strvyl(ic,1,kc)-acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)-fornow)
        bcl2yl(ic,1,kc) = strvyl(ic,1,kc)  &
            *(bcl2yl(ic,1,kc)-bcl5yl(ic,1,kc)*ova2yl(ic,kc))
        bcl3yl(ic,1,kc) = strvyl(ic,1,kc)*bcl3yl(ic,1,kc)
        bcl4yl(ic,1,kc) = strvyl(ic,1,kc)*bcl4yl(ic,1,kc)
        bcl5yl(ic,1,kc) = half*(strvyl(ic,1,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y UNCHANGED
        bcl3yl(ic,1,kc) = -dudtyl(ic,kc) - bcl3yl(ic,1,kc)
        bcl4yl(ic,1,kc) = -dwdtyl(ic,kc) - bcl4yl(ic,1,kc)
        bcl5yl(ic,1,kc) = bcl1yl(ic,1,kc)  &
            - strdyl(ic,1,kc)*acouyl(ic,kc)*dvdtyl(ic,kc) - bcl5yl(ic,1,kc)
        bcl2yl(ic,1,kc) = gam1yl(ic,kc)*ova2yl(ic,kc)  &
            *(bcl1yl(ic,1,kc)+bcl5yl(ic,1,kc))  &
            + strdyl(ic,1,kc)*(dtdtyl(ic,kc)/strtyl(ic,1,kc)  &
            - sorpyl(ic,kc)/strpyl(ic,1,kc)) - bcl2yl(ic,1,kc)
        
      END DO
    END DO
    
!         LYY
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYY
          bclyyl(ic,kc,ispec) = strvyl(ic,1,kc)*bclyyl(ic,kc,ispec)
          
!               UPDATE L2Y
          bcl2yl(ic,1,kc) = bcl2yl(ic,1,kc) + (rateyl(ic,kc,ispec)  &
              - strdyl(ic,1,kc)*bclyyl(ic,kc,ispec)) *rgspec(ispec)/strryl(ic,1,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,1,kc)  &
            - bcl5yl(ic,1,kc)*ova2yl(ic,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)  &
              - (bcl2yl(ic,1,kc)+bcl5yl(ic,1,kc)*ova2yl(ic,kc))*stryyl(ic,kc,ispec)
          
        END DO
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
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        strhyr(ic,kc,ispec) = strhyr(ic,kc,ispec)  &
            - strgyr(ic,1,kc)*strtyr(ic,1,kc)*rgspec(ispec)/strryr(ic,1,kc)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      gam1yr(ic,kc) = strgyr(ic,1,kc) - strryr(ic,1,kc)
      streyr(ic,1,kc) = streyr(ic,1,kc) - gam1yr(ic,kc)*strtyr(ic,1,kc)
      
      gam1yr(ic,kc) = strryr(ic,1,kc)/gam1yr(ic,kc)
      ovgmyr(ic,kc) = one/gam1yr(ic,kc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      fornow = strgyr(ic,1,kc)*gam1yr(ic,kc)*strtyr(ic,1,kc)
      acouyr(ic,kc) = SQRT(fornow)
      ova2yr(ic,kc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbcyr == nsbco1)THEN
    
!         OUTFLOW BC No 1
!         SUBSONIC NON-REFLECTING OUTFLOW
!         WITH OPTION TO SET PRESSURE AT INFINITY
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyr(ic,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          sorpyr(ic,kc) = sorpyr(ic,kc)  &
              + strhyr(ic,kc,ispec)*rateyr(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyr(ic,kc) = -sorpyr(ic,kc)*gam1yr(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L1Y AS REQUIRED
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L1Y
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)-strdyr(ic,1,kc)*acouyr(ic,kc)*bcl1yr(ic,1,kc))
        
!             SUBTRACT FROM NEW VALUE OF L1Y
        bcl1yr(ic,1,kc)= half*sorpyr(ic,kc)  &
            + cobcyr*acouyr(ic,kc)*(strpyr(ic,1,kc)-pinfyr) - bcl1yr(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,kc)
        
        urhs(ic,jstol,kc) = urhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,kc)*struyr(ic,1,kc)
        
        vrhs(ic,jstol,kc) = vrhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,kc)*(strvyr(ic,1,kc)-acouyr(ic,kc))
        
        wrhs(ic,jstol,kc) = wrhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,kc)*strwyr(ic,1,kc)
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*(ova2yr(ic,kc)*streyr(ic,1,kc)  &
            - strvyr(ic,1,kc)/acouyr(ic,kc) + ovgmyr(ic,kc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)  &
              - bcl1yr(ic,1,kc)*ova2yr(ic,kc)*stryyr(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
  IF(nsbcyr == nsbci1)THEN
    
!         INFLOW BC No 1
!         SUBSONIC NON-REFLECTING LAMINAR INFLOW
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyr(ic,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          sorpyr(ic,kc) = sorpyr(ic,kc)  &
              + strhyr(ic,kc,ispec)*rateyr(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyr(ic,kc) = -sorpyr(ic,kc)*gam1yr(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Y-L4Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyr(ic,1,kc)*acouyr(ic,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl2yr(ic,1,kc) = strvyr(ic,1,kc)  &
            *(bcl2yr(ic,1,kc)-bcl5yr(ic,1,kc)*ova2yr(ic,kc))
        bcl3yr(ic,1,kc) = strvyr(ic,1,kc)*bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = strvyr(ic,1,kc)*bcl4yr(ic,1,kc)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
!             L5Y UNCHANGED
        bcl1yr(ic,1,kc) = half*sorpyr(ic,kc)  &
            + cobcyr*acouyr(ic,kc)*(strpyr(ic,1,kc)-pinfyr) - bcl1yr(ic,1,kc)
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
          bclyyr(ic,kc,ispec) = strvyr(ic,1,kc)*bclyyr(ic,kc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
          bclyyr(ic,kc,ispec) = rateyr(ic,kc,ispec)/strdyr(ic,1,kc)  &
              - bclyyr(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,kc)  &
            - bcl2yr(ic,1,kc)
        
        urhs(ic,jstol,kc) = urhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,kc)*struyr(ic,1,kc)  &
            - bcl2yr(ic,1,kc)*struyr(ic,1,kc) - bcl3yr(ic,1,kc)*strdyr(ic,1,kc)
        
        vrhs(ic,jstol,kc) = vrhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,kc)*(strvyr(ic,1,kc)-acouyr(ic,kc))  &
            - bcl2yr(ic,1,kc)*strvyr(ic,1,kc)
        
        wrhs(ic,jstol,kc) = wrhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*ova2yr(ic,kc)*strwyr(ic,1,kc)  &
            - bcl2yr(ic,1,kc)*strwyr(ic,1,kc) - bcl4yr(ic,1,kc)*strdyr(ic,1,kc)
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*(ova2yr(ic,kc)*streyr(ic,1,kc)  &
            + strvyr(ic,1,kc)/acouyr(ic,kc) + ovgmyr(ic,kc))  &
            - bcl2yr(ic,1,kc)*streyr(ic,1,kc)  &
            - bcl3yr(ic,1,kc)*strdyr(ic,1,kc)*struyr(ic,1,kc)  &
            - bcl4yr(ic,1,kc)*strdyr(ic,1,kc)*strwyr(ic,1,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = bclyyr(ic,kc,ispec)*strdyr(ic,1,kc)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc) - fornow*strhyr(ic,kc,ispec)
          
          yrhs(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)  &
              - (bcl2yr(ic,1,kc)+bcl1yr(ic,1,kc)*ova2yr(ic,kc))*stryyr(ic,kc,ispec)  &
              - fornow
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyr == nsbci2)THEN
    
!         INFLOW BOUNDARY CONDITION No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sydtyr(ic,kc) = zero
        sorpyr(ic,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          sydtyr(ic,kc) = sydtyr(ic,kc) + dydtyr(ic,kc,ispec)*rgspec(ispec)
          sorpyr(ic,kc) = sorpyr(ic,kc)  &
              + strhyr(ic,kc,ispec)*rateyr(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sydtyr(ic,kc) = sydtyr(ic,kc)/strryr(ic,1,kc)
        sorpyr(ic,kc) = -sorpyr(ic,kc)*gam1yr(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Y,L2Y,L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyr(ic,1,kc)*acouyr(ic,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl2yr(ic,1,kc) = strvyr(ic,1,kc)  &
            *(bcl2yr(ic,1,kc)-bcl5yr(ic,1,kc)*ova2yr(ic,kc))
        bcl5yr(ic,1,kc) = half*(strvyr(ic,1,kc)+acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Y UNCHANGED
        bcl1yr(ic,1,kc) = bcl5yr(ic,1,kc)  &
            + strdyr(ic,1,kc)*acouyr(ic,kc)*dvdtyr(ic,kc) - bcl1yr(ic,1,kc)
        bcl2yr(ic,1,kc) = gam1yr(ic,kc)*ova2yr(ic,kc)  &
            *(bcl1yr(ic,1,kc)+bcl5yr(ic,1,kc))  &
            + strdyr(ic,1,kc)*(dtdtyr(ic,kc)/strtyr(ic,1,kc)  &
            - sorpyr(ic,kc)/strpyr(ic,1,kc) + sydtyr(ic,kc))  &
            - bcl2yr(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,kc)  &
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
        fornow = strdyr(ic,1,kc)*acouyr(ic,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl2yr(ic,1,kc) = strvyr(ic,1,kc)  &
            *(bcl2yr(ic,1,kc)-bcl5yr(ic,1,kc)*ova2yr(ic,kc))
        bcl3yr(ic,1,kc) = strvyr(ic,1,kc)*bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = strvyr(ic,1,kc)*bcl4yr(ic,1,kc)
        bcl5yr(ic,1,kc) = half*(strvyr(ic,1,kc)+acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Y UNCHANGED
        fornow = bcl5yr(ic,1,kc) + strdyr(ic,1,kc)*acouyr(ic,kc)*dvdtyr(ic,kc)
        bcl1yr(ic,1,kc) = fornow - bcl1yr(ic,1,kc)
        bcl2yr(ic,1,kc) = -dddtyr(ic,kc)  &
            - ova2yr(ic,kc)*(bcl1yr(ic,1,kc)+fornow) - bcl2yr(ic,1,kc)
        bcl3yr(ic,1,kc) = -dudtyr(ic,kc) - bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = -dwdtyr(ic,kc) - bcl4yr(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*(ova2yr(ic,kc)*streyr(ic,1,kc)  &
            + strvyr(ic,1,kc)/acouyr(ic,kc) + ovgmyr(ic,kc))  &
            - bcl2yr(ic,1,kc)*streyr(ic,1,kc)  &
            - bcl3yr(ic,1,kc)*strdyr(ic,1,kc)*struyr(ic,1,kc)  &
            - bcl4yr(ic,1,kc)*strdyr(ic,1,kc)*strwyr(ic,1,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          bclyyr(ic,kc,ispec) = rateyr(ic,kc,ispec)/strdyr(ic,1,kc)  &
              - dydtyr(ic,kc,ispec) - strvyr(ic,1,kc)*bclyyr(ic,kc,ispec)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
              - bclyyr(ic,kc,ispec)*strdyr(ic,1,kc)*strhyr(ic,kc,ispec)
          
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
        fornow = strdyr(ic,1,kc)*acouyr(ic,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl3yr(ic,1,kc) = strvyr(ic,1,kc)*bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = strvyr(ic,1,kc)*bcl4yr(ic,1,kc)
        bcl5yr(ic,1,kc) = half*(strvyr(ic,1,kc)+acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L2Y,L5Y UNCHANGED
        bcl1yr(ic,1,kc) = bcl5yr(ic,1,kc)  &
            + strdyr(ic,1,kc)*acouyr(ic,kc)*dvdtyr(ic,kc) - bcl1yr(ic,1,kc)
        bcl3yr(ic,1,kc) = -dudtyr(ic,kc) - bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = -dwdtyr(ic,kc) - bcl4yr(ic,1,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,kc)
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,1,kc)*(ova2yr(ic,kc)*streyr(ic,1,kc)  &
            + strvyr(ic,1,kc)/acouyr(ic,kc) + ovgmyr(ic,kc))  &
            - bcl3yr(ic,1,kc)*strdyr(ic,1,kc)*struyr(ic,1,kc)  &
            - bcl4yr(ic,1,kc)*strdyr(ic,1,kc)*strwyr(ic,1,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)  &
              - bcl1yr(ic,1,kc)*ova2yr(ic,kc)*stryyr(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbcyr == nsbcw2)THEN
    
!         WALL BOUNDARY CONDITION No 2
!         NO-SLIP WALL - ISOTHERMAL
    
!         VELOCITY AND TEMPERATURE IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyr(ic,kc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          sorpyr(ic,kc) = sorpyr(ic,kc)  &
              + strhyr(ic,kc,ispec)*rateyr(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        sorpyr(ic,kc) = -sorpyr(ic,kc)*gam1yr(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Y-L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyr(ic,1,kc)*acouyr(ic,kc)*bcl1yr(ic,1,kc)
        bcl1yr(ic,1,kc) = half*(strvyr(ic,1,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)-fornow)
        bcl2yr(ic,1,kc) = strvyr(ic,1,kc)  &
            *(bcl2yr(ic,1,kc)-bcl5yr(ic,1,kc)*ova2yr(ic,kc))
        bcl3yr(ic,1,kc) = strvyr(ic,1,kc)*bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = strvyr(ic,1,kc)*bcl4yr(ic,1,kc)
        bcl5yr(ic,1,kc) = half*(strvyr(ic,1,kc)+acouyr(ic,kc))  &
            *(bcl5yr(ic,1,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Y UNCHANGED
        bcl1yr(ic,1,kc) = bcl5yr(ic,1,kc)  &
            + strdyr(ic,1,kc)*acouyr(ic,kc)*dvdtyr(ic,kc) - bcl1yr(ic,1,kc)
        bcl3yr(ic,1,kc) = -dudtyr(ic,kc) - bcl3yr(ic,1,kc)
        bcl4yr(ic,1,kc) = -dwdtyr(ic,kc) - bcl4yr(ic,1,kc)
        bcl2yr(ic,1,kc) = gam1yr(ic,kc)*ova2yr(ic,kc)  &
            *(bcl1yr(ic,1,kc)+bcl5yr(ic,1,kc))  &
            + strdyr(ic,1,kc)*(dtdtyr(ic,kc)/strtyr(ic,1,kc)  &
            - sorpyr(ic,kc)/strpyr(ic,1,kc)) - bcl2yr(ic,1,kc)
        
      END DO
    END DO
    
!         LYY
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYY
          bclyyr(ic,kc,ispec) = strvyr(ic,1,kc)*bclyyr(ic,kc,ispec)
          
!               UPDATE L2Y
          bcl2yr(ic,1,kc) = bcl2yr(ic,1,kc) + (rateyr(ic,kc,ispec)  &
              - strdyr(ic,1,kc)*bclyyr(ic,kc,ispec)) *rgspec(ispec)/strryr(ic,1,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,1,kc)*ova2yr(ic,kc)  &
            - bcl2yr(ic,1,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)  &
              - (bcl2yr(ic,1,kc)+bcl1yr(ic,1,kc)*ova2yr(ic,kc))*stryyr(ic,kc,ispec)
          
        END DO
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
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        strhzl(ic,jc,ispec) = strhzl(ic,jc,ispec)  &
            - strgzl(ic,jc,1)*strtzl(ic,jc,1)*rgspec(ispec)/strrzl(ic,jc,1)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      gam1zl(ic,jc) = strgzl(ic,jc,1) - strrzl(ic,jc,1)
      strezl(ic,jc,1) = strezl(ic,jc,1) - gam1zl(ic,jc)*strtzl(ic,jc,1)
      
      gam1zl(ic,jc) = strrzl(ic,jc,1)/gam1zl(ic,jc)
      ovgmzl(ic,jc) = one/gam1zl(ic,jc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = strgzl(ic,jc,1)*gam1zl(ic,jc)*strtzl(ic,jc,1)
      acouzl(ic,jc) = SQRT(fornow)
      ova2zl(ic,jc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbczl == nsbco1)THEN
    
!         OUTFLOW BC No 1
!         SUBSONIC NON-REFLECTING OUTFLOW
!         WITH OPTION TO SET PRESSURE AT INFINITY
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzl(ic,jc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          sorpzl(ic,jc) = sorpzl(ic,jc)  &
              + strhzl(ic,jc,ispec)*ratezl(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzl(ic,jc) = -sorpzl(ic,jc)*gam1zl(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L5Z AS REQUIRED
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L5Z
        bcl5zl(ic,jc,1) = half*(strwzl(ic,jc,1)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)+strdzl(ic,jc,1)*acouzl(ic,jc)*bcl1zl(ic,jc,1))
        
!             SUBTRACT FROM NEW VALUE OF L5Z
        bcl5zl(ic,jc,1)= half*sorpzl(ic,jc)  &
            + cobczl*acouzl(ic,jc)*(strpzl(ic,jc,1)-pinfzl) - bcl5zl(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl5zl(ic,jc,1)*ova2zl(ic,jc)
        
        urhs(ic,jc,kstal) = urhs(ic,jc,kstal)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc)*struzl(ic,jc,1)
        
        vrhs(ic,jc,kstal) = vrhs(ic,jc,kstal)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc)*strvzl(ic,jc,1)
        
        wrhs(ic,jc,kstal) = wrhs(ic,jc,kstal)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc)*(strwzl(ic,jc,1)+acouzl(ic,jc))
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
            - bcl5zl(ic,jc,1)*(ova2zl(ic,jc)*strezl(ic,jc,1)  &
            + strwzl(ic,jc,1)/acouzl(ic,jc) + ovgmzl(ic,jc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)  &
              - bcl5zl(ic,jc,1)*ova2zl(ic,jc)*stryzl(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
  IF(nsbczl == nsbci1)THEN
    
!         INFLOW BC No 1
!         SUBSONIC NON-REFLECTING LAMINAR INFLOW
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzl(ic,jc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          sorpzl(ic,jc) = sorpzl(ic,jc)  &
              + strhzl(ic,jc,ispec)*ratezl(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzl(ic,jc) = -sorpzl(ic,jc)*gam1zl(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L2Z-L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzl(ic,jc,1)*acouzl(ic,jc)*bcl1zl(ic,jc,1)
        bcl2zl(ic,jc,1) = strwzl(ic,jc,1)  &
            *(bcl2zl(ic,jc,1)-bcl5zl(ic,jc,1)*ova2zl(ic,jc))
        bcl3zl(ic,jc,1) = strwzl(ic,jc,1)*bcl3zl(ic,jc,1)
        bcl4zl(ic,jc,1) = strwzl(ic,jc,1)*bcl4zl(ic,jc,1)
        bcl5zl(ic,jc,1) = half*(strwzl(ic,jc,1)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
!             L1Z UNCHANGED
        bcl2zl(ic,jc,1) = -bcl2zl(ic,jc,1)
        bcl3zl(ic,jc,1) = -bcl3zl(ic,jc,1)
        bcl4zl(ic,jc,1) = -bcl4zl(ic,jc,1)
        bcl5zl(ic,jc,1) = half*sorpzl(ic,jc)  &
            + cobczl*acouzl(ic,jc)*(strpzl(ic,jc,1)-pinfzl) - bcl5zl(ic,jc,1)
        
      END DO
    END DO
    
!         LYZ
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
!               OLD VALUE OF L's
          bclyzl(ic,jc,ispec) = strwzl(ic,jc,1)*bclyzl(ic,jc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
          bclyzl(ic,jc,ispec) = ratezl(ic,jc,ispec)/strdzl(ic,jc,1)  &
              - bclyzl(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc)
        
        urhs(ic,jc,kstal) = urhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*struzl(ic,jc,1)  &
            - bcl3zl(ic,jc,1)*strdzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc)*struzl(ic,jc,1)
        
        vrhs(ic,jc,kstal) = vrhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*strvzl(ic,jc,1)  &
            - bcl4zl(ic,jc,1)*strdzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc)*strvzl(ic,jc,1)
        
        wrhs(ic,jc,kstal) = wrhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*strwzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc)*(strwzl(ic,jc,1)+acouzl(ic,jc))
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*strezl(ic,jc,1)  &
            - bcl3zl(ic,jc,1)*strdzl(ic,jc,1)*struzl(ic,jc,1)  &
            - bcl4zl(ic,jc,1)*strdzl(ic,jc,1)*strvzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*(ova2zl(ic,jc)*strezl(ic,jc,1)  &
            + strwzl(ic,jc,1)/acouzl(ic,jc) + ovgmzl(ic,jc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = bclyzl(ic,jc,ispec)*strdzl(ic,jc,1)
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - fornow*strhzl(ic,jc,ispec)
          
          yrhs(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)  &
              - (bcl2zl(ic,jc,1)+bcl5zl(ic,jc,1)*ova2zl(ic,jc))*stryzl(ic,jc,ispec)  &
              - fornow
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczl == nsbci2)THEN
    
!         INFLOW BOUNDARY CONDITION No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sydtzl(ic,jc) = zero
        sorpzl(ic,jc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          sydtzl(ic,jc) = sydtzl(ic,jc) + dydtzl(ic,jc,ispec)*rgspec(ispec)
          sorpzl(ic,jc) = sorpzl(ic,jc)  &
              + strhzl(ic,jc,ispec)*ratezl(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sydtzl(ic,jc) = sydtzl(ic,jc)/strrzl(ic,jc,1)
        sorpzl(ic,jc) = -sorpzl(ic,jc)*gam1zl(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Z,L2Z,L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzl(ic,jc,1)*acouzl(ic,jc)*bcl1zl(ic,jc,1)
        bcl1zl(ic,jc,1) = half*(strwzl(ic,jc,1)-acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)-fornow)
        bcl2zl(ic,jc,1) = strwzl(ic,jc,1)  &
            *(bcl2zl(ic,jc,1)-bcl5zl(ic,jc,1)*ova2zl(ic,jc))
        bcl5zl(ic,jc,1) = half*(strwzl(ic,jc,1)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Z UNCHANGED
        bcl5zl(ic,jc,1) = bcl1zl(ic,jc,1)  &
            - strdzl(ic,jc,1)*acouzl(ic,jc)*dwdtzl(ic,jc) - bcl5zl(ic,jc,1)
        bcl2zl(ic,jc,1) = gam1zl(ic,jc)*ova2zl(ic,jc)  &
            *(bcl1zl(ic,jc,1)+bcl5zl(ic,jc,1))  &
            + strdzl(ic,jc,1)*(dtdtzl(ic,jc)/strtzl(ic,jc,1)  &
            - sorpzl(ic,jc)/strpzl(ic,jc,1) + sydtzl(ic,jc))  &
            - bcl2zl(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc)
        
      END DO
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczl == nsbci3)THEN
    
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
        fornow = strdzl(ic,jc,1)*acouzl(ic,jc)*bcl1zl(ic,jc,1)
        bcl1zl(ic,jc,1) = half*(strwzl(ic,jc,1)-acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)-fornow)
        bcl2zl(ic,jc,1) = strwzl(ic,jc,1)  &
            *(bcl2zl(ic,jc,1)-bcl5zl(ic,jc,1)*ova2zl(ic,jc))
        bcl3zl(ic,jc,1) = strwzl(ic,jc,1)*bcl3zl(ic,jc,1)
        bcl4zl(ic,jc,1) = strwzl(ic,jc,1)*bcl4zl(ic,jc,1)
        bcl5zl(ic,jc,1) = half*(strwzl(ic,jc,1)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Z UNCHANGED
        fornow = bcl1zl(ic,jc,1) - strdzl(ic,jc,1)*acouzl(ic,jc)*dwdtzl(ic,jc)
        bcl2zl(ic,jc,1) = -dddtzl(ic,jc)  &
            - ova2zl(ic,jc)*(bcl1zl(ic,jc,1)+fornow) - bcl2zl(ic,jc,1)
        bcl3zl(ic,jc,1) = -dudtzl(ic,jc) - bcl3zl(ic,jc,1)
        bcl4zl(ic,jc,1) = -dvdtzl(ic,jc) - bcl4zl(ic,jc,1)
        bcl5zl(ic,jc,1) = fornow - bcl5zl(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - bcl2zl(ic,jc,1)*strezl(ic,jc,1)  &
            - bcl3zl(ic,jc,1)*strdzl(ic,jc,1)*struzl(ic,jc,1)  &
            - bcl4zl(ic,jc,1)*strdzl(ic,jc,1)*strvzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*(ova2zl(ic,jc)*strezl(ic,jc,1)  &
            + strwzl(ic,jc,1)/acouzl(ic,jc) + ovgmzl(ic,jc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          bclyzl(ic,jc,ispec) = ratezl(ic,jc,ispec)/strdzl(ic,jc,1)  &
              - dydtzl(ic,jc,ispec) - strwzl(ic,jc,1)*bclyzl(ic,jc,ispec)
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
              - bclyzl(ic,jc,ispec)*strdzl(ic,jc,1)*strhzl(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       WALL BOUNDARY CONDITIONS
!       ------------------------
  
  IF(nsbczl == nsbcw1)THEN
    
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
        fornow = strdzl(ic,jc,1)*acouzl(ic,jc)*bcl1zl(ic,jc,1)
        bcl1zl(ic,jc,1) = half*(strwzl(ic,jc,1)-acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)-fornow)
        bcl3zl(ic,jc,1) = strwzl(ic,jc,1)*bcl3zl(ic,jc,1)
        bcl4zl(ic,jc,1) = strwzl(ic,jc,1)*bcl4zl(ic,jc,1)
        bcl5zl(ic,jc,1) = half*(strwzl(ic,jc,1)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Z,L2Z UNCHANGED
        bcl3zl(ic,jc,1) = -dudtzl(ic,jc) - bcl3zl(ic,jc,1)
        bcl4zl(ic,jc,1) = -dvdtzl(ic,jc) - bcl4zl(ic,jc,1)
        bcl5zl(ic,jc,1) = bcl1zl(ic,jc,1)  &
            - strdzl(ic,jc,1)*acouzl(ic,jc)*dwdtzl(ic,jc) - bcl5zl(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl5zl(ic,jc,1)*ova2zl(ic,jc)
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
            - bcl3zl(ic,jc,1)*strdzl(ic,jc,1)*struzl(ic,jc,1)  &
            - bcl4zl(ic,jc,1)*strdzl(ic,jc,1)*strvzl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*(ova2zl(ic,jc)*strezl(ic,jc,1)  &
            + strwzl(ic,jc,1)/acouzl(ic,jc) + ovgmzl(ic,jc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)  &
              - bcl5zl(ic,jc,1)*ova2zl(ic,jc)*stryzl(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczl == nsbcw2)THEN
    
!         WALL BOUNDARY CONDITION No 2
!         NO-SLIP WALL - ISOTHERMAL
    
!         VELOCITY AND TEMPERATURE IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzl(ic,jc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          sorpzl(ic,jc) = sorpzl(ic,jc)  &
              + strhzl(ic,jc,ispec)*ratezl(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzl(ic,jc) = -sorpzl(ic,jc)*gam1zl(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Z-L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzl(ic,jc,1)*acouzl(ic,jc)*bcl1zl(ic,jc,1)
        bcl1zl(ic,jc,1) = half*(strwzl(ic,jc,1)-acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)-fornow)
        bcl2zl(ic,jc,1) = strwzl(ic,jc,1)  &
            *(bcl2zl(ic,jc,1)-bcl5zl(ic,jc,1)*ova2zl(ic,jc))
        bcl3zl(ic,jc,1) = strwzl(ic,jc,1)*bcl3zl(ic,jc,1)
        bcl4zl(ic,jc,1) = strwzl(ic,jc,1)*bcl4zl(ic,jc,1)
        bcl5zl(ic,jc,1) = half*(strwzl(ic,jc,1)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y UNCHANGED
        bcl3zl(ic,jc,1) = -dudtzl(ic,jc) - bcl3zl(ic,jc,1)
        bcl4zl(ic,jc,1) = -dvdtzl(ic,jc) - bcl4zl(ic,jc,1)
        bcl5zl(ic,jc,1) = bcl1zl(ic,jc,1)  &
            - strdzl(ic,jc,1)*acouzl(ic,jc)*dwdtzl(ic,jc) - bcl5zl(ic,jc,1)
        bcl2zl(ic,jc,1) = gam1zl(ic,jc)*ova2zl(ic,jc)  &
            *(bcl1zl(ic,jc,1)+bcl5zl(ic,jc,1))  &
            + strdzl(ic,jc,1)*(dtdtzl(ic,jc)/strtzl(ic,jc,1)  &
            - sorpzl(ic,jc)/strpzl(ic,jc,1)) - bcl2zl(ic,jc,1)
        
      END DO
    END DO
    
!         LYZ
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYZ
          bclyzl(ic,jc,ispec) = strwzl(ic,jc,1)*bclyzl(ic,jc,ispec)
          
!               UPDATE L2Z
          bcl2zl(ic,jc,1) = bcl2zl(ic,jc,1) + (ratezl(ic,jc,ispec)  &
              - strdzl(ic,jc,1)*bclyzl(ic,jc,ispec)) *rgspec(ispec)/strrzl(ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc,1)  &
            - bcl5zl(ic,jc,1)*ova2zl(ic,jc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)  &
              - (bcl2zl(ic,jc,1)+bcl5zl(ic,jc,1)*ova2zl(ic,jc))*stryzl(ic,jc,ispec)
          
        END DO
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
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        strhzr(ic,jc,ispec) = strhzr(ic,jc,ispec)  &
            - strgzr(ic,jc,1)*strtzr(ic,jc,1)*rgspec(ispec)/strrzr(ic,jc,1)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      gam1zr(ic,jc) = strgzr(ic,jc,1) - strrzr(ic,jc,1)
      strezr(ic,jc,1) = strezr(ic,jc,1) - gam1zr(ic,jc)*strtzr(ic,jc,1)
      
      gam1zr(ic,jc) = strrzr(ic,jc,1)/gam1zr(ic,jc)
      ovgmzr(ic,jc) = one/gam1zr(ic,jc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = strgzr(ic,jc,1)*gam1zr(ic,jc)*strtzr(ic,jc,1)
      acouzr(ic,jc) = SQRT(fornow)
      ova2zr(ic,jc) = one/fornow
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbczr == nsbco1)THEN
    
!         OUTFLOW BC No 1
!         SUBSONIC NON-REFLECTING OUTFLOW
!         WITH OPTION TO SET PRESSURE AT INFINITY
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzr(ic,jc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          sorpzr(ic,jc) = sorpzr(ic,jc)  &
              + strhzr(ic,jc,ispec)*ratezr(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzr(ic,jc) = -sorpzr(ic,jc)*gam1zr(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L1Z AS REQUIRED
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L1Z
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc,1)-strdzr(ic,jc,1)*acouzr(ic,jc)*bcl1zr(ic,jc,1))
        
!             SUBTRACT FROM NEW VALUE OF L1Z
        bcl1zr(ic,jc,1)= half*sorpzr(ic,jc)  &
            + cobczr*acouzr(ic,jc)*(strpzr(ic,jc,1)-pinfzr) - bcl1zr(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc)
        
        urhs(ic,jc,kstol) = urhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc)*struzr(ic,jc,1)
        
        vrhs(ic,jc,kstol) = vrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc)*strvzr(ic,jc,1)
        
        wrhs(ic,jc,kstol) = wrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc)*(strwzr(ic,jc,1)-acouzr(ic,jc))
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*(ova2zr(ic,jc)*strezr(ic,jc,1)  &
            - strwzr(ic,jc,1)/acouzr(ic,jc) + ovgmzr(ic,jc))
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)  &
              - bcl1zr(ic,jc,1)*ova2zr(ic,jc)*stryzr(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
!       INFLOW BOUNDARY CONDITIONS
!       --------------------------
  
  IF(nsbczr == nsbci1)THEN
    
!         INFLOW BC No 1
!         SUBSONIC NON-REFLECTING LAMINAR INFLOW
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzr(ic,jc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          sorpzr(ic,jc) = sorpzr(ic,jc)  &
              + strhzr(ic,jc,ispec)*ratezr(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzr(ic,jc) = -sorpzr(ic,jc)*gam1zr(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Y-L4Y
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl2zr(ic,jc,1) = strwzr(ic,jc,1)  &
            *(bcl2zr(ic,jc,1)-bcl5zr(ic,jc,1)*ova2zr(ic,jc))
        bcl3zr(ic,jc,1) = strwzr(ic,jc,1)*bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = strwzr(ic,jc,1)*bcl4zr(ic,jc,1)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
!             L5Z UNCHANGED
        bcl1zr(ic,jc,1) = half*sorpzr(ic,jc)  &
            + cobczr*acouzr(ic,jc)*(strpzr(ic,jc,1)-pinfzr) - bcl1zr(ic,jc,1)
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
          bclyzr(ic,jc,ispec) = strwzr(ic,jc,1)*bclyzr(ic,jc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
          bclyzr(ic,jc,ispec) = ratezr(ic,jc,ispec)/strdzr(ic,jc,1)  &
              - bclyzr(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc)  &
            - bcl2zr(ic,jc,1)
        
        urhs(ic,jc,kstol) = urhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc)*struzr(ic,jc,1)  &
            - bcl2zr(ic,jc,1)*struzr(ic,jc,1) - bcl3zr(ic,jc,1)*strdzr(ic,jc,1)
        
        vrhs(ic,jc,kstol) = vrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc)*strvzr(ic,jc,1)  &
            - bcl2zr(ic,jc,1)*strvzr(ic,jc,1) - bcl4zr(ic,jc,1)*strdzr(ic,jc,1)
        
        wrhs(ic,jc,kstol) = wrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*ova2zr(ic,jc)*(strwzr(ic,jc,1)-acouzr(ic,jc))  &
            - bcl2zr(ic,jc,1)*strwzr(ic,jc,1)
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*(ova2zr(ic,jc)*strezr(ic,jc,1)  &
            + strwzr(ic,jc,1)/acouzr(ic,jc) + ovgmzr(ic,jc))  &
            - bcl2zr(ic,jc,1)*strezr(ic,jc,1)  &
            - bcl3zr(ic,jc,1)*strdzr(ic,jc,1)*struzr(ic,jc,1)  &
            - bcl4zr(ic,jc,1)*strdzr(ic,jc,1)*strvzr(ic,jc,1)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = bclyzr(ic,jc,ispec)*strdzr(ic,jc,1)
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol) - fornow*strhzr(ic,jc,ispec)
          
          yrhs(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)  &
              - (bcl2zr(ic,jc,1)+bcl5zr(ic,jc,1)*ova2zr(ic,jc))*stryzr(ic,jc,ispec)  &
              - fornow
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczr == nsbci2)THEN
    
!         INFLOW BOUNDARY CONDITION No 2
!         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sydtzr(ic,jc) = zero
        sorpzr(ic,jc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          sydtzr(ic,jc) = sydtzr(ic,jc) + dydtzr(ic,jc,ispec)*rgspec(ispec)
          sorpzr(ic,jc) = sorpzr(ic,jc)  &
              + strhzr(ic,jc,ispec)*ratezr(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sydtzr(ic,jc) = sydtzr(ic,jc)/strrzr(ic,jc,1)
        sorpzr(ic,jc) = -sorpzr(ic,jc)*gam1zr(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Z,L2Z,L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl2zr(ic,jc,1) = strwzr(ic,jc,1)  &
            *(bcl2zr(ic,jc,1)-bcl5zr(ic,jc,1)*ova2zr(ic,jc))
        bcl5zr(ic,jc,1) = half*(strwzr(ic,jc,1)+acouzr(ic,jc))  &
            *(bcl5zr(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Z UNCHANGED
        bcl1zr(ic,jc,1) = bcl5zr(ic,jc,1)  &
            + strdzr(ic,jc,1)*acouzr(ic,jc)*dwdtzr(ic,jc) - bcl1zr(ic,jc,1)
        bcl2zr(ic,jc,1) = gam1zr(ic,jc)*ova2zr(ic,jc)  &
            *(bcl1zr(ic,jc,1)+bcl5zr(ic,jc,1))  &
            + strdzr(ic,jc,1)*(dtdtzr(ic,jc)/strtzr(ic,jc,1)  &
            - sorpzr(ic,jc)/strpzr(ic,jc,1) + sydtzr(ic,jc))  &
            - bcl2zr(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc)  &
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
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl2zr(ic,jc,1) = strwzr(ic,jc,1)  &
            *(bcl2zr(ic,jc,1)-bcl5zr(ic,jc,1)*ova2zr(ic,jc))
        bcl3zr(ic,jc,1) = strwzr(ic,jc,1)*bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = strwzr(ic,jc,1)*bcl4zr(ic,jc,1)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Z UNCHANGED
        fornow = bcl5zr(ic,jc,1) + strdzr(ic,jc,1)*acouzr(ic,jc)*dwdtzr(ic,jc)
        bcl1zr(ic,jc,1) = fornow - bcl1zr(ic,jc,1)
        bcl2zr(ic,jc,1) = -dddtzr(ic,jc)  &
            - ova2zr(ic,jc)*(bcl1zr(ic,jc,1)+fornow) - bcl2zr(ic,jc,1)
        bcl3zr(ic,jc,1) = -dudtzr(ic,jc) - bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = -dvdtzr(ic,jc) - bcl4zr(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*(ova2zr(ic,jc)*strezr(ic,jc,1)  &
            + strwzr(ic,jc,1)/acouzr(ic,jc) + ovgmzr(ic,jc))  &
            - bcl2zr(ic,jc,1)*strezr(ic,jc,1)  &
            - bcl3zr(ic,jc,1)*strdzr(ic,jc,1)*struzr(ic,jc,1)  &
            - bcl4zr(ic,jc,1)*strdzr(ic,jc,1)*strvzr(ic,jc,1)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          bclyzr(ic,jc,ispec) = ratezr(ic,jc,ispec)/strdzr(ic,jc,1)  &
              - dydtzr(ic,jc,ispec) - strwzr(ic,jc,1)*bclyzr(ic,jc,ispec)
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
              - bclyzr(ic,jc,ispec)*strdzr(ic,jc,1)*strhzr(ic,jc,ispec)
          
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
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl3zr(ic,jc,1) = strwzr(ic,jc,1)*bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = strwzr(ic,jc,1)*bcl4zr(ic,jc,1)
        bcl5zr(ic,jc,1) = half*(strwzr(ic,jc,1)+acouzr(ic,jc))  &
            *(bcl5zr(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L2Z,L5Z UNCHANGED
        bcl1zr(ic,jc,1) = bcl5zr(ic,jc,1)  &
            + strdzr(ic,jc,1)*acouzr(ic,jc)*dwdtzr(ic,jc) - bcl1zr(ic,jc,1)
        bcl3zr(ic,jc,1) = -dudtzr(ic,jc) - bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = -dvdtzr(ic,jc) - bcl4zr(ic,jc,1)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc)
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc,1)*(ova2zr(ic,jc)*strezr(ic,jc,1)  &
            + strwzr(ic,jc,1)/acouzr(ic,jc) + ovgmzr(ic,jc))  &
            - bcl3zr(ic,jc,1)*strdzr(ic,jc,1)*struzr(ic,jc,1)  &
            - bcl4zr(ic,jc,1)*strdzr(ic,jc,1)*strvzr(ic,jc,1)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)  &
              - bcl1zr(ic,jc,1)*ova2zr(ic,jc)*stryzr(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
  IF(nsbczr == nsbcw2)THEN
    
!         WALL BOUNDARY CONDITION No 2
!         NO-SLIP WALL - ISOTHERMAL
    
!         VELOCITY AND TEMPERATURE IMPOSED
!         AS FUNCTIONS OF TIME
!         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
!         SET IN SUBROUTINE BOUNDT
    
!         PRECOMPUTE CHEMISTRY TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzr(ic,jc) = zero
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          sorpzr(ic,jc) = sorpzr(ic,jc)  &
              + strhzr(ic,jc,ispec)*ratezr(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        sorpzr(ic,jc) = -sorpzr(ic,jc)*gam1zr(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Z-L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzr(ic,jc,1)*acouzr(ic,jc)*bcl1zr(ic,jc,1)
        bcl1zr(ic,jc,1) = half*(strwzr(ic,jc,1)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc,1)-fornow)
        bcl2zr(ic,jc,1) = strwzr(ic,jc,1)  &
            *(bcl2zr(ic,jc,1)-bcl5zr(ic,jc,1)*ova2zr(ic,jc))
        bcl3zr(ic,jc,1) = strwzr(ic,jc,1)*bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = strwzr(ic,jc,1)*bcl4zr(ic,jc,1)
        bcl5zr(ic,jc,1) = half*(strwzr(ic,jc,1)+acouzr(ic,jc))  &
            *(bcl5zr(ic,jc,1)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Z UNCHANGED
        bcl1zr(ic,jc,1) = bcl5zr(ic,jc,1)  &
            + strdzr(ic,jc,1)*acouzr(ic,jc)*dwdtzr(ic,jc) - bcl1zr(ic,jc,1)
        bcl3zr(ic,jc,1) = -dudtzr(ic,jc) - bcl3zr(ic,jc,1)
        bcl4zr(ic,jc,1) = -dvdtzr(ic,jc) - bcl4zr(ic,jc,1)
        bcl2zr(ic,jc,1) = gam1zr(ic,jc)*ova2zr(ic,jc)  &
            *(bcl1zr(ic,jc,1)+bcl5zr(ic,jc,1))  &
            + strdzr(ic,jc,1)*(dtdtzr(ic,jc)/strtzr(ic,jc,1)  &
            - sorpzr(ic,jc)/strpzr(ic,jc,1)) - bcl2zr(ic,jc,1)
        
      END DO
    END DO
    
!         LYZ
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYZ
          bclyzr(ic,jc,ispec) = strwzr(ic,jc,1)*bclyzr(ic,jc,ispec)
          
!               UPDATE L2Z
          bcl2zr(ic,jc,1) = bcl2zr(ic,jc,1) + (ratezr(ic,jc,ispec)  &
              - strdzr(ic,jc,1)*bclyzr(ic,jc,ispec)) *rgspec(ispec)/strrzr(ic,jc,1)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc,1)*ova2zr(ic,jc)  &
            - bcl2zr(ic,jc,1)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)  &
              - (bcl2zr(ic,jc,1)+bcl1zr(ic,jc,1)*ova2zr(ic,jc))*stryzr(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
  END IF
  
!       =======================================================================
  
END IF
!     Z-DIRECTION RIGHT-HAND END

!     =========================================================================


RETURN
END SUBROUTINE bounds
