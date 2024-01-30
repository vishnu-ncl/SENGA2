SUBROUTINE bounds

! Code converted using TO_F90 by Alan Miller
! Date: 2024-01-30  Time: 13:43:22

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

! VISHNU MOHAN (VM), KHALIL ABO AMSHA (KAA) AND NILANJAN CHAKRABORTY (NC)
! IMPLEMENTED THE TRANSVERSE DERIVATIVE TERMS FOR BOUNDARY CONDITION AS PER METHOD
! OF YOO AND IM


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
REAL(kind=8) :: fornow
INTEGER :: ic,jc,kc
INTEGER :: ispec
!     THESE ARE INTRODUCED FOR LODATO'S BC EQN. 3.80- NC
DOUBLE PRECISION :: tt1xl(nysize,nzsize), tt2xl(nysize,nzsize),  &
    tt3xl(nysize,nzsize), tt4xl(nysize,nzsize),  &
    tt5xl(nysize,nzsize), tt6xl(nysize,nzsize,nspec)
DOUBLE PRECISION :: tt1xr(nysize,nzsize), tt2xr(nysize,nzsize),  &
    tt3xr(nysize,nzsize), tt4xr(nysize,nzsize),  &
    tt5xr(nysize,nzsize), tt6xr(nysize,nzsize,nspec)
!UMAIRS CORRECTION HERE
DOUBLE PRECISION :: tt1yl(nxsize,nzsize), tt2yl(nxsize,nzsize),  &
    tt3yl(nxsize,nzsize), tt4yl(nxsize,nzsize),  &
    tt5yl(nxsize,nzsize), tt6yl(nxsize,nzsize,nspec)
DOUBLE PRECISION :: tt1yr(nxsize,nzsize), tt2yr(nxsize,nzsize),  &
    tt3yr(nxsize,nzsize), tt4yr(nxsize,nzsize),  &
    tt5yr(nxsize,nzsize), tt6yr(nxsize,nzsize,nspec)
DOUBLE PRECISION :: tt1zl(nxsize,nysize), tt2zl(nxsize,nysize),  &
    tt3zl(nxsize,nysize), tt4zl(nxsize,nysize),  &
    tt5zl(nxsize,nysize), tt6zl(nxsize,nysize,nspec)
DOUBLE PRECISION :: tt1zr(nxsize,nysize), tt2zr(nxsize,nysize),  &
    tt3zr(nxsize,nysize), tt4zr(nxsize,nysize),  &
    tt5zr(nxsize,nysize), tt6zr(nxsize,nysize,nspec)
DOUBLE PRECISION :: bet
INTEGER :: flag_bet_xl,flag_bet_xr,flag_bet_yl,flag_bet_yr,  &
    flag_bet_zl,flag_bet_zr
!VM: SWITCH TO SET BET=LOCAL_MACH
!FLAG_BET=0 POINTS TO GLOBAL VALUE OF BET
!FLAG_BET=1 POINTS TO LOCAL VALUE OF BET

INTEGER :: flag_pio_xl,flag_pio_xr,flag_pio_yl,flag_pio_yr,  &
    flag_pio_zl,flag_pio_zr
!VM: SWITCH TO TURN POINTWISE INFLOW-OUTFLOW ON/OFF
!1 IS ON, 0 IS OFF
bet=0.050D0

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
            - strgxl(jc,kc)*strtxl(jc,kc)*rgspec(ispec)/strrxl(jc,kc)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      gam1xl(jc,kc) = strgxl(jc,kc) - strrxl(jc,kc)
      strexl(jc,kc) = strexl(jc,kc) - gam1xl(jc,kc)*strtxl(jc,kc)
      
      gam1xl(jc,kc) = strrxl(jc,kc)/gam1xl(jc,kc)
      ovgmxl(jc,kc) = one/gam1xl(jc,kc)
      
      
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      fornow = strgxl(jc,kc)*gam1xl(jc,kc)*strtxl(jc,kc)
      acouxl(jc,kc) = SQRT(fornow)
      ova2xl(jc,kc) = one/fornow
      
!     THESE ARE INTRODUCED FOR LODATO'S BC- NC
      
      tt1xl(jc,kc)=t51bxl(jc,kc)+ t52bxl(jc,kc)*(gam1xl(jc,kc)+1.0)-  &
          strdxl(jc,kc)* acouxl(jc,kc)*t2bxl(jc,kc)
      
      tt2xl(jc,kc)=acouxl(jc,kc)*acouxl(jc,kc)*  &
          t1bxl(jc,kc)-t51bxl(jc,kc)-(gam1xl(jc,kc)+1.0)* t52bxl(jc,kc)
      
      tt3xl(jc,kc)=t3bxl(jc,kc)
      tt4xl(jc,kc)=t4bxl(jc,kc)
      
      tt5xl(jc,kc)=t51bxl(jc,kc)+ t52bxl(jc,kc)*(gam1xl(jc,kc)+1.0)+  &
          strdxl(jc,kc)* acouxl(jc,kc)*t2bxl(jc,kc)
      
      DO is=1,nspec
        tt6xl(jc,kc,is)=t6bxl(jc,kc,is)
      END DO
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbcxl == nsbco1)THEN
    flag_pio_xl=nxlprm(1)!0
    flag_bet_xl=nxlprm(2)!1
    
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
        bcl2xl(jc,kc) = struxl(jc,kc)  &
            *(bcl2xl(jc,kc)-bcl5xl(jc,kc)*ova2xl(jc,kc))
        bcl3xl(jc,kc) = struxl(jc,kc)*bcl3xl(jc,kc)
        bcl4xl(jc,kc) = struxl(jc,kc)*bcl4xl(jc,kc)
!             OLD VALUE OF L5X
        bcl5xl(jc,kc) = half*(struxl(jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)+strdxl(jc,kc)*acouxl(jc,kc)*bcl1xl(jc,kc))
        
!             SUBTRACT FROM NEW VALUE OF L5X
!             TERM INTRODUCED FOR LODATO'S BC- NC
        IF(flag_bet_xl==1) THEN
          bet=struxl(jc,kc)*struxl(jc,kc)+ strvxl(jc,kc)*strvxl(jc,kc)+  &
              strwxl(jc,kc)*strwxl(jc,kc)
          bet=SQRT(bet)/acouxl(jc,kc)
        END IF
        IF((struxl(jc,kc) > zero).AND.(flag_pio_xl==1))THEN
          
!             SUBTRACT FROM NEW VALUE OF L5X
          bcl2xl(jc,kc) = -bcl2xl(jc,kc)  &
              -ova2xl(jc,kc)*sorpxl(jc,kc)
          
          bcl3xl(jc,kc) = 0.1*(strvxl(jc,kc)-0.0D0)  &
              -bcl3xl(jc,kc)
          
          bcl4xl(jc,kc) = 0.1*(strwxl(jc,kc)-0.0D0)  &
              -bcl4xl(jc,kc)
          bcl5xl(jc,kc)= half*sorpxl(jc,kc)  &
              + cobcxl*acouxl(jc,kc)*(strpxl(jc,kc)-pinfxl)  &
              +0.5*(1.0-bet)*tt5xl(jc,kc)- bcl5xl(jc,kc)  &
              + 100.0*strdxl(jc,kc)*(one/xgdlen)  &
              *(acouxl(jc,kc)**2.0-struxl(jc,kc)**2.0) *(struxl(jc,kc)-0.0)
        ELSE
          
          bcl5xl(jc,kc)= half*sorpxl(jc,kc)  &
              + cobcxl*acouxl(jc,kc)*(strpxl(jc,kc)-pinfxl)  &
              +0.5*(1.0-bet)*tt5xl(jc,kc)- bcl5xl(jc,kc)
        END IF
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        IF((struxl(jc,kc) > zero).AND.(flag_pio_xl==1))THEN
          
          drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(jc,kc)  &
              - bcl5xl(jc,kc)*ova2xl(jc,kc)
          
          urhs(istal,jc,kc) = urhs(istal,jc,kc)  &
              - bcl2xl(jc,kc)*struxl(jc,kc)  &
              - bcl5xl(jc,kc)*ova2xl(jc,kc)*(struxl(jc,kc)+acouxl(jc,kc))
          
          vrhs(istal,jc,kc) = vrhs(istal,jc,kc)  &
              - bcl2xl(jc,kc)*strvxl(jc,kc) - bcl3xl(jc,kc)*strdxl(jc,kc)  &
              - bcl5xl(jc,kc)*ova2xl(jc,kc)*strvxl(jc,kc)
          
          wrhs(istal,jc,kc) = wrhs(istal,jc,kc)  &
              - bcl2xl(jc,kc)*strwxl(jc,kc) - bcl4xl(jc,kc)*strdxl(jc,kc)  &
              - bcl5xl(jc,kc)*ova2xl(jc,kc)*strwxl(jc,kc)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
              - bcl2xl(jc,kc)*strexl(jc,kc)  &
              - bcl3xl(jc,kc)*strdxl(jc,kc)*strvxl(jc,kc)  &
              - bcl4xl(jc,kc)*strdxl(jc,kc)*strwxl(jc,kc)  &
              - bcl5xl(jc,kc)*(ova2xl(jc,kc)*strexl(jc,kc)  &
              + struxl(jc,kc)/acouxl(jc,kc) + ovgmxl(jc,kc))
          
        ELSE
          drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl5xl(jc,kc)*ova2xl(jc,kc)
          
          urhs(istal,jc,kc) = urhs(istal,jc,kc)  &
              - bcl5xl(jc,kc)*ova2xl(jc,kc)*(struxl(jc,kc)+acouxl(jc,kc))
          
          vrhs(istal,jc,kc) = vrhs(istal,jc,kc)  &
              - bcl5xl(jc,kc)*ova2xl(jc,kc)*strvxl(jc,kc)
          
          wrhs(istal,jc,kc) = wrhs(istal,jc,kc)  &
              - bcl5xl(jc,kc)*ova2xl(jc,kc)*strwxl(jc,kc)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
              - bcl5xl(jc,kc)*(ova2xl(jc,kc)*strexl(jc,kc)  &
              + struxl(jc,kc)/acouxl(jc,kc) + ovgmxl(jc,kc))
        END IF
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          IF((struxl(jc,kc) > zero).AND.(flag_pio_xl==1))THEN
            fornow = bclyxl(jc,kc,ispec)*strdxl(jc,kc)
            
            erhs(istal,jc,kc) = erhs(istal,jc,kc) - fornow*strhxl(jc,kc,ispec)
            
            yrhs(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)  &
                - (bcl2xl(jc,kc)+bcl5xl(jc,kc)*ova2xl(jc,kc))*stryxl(jc,kc,ispec)  &
                - fornow
            
          ELSE
            yrhs(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)  &
                - bcl5xl(jc,kc)*ova2xl(jc,kc)*stryxl(jc,kc,ispec)
          END IF
          
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
        fornow = strdxl(jc,kc)*acouxl(jc,kc)*bcl1xl(jc,kc)
        bcl2xl(jc,kc) = struxl(jc,kc)  &
            *(bcl2xl(jc,kc)-bcl5xl(jc,kc)*ova2xl(jc,kc))
        bcl3xl(jc,kc) = struxl(jc,kc)*bcl3xl(jc,kc)
        bcl4xl(jc,kc) = struxl(jc,kc)*bcl4xl(jc,kc)
        bcl5xl(jc,kc) = half*(struxl(jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
!             L1X UNCHANGED
        bcl2xl(jc,kc) = -bcl2xl(jc,kc)
        bcl3xl(jc,kc) = -bcl3xl(jc,kc)
        bcl4xl(jc,kc) = -bcl4xl(jc,kc)
        bcl5xl(jc,kc) = half*sorpxl(jc,kc)  &
            + cobcxl*acouxl(jc,kc)*(strpxl(jc,kc)-pinfxl) - bcl5xl(jc,kc)
        
      END DO
    END DO
    
!         LYX
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
!               OLD VALUE OF L's
          bclyxl(jc,kc,ispec) = struxl(jc,kc)*bclyxl(jc,kc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
          bclyxl(jc,kc,ispec) = ratexl(jc,kc,ispec)/strdxl(jc,kc)  &
              - bclyxl(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(jc,kc)  &
            - bcl5xl(jc,kc)*ova2xl(jc,kc)
        
        urhs(istal,jc,kc) = urhs(istal,jc,kc) - bcl2xl(jc,kc)*struxl(jc,kc)  &
            - bcl5xl(jc,kc)*ova2xl(jc,kc)*(struxl(jc,kc)+acouxl(jc,kc))
        
        vrhs(istal,jc,kc) = vrhs(istal,jc,kc) - bcl2xl(jc,kc)*strvxl(jc,kc)  &
            - bcl3xl(jc,kc)*strdxl(jc,kc)  &
            - bcl5xl(jc,kc)*ova2xl(jc,kc)*strvxl(jc,kc)
        
        wrhs(istal,jc,kc) = wrhs(istal,jc,kc) - bcl2xl(jc,kc)*strwxl(jc,kc)  &
            - bcl4xl(jc,kc)*strdxl(jc,kc)  &
            - bcl5xl(jc,kc)*ova2xl(jc,kc)*strwxl(jc,kc)
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc) - bcl2xl(jc,kc)*strexl(jc,kc)  &
            - bcl3xl(jc,kc)*strdxl(jc,kc)*strvxl(jc,kc)  &
            - bcl4xl(jc,kc)*strdxl(jc,kc)*strwxl(jc,kc)  &
            - bcl5xl(jc,kc)*(ova2xl(jc,kc)*strexl(jc,kc)  &
            + struxl(jc,kc)/acouxl(jc,kc) + ovgmxl(jc,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = bclyxl(jc,kc,ispec)*strdxl(jc,kc)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc) - fornow*strhxl(jc,kc,ispec)
          
          yrhs(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)  &
              - (bcl2xl(jc,kc)+bcl5xl(jc,kc)*ova2xl(jc,kc))*stryxl(jc,kc,ispec)  &
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
        
        sydtxl(jc,kc) = sydtxl(jc,kc)/strrxl(jc,kc)
        sorpxl(jc,kc) = -sorpxl(jc,kc)*gam1xl(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1X,L2X,L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxl(jc,kc)*acouxl(jc,kc)*bcl1xl(jc,kc)
        bcl1xl(jc,kc) = half*(struxl(jc,kc)-acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)-fornow)
        bcl2xl(jc,kc) = struxl(jc,kc)  &
            *(bcl2xl(jc,kc)-bcl5xl(jc,kc)*ova2xl(jc,kc))
        bcl5xl(jc,kc) = half*(struxl(jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1X UNCHANGED
        bcl5xl(jc,kc) = bcl1xl(jc,kc)  &
            - strdxl(jc,kc)*acouxl(jc,kc)*dudtxl(jc,kc) - bcl5xl(jc,kc)
        bcl2xl(jc,kc) = gam1xl(jc,kc)*ova2xl(jc,kc)  &
            *(bcl1xl(jc,kc)+bcl5xl(jc,kc))  &
            + strdxl(jc,kc)*(dtdtxl(jc,kc)/strtxl(jc,kc)  &
            - sorpxl(jc,kc)/strpxl(jc,kc) + sydtxl(jc,kc))  &
            - bcl2xl(jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(jc,kc)  &
            - bcl5xl(jc,kc)*ova2xl(jc,kc)
        
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
        fornow = strdxl(jc,kc)*acouxl(jc,kc)*bcl1xl(jc,kc)
        bcl1xl(jc,kc) = half*(struxl(jc,kc)-acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)-fornow)
        bcl2xl(jc,kc) = struxl(jc,kc)  &
            *(bcl2xl(jc,kc)-bcl5xl(jc,kc)*ova2xl(jc,kc))
        bcl3xl(jc,kc) = struxl(jc,kc)*bcl3xl(jc,kc)
        bcl4xl(jc,kc) = struxl(jc,kc)*bcl4xl(jc,kc)
        bcl5xl(jc,kc) = half*(struxl(jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1X UNCHANGED
        fornow = bcl1xl(jc,kc) - strdxl(jc,kc)*acouxl(jc,kc)*dudtxl(jc,kc)
        bcl2xl(jc,kc) = -dddtxl(jc,kc)  &
            - ova2xl(jc,kc)*(bcl1xl(jc,kc)+fornow) - bcl2xl(jc,kc)
        bcl3xl(jc,kc) = -dvdtxl(jc,kc) - bcl3xl(jc,kc)
        bcl4xl(jc,kc) = -dwdtxl(jc,kc) - bcl4xl(jc,kc)
        bcl5xl(jc,kc) = fornow - bcl5xl(jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc) - bcl2xl(jc,kc)*strexl(jc,kc)  &
            - bcl3xl(jc,kc)*strdxl(jc,kc)*strvxl(jc,kc)  &
            - bcl4xl(jc,kc)*strdxl(jc,kc)*strwxl(jc,kc)  &
            - bcl5xl(jc,kc)*(ova2xl(jc,kc)*strexl(jc,kc)  &
            + struxl(jc,kc)/acouxl(jc,kc) + ovgmxl(jc,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          bclyxl(jc,kc,ispec) = ratexl(jc,kc,ispec)/strdxl(jc,kc)  &
              - dydtxl(jc,kc,ispec) - struxl(jc,kc)*bclyxl(jc,kc,ispec)
          
          erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
              - bclyxl(jc,kc,ispec)*strdxl(jc,kc)*strhxl(jc,kc,ispec)
          
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
        fornow = strdxl(jc,kc)*acouxl(jc,kc)*bcl1xl(jc,kc)
        bcl1xl(jc,kc) = half*(struxl(jc,kc)-acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)-fornow)
        bcl3xl(jc,kc) = struxl(jc,kc)*bcl3xl(jc,kc)
        bcl4xl(jc,kc) = struxl(jc,kc)*bcl4xl(jc,kc)
        bcl5xl(jc,kc) = half*(struxl(jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1X,L2X UNCHANGED
        bcl3xl(jc,kc) = -dvdtxl(jc,kc) - bcl3xl(jc,kc)
        bcl4xl(jc,kc) = -dwdtxl(jc,kc) - bcl4xl(jc,kc)
        bcl5xl(jc,kc) = bcl1xl(jc,kc)  &
            - strdxl(jc,kc)*acouxl(jc,kc)*dudtxl(jc,kc) - bcl5xl(jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl5xl(jc,kc)*ova2xl(jc,kc)
        
        erhs(istal,jc,kc) = erhs(istal,jc,kc)  &
            - bcl3xl(jc,kc)*strdxl(jc,kc)*strvxl(jc,kc)  &
            - bcl4xl(jc,kc)*strdxl(jc,kc)*strwxl(jc,kc)  &
            - bcl5xl(jc,kc)*(ova2xl(jc,kc)*strexl(jc,kc)  &
            + struxl(jc,kc)/acouxl(jc,kc) + ovgmxl(jc,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)  &
              - bcl5xl(jc,kc)*ova2xl(jc,kc)*stryxl(jc,kc,ispec)
          
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
        fornow = strdxl(jc,kc)*acouxl(jc,kc)*bcl1xl(jc,kc)
        bcl1xl(jc,kc) = half*(struxl(jc,kc)-acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)-fornow)
        bcl2xl(jc,kc) = struxl(jc,kc)  &
            *(bcl2xl(jc,kc)-bcl5xl(jc,kc)*ova2xl(jc,kc))
        bcl3xl(jc,kc) = struxl(jc,kc)*bcl3xl(jc,kc)
        bcl4xl(jc,kc) = struxl(jc,kc)*bcl4xl(jc,kc)
        bcl5xl(jc,kc) = half*(struxl(jc,kc)+acouxl(jc,kc))  &
            *(bcl5xl(jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1X UNCHANGED
        bcl3xl(jc,kc) = -dvdtxl(jc,kc) - bcl3xl(jc,kc)
        bcl4xl(jc,kc) = -dwdtxl(jc,kc) - bcl4xl(jc,kc)
        bcl5xl(jc,kc) = bcl1xl(jc,kc)  &
            - strdxl(jc,kc)*acouxl(jc,kc)*dudtxl(jc,kc) - bcl5xl(jc,kc)
        bcl2xl(jc,kc) = gam1xl(jc,kc)*ova2xl(jc,kc)  &
            *(bcl1xl(jc,kc)+bcl5xl(jc,kc))  &
            + strdxl(jc,kc)*(dtdtxl(jc,kc)/strtxl(jc,kc)  &
            - sorpxl(jc,kc)/strpxl(jc,kc)) - bcl2xl(jc,kc)
        
      END DO
    END DO
    
!         LYX
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
!               OLD VALUE OF LYX
          bclyxl(jc,kc,ispec) = struxl(jc,kc)*bclyxl(jc,kc,ispec)
          
!               UPDATE L2X
          bcl2xl(jc,kc) = bcl2xl(jc,kc) + (ratexl(jc,kc,ispec)  &
              - strdxl(jc,kc)*bclyxl(jc,kc,ispec)) *rgspec(ispec)/strrxl(jc,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istal,jc,kc) = drhs(istal,jc,kc) - bcl2xl(jc,kc)  &
            - bcl5xl(jc,kc)*ova2xl(jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istal,jc,kc,ispec) = yrhs(istal,jc,kc,ispec)  &
              - (bcl2xl(jc,kc)+bcl5xl(jc,kc)*ova2xl(jc,kc))*stryxl(jc,kc,ispec)
          
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
            - strgxr(jc,kc)*strtxr(jc,kc)*rgspec(ispec)/strrxr(jc,kc)
        
      END DO
      
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      gam1xr(jc,kc) = strgxr(jc,kc) - strrxr(jc,kc)
      strexr(jc,kc) = strexr(jc,kc) - gam1xr(jc,kc)*strtxr(jc,kc)
      
      gam1xr(jc,kc) = strrxr(jc,kc)/gam1xr(jc,kc)
      ovgmxr(jc,kc) = one/gam1xr(jc,kc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
      fornow = strgxr(jc,kc)*gam1xr(jc,kc)*strtxr(jc,kc)
      acouxr(jc,kc) = SQRT(fornow)
      ova2xr(jc,kc) = one/fornow
      
!     THESE ARE INTRODUCED FOR LODATO'S BC- NC
      
      tt1xr(jc,kc)=t51bxr(jc,kc)+ t52bxr(jc,kc)*(gam1xr(jc,kc)+1.0)-  &
          strdxr(jc,kc)* acouxr(jc,kc)*t2bxr(jc,kc)
      
      tt2xr(jc,kc)=acouxr(jc,kc)*acouxr(jc,kc)*  &
          t1bxr(jc,kc)-t51bxr(jc,kc)-(gam1xr(jc,kc)+1.0)* t52bxr(jc,kc)
      
      tt3xr(jc,kc)=t3bxr(jc,kc)
      tt4xr(jc,kc)=t4bxr(jc,kc)
      
      tt5xr(jc,kc)=t51bxr(jc,kc)+ t52bxr(jc,kc)*(gam1xr(jc,kc)+1.0)+  &
          strdxr(jc,kc)* acouxr(jc,kc)*t2bxr(jc,kc)
      
      DO is=1,nspec
        tt6xr(jc,kc,is)=t6bxr(jc,kc,is)
      END DO
      
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbcxr == nsbco1)THEN
    flag_pio_xr=nxrprm(1)!0
    flag_bet_xr=nxrprm(2)!1
    
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
        bcl2xr(jc,kc) = struxr(jc,kc)*(bcl2xr(jc,kc)-ova2xr(jc,kc)  &
            *bcl5xr(jc,kc))
        
        bcl3xr(jc,kc) = struxr(jc,kc)*bcl3xr(jc,kc)
        bcl4xr(jc,kc) = struxr(jc,kc)*bcl4xr(jc,kc)
        
!             OLD VALUE OF L1X
        bcl1xr(jc,kc) = half*(struxr(jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)-strdxr(jc,kc)*acouxr(jc,kc)*bcl1xr(jc,kc))
        
!             SUBTRACT FROM NEW VALUE OF L1X
!     TERMS INTRODUCED FOR LODATO'S BC- NC
        IF(flag_bet_xr==1) THEN
          bet=struxr(jc,kc)*struxr(jc,kc)+ strvxr(jc,kc)*strvxr(jc,kc)+  &
              strwxr(jc,kc)*strwxr(jc,kc)
          bet=SQRT(bet)/acouxr(jc,kc)
        END IF
        IF((struxr(jc,kc) < zero).AND.(flag_pio_xr==1))THEN
          bcl2xr(jc,kc)=-ova2xr(jc,kc)*sorpxr(jc,kc) -bcl2xr(jc,kc)
          
          bcl3xr(jc,kc)=0.1*(strvxr(jc,kc)-0.0D0)-bcl3xr(jc,kc)
          
          bcl4xr(jc,kc)=0.1*(strwxr(jc,kc)-0.0D0)-bcl4xr(jc,kc)
          bcl1xr(jc,kc)= half*sorpxr(jc,kc)  &
              + cobcxr*acouxr(jc,kc)*(strpxr(jc,kc)-pinfxr)  &
              +0.5*(1.0-bet)*tt1xr(jc,kc)- bcl1xr(jc,kc)  &
              -100.0*strdxr(jc,kc)*(one/xgdlen)  &
              *(acouxr(jc,kc)**2.0-struxr(jc,kc)**2.0) *(struxr(jc,kc)-0.0D0)
        ELSE
          bcl1xr(jc,kc)= half*sorpxr(jc,kc)  &
              + cobcxr*acouxr(jc,kc)*(strpxr(jc,kc)-pinfxr)  &
              +0.5*(1.0-bet)*tt1xr(jc,kc)- bcl1xr(jc,kc)
        END IF
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        IF((struxr(jc,kc) < zero).AND.(flag_pio_xr==1))THEN
          drhs(istol,jc,kc) = drhs(istol,jc,kc)  &
              - bcl1xr(jc,kc)*ova2xr(jc,kc) - bcl2xr(jc,kc)
          
          urhs(istol,jc,kc) = urhs(istol,jc,kc)  &
              - bcl1xr(jc,kc)*ova2xr(jc,kc)*(struxr(jc,kc)-acouxr(jc,kc))  &
              - bcl2xr(jc,kc)*struxr(jc,kc)
          
          vrhs(istol,jc,kc) = vrhs(istol,jc,kc)  &
              - bcl1xr(jc,kc)*ova2xr(jc,kc)*strvxr(jc,kc)  &
              - bcl2xr(jc,kc)*strvxr(jc,kc) - bcl3xr(jc,kc)*strdxr(jc,kc)
          
          wrhs(istol,jc,kc) = wrhs(istol,jc,kc)  &
              - bcl1xr(jc,kc)*ova2xr(jc,kc)*strwxr(jc,kc)  &
              - bcl2xr(jc,kc)*strwxr(jc,kc) - bcl4xr(jc,kc)*strdxr(jc,kc)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
              - bcl1xr(jc,kc)*(ova2xr(jc,kc)*strexr(jc,kc)  &
              - struxr(jc,kc)/acouxr(jc,kc) + ovgmxr(jc,kc))  &
              - bcl2xr(jc,kc)*strexr(jc,kc)  &
              - bcl3xr(jc,kc)*strdxr(jc,kc)*strvxr(jc,kc)  &
              - bcl4xr(jc,kc)*strdxr(jc,kc)*strwxr(jc,kc)
          
          
        ELSE
          
          drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(jc,kc)*ova2xr(jc,kc)
          
          urhs(istol,jc,kc) = urhs(istol,jc,kc)  &
              - bcl1xr(jc,kc)*ova2xr(jc,kc)*(struxr(jc,kc)-acouxr(jc,kc))
          
          vrhs(istol,jc,kc) = vrhs(istol,jc,kc)  &
              - bcl1xr(jc,kc)*ova2xr(jc,kc)*strvxr(jc,kc)
          
          wrhs(istol,jc,kc) = wrhs(istol,jc,kc)  &
              - bcl1xr(jc,kc)*ova2xr(jc,kc)*strwxr(jc,kc)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
              - bcl1xr(jc,kc)*(ova2xr(jc,kc)*strexr(jc,kc)  &
              - struxr(jc,kc)/acouxr(jc,kc) + ovgmxr(jc,kc))
        END IF
        
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          IF((struxr(jc,kc) < zero).AND.(flag_pio_xr==1))THEN
            
            fornow = bclyxr(jc,kc,ispec)*strdxr(jc,kc)
            
            erhs(istol,jc,kc) = erhs(istol,jc,kc) - fornow*strhxr(jc,kc,ispec)
            
            yrhs(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)  &
                - (bcl2xr(jc,kc)+bcl1xr(jc,kc)*ova2xr(jc,kc))*stryxr(jc,kc,ispec)  &
                - fornow
            
          ELSE
            yrhs(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)  &
                - bcl1xr(jc,kc)*ova2xr(jc,kc)*stryxr(jc,kc,ispec)
          END IF
          
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
        fornow = strdxr(jc,kc)*acouxr(jc,kc)*bcl1xr(jc,kc)
        bcl1xr(jc,kc) = half*(struxr(jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)-fornow)
        bcl2xr(jc,kc) = struxr(jc,kc)  &
            *(bcl2xr(jc,kc)-bcl5xr(jc,kc)*ova2xr(jc,kc))
        bcl3xr(jc,kc) = struxr(jc,kc)*bcl3xr(jc,kc)
        bcl4xr(jc,kc) = struxr(jc,kc)*bcl4xr(jc,kc)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
!             L5X UNCHANGED
        bcl1xr(jc,kc) = half*sorpxr(jc,kc)  &
            + cobcxr*acouxr(jc,kc)*(strpxr(jc,kc)-pinfxr) - bcl1xr(jc,kc)
        bcl2xr(jc,kc) = -bcl2xr(jc,kc)
        bcl3xr(jc,kc) = -bcl3xr(jc,kc)
        bcl4xr(jc,kc) = -bcl4xr(jc,kc)
        
      END DO
    END DO
    
!         LYX
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
!               OLD VALUE OF L's
          bclyxr(jc,kc,ispec) = struxr(jc,kc)*bclyxr(jc,kc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
          bclyxr(jc,kc,ispec) = ratexr(jc,kc,ispec)/strdxr(jc,kc)  &
              - bclyxr(jc,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(jc,kc)*ova2xr(jc,kc)  &
            - bcl2xr(jc,kc)
        
        urhs(istol,jc,kc) = urhs(istol,jc,kc)  &
            - bcl1xr(jc,kc)*ova2xr(jc,kc)*(struxr(jc,kc)-acouxr(jc,kc))  &
            - bcl2xr(jc,kc)*struxr(jc,kc)
        
        vrhs(istol,jc,kc) = vrhs(istol,jc,kc)  &
            - bcl1xr(jc,kc)*ova2xr(jc,kc)*strvxr(jc,kc)  &
            - bcl2xr(jc,kc)*strvxr(jc,kc) - bcl3xr(jc,kc)*strdxr(jc,kc)
        
        wrhs(istol,jc,kc) = wrhs(istol,jc,kc)  &
            - bcl1xr(jc,kc)*ova2xr(jc,kc)*strwxr(jc,kc)  &
            - bcl2xr(jc,kc)*strwxr(jc,kc) - bcl4xr(jc,kc)*strdxr(jc,kc)
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(jc,kc)*(ova2xr(jc,kc)*strexr(jc,kc)  &
            - struxr(jc,kc)/acouxr(jc,kc) + ovgmxr(jc,kc))  &
            - bcl2xr(jc,kc)*strexr(jc,kc)  &
            - bcl3xr(jc,kc)*strdxr(jc,kc)*strvxr(jc,kc)  &
            - bcl4xr(jc,kc)*strdxr(jc,kc)*strwxr(jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          fornow = bclyxr(jc,kc,ispec)*strdxr(jc,kc)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc) - fornow*strhxr(jc,kc,ispec)
          
          yrhs(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)  &
              - (bcl2xr(jc,kc)+bcl1xr(jc,kc)*ova2xr(jc,kc))*stryxr(jc,kc,ispec)  &
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
        
        sydtxr(jc,kc) = sydtxr(jc,kc)/strrxr(jc,kc)
        sorpxr(jc,kc) = -sorpxr(jc,kc)*gam1xr(jc,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1X,L2X,L5X
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
!             OLD VALUE OF L's
        fornow = strdxr(jc,kc)*acouxr(jc,kc)*bcl1xr(jc,kc)
        bcl1xr(jc,kc) = half*(struxr(jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)-fornow)
        bcl2xr(jc,kc) = struxr(jc,kc)  &
            *(bcl2xr(jc,kc)-bcl5xr(jc,kc)*ova2xr(jc,kc))
        bcl5xr(jc,kc) = half*(struxr(jc,kc)+acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5X UNCHANGED
        bcl1xr(jc,kc) = bcl5xr(jc,kc)  &
            + strdxr(jc,kc)*acouxr(jc,kc)*dudtxr(jc,kc) - bcl1xr(jc,kc)
        bcl2xr(jc,kc) = gam1xr(jc,kc)*ova2xr(jc,kc)  &
            *(bcl1xr(jc,kc)+bcl5xr(jc,kc))  &
            + strdxr(jc,kc)*(dtdtxr(jc,kc)/strtxr(jc,kc)  &
            - sorpxr(jc,kc)/strpxr(jc,kc) + sydtxr(jc,kc))  &
            - bcl2xr(jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl2xr(jc,kc)  &
            - bcl1xr(jc,kc)*ova2xr(jc,kc)
        
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
        fornow = strdxr(jc,kc)*acouxr(jc,kc)*bcl1xr(jc,kc)
        bcl1xr(jc,kc) = half*(struxr(jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)-fornow)
        bcl2xr(jc,kc) = struxr(jc,kc)  &
            *(bcl2xr(jc,kc)-bcl5xr(jc,kc)*ova2xr(jc,kc))
        bcl3xr(jc,kc) = struxr(jc,kc)*bcl3xr(jc,kc)
        bcl4xr(jc,kc) = struxr(jc,kc)*bcl4xr(jc,kc)
        bcl5xr(jc,kc) = half*(struxr(jc,kc)+acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5X UNCHANGED
        fornow = bcl5xr(jc,kc) + strdxr(jc,kc)*acouxr(jc,kc)*dudtxr(jc,kc)
        bcl1xr(jc,kc) = fornow - bcl1xr(jc,kc)
        bcl2xr(jc,kc) = -dddtxr(jc,kc)  &
            - ova2xr(jc,kc)*(bcl5xr(jc,kc)+fornow) - bcl2xr(jc,kc)
        bcl3xr(jc,kc) = -dvdtxr(jc,kc) - bcl3xr(jc,kc)
        bcl4xr(jc,kc) = -dwdtxr(jc,kc) - bcl4xr(jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(jc,kc)*(ova2xr(jc,kc)*strexr(jc,kc)  &
            - struxr(jc,kc)/acouxr(jc,kc) + ovgmxr(jc,kc))  &
            - bcl2xr(jc,kc)*strexr(jc,kc)  &
            - bcl3xr(jc,kc)*strdxr(jc,kc)*strvxr(jc,kc)  &
            - bcl4xr(jc,kc)*strdxr(jc,kc)*strwxr(jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          bclyxr(jc,kc,ispec) = ratexr(jc,kc,ispec)/strdxr(jc,kc)  &
              - dydtxr(jc,kc,ispec) - struxr(jc,kc)*bclyxr(jc,kc,ispec)
          
          erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
              - bclyxr(jc,kc,ispec)*strdxr(jc,kc)*strhxr(jc,kc,ispec)
          
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
        fornow = strdxr(jc,kc)*acouxr(jc,kc)*bcl1xr(jc,kc)
        bcl1xr(jc,kc) = half*(struxr(jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)-fornow)
        bcl3xr(jc,kc) = struxr(jc,kc)*bcl3xr(jc,kc)
        bcl4xr(jc,kc) = struxr(jc,kc)*bcl4xr(jc,kc)
        bcl5xr(jc,kc) = half*(struxr(jc,kc)+acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L2X,L5X UNCHANGED
        bcl1xr(jc,kc) = bcl5xr(jc,kc)  &
            + strdxr(jc,kc)*acouxr(jc,kc)*dudtxr(jc,kc) - bcl1xr(jc,kc)
        bcl3xr(jc,kc) = -dvdtxr(jc,kc) - bcl3xr(jc,kc)
        bcl4xr(jc,kc) = -dwdtxr(jc,kc) - bcl4xr(jc,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(jc,kc)*ova2xr(jc,kc)
        
        erhs(istol,jc,kc) = erhs(istol,jc,kc)  &
            - bcl1xr(jc,kc)*(ova2xr(jc,kc)*strexr(jc,kc)  &
            + struxr(jc,kc)/acouxr(jc,kc) + ovgmxr(jc,kc))  &
            - bcl3xr(jc,kc)*strdxr(jc,kc)*strvxr(jc,kc)  &
            - bcl4xr(jc,kc)*strdxr(jc,kc)*strwxr(jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)  &
              - bcl1xr(jc,kc)*ova2xr(jc,kc)*stryxr(jc,kc,ispec)
          
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
        fornow = strdxr(jc,kc)*acouxr(jc,kc)*bcl1xr(jc,kc)
        bcl1xr(jc,kc) = half*(struxr(jc,kc)-acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)-fornow)
        bcl2xr(jc,kc) = struxr(jc,kc)  &
            *(bcl2xr(jc,kc)-bcl5xr(jc,kc)*ova2xr(jc,kc))
        bcl3xr(jc,kc) = struxr(jc,kc)*bcl3xr(jc,kc)
        bcl4xr(jc,kc) = struxr(jc,kc)*bcl4xr(jc,kc)
        bcl5xr(jc,kc) = half*(struxr(jc,kc)+acouxr(jc,kc))  &
            *(bcl5xr(jc,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5X UNCHANGED
        bcl1xr(jc,kc) = bcl5xr(jc,kc)  &
            + strdxr(jc,kc)*acouxr(jc,kc)*dudtxr(jc,kc) - bcl1xr(jc,kc)
        bcl3xr(jc,kc) = -dvdtxr(jc,kc) - bcl3xr(jc,kc)
        bcl4xr(jc,kc) = -dwdtxr(jc,kc) - bcl4xr(jc,kc)
        bcl2xr(jc,kc) = gam1xr(jc,kc)*ova2xr(jc,kc)  &
            *(bcl1xr(jc,kc)+bcl5xr(jc,kc))  &
            + strdxr(jc,kc)*(dtdtxr(jc,kc)/strtxr(jc,kc)  &
            - sorpxr(jc,kc)/strpxr(jc,kc)) - bcl2xr(jc,kc)
        
      END DO
    END DO
    
!         LYX
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
!               OLD VALUE OF LYX
          bclyxr(jc,kc,ispec) = struxr(jc,kc)*bclyxr(jc,kc,ispec)
          
!               UPDATE L2X
          bcl2xr(jc,kc) = bcl2xr(jc,kc) + (ratexr(jc,kc,ispec)  &
              - strdxr(jc,kc)*bclyxr(jc,kc,ispec)) *rgspec(ispec)/strrxr(jc,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        
        drhs(istol,jc,kc) = drhs(istol,jc,kc) - bcl1xr(jc,kc)*ova2xr(jc,kc)  &
            - bcl2xr(jc,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO jc = jstal,jstol
          
          yrhs(istol,jc,kc,ispec) = yrhs(istol,jc,kc,ispec)  &
              - (bcl2xr(jc,kc)+bcl1xr(jc,kc)*ova2xr(jc,kc))*stryxr(jc,kc,ispec)
          
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
            - strgyl(ic,kc)*strtyl(ic,kc)*rgspec(ispec)/strryl(ic,kc)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      gam1yl(ic,kc) = strgyl(ic,kc) - strryl(ic,kc)
      streyl(ic,kc) = streyl(ic,kc) - gam1yl(ic,kc)*strtyl(ic,kc)
      
      gam1yl(ic,kc) = strryl(ic,kc)/gam1yl(ic,kc)
      ovgmyl(ic,kc) = one/gam1yl(ic,kc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      fornow = strgyl(ic,kc)*gam1yl(ic,kc)*strtyl(ic,kc)
      acouyl(ic,kc) = SQRT(fornow)
      ova2yl(ic,kc) = one/fornow
      
!     THESE ARE INTRODUCED FOR LODATO'S BC- NC
      
      tt1yl(ic,kc)=t51byl(ic,kc)+ t52byl(ic,kc)*(gam1yl(ic,kc)+1.0)-  &
          strdyl(ic,kc)* acouyl(ic,kc)*t2byl(ic,kc)
      
      tt2yl(ic,kc)=acouyl(ic,kc)*acouyl(ic,kc)*  &
          t1byl(ic,kc)-t51byl(ic,kc)-(gam1yl(ic,kc)+1.0)* t52byl(ic,kc)
      
      tt3yl(ic,kc)=t3byl(ic,kc)
      tt4yl(ic,kc)=t4byl(ic,kc)
      
      tt5yl(ic,kc)=t51byl(ic,kc)+ t52byl(ic,kc)*(gam1yl(ic,kc)+1.0)+  &
          strdyl(ic,kc)* acouyl(ic,kc)*t2byl(ic,kc)
      
      DO is=1,nspec
        tt6yl(ic,kc,is)=t6byl(ic,kc,is)
      END DO
      
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbcyl == nsbco1)THEN
    flag_pio_yl=nylprm(1)!0
    flag_bet_yl=nylprm(2)!1
    
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
        bcl2yl(ic,kc) = strvyl(ic,kc)  &
            *(bcl2yl(ic,kc)-bcl5yl(ic,kc)*ova2yl(ic,kc))
        bcl3yl(ic,kc) = strvyl(ic,kc)*bcl3yl(ic,kc)
        bcl4yl(ic,kc) = strvyl(ic,kc)*bcl4yl(ic,kc)
!             OLD VALUE OF L5Y
        bcl5yl(ic,kc) = half*(strvyl(ic,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)+strdyl(ic,kc)*acouyl(ic,kc)*bcl1yl(ic,kc))
        
!             SUBTRACT FROM NEW VALUE OF L5Y
!             TERMS INTRODUCED FOR LODATO'S BC- NC
        IF(flag_bet_yl==1) THEN
          bet=struyl(ic,kc)*struyl(ic,kc)+ strvyl(ic,kc)*strvyl(ic,kc)+  &
              strwyl(ic,kc)*strwyl(ic,kc)
          bet=SQRT(bet)/acouyl(ic,kc)
        END IF
        IF((strvyl(ic,kc) > zero).AND.(flag_pio_yl==1))THEN
          bcl2yl(ic,kc) = -bcl2yl(ic,kc)  &
              -ova2yl(ic,kc)*sorpyl(ic,kc)
          
          bcl3yl(ic,kc) = 0.1*(struyl(ic,kc)-0.0D0)  &
              -bcl3yl(ic,kc)
          
          bcl4yl(ic,kc) = 0.1*(strwyl(ic,kc)-0.0D0)  &
              -bcl4yl(ic,kc)
          bcl5yl(ic,kc)= half*sorpyl(ic,kc)  &
              + cobcyl*acouyl(ic,kc)*(strpyl(ic,kc)-pinfyl)  &
              +0.5*(1.0-bet)*tt5yl(ic,kc)- bcl5yl(ic,kc)  &
              + 100.0*strdyl(ic,kc)*(one/ygdlen) &!ASK NC: CHANGED TO YGDLEN  &
              *(acouyl(ic,kc)**2.0-strvyl(ic,kc)**2.0)*(strvyl(ic,kc)-0.0)
          
        ELSE
          bcl5yl(ic,kc)= half*sorpyl(ic,kc)  &
              + cobcyl*acouyl(ic,kc)*(strpyl(ic,kc)-pinfyl)  &
              +0.5*(1.0-bet)*tt5yl(ic,kc)- bcl5yl(ic,kc)
        END IF
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        IF((strvyl(ic,kc) > zero).AND.(flag_pio_yl==1))THEN
          
          drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,kc)  &
              - bcl5yl(ic,kc)*ova2yl(ic,kc)
          
          urhs(ic,jstal,kc) = urhs(ic,jstal,kc)  &
              - bcl2yl(ic,kc)*struyl(ic,kc) - bcl3yl(ic,kc)*strdyl(ic,kc)  &
              - bcl5yl(ic,kc)*ova2yl(ic,kc)*struyl(ic,kc)
          
          vrhs(ic,jstal,kc) = vrhs(ic,jstal,kc)  &
              - bcl2yl(ic,kc)*strvyl(ic,kc)  &
              - bcl5yl(ic,kc)*ova2yl(ic,kc)*(strvyl(ic,kc)+acouyl(ic,kc))
          
          wrhs(ic,jstal,kc) = wrhs(ic,jstal,kc)  &
              - bcl2yl(ic,kc)*strwyl(ic,kc) - bcl4yl(ic,kc)*strdyl(ic,kc)  &
              - bcl5yl(ic,kc)*ova2yl(ic,kc)*strwyl(ic,kc)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
              - bcl2yl(ic,kc)*streyl(ic,kc)  &
              - bcl3yl(ic,kc)*strdyl(ic,kc)*struyl(ic,kc)  &
              - bcl4yl(ic,kc)*strdyl(ic,kc)*strwyl(ic,kc)  &
              - bcl5yl(ic,kc)*(ova2yl(ic,kc)*streyl(ic,kc)  &
              + strvyl(ic,kc)/acouyl(ic,kc) + ovgmyl(ic,kc))
          
        ELSE
          drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl5yl(ic,kc)*ova2yl(ic,kc)
          
          urhs(ic,jstal,kc) = urhs(ic,jstal,kc)  &
              - bcl5yl(ic,kc)*ova2yl(ic,kc)*struyl(ic,kc)
          
          vrhs(ic,jstal,kc) = vrhs(ic,jstal,kc)  &
              - bcl5yl(ic,kc)*ova2yl(ic,kc)*(strvyl(ic,kc)+acouyl(ic,kc))
          
          wrhs(ic,jstal,kc) = wrhs(ic,jstal,kc)  &
              - bcl5yl(ic,kc)*ova2yl(ic,kc)*strwyl(ic,kc)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
              - bcl5yl(ic,kc)*(ova2yl(ic,kc)*streyl(ic,kc)  &
              + strvyl(ic,kc)/acouyl(ic,kc) + ovgmyl(ic,kc))
        END IF
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          IF((strvyl(ic,kc) > zero).AND.(flag_pio_yl==1))THEN
            
            fornow = bclyyl(ic,kc,ispec)*strdyl(ic,kc)
            
            erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - fornow*strhyl(ic,kc,ispec)
            
            yrhs(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)  &
                - (bcl2yl(ic,kc)+bcl5yl(ic,kc)*ova2yl(ic,kc))*stryyl(ic,kc,ispec)  &
                - fornow
            
            
          ELSE
            yrhs(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)  &
                - bcl5yl(ic,kc)*ova2yl(ic,kc)*stryyl(ic,kc,ispec)
          END IF
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
        fornow = strdyl(ic,kc)*acouyl(ic,kc)*bcl1yl(ic,kc)
        bcl2yl(ic,kc) = strvyl(ic,kc)  &
            *(bcl2yl(ic,kc)-bcl5yl(ic,kc)*ova2yl(ic,kc))
        bcl3yl(ic,kc) = strvyl(ic,kc)*bcl3yl(ic,kc)
        bcl4yl(ic,kc) = strvyl(ic,kc)*bcl4yl(ic,kc)
        bcl5yl(ic,kc) = half*(strvyl(ic,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
!             L1Y UNCHANGED
        bcl2yl(ic,kc) = -bcl2yl(ic,kc)
        bcl3yl(ic,kc) = -bcl3yl(ic,kc)
        bcl4yl(ic,kc) = -bcl4yl(ic,kc)
        bcl5yl(ic,kc) = half*sorpyl(ic,kc)  &
            + cobcyl*acouyl(ic,kc)*(strpyl(ic,kc)-pinfyl) - bcl5yl(ic,kc)
        
      END DO
    END DO
    
!         LYY
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
!               OLD VALUE OF L's
          bclyyl(ic,kc,ispec) = strvyl(ic,kc)*bclyyl(ic,kc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
          bclyyl(ic,kc,ispec) = rateyl(ic,kc,ispec)/strdyl(ic,kc)  &
              - bclyyl(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,kc)  &
            - bcl5yl(ic,kc)*ova2yl(ic,kc)
        
        urhs(ic,jstal,kc) = urhs(ic,jstal,kc) - bcl2yl(ic,kc)*struyl(ic,kc)  &
            - bcl3yl(ic,kc)*strdyl(ic,kc)  &
            - bcl5yl(ic,kc)*ova2yl(ic,kc)*struyl(ic,kc)
        
        vrhs(ic,jstal,kc) = vrhs(ic,jstal,kc) - bcl2yl(ic,kc)*strvyl(ic,kc)  &
            - bcl5yl(ic,kc)*ova2yl(ic,kc)*(strvyl(ic,kc)+acouyl(ic,kc))
        
        wrhs(ic,jstal,kc) = wrhs(ic,jstal,kc) - bcl2yl(ic,kc)*strwyl(ic,kc)  &
            - bcl4yl(ic,kc)*strdyl(ic,kc)  &
            - bcl5yl(ic,kc)*ova2yl(ic,kc)*strwyl(ic,kc)
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - bcl2yl(ic,kc)*streyl(ic,kc)  &
            - bcl3yl(ic,kc)*strdyl(ic,kc)*struyl(ic,kc)  &
            - bcl4yl(ic,kc)*strdyl(ic,kc)*strwyl(ic,kc)  &
            - bcl5yl(ic,kc)*(ova2yl(ic,kc)*streyl(ic,kc)  &
            + strvyl(ic,kc)/acouyl(ic,kc) + ovgmyl(ic,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = bclyyl(ic,kc,ispec)*strdyl(ic,kc)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - fornow*strhyl(ic,kc,ispec)
          
          yrhs(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)  &
              - (bcl2yl(ic,kc)+bcl5yl(ic,kc)*ova2yl(ic,kc))*stryyl(ic,kc,ispec)  &
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
        
        sydtyl(ic,kc) = sydtyl(ic,kc)/strryl(ic,kc)
        sorpyl(ic,kc) = -sorpyl(ic,kc)*gam1yl(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Y,L2Y,L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyl(ic,kc)*acouyl(ic,kc)*bcl1yl(ic,kc)
        bcl1yl(ic,kc) = half*(strvyl(ic,kc)-acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)-fornow)
        bcl2yl(ic,kc) = strvyl(ic,kc)  &
            *(bcl2yl(ic,kc)-bcl5yl(ic,kc)*ova2yl(ic,kc))
        bcl5yl(ic,kc) = half*(strvyl(ic,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y UNCHANGED
        bcl5yl(ic,kc) = bcl1yl(ic,kc)  &
            - strdyl(ic,kc)*acouyl(ic,kc)*dvdtyl(ic,kc) - bcl5yl(ic,kc)
        bcl2yl(ic,kc) = gam1yl(ic,kc)*ova2yl(ic,kc)  &
            *(bcl1yl(ic,kc)+bcl5yl(ic,kc))  &
            + strdyl(ic,kc)*(dtdtyl(ic,kc)/strtyl(ic,kc)  &
            - sorpyl(ic,kc)/strpyl(ic,kc) + sydtyl(ic,kc))  &
            - bcl2yl(ic,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,kc)  &
            - bcl5yl(ic,kc)*ova2yl(ic,kc)
        
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
        fornow = strdyl(ic,kc)*acouyl(ic,kc)*bcl1yl(ic,kc)
        bcl1yl(ic,kc) = half*(strvyl(ic,kc)-acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)-fornow)
        bcl2yl(ic,kc) = strvyl(ic,kc)  &
            *(bcl2yl(ic,kc)-bcl5yl(ic,kc)*ova2yl(ic,kc))
        bcl3yl(ic,kc) = strvyl(ic,kc)*bcl3yl(ic,kc)
        bcl4yl(ic,kc) = strvyl(ic,kc)*bcl4yl(ic,kc)
        bcl5yl(ic,kc) = half*(strvyl(ic,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y UNCHANGED
        fornow = bcl1yl(ic,kc) - strdyl(ic,kc)*acouyl(ic,kc)*dvdtyl(ic,kc)
        bcl2yl(ic,kc) = -dddtyl(ic,kc)  &
            - ova2yl(ic,kc)*(bcl1yl(ic,kc)+fornow) - bcl2yl(ic,kc)
        bcl3yl(ic,kc) = -dudtyl(ic,kc) - bcl3yl(ic,kc)
        bcl4yl(ic,kc) = -dwdtyl(ic,kc) - bcl4yl(ic,kc)
        bcl5yl(ic,kc) = fornow - bcl5yl(ic,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc) - bcl2yl(ic,kc)*streyl(ic,kc)  &
            - bcl3yl(ic,kc)*strdyl(ic,kc)*struyl(ic,kc)  &
            - bcl4yl(ic,kc)*strdyl(ic,kc)*strwyl(ic,kc)  &
            - bcl5yl(ic,kc)*(ova2yl(ic,kc)*streyl(ic,kc)  &
            + strvyl(ic,kc)/acouyl(ic,kc) + ovgmyl(ic,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          bclyyl(ic,kc,ispec) = rateyl(ic,kc,ispec)/strdyl(ic,kc)  &
              - dydtyl(ic,kc,ispec) - strvyl(ic,kc)*bclyyl(ic,kc,ispec)
          
          erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
              - bclyyl(ic,kc,ispec)*strdyl(ic,kc)*strhyl(ic,kc,ispec)
          
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
        fornow = strdyl(ic,kc)*acouyl(ic,kc)*bcl1yl(ic,kc)
        bcl1yl(ic,kc) = half*(strvyl(ic,kc)-acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)-fornow)
        bcl3yl(ic,kc) = strvyl(ic,kc)*bcl3yl(ic,kc)
        bcl4yl(ic,kc) = strvyl(ic,kc)*bcl4yl(ic,kc)
        bcl5yl(ic,kc) = half*(strvyl(ic,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y,L2Y UNCHANGED
        bcl3yl(ic,kc) = -dudtyl(ic,kc) - bcl3yl(ic,kc)
        bcl4yl(ic,kc) = -dwdtyl(ic,kc) - bcl4yl(ic,kc)
        bcl5yl(ic,kc) = bcl1yl(ic,kc)  &
            - strdyl(ic,kc)*acouyl(ic,kc)*dvdtyl(ic,kc) - bcl5yl(ic,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl5yl(ic,kc)*ova2yl(ic,kc)
        
        erhs(ic,jstal,kc) = erhs(ic,jstal,kc)  &
            - bcl3yl(ic,kc)*strdyl(ic,kc)*struyl(ic,kc)  &
            - bcl4yl(ic,kc)*strdyl(ic,kc)*strwyl(ic,kc)  &
            - bcl5yl(ic,kc)*(ova2yl(ic,kc)*streyl(ic,kc)  &
            + strvyl(ic,kc)/acouyl(ic,kc) + ovgmyl(ic,kc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)  &
              - bcl5yl(ic,kc)*ova2yl(ic,kc)*stryyl(ic,kc,ispec)
          
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
        fornow = strdyl(ic,kc)*acouyl(ic,kc)*bcl1yl(ic,kc)
        bcl1yl(ic,kc) = half*(strvyl(ic,kc)-acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)-fornow)
        bcl2yl(ic,kc) = strvyl(ic,kc)  &
            *(bcl2yl(ic,kc)-bcl5yl(ic,kc)*ova2yl(ic,kc))
        bcl3yl(ic,kc) = strvyl(ic,kc)*bcl3yl(ic,kc)
        bcl4yl(ic,kc) = strvyl(ic,kc)*bcl4yl(ic,kc)
        bcl5yl(ic,kc) = half*(strvyl(ic,kc)+acouyl(ic,kc))  &
            *(bcl5yl(ic,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y UNCHANGED
        bcl3yl(ic,kc) = -dudtyl(ic,kc) - bcl3yl(ic,kc)
        bcl4yl(ic,kc) = -dwdtyl(ic,kc) - bcl4yl(ic,kc)
        bcl5yl(ic,kc) = bcl1yl(ic,kc)  &
            - strdyl(ic,kc)*acouyl(ic,kc)*dvdtyl(ic,kc) - bcl5yl(ic,kc)
        bcl2yl(ic,kc) = gam1yl(ic,kc)*ova2yl(ic,kc)  &
            *(bcl1yl(ic,kc)+bcl5yl(ic,kc))  &
            + strdyl(ic,kc)*(dtdtyl(ic,kc)/strtyl(ic,kc)  &
            - sorpyl(ic,kc)/strpyl(ic,kc)) - bcl2yl(ic,kc)
        
      END DO
    END DO
    
!         LYY
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYY
          bclyyl(ic,kc,ispec) = strvyl(ic,kc)*bclyyl(ic,kc,ispec)
          
!               UPDATE L2Y
          bcl2yl(ic,kc) = bcl2yl(ic,kc) + (rateyl(ic,kc,ispec)  &
              - strdyl(ic,kc)*bclyyl(ic,kc,ispec)) *rgspec(ispec)/strryl(ic,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstal,kc) = drhs(ic,jstal,kc) - bcl2yl(ic,kc)  &
            - bcl5yl(ic,kc)*ova2yl(ic,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstal,kc,ispec) = yrhs(ic,jstal,kc,ispec)  &
              - (bcl2yl(ic,kc)+bcl5yl(ic,kc)*ova2yl(ic,kc))*stryyl(ic,kc,ispec)
          
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
            - strgyr(ic,kc)*strtyr(ic,kc)*rgspec(ispec)/strryr(ic,kc)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      gam1yr(ic,kc) = strgyr(ic,kc) - strryr(ic,kc)
      streyr(ic,kc) = streyr(ic,kc) - gam1yr(ic,kc)*strtyr(ic,kc)
      
      gam1yr(ic,kc) = strryr(ic,kc)/gam1yr(ic,kc)
      ovgmyr(ic,kc) = one/gam1yr(ic,kc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO kc = kstal,kstol
    DO ic = istal,istol
      
      fornow = strgyr(ic,kc)*gam1yr(ic,kc)*strtyr(ic,kc)
      acouyr(ic,kc) = SQRT(fornow)
      ova2yr(ic,kc) = one/fornow
      
!     THESE ARE INTRODUCED FOR LODATO'S BC- NC
      
      tt1yr(ic,kc)=t51byr(ic,kc)+ t52byr(ic,kc)*(gam1yr(ic,kc)+1.0)-  &
          strdyr(ic,kc)* acouyr(ic,kc)*t2byr(ic,kc)
      
      tt2yr(ic,kc)=acouyr(ic,kc)*acouyr(ic,kc)*  &
          t1byr(ic,kc)-t51byr(ic,kc)-(gam1yr(ic,kc)+1.0)* t52byr(ic,kc)
      
      tt3yr(ic,kc)=t3byr(ic,kc)
      tt4yr(ic,kc)=t4byr(ic,kc)
      
      tt5yr(ic,kc)=t51byr(ic,kc)+ t52byr(ic,kc)*(gam1yr(ic,kc)+1.0)+  &
          strdyr(ic,kc)* acouyr(ic,kc)*t2byr(ic,kc)
      
      DO is=1,nspec
        tt6yr(ic,kc,is)=t6byr(ic,kc,is)
      END DO
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbcyr == nsbco1)THEN
    flag_pio_yr=nyrprm(1)!0
    flag_bet_yr=nyrprm(2)!1
    
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
        
        bcl2yr(ic,kc) = strvyr(ic,kc)  &
            *(bcl2yr(ic,kc)-bcl5yr(ic,kc)*ova2yr(ic,kc))
        bcl3yr(ic,kc) = strvyr(ic,kc)*bcl3yr(ic,kc)
        bcl4yr(ic,kc) = strvyr(ic,kc)*bcl4yr(ic,kc)
!             OLD VALUE OF L1Y
        bcl1yr(ic,kc) = half*(strvyr(ic,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)-strdyr(ic,kc)*acouyr(ic,kc)*bcl1yr(ic,kc))
        
!             SUBTRACT FROM NEW VALUE OF L1Y
!     TERMS INTRODUCED FOR LODATO'S BC- NC
        IF(flag_bet_yr==1) THEN
          bet=struyr(ic,kc)*struyr(ic,kc)+ strvyr(ic,kc)*strvyr(ic,kc)+  &
              strwyr(ic,kc)*strwyr(ic,kc)
          bet=SQRT(bet)/acouyr(ic,kc)
        END IF
        
        IF((strvyr(ic,kc) < zero).AND.(flag_pio_yr==1))THEN
          bcl2yr(ic,kc) = -bcl2yr(ic,kc)  &
              -ova2yr(ic,kc)*sorpyr(ic,kc)
          
          bcl3yr(ic,kc) = 0.1*(struyr(ic,kc)-0.0D0)  &
              -bcl3yr(ic,kc)
          
          bcl4yr(ic,kc) = 0.1*(strwyr(ic,kc)-0.0D0)  &
              -bcl4yr(ic,kc)
          bcl1yr(ic,kc)= half*sorpyr(ic,kc)  &
              + cobcyr*acouyr(ic,kc)*(strpyr(ic,kc)-pinfyr)  &
              +0.5*(1.0-bet)*tt1yr(ic,kc)- bcl1yr(ic,kc) &
!              17 JAN 2024 VM: change making negative
!     +                     + 100.0*STRDYR(IC,KC)*(ONE/YGDLEN)!ASK NC: CHANGED TO YGDLEN  &
          - 100.0*strdyr(ic,kc)*(one/ygdlen)  &!ASK NC: CHANGED TO YGDLEN
              *(acouyr(ic,kc)**2.0-strvyr(ic,kc)**2.0) *(strvyr(ic,kc)-0.0)
        ELSE
          bcl1yr(ic,kc)= half*sorpyr(ic,kc)  &
              + cobcyr*acouyr(ic,kc)*(strpyr(ic,kc)-pinfyr)  &
              +0.5*(1.0-bet)*tt1yr(ic,kc)- bcl1yr(ic,kc)
        END IF
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        IF((strvyr(ic,kc) < zero).AND.(flag_pio_yr==1))THEN
          
          drhs(ic,jstol,kc) = drhs(ic,jstol,kc)  &
              - bcl1yr(ic,kc)*ova2yr(ic,kc) - bcl2yr(ic,kc)
          
          urhs(ic,jstol,kc) = urhs(ic,jstol,kc)  &
              - bcl1yr(ic,kc)*ova2yr(ic,kc)*struyr(ic,kc)  &
              - bcl2yr(ic,kc)*struyr(ic,kc) - bcl3yr(ic,kc)*strdyr(ic,kc)
          
          vrhs(ic,jstol,kc) = vrhs(ic,jstol,kc)  &
              - bcl1yr(ic,kc)*ova2yr(ic,kc)*(strvyr(ic,kc)-acouyr(ic,kc))  &
              - bcl2yr(ic,kc)*strvyr(ic,kc)
          
          wrhs(ic,jstol,kc) = wrhs(ic,jstol,kc)  &
              - bcl1yr(ic,kc)*ova2yr(ic,kc)*strwyr(ic,kc)  &
              - bcl2yr(ic,kc)*strwyr(ic,kc) - bcl4yr(ic,kc)*strdyr(ic,kc)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
              - bcl1yr(ic,kc)*(ova2yr(ic,kc)*streyr(ic,kc) &
!              16 JAN 2024: VM+NC sign change
!     +                                     + STRVYR(IC,KC)/ACOUYR(IC,KC)  &
          - strvyr(ic,kc)/acouyr(ic,kc) + ovgmyr(ic,kc))  &
              - bcl2yr(ic,kc)*streyr(ic,kc)  &
              - bcl3yr(ic,kc)*strdyr(ic,kc)*struyr(ic,kc)  &
              - bcl4yr(ic,kc)*strdyr(ic,kc)*strwyr(ic,kc)
          
        ELSE
          drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,kc)*ova2yr(ic,kc)
          
          urhs(ic,jstol,kc) = urhs(ic,jstol,kc)  &
              - bcl1yr(ic,kc)*ova2yr(ic,kc)*struyr(ic,kc)
          
          vrhs(ic,jstol,kc) = vrhs(ic,jstol,kc)  &
              - bcl1yr(ic,kc)*ova2yr(ic,kc)*(strvyr(ic,kc)-acouyr(ic,kc))
          
          wrhs(ic,jstol,kc) = wrhs(ic,jstol,kc)  &
              - bcl1yr(ic,kc)*ova2yr(ic,kc)*strwyr(ic,kc)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
              - bcl1yr(ic,kc)*(ova2yr(ic,kc)*streyr(ic,kc)  &
              - strvyr(ic,kc)/acouyr(ic,kc) + ovgmyr(ic,kc))
        END IF
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          IF((strvyr(ic,kc) < zero).AND.(flag_pio_yr==1))THEN
            
            fornow = bclyyr(ic,kc,ispec)*strdyr(ic,kc)
            
            erhs(ic,jstol,kc) = erhs(ic,jstol,kc) - fornow*strhyr(ic,kc,ispec)
            
            yrhs(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)  &
                - (bcl2yr(ic,kc)+bcl1yr(ic,kc)*ova2yr(ic,kc))*stryyr(ic,kc,ispec)  &
                - fornow
            
          ELSE
            yrhs(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)  &
                - bcl1yr(ic,kc)*ova2yr(ic,kc)*stryyr(ic,kc,ispec)
          END IF
          
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
        fornow = strdyr(ic,kc)*acouyr(ic,kc)*bcl1yr(ic,kc)
        bcl1yr(ic,kc) = half*(strvyr(ic,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)-fornow)
        bcl2yr(ic,kc) = strvyr(ic,kc)  &
            *(bcl2yr(ic,kc)-bcl5yr(ic,kc)*ova2yr(ic,kc))
        bcl3yr(ic,kc) = strvyr(ic,kc)*bcl3yr(ic,kc)
        bcl4yr(ic,kc) = strvyr(ic,kc)*bcl4yr(ic,kc)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
!             L5Y UNCHANGED
        bcl1yr(ic,kc) = half*sorpyr(ic,kc)  &
            + cobcyr*acouyr(ic,kc)*(strpyr(ic,kc)-pinfyr) - bcl1yr(ic,kc)
        bcl2yr(ic,kc) = -bcl2yr(ic,kc)
        bcl3yr(ic,kc) = -bcl3yr(ic,kc)
        bcl4yr(ic,kc) = -bcl4yr(ic,kc)
        
      END DO
    END DO
    
!         LYY
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
!               OLD VALUE OF L's
          bclyyr(ic,kc,ispec) = strvyr(ic,kc)*bclyyr(ic,kc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
          bclyyr(ic,kc,ispec) = rateyr(ic,kc,ispec)/strdyr(ic,kc)  &
              - bclyyr(ic,kc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,kc)*ova2yr(ic,kc)  &
            - bcl2yr(ic,kc)
        
        urhs(ic,jstol,kc) = urhs(ic,jstol,kc)  &
            - bcl1yr(ic,kc)*ova2yr(ic,kc)*struyr(ic,kc)  &
            - bcl2yr(ic,kc)*struyr(ic,kc) - bcl3yr(ic,kc)*strdyr(ic,kc)
        
        vrhs(ic,jstol,kc) = vrhs(ic,jstol,kc)  &
            - bcl1yr(ic,kc)*ova2yr(ic,kc)*(strvyr(ic,kc)-acouyr(ic,kc))  &
            - bcl2yr(ic,kc)*strvyr(ic,kc)
        
        wrhs(ic,jstol,kc) = wrhs(ic,jstol,kc)  &
            - bcl1yr(ic,kc)*ova2yr(ic,kc)*strwyr(ic,kc)  &
            - bcl2yr(ic,kc)*strwyr(ic,kc) - bcl4yr(ic,kc)*strdyr(ic,kc)
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,kc)*(ova2yr(ic,kc)*streyr(ic,kc) &
!                     16 JAN 2024: VM+NC sign change
!     +                                     + STRVYR(IC,KC)/ACOUYR(IC,KC)  &
        - strvyr(ic,kc)/acouyr(ic,kc) + ovgmyr(ic,kc))  &
            - bcl2yr(ic,kc)*streyr(ic,kc)  &
            - bcl3yr(ic,kc)*strdyr(ic,kc)*struyr(ic,kc)  &
            - bcl4yr(ic,kc)*strdyr(ic,kc)*strwyr(ic,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          fornow = bclyyr(ic,kc,ispec)*strdyr(ic,kc)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc) - fornow*strhyr(ic,kc,ispec)
          
          yrhs(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)  &
              - (bcl2yr(ic,kc)+bcl1yr(ic,kc)*ova2yr(ic,kc))*stryyr(ic,kc,ispec)  &
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
        
        sydtyr(ic,kc) = sydtyr(ic,kc)/strryr(ic,kc)
        sorpyr(ic,kc) = -sorpyr(ic,kc)*gam1yr(ic,kc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Y,L2Y,L5Y
    DO kc = kstal,kstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdyr(ic,kc)*acouyr(ic,kc)*bcl1yr(ic,kc)
        bcl1yr(ic,kc) = half*(strvyr(ic,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)-fornow)
        bcl2yr(ic,kc) = strvyr(ic,kc)  &
            *(bcl2yr(ic,kc)-bcl5yr(ic,kc)*ova2yr(ic,kc))
        bcl5yr(ic,kc) = half*(strvyr(ic,kc)+acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Y UNCHANGED
        bcl1yr(ic,kc) = bcl5yr(ic,kc)  &
            + strdyr(ic,kc)*acouyr(ic,kc)*dvdtyr(ic,kc) - bcl1yr(ic,kc)
        bcl2yr(ic,kc) = gam1yr(ic,kc)*ova2yr(ic,kc)  &
            *(bcl1yr(ic,kc)+bcl5yr(ic,kc))  &
            + strdyr(ic,kc)*(dtdtyr(ic,kc)/strtyr(ic,kc)  &
            - sorpyr(ic,kc)/strpyr(ic,kc) + sydtyr(ic,kc))  &
            - bcl2yr(ic,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,kc)*ova2yr(ic,kc)  &
            - bcl2yr(ic,kc)
        
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
        fornow = strdyr(ic,kc)*acouyr(ic,kc)*bcl1yr(ic,kc)
        bcl1yr(ic,kc) = half*(strvyr(ic,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)-fornow)
        bcl2yr(ic,kc) = strvyr(ic,kc)  &
            *(bcl2yr(ic,kc)-bcl5yr(ic,kc)*ova2yr(ic,kc))
        bcl3yr(ic,kc) = strvyr(ic,kc)*bcl3yr(ic,kc)
        bcl4yr(ic,kc) = strvyr(ic,kc)*bcl4yr(ic,kc)
        bcl5yr(ic,kc) = half*(strvyr(ic,kc)+acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Y UNCHANGED
        fornow = bcl5yr(ic,kc) + strdyr(ic,kc)*acouyr(ic,kc)*dvdtyr(ic,kc)
        bcl1yr(ic,kc) = fornow - bcl1yr(ic,kc)
        bcl2yr(ic,kc) = -dddtyr(ic,kc)  &
            - ova2yr(ic,kc)*(bcl1yr(ic,kc)+fornow) - bcl2yr(ic,kc)
        bcl3yr(ic,kc) = -dudtyr(ic,kc) - bcl3yr(ic,kc)
        bcl4yr(ic,kc) = -dwdtyr(ic,kc) - bcl4yr(ic,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,kc)*(ova2yr(ic,kc)*streyr(ic,kc)  &
            + strvyr(ic,kc)/acouyr(ic,kc) + ovgmyr(ic,kc))  &
            - bcl2yr(ic,kc)*streyr(ic,kc)  &
            - bcl3yr(ic,kc)*strdyr(ic,kc)*struyr(ic,kc)  &
            - bcl4yr(ic,kc)*strdyr(ic,kc)*strwyr(ic,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          bclyyr(ic,kc,ispec) = rateyr(ic,kc,ispec)/strdyr(ic,kc)  &
              - dydtyr(ic,kc,ispec) - strvyr(ic,kc)*bclyyr(ic,kc,ispec)
          
          erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
              - bclyyr(ic,kc,ispec)*strdyr(ic,kc)*strhyr(ic,kc,ispec)
          
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
        fornow = strdyr(ic,kc)*acouyr(ic,kc)*bcl1yr(ic,kc)
        bcl1yr(ic,kc) = half*(strvyr(ic,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)-fornow)
        bcl3yr(ic,kc) = strvyr(ic,kc)*bcl3yr(ic,kc)
        bcl4yr(ic,kc) = strvyr(ic,kc)*bcl4yr(ic,kc)
        bcl5yr(ic,kc) = half*(strvyr(ic,kc)+acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L2Y,L5Y UNCHANGED
        bcl1yr(ic,kc) = bcl5yr(ic,kc)  &
            + strdyr(ic,kc)*acouyr(ic,kc)*dvdtyr(ic,kc) - bcl1yr(ic,kc)
        bcl3yr(ic,kc) = -dudtyr(ic,kc) - bcl3yr(ic,kc)
        bcl4yr(ic,kc) = -dwdtyr(ic,kc) - bcl4yr(ic,kc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,kc)*ova2yr(ic,kc)
        
        erhs(ic,jstol,kc) = erhs(ic,jstol,kc)  &
            - bcl1yr(ic,kc)*(ova2yr(ic,kc)*streyr(ic,kc)  &
            + strvyr(ic,kc)/acouyr(ic,kc) + ovgmyr(ic,kc))  &
            - bcl3yr(ic,kc)*strdyr(ic,kc)*struyr(ic,kc)  &
            - bcl4yr(ic,kc)*strdyr(ic,kc)*strwyr(ic,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)  &
              - bcl1yr(ic,kc)*ova2yr(ic,kc)*stryyr(ic,kc,ispec)
          
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
        fornow = strdyr(ic,kc)*acouyr(ic,kc)*bcl1yr(ic,kc)
        bcl1yr(ic,kc) = half*(strvyr(ic,kc)-acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)-fornow)
        bcl2yr(ic,kc) = strvyr(ic,kc)  &
            *(bcl2yr(ic,kc)-bcl5yr(ic,kc)*ova2yr(ic,kc))
        bcl3yr(ic,kc) = strvyr(ic,kc)*bcl3yr(ic,kc)
        bcl4yr(ic,kc) = strvyr(ic,kc)*bcl4yr(ic,kc)
        bcl5yr(ic,kc) = half*(strvyr(ic,kc)+acouyr(ic,kc))  &
            *(bcl5yr(ic,kc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Y UNCHANGED
        bcl1yr(ic,kc) = bcl5yr(ic,kc)  &
            + strdyr(ic,kc)*acouyr(ic,kc)*dvdtyr(ic,kc) - bcl1yr(ic,kc)
        bcl3yr(ic,kc) = -dudtyr(ic,kc) - bcl3yr(ic,kc)
        bcl4yr(ic,kc) = -dwdtyr(ic,kc) - bcl4yr(ic,kc)
        bcl2yr(ic,kc) = gam1yr(ic,kc)*ova2yr(ic,kc)  &
            *(bcl1yr(ic,kc)+bcl5yr(ic,kc))  &
            + strdyr(ic,kc)*(dtdtyr(ic,kc)/strtyr(ic,kc)  &
            - sorpyr(ic,kc)/strpyr(ic,kc)) - bcl2yr(ic,kc)
        
      END DO
    END DO
    
!         LYY
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYY
          bclyyr(ic,kc,ispec) = strvyr(ic,kc)*bclyyr(ic,kc,ispec)
          
!               UPDATE L2Y
          bcl2yr(ic,kc) = bcl2yr(ic,kc) + (rateyr(ic,kc,ispec)  &
              - strdyr(ic,kc)*bclyyr(ic,kc,ispec)) *rgspec(ispec)/strryr(ic,kc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO kc = kstal,kstol
      DO ic = istal,istol
        
        drhs(ic,jstol,kc) = drhs(ic,jstol,kc) - bcl1yr(ic,kc)*ova2yr(ic,kc)  &
            - bcl2yr(ic,kc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO kc = kstal,kstol
        DO ic = istal,istol
          
          yrhs(ic,jstol,kc,ispec) = yrhs(ic,jstol,kc,ispec)  &
              - (bcl2yr(ic,kc)+bcl1yr(ic,kc)*ova2yr(ic,kc))*stryyr(ic,kc,ispec)
          
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
            - strgzl(ic,jc)*strtzl(ic,jc)*rgspec(ispec)/strrzl(ic,jc)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      gam1zl(ic,jc) = strgzl(ic,jc) - strrzl(ic,jc)
      strezl(ic,jc) = strezl(ic,jc) - gam1zl(ic,jc)*strtzl(ic,jc)
      
      gam1zl(ic,jc) = strrzl(ic,jc)/gam1zl(ic,jc)
      ovgmzl(ic,jc) = one/gam1zl(ic,jc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = strgzl(ic,jc)*gam1zl(ic,jc)*strtzl(ic,jc)
      acouzl(ic,jc) = SQRT(fornow)
      ova2zl(ic,jc) = one/fornow
      
!     THESE ARE INTRODUCED FOR LODATO'S BC- NC
      
      tt1zl(ic,jc)=t51bzl(ic,jc)+ t52bzl(ic,jc)*(gam1zl(ic,jc)+1.0)-  &
          strdzl(ic,jc)* acouzl(ic,jc)*t2bzl(ic,jc)
      
      tt2zl(ic,jc)=acouzl(ic,jc)*acouzl(ic,jc)*  &
          t1bzl(ic,jc)-t51bzl(ic,jc)-(gam1zl(ic,jc)+1.0)* t52bzl(ic,jc)
      
      tt3zl(ic,jc)=t3bzl(ic,jc)
      tt4zl(ic,jc)=t4bzl(ic,jc)
      
      tt5zl(ic,jc)=t51bzl(ic,jc)+ t52bzl(ic,jc)*(gam1zl(ic,jc)+1.0)+  &
          strdzl(ic,jc)* acouzl(ic,jc)*t2bzl(ic,jc)
      
      DO is=1,nspec
        tt6zl(ic,jc,is)=t6bzl(ic,jc,is)
      END DO
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbczl == nsbco1)THEN
    flag_pio_zl=nzlprm(1)!0
    flag_bet_zl=nzlprm(2)!1
    
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
        bcl2zl(ic,jc) = strwzl(ic,jc)  &
            *(bcl2zl(ic,jc)-bcl5zl(ic,jc)*ova2zl(ic,jc))
        bcl3zl(ic,jc) = strwzl(ic,jc)*bcl3zl(ic,jc)
        bcl4zl(ic,jc) = strwzl(ic,jc)*bcl4zl(ic,jc)
        
!             OLD VALUE OF L5Z
        bcl5zl(ic,jc) = half*(strwzl(ic,jc)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)+strdzl(ic,jc)*acouzl(ic,jc)*bcl1zl(ic,jc))
        
!             SUBTRACT FROM NEW VALUE OF L5Z
!             TERMS INTRODUCED FOR LODATO'S BC- NC
        IF(flag_bet_zl==1) THEN
          bet=struzl(ic,jc)*struzl(ic,jc)+ strvzl(ic,jc)*strvzl(ic,jc)+  &
              strwzl(ic,jc)*strwzl(ic,jc)
          bet=SQRT(bet)/acouzl(ic,jc)
        END IF
        IF((strwzl(ic,jc) > zero).AND.(flag_pio_zl==1))THEN
          bcl2zl(ic,jc) = -bcl2zl(ic,jc)  &
              -ova2zl(ic,jc)*sorpzl(ic,jc)
          
          bcl3zl(ic,jc) = 0.1*(struzl(ic,jc)-0.0D0)  &
              -bcl3zl(ic,jc)
          
          bcl4zl(ic,jc) = 0.1*(strvzl(ic,jc)-0.0D0)  &
              -bcl4zl(ic,jc)
          bcl5zl(ic,jc)= half*sorpzl(ic,jc)  &
              + cobczl*acouzl(ic,jc)*(strpzl(ic,jc)-pinfzl)  &
              +0.5*(1.0-bet)*tt5zl(ic,jc)- bcl5zl(ic,jc)  &
              + 100.0*strdzl(ic,jc)*(one/zgdlen)  &
              *(acouzl(ic,jc)**2.0-strwzl(ic,jc)**2.0) *(strwzl(ic,jc)-0.0)
        ELSE
          bcl5zl(ic,jc)= half*sorpzl(ic,jc)  &
              + cobczl*acouzl(ic,jc)*(strpzl(ic,jc)-pinfzl)  &
              +0.5*(1.0-bet)*tt5zl(ic,jc)- bcl5zl(ic,jc)
        END IF
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        IF((strwzl(ic,jc) > zero).AND.(flag_pio_zl==1))THEN
          drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc)  &
              - bcl5zl(ic,jc)*ova2zl(ic,jc)
          
          urhs(ic,jc,kstal) = urhs(ic,jc,kstal)  &
              - bcl2zl(ic,jc)*struzl(ic,jc) - bcl3zl(ic,jc)*strdzl(ic,jc)  &
              - bcl5zl(ic,jc)*ova2zl(ic,jc)*struzl(ic,jc)
          
          vrhs(ic,jc,kstal) = vrhs(ic,jc,kstal)  &
              - bcl2zl(ic,jc)*strvzl(ic,jc) - bcl4zl(ic,jc)*strdzl(ic,jc)  &
              - bcl5zl(ic,jc)*ova2zl(ic,jc)*strvzl(ic,jc)
          
          wrhs(ic,jc,kstal) = wrhs(ic,jc,kstal)  &
              - bcl2zl(ic,jc)*strwzl(ic,jc)  &
              - bcl5zl(ic,jc)*ova2zl(ic,jc)*(strwzl(ic,jc)+acouzl(ic,jc))
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
              - bcl2zl(ic,jc)*strezl(ic,jc)  &
              - bcl3zl(ic,jc)*strdzl(ic,jc)*struzl(ic,jc)  &
              - bcl4zl(ic,jc)*strdzl(ic,jc)*strvzl(ic,jc)  &
              - bcl5zl(ic,jc)*(ova2zl(ic,jc)*strezl(ic,jc)  &
              + strwzl(ic,jc)/acouzl(ic,jc) + ovgmzl(ic,jc))
        ELSE
          drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl5zl(ic,jc)*ova2zl(ic,jc)
          
          urhs(ic,jc,kstal) = urhs(ic,jc,kstal)  &
              - bcl5zl(ic,jc)*ova2zl(ic,jc)*struzl(ic,jc)
          
          vrhs(ic,jc,kstal) = vrhs(ic,jc,kstal)  &
              - bcl5zl(ic,jc)*ova2zl(ic,jc)*strvzl(ic,jc)
          
          wrhs(ic,jc,kstal) = wrhs(ic,jc,kstal)  &
              - bcl5zl(ic,jc)*ova2zl(ic,jc)*(strwzl(ic,jc)+acouzl(ic,jc))
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
              - bcl5zl(ic,jc)*(ova2zl(ic,jc)*strezl(ic,jc)  &
              + strwzl(ic,jc)/acouzl(ic,jc) + ovgmzl(ic,jc))
        END IF
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          IF((strwzl(ic,jc) > zero).AND.(flag_pio_zl==1))THEN
            fornow = bclyzl(ic,jc,ispec)*strdzl(ic,jc)
            
            erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - fornow*strhzl(ic,jc,ispec)
            
            yrhs(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)  &
                - (bcl2zl(ic,jc)+bcl5zl(ic,jc)*ova2zl(ic,jc))*stryzl(ic,jc,ispec)  &
                - fornow
          ELSE
            yrhs(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)  &
                - bcl5zl(ic,jc)*ova2zl(ic,jc)*stryzl(ic,jc,ispec)
          END IF
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
        fornow = strdzl(ic,jc)*acouzl(ic,jc)*bcl1zl(ic,jc)
        bcl2zl(ic,jc) = strwzl(ic,jc)  &
            *(bcl2zl(ic,jc)-bcl5zl(ic,jc)*ova2zl(ic,jc))
        bcl3zl(ic,jc) = strwzl(ic,jc)*bcl3zl(ic,jc)
        bcl4zl(ic,jc) = strwzl(ic,jc)*bcl4zl(ic,jc)
        bcl5zl(ic,jc) = half*(strwzl(ic,jc)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
!             L1Z UNCHANGED
        bcl2zl(ic,jc) = -bcl2zl(ic,jc)
        bcl3zl(ic,jc) = -bcl3zl(ic,jc)
        bcl4zl(ic,jc) = -bcl4zl(ic,jc)
        bcl5zl(ic,jc) = half*sorpzl(ic,jc)  &
            + cobczl*acouzl(ic,jc)*(strpzl(ic,jc)-pinfzl) - bcl5zl(ic,jc)
        
      END DO
    END DO
    
!         LYZ
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
!               OLD VALUE OF L's
          bclyzl(ic,jc,ispec) = strwzl(ic,jc)*bclyzl(ic,jc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
          bclyzl(ic,jc,ispec) = ratezl(ic,jc,ispec)/strdzl(ic,jc)  &
              - bclyzl(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc)  &
            - bcl5zl(ic,jc)*ova2zl(ic,jc)
        
        urhs(ic,jc,kstal) = urhs(ic,jc,kstal) - bcl2zl(ic,jc)*struzl(ic,jc)  &
            - bcl3zl(ic,jc)*strdzl(ic,jc)  &
            - bcl5zl(ic,jc)*ova2zl(ic,jc)*struzl(ic,jc)
        
        vrhs(ic,jc,kstal) = vrhs(ic,jc,kstal) - bcl2zl(ic,jc)*strvzl(ic,jc)  &
            - bcl4zl(ic,jc)*strdzl(ic,jc)  &
            - bcl5zl(ic,jc)*ova2zl(ic,jc)*strvzl(ic,jc)
        
        wrhs(ic,jc,kstal) = wrhs(ic,jc,kstal) - bcl2zl(ic,jc)*strwzl(ic,jc)  &
            - bcl5zl(ic,jc)*ova2zl(ic,jc)*(strwzl(ic,jc)+acouzl(ic,jc))
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - bcl2zl(ic,jc)*strezl(ic,jc)  &
            - bcl3zl(ic,jc)*strdzl(ic,jc)*struzl(ic,jc)  &
            - bcl4zl(ic,jc)*strdzl(ic,jc)*strvzl(ic,jc)  &
            - bcl5zl(ic,jc)*(ova2zl(ic,jc)*strezl(ic,jc)  &
            + strwzl(ic,jc)/acouzl(ic,jc) + ovgmzl(ic,jc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = bclyzl(ic,jc,ispec)*strdzl(ic,jc)
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - fornow*strhzl(ic,jc,ispec)
          
          yrhs(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)  &
              - (bcl2zl(ic,jc)+bcl5zl(ic,jc)*ova2zl(ic,jc))*stryzl(ic,jc,ispec)  &
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
        
        sydtzl(ic,jc) = sydtzl(ic,jc)/strrzl(ic,jc)
        sorpzl(ic,jc) = -sorpzl(ic,jc)*gam1zl(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Z,L2Z,L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzl(ic,jc)*acouzl(ic,jc)*bcl1zl(ic,jc)
        bcl1zl(ic,jc) = half*(strwzl(ic,jc)-acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)-fornow)
        bcl2zl(ic,jc) = strwzl(ic,jc)  &
            *(bcl2zl(ic,jc)-bcl5zl(ic,jc)*ova2zl(ic,jc))
        bcl5zl(ic,jc) = half*(strwzl(ic,jc)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Z UNCHANGED
        bcl5zl(ic,jc) = bcl1zl(ic,jc)  &
            - strdzl(ic,jc)*acouzl(ic,jc)*dwdtzl(ic,jc) - bcl5zl(ic,jc)
        bcl2zl(ic,jc) = gam1zl(ic,jc)*ova2zl(ic,jc)  &
            *(bcl1zl(ic,jc)+bcl5zl(ic,jc))  &
            + strdzl(ic,jc)*(dtdtzl(ic,jc)/strtzl(ic,jc)  &
            - sorpzl(ic,jc)/strpzl(ic,jc) + sydtzl(ic,jc))  &
            - bcl2zl(ic,jc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc)  &
            - bcl5zl(ic,jc)*ova2zl(ic,jc)
        
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
        fornow = strdzl(ic,jc)*acouzl(ic,jc)*bcl1zl(ic,jc)
        bcl1zl(ic,jc) = half*(strwzl(ic,jc)-acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)-fornow)
        bcl2zl(ic,jc) = strwzl(ic,jc)  &
            *(bcl2zl(ic,jc)-bcl5zl(ic,jc)*ova2zl(ic,jc))
        bcl3zl(ic,jc) = strwzl(ic,jc)*bcl3zl(ic,jc)
        bcl4zl(ic,jc) = strwzl(ic,jc)*bcl4zl(ic,jc)
        bcl5zl(ic,jc) = half*(strwzl(ic,jc)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Z UNCHANGED
        fornow = bcl1zl(ic,jc) - strdzl(ic,jc)*acouzl(ic,jc)*dwdtzl(ic,jc)
        bcl2zl(ic,jc) = -dddtzl(ic,jc)  &
            - ova2zl(ic,jc)*(bcl1zl(ic,jc)+fornow) - bcl2zl(ic,jc)
        bcl3zl(ic,jc) = -dudtzl(ic,jc) - bcl3zl(ic,jc)
        bcl4zl(ic,jc) = -dvdtzl(ic,jc) - bcl4zl(ic,jc)
        bcl5zl(ic,jc) = fornow - bcl5zl(ic,jc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal) - bcl2zl(ic,jc)*strezl(ic,jc)  &
            - bcl3zl(ic,jc)*strdzl(ic,jc)*struzl(ic,jc)  &
            - bcl4zl(ic,jc)*strdzl(ic,jc)*strvzl(ic,jc)  &
            - bcl5zl(ic,jc)*(ova2zl(ic,jc)*strezl(ic,jc)  &
            + strwzl(ic,jc)/acouzl(ic,jc) + ovgmzl(ic,jc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          bclyzl(ic,jc,ispec) = ratezl(ic,jc,ispec)/strdzl(ic,jc)  &
              - dydtzl(ic,jc,ispec) - strwzl(ic,jc)*bclyzl(ic,jc,ispec)
          
          erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
              - bclyzl(ic,jc,ispec)*strdzl(ic,jc)*strhzl(ic,jc,ispec)
          
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
        fornow = strdzl(ic,jc)*acouzl(ic,jc)*bcl1zl(ic,jc)
        bcl1zl(ic,jc) = half*(strwzl(ic,jc)-acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)-fornow)
        bcl3zl(ic,jc) = strwzl(ic,jc)*bcl3zl(ic,jc)
        bcl4zl(ic,jc) = strwzl(ic,jc)*bcl4zl(ic,jc)
        bcl5zl(ic,jc) = half*(strwzl(ic,jc)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Z,L2Z UNCHANGED
        bcl3zl(ic,jc) = -dudtzl(ic,jc) - bcl3zl(ic,jc)
        bcl4zl(ic,jc) = -dvdtzl(ic,jc) - bcl4zl(ic,jc)
        bcl5zl(ic,jc) = bcl1zl(ic,jc)  &
            - strdzl(ic,jc)*acouzl(ic,jc)*dwdtzl(ic,jc) - bcl5zl(ic,jc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl5zl(ic,jc)*ova2zl(ic,jc)
        
        erhs(ic,jc,kstal) = erhs(ic,jc,kstal)  &
            - bcl3zl(ic,jc)*strdzl(ic,jc)*struzl(ic,jc)  &
            - bcl4zl(ic,jc)*strdzl(ic,jc)*strvzl(ic,jc)  &
            - bcl5zl(ic,jc)*(ova2zl(ic,jc)*strezl(ic,jc)  &
            + strwzl(ic,jc)/acouzl(ic,jc) + ovgmzl(ic,jc))
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)  &
              - bcl5zl(ic,jc)*ova2zl(ic,jc)*stryzl(ic,jc,ispec)
          
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
        fornow = strdzl(ic,jc)*acouzl(ic,jc)*bcl1zl(ic,jc)
        bcl1zl(ic,jc) = half*(strwzl(ic,jc)-acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)-fornow)
        bcl2zl(ic,jc) = strwzl(ic,jc)  &
            *(bcl2zl(ic,jc)-bcl5zl(ic,jc)*ova2zl(ic,jc))
        bcl3zl(ic,jc) = strwzl(ic,jc)*bcl3zl(ic,jc)
        bcl4zl(ic,jc) = strwzl(ic,jc)*bcl4zl(ic,jc)
        bcl5zl(ic,jc) = half*(strwzl(ic,jc)+acouzl(ic,jc))  &
            *(bcl5zl(ic,jc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L1Y UNCHANGED
        bcl3zl(ic,jc) = -dudtzl(ic,jc) - bcl3zl(ic,jc)
        bcl4zl(ic,jc) = -dvdtzl(ic,jc) - bcl4zl(ic,jc)
        bcl5zl(ic,jc) = bcl1zl(ic,jc)  &
            - strdzl(ic,jc)*acouzl(ic,jc)*dwdtzl(ic,jc) - bcl5zl(ic,jc)
        bcl2zl(ic,jc) = gam1zl(ic,jc)*ova2zl(ic,jc)  &
            *(bcl1zl(ic,jc)+bcl5zl(ic,jc))  &
            + strdzl(ic,jc)*(dtdtzl(ic,jc)/strtzl(ic,jc)  &
            - sorpzl(ic,jc)/strpzl(ic,jc)) - bcl2zl(ic,jc)
        
      END DO
    END DO
    
!         LYZ
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYZ
          bclyzl(ic,jc,ispec) = strwzl(ic,jc)*bclyzl(ic,jc,ispec)
          
!               UPDATE L2Z
          bcl2zl(ic,jc) = bcl2zl(ic,jc) + (ratezl(ic,jc,ispec)  &
              - strdzl(ic,jc)*bclyzl(ic,jc,ispec)) *rgspec(ispec)/strrzl(ic,jc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstal) = drhs(ic,jc,kstal) - bcl2zl(ic,jc)  &
            - bcl5zl(ic,jc)*ova2zl(ic,jc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstal,ispec) = yrhs(ic,jc,kstal,ispec)  &
              - (bcl2zl(ic,jc)+bcl5zl(ic,jc)*ova2zl(ic,jc))*stryzl(ic,jc,ispec)
          
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
            - strgzr(ic,jc)*strtzr(ic,jc)*rgspec(ispec)/strrzr(ic,jc)
        
      END DO
    END DO
    
  END DO
  
!       REDUCED INTERNAL ENERGY
!       -----------------------
!       GAMMA-1, 1/(GAMMA-1)
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      gam1zr(ic,jc) = strgzr(ic,jc) - strrzr(ic,jc)
      strezr(ic,jc) = strezr(ic,jc) - gam1zr(ic,jc)*strtzr(ic,jc)
      
      gam1zr(ic,jc) = strrzr(ic,jc)/gam1zr(ic,jc)
      ovgmzr(ic,jc) = one/gam1zr(ic,jc)
      
    END DO
  END DO
  
!       SPEED OF SOUND
!       --------------
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fornow = strgzr(ic,jc)*gam1zr(ic,jc)*strtzr(ic,jc)
      acouzr(ic,jc) = SQRT(fornow)
      ova2zr(ic,jc) = one/fornow
      
!     THESE ARE INTRODUCED FOR LODATO'S BC- NC
      
      tt1zr(ic,jc)=t51bzr(ic,jc)+ t52bzr(ic,jc)*(gam1zr(ic,jc)+1.0)-  &
          strdzr(ic,jc)* acouzr(ic,jc)*t2bzr(ic,jc)
      
      tt2zr(ic,jc)=acouzr(ic,jc)*acouzr(ic,jc)*  &
          t1bzr(ic,jc)-t51bzr(ic,jc)-(gam1zr(ic,jc)+1.0)* t52bzr(ic,jc)
      
      tt3zr(ic,jc)=t3bzr(ic,jc)
      tt4zr(ic,jc)=t4bzr(ic,jc)
      
      tt5zr(ic,jc)=t51bzr(ic,jc)+ t52bzr(ic,jc)*(gam1zr(ic,jc)+1.0)+  &
          strdzr(ic,jc)* acouzr(ic,jc)*t2bzr(ic,jc)
      
      DO is=1,nspec
        tt6zr(ic,jc,is)=t6bzr(ic,jc,is)
      END DO
      
    END DO
  END DO
  
!       =======================================================================
  
!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------
  
  IF(nsbczr == nsbco1)THEN
    flag_pio_zr=nzrprm(1)!0
    flag_bet_zr=nzrprm(2)!1
    
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
        
        bcl2zr(ic,jc) = strwzr(ic,jc)  &
            *(bcl2zr(ic,jc)-bcl5zr(ic,jc)*ova2zr(ic,jc))
        bcl3zr(ic,jc) = strwzr(ic,jc)*bcl3zr(ic,jc)
        bcl4zr(ic,jc) = strwzr(ic,jc)*bcl4zr(ic,jc)
!             OLD VALUE OF L1Z
        bcl1zr(ic,jc) = half*(strwzr(ic,jc)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc)-strdzr(ic,jc)*acouzr(ic,jc)*bcl1zr(ic,jc))
        
!             SUBTRACT FROM NEW VALUE OF L1Z
!             TERMS INTRODUCED FOR LODATO'S BC- NC
        IF(flag_bet_zr==1) THEN
          bet=struzr(ic,jc)*struzr(ic,jc)+ strvzr(ic,jc)*strvzr(ic,jc)+  &
              strwzr(ic,jc)*strwzr(ic,jc)
          bet=SQRT(bet)/acouzr(ic,jc)
        END IF
        IF((strwzr(ic,jc) < zero).AND.(flag_pio_zr==1))THEN
          bcl2zr(ic,jc) = -bcl2zr(ic,jc)  &
              -ova2zr(ic,jc)*sorpzr(ic,jc)
          
          bcl3zr(ic,jc) = 0.1*(struzr(ic,jc)-0.0D0)  &
              -bcl3zr(ic,jc)
          
          bcl4zr(ic,jc) = 0.1*(strvzr(ic,jc)-0.0D0)  &
              -bcl4zr(ic,jc)
          bcl1zr(ic,jc)= half*sorpzr(ic,jc)  &
              + cobczr*acouzr(ic,jc)*(strpzr(ic,jc)-pinfzr)  &
              + 0.5*(1.0-bet)*tt1zr(ic,jc)- bcl1zr(ic,jc) &
!              19 JAN 2024 VM: change making negative
!     +                     + 100.0*STRDZR(IC,JC)*(ONE/ZGDLEN)  &
          - 100.0*strdzr(ic,jc)*(one/zgdlen)  &
              *(acouzr(ic,jc)**2.0-strwzr(ic,jc)**2.0) *(strwzr(ic,jc)-0.0)
        ELSE
          bcl1zr(ic,jc)= half*sorpzr(ic,jc)  &
              + cobczr*acouzr(ic,jc)*(strpzr(ic,jc)-pinfzr)  &
              + 0.5*(1.0-bet)*tt1zr(ic,jc)- bcl1zr(ic,jc)
        END IF
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        IF((strwzr(ic,jc) < zero).AND.(flag_pio_zr==1))THEN
          drhs(ic,jc,kstol) = drhs(ic,jc,kstol)  &
              - bcl1zr(ic,jc)*ova2zr(ic,jc) - bcl2zr(ic,jc)
          
          urhs(ic,jc,kstol) = urhs(ic,jc,kstol)  &
              - bcl1zr(ic,jc)*ova2zr(ic,jc)*struzr(ic,jc)  &
              - bcl2zr(ic,jc)*struzr(ic,jc) - bcl3zr(ic,jc)*strdzr(ic,jc)
          
          vrhs(ic,jc,kstol) = vrhs(ic,jc,kstol)  &
              - bcl1zr(ic,jc)*ova2zr(ic,jc)*strvzr(ic,jc)  &
              - bcl2zr(ic,jc)*strvzr(ic,jc) - bcl4zr(ic,jc)*strdzr(ic,jc)
          
          wrhs(ic,jc,kstol) = wrhs(ic,jc,kstol)  &
              - bcl1zr(ic,jc)*ova2zr(ic,jc)*(strwzr(ic,jc)-acouzr(ic,jc))  &
              - bcl2zr(ic,jc)*strwzr(ic,jc)
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
              - bcl1zr(ic,jc)*(ova2zr(ic,jc)*strezr(ic,jc) &
!              16 JAN 2024: VM+NC sign change
!     +                                     + STRWZR(IC,JC)/ACOUZR(IC,JC)  &
          - strwzr(ic,jc)/acouzr(ic,jc) + ovgmzr(ic,jc))  &
              - bcl2zr(ic,jc)*strezr(ic,jc)  &
              - bcl3zr(ic,jc)*strdzr(ic,jc)*struzr(ic,jc)  &
              - bcl4zr(ic,jc)*strdzr(ic,jc)*strvzr(ic,jc)
        ELSE
          drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc)*ova2zr(ic,jc)
          
          urhs(ic,jc,kstol) = urhs(ic,jc,kstol)  &
              - bcl1zr(ic,jc)*ova2zr(ic,jc)*struzr(ic,jc)
          
          vrhs(ic,jc,kstol) = vrhs(ic,jc,kstol)  &
              - bcl1zr(ic,jc)*ova2zr(ic,jc)*strvzr(ic,jc)
          
          wrhs(ic,jc,kstol) = wrhs(ic,jc,kstol)  &
              - bcl1zr(ic,jc)*ova2zr(ic,jc)*(strwzr(ic,jc)-acouzr(ic,jc))
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
              - bcl1zr(ic,jc)*(ova2zr(ic,jc)*strezr(ic,jc)  &
              - strwzr(ic,jc)/acouzr(ic,jc) + ovgmzr(ic,jc))
        END IF
      END DO
    END DO
    
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          IF((strwzr(ic,jc) < zero).AND.(flag_pio_zr==1))THEN
            fornow = bclyzr(ic,jc,ispec)*strdzr(ic,jc)
            
            erhs(ic,jc,kstol) = erhs(ic,jc,kstol) - fornow*strhzr(ic,jc,ispec)
            
            yrhs(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)  &
                - (bcl2zr(ic,jc)+bcl5zr(ic,jc)*ova2zr(ic,jc))*stryzr(ic,jc,ispec)  &
                - fornow
            
          ELSE
            yrhs(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)  &
                - bcl1zr(ic,jc)*ova2zr(ic,jc)*stryzr(ic,jc,ispec)
          END IF
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
        fornow = strdzr(ic,jc)*acouzr(ic,jc)*bcl1zr(ic,jc)
        bcl1zr(ic,jc) = half*(strwzr(ic,jc)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc)-fornow)
        bcl2zr(ic,jc) = strwzr(ic,jc)  &
            *(bcl2zr(ic,jc)-bcl5zr(ic,jc)*ova2zr(ic,jc))
        bcl3zr(ic,jc) = strwzr(ic,jc)*bcl3zr(ic,jc)
        bcl4zr(ic,jc) = strwzr(ic,jc)*bcl4zr(ic,jc)
        
!             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
!             L5Z UNCHANGED
        bcl1zr(ic,jc) = half*sorpzr(ic,jc)  &
            + cobczr*acouzr(ic,jc)*(strpzr(ic,jc)-pinfzr) - bcl1zr(ic,jc)
        bcl2zr(ic,jc) = -bcl2zr(ic,jc)
        bcl3zr(ic,jc) = -bcl3zr(ic,jc)
        bcl4zr(ic,jc) = -bcl4zr(ic,jc)
        
      END DO
    END DO
    
!         LYZ
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
!               OLD VALUE OF L's
          bclyzr(ic,jc,ispec) = strwzr(ic,jc)*bclyzr(ic,jc,ispec)
          
!               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
          bclyzr(ic,jc,ispec) = ratezr(ic,jc,ispec)/strdzr(ic,jc)  &
              - bclyzr(ic,jc,ispec)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc)*ova2zr(ic,jc)  &
            - bcl2zr(ic,jc)
        
        urhs(ic,jc,kstol) = urhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc)*ova2zr(ic,jc)*struzr(ic,jc)  &
            - bcl2zr(ic,jc)*struzr(ic,jc) - bcl3zr(ic,jc)*strdzr(ic,jc)
        
        vrhs(ic,jc,kstol) = vrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc)*ova2zr(ic,jc)*strvzr(ic,jc)  &
            - bcl2zr(ic,jc)*strvzr(ic,jc) - bcl4zr(ic,jc)*strdzr(ic,jc)
        
        wrhs(ic,jc,kstol) = wrhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc)*ova2zr(ic,jc)*(strwzr(ic,jc)-acouzr(ic,jc))  &
            - bcl2zr(ic,jc)*strwzr(ic,jc)
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc)*(ova2zr(ic,jc)*strezr(ic,jc) &
!              16 JAN 2024: VM+NC sign change
!     +                                     + STRWZR(IC,JC)/ACOUZR(IC,JC)  &
        - strwzr(ic,jc)/acouzr(ic,jc) + ovgmzr(ic,jc))  &
            - bcl2zr(ic,jc)*strezr(ic,jc)  &
            - bcl3zr(ic,jc)*strdzr(ic,jc)*struzr(ic,jc)  &
            - bcl4zr(ic,jc)*strdzr(ic,jc)*strvzr(ic,jc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          fornow = bclyzr(ic,jc,ispec)*strdzr(ic,jc)
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol) - fornow*strhzr(ic,jc,ispec)
          
          yrhs(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)  &
              - (bcl2zr(ic,jc)+bcl5zr(ic,jc)*ova2zr(ic,jc))*stryzr(ic,jc,ispec)  &
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
        
        sydtzr(ic,jc) = sydtzr(ic,jc)/strrzr(ic,jc)
        sorpzr(ic,jc) = -sorpzr(ic,jc)*gam1zr(ic,jc)
        
      END DO
    END DO
    
!         SPECIFY L's AS REQUIRED
!         L1Z,L2Z,L5Z
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             OLD VALUE OF L's
        fornow = strdzr(ic,jc)*acouzr(ic,jc)*bcl1zr(ic,jc)
        bcl1zr(ic,jc) = half*(strwzr(ic,jc)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc)-fornow)
        bcl2zr(ic,jc) = strwzr(ic,jc)  &
            *(bcl2zr(ic,jc)-bcl5zr(ic,jc)*ova2zr(ic,jc))
        bcl5zr(ic,jc) = half*(strwzr(ic,jc)+acouzr(ic,jc))  &
            *(bcl5zr(ic,jc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Z UNCHANGED
        bcl1zr(ic,jc) = bcl5zr(ic,jc)  &
            + strdzr(ic,jc)*acouzr(ic,jc)*dwdtzr(ic,jc) - bcl1zr(ic,jc)
        bcl2zr(ic,jc) = gam1zr(ic,jc)*ova2zr(ic,jc)  &
            *(bcl1zr(ic,jc)+bcl5zr(ic,jc))  &
            + strdzr(ic,jc)*(dtdtzr(ic,jc)/strtzr(ic,jc)  &
            - sorpzr(ic,jc)/strpzr(ic,jc) + sydtzr(ic,jc))  &
            - bcl2zr(ic,jc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc)*ova2zr(ic,jc)  &
            - bcl2zr(ic,jc)
        
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
        fornow = strdzr(ic,jc)*acouzr(ic,jc)*bcl1zr(ic,jc)
        bcl1zr(ic,jc) = half*(strwzr(ic,jc)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc)-fornow)
        bcl2zr(ic,jc) = strwzr(ic,jc)  &
            *(bcl2zr(ic,jc)-bcl5zr(ic,jc)*ova2zr(ic,jc))
        bcl3zr(ic,jc) = strwzr(ic,jc)*bcl3zr(ic,jc)
        bcl4zr(ic,jc) = strwzr(ic,jc)*bcl4zr(ic,jc)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Z UNCHANGED
        fornow = bcl5zr(ic,jc) + strdzr(ic,jc)*acouzr(ic,jc)*dwdtzr(ic,jc)
        bcl1zr(ic,jc) = fornow - bcl1zr(ic,jc)
        bcl2zr(ic,jc) = -dddtzr(ic,jc)  &
            - ova2zr(ic,jc)*(bcl1zr(ic,jc)+fornow) - bcl2zr(ic,jc)
        bcl3zr(ic,jc) = -dudtzr(ic,jc) - bcl3zr(ic,jc)
        bcl4zr(ic,jc) = -dvdtzr(ic,jc) - bcl4zr(ic,jc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc)*(ova2zr(ic,jc)*strezr(ic,jc)  &
            + strwzr(ic,jc)/acouzr(ic,jc) + ovgmzr(ic,jc))  &
            - bcl2zr(ic,jc)*strezr(ic,jc)  &
            - bcl3zr(ic,jc)*strdzr(ic,jc)*struzr(ic,jc)  &
            - bcl4zr(ic,jc)*strdzr(ic,jc)*strvzr(ic,jc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          bclyzr(ic,jc,ispec) = ratezr(ic,jc,ispec)/strdzr(ic,jc)  &
              - dydtzr(ic,jc,ispec) - strwzr(ic,jc)*bclyzr(ic,jc,ispec)
          
          erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
              - bclyzr(ic,jc,ispec)*strdzr(ic,jc)*strhzr(ic,jc,ispec)
          
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
        fornow = strdzr(ic,jc)*acouzr(ic,jc)*bcl1zr(ic,jc)
        bcl1zr(ic,jc) = half*(strwzr(ic,jc)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc)-fornow)
        bcl3zr(ic,jc) = strwzr(ic,jc)*bcl3zr(ic,jc)
        bcl4zr(ic,jc) = strwzr(ic,jc)*bcl4zr(ic,jc)
        bcl5zr(ic,jc) = half*(strwzr(ic,jc)+acouzr(ic,jc))  &
            *(bcl5zr(ic,jc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L2Z,L5Z UNCHANGED
        bcl1zr(ic,jc) = bcl5zr(ic,jc)  &
            + strdzr(ic,jc)*acouzr(ic,jc)*dwdtzr(ic,jc) - bcl1zr(ic,jc)
        bcl3zr(ic,jc) = -dudtzr(ic,jc) - bcl3zr(ic,jc)
        bcl4zr(ic,jc) = -dvdtzr(ic,jc) - bcl4zr(ic,jc)
        
      END DO
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc)*ova2zr(ic,jc)
        
        erhs(ic,jc,kstol) = erhs(ic,jc,kstol)  &
            - bcl1zr(ic,jc)*(ova2zr(ic,jc)*strezr(ic,jc)  &
            + strwzr(ic,jc)/acouzr(ic,jc) + ovgmzr(ic,jc))  &
            - bcl3zr(ic,jc)*strdzr(ic,jc)*struzr(ic,jc)  &
            - bcl4zr(ic,jc)*strdzr(ic,jc)*strvzr(ic,jc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)  &
              - bcl1zr(ic,jc)*ova2zr(ic,jc)*stryzr(ic,jc,ispec)
          
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
        fornow = strdzr(ic,jc)*acouzr(ic,jc)*bcl1zr(ic,jc)
        bcl1zr(ic,jc) = half*(strwzr(ic,jc)-acouzr(ic,jc))  &
            *(bcl5zr(ic,jc)-fornow)
        bcl2zr(ic,jc) = strwzr(ic,jc)  &
            *(bcl2zr(ic,jc)-bcl5zr(ic,jc)*ova2zr(ic,jc))
        bcl3zr(ic,jc) = strwzr(ic,jc)*bcl3zr(ic,jc)
        bcl4zr(ic,jc) = strwzr(ic,jc)*bcl4zr(ic,jc)
        bcl5zr(ic,jc) = half*(strwzr(ic,jc)+acouzr(ic,jc))  &
            *(bcl5zr(ic,jc)+fornow)
        
!             SUBTRACT FROM NEW VALUE OF L's
!             L5Z UNCHANGED
        bcl1zr(ic,jc) = bcl5zr(ic,jc)  &
            + strdzr(ic,jc)*acouzr(ic,jc)*dwdtzr(ic,jc) - bcl1zr(ic,jc)
        bcl3zr(ic,jc) = -dudtzr(ic,jc) - bcl3zr(ic,jc)
        bcl4zr(ic,jc) = -dvdtzr(ic,jc) - bcl4zr(ic,jc)
        bcl2zr(ic,jc) = gam1zr(ic,jc)*ova2zr(ic,jc)  &
            *(bcl1zr(ic,jc)+bcl5zr(ic,jc))  &
            + strdzr(ic,jc)*(dtdtzr(ic,jc)/strtzr(ic,jc)  &
            - sorpzr(ic,jc)/strpzr(ic,jc)) - bcl2zr(ic,jc)
        
      END DO
    END DO
    
!         LYZ
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
!               OLD VALUE OF LYZ
          bclyzr(ic,jc,ispec) = strwzr(ic,jc)*bclyzr(ic,jc,ispec)
          
!               UPDATE L2Z
          bcl2zr(ic,jc) = bcl2zr(ic,jc) + (ratezr(ic,jc,ispec)  &
              - strdzr(ic,jc)*bclyzr(ic,jc,ispec)) *rgspec(ispec)/strrzr(ic,jc)
          
        END DO
      END DO
      
    END DO
    
!         ADD TO CONSERVATIVE SOURCE TERMS
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        drhs(ic,jc,kstol) = drhs(ic,jc,kstol) - bcl1zr(ic,jc)*ova2zr(ic,jc)  &
            - bcl2zr(ic,jc)
        
      END DO
    END DO
    
    DO ispec = 1,nspec
      
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yrhs(ic,jc,kstol,ispec) = yrhs(ic,jc,kstol,ispec)  &
              - (bcl2zr(ic,jc)+bcl1zr(ic,jc)*ova2zr(ic,jc))*stryzr(ic,jc,ispec)
          
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
