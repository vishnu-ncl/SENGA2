SUBROUTINE temper
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:26:27

!     *************************************************************************

!     TEMPER
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     16-NOV-2002:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     COMPUTES TEMPERATURE AND PRESSURE

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     PARAMETERS
!     ==========
real(kind=8):: toltmp
PARAMETER(toltmp = 1.0E-10)
INTEGER :: ntitrs
PARAMETER(ntitrs = 100)


!     LOCAL DATA
!     ==========
real(kind=8):: tcoeff(0:nctmax),tderiv(1:nctmax)
real(kind=8):: ukuk
real(kind=8):: tempor,tfpoly,tdpoly,deltmp,cpfory
INTEGER :: ic,jc,kc,ispec,itint,icp,ititrs
INTEGER :: iindex,ipower,icoef1,icoef2


!     BEGIN
!     =====

!     =========================================================================

!     TEMPERATURE AND PRESSURE
!     ------------------------

!     TEMPERATURE AND PRESSURE ARE PARALLEL

DO kc = kstalt,kstolt
  DO jc = jstalt,jstolt
    DO ic = istalt,istolt
      
!           ===================================================================
      
!           KINETIC ENERGY TERM
      ukuk = (urhs(ic,jc,kc)*urhs(ic,jc,kc)  &
          +  vrhs(ic,jc,kc)*vrhs(ic,jc,kc)  &
          +  wrhs(ic,jc,kc)*wrhs(ic,jc,kc))/drhs(ic,jc,kc)
      
!           ===================================================================
      
!           INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
!           AND ITS DERIVATIVE
      tcoeff(0) = half*ukuk - erhs(ic,jc,kc)
      DO icp = 1, nctmax
        tcoeff(icp) = zero
        tderiv(icp) = zero
      END DO
      
!           ===================================================================
      
!           USE STORE7 TO ACCUMULATE MIXTURE SPECIFIC GAS CONSTANT
!           INITIALISE STORE7
      store7(ic,jc,kc) = zero
      
!           ===================================================================
      
!           RUN THROUGH ALL SPECIES
      DO ispec = 1,nspec
        
!             =================================================================
        
!             LOCATE TEMPERATURE IN AN INTERVAL
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1
        icoef2 = ntbase**ipower
        icoef1 = icoef2*ntbase
        itint = 1 + MOD(itndex(ic,jc,kc,iindex),icoef1)/icoef2
        
!             =================================================================
        
!             CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
        tcoeff(0) = tcoeff(0) + yrhs(ic,jc,kc,ispec)*  &
            amascp(ncenth(itint,ispec),itint,ispec)
        tcoeff(1) = tcoeff(1) + yrhs(ic,jc,kc,ispec)*amasct(1,itint,ispec)
        tderiv(1) = tcoeff(1)
        DO icp = 2, ncpoly(itint,ispec)
          tcoeff(icp) = tcoeff(icp) + yrhs(ic,jc,kc,ispec)*  &
              amasct(icp,itint,ispec)
          tderiv(icp) = tderiv(icp) + yrhs(ic,jc,kc,ispec)*  &
              amascp(icp,itint,ispec)
        END DO
        
!             =================================================================
        
!             USE STORE7
!             TO ACCUMULATE (DENSITY TIMES) MIXTURE SPECIFIC GAS CONSTANT
        store7(ic,jc,kc) = store7(ic,jc,kc)  &
            + yrhs(ic,jc,kc,ispec)*rgspec(ispec)
        
!             =================================================================
        
      END DO
!           END OF RUN THROUGH ALL SPECIES
      
!           ===================================================================
      
!           SOLVE FOR TEMPERATURE
!           USING NEWTON-RAPHSON
      tempor = trun(ic,jc,kc)
      ititrs = 1
      1000        CONTINUE
      
!             EVALUATE TEMPERATURE POLYNOMIAL AND ITS DERIVATIVE
      tfpoly = tcoeff(nctmax)
      tdpoly = tderiv(nctmax)
      DO icp = nctmm1,1,-1
        tfpoly = tcoeff(icp) + tfpoly*tempor
        tdpoly = tderiv(icp) + tdpoly*tempor
      END DO
      tfpoly = tcoeff(0) + tfpoly*tempor
      
!             EVALUATE TEMPERATURE CORRECTION
      deltmp = -tfpoly/tdpoly
      
!             CHECK FOR CONVERGENCE
      IF(ABS(deltmp) > toltmp)THEN
        IF(ititrs < ntitrs)THEN
          tempor = tempor + deltmp
          ititrs = ititrs + 1
          GO TO 1000
        ELSE
          WRITE(6,*) 'Fatal: TEMPER: T iteration failed to converge'
          WRITE(6,*)'processor:',iproc
          WRITE(6,*)'at point:',ic,jc,kc
          WRITE(6,*)'with values:',tempor,deltmp
          WRITE(6,*)drhs(ic,jc,kc)
          WRITE(6,*)urhs(ic,jc,kc)
          WRITE(6,*)vrhs(ic,jc,kc)
          WRITE(6,*)wrhs(ic,jc,kc)
          WRITE(6,*)erhs(ic,jc,kc)
          STOP
        END IF
      END IF
      
!           END OF LOOP 1000
      
!           ===================================================================
      
!           SET THE NEW TEMPERATURE
      trun(ic,jc,kc) = tempor
      
!           ===================================================================
      
!           FOR ALL SPECIES RELOCATE TEMPERATURE IN AN INTERVAL
!           EVALUATE MIXTURE SPECIFIC HEAT CP
      DO iindex = 1,nintmx
        itndex(ic,jc,kc,iindex) = 0
      END DO
      transp(ic,jc,kc) = zero
      DO ispec = 1,nspec
        
        itint = 1
        1100          CONTINUE
        IF(trun(ic,jc,kc) > tinthi(itint,ispec))THEN
          IF(itint < ntint(ispec))THEN
            itint = itint + 1
            GO TO 1100
          END IF
        END IF
!             END OF LOOP 1100
        
!             SET THE TEMPERATURE INTERVAL INDEX
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1
        itndex(ic,jc,kc,iindex) = itndex(ic,jc,kc,iindex)  &
            +(itint-1)*ntbase**ipower
        
!             =================================================================
        
!             EVALUATE MIXTURE SPECIFIC HEAT CP
        cpfory = amascp(ncpoly(itint,ispec),itint,ispec)
        DO icp = ncpom1(itint,ispec),1,-1
          cpfory = cpfory*trun(ic,jc,kc) + amascp(icp,itint,ispec)
        END DO
        transp(ic,jc,kc) = transp(ic,jc,kc) + yrhs(ic,jc,kc,ispec)*cpfory
        
      END DO
      transp(ic,jc,kc) = transp(ic,jc,kc)/drhs(ic,jc,kc)
      
!           ===================================================================
      
!           EVALUATE MIXTURE PRESSURE
      prun(ic,jc,kc) = trun(ic,jc,kc)*store7(ic,jc,kc)
      
!           ===================================================================
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE temper
