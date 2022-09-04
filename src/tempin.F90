SUBROUTINE tempin
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-04  Time: 20:55:58

!     *************************************************************************

!     TEMPIN
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     15-MAY-2004:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     INITIALISES TEMPERATURE AND PRESSURE
!     USES BISECTION METHOD FOR ROBUSTNESS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     PARAMETERS
!     ==========
REAL(KIND=8) :: toltmp
PARAMETER(toltmp = 0.00010_8)
REAL(KIND=8) :: tininc
PARAMETER(tininc = 50.0_8)
REAL(KIND=8) :: tlimlo,tlimhi
PARAMETER(tlimlo = 200.0_8, tlimhi = 3000.0_8)


!     LOCAL DATA
!     ==========
REAL(KIND=8) :: tcoeff(0:nctmax)
REAL(KIND=8) :: ukuk
REAL(KIND=8) :: tempor,tupper,tlower,tresid,tuk2me,cpfory
INTEGER :: ic,jc,kc,ispec,itint,icp
INTEGER :: iindex,ipower
LOGICAL :: fnconv


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
      
      tuk2me = half*ukuk - erhs(ic,jc,kc)
      
!           ===================================================================
      
!           SOLVE FOR TEMPERATURE
!           =====================
!           USING BISECTION
      tlower = tlimlo
      tupper = tlimhi
      
!           SET NON-CONVERGENCE FLAG
      fnconv = .true.
      
!           SET TEMPERATURE FROM INITIAL GUESS
      tempor = trin
      
!           ===================================================================
      
!           FIRST BRACKET THE ROOT
!           ======================
      
!           -------------------------------------------------------------------
      
!           INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
      tcoeff(0) = tuk2me
      DO icp = 1, nctmax
        tcoeff(icp) = zero
      END DO
      
!           FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
      DO ispec = 1,nspec
        
        itint = 1
        1100          CONTINUE
        IF(tempor > tinthi(itint,ispec))THEN
          IF(itint < ntint(ispec))THEN
            itint = itint + 1
            GO TO 1100
          END IF
        END IF
!             END OF LOOP 1100
        
!             CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
        tcoeff(0) = tcoeff(0) + yrhs(ic,jc,kc,ispec)*  &
            amascp(ncenth(itint,ispec),itint,ispec)
        
        tcoeff(1) = tcoeff(1) + yrhs(ic,jc,kc,ispec)*amasct(1,itint,ispec)
        DO icp = 2, ncpoly(itint,ispec)
          tcoeff(icp) = tcoeff(icp)  &
              + yrhs(ic,jc,kc,ispec)*amasct(icp,itint,ispec)
        END DO
        
      END DO
!           END OF RUN THROUGH ALL SPECIES
      
!           -------------------------------------------------------------------
      
!           EVALUATE TEMPERATURE RESIDUAL
      tresid = tcoeff(nctmax)
      DO icp = nctmm1,1,-1
        tresid = tcoeff(icp) + tresid*tempor
      END DO
      tresid = tcoeff(0) + tresid*tempor
      
!           -------------------------------------------------------------------
      
!           CHECK INITIAL GUESS FOR CONVERGENCE
      IF(ABS(tresid) < toltmp)THEN
        
!             -----------------------------------------------------------------
        
!             CONVERGED ON FIRST PASS
        fnconv = .false.
        
!             -----------------------------------------------------------------
        
      ELSE IF(tresid < zero)THEN
        
!             -----------------------------------------------------------------
        
!             INITIAL GUESS IS TOO LOW
        2000          CONTINUE
        
!               SET INITIAL GUESS AS LOWER LIMIT AND TRY AGAIN
        tlower = tempor
        tempor = tempor + tininc
        
!               ---------------------------------------------------------------
        
!               INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
        tcoeff(0) = tuk2me
        DO icp = 1, nctmax
          tcoeff(icp) = zero
        END DO
        
!               FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
        DO ispec = 1,nspec
          
          itint = 1
          2100              CONTINUE
          IF(tempor > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 2100
            END IF
          END IF
!                 END OF LOOP 2100
          
!                 CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
          tcoeff(0) = tcoeff(0) + yrhs(ic,jc,kc,ispec)*  &
              amascp(ncenth(itint,ispec),itint,ispec)
          tcoeff(1) = tcoeff(1) + yrhs(ic,jc,kc,ispec)*amasct(1,itint,ispec)
          DO icp = 2, ncpoly(itint,ispec)
            tcoeff(icp) = tcoeff(icp)  &
                + yrhs(ic,jc,kc,ispec)*amasct(icp,itint,ispec)
          END DO
          
        END DO
!               END OF RUN THROUGH ALL SPECIES
        
!               ---------------------------------------------------------------
        
!               EVALUATE TEMPERATURE RESIDUAL
        tresid = tcoeff(nctmax)
        DO icp = nctmm1,1,-1
          tresid = tcoeff(icp) + tresid*tempor
        END DO
        tresid = tcoeff(0) + tresid*tempor
        
!               ---------------------------------------------------------------
        
!               CHECK NEW GUESS FOR CONVERGENCE
        IF(ABS(tresid) < toltmp)THEN
          
!                 -------------------------------------------------------------
          
!                 NEW GUESS HAS CONVERGED
          fnconv = .false.
          
!                 -------------------------------------------------------------
          
        ELSE IF(tresid < zero)THEN
          
!                 -------------------------------------------------------------
          
!                 NEW GUESS IS STILL TOO LOW: GO ROUND AGAIN
          IF(tempor < tlimhi)THEN
            GO TO 2000
          ELSE
            WRITE(6,*) 'Fatal: TEMPIN: T upper bracket failed to converge'
            WRITE(6,*)'processor:',iproc
            WRITE(6,*)'at point:',ic,jc,kc
            WRITE(6,*)'with values:',tempor,tresid
            WRITE(6,*)drhs(ic,jc,kc)
            WRITE(6,*)urhs(ic,jc,kc)
            WRITE(6,*)vrhs(ic,jc,kc)
            WRITE(6,*)wrhs(ic,jc,kc)
            WRITE(6,*)erhs(ic,jc,kc)
            DO ispec = 1, nspec
              WRITE(6,*)yrhs(ic,jc,kc,ispec)
            END DO
            STOP
          END IF
          
!                 -------------------------------------------------------------
          
        ELSE IF(tresid > zero)THEN
          
!                 -------------------------------------------------------------
          
!                 ROOT IS BRACKETED
          tupper = tempor
          
!                 -------------------------------------------------------------
          
        END IF
        
!               ---------------------------------------------------------------
        
!             END OF LOOP 2000
        
!             -----------------------------------------------------------------
        
      ELSE IF(tresid > zero)THEN
        
!             -----------------------------------------------------------------
        
!             INITIAL GUESS IS TOO HIGH
        3000          CONTINUE
        
!               SET INITIAL GUESS AS UPPER LIMIT AND TRY AGAIN
        tupper = tempor
        tempor = tempor - tininc
        
!               ---------------------------------------------------------------
        
!               INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
        tcoeff(0) = tuk2me
        DO icp = 1, nctmax
          tcoeff(icp) = zero
        END DO
        
!               FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
        DO ispec = 1,nspec
          
          itint = 1
          3100              CONTINUE
          IF(tempor > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 3100
            END IF
          END IF
!                 END OF LOOP 3100
          
!                 CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
          tcoeff(0) = tcoeff(0) + yrhs(ic,jc,kc,ispec)*  &
              amascp(ncenth(itint,ispec),itint,ispec)
          tcoeff(1) = tcoeff(1) + yrhs(ic,jc,kc,ispec)*amasct(1,itint,ispec)
          DO icp = 2, ncpoly(itint,ispec)
            tcoeff(icp) = tcoeff(icp)  &
                + yrhs(ic,jc,kc,ispec)*amasct(icp,itint,ispec)
          END DO
          
        END DO
!               END OF RUN THROUGH ALL SPECIES
        
!               ---------------------------------------------------------------
        
!               EVALUATE TEMPERATURE RESIDUAL
        tresid = tcoeff(nctmax)
        DO icp = nctmm1,1,-1
          tresid = tcoeff(icp) + tresid*tempor
        END DO
        tresid = tcoeff(0) + tresid*tempor
        
!               ---------------------------------------------------------------
        
!               CHECK NEW GUESS FOR CONVERGENCE
        IF(ABS(tresid) < toltmp)THEN
          
!                 -------------------------------------------------------------
          
!                 NEW GUESS HAS CONVERGED
          fnconv = .false.
          
!                 -------------------------------------------------------------
          
        ELSE IF(tresid > zero)THEN
          
!                 -------------------------------------------------------------
          
!                 NEW GUESS IS STILL TOO HIGH: GO ROUND AGAIN
          IF(tempor > tlimlo)THEN
            GO TO 3000
          ELSE
            WRITE(6,*) 'Fatal: TEMPIN: T lower bracket failed to converge'
            WRITE(6,*)'processor:',iproc
            WRITE(6,*)'at point:',ic,jc,kc
            WRITE(6,*)'with values:',tempor,tresid
            WRITE(6,*)drhs(ic,jc,kc)
            WRITE(6,*)urhs(ic,jc,kc)
            WRITE(6,*)vrhs(ic,jc,kc)
            WRITE(6,*)wrhs(ic,jc,kc)
            WRITE(6,*)erhs(ic,jc,kc)
            DO ispec = 1, nspec
              WRITE(6,*)yrhs(ic,jc,kc,ispec)
            END DO
            STOP
          END IF
          
!                 -------------------------------------------------------------
          
        ELSE IF(tresid < zero)THEN
          
!                 -------------------------------------------------------------
          
!                 ROOT IS BRACKETED
          tlower = tempor
          
!                 -------------------------------------------------------------
          
        END IF
        
!               ---------------------------------------------------------------
        
!             END OF LOOP 3000
        
!             -----------------------------------------------------------------
        
      END IF
!           END OF CHECK INITIAL GUESS FOR CONVERGENCE
      
!           ===================================================================
      
!           ROOT IS BRACKETED
!           =================
!           NOW USE BISECTION TO REFINE THE ROOT
      
      IF(fnconv)THEN
        
        4000          CONTINUE
        
!               BISECT
        tempor = half*(tlower+tupper)
        
!               ---------------------------------------------------------------
        
!               INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
        tcoeff(0) = tuk2me
        DO icp = 1, nctmax
          tcoeff(icp) = zero
        END DO
        
!               FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
        DO ispec = 1,nspec
          
          itint = 1
          4100              CONTINUE
          IF(tempor > tinthi(itint,ispec))THEN
            IF(itint < ntint(ispec))THEN
              itint = itint + 1
              GO TO 4100
            END IF
          END IF
!                 END OF LOOP 4100
          
!                 CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
          tcoeff(0) = tcoeff(0) + yrhs(ic,jc,kc,ispec)*  &
              amascp(ncenth(itint,ispec),itint,ispec)
          tcoeff(1) = tcoeff(1) + yrhs(ic,jc,kc,ispec)*amasct(1,itint,ispec)
          DO icp = 2, ncpoly(itint,ispec)
            tcoeff(icp) = tcoeff(icp)  &
                + yrhs(ic,jc,kc,ispec)*amasct(icp,itint,ispec)
          END DO
          
        END DO
!               END OF RUN THROUGH ALL SPECIES
        
!               ---------------------------------------------------------------
        
!               EVALUATE TEMPERATURE RESIDUAL
        tresid = tcoeff(nctmax)
        DO icp = nctmm1,1,-1
          tresid = tcoeff(icp) + tresid*tempor
        END DO
        tresid = tcoeff(0) + tresid*tempor
        
!               ---------------------------------------------------------------
        
        IF(ABS(tresid) < toltmp)THEN
          
!                 CONVERGED
          trun(ic,jc,kc) = tempor
          
        ELSE IF(tresid < zero)THEN
          
          tlower = tempor
          GO TO 4000
          
        ELSE IF(tresid > zero)THEN
          
          tupper = tempor
          GO TO 4000
          
        END IF
        
!             -----------------------------------------------------------------
!             END OF LOOP 4000
        
      END IF
!           END OF BISECTION
      
!           ===================================================================
      
!           SET THE NEW TEMPERATURE
      trun(ic,jc,kc) = tempor
      
!           ===================================================================
      
!           CONSTRUCT THE TEMPERATURE INTERVAL INDEX
!           EVALUATE PRESSURE
!           EVALUATE MIXTURE SPECIFIC HEAT CP
      DO iindex = 1,nintmx
        itndex(ic,jc,kc,iindex) = 0
      END DO
      store7(ic,jc,kc) = zero
      transp(ic,jc,kc) = zero
      DO ispec = 1,nspec
        
        itint = 1
        5100          CONTINUE
        IF(trun(ic,jc,kc) > tinthi(itint,ispec))THEN
          IF(itint < ntint(ispec))THEN
            itint = itint + 1
            GO TO 5100
          END IF
        END IF
!             END OF LOOP 5100
        
!             SET THE TEMPERATURE INDEX
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1
        itndex(ic,jc,kc,iindex) = itndex(ic,jc,kc,iindex)  &
            + (itint-1)*ntbase**ipower
        
!             =================================================================
        
!             EVALUATE MIXTURE SPECIFIC HEAT CP
        cpfory = amascp(ncpoly(itint,ispec),itint,ispec)
        DO icp = ncpom1(itint,ispec),1,-1
          cpfory = cpfory*trun(ic,jc,kc) + amascp(icp,itint,ispec)
        END DO
        transp(ic,jc,kc) = transp(ic,jc,kc) + yrhs(ic,jc,kc,ispec)*cpfory
        
!             =================================================================
        
!             EVALUATE (DENSITY TIMES) MIXTURE GAS CONSTANT FOR PRESSURE
        store7(ic,jc,kc) = store7(ic,jc,kc)  &
            + yrhs(ic,jc,kc,ispec)*rgspec(ispec)
        
!             =================================================================
        
      END DO
      transp(ic,jc,kc) = transp(ic,jc,kc)/drhs(ic,jc,kc)
      
!           ===================================================================
      
!           EVALUATE PRESSURE
      prun(ic,jc,kc) = trun(ic,jc,kc)*store7(ic,jc,kc)
      
!           ===================================================================
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE tempin
