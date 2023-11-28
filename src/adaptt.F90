SUBROUTINE adaptt
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:58:13

!     *************************************************************************

!     ADAPTT
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     19-JAN-2003:  CREATED
!     23-AUG-2009:  RSC REVISE ERROR NORM EVALUATION
!     08-AUG-2012:  RSC EVALUATE ALL SPECIES
!     09-AUG-2012   RSC/RACG USE GLOBAL ERROR

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     COMPUTES NEW TIMESTEP FOR ERK SCHEME

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
REAL(kind=8) :: erytot(nspcmx)
REAL(kind=8) :: erdtot,erutot,ervtot,erwtot,eretot
REAL(kind=8) :: errmax,tratio,tstold
!     RSC/RACG 09-AUG-2012 USE GLOBAL ERROR
!      REAL(kind=8) TSTLOC
REAL(kind=8) :: errloc
REAL(kind=8) :: fornow
INTEGER :: ic,jc,kc,ispec


!     BEGIN
!     =====

!     =========================================================================

!     CHECK ADAPTION FLAG
!     -------------------
IF(fladpt)THEN
  
!       =======================================================================
  
!       INITIALISE THE ERROR NORM TOTALS
!       --------------------------------
  erdtot = zero
  erutot = zero
  ervtot = zero
  erwtot = zero
  eretot = zero
!       RSC 08-AUG-2012 EVALUATE ALL SPECIES
!        DO ISPEC = 1,NSPM1
  DO ispec = 1,nspec
    erytot(ispec) = zero
  END DO
  
!       =======================================================================
  
!       ERK ERROR EVALUATION
!       --------------------
!       USING ERK ERROR ARRAYS
!       RSC 23 AUG-2009 REVISE ERROR NORM EVALUATION
  
!       EVALUATE ERROR NORMS
  DO kc = kstald,kstold
    DO jc = jstald,jstold
      DO ic = istald,istold
        
!              FORNOW = ABS(DERR(IC,JC,KC))
        fornow = ABS(derr(ic,jc,kc))/(ABS(drun(ic,jc,kc))+erdnrm)
        IF(fornow > erdtot)erdtot = fornow
        
      END DO
    END DO
  END DO
  
  DO kc = kstalu,kstolu
    DO jc = jstalu,jstolu
      DO ic = istalu,istolu
        
!              FORNOW = ABS(UERR(IC,JC,KC))
        fornow = ABS(uerr(ic,jc,kc))/(ABS(urun(ic,jc,kc))+erunrm)
        IF(fornow > erutot)erutot = fornow
        
      END DO
    END DO
  END DO
  
  DO kc = kstalv,kstolv
    DO jc = jstalv,jstolv
      DO ic = istalv,istolv
        
!              FORNOW = ABS(VERR(IC,JC,KC))
        fornow = ABS(verr(ic,jc,kc))/(ABS(vrun(ic,jc,kc))+ervnrm)
        IF(fornow > ervtot)ervtot = fornow
        
      END DO
    END DO
  END DO
  
  DO kc = kstalw,kstolw
    DO jc = jstalw,jstolw
      DO ic = istalw,istolw
        
!              FORNOW = ABS(WERR(IC,JC,KC))
        fornow = ABS(werr(ic,jc,kc))/(ABS(wrun(ic,jc,kc))+erwnrm)
        IF(fornow > erwtot)erwtot = fornow
        
      END DO
    END DO
  END DO
  
  DO kc = kstale,kstole
    DO jc = jstale,jstole
      DO ic = istale,istole
        
!              FORNOW = ABS(EERR(IC,JC,KC))
        fornow = ABS(eerr(ic,jc,kc))/(ABS(erun(ic,jc,kc))+erenrm)
        IF(fornow > eretot)eretot = fornow
        
      END DO
    END DO
  END DO
  
!       RSC 08-AUG-2012 EVALUATE ALL SPECIES
!        DO ISPEC = 1,NSPM1
  DO ispec = 1,nspec
    
    DO kc = kstaly,kstoly
      DO jc = jstaly,jstoly
        DO ic = istaly,istoly
          
!                FORNOW = ABS(YERR(IC,JC,KC,ISPEC))
          fornow = ABS(yerr(ic,jc,kc,ispec))  &
              /(ABS(yrun(ic,jc,kc,ispec))+erynrm(ispec))
          IF(fornow > erytot(ispec))erytot(ispec) = fornow
          
        END DO
      END DO
    END DO
    
  END DO
  
!       =======================================================================
  
!C       NORMALISE THE ERROR NORMS
!C       -------------------------
!        ERDTOT = ABS(ERDTOT)*ERDNRM
!        ERUTOT = ABS(ERUTOT)*ERUNRM
!        ERVTOT = ABS(ERVTOT)*ERVNRM
!        ERWTOT = ABS(ERWTOT)*ERWNRM
!        ERETOT = ABS(ERETOT)*ERENRM
!        DO ISPEC = 1,NSPM1
!          ERYTOT(ISPEC) = ABS(ERYTOT(ISPEC))*ERYNRM(ISPEC)
!        ENDDO
  
!       =======================================================================
  
!       FIND THE MAXIMUM
!       ----------------
  errmax = zero
  IF(erdtot > errmax)THEN
    errmax = erdtot
    inderr = -4
  END IF
  IF(erutot > errmax)THEN
    errmax = erutot
    inderr = -1
  END IF
  IF(ervtot > errmax)THEN
    errmax = ervtot
    inderr = -2
  END IF
  IF(erwtot > errmax)THEN
    errmax = erwtot
    inderr = -3
  END IF
  IF(eretot > errmax)THEN
    errmax = eretot
    inderr = 0
  END IF
!       RSC 08-AUG-2012 EVALUATE ALL SPECIES
!        DO ISPEC = 1,NSPM1
  DO ispec = 1,nspec
    IF(erytot(ispec) > errmax)THEN
      errmax = erytot(ispec)
      inderr = ispec
    END IF
  END DO
  
!       =======================================================================
  
!       FIND THE LARGEST GLOBAL ERROR
!       -----------------------------
!       RSC/RACG 09-AUG-2012 USE GLOBAL ERROR
  errloc = errmax
  CALL p_gmax(errloc,errmax)
  
!       =======================================================================
  
!       EVALUATE THE NEW TIME STEP
!       --------------------------
!       ZERO CHECK
  IF(errmax < errlow)errmax = errlow
  
!       -----------------------------------------------------------------------
  
!C       I-CONTROLLER
!        TRATIO = CTMULT*EXP(CTALPH*LOG(ERRTOL/ERRMAX))
  
!       -----------------------------------------------------------------------
  
!C       PI-CONTROLLER
!        TRATIO = CTMULT*EXP(CTALPH*LOG(ERRTOL/ERRMAX)
!     +                     +CTBETA*LOG(ERROLD/ERRTOL))
!        ERROLD = ERRMAX
  
!       -----------------------------------------------------------------------
  
!       PID-CONTROLLER
  tratio = ctmult*EXP(ctalph*LOG(errtol/errmax) +ctbeta*LOG(errold/errtol)  &
      +ctgama*LOG(errtol/errldr))
  errldr = errold
  errold = errmax
  
!       -----------------------------------------------------------------------
  
!       LIMIT CHANGES TO TIME STEP
  IF(tratio > trmax)tratio = trmax
  IF(tratio < trmin)tratio = trmin
  
!       -----------------------------------------------------------------------
  
!       SAVE THE OLD TIME STEP
  tstold = tstep
  
!       SET THE NEW TIME STEP
  tstep = tstep*tratio
  
!       LIMIT THE TIME STEP
  IF(tstep > tsmax)tstep = tsmax
  IF(tstep < tsmin)tstep = tsmin
  
!       =======================================================================
  
!       PARALLEL TRANSFER TO SET NEW GLOBAL TIME STEP
!       ---------------------------------------------
!       NEW TIME STEP IS THE GLOBAL MINIMUM OVER ALL PROCESSORS
!       RSC/RACG 09-AUG-2012 USE GLOBAL ERROR
!        TSTLOC = TSTEP
!        CALL P_GMIN(TSTLOC,TSTEP)
  
!       =======================================================================
  
!       UPDATE THE TIME ADVANCEMENT COEFFICIENTS
!       AND THE RK SUBSTEP TIME LEVELS
  tratio = tstep/tstold
  DO irkstp = 1, nrkstp
    rklhs(irkstp) = rklhs(irkstp)*tratio
    rkrhs(irkstp) = rkrhs(irkstp)*tratio
    rkerr(irkstp) = rkerr(irkstp)*tratio
    rktim(irkstp) = rktim(irkstp)*tratio
  END DO
  
!       =======================================================================
  
!       (RE)INITIALISE ERK ERROR ARRAYS
!       -------------------------------
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
        derr(ic,jc,kc) = zero
        uerr(ic,jc,kc) = zero
        verr(ic,jc,kc) = zero
        werr(ic,jc,kc) = zero
        eerr(ic,jc,kc) = zero
        
      END DO
    END DO
  END DO
!       RSC 08-AUG-2012 EVALUATE ALL SPECIES
!        DO ISPEC = 1,NSPM1
  DO ispec = 1,nspec
    
    DO kc = kstal,kstol
      DO jc = jstal,jstol
        DO ic = istal,istol
          
          yerr(ic,jc,kc,ispec) = zero
          
        END DO
      END DO
    END DO
    
  END DO
  
!       =======================================================================
  
!       (RE)INITIALISE ERK SUBSTEP ERROR NORMS
!       --------------------------------------
  DO irkstp = 1,nrkstp
    
    erdrhs(irkstp) = zero
    erurhs(irkstp) = zero
    ervrhs(irkstp) = zero
    erwrhs(irkstp) = zero
    ererhs(irkstp) = zero
!         RSC 08-AUG-2012 EVALUATE ALL SPECIES
!          DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec
      eryrhs(ispec,irkstp) = zero
    END DO
    
  END DO
  
!       =======================================================================
  
END IF
!     ADAPTION FLAG

!     =========================================================================


RETURN
END SUBROUTINE adaptt
