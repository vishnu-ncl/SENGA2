      SUBROUTINE DFMSTR
 
C     *************************************************************************
C
C     DFMSTR
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     10-NOV-2013:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIALISES MESH STRETCHING
C     IN ONE DIMENSION ONLY
C     WITH SPECIFIED MAPPING FUNCTION FOR STRETCH RATE: TANH
C     ONE-SIDED (LEFT)
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION DHATDG(NSZMAX),D2HDG2(NSZMAX)
      DOUBLE PRECISION GCMFAC,GCMHOF,GCMFIX,GCMCOF,GCMOFF,DELTAG
      DOUBLE PRECISION FORNOW
      INTEGER IGSTAL,IGSTOL,IGOFST
      INTEGER IC,ICPROC


C     BEGIN
C     =====
      
C     =========================================================================

C     -------------------------------------------------------------------------

C     SPECIFY THE DIMENSION FOR COORDINATE STRETCHING
C     Y-DIRECTION
      IGOFST = 0
      DO ICPROC = 0, IYPROC-1
        IGOFST = IGOFST + NPMAPY(ICPROC)
      ENDDO
      IGSTAL = IGOFST + 1
      IGSTOL = IGOFST + NPMAPY(IYPROC)
      DELTAG = DELTAY

C     SPECIFY THE SCALE FACTOR FOR STRETCH RATE
      GCMFAC = 1.0D2

C     SPECIFY THE OFFSET FOR STRETCH RATE AT THE ORIGIN 
      GCMHOF = 5.0D-1 

C     SPECIFY THE FIXED POINT FOR NORMALISATION
      GCMFIX = YGDLEN

C     SPECIFY THE STRETCHING FUNCTION FOR MESH SPACING
C     HYPERBOLIC TANGENT => TANH

C     HENCE:
C     INTEGRAL OF STRETCHING FUNCTION FOR STRETCHED COORDINATE LOCATION
C     HYPERBOLIC TANGENT => LN(COSH)

C     DERIVATIVE OF STRETCHING FUNCTION FOR SECOND DERIVATIVE EVALUATION
C     HYPERBOLIC TANGENT => SECH^2

C     -------------------------------------------------------------------------

C     NORMALISE THE STRETCHED COORDINATES
C     NOTE THAT THE NORMALISATION IS INHERENTLY GLOBAL
      GCMOFF = LOG(COSH(GCMHOF))
      GCMCOF = LOG(COSH(GCMFAC*GCMFIX + GCMHOF)) - GCMOFF
      GCMCOF = GCMFAC*GCMFIX/GCMCOF

C     REGULAR COORDINATES
      DO IC = IGSTAL,IGSTOL
        GCMREG(IC) = REAL(IC-1)*DELTAG
      ENDDO     

C     SPECIFIED STRETCHING FUNCTION FOR MESH SPACING
C     HYPERBOLIC TANGENT
      DO IC = IGSTAL,IGSTOL
        DHATDG(IC) = GCMCOF*TANH(GCMFAC*GCMREG(IC) + GCMHOF)
        DGDHAT(IC) = ONE/DHATDG(IC)
        DGDHSQ(IC) = DGDHAT(IC)*DGDHAT(IC)
      ENDDO     

C     INTEGRAL OF STRETCHING FUNCTION FOR STRETCHED COORDINATE LOCATION
C     HYPERBOLIC TANGENT => LN(COSH)
      DO IC = IGSTAL,IGSTOL
        FORNOW = LOG(COSH(GCMFAC*GCMREG(IC) + GCMHOF))- GCMOFF
        GCMSTR(IC) = FORNOW*GCMCOF/GCMFAC
      ENDDO     

C     DERIVATIVE OF STRETCHING FUNCTION FOR SECOND DERIVATIVE EVALUATION
C     HYPERBOLIC TANGENT => SECH^2
      DO IC = IGSTAL,IGSTOL
        FORNOW = COSH(GCMFAC*GCMREG(IC) + GCMHOF)
        FORNOW = FORNOW*FORNOW
        D2HDG2(IC) = GCMCOF*GCMFAC/FORNOW
      ENDDO     

C     HYPERBOLIC COTANGENT => COSECH^2
      DO IC = IGSTAL,IGSTOL
        FORNOW = SINH(GCMFAC*GCMREG(IC) + GCMHOF)
        FORNOW = -GCMCOF*FORNOW*FORNOW
        D2GDH2(IC) = GCMFAC/FORNOW
      ENDDO     

C     DIAGNOSTICS
      WRITE(6,'(2I5)')IGSTAL,IGSTOL
      WRITE(6,'(5(1PE15.7))')GCMFIX,GCMFAC,GCMHOF,GCMOFF,GCMCOF
C      DO IC = IGSTAL,IGSTOL
C        WRITE(6,'(6(1PE12.4))')GCMREG(IC),GCMSTR(IC),
C     +                         DHATDG(IC),D2HDG2(IC),
C     +                         DGDHAT(IC),D2GDH2(IC)
C      ENDDO     

      OPEN(UNIT=9,FILE='meshpoints.res',FORM='FORMATTED')
      DO IC = IGSTAL,IGSTOL
        WRITE(9,'(2(1PE12.4))')GCMREG(IC),GCMSTR(IC)
      ENDDO     
      WRITE(9,*)
      CLOSE(9)

      OPEN(UNIT=9,FILE='dhatdg.res',FORM='FORMATTED')
      DO IC = IGSTAL,IGSTOL
        WRITE(9,'(2(1PE12.4))')GCMREG(IC),DHATDG(IC)
      ENDDO     
      WRITE(9,*)
      CLOSE(9)

      OPEN(UNIT=9,FILE='strpoints.res',FORM='FORMATTED')
      DO IC = IGSTAL,IGSTOL
        WRITE(9,'(2(1PE12.4))')GCMSTR(IC),GCMREG(IC)
      ENDDO     
      WRITE(9,*)
      CLOSE(9)

      OPEN(UNIT=9,FILE='dgdhat.res',FORM='FORMATTED')
      DO IC = IGSTAL+1,IGSTOL
        WRITE(9,'(2(1PE12.4))')GCMREG(IC),DGDHAT(IC)
      ENDDO     
      WRITE(9,*)
      CLOSE(9)

      OPEN(UNIT=9,FILE='d2gdh2.res',FORM='FORMATTED')
      DO IC = IGSTAL+1,IGSTOL
        WRITE(9,'(2(1PE12.4))')GCMREG(IC),D2GDH2(IC)
      ENDDO     
      WRITE(9,*)
      CLOSE(9)

C     =========================================================================


      RETURN
      END
