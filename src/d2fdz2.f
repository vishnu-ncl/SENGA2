      SUBROUTINE D2FDZ2(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     D2FDZ2
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     11-APR-2003:  RSC MODIFIED FOR SENGA2
C     10-OCT-2004:  RSC NULL VERSION
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES SECOND Z-DERIVATIVE OF SPECIFIED FUNCTION
C     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
C     EXPLICIT 8TH,6TH,4TH,2ND ORDER END CONDITIONS
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      DOUBLE PRECISION FUNCTN(NXBIGL:NXBIGR,NYBIGL:NYBIGR,NZBIGL:NZBIGR)
      DOUBLE PRECISION FDERIV(NXSIZE,NYSIZE,NZSIZE)


C     LOCAL DATA
C     ==========
      INTEGER IC,JC,KC


C     BEGIN
C     =====

C     =========================================================================

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FDERIV(IC,JC,KC) = ZERO

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
