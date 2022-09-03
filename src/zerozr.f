      SUBROUTINE ZEROZR(FARRAY)
 
C     *************************************************************************
C
C     ZEROZR
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     06-JUL-2003:  RSC MODIFIED FOR SENGA2
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     ZEROS THE Z-WISE RIGHT END ELEMENTS OF THE ARRAY
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      USE com_senga
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      DOUBLE PRECISION FARRAY(NXSIZE,NYSIZE,NZSIZE)


C     LOCAL DATA
C     ==========
      INTEGER IC,JC


C     BEGIN
C     =====

C     =========================================================================

      DO JC = JSTAL,JSTOL
        DO IC = ISTAL,ISTOL

          FARRAY(IC,JC,KSTOL) = ZERO

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
