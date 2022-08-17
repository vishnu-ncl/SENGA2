      SUBROUTINE ZEROYL(FARRAY)
 
C     *************************************************************************
C
C     ZEROYL
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
C     ZEROS THE Y-WISE LEFT END ELEMENTS OF THE ARRAY
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      DOUBLE PRECISION FARRAY(NXSIZE,NYSIZE,NZSIZE)


C     LOCAL DATA
C     ==========
      INTEGER IC,KC


C     BEGIN
C     =====

C     =========================================================================

      DO KC = KSTAL,KSTOL
        DO IC = ISTAL,ISTOL

          FARRAY(IC,JSTAL,KC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
