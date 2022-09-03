      SUBROUTINE BCDTYL
 
C     *************************************************************************
C
C     BCDTYL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     26-OCT-2013:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR DENSITY
C     AND ITS TIME DERIVATIVE
C
C     Y-DIRECTION LEFT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      USE com_senga
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER IC,KC


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     EVALUATE AND RETURN STRDYL,DDDTYL

      DO KC = KSTAL,KSTOL
        DO IC = ISTAL,ISTOL

C         SET DENSITY TO CONSTANT (INITIAL) VALUE
          STRDYL(IC,KC) = DRIN

C         SET DENSITY TIME DERIVATIVE TO ZERO
          DDDTYL(IC,KC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
