      SUBROUTINE FINCOM
 
C     *************************************************************************
C
C     FINCOM
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     15-JAN-2003:  CREATED
C     08-AUG-2012:  RSC EVALUATE ALL SPECIES
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     COMPUTES FINAL SOLUTION VALUES IN ERK SCHEME
C     BY DOING A LINEAR COMBINATION OF LEFT- AND RIGHT-HAND SIDES
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      USE com_senga
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER IC,JC,KC,ISPEC

C     -------------------------------------------------------------------------

C     BEGIN
C     =====

C     =========================================================================

C     FINAL ERK SUBSTEP
C     =================

C     -------------------------------------------------------------------------
C     NOTE: ALL ERK ERROR ARRAYS ARE INITIALISED TO ZERO IN SUBROUTINE ADAPTT
C     -------------------------------------------------------------------------

C     DENSITY
C     ----------
      DO KC = KSTALD,KSTOLD
        DO JC = JSTALD,JSTOLD
          DO IC = ISTALD,ISTOLD

            DERR(IC,JC,KC) = DERR(IC,JC,KC)
     +                     + RKERR(NRKSTP)*DRHS(IC,JC,KC)

            DRUN(IC,JC,KC) = DRUN(IC,JC,KC)
     +                     + RKLHS(NRKSTP)*DRHS(IC,JC,KC)
            DRHS(IC,JC,KC) = DRUN(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     -------------------------------------------------------------------------

C     U VELOCITY
C     ----------
      DO KC = KSTALU,KSTOLU
        DO JC = JSTALU,JSTOLU
          DO IC = ISTALU,ISTOLU

            UERR(IC,JC,KC) = UERR(IC,JC,KC)
     +                     + RKERR(NRKSTP)*URHS(IC,JC,KC)

            URUN(IC,JC,KC) = URUN(IC,JC,KC)
     +                     + RKLHS(NRKSTP)*URHS(IC,JC,KC)
            URHS(IC,JC,KC) = URUN(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     -------------------------------------------------------------------------

C     V-VELOCITY
C     ----------
      DO KC = KSTALV,KSTOLV
        DO JC = JSTALV,JSTOLV
          DO IC = ISTALV,ISTOLV

            VERR(IC,JC,KC) = VERR(IC,JC,KC)
     +                     + RKERR(NRKSTP)*VRHS(IC,JC,KC)

            VRUN(IC,JC,KC) = VRUN(IC,JC,KC)
     +                     + RKLHS(NRKSTP)*VRHS(IC,JC,KC)
            VRHS(IC,JC,KC) = VRUN(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     -------------------------------------------------------------------------

C     W-VELOCITY
C     ----------
      DO KC = KSTALW,KSTOLW
        DO JC = JSTALW,JSTOLW
          DO IC = ISTALW,ISTOLW

            WERR(IC,JC,KC) = WERR(IC,JC,KC)
     +                     + RKERR(NRKSTP)*WRHS(IC,JC,KC)

            WRUN(IC,JC,KC) = WRUN(IC,JC,KC)
     +                     + RKLHS(NRKSTP)*WRHS(IC,JC,KC)
            WRHS(IC,JC,KC) = WRUN(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
      
C     -------------------------------------------------------------------------

C     STAGNATION INTERNAL ENERGY
C     --------------------------
      DO KC = KSTALE,KSTOLE
        DO JC = JSTALE,JSTOLE
          DO IC = ISTALE,ISTOLE

            EERR(IC,JC,KC) = EERR(IC,JC,KC)
     +                     + RKERR(NRKSTP)*ERHS(IC,JC,KC)

            ERUN(IC,JC,KC) = ERUN(IC,JC,KC)
     +                     + RKLHS(NRKSTP)*ERHS(IC,JC,KC)
            ERHS(IC,JC,KC) = ERUN(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
 
C     -------------------------------------------------------------------------

C     SPECIES MASS FRACTIONS
C     ----------------------
C     RSC 08-AUG-2012 EVALUATE ALL SPECIES
C      DO ISPEC = 1,NSPM1
      DO ISPEC = 1,NSPEC

        DO KC = KSTALY,KSTOLY
          DO JC = JSTALY,JSTOLY
            DO IC = ISTALY,ISTOLY

              YERR(IC,JC,KC,ISPEC) = YERR(IC,JC,KC,ISPEC)
     +                             + RKERR(NRKSTP)*YRHS(IC,JC,KC,ISPEC)

              YRUN(IC,JC,KC,ISPEC) = YRUN(IC,JC,KC,ISPEC)
     +                             + RKLHS(NRKSTP)*YRHS(IC,JC,KC,ISPEC)
              YRHS(IC,JC,KC,ISPEC) = YRUN(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     -------------------------------------------------------------------------

CC     NTH SPECIES
C      DO KC = KSTALY,KSTOLY
C        DO JC = JSTALY,JSTOLY
C          DO IC = ISTALY,ISTOLY
C 
C            YRUN(IC,JC,KC,NSPEC) = ZERO
C
C          ENDDO
C        ENDDO
C      ENDDO
C
C      DO ISPEC = 1,NSPM1
C        DO KC = KSTALY,KSTOLY
C          DO JC = JSTALY,JSTOLY
C            DO IC = ISTALY,ISTOLY
C 
C              YRUN(IC,JC,KC,NSPEC) = YRUN(IC,JC,KC,NSPEC)
C     +                             + YRUN(IC,JC,KC,ISPEC)
C
C            ENDDO
C          ENDDO
C        ENDDO
C      ENDDO
C
C      DO KC = KSTALY,KSTOLY
C        DO JC = JSTALY,JSTOLY
C          DO IC = ISTALY,ISTOLY
C 
C            YRUN(IC,JC,KC,NSPEC)
C     +        = DRUN(IC,JC,KC)*(ONE-YRUN(IC,JC,KC,NSPEC)/DRUN(IC,JC,KC))
C
C            YRHS(IC,JC,KC,NSPEC) = YRUN(IC,JC,KC,NSPEC)
C
C          ENDDO
C        ENDDO
C      ENDDO

C     =========================================================================


      RETURN
      END
