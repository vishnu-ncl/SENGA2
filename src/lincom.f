      SUBROUTINE LINCOM
 
C     *************************************************************************
C
C     LINCOM
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
C     COMPUTES INTERMEDIATE SOLUTION VALUES IN ERK SCHEME
C     BY DOING LINEAR COMBINATIONS OF LEFT- AND RIGHT-HAND SIDES
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FORNOW
      INTEGER IC,JC,KC,ISPEC


C     BEGIN
C     =====

C     =========================================================================

C     ERK SUBSTEP
C     ===========

C     -------------------------------------------------------------------------
C     NOTE: ALL ERK ERROR ARRAYS ARE INITIALISED TO ZERO IN SUBROUTINE ADAPTT
C     -------------------------------------------------------------------------

C     DENSITY
C     -------
      DO KC = KSTALD,KSTOLD
        DO JC = JSTALD,JSTOLD
          DO IC = ISTALD,ISTOLD

            DERR(IC,JC,KC) = DERR(IC,JC,KC)
     +                              + RKERR(IRKSTP)*DRHS(IC,JC,KC)

            FORNOW = DRUN(IC,JC,KC)
            DRUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*DRHS(IC,JC,KC)
            DRHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*DRHS(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO
                       
C     -------------------------------------------------------------------------

C     U-VELOCITY
C     ----------
      DO KC = KSTALU,KSTOLU
        DO JC = JSTALU,JSTOLU
          DO IC = ISTALU,ISTOLU

            UERR(IC,JC,KC) = UERR(IC,JC,KC)
     +                              + RKERR(IRKSTP)*URHS(IC,JC,KC)

            FORNOW = URUN(IC,JC,KC)
            URUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*URHS(IC,JC,KC)
            URHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*URHS(IC,JC,KC)
            
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
     +                              + RKERR(IRKSTP)*VRHS(IC,JC,KC)

            FORNOW = VRUN(IC,JC,KC)
            VRUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*VRHS(IC,JC,KC)
            VRHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*VRHS(IC,JC,KC)

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
     +                              + RKERR(IRKSTP)*WRHS(IC,JC,KC)

            FORNOW = WRUN(IC,JC,KC)
            WRUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*WRHS(IC,JC,KC)
            WRHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*WRHS(IC,JC,KC)

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
     +                              + RKERR(IRKSTP)*ERHS(IC,JC,KC)

            FORNOW = ERUN(IC,JC,KC)
            ERUN(IC,JC,KC) = FORNOW + RKLHS(IRKSTP)*ERHS(IC,JC,KC)
            ERHS(IC,JC,KC) = FORNOW + RKRHS(IRKSTP)*ERHS(IC,JC,KC)

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
     +                             + RKERR(IRKSTP)*YRHS(IC,JC,KC,ISPEC)

              FORNOW = YRUN(IC,JC,KC,ISPEC)
              YRUN(IC,JC,KC,ISPEC) = FORNOW
     +                             + RKLHS(IRKSTP)*YRHS(IC,JC,KC,ISPEC)
              YRHS(IC,JC,KC,ISPEC) = FORNOW
     +                             + RKRHS(IRKSTP)*YRHS(IC,JC,KC,ISPEC)

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
C            YRHS(IC,JC,KC,NSPEC) = ZERO
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
C              YRHS(IC,JC,KC,NSPEC) = YRHS(IC,JC,KC,NSPEC)
C     +                             + YRHS(IC,JC,KC,ISPEC)
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
C            YRHS(IC,JC,KC,NSPEC)
C     +        = DRHS(IC,JC,KC)*(ONE-YRHS(IC,JC,KC,NSPEC)/DRHS(IC,JC,KC))
C
C          ENDDO
C        ENDDO
C      ENDDO

C     =========================================================================


      RETURN
      END
