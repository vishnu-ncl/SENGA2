      SUBROUTINE BOUNDS 
 
C     *************************************************************************
C
C     BOUNDS
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     13-JUL-2003:  RSC MODIFIED FOR SENGA2
C     08-AUG-2012:  RSC EVALUATE ALL SPECIES
C     26-OCT-2013:  RSC ACTIVATE ALL BCS ON ALL SIDES
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     COMPUTES CHARACTERISTIC BOUNDARY CONDITIONS FOR ALL PDES
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
      INTEGER IC,JC,KC
      INTEGER ISPEC


C     BEGIN
C     =====

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================


C     X-DIRECTION LEFT-HAND END
C     -------------------------
      IF(FXLCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUXL = PRIMITIVE U-VELOCITY COMPONENT
C       STRVXL = PRIMITIVE V-VELOCITY COMPONENT
C       STRWXL = PRIMITIVE W-VELOCITY COMPONENT
C       STRPXL = PRESSURE
C       STRDXL = DENSITY
C       STRTXL = TEMPERATURE
C       STREXL = INTERNAL ENERGY
C       STRGXL = MIXTURE CP
C       STRRXL = MIXTURE SPECIFIC GAS CONSTANT
C       STRYXL(ISPEC) = SPECIES MASS FRACTION
C       RATEXL(ISPEC) = SPECIES REACTION RATE
C       STRHXL(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1XL = DUDX
C       BCL2XL = DRHODX
C       BCL3XL = DVDX
C       BCL4XL = DWDX
C       BCL5XL = DPDX
C       BCLYXL(ISPEC) = DYDX

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRHXL(JC,KC,ISPEC) = STRHXL(JC,KC,ISPEC)
     +         - STRGXL(JC,KC)*STRTXL(JC,KC)*RGSPEC(ISPEC)/STRRXL(JC,KC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            GAM1XL(JC,KC) = STRGXL(JC,KC) - STRRXL(JC,KC)
            STREXL(JC,KC) = STREXL(JC,KC) - GAM1XL(JC,KC)*STRTXL(JC,KC)

            GAM1XL(JC,KC) = STRRXL(JC,KC)/GAM1XL(JC,KC)
            OVGMXL(JC,KC) = ONE/GAM1XL(JC,KC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            FORNOW = STRGXL(JC,KC)*GAM1XL(JC,KC)*STRTXL(JC,KC)
            ACOUXL(JC,KC) = SQRT(FORNOW)
            OVA2XL(JC,KC) = ONE/FORNOW

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCXL.EQ.NSBCO1)THEN

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXL(JC,KC) = SORPXL(JC,KC)
     +        + STRHXL(JC,KC,ISPEC)*RATEXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = -SORPXL(JC,KC)*GAM1XL(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L5X AS REQUIRED
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

C             OLD VALUE OF L5X
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +        *(BCL5XL(JC,KC)+STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC))

C             SUBTRACT FROM NEW VALUE OF L5X
              BCL5XL(JC,KC)= HALF*SORPXL(JC,KC)
     +                     + COBCXL*ACOUXL(JC,KC)*(STRPXL(JC,KC)-PINFXL)
     +                     - BCL5XL(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

              URHS(ISTAL,JC,KC) = URHS(ISTAL,JC,KC)
     +   - BCL5XL(JC,KC)*OVA2XL(JC,KC)*(STRUXL(JC,KC)+ACOUXL(JC,KC))

              VRHS(ISTAL,JC,KC) = VRHS(ISTAL,JC,KC)
     +                   - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRVXL(JC,KC)

              WRHS(ISTAL,JC,KC) = WRHS(ISTAL,JC,KC)
     +                   - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRWXL(JC,KC)

              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +               - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
     +                              + STRUXL(JC,KC)/ACOUXL(JC,KC)
     +                              + OVGMXL(JC,KC))

            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     +                 - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRYXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCXL.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXL(JC,KC) = SORPXL(JC,KC)
     +        + STRHXL(JC,KC,ISPEC)*RATEXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = -SORPXL(JC,KC)*GAM1XL(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L2X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                      *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
C             L1X UNCHANGED
              BCL2XL(JC,KC) = -BCL2XL(JC,KC)
              BCL3XL(JC,KC) = -BCL3XL(JC,KC)
              BCL4XL(JC,KC) = -BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*SORPXL(JC,KC)
     +                     + COBCXL*ACOUXL(JC,KC)*(STRPXL(JC,KC)-PINFXL)
     +                      - BCL5XL(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF L's
                BCLYXL(JC,KC,ISPEC) = STRUXL(JC,KC)*BCLYXL(JC,KC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
                BCLYXL(JC,KC,ISPEC) = RATEXL(JC,KC,ISPEC)/STRDXL(JC,KC)
     +                              - BCLYXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

              URHS(ISTAL,JC,KC) = URHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STRUXL(JC,KC)
     +      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*(STRUXL(JC,KC)+ACOUXL(JC,KC))

              VRHS(ISTAL,JC,KC) = VRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STRVXL(JC,KC)
     +                          - BCL3XL(JC,KC)*STRDXL(JC,KC)
     +                      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRVXL(JC,KC)

              WRHS(ISTAL,JC,KC) = WRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STRWXL(JC,KC)
     +                          - BCL4XL(JC,KC)*STRDXL(JC,KC)
     +                      - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRWXL(JC,KC)

              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STREXL(JC,KC)
     +                      - BCL3XL(JC,KC)*STRDXL(JC,KC)*STRVXL(JC,KC)
     +                      - BCL4XL(JC,KC)*STRDXL(JC,KC)*STRWXL(JC,KC)
     +                      - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
     +                                     + STRUXL(JC,KC)/ACOUXL(JC,KC)
     +                                     + OVGMXL(JC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = BCLYXL(JC,KC,ISPEC)*STRDXL(JC,KC)

                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                            - FORNOW*STRHXL(JC,KC,ISPEC)

                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     + - (BCL2XL(JC,KC)+BCL5XL(JC,KC)*OVA2XL(JC,KC))*STRYXL(JC,KC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXL.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SYDTXL(JC,KC) = ZERO
              SORPXL(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SYDTXL(JC,KC) = SYDTXL(JC,KC)
     +                        + DYDTXL(JC,KC,ISPEC)*RGSPEC(ISPEC)
                SORPXL(JC,KC) = SORPXL(JC,KC)
     +                        + STRHXL(JC,KC,ISPEC)*RATEXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SYDTXL(JC,KC) = SYDTXL(JC,KC)/STRRXL(JC,KC)
              SORPXL(JC,KC) = -SORPXL(JC,KC)*GAM1XL(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X,L2X,L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL1XL(JC,KC) = HALF*(STRUXL(JC,KC)-ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)-FORNOW)
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                      *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1X UNCHANGED
              BCL5XL(JC,KC) = BCL1XL(JC,KC)
     +                      - STRDXL(JC,KC)*ACOUXL(JC,KC)*DUDTXL(JC,KC)
     +                      - BCL5XL(JC,KC)
              BCL2XL(JC,KC) = GAM1XL(JC,KC)*OVA2XL(JC,KC)
     +                       *(BCL1XL(JC,KC)+BCL5XL(JC,KC))
     +                      + STRDXL(JC,KC)*(DTDTXL(JC,KC)/STRTXL(JC,KC)
     +                      - SORPXL(JC,KC)/STRPXL(JC,KC)
     +                      + SYDTXL(JC,KC))
     +                      - BCL2XL(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCXL.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL1XL(JC,KC) = HALF*(STRUXL(JC,KC)-ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)-FORNOW)
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                    *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1X UNCHANGED
              FORNOW = BCL1XL(JC,KC)
     +               - STRDXL(JC,KC)*ACOUXL(JC,KC)*DUDTXL(JC,KC)
              BCL2XL(JC,KC) = -DDDTXL(JC,KC)
     +                      - OVA2XL(JC,KC)*(BCL1XL(JC,KC)+FORNOW)
     +                      - BCL2XL(JC,KC)
              BCL3XL(JC,KC) = -DVDTXL(JC,KC) - BCL3XL(JC,KC)
              BCL4XL(JC,KC) = -DWDTXL(JC,KC) - BCL4XL(JC,KC)
              BCL5XL(JC,KC) = FORNOW - BCL5XL(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)*STREXL(JC,KC)
     +                  - BCL3XL(JC,KC)*STRDXL(JC,KC)*STRVXL(JC,KC)
     +                  - BCL4XL(JC,KC)*STRDXL(JC,KC)*STRWXL(JC,KC)
     +                  - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
     +                                 + STRUXL(JC,KC)/ACOUXL(JC,KC)
     +                                 + OVGMXL(JC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                BCLYXL(JC,KC,ISPEC) = RATEXL(JC,KC,ISPEC)/STRDXL(JC,KC)
     +                              - DYDTXL(JC,KC,ISPEC)
     +                              - STRUXL(JC,KC)*BCLYXL(JC,KC,ISPEC)

                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +          - BCLYXL(JC,KC,ISPEC)*STRDXL(JC,KC)*STRHXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

        IF(NSBCXL.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 1
C         NO-SLIP WALL - ADIABATIC

C         ALL VELOCITY COMPONENTS IMPOSED
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1X,L3X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL1XL(JC,KC) = HALF*(STRUXL(JC,KC)-ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)-FORNOW)
              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1X,L2X UNCHANGED
              BCL3XL(JC,KC) = -DVDTXL(JC,KC) - BCL3XL(JC,KC)
              BCL4XL(JC,KC) = -DWDTXL(JC,KC) - BCL4XL(JC,KC)
              BCL5XL(JC,KC) = BCL1XL(JC,KC)
     +                      - STRDXL(JC,KC)*ACOUXL(JC,KC)*DUDTXL(JC,KC)
     +                      - BCL5XL(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

              ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +                  - BCL3XL(JC,KC)*STRDXL(JC,KC)*STRVXL(JC,KC)
     +                  - BCL4XL(JC,KC)*STRDXL(JC,KC)*STRWXL(JC,KC)
     +                  - BCL5XL(JC,KC)*(OVA2XL(JC,KC)*STREXL(JC,KC)
     +                                 + STRUXL(JC,KC)/ACOUXL(JC,KC)
     +                                 + OVGMXL(JC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     +                 - BCL5XL(JC,KC)*OVA2XL(JC,KC)*STRYXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCXL.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXL(JC,KC) = SORPXL(JC,KC)
     +                        + STRHXL(JC,KC,ISPEC)*RATEXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXL(JC,KC) = -SORPXL(JC,KC)*GAM1XL(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXL(JC,KC)*ACOUXL(JC,KC)*BCL1XL(JC,KC)
              BCL1XL(JC,KC) = HALF*(STRUXL(JC,KC)-ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)-FORNOW)
              BCL2XL(JC,KC) = STRUXL(JC,KC)
     +                    *(BCL2XL(JC,KC)-BCL5XL(JC,KC)*OVA2XL(JC,KC))
              BCL3XL(JC,KC) = STRUXL(JC,KC)*BCL3XL(JC,KC)
              BCL4XL(JC,KC) = STRUXL(JC,KC)*BCL4XL(JC,KC)
              BCL5XL(JC,KC) = HALF*(STRUXL(JC,KC)+ACOUXL(JC,KC))
     +                      *(BCL5XL(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1X UNCHANGED
              BCL3XL(JC,KC) = -DVDTXL(JC,KC) - BCL3XL(JC,KC)
              BCL4XL(JC,KC) = -DWDTXL(JC,KC) - BCL4XL(JC,KC)
              BCL5XL(JC,KC) = BCL1XL(JC,KC)
     +                      - STRDXL(JC,KC)*ACOUXL(JC,KC)*DUDTXL(JC,KC)
     +                      - BCL5XL(JC,KC)
              BCL2XL(JC,KC) = GAM1XL(JC,KC)*OVA2XL(JC,KC)
     +                       *(BCL1XL(JC,KC)+BCL5XL(JC,KC))
     +                      + STRDXL(JC,KC)*(DTDTXL(JC,KC)/STRTXL(JC,KC)
     +                      - SORPXL(JC,KC)/STRPXL(JC,KC))
     +                      - BCL2XL(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF LYX
                BCLYXL(JC,KC,ISPEC) = STRUXL(JC,KC)*BCLYXL(JC,KC,ISPEC)

C               UPDATE L2X
                BCL2XL(JC,KC) = BCL2XL(JC,KC)
     +                        + (RATEXL(JC,KC,ISPEC)
     +                         - STRDXL(JC,KC)*BCLYXL(JC,KC,ISPEC))
     +                          *RGSPEC(ISPEC)/STRRXL(JC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)
     +                          - BCL2XL(JC,KC)
     +                          - BCL5XL(JC,KC)*OVA2XL(JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTAL,JC,KC,ISPEC) = YRHS(ISTAL,JC,KC,ISPEC)
     + - (BCL2XL(JC,KC)+BCL5XL(JC,KC)*OVA2XL(JC,KC))*STRYXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF  
C     X-DIRECTION LEFT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     X-DIRECTION RIGHT-HAND END
C     --------------------------
      IF(FXRCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUXR = PRIMITIVE U-VELOCITY COMPONENT
C       STRVXR = PRIMITIVE V-VELOCITY COMPONENT
C       STRWXR = PRIMITIVE W-VELOCITY COMPONENT
C       STRPXR = PRESSURE
C       STRDXR = DENSITY
C       STRTXR = TEMPERATURE
C       STREXR = INTERNAL ENERGY
C       STRGXR = MIXTURE CP
C       STRRXR = MIXTURE SPECIFIC GAS CONSTANT
C       STRYXR(ISPEC) = SPECIES MASS FRACTION
C       RATEXR(ISPEC) = SPECIES REACTION RATE
C       STRHXR(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1XR = DUDX
C       BCL2XR = DRHODX
C       BCL3XR = DVDX
C       BCL4XR = DWDX
C       BCL5XR = DPDX
C       BCLYXR(ISPEC) = DYDX

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              STRHXR(JC,KC,ISPEC) = STRHXR(JC,KC,ISPEC)
     +         - STRGXR(JC,KC)*STRTXR(JC,KC)*RGSPEC(ISPEC)/STRRXR(JC,KC)

            ENDDO

          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            GAM1XR(JC,KC) = STRGXR(JC,KC) - STRRXR(JC,KC)
            STREXR(JC,KC) = STREXR(JC,KC) - GAM1XR(JC,KC)*STRTXR(JC,KC)

            GAM1XR(JC,KC) = STRRXR(JC,KC)/GAM1XR(JC,KC)
            OVGMXR(JC,KC) = ONE/GAM1XR(JC,KC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            FORNOW = STRGXR(JC,KC)*GAM1XR(JC,KC)*STRTXR(JC,KC)
            ACOUXR(JC,KC) = SQRT(FORNOW)
            OVA2XR(JC,KC) = ONE/FORNOW

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCXR.EQ.NSBCO1)THEN

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXR(JC,KC) = SORPXR(JC,KC)
     +        + STRHXR(JC,KC,ISPEC)*RATEXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = -SORPXR(JC,KC)*GAM1XR(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L1X AS REQUIRED
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

C             OLD VALUE OF L1X
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +        *(BCL5XR(JC,KC)-STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC))

C             SUBTRACT FROM NEW VALUE OF L1X
              BCL1XR(JC,KC)= HALF*SORPXR(JC,KC)
     +                     + COBCXR*ACOUXR(JC,KC)*(STRPXR(JC,KC)-PINFXR)
     +                     - BCL1XR(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)

              URHS(ISTOL,JC,KC) = URHS(ISTOL,JC,KC)
     +   - BCL1XR(JC,KC)*OVA2XR(JC,KC)*(STRUXR(JC,KC)-ACOUXR(JC,KC))

              VRHS(ISTOL,JC,KC) = VRHS(ISTOL,JC,KC)
     +                   - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRVXR(JC,KC)

              WRHS(ISTOL,JC,KC) = WRHS(ISTOL,JC,KC)
     +                   - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRWXR(JC,KC)

              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +               - BCL1XR(JC,KC)*(OVA2XR(JC,KC)*STREXR(JC,KC)
     +                              - STRUXR(JC,KC)/ACOUXR(JC,KC)
     +                              + OVGMXR(JC,KC))

            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
     +                 - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRYXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCXR.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXR(JC,KC) = SORPXR(JC,KC)
     +        + STRHXR(JC,KC,ISPEC)*RATEXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = -SORPXR(JC,KC)*GAM1XR(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X-L4X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL2XR(JC,KC) = STRUXR(JC,KC)
     +                      *(BCL2XR(JC,KC)-BCL5XR(JC,KC)*OVA2XR(JC,KC))
              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
C             L5X UNCHANGED
              BCL1XR(JC,KC) = HALF*SORPXR(JC,KC)
     +                     + COBCXR*ACOUXR(JC,KC)*(STRPXR(JC,KC)-PINFXR)
     +                      - BCL1XR(JC,KC)
              BCL2XR(JC,KC) = -BCL2XR(JC,KC)
              BCL3XR(JC,KC) = -BCL3XR(JC,KC)
              BCL4XR(JC,KC) = -BCL4XR(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF L's
                BCLYXR(JC,KC,ISPEC) = STRUXR(JC,KC)*BCLYXR(JC,KC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
                BCLYXR(JC,KC,ISPEC) = RATEXR(JC,KC,ISPEC)/STRDXR(JC,KC)
     +                              - BCLYXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)
     +                          - BCL2XR(JC,KC)

              URHS(ISTOL,JC,KC) = URHS(ISTOL,JC,KC)
     +      - BCL1XR(JC,KC)*OVA2XR(JC,KC)*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                          - BCL2XR(JC,KC)*STRUXR(JC,KC)

              VRHS(ISTOL,JC,KC) = VRHS(ISTOL,JC,KC)
     +                      - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRVXR(JC,KC)
     +                          - BCL2XR(JC,KC)*STRVXR(JC,KC)
     +                          - BCL3XR(JC,KC)*STRDXR(JC,KC)

              WRHS(ISTOL,JC,KC) = WRHS(ISTOL,JC,KC)
     +                      - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRWXR(JC,KC)
     +                          - BCL2XR(JC,KC)*STRWXR(JC,KC)
     +                          - BCL4XR(JC,KC)*STRDXR(JC,KC)

              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                      - BCL1XR(JC,KC)*(OVA2XR(JC,KC)*STREXR(JC,KC)
     +                                     - STRUXR(JC,KC)/ACOUXR(JC,KC)
     +                                     + OVGMXR(JC,KC))
     +                      - BCL2XR(JC,KC)*STREXR(JC,KC)
     +                      - BCL3XR(JC,KC)*STRDXR(JC,KC)*STRVXR(JC,KC)
     +                      - BCL4XR(JC,KC)*STRDXR(JC,KC)*STRWXR(JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                FORNOW = BCLYXR(JC,KC,ISPEC)*STRDXR(JC,KC)

                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                            - FORNOW*STRHXR(JC,KC,ISPEC)

                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
     + - (BCL2XR(JC,KC)+BCL1XR(JC,KC)*OVA2XR(JC,KC))*STRYXR(JC,KC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXR.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SYDTXR(JC,KC) = ZERO
              SORPXR(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SYDTXR(JC,KC) = SYDTXR(JC,KC)
     +                        + DYDTXR(JC,KC,ISPEC)*RGSPEC(ISPEC)
                SORPXR(JC,KC) = SORPXR(JC,KC)
     +                        + STRHXR(JC,KC,ISPEC)*RATEXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SYDTXR(JC,KC) = SYDTXR(JC,KC)/STRRXR(JC,KC)
              SORPXR(JC,KC) = -SORPXR(JC,KC)*GAM1XR(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X,L2X,L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL2XR(JC,KC) = STRUXR(JC,KC)
     +                      *(BCL2XR(JC,KC)-BCL5XR(JC,KC)*OVA2XR(JC,KC))
              BCL5XR(JC,KC) = HALF*(STRUXR(JC,KC)+ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5X UNCHANGED
              BCL1XR(JC,KC) = BCL5XR(JC,KC)
     +                      + STRDXR(JC,KC)*ACOUXR(JC,KC)*DUDTXR(JC,KC)
     +                      - BCL1XR(JC,KC)
              BCL2XR(JC,KC) = GAM1XR(JC,KC)*OVA2XR(JC,KC)
     +                       *(BCL1XR(JC,KC)+BCL5XR(JC,KC))
     +                      + STRDXR(JC,KC)*(DTDTXR(JC,KC)/STRTXR(JC,KC)
     +                      - SORPXR(JC,KC)/STRPXR(JC,KC)
     +                      + SYDTXR(JC,KC))
     +                      - BCL2XR(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL2XR(JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCXR.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL2XR(JC,KC) = STRUXR(JC,KC)
     +                    *(BCL2XR(JC,KC)-BCL5XR(JC,KC)*OVA2XR(JC,KC))
              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)
              BCL5XR(JC,KC) = HALF*(STRUXR(JC,KC)+ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5X UNCHANGED
              FORNOW = BCL5XR(JC,KC)
     +               + STRDXR(JC,KC)*ACOUXR(JC,KC)*DUDTXR(JC,KC)
              BCL1XR(JC,KC) = FORNOW - BCL1XR(JC,KC)
              BCL2XR(JC,KC) = -DDDTXR(JC,KC)
     +                      - OVA2XR(JC,KC)*(BCL5XR(JC,KC)+FORNOW)
     +                      - BCL2XR(JC,KC)
              BCL3XR(JC,KC) = -DVDTXR(JC,KC) - BCL3XR(JC,KC)
              BCL4XR(JC,KC) = -DWDTXR(JC,KC) - BCL4XR(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                  - BCL1XR(JC,KC)*(OVA2XR(JC,KC)*STREXR(JC,KC)
     +                                 - STRUXR(JC,KC)/ACOUXR(JC,KC)
     +                                 + OVGMXR(JC,KC))
     +                  - BCL2XR(JC,KC)*STREXR(JC,KC)
     +                  - BCL3XR(JC,KC)*STRDXR(JC,KC)*STRVXR(JC,KC)
     +                  - BCL4XR(JC,KC)*STRDXR(JC,KC)*STRWXR(JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                BCLYXR(JC,KC,ISPEC) = RATEXR(JC,KC,ISPEC)/STRDXR(JC,KC)
     +                              - DYDTXR(JC,KC,ISPEC)
     +                              - STRUXR(JC,KC)*BCLYXR(JC,KC,ISPEC)

                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +          - BCLYXR(JC,KC,ISPEC)*STRDXR(JC,KC)*STRHXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

        IF(NSBCXR.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 1
C         NO-SLIP WALL - ADIABATIC

C         ALL VELOCITY COMPONENTS IMPOSED
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1X,L3X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)
              BCL5XR(JC,KC) = HALF*(STRUXR(JC,KC)+ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L2X,L5X UNCHANGED
              BCL1XR(JC,KC) = BCL5XR(JC,KC)
     +                      + STRDXR(JC,KC)*ACOUXR(JC,KC)*DUDTXR(JC,KC)
     +                      - BCL1XR(JC,KC)
              BCL3XR(JC,KC) = -DVDTXR(JC,KC) - BCL3XR(JC,KC)
              BCL4XR(JC,KC) = -DWDTXR(JC,KC) - BCL4XR(JC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)

              ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +                  - BCL1XR(JC,KC)*(OVA2XR(JC,KC)*STREXR(JC,KC)
     +                                 + STRUXR(JC,KC)/ACOUXR(JC,KC)
     +                                 + OVGMXR(JC,KC))
     +                  - BCL3XR(JC,KC)*STRDXR(JC,KC)*STRVXR(JC,KC)
     +                  - BCL4XR(JC,KC)*STRDXR(JC,KC)*STRWXR(JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
     +                 - BCL1XR(JC,KC)*OVA2XR(JC,KC)*STRYXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCXR.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                SORPXR(JC,KC) = SORPXR(JC,KC)
     +                        + STRHXR(JC,KC,ISPEC)*RATEXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              SORPXR(JC,KC) = -SORPXR(JC,KC)*GAM1XR(JC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1X-L5X
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDXR(JC,KC)*ACOUXR(JC,KC)*BCL1XR(JC,KC)
              BCL1XR(JC,KC) = HALF*(STRUXR(JC,KC)-ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)-FORNOW)
              BCL2XR(JC,KC) = STRUXR(JC,KC)
     +                    *(BCL2XR(JC,KC)-BCL5XR(JC,KC)*OVA2XR(JC,KC))
              BCL3XR(JC,KC) = STRUXR(JC,KC)*BCL3XR(JC,KC)
              BCL4XR(JC,KC) = STRUXR(JC,KC)*BCL4XR(JC,KC)
              BCL5XR(JC,KC) = HALF*(STRUXR(JC,KC)+ACOUXR(JC,KC))
     +                      *(BCL5XR(JC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5X UNCHANGED
              BCL1XR(JC,KC) = BCL5XR(JC,KC)
     +                      + STRDXR(JC,KC)*ACOUXR(JC,KC)*DUDTXR(JC,KC)
     +                      - BCL1XR(JC,KC)
              BCL3XR(JC,KC) = -DVDTXR(JC,KC) - BCL3XR(JC,KC)
              BCL4XR(JC,KC) = -DWDTXR(JC,KC) - BCL4XR(JC,KC)
              BCL2XR(JC,KC) = GAM1XR(JC,KC)*OVA2XR(JC,KC)
     +                       *(BCL1XR(JC,KC)+BCL5XR(JC,KC))
     +                      + STRDXR(JC,KC)*(DTDTXR(JC,KC)/STRTXR(JC,KC)
     +                      - SORPXR(JC,KC)/STRPXR(JC,KC))
     +                      - BCL2XR(JC,KC)

            ENDDO
          ENDDO

C         LYX
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

C               OLD VALUE OF LYX
                BCLYXR(JC,KC,ISPEC) = STRUXR(JC,KC)*BCLYXR(JC,KC,ISPEC)

C               UPDATE L2X
                BCL2XR(JC,KC) = BCL2XR(JC,KC)
     +                        + (RATEXR(JC,KC,ISPEC)
     +                         - STRDXR(JC,KC)*BCLYXR(JC,KC,ISPEC))
     +                          *RGSPEC(ISPEC)/STRRXR(JC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)
     +                          - BCL1XR(JC,KC)*OVA2XR(JC,KC)
     +                          - BCL2XR(JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTOL,JC,KC,ISPEC) = YRHS(ISTOL,JC,KC,ISPEC)
     + - (BCL2XR(JC,KC)+BCL1XR(JC,KC)*OVA2XR(JC,KC))*STRYXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     X-DIRECTION RIGHT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     Y-DIRECTION LEFT-HAND END
C     -------------------------
      IF(FYLCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUYL = PRIMITIVE U-VELOCITY COMPONENT
C       STRVYL = PRIMITIVE V-VELOCITY COMPONENT
C       STRWYL = PRIMITIVE W-VELOCITY COMPONENT
C       STRPYL = PRESSURE
C       STRDYL = DENSITY
C       STRTYL = TEMPERATURE
C       STREYL = INTERNAL ENERGY
C       STRGYL = MIXTURE CP
C       STRRYL = MIXTURE SPECIFIC GAS CONSTANT
C       STRYYL(ISPEC) = SPECIES MASS FRACTION
C       RATEYL(ISPEC) = SPECIES REACTION RATE
C       STRHYL(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1YL = DVDY
C       BCL2YL = DRHODY
C       BCL3YL = DUDY
C       BCL4YL = DWDY
C       BCL5YL = DPDY
C       BCLYYL(ISPEC) = DYDY

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRHYL(IC,KC,ISPEC) = STRHYL(IC,KC,ISPEC)
     +         - STRGYL(IC,KC)*STRTYL(IC,KC)*RGSPEC(ISPEC)/STRRYL(IC,KC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            GAM1YL(IC,KC) = STRGYL(IC,KC) - STRRYL(IC,KC)
            STREYL(IC,KC) = STREYL(IC,KC) - GAM1YL(IC,KC)*STRTYL(IC,KC)

            GAM1YL(IC,KC) = STRRYL(IC,KC)/GAM1YL(IC,KC)
            OVGMYL(IC,KC) = ONE/GAM1YL(IC,KC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = STRGYL(IC,KC)*GAM1YL(IC,KC)*STRTYL(IC,KC)
            ACOUYL(IC,KC) = SQRT(FORNOW)
            OVA2YL(IC,KC) = ONE/FORNOW

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCYL.EQ.NSBCO1)THEN

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1, NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYL(IC,KC) = SORPYL(IC,KC)
     +        + STRHYL(IC,KC,ISPEC)*RATEYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = -SORPYL(IC,KC)*GAM1YL(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L5Y AS REQUIRED
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L5Y
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +        *(BCL5YL(IC,KC)+STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC))

C             SUBTRACT FROM NEW VALUE OF L5Y
              BCL5YL(IC,KC)= HALF*SORPYL(IC,KC)
     +                     + COBCYL*ACOUYL(IC,KC)*(STRPYL(IC,KC)-PINFYL)
     +                     - BCL5YL(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)

              URHS(IC,JSTAL,KC) = URHS(IC,JSTAL,KC)
     +                   - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRUYL(IC,KC)

              VRHS(IC,JSTAL,KC) = VRHS(IC,JSTAL,KC)
     +   - BCL5YL(IC,KC)*OVA2YL(IC,KC)*(STRVYL(IC,KC)+ACOUYL(IC,KC))

              WRHS(IC,JSTAL,KC) = WRHS(IC,JSTAL,KC)
     +                   - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRWYL(IC,KC)

              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +               - BCL5YL(IC,KC)*(OVA2YL(IC,KC)*STREYL(IC,KC)
     +                              + STRVYL(IC,KC)/ACOUYL(IC,KC)
     +                              + OVGMYL(IC,KC))

            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
     +                 - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRYYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCYL.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYL(IC,KC) = SORPYL(IC,KC)
     +        + STRHYL(IC,KC,ISPEC)*RATEYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = -SORPYL(IC,KC)*GAM1YL(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L2Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                      *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
C             L1Y UNCHANGED
              BCL2YL(IC,KC) = -BCL2YL(IC,KC)
              BCL3YL(IC,KC) = -BCL3YL(IC,KC)
              BCL4YL(IC,KC) = -BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*SORPYL(IC,KC)
     +                     + COBCYL*ACOUYL(IC,KC)*(STRPYL(IC,KC)-PINFYL)
     +                      - BCL5YL(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF L's
                BCLYYL(IC,KC,ISPEC) = STRVYL(IC,KC)*BCLYYL(IC,KC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
                BCLYYL(IC,KC,ISPEC) = RATEYL(IC,KC,ISPEC)/STRDYL(IC,KC)
     +                              - BCLYYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)

              URHS(IC,JSTAL,KC) = URHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STRUYL(IC,KC)
     +                          - BCL3YL(IC,KC)*STRDYL(IC,KC)
     +                      - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRUYL(IC,KC)

              VRHS(IC,JSTAL,KC) = VRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STRVYL(IC,KC)
     +      - BCL5YL(IC,KC)*OVA2YL(IC,KC)*(STRVYL(IC,KC)+ACOUYL(IC,KC))

              WRHS(IC,JSTAL,KC) = WRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STRWYL(IC,KC)
     +                          - BCL4YL(IC,KC)*STRDYL(IC,KC)
     +                      - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRWYL(IC,KC)

              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STREYL(IC,KC)
     +                      - BCL3YL(IC,KC)*STRDYL(IC,KC)*STRUYL(IC,KC)
     +                      - BCL4YL(IC,KC)*STRDYL(IC,KC)*STRWYL(IC,KC)
     +                      - BCL5YL(IC,KC)*(OVA2YL(IC,KC)*STREYL(IC,KC)
     +                                     + STRVYL(IC,KC)/ACOUYL(IC,KC)
     +                                     + OVGMYL(IC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = BCLYYL(IC,KC,ISPEC)*STRDYL(IC,KC)

                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                            - FORNOW*STRHYL(IC,KC,ISPEC)

                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
     + - (BCL2YL(IC,KC)+BCL5YL(IC,KC)*OVA2YL(IC,KC))*STRYYL(IC,KC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYL.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SYDTYL(IC,KC) = ZERO
              SORPYL(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SYDTYL(IC,KC) = SYDTYL(IC,KC)
     +                        + DYDTYL(IC,KC,ISPEC)*RGSPEC(ISPEC)
                SORPYL(IC,KC) = SORPYL(IC,KC)
     +                        + STRHYL(IC,KC,ISPEC)*RATEYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SYDTYL(IC,KC) = SYDTYL(IC,KC)/STRRYL(IC,KC)
              SORPYL(IC,KC) = -SORPYL(IC,KC)*GAM1YL(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y,L2Y,L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL1YL(IC,KC) = HALF*(STRVYL(IC,KC)-ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)-FORNOW)
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                      *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              BCL5YL(IC,KC) = BCL1YL(IC,KC)
     +                      - STRDYL(IC,KC)*ACOUYL(IC,KC)*DVDTYL(IC,KC)
     +                      - BCL5YL(IC,KC)
              BCL2YL(IC,KC) = GAM1YL(IC,KC)*OVA2YL(IC,KC)
     +                       *(BCL1YL(IC,KC)+BCL5YL(IC,KC))
     +                      + STRDYL(IC,KC)*(DTDTYL(IC,KC)/STRTYL(IC,KC)
     +                      - SORPYL(IC,KC)/STRPYL(IC,KC)
     +                      + SYDTYL(IC,KC))
     +                      - BCL2YL(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCYL.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL1YL(IC,KC) = HALF*(STRVYL(IC,KC)-ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)-FORNOW)
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                    *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              FORNOW = BCL1YL(IC,KC)
     +               - STRDYL(IC,KC)*ACOUYL(IC,KC)*DVDTYL(IC,KC)
              BCL2YL(IC,KC) = -DDDTYL(IC,KC)
     +                      - OVA2YL(IC,KC)*(BCL1YL(IC,KC)+FORNOW)
     +                      - BCL2YL(IC,KC)
              BCL3YL(IC,KC) = -DUDTYL(IC,KC) - BCL3YL(IC,KC)
              BCL4YL(IC,KC) = -DWDTYL(IC,KC) - BCL4YL(IC,KC)
              BCL5YL(IC,KC) = FORNOW - BCL5YL(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)*STREYL(IC,KC)
     +                  - BCL3YL(IC,KC)*STRDYL(IC,KC)*STRUYL(IC,KC)
     +                  - BCL4YL(IC,KC)*STRDYL(IC,KC)*STRWYL(IC,KC)
     +                  - BCL5YL(IC,KC)*(OVA2YL(IC,KC)*STREYL(IC,KC)
     +                                 + STRVYL(IC,KC)/ACOUYL(IC,KC)
     +                                 + OVGMYL(IC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                BCLYYL(IC,KC,ISPEC) = RATEYL(IC,KC,ISPEC)/STRDYL(IC,KC)
     +                              - DYDTYL(IC,KC,ISPEC)
     +                              - STRVYL(IC,KC)*BCLYYL(IC,KC,ISPEC)

                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +          - BCLYYL(IC,KC,ISPEC)*STRDYL(IC,KC)*STRHYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

        IF(NSBCYL.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 1
C         NO-SLIP WALL - ADIABATIC

C         ALL VELOCITY COMPONENTS IMPOSED
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Y,L3Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL1YL(IC,KC) = HALF*(STRVYL(IC,KC)-ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)-FORNOW)
              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y,L2Y UNCHANGED
              BCL3YL(IC,KC) = -DUDTYL(IC,KC) - BCL3YL(IC,KC)
              BCL4YL(IC,KC) = -DWDTYL(IC,KC) - BCL4YL(IC,KC)
              BCL5YL(IC,KC) = BCL1YL(IC,KC)
     +                      - STRDYL(IC,KC)*ACOUYL(IC,KC)*DVDTYL(IC,KC)
     +                      - BCL5YL(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)

              ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +                  - BCL3YL(IC,KC)*STRDYL(IC,KC)*STRUYL(IC,KC)
     +                  - BCL4YL(IC,KC)*STRDYL(IC,KC)*STRWYL(IC,KC)
     +                  - BCL5YL(IC,KC)*(OVA2YL(IC,KC)*STREYL(IC,KC)
     +                                 + STRVYL(IC,KC)/ACOUYL(IC,KC)
     +                                 + OVGMYL(IC,KC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
     +                 - BCL5YL(IC,KC)*OVA2YL(IC,KC)*STRYYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  


C       =======================================================================

        IF(NSBCYL.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYL(IC,KC) = SORPYL(IC,KC)
     +                        + STRHYL(IC,KC,ISPEC)*RATEYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYL(IC,KC) = -SORPYL(IC,KC)*GAM1YL(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYL(IC,KC)*ACOUYL(IC,KC)*BCL1YL(IC,KC)
              BCL1YL(IC,KC) = HALF*(STRVYL(IC,KC)-ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)-FORNOW)
              BCL2YL(IC,KC) = STRVYL(IC,KC)
     +                    *(BCL2YL(IC,KC)-BCL5YL(IC,KC)*OVA2YL(IC,KC))
              BCL3YL(IC,KC) = STRVYL(IC,KC)*BCL3YL(IC,KC)
              BCL4YL(IC,KC) = STRVYL(IC,KC)*BCL4YL(IC,KC)
              BCL5YL(IC,KC) = HALF*(STRVYL(IC,KC)+ACOUYL(IC,KC))
     +                      *(BCL5YL(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              BCL3YL(IC,KC) = -DUDTYL(IC,KC) - BCL3YL(IC,KC)
              BCL4YL(IC,KC) = -DWDTYL(IC,KC) - BCL4YL(IC,KC)
              BCL5YL(IC,KC) = BCL1YL(IC,KC)
     +                      - STRDYL(IC,KC)*ACOUYL(IC,KC)*DVDTYL(IC,KC)
     +                      - BCL5YL(IC,KC)
              BCL2YL(IC,KC) = GAM1YL(IC,KC)*OVA2YL(IC,KC)
     +                       *(BCL1YL(IC,KC)+BCL5YL(IC,KC))
     +                      + STRDYL(IC,KC)*(DTDTYL(IC,KC)/STRTYL(IC,KC)
     +                      - SORPYL(IC,KC)/STRPYL(IC,KC))
     +                      - BCL2YL(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYY
                BCLYYL(IC,KC,ISPEC) = STRVYL(IC,KC)*BCLYYL(IC,KC,ISPEC)

C               UPDATE L2Y
                BCL2YL(IC,KC) = BCL2YL(IC,KC)
     +                        + (RATEYL(IC,KC,ISPEC)
     +                         - STRDYL(IC,KC)*BCLYYL(IC,KC,ISPEC))
     +                          *RGSPEC(ISPEC)/STRRYL(IC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)
     +                          - BCL2YL(IC,KC)
     +                          - BCL5YL(IC,KC)*OVA2YL(IC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTAL,KC,ISPEC) = YRHS(IC,JSTAL,KC,ISPEC)
     + - (BCL2YL(IC,KC)+BCL5YL(IC,KC)*OVA2YL(IC,KC))*STRYYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     Y-DIRECTION LEFT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     Y-DIRECTION RIGHT-HAND END
C     --------------------------
      IF(FYRCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUYR = PRIMITIVE U-VELOCITY COMPONENT
C       STRVYR = PRIMITIVE V-VELOCITY COMPONENT
C       STRWYR = PRIMITIVE W-VELOCITY COMPONENT
C       STRPYR = PRESSURE
C       STRDYR = DENSITY
C       STRTYR = TEMPERATURE
C       STREYR = INTERNAL ENERGY
C       STRGYR = MIXTURE CP
C       STRRYR = MIXTURE SPECIFIC GAS CONSTANT
C       STRYYR(ISPEC) = SPECIES MASS FRACTION
C       RATEYR(ISPEC) = SPECIES REACTION RATE
C       STRHYR(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1YR = DVDY
C       BCL2YR = DRHODY
C       BCL3YR = DUDY
C       BCL4YR = DWDY
C       BCL5YR = DPDY
C       BCLYYR(ISPEC) = DYDY

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              STRHYR(IC,KC,ISPEC) = STRHYR(IC,KC,ISPEC)
     +         - STRGYR(IC,KC)*STRTYR(IC,KC)*RGSPEC(ISPEC)/STRRYR(IC,KC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            GAM1YR(IC,KC) = STRGYR(IC,KC) - STRRYR(IC,KC)
            STREYR(IC,KC) = STREYR(IC,KC) - GAM1YR(IC,KC)*STRTYR(IC,KC)

            GAM1YR(IC,KC) = STRRYR(IC,KC)/GAM1YR(IC,KC)
            OVGMYR(IC,KC) = ONE/GAM1YR(IC,KC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = STRGYR(IC,KC)*GAM1YR(IC,KC)*STRTYR(IC,KC)
            ACOUYR(IC,KC) = SQRT(FORNOW)
            OVA2YR(IC,KC) = ONE/FORNOW

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCYR.EQ.NSBCO1)THEN

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYR(IC,KC) = SORPYR(IC,KC)
     +        + STRHYR(IC,KC,ISPEC)*RATEYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = -SORPYR(IC,KC)*GAM1YR(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L1Y AS REQUIRED
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L1Y
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +        *(BCL5YR(IC,KC)-STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC))

C             SUBTRACT FROM NEW VALUE OF L1Y
              BCL1YR(IC,KC)= HALF*SORPYR(IC,KC)
     +                     + COBCYR*ACOUYR(IC,KC)*(STRPYR(IC,KC)-PINFYR)
     +                     - BCL1YR(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)

              URHS(IC,JSTOL,KC) = URHS(IC,JSTOL,KC)
     +                   - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRUYR(IC,KC)

              VRHS(IC,JSTOL,KC) = VRHS(IC,JSTOL,KC)
     +   - BCL1YR(IC,KC)*OVA2YR(IC,KC)*(STRVYR(IC,KC)-ACOUYR(IC,KC))

              WRHS(IC,JSTOL,KC) = WRHS(IC,JSTOL,KC)
     +                   - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRWYR(IC,KC)

              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +               - BCL1YR(IC,KC)*(OVA2YR(IC,KC)*STREYR(IC,KC)
     +                              - STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                              + OVGMYR(IC,KC))

            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
     +                 - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRYYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCYR.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYR(IC,KC) = SORPYR(IC,KC)
     +        + STRHYR(IC,KC,ISPEC)*RATEYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = -SORPYR(IC,KC)*GAM1YR(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L4Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                      *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
C             L5Y UNCHANGED
              BCL1YR(IC,KC) = HALF*SORPYR(IC,KC)
     +                     + COBCYR*ACOUYR(IC,KC)*(STRPYR(IC,KC)-PINFYR)
     +                      - BCL1YR(IC,KC)
              BCL2YR(IC,KC) = -BCL2YR(IC,KC)
              BCL3YR(IC,KC) = -BCL3YR(IC,KC)
              BCL4YR(IC,KC) = -BCL4YR(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF L's
                BCLYYR(IC,KC,ISPEC) = STRVYR(IC,KC)*BCLYYR(IC,KC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
                BCLYYR(IC,KC,ISPEC) = RATEYR(IC,KC,ISPEC)/STRDYR(IC,KC)
     +                              - BCLYYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)
     +                          - BCL2YR(IC,KC)

              URHS(IC,JSTOL,KC) = URHS(IC,JSTOL,KC)
     +                      - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRUYR(IC,KC)
     +                          - BCL2YR(IC,KC)*STRUYR(IC,KC)
     +                          - BCL3YR(IC,KC)*STRDYR(IC,KC)

              VRHS(IC,JSTOL,KC) = VRHS(IC,JSTOL,KC)
     +      - BCL1YR(IC,KC)*OVA2YR(IC,KC)*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                          - BCL2YR(IC,KC)*STRVYR(IC,KC)

              WRHS(IC,JSTOL,KC) = WRHS(IC,JSTOL,KC)
     +                      - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRWYR(IC,KC)
     +                          - BCL2YR(IC,KC)*STRWYR(IC,KC)
     +                          - BCL4YR(IC,KC)*STRDYR(IC,KC)

              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                      - BCL1YR(IC,KC)*(OVA2YR(IC,KC)*STREYR(IC,KC)
     +                                     + STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                                     + OVGMYR(IC,KC))
     +                          - BCL2YR(IC,KC)*STREYR(IC,KC)
     +                      - BCL3YR(IC,KC)*STRDYR(IC,KC)*STRUYR(IC,KC)
     +                      - BCL4YR(IC,KC)*STRDYR(IC,KC)*STRWYR(IC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = BCLYYR(IC,KC,ISPEC)*STRDYR(IC,KC)

                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                            - FORNOW*STRHYR(IC,KC,ISPEC)

                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
     + - (BCL2YR(IC,KC)+BCL1YR(IC,KC)*OVA2YR(IC,KC))*STRYYR(IC,KC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYR.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SYDTYR(IC,KC) = ZERO
              SORPYR(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SYDTYR(IC,KC) = SYDTYR(IC,KC)
     +                        + DYDTYR(IC,KC,ISPEC)*RGSPEC(ISPEC)
                SORPYR(IC,KC) = SORPYR(IC,KC)
     +                        + STRHYR(IC,KC,ISPEC)*RATEYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SYDTYR(IC,KC) = SYDTYR(IC,KC)/STRRYR(IC,KC)
              SORPYR(IC,KC) = -SORPYR(IC,KC)*GAM1YR(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y,L2Y,L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                      *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL5YR(IC,KC) = HALF*(STRVYR(IC,KC)+ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Y UNCHANGED
              BCL1YR(IC,KC) = BCL5YR(IC,KC)
     +                      + STRDYR(IC,KC)*ACOUYR(IC,KC)*DVDTYR(IC,KC)
     +                      - BCL1YR(IC,KC)
              BCL2YR(IC,KC) = GAM1YR(IC,KC)*OVA2YR(IC,KC)
     +                       *(BCL1YR(IC,KC)+BCL5YR(IC,KC))
     +                      + STRDYR(IC,KC)*(DTDTYR(IC,KC)/STRTYR(IC,KC)
     +                      - SORPYR(IC,KC)/STRPYR(IC,KC)
     +                      + SYDTYR(IC,KC))
     +                      - BCL2YR(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)
     +                          - BCL2YR(IC,KC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCYR.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                    *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)
              BCL5YR(IC,KC) = HALF*(STRVYR(IC,KC)+ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Y UNCHANGED
              FORNOW = BCL5YR(IC,KC)
     +               + STRDYR(IC,KC)*ACOUYR(IC,KC)*DVDTYR(IC,KC)
              BCL1YR(IC,KC) = FORNOW - BCL1YR(IC,KC)
              BCL2YR(IC,KC) = -DDDTYR(IC,KC)
     +                      - OVA2YR(IC,KC)*(BCL1YR(IC,KC)+FORNOW)
     +                      - BCL2YR(IC,KC)
              BCL3YR(IC,KC) = -DUDTYR(IC,KC) - BCL3YR(IC,KC)
              BCL4YR(IC,KC) = -DWDTYR(IC,KC) - BCL4YR(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                  - BCL1YR(IC,KC)*(OVA2YR(IC,KC)*STREYR(IC,KC)
     +                                 + STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                                 + OVGMYR(IC,KC))
     +                          - BCL2YR(IC,KC)*STREYR(IC,KC)
     +                  - BCL3YR(IC,KC)*STRDYR(IC,KC)*STRUYR(IC,KC)
     +                  - BCL4YR(IC,KC)*STRDYR(IC,KC)*STRWYR(IC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                BCLYYR(IC,KC,ISPEC) = RATEYR(IC,KC,ISPEC)/STRDYR(IC,KC)
     +                              - DYDTYR(IC,KC,ISPEC)
     +                              - STRVYR(IC,KC)*BCLYYR(IC,KC,ISPEC)

                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +          - BCLYYR(IC,KC,ISPEC)*STRDYR(IC,KC)*STRHYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

        IF(NSBCYR.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 1
C         NO-SLIP WALL - ADIABATIC

C         ALL VELOCITY COMPONENTS IMPOSED
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Y,L3Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)
              BCL5YR(IC,KC) = HALF*(STRVYR(IC,KC)+ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L2Y,L5Y UNCHANGED
              BCL1YR(IC,KC) = BCL5YR(IC,KC)
     +                      + STRDYR(IC,KC)*ACOUYR(IC,KC)*DVDTYR(IC,KC)
     +                      - BCL1YR(IC,KC)
              BCL3YR(IC,KC) = -DUDTYR(IC,KC) - BCL3YR(IC,KC)
              BCL4YR(IC,KC) = -DWDTYR(IC,KC) - BCL4YR(IC,KC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)

              ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +                  - BCL1YR(IC,KC)*(OVA2YR(IC,KC)*STREYR(IC,KC)
     +                                 + STRVYR(IC,KC)/ACOUYR(IC,KC)
     +                                 + OVGMYR(IC,KC))
     +                  - BCL3YR(IC,KC)*STRDYR(IC,KC)*STRUYR(IC,KC)
     +                  - BCL4YR(IC,KC)*STRDYR(IC,KC)*STRWYR(IC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
     +                 - BCL1YR(IC,KC)*OVA2YR(IC,KC)*STRYYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCYR.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                SORPYR(IC,KC) = SORPYR(IC,KC)
     +                        + STRHYR(IC,KC,ISPEC)*RATEYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              SORPYR(IC,KC) = -SORPYR(IC,KC)*GAM1YR(IC,KC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L5Y
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDYR(IC,KC)*ACOUYR(IC,KC)*BCL1YR(IC,KC)
              BCL1YR(IC,KC) = HALF*(STRVYR(IC,KC)-ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)-FORNOW)
              BCL2YR(IC,KC) = STRVYR(IC,KC)
     +                    *(BCL2YR(IC,KC)-BCL5YR(IC,KC)*OVA2YR(IC,KC))
              BCL3YR(IC,KC) = STRVYR(IC,KC)*BCL3YR(IC,KC)
              BCL4YR(IC,KC) = STRVYR(IC,KC)*BCL4YR(IC,KC)
              BCL5YR(IC,KC) = HALF*(STRVYR(IC,KC)+ACOUYR(IC,KC))
     +                      *(BCL5YR(IC,KC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Y UNCHANGED
              BCL1YR(IC,KC) = BCL5YR(IC,KC)
     +                      + STRDYR(IC,KC)*ACOUYR(IC,KC)*DVDTYR(IC,KC)
     +                      - BCL1YR(IC,KC)
              BCL3YR(IC,KC) = -DUDTYR(IC,KC) - BCL3YR(IC,KC)
              BCL4YR(IC,KC) = -DWDTYR(IC,KC) - BCL4YR(IC,KC)
              BCL2YR(IC,KC) = GAM1YR(IC,KC)*OVA2YR(IC,KC)
     +                       *(BCL1YR(IC,KC)+BCL5YR(IC,KC))
     +                      + STRDYR(IC,KC)*(DTDTYR(IC,KC)/STRTYR(IC,KC)
     +                      - SORPYR(IC,KC)/STRPYR(IC,KC))
     +                      - BCL2YR(IC,KC)

            ENDDO
          ENDDO

C         LYY
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYY
                BCLYYR(IC,KC,ISPEC) = STRVYR(IC,KC)*BCLYYR(IC,KC,ISPEC)

C               UPDATE L2Y
                BCL2YR(IC,KC) = BCL2YR(IC,KC)
     +                        + (RATEYR(IC,KC,ISPEC)
     +                         - STRDYR(IC,KC)*BCLYYR(IC,KC,ISPEC))
     +                          *RGSPEC(ISPEC)/STRRYR(IC,KC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)
     +                          - BCL1YR(IC,KC)*OVA2YR(IC,KC)
     +                          - BCL2YR(IC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTOL,KC,ISPEC) = YRHS(IC,JSTOL,KC,ISPEC)
     + - (BCL2YR(IC,KC)+BCL1YR(IC,KC)*OVA2YR(IC,KC))*STRYYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     Y-DIRECTION RIGHT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     Z-DIRECTION LEFT-HAND END
C     -------------------------
      IF(FZLCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUZL = PRIMITIVE U-VELOCITY COMPONENT
C       STRVZL = PRIMITIVE V-VELOCITY COMPONENT
C       STRWZL = PRIMITIVE W-VELOCITY COMPONENT
C       STRPZL = PRESSURE
C       STRDZL = DENSITY
C       STRTZL = TEMPERATURE
C       STREZL = INTERNAL ENERGY
C       STRGZL = MIXTURE CP
C       STRRZL = MIXTURE SPECIFIC GAS CONSTANT
C       STRYZL(ISPEC) = SPECIES MASS FRACTION
C       RATEZL(ISPEC) = SPECIES REACTION RATE
C       STRHZL(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1ZL = DWDZ
C       BCL2ZL = DRHODZ
C       BCL3ZL = DUDZ
C       BCL4ZL = DVDZ
C       BCL5ZL = DPDZ
C       BCLYZL(ISPEC) = DYDZ

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRHZL(IC,JC,ISPEC) = STRHZL(IC,JC,ISPEC)
     +         - STRGZL(IC,JC)*STRTZL(IC,JC)*RGSPEC(ISPEC)/STRRZL(IC,JC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            GAM1ZL(IC,JC) = STRGZL(IC,JC) - STRRZL(IC,JC)
            STREZL(IC,JC) = STREZL(IC,JC) - GAM1ZL(IC,JC)*STRTZL(IC,JC)

            GAM1ZL(IC,JC) = STRRZL(IC,JC)/GAM1ZL(IC,JC)
            OVGMZL(IC,JC) = ONE/GAM1ZL(IC,JC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = STRGZL(IC,JC)*GAM1ZL(IC,JC)*STRTZL(IC,JC)
            ACOUZL(IC,JC) = SQRT(FORNOW)
            OVA2ZL(IC,JC) = ONE/FORNOW

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCZL.EQ.NSBCO1)THEN

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZL(IC,JC) = SORPZL(IC,JC)
     +        + STRHZL(IC,JC,ISPEC)*RATEZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = -SORPZL(IC,JC)*GAM1ZL(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L5Z AS REQUIRED
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L5Z
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +        *(BCL5ZL(IC,JC)+STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC))

C             SUBTRACT FROM NEW VALUE OF L5Z
              BCL5ZL(IC,JC)= HALF*SORPZL(IC,JC)
     +                     + COBCZL*ACOUZL(IC,JC)*(STRPZL(IC,JC)-PINFZL)
     +                     - BCL5ZL(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)

              URHS(IC,JC,KSTAL) = URHS(IC,JC,KSTAL)
     +                   - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRUZL(IC,JC)

              VRHS(IC,JC,KSTAL) = VRHS(IC,JC,KSTAL)
     +                   - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRVZL(IC,JC)

              WRHS(IC,JC,KSTAL) = WRHS(IC,JC,KSTAL)
     +   - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*(STRWZL(IC,JC)+ACOUZL(IC,JC))

              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +               - BCL5ZL(IC,JC)*(OVA2ZL(IC,JC)*STREZL(IC,JC)
     +                              + STRWZL(IC,JC)/ACOUZL(IC,JC)
     +                              + OVGMZL(IC,JC))

            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
     +                 - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRYZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCZL.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZL(IC,JC) = SORPZL(IC,JC)
     +        + STRHZL(IC,JC,ISPEC)*RATEZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = -SORPZL(IC,JC)*GAM1ZL(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L2Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                      *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
C             L1Z UNCHANGED
              BCL2ZL(IC,JC) = -BCL2ZL(IC,JC)
              BCL3ZL(IC,JC) = -BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = -BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*SORPZL(IC,JC)
     +                     + COBCZL*ACOUZL(IC,JC)*(STRPZL(IC,JC)-PINFZL)
     +                      - BCL5ZL(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF L's
                BCLYZL(IC,JC,ISPEC) = STRWZL(IC,JC)*BCLYZL(IC,JC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
                BCLYZL(IC,JC,ISPEC) = RATEZL(IC,JC,ISPEC)/STRDZL(IC,JC)
     +                              - BCLYZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)

              URHS(IC,JC,KSTAL) = URHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STRUZL(IC,JC)
     +                          - BCL3ZL(IC,JC)*STRDZL(IC,JC)
     +                      - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRUZL(IC,JC)

              VRHS(IC,JC,KSTAL) = VRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STRVZL(IC,JC)
     +                          - BCL4ZL(IC,JC)*STRDZL(IC,JC)
     +                      - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRVZL(IC,JC)

              WRHS(IC,JC,KSTAL) = WRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STRWZL(IC,JC)
     +      - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*(STRWZL(IC,JC)+ACOUZL(IC,JC))

              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STREZL(IC,JC)
     +                      - BCL3ZL(IC,JC)*STRDZL(IC,JC)*STRUZL(IC,JC)
     +                      - BCL4ZL(IC,JC)*STRDZL(IC,JC)*STRVZL(IC,JC)
     +                      - BCL5ZL(IC,JC)*(OVA2ZL(IC,JC)*STREZL(IC,JC)
     +                                     + STRWZL(IC,JC)/ACOUZL(IC,JC)
     +                                     + OVGMZL(IC,JC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = BCLYZL(IC,JC,ISPEC)*STRDZL(IC,JC)

                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                            - FORNOW*STRHZL(IC,JC,ISPEC)

                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
     + - (BCL2ZL(IC,JC)+BCL5ZL(IC,JC)*OVA2ZL(IC,JC))*STRYZL(IC,JC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZL.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SYDTZL(IC,JC) = ZERO
              SORPZL(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SYDTZL(IC,JC) = SYDTZL(IC,JC)
     +                        + DYDTZL(IC,JC,ISPEC)*RGSPEC(ISPEC)
                SORPZL(IC,JC) = SORPZL(IC,JC)
     +                        + STRHZL(IC,JC,ISPEC)*RATEZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SYDTZL(IC,JC) = SYDTZL(IC,JC)/STRRZL(IC,JC)
              SORPZL(IC,JC) = -SORPZL(IC,JC)*GAM1ZL(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z,L2Z,L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL1ZL(IC,JC) = HALF*(STRWZL(IC,JC)-ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)-FORNOW)
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                      *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Z UNCHANGED
              BCL5ZL(IC,JC) = BCL1ZL(IC,JC)
     +                      - STRDZL(IC,JC)*ACOUZL(IC,JC)*DWDTZL(IC,JC)
     +                      - BCL5ZL(IC,JC)
              BCL2ZL(IC,JC) = GAM1ZL(IC,JC)*OVA2ZL(IC,JC)
     +                       *(BCL1ZL(IC,JC)+BCL5ZL(IC,JC))
     +                      + STRDZL(IC,JC)*(DTDTZL(IC,JC)/STRTZL(IC,JC)
     +                      - SORPZL(IC,JC)/STRPZL(IC,JC)
     +                      + SYDTZL(IC,JC))
     +                      - BCL2ZL(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCZL.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL1ZL(IC,JC) = HALF*(STRWZL(IC,JC)-ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)-FORNOW)
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                    *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Z UNCHANGED
              FORNOW = BCL1ZL(IC,JC)
     +               - STRDZL(IC,JC)*ACOUZL(IC,JC)*DWDTZL(IC,JC)
              BCL2ZL(IC,JC) = -DDDTZL(IC,JC)
     +                      - OVA2ZL(IC,JC)*(BCL1ZL(IC,JC)+FORNOW)
     +                      - BCL2ZL(IC,JC)
              BCL3ZL(IC,JC) = -DUDTZL(IC,JC) - BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = -DVDTZL(IC,JC) - BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = FORNOW - BCL5ZL(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)*STREZL(IC,JC)
     +                  - BCL3ZL(IC,JC)*STRDZL(IC,JC)*STRUZL(IC,JC)
     +                  - BCL4ZL(IC,JC)*STRDZL(IC,JC)*STRVZL(IC,JC)
     +                  - BCL5ZL(IC,JC)*(OVA2ZL(IC,JC)*STREZL(IC,JC)
     +                                 + STRWZL(IC,JC)/ACOUZL(IC,JC)
     +                                 + OVGMZL(IC,JC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                BCLYZL(IC,JC,ISPEC) = RATEZL(IC,JC,ISPEC)/STRDZL(IC,JC)
     +                              - DYDTZL(IC,JC,ISPEC)
     +                              - STRWZL(IC,JC)*BCLYZL(IC,JC,ISPEC)

                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +          - BCLYZL(IC,JC,ISPEC)*STRDZL(IC,JC)*STRHZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

        IF(NSBCZL.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 1
C         NO-SLIP WALL - ADIABATIC

C         ALL VELOCITY COMPONENTS IMPOSED
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Z,L3Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL1ZL(IC,JC) = HALF*(STRWZL(IC,JC)-ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)-FORNOW)
              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Z,L2Z UNCHANGED
              BCL3ZL(IC,JC) = -DUDTZL(IC,JC) - BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = -DVDTZL(IC,JC) - BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = BCL1ZL(IC,JC)
     +                      - STRDZL(IC,JC)*ACOUZL(IC,JC)*DWDTZL(IC,JC)
     +                      - BCL5ZL(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)

              ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +                  - BCL3ZL(IC,JC)*STRDZL(IC,JC)*STRUZL(IC,JC)
     +                  - BCL4ZL(IC,JC)*STRDZL(IC,JC)*STRVZL(IC,JC)
     +                  - BCL5ZL(IC,JC)*(OVA2ZL(IC,JC)*STREZL(IC,JC)
     +                                 + STRWZL(IC,JC)/ACOUZL(IC,JC)
     +                                 + OVGMZL(IC,JC))

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
     +                 - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)*STRYZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCZL.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZL(IC,JC) = SORPZL(IC,JC)
     +                        + STRHZL(IC,JC,ISPEC)*RATEZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZL(IC,JC) = -SORPZL(IC,JC)*GAM1ZL(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZL(IC,JC)*ACOUZL(IC,JC)*BCL1ZL(IC,JC)
              BCL1ZL(IC,JC) = HALF*(STRWZL(IC,JC)-ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)-FORNOW)
              BCL2ZL(IC,JC) = STRWZL(IC,JC)
     +                    *(BCL2ZL(IC,JC)-BCL5ZL(IC,JC)*OVA2ZL(IC,JC))
              BCL3ZL(IC,JC) = STRWZL(IC,JC)*BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = STRWZL(IC,JC)*BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = HALF*(STRWZL(IC,JC)+ACOUZL(IC,JC))
     +                      *(BCL5ZL(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L1Y UNCHANGED
              BCL3ZL(IC,JC) = -DUDTZL(IC,JC) - BCL3ZL(IC,JC)
              BCL4ZL(IC,JC) = -DVDTZL(IC,JC) - BCL4ZL(IC,JC)
              BCL5ZL(IC,JC) = BCL1ZL(IC,JC)
     +                      - STRDZL(IC,JC)*ACOUZL(IC,JC)*DWDTZL(IC,JC)
     +                      - BCL5ZL(IC,JC)
              BCL2ZL(IC,JC) = GAM1ZL(IC,JC)*OVA2ZL(IC,JC)
     +                       *(BCL1ZL(IC,JC)+BCL5ZL(IC,JC))
     +                      + STRDZL(IC,JC)*(DTDTZL(IC,JC)/STRTZL(IC,JC)
     +                      - SORPZL(IC,JC)/STRPZL(IC,JC))
     +                      - BCL2ZL(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYZ
                BCLYZL(IC,JC,ISPEC) = STRWZL(IC,JC)*BCLYZL(IC,JC,ISPEC)

C               UPDATE L2Z
                BCL2ZL(IC,JC) = BCL2ZL(IC,JC)
     +                        + (RATEZL(IC,JC,ISPEC)
     +                         - STRDZL(IC,JC)*BCLYZL(IC,JC,ISPEC))
     +                          *RGSPEC(ISPEC)/STRRZL(IC,JC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)
     +                          - BCL2ZL(IC,JC)
     +                          - BCL5ZL(IC,JC)*OVA2ZL(IC,JC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTAL,ISPEC) = YRHS(IC,JC,KSTAL,ISPEC)
     + - (BCL2ZL(IC,JC)+BCL5ZL(IC,JC)*OVA2ZL(IC,JC))*STRYZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     Z-DIRECTION LEFT-HAND END

C     =========================================================================
C     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     =========================================================================

C     Z-DIRECTION RIGHT-HAND END
C     --------------------------
      IF(FZRCNV)THEN

C       =======================================================================

C       STR ARRAYS CONTAIN STORED VALUES
C       STRUZR = PRIMITIVE U-VELOCITY COMPONENT
C       STRVZR = PRIMITIVE V-VELOCITY COMPONENT
C       STRWZR = PRIMITIVE W-VELOCITY COMPONENT
C       STRPZR = PRESSURE
C       STRDZR = DENSITY
C       STRTZR = TEMPERATURE
C       STREZR = INTERNAL ENERGY
C       STRGZR = MIXTURE CP
C       STRRZR = MIXTURE SPECIFIC GAS CONSTANT
C       STRYZR(ISPEC) = SPECIES MASS FRACTION
C       RATEZR(ISPEC) = SPECIES REACTION RATE
C       STRHZR(ISPEC) = SPECIES ENTHALPY

C       BCL ARRAYS CONTAIN FIRST DERIVATIVES
C       BCL1ZR = DWDR
C       BCL2ZR = DRHODR
C       BCL3ZR = DUDR
C       BCL4ZR = DVDR
C       BCL5ZR = DPDR
C       BCLYZR(ISPEC) = DYDR

C       =======================================================================

C       REDUCED SPECIES ENTHALPY
C       ------------------------
        DO ISPEC = 1,NSPEC

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              STRHZR(IC,JC,ISPEC) = STRHZR(IC,JC,ISPEC)
     +         - STRGZR(IC,JC)*STRTZR(IC,JC)*RGSPEC(ISPEC)/STRRZR(IC,JC)

            ENDDO
          ENDDO

        ENDDO

C       REDUCED INTERNAL ENERGY
C       -----------------------
C       GAMMA-1, 1/(GAMMA-1)
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            GAM1ZR(IC,JC) = STRGZR(IC,JC) - STRRZR(IC,JC)
            STREZR(IC,JC) = STREZR(IC,JC) - GAM1ZR(IC,JC)*STRTZR(IC,JC)

            GAM1ZR(IC,JC) = STRRZR(IC,JC)/GAM1ZR(IC,JC)
            OVGMZR(IC,JC) = ONE/GAM1ZR(IC,JC)

          ENDDO
        ENDDO

C       SPEED OF SOUND
C       --------------
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = STRGZR(IC,JC)*GAM1ZR(IC,JC)*STRTZR(IC,JC)
            ACOUZR(IC,JC) = SQRT(FORNOW)
            OVA2ZR(IC,JC) = ONE/FORNOW

          ENDDO
        ENDDO

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

        IF(NSBCZR.EQ.NSBCO1)THEN

C         OUTFLOW BC No 1
C         SUBSONIC NON-REFLECTING OUTFLOW
C         WITH OPTION TO SET PRESSURE AT INFINITY

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZR(IC,JC) = SORPZR(IC,JC)
     +        + STRHZR(IC,JC,ISPEC)*RATEZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = -SORPZR(IC,JC)*GAM1ZR(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L1Z AS REQUIRED
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L1Z
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +        *(BCL5ZR(IC,JC)-STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC))

C             SUBTRACT FROM NEW VALUE OF L1Z
              BCL1ZR(IC,JC)= HALF*SORPZR(IC,JC)
     +                     + COBCZR*ACOUZR(IC,JC)*(STRPZR(IC,JC)-PINFZR)
     +                     - BCL1ZR(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)

              URHS(IC,JC,KSTOL) = URHS(IC,JC,KSTOL)
     +                   - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRUZR(IC,JC)

              VRHS(IC,JC,KSTOL) = VRHS(IC,JC,KSTOL)
     +                   - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRVZR(IC,JC)

              WRHS(IC,JC,KSTOL) = WRHS(IC,JC,KSTOL)
     +   - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*(STRWZR(IC,JC)-ACOUZR(IC,JC))

              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +               - BCL1ZR(IC,JC)*(OVA2ZR(IC,JC)*STREZR(IC,JC)
     +                              - STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                              + OVGMZR(IC,JC))

            ENDDO
          ENDDO

C         RSC 08-AUG-2012 EVALUATE ALL SPECIES
C          DO ISPEC = 1,NSPM1
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
     +                 - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRYZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

        IF(NSBCZR.EQ.NSBCI1)THEN

C         INFLOW BC No 1
C         SUBSONIC NON-REFLECTING LAMINAR INFLOW

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZR(IC,JC) = SORPZR(IC,JC)
     +        + STRHZR(IC,JC,ISPEC)*RATEZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = -SORPZR(IC,JC)*GAM1ZR(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Y-L4Y
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                      *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)

C             SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
C             L5Z UNCHANGED
              BCL1ZR(IC,JC) = HALF*SORPZR(IC,JC)
     +                     + COBCZR*ACOUZR(IC,JC)*(STRPZR(IC,JC)-PINFZR)
     +                      - BCL1ZR(IC,JC)
              BCL2ZR(IC,JC) = -BCL2ZR(IC,JC)
              BCL3ZR(IC,JC) = -BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = -BCL4ZR(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF L's
                BCLYZR(IC,JC,ISPEC) = STRWZR(IC,JC)*BCLYZR(IC,JC,ISPEC)

C               SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
                BCLYZR(IC,JC,ISPEC) = RATEZR(IC,JC,ISPEC)/STRDZR(IC,JC)
     +                              - BCLYZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                          - BCL2ZR(IC,JC)

              URHS(IC,JC,KSTOL) = URHS(IC,JC,KSTOL)
     +                      - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRUZR(IC,JC)
     +                          - BCL2ZR(IC,JC)*STRUZR(IC,JC)
     +                          - BCL3ZR(IC,JC)*STRDZR(IC,JC)

              VRHS(IC,JC,KSTOL) = VRHS(IC,JC,KSTOL)
     +                      - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRVZR(IC,JC)
     +                          - BCL2ZR(IC,JC)*STRVZR(IC,JC)
     +                          - BCL4ZR(IC,JC)*STRDZR(IC,JC)

              WRHS(IC,JC,KSTOL) = WRHS(IC,JC,KSTOL)
     +      - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                          - BCL2ZR(IC,JC)*STRWZR(IC,JC)

              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                      - BCL1ZR(IC,JC)*(OVA2ZR(IC,JC)*STREZR(IC,JC)
     +                                     + STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                                     + OVGMZR(IC,JC))
     +                          - BCL2ZR(IC,JC)*STREZR(IC,JC)
     +                      - BCL3ZR(IC,JC)*STRDZR(IC,JC)*STRUZR(IC,JC)
     +                      - BCL4ZR(IC,JC)*STRDZR(IC,JC)*STRVZR(IC,JC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                FORNOW = BCLYZR(IC,JC,ISPEC)*STRDZR(IC,JC)

                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                            - FORNOW*STRHZR(IC,JC,ISPEC)

                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
     + - (BCL2ZR(IC,JC)+BCL5ZR(IC,JC)*OVA2ZR(IC,JC))*STRYZR(IC,JC,ISPEC)
     +                                  - FORNOW

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZR.EQ.NSBCI2)THEN

C         INFLOW BOUNDARY CONDITION No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         VELOCITY, TEMPERATURE AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SYDTZR(IC,JC) = ZERO
              SORPZR(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SYDTZR(IC,JC) = SYDTZR(IC,JC)
     +                        + DYDTZR(IC,JC,ISPEC)*RGSPEC(ISPEC)
                SORPZR(IC,JC) = SORPZR(IC,JC)
     +                        + STRHZR(IC,JC,ISPEC)*RATEZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SYDTZR(IC,JC) = SYDTZR(IC,JC)/STRRZR(IC,JC)
              SORPZR(IC,JC) = -SORPZR(IC,JC)*GAM1ZR(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z,L2Z,L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                      *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL5ZR(IC,JC) = HALF*(STRWZR(IC,JC)+ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Z UNCHANGED
              BCL1ZR(IC,JC) = BCL5ZR(IC,JC)
     +                      + STRDZR(IC,JC)*ACOUZR(IC,JC)*DWDTZR(IC,JC)
     +                      - BCL1ZR(IC,JC)
              BCL2ZR(IC,JC) = GAM1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                       *(BCL1ZR(IC,JC)+BCL5ZR(IC,JC))
     +                      + STRDZR(IC,JC)*(DTDTZR(IC,JC)/STRTZR(IC,JC)
     +                      - SORPZR(IC,JC)/STRPZR(IC,JC)
     +                      + SYDTZR(IC,JC))
     +                      - BCL2ZR(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                          - BCL2ZR(IC,JC)

            ENDDO
          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCZR.EQ.NSBCI3)THEN

C         INFLOW BOUNDARY CONDITION No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         VELOCITY, DENSITY AND MASS FRACTIONS IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                    *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Z UNCHANGED
              FORNOW = BCL5ZR(IC,JC)
     +               + STRDZR(IC,JC)*ACOUZR(IC,JC)*DWDTZR(IC,JC)
              BCL1ZR(IC,JC) = FORNOW - BCL1ZR(IC,JC)
              BCL2ZR(IC,JC) = -DDDTZR(IC,JC)
     +                      - OVA2ZR(IC,JC)*(BCL1ZR(IC,JC)+FORNOW)
     +                      - BCL2ZR(IC,JC)
              BCL3ZR(IC,JC) = -DUDTZR(IC,JC) - BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = -DVDTZR(IC,JC) - BCL4ZR(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                  - BCL1ZR(IC,JC)*(OVA2ZR(IC,JC)*STREZR(IC,JC)
     +                                 + STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                                 + OVGMZR(IC,JC))
     +                          - BCL2ZR(IC,JC)*STREZR(IC,JC)
     +                  - BCL3ZR(IC,JC)*STRDZR(IC,JC)*STRUZR(IC,JC)
     +                  - BCL4ZR(IC,JC)*STRDZR(IC,JC)*STRVZR(IC,JC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                BCLYZR(IC,JC,ISPEC) = RATEZR(IC,JC,ISPEC)/STRDZR(IC,JC)
     +                              - DYDTZR(IC,JC,ISPEC)
     +                              - STRWZR(IC,JC)*BCLYZR(IC,JC,ISPEC)

                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +          - BCLYZR(IC,JC,ISPEC)*STRDZR(IC,JC)*STRHZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

C       WALL BOUNDARY CONDITIONS
C       ------------------------

        IF(NSBCZR.EQ.NSBCW1)THEN

C         WALL BOUNDARY CONDITION No 1
C         NO-SLIP WALL - ADIABATIC

C         ALL VELOCITY COMPONENTS IMPOSED
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         SPECIFY L's AS REQUIRED
C         L1Z,L3Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)
              BCL5ZR(IC,JC) = HALF*(STRWZR(IC,JC)+ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L2Z,L5Z UNCHANGED
              BCL1ZR(IC,JC) = BCL5ZR(IC,JC)
     +                      + STRDZR(IC,JC)*ACOUZR(IC,JC)*DWDTZR(IC,JC)
     +                      - BCL1ZR(IC,JC)
              BCL3ZR(IC,JC) = -DUDTZR(IC,JC) - BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = -DVDTZR(IC,JC) - BCL4ZR(IC,JC)

            ENDDO
          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)

              ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +                  - BCL1ZR(IC,JC)*(OVA2ZR(IC,JC)*STREZR(IC,JC)
     +                                 + STRWZR(IC,JC)/ACOUZR(IC,JC)
     +                                 + OVGMZR(IC,JC))
     +                  - BCL3ZR(IC,JC)*STRDZR(IC,JC)*STRUZR(IC,JC)
     +                  - BCL4ZR(IC,JC)*STRDZR(IC,JC)*STRVZR(IC,JC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
     +                 - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)*STRYZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

        IF(NSBCZR.EQ.NSBCW2)THEN

C         WALL BOUNDARY CONDITION No 2
C         NO-SLIP WALL - ISOTHERMAL

C         VELOCITY AND TEMPERATURE IMPOSED
C         AS FUNCTIONS OF TIME
C         VALUES AND TIME DERIVATIVES OF PRIMITIVE VARIABLES
C         SET IN SUBROUTINE BOUNDT

C         PRECOMPUTE CHEMISTRY TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = ZERO

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                SORPZR(IC,JC) = SORPZR(IC,JC)
     +                        + STRHZR(IC,JC,ISPEC)*RATEZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              SORPZR(IC,JC) = -SORPZR(IC,JC)*GAM1ZR(IC,JC)

            ENDDO
          ENDDO

C         SPECIFY L's AS REQUIRED
C         L1Z-L5Z
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL
           
C             OLD VALUE OF L's
              FORNOW = STRDZR(IC,JC)*ACOUZR(IC,JC)*BCL1ZR(IC,JC)
              BCL1ZR(IC,JC) = HALF*(STRWZR(IC,JC)-ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)-FORNOW)
              BCL2ZR(IC,JC) = STRWZR(IC,JC)
     +                    *(BCL2ZR(IC,JC)-BCL5ZR(IC,JC)*OVA2ZR(IC,JC))
              BCL3ZR(IC,JC) = STRWZR(IC,JC)*BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = STRWZR(IC,JC)*BCL4ZR(IC,JC)
              BCL5ZR(IC,JC) = HALF*(STRWZR(IC,JC)+ACOUZR(IC,JC))
     +                      *(BCL5ZR(IC,JC)+FORNOW)

C             SUBTRACT FROM NEW VALUE OF L's
C             L5Z UNCHANGED
              BCL1ZR(IC,JC) = BCL5ZR(IC,JC)
     +                      + STRDZR(IC,JC)*ACOUZR(IC,JC)*DWDTZR(IC,JC)
     +                      - BCL1ZR(IC,JC)
              BCL3ZR(IC,JC) = -DUDTZR(IC,JC) - BCL3ZR(IC,JC)
              BCL4ZR(IC,JC) = -DVDTZR(IC,JC) - BCL4ZR(IC,JC)
              BCL2ZR(IC,JC) = GAM1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                       *(BCL1ZR(IC,JC)+BCL5ZR(IC,JC))
     +                      + STRDZR(IC,JC)*(DTDTZR(IC,JC)/STRTZR(IC,JC)
     +                      - SORPZR(IC,JC)/STRPZR(IC,JC))
     +                      - BCL2ZR(IC,JC)

            ENDDO
          ENDDO

C         LYZ
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

C               OLD VALUE OF LYZ
                BCLYZR(IC,JC,ISPEC) = STRWZR(IC,JC)*BCLYZR(IC,JC,ISPEC)

C               UPDATE L2Z
                BCL2ZR(IC,JC) = BCL2ZR(IC,JC)
     +                        + (RATEZR(IC,JC,ISPEC)
     +                         - STRDZR(IC,JC)*BCLYZR(IC,JC,ISPEC))
     +                          *RGSPEC(ISPEC)/STRRZR(IC,JC)

              ENDDO
            ENDDO

          ENDDO

C         ADD TO CONSERVATIVE SOURCE TERMS
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)
     +                          - BCL1ZR(IC,JC)*OVA2ZR(IC,JC)
     +                          - BCL2ZR(IC,JC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTOL,ISPEC) = YRHS(IC,JC,KSTOL,ISPEC)
     + - (BCL2ZR(IC,JC)+BCL1ZR(IC,JC)*OVA2ZR(IC,JC))*STRYZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF  

C       =======================================================================

      ENDIF
C     Z-DIRECTION RIGHT-HAND END

C     =========================================================================


      RETURN
      END
