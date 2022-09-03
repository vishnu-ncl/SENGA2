      SUBROUTINE BOUNDT
 
C     *************************************************************************
C
C     BOUNDT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     29-SEP-2003:  CREATED
C     10-MAY-2015:  RSC WALL BCS UPDATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     APPLIES BOUNDARY CONDITIONS TO PRIMITIVE VARIABLES
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      USE com_senga
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FORNOW
      INTEGER JC,KC
      INTEGER ISPEC
      INTEGER IINDEX,IPOWER,ICOEF1,ICOEF2
      INTEGER ITINT,ICP


C     BEGIN
C     =====

C     =========================================================================
CKA   FIX INFLOW BC
CKA   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
      BTIME  = RKTIM(IRKSTP)
      FUPELC = .FALSE.

C     X-DIRECTION LEFT-HAND END
C     -------------------------

C     GLOBAL BC SUPPORT
C     TURBULENT INFLOW VELOCITY FIELD
      IF(FXLTRB)CALL BCUTXL

C     LOCAL BC SUPPORT
      IF(FXLCNV)THEN

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

C       OUTFLOW BC No 1
C       SUBSONIC NON-REFLECTING OUTFLOW
C       WITH OPTION TO SET PRESSURE AT INFINITY
C       REQUIRES NO ACTION HERE

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

C       INFLOW BC No 1
C       SUBSONIC NON-REFLECTING LAMINAR INFLOW
C       REQUIRES NO ACTION HERE

C       =======================================================================

        IF(NSBCXL.EQ.NSBCI2)THEN

C         INFLOW BC No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTXL

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTXL

C         SET TEMPERATURE INTERVAL INDEX
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(ISTAL,JC,KC,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
1000            CONTINUE
                  IF(STRTXL(JC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 1000
                    ENDIF
                  ENDIF
C               END OF LOOP 1000

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(ISTAL,JC,KC,IINDEX) = ITNDEX(ISTAL,JC,KC,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              URHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*STRUXL(JC,KC)
              VRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*STRVXL(JC,KC)
              WRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*STRWXL(JC,KC)
              ERHS(ISTAL,JC,KC) = HALF*(STRUXL(JC,KC)*STRUXL(JC,KC)
     +                                + STRVXL(JC,KC)*STRVXL(JC,KC)
     +                                + STRWXL(JC,KC)*STRWXL(JC,KC))
              ERHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*ERHS(ISTAL,JC,KC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTXL

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                ITINT = 1 +MOD(ITNDEX(ISTAL,JC,KC,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTXL(JC,KC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTXL(JC,KC)

                YRHS(ISTAL,JC,KC,ISPEC)
     +                           = DRHS(ISTAL,JC,KC)*STRYXL(JC,KC,ISPEC)

                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTXL(JC,KC))*YRHS(ISTAL,JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXL.EQ.NSBCI3)THEN 

C         INFLOW BC No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         SET DENSITY AND TIME DERIVATIVE
          CALL BCDTXL

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTXL

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTAL,JC,KC) = STRDXL(JC,KC)

              URHS(ISTAL,JC,KC) = STRDXL(JC,KC)*STRUXL(JC,KC)
              VRHS(ISTAL,JC,KC) = STRDXL(JC,KC)*STRVXL(JC,KC)
              WRHS(ISTAL,JC,KC) = STRDXL(JC,KC)*STRWXL(JC,KC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTXL

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTAL,JC,KC,ISPEC)
     +                    = STRDXL(JC,KC)*STRYXL(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXL.EQ.NSBCW1)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ADIABATIC

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTXL

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              URHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*STRUXL(JC,KC)
              VRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*STRVXL(JC,KC)
              WRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*STRWXL(JC,KC)

            ENDDO
          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXL.EQ.NSBCW2)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ISOTHERMAL

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTXL

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTXL

C         SET TEMPERATURE INTERVAL INDEX
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(ISTAL,JC,KC,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
1100            CONTINUE
                  IF(STRTXL(JC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 1100
                    ENDIF
                  ENDIF
C               END OF LOOP 1100

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(ISTAL,JC,KC,IINDEX) = ITNDEX(ISTAL,JC,KC,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              URHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*STRUXL(JC,KC)
              VRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*STRVXL(JC,KC)
              WRHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*STRWXL(JC,KC)
              ERHS(ISTAL,JC,KC) = HALF*(STRUXL(JC,KC)*STRUXL(JC,KC)
     +                                + STRVXL(JC,KC)*STRVXL(JC,KC)
     +                                + STRWXL(JC,KC)*STRWXL(JC,KC))
              ERHS(ISTAL,JC,KC) = DRHS(ISTAL,JC,KC)*ERHS(ISTAL,JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                ITINT = 1 +MOD(ITNDEX(ISTAL,JC,KC,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTXL(JC,KC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTXL(JC,KC)

                ERHS(ISTAL,JC,KC) = ERHS(ISTAL,JC,KC)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTXL(JC,KC))*YRHS(ISTAL,JC,KC,ISPEC)

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

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

C       OUTFLOW BC No 1
C       SUBSONIC NON-REFLECTING OUTFLOW
C       WITH OPTION TO SET PRESSURE AT INFINITY
C       REQUIRES NO ACTION HERE

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

C       INFLOW BC No 1
C       SUBSONIC NON-REFLECTING LAMINAR INFLOW
C       REQUIRES NO ACTION HERE

C       =======================================================================

        IF(NSBCXR.EQ.NSBCI2)THEN

C         INFLOW BC No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTXR

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTXR

C         SET TEMPERATURE INTERVAL INDEX
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(ISTOL,JC,KC,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
1500            CONTINUE
                  IF(STRTXR(JC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 1500
                    ENDIF
                  ENDIF
C               END OF LOOP 1500

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(ISTOL,JC,KC,IINDEX) = ITNDEX(ISTOL,JC,KC,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              URHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*STRUXR(JC,KC)
              VRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*STRVXR(JC,KC)
              WRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*STRWXR(JC,KC)
              ERHS(ISTOL,JC,KC) = HALF*(STRUXR(JC,KC)*STRUXR(JC,KC)
     +                                + STRVXR(JC,KC)*STRVXR(JC,KC)
     +                                + STRWXR(JC,KC)*STRWXR(JC,KC))
              ERHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*ERHS(ISTOL,JC,KC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTXR

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                ITINT = 1 +MOD(ITNDEX(ISTOL,JC,KC,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTXR(JC,KC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTXR(JC,KC)

                YRHS(ISTOL,JC,KC,ISPEC)
     +                           = DRHS(ISTOL,JC,KC)*STRYXR(JC,KC,ISPEC)

                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTXR(JC,KC))*YRHS(ISTOL,JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXR.EQ.NSBCI3)THEN 

C         INFLOW BC No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         SET DENSITY AND TIME DERIVATIVE
          CALL BCDTXR

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTXR

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DRHS(ISTOL,JC,KC) = STRDXR(JC,KC)

              URHS(ISTOL,JC,KC) = STRDXR(JC,KC)*STRUXR(JC,KC)
              VRHS(ISTOL,JC,KC) = STRDXR(JC,KC)*STRVXR(JC,KC)
              WRHS(ISTOL,JC,KC) = STRDXR(JC,KC)*STRWXR(JC,KC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTXR

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                YRHS(ISTOL,JC,KC,ISPEC)
     +                    = STRDXR(JC,KC)*STRYXR(JC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXR.EQ.NSBCW1)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ADIABATIC

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTXR

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              URHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*STRUXR(JC,KC)
              VRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*STRVXR(JC,KC)
              WRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*STRWXR(JC,KC)

            ENDDO
          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCXR.EQ.NSBCW2)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ISOTHERMAL

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTXR

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTXR

C         SET TEMPERATURE INTERVAL INDEX
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(ISTOL,JC,KC,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
1600            CONTINUE
                  IF(STRTXR(JC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 1600
                    ENDIF
                  ENDIF
C               END OF LOOP 1600

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(ISTOL,JC,KC,IINDEX) = ITNDEX(ISTOL,JC,KC,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL

              URHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*STRUXR(JC,KC)
              VRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*STRVXR(JC,KC)
              WRHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*STRWXR(JC,KC)
              ERHS(ISTOL,JC,KC) = HALF*(STRUXR(JC,KC)*STRUXR(JC,KC)
     +                                + STRVXR(JC,KC)*STRVXR(JC,KC)
     +                                + STRWXR(JC,KC)*STRWXR(JC,KC))
              ERHS(ISTOL,JC,KC) = DRHS(ISTOL,JC,KC)*ERHS(ISTOL,JC,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO KC = KSTAL,KSTOL
              DO JC = JSTAL,JSTOL

                ITINT = 1 +MOD(ITNDEX(ISTOL,JC,KC,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTXR(JC,KC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTXR(JC,KC)

                ERHS(ISTOL,JC,KC) = ERHS(ISTOL,JC,KC)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTXR(JC,KC))*YRHS(ISTOL,JC,KC,ISPEC)

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

C     GLOBAL BC SUPPORT
C     TURBULENT INFLOW VELOCITY FIELD
      IF(FYLTRB)CALL BCUTYL

C     LOCAL BC SUPPORT
      IF(FYLCNV)THEN

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

C       OUTFLOW BC No 1
C       SUBSONIC NON-REFLECTING OUTFLOW
C       WITH OPTION TO SET PRESSURE AT INFINITY
C       REQUIRES NO ACTION HERE

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

C       INFLOW BC No 1
C       SUBSONIC NON-REFLECTING LAMINAR INFLOW
C       REQUIRES NO ACTION HERE

C       =======================================================================

        IF(NSBCYL.EQ.NSBCI2)THEN

C         INFLOW BC No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTYL

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTYL

C         SET TEMPERATURE INTERVAL INDEX
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(IC,JSTAL,KC,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
2000            CONTINUE
                  IF(STRTYL(IC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 2000
                    ENDIF
                  ENDIF
C               END OF LOOP 2000

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(IC,JSTAL,KC,IINDEX) = ITNDEX(IC,JSTAL,KC,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*STRUYL(IC,KC)
              VRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*STRVYL(IC,KC)
              WRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*STRWYL(IC,KC)
              ERHS(IC,JSTAL,KC) = HALF*(STRUYL(IC,KC)*STRUYL(IC,KC)
     +                                + STRVYL(IC,KC)*STRVYL(IC,KC)
     +                                + STRWYL(IC,KC)*STRWYL(IC,KC))
              ERHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*ERHS(IC,JSTAL,KC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTYL

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                ITINT = 1 +MOD(ITNDEX(IC,JSTAL,KC,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTYL(IC,KC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTYL(IC,KC)

                YRHS(IC,JSTAL,KC,ISPEC)
     +                           = DRHS(IC,JSTAL,KC)*STRYYL(IC,KC,ISPEC)

                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTYL(IC,KC))*YRHS(IC,JSTAL,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYL.EQ.NSBCI3)THEN 

C         INFLOW BC No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         SET DENSITY AND TIME DERIVATIVE
          CALL BCDTYL

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTYL

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTAL,KC) = STRDYL(IC,KC)

              URHS(IC,JSTAL,KC) = STRDYL(IC,KC)*STRUYL(IC,KC)
              VRHS(IC,JSTAL,KC) = STRDYL(IC,KC)*STRVYL(IC,KC)
              WRHS(IC,JSTAL,KC) = STRDYL(IC,KC)*STRWYL(IC,KC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTYL

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTAL,KC,ISPEC)
     +                    = STRDYL(IC,KC)*STRYYL(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYL.EQ.NSBCW1)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ADIABATIC

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTYL

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*STRUYL(IC,KC)
              VRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*STRVYL(IC,KC)
              WRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*STRWYL(IC,KC)

            ENDDO
          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYL.EQ.NSBCW2)THEN 

C         WALL BC No 2
C         NO-SLIP WALL - ISOTHERMAL

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTYL

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTYL

C         SET TEMPERATURE INTERVAL INDEX
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(IC,JSTAL,KC,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
2100            CONTINUE
                  IF(STRTYL(IC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 2100
                    ENDIF
                  ENDIF
C               END OF LOOP 2100

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(IC,JSTAL,KC,IINDEX) = ITNDEX(IC,JSTAL,KC,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*STRUYL(IC,KC)
              VRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*STRVYL(IC,KC)
              WRHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*STRWYL(IC,KC)
              ERHS(IC,JSTAL,KC) = HALF*(STRUYL(IC,KC)*STRUYL(IC,KC)
     +                                + STRVYL(IC,KC)*STRVYL(IC,KC)
     +                                + STRWYL(IC,KC)*STRWYL(IC,KC))
              ERHS(IC,JSTAL,KC) = DRHS(IC,JSTAL,KC)*ERHS(IC,JSTAL,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                ITINT = 1 +MOD(ITNDEX(IC,JSTAL,KC,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTYL(IC,KC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTYL(IC,KC)

                ERHS(IC,JSTAL,KC) = ERHS(IC,JSTAL,KC)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTYL(IC,KC))*YRHS(IC,JSTAL,KC,ISPEC)

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

C     GLOBAL BC SUPPORT
C     TURBULENT INFLOW VELOCITY FIELD
      IF(FYRTRB)CALL BCUTYR

C     LOCAL BC SUPPORT
      IF(FYRCNV)THEN

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

C       OUTFLOW BC No 1
C       SUBSONIC NON-REFLECTING OUTFLOW
C       WITH OPTION TO SET PRESSURE AT INFINITY
C       REQUIRES NO ACTION HERE

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

C       INFLOW BC No 1
C       SUBSONIC NON-REFLECTING LAMINAR INFLOW
C       REQUIRES NO ACTION HERE

C       =======================================================================

        IF(NSBCYR.EQ.NSBCI2)THEN

C         INFLOW BC No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTYR

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTYR

C         SET TEMPERATURE INTERVAL INDEX
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(IC,JSTOL,KC,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
2500            CONTINUE
                  IF(STRTYR(IC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 2500
                    ENDIF
                  ENDIF
C               END OF LOOP 2500

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(IC,JSTOL,KC,IINDEX) = ITNDEX(IC,JSTOL,KC,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*STRUYR(IC,KC)
              VRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*STRVYR(IC,KC)
              WRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*STRWYR(IC,KC)
              ERHS(IC,JSTOL,KC) = HALF*(STRUYR(IC,KC)*STRUYR(IC,KC)
     +                                + STRVYR(IC,KC)*STRVYR(IC,KC)
     +                                + STRWYR(IC,KC)*STRWYR(IC,KC))
              ERHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*ERHS(IC,JSTOL,KC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTYR

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                ITINT = 1 +MOD(ITNDEX(IC,JSTOL,KC,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTYR(IC,KC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTYR(IC,KC)

                YRHS(IC,JSTOL,KC,ISPEC)
     +                           = DRHS(IC,JSTOL,KC)*STRYYR(IC,KC,ISPEC)

                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTYR(IC,KC))*YRHS(IC,JSTOL,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYR.EQ.NSBCI3)THEN 

C         INFLOW BC No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         SET DENSITY AND TIME DERIVATIVE
          CALL BCDTYR

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTYR

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JSTOL,KC) = STRDYR(IC,KC)

              URHS(IC,JSTOL,KC) = STRDYR(IC,KC)*STRUYR(IC,KC)
              VRHS(IC,JSTOL,KC) = STRDYR(IC,KC)*STRVYR(IC,KC)
              WRHS(IC,JSTOL,KC) = STRDYR(IC,KC)*STRWYR(IC,KC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTYR

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JSTOL,KC,ISPEC)
     +                    = STRDYR(IC,KC)*STRYYR(IC,KC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYR.EQ.NSBCW1)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ADIABATIC

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTYR

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*STRUYR(IC,KC)
              VRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*STRVYR(IC,KC)
              WRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*STRWYR(IC,KC)

            ENDDO
          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCYR.EQ.NSBCW2)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ISOTHERMAL

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTYR

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTYR

C         SET TEMPERATURE INTERVAL INDEX
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(IC,JSTOL,KC,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
2600            CONTINUE
                  IF(STRTYR(IC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 2600
                    ENDIF
                  ENDIF
C               END OF LOOP 2600

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(IC,JSTOL,KC,IINDEX) = ITNDEX(IC,JSTOL,KC,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO KC = KSTAL,KSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*STRUYR(IC,KC)
              VRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*STRVYR(IC,KC)
              WRHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*STRWYR(IC,KC)
              ERHS(IC,JSTOL,KC) = HALF*(STRUYR(IC,KC)*STRUYR(IC,KC)
     +                                + STRVYR(IC,KC)*STRVYR(IC,KC)
     +                                + STRWYR(IC,KC)*STRWYR(IC,KC))
              ERHS(IC,JSTOL,KC) = DRHS(IC,JSTOL,KC)*ERHS(IC,JSTOL,KC)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO KC = KSTAL,KSTOL
              DO IC = ISTAL,ISTOL

                ITINT = 1 +MOD(ITNDEX(IC,JSTOL,KC,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTYR(IC,KC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTYR(IC,KC)

                ERHS(IC,JSTOL,KC) = ERHS(IC,JSTOL,KC)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTYR(IC,KC))*YRHS(IC,JSTOL,KC,ISPEC)

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

C     GLOBAL BC SUPPORT
C     TURBULENT INFLOW VELOCITY FIELD
      IF(FZLTRB)CALL BCUTZL

C     LOCAL BC SUPPORT
      IF(FZLCNV)THEN

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

C       OUTFLOW BC No 1
C       SUBSONIC NON-REFLECTING OUTFLOW
C       WITH OPTION TO SET PRESSURE AT INFINITY
C       REQUIRES NO ACTION HERE

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

C       INFLOW BC No 1
C       SUBSONIC NON-REFLECTING LAMINAR INFLOW
C       REQUIRES NO ACTION HERE

C       =======================================================================

        IF(NSBCZL.EQ.NSBCI2)THEN

C         INFLOW BC No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTZL

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTZL

C         SET TEMPERATURE INTERVAL INDEX
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(IC,JC,KSTAL,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
3000            CONTINUE
                  IF(STRTZL(IC,JC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 3000
                    ENDIF
                  ENDIF
C               END OF LOOP 3000

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(IC,JC,KSTAL,IINDEX) = ITNDEX(IC,JC,KSTAL,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*STRUZL(IC,JC)
              VRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*STRVZL(IC,JC)
              WRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*STRWZL(IC,JC)
              ERHS(IC,JC,KSTAL) = HALF*(STRUZL(IC,JC)*STRUZL(IC,JC)
     +                                + STRVZL(IC,JC)*STRVZL(IC,JC)
     +                                + STRWZL(IC,JC)*STRWZL(IC,JC))
              ERHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*ERHS(IC,JC,KSTAL)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTZL

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                ITINT = 1 +MOD(ITNDEX(IC,JC,KSTAL,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTZL(IC,JC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTZL(IC,JC)

                YRHS(IC,JC,KSTAL,ISPEC)
     +                           = DRHS(IC,JC,KSTAL)*STRYZL(IC,JC,ISPEC)

                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTZL(IC,JC))*YRHS(IC,JC,KSTAL,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZL.EQ.NSBCI3)THEN 

C         INFLOW BC No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         SET DENSITY AND TIME DERIVATIVE
          CALL BCDTZL

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTZL

C         CONSERVATIVE VARIABLES
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTAL) = STRDZL(IC,JC)

              URHS(IC,JC,KSTAL) = STRDZL(IC,JC)*STRUZL(IC,JC)
              VRHS(IC,JC,KSTAL) = STRDZL(IC,JC)*STRVZL(IC,JC)
              WRHS(IC,JC,KSTAL) = STRDZL(IC,JC)*STRWZL(IC,JC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTZL

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTAL,ISPEC)
     +                    = STRDZL(IC,JC)*STRYZL(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZL.EQ.NSBCW1)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ADIABATIC

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTZL

C         CONSERVATIVE VARIABLES
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*STRUZL(IC,JC)
              VRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*STRVZL(IC,JC)
              WRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*STRWZL(IC,JC)

            ENDDO
          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZL.EQ.NSBCW2)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ISOTHERMAL

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTZL

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTZL

C         SET TEMPERATURE INTERVAL INDEX
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(IC,JC,KSTAL,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
3100            CONTINUE
                  IF(STRTZL(IC,JC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 3100
                    ENDIF
                  ENDIF
C               END OF LOOP 3100

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(IC,JC,KSTAL,IINDEX) = ITNDEX(IC,JC,KSTAL,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*STRUZL(IC,JC)
              VRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*STRVZL(IC,JC)
              WRHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*STRWZL(IC,JC)
              ERHS(IC,JC,KSTAL) = HALF*(STRUZL(IC,JC)*STRUZL(IC,JC)
     +                                + STRVZL(IC,JC)*STRVZL(IC,JC)
     +                                + STRWZL(IC,JC)*STRWZL(IC,JC))
              ERHS(IC,JC,KSTAL) = DRHS(IC,JC,KSTAL)*ERHS(IC,JC,KSTAL)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                ITINT = 1 +MOD(ITNDEX(IC,JC,KSTAL,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTZL(IC,JC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTZL(IC,JC)

                ERHS(IC,JC,KSTAL) = ERHS(IC,JC,KSTAL)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTZL(IC,JC))*YRHS(IC,JC,KSTAL,ISPEC)

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

C     GLOBAL BC SUPPORT
C     TURBULENT INFLOW VELOCITY FIELD
      IF(FZRTRB)CALL BCUTZR

C     LOCAL BC SUPPORT
      IF(FZRCNV)THEN

C       =======================================================================

C       OUTFLOW BOUNDARY CONDITIONS
C       ---------------------------

C       OUTFLOW BC No 1
C       SUBSONIC NON-REFLECTING OUTFLOW
C       WITH OPTION TO SET PRESSURE AT INFINITY
C       REQUIRES NO ACTION HERE

C       =======================================================================

C       INFLOW BOUNDARY CONDITIONS
C       --------------------------

C       INFLOW BC No 1
C       SUBSONIC NON-REFLECTING LAMINAR INFLOW
C       REQUIRES NO ACTION HERE

C       =======================================================================

        IF(NSBCZR.EQ.NSBCI2)THEN

C         INFLOW BC No 2
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTZR

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTZR

C         SET TEMPERATURE INTERVAL INDEX
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(IC,JC,KSTOL,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
3500            CONTINUE
                  IF(STRTZR(IC,JC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 3500
                    ENDIF
                  ENDIF
C               END OF LOOP 3500

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(IC,JC,KSTOL,IINDEX) = ITNDEX(IC,JC,KSTOL,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*STRUZR(IC,JC)
              VRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*STRVZR(IC,JC)
              WRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*STRWZR(IC,JC)
              ERHS(IC,JC,KSTOL) = HALF*(STRUZR(IC,JC)*STRUZR(IC,JC)
     +                                + STRVZR(IC,JC)*STRVZR(IC,JC)
     +                                + STRWZR(IC,JC)*STRWZR(IC,JC))
              ERHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*ERHS(IC,JC,KSTOL)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTZR

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                ITINT = 1 +MOD(ITNDEX(IC,JC,KSTOL,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTZR(IC,JC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTZR(IC,JC)

                YRHS(IC,JC,KSTOL,ISPEC)
     +                           = DRHS(IC,JC,KSTOL)*STRYZR(IC,JC,ISPEC)

                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTZR(IC,JC))*YRHS(IC,JC,KSTOL,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZR.EQ.NSBCI3)THEN 

C         INFLOW BC No 3
C         SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

C         SET DENSITY AND TIME DERIVATIVE
          CALL BCDTZR

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTZR

C         CONSERVATIVE VARIABLES
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DRHS(IC,JC,KSTOL) = STRDZR(IC,JC)

              URHS(IC,JC,KSTOL) = STRDZR(IC,JC)*STRUZR(IC,JC)
              VRHS(IC,JC,KSTOL) = STRDZR(IC,JC)*STRVZR(IC,JC)
              WRHS(IC,JC,KSTOL) = STRDZR(IC,JC)*STRWZR(IC,JC)

            ENDDO
          ENDDO

C         SET MASS FRACTIONS AND TIME DERIVATIVES
          CALL BCYTZR

C         CONSERVATIVE VARIABLES
          DO ISPEC = 1,NSPEC

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                YRHS(IC,JC,KSTOL,ISPEC)
     +                    = STRDZR(IC,JC)*STRYZR(IC,JC,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZR.EQ.NSBCW1)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ADIABATIC

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTZR

C         CONSERVATIVE VARIABLES
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*STRUZR(IC,JC)
              VRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*STRVZR(IC,JC)
              WRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*STRWZR(IC,JC)

            ENDDO
          ENDDO

        ENDIF

C       =======================================================================

        IF(NSBCZR.EQ.NSBCW2)THEN 

C         WALL BC No 1
C         NO-SLIP WALL - ISOTHERMAL

C         SET VELOCITY COMPONENTS AND TIME DERIVATIVES
          CALL BCUTZR

C         SET TEMPERATURE AND TIME DERIVATIVE
          CALL BCTTZR

C         SET TEMPERATURE INTERVAL INDEX
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              DO IINDEX = 1,NINTMX
                ITNDEX(IC,JC,KSTOL,IINDEX) = 0
              ENDDO

              DO ISPEC = 1,NSPEC

                ITINT = 1
3600            CONTINUE
                  IF(STRTZR(IC,JC).GT.TINTHI(ITINT,ISPEC))THEN
                    IF(ITINT.LT.NTINT(ISPEC))THEN
                      ITINT = ITINT + 1
                      GOTO 3600
                    ENDIF
                  ENDIF
C               END OF LOOP 3600

C               SET THE TEMPERATURE INTERVAL INDEX
                IINDEX = 1 + (ISPEC-1)/NSPIMX
                IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
                ITNDEX(IC,JC,KSTOL,IINDEX) = ITNDEX(IC,JC,KSTOL,IINDEX)
     +                                     + (ITINT-1)*NTBASE**IPOWER

              ENDDO

            ENDDO
          ENDDO

C         CONSERVATIVE VARIABLES
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              URHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*STRUZR(IC,JC)
              VRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*STRVZR(IC,JC)
              WRHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*STRWZR(IC,JC)
              ERHS(IC,JC,KSTOL) = HALF*(STRUZR(IC,JC)*STRUZR(IC,JC)
     +                                + STRVZR(IC,JC)*STRVZR(IC,JC)
     +                                + STRWZR(IC,JC)*STRWZR(IC,JC))
              ERHS(IC,JC,KSTOL) = DRHS(IC,JC,KSTOL)*ERHS(IC,JC,KSTOL)

            ENDDO
          ENDDO

          DO ISPEC = 1,NSPEC

C           TEMPERATURE INTERVAL INDEXING 
            IINDEX = 1 + (ISPEC-1)/NSPIMX
            IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
            ICOEF2 = NTBASE**IPOWER
            ICOEF1 = ICOEF2*NTBASE

            DO JC = JSTAL,JSTOL
              DO IC = ISTAL,ISTOL

                ITINT = 1 +MOD(ITNDEX(IC,JC,KSTOL,IINDEX),ICOEF1)/ICOEF2
                FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
                DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                  FORNOW = FORNOW*STRTZR(IC,JC)
     +                   + AMASCH(ICP,ITINT,ISPEC)
                ENDDO
                FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +                 + FORNOW*STRTZR(IC,JC)

                ERHS(IC,JC,KSTOL) = ERHS(IC,JC,KSTOL)
     +    + (FORNOW-RGSPEC(ISPEC)*STRTZR(IC,JC))*YRHS(IC,JC,KSTOL,ISPEC)

              ENDDO
            ENDDO

          ENDDO

        ENDIF

      ENDIF
C     Z-DIRECTION RIGHT-HAND END

C     =========================================================================


      RETURN
      END
