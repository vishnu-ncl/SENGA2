SUBROUTINE tempin_kernel_main(trun,drhs,urhs,vrhs,wrhs,erhs,yrhs,amascp,amasct,ncpoly,ncenth,tinthi,ntint,trin,nspec,iproc,idx)

    real(kind=8), dimension(1), intent(in) :: drhs,urhs,vrhs,wrhs,erhs,yrhs
    real(kind=8), dimension(1) :: trun

    real(kind=8), intent(in) :: trin
    integer(kind=4), dimension(3), intent(in) :: idx

    integer(kind=4), intent(in) :: nspec,iproc

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amascp,amasct
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncenth
    real(kind=8), dimension(ntinmx,nspcmx), intent(in) :: tinthi
    integer(kind=4), dimension(nspcmx), intent(in) :: ntint

!   PARAMETERS
!   ==========
    real(kind=8), parameter :: toltmp=0.00010_8
    real(kind=8), parameter :: tininc=50.0_8
    real(kind=8), parameter :: tlimlo=200.0_8
    real(kind=8), parameter :: tlimhi=3000.0_8

!   LOCAL DATA
!   ==========
    real(kind=8) :: tcoeff(0:nctmax)
    real(kind=8) :: ukuk
    real(kind=8) :: tempor,tupper,tlower,tresid,tuk2me,cpfory
    integer(kind=4) :: ic,jc,kc,ispec,itint,icp
    integer(kind=4) :: iindex,ipower
    logical :: fnconv

    ic = idx(1)
    jc = idx(2)
    kc = idx(3)

!   KINETIC ENERGY TERM
    ukuk = (urhs(OPS_ACC3(0,0,0))*urhs(OPS_ACC3(0,0,0))  &
          + vrhs(OPS_ACC4(0,0,0))*vrhs(OPS_ACC4(0,0,0))  &
          + wrhs(OPS_ACC5(0,0,0))*wrhs(OPS_ACC5(0,0,0)))/drhs(OPS_ACC2(0,0,0))

    tuk2me = 0.5_8*ukuk - erhs(OPS_ACC6(0,0,0))

!   ===================================================================

!   SOLVE FOR TEMPERATURE
!   =====================
!   USING BISECTION
    tlower = tlimlo
    tupper = tlimhi

!   SET NON-CONVERGENCE FLAG
    fnconv = .true.

!   SET TEMPERATURE FROM INITIAL GUESS
    tempor = trin

!   ===================================================================

!   FIRST BRACKET THE ROOT
!   ======================

!   -------------------------------------------------------------------

!   INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
    tcoeff(0) = tuk2me
    DO icp = 1, nctmax
        tcoeff(icp) = 0.0_8
    END DO

!   FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
    DO ispec = 1, nspec

        itint = 1
        DO WHILE (tempor > tinthi(itint,ispec) .and. itint < ntint(ispec))
            itint = itint + 1
        END DO

!       CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
        tcoeff(0) = tcoeff(0) + yrhs(OPS_ACC_MD7(ispec,0,0,0))*  &
                                amascp(ncenth(itint,ispec),itint,ispec)

        tcoeff(1) = tcoeff(1) + yrhs(OPS_ACC_MD7(ispec,0,0,0))*amasct(1,itint,ispec)
        DO icp = 2, ncpoly(itint,ispec)
            tcoeff(icp) = tcoeff(icp)  &
                        + yrhs(OPS_ACC_MD7(ispec,0,0,0))*amasct(icp,itint,ispec)
        END DO

    END DO
!   END OF RUN THROUGH ALL SPECIES

!   -------------------------------------------------------------------

!   EVALUATE TEMPERATURE RESIDUAL
    tresid = tcoeff(nctmax)
    DO icp = nctmm1,1,-1
        tresid = tcoeff(icp) + tresid*tempor
    END DO
    tresid = tcoeff(0) + tresid*tempor

!   -------------------------------------------------------------------

!   CHECK INITIAL GUESS FOR CONVERGENCE
    IF(ABS(tresid) < toltmp) THEN
!   -----------------------------------------------------------------
!       CONVERGED ON FIRST PASS
        fnconv = .false.
!   -----------------------------------------------------------------
    ELSE IF(tresid < 0.0_8)THEN
!   -----------------------------------------------------------------
!       INITIAL GUESS IS TOO LOW 2000
        DO

!           SET INITIAL GUESS AS LOWER LIMIT AND TRY AGAIN
            tlower = tempor
            tempor = tempor + tininc

!           ---------------------------------------------------------------
!           INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
            tcoeff(0) = tuk2me
            DO icp = 1, nctmax
                tcoeff(icp) = 0.0_8
            END DO

!           FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
            DO ispec = 1, nspec

                itint = 1
                DO WHILE (tempor > tinthi(itint,ispec) .and. itint < ntint(ispec))
                    itint = itint + 1
                END DO

!               CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
                tcoeff(0) = tcoeff(0) + yrhs(OPS_ACC_MD7(ispec,0,0,0))*  &
                                        amascp(ncenth(itint,ispec),itint,ispec)
                tcoeff(1) = tcoeff(1) + yrhs(OPS_ACC_MD7(ispec,0,0,0))*amasct(1,itint,ispec)

                DO icp = 2, ncpoly(itint,ispec)
                    tcoeff(icp) = tcoeff(icp)  &
                                + yrhs(OPS_ACC_MD7(ispec,0,0,0))*amasct(icp,itint,ispec)
                END DO

            END DO
!           END OF RUN THROUGH ALL SPECIES

!           ---------------------------------------------------------------

!           EVALUATE TEMPERATURE RESIDUAL
            tresid = tcoeff(nctmax)
            DO icp = nctmm1,1,-1
                tresid = tcoeff(icp) + tresid*tempor
            END DO
            tresid = tcoeff(0) + tresid*tempor

!           ---------------------------------------------------------------
!           CHECK NEW GUESS FOR CONVERGENCE
            IF(ABS(tresid) < toltmp) THEN
!           -------------------------------------------------------------
!               NEW GUESS HAS CONVERGED
                fnconv = .false.
                EXIT
!           -------------------------------------------------------------
            ELSE IF(tresid < 0.0_8) THEN
!           -------------------------------------------------------------
!               NEW GUESS IS STILL TOO LOW: GO ROUND AGAIN
                IF(tempor >= tlimhi) THEN
                    print *,'Fatal: TEMPIN: T upper bracket failed to converge'
                    print *,'processor:',iproc
                    print *,'at point:',ic,jc,kc
                    print *,'with values:',tempor,tresid
                    print *,drhs(OPS_ACC2(0,0,0))
                    print *,urhs(OPS_ACC3(0,0,0))
                    print *,vrhs(OPS_ACC4(0,0,0))
                    print *,wrhs(OPS_ACC5(0,0,0))
                    print *,erhs(OPS_ACC6(0,0,0))
                    DO ispec = 1, nspec
                        print *,yrhs(OPS_ACC_MD7(ispec,0,0,0))
                    END DO
                    STOP
                END IF
!           -------------------------------------------------------------
            ELSE IF(tresid > 0.0_8) THEN
!           -------------------------------------------------------------
!               ROOT IS BRACKETED
                tupper = tempor
                EXIT
!           -------------------------------------------------------------
            END IF
!       ---------------------------------------------------------------
        END DO
!       END OF LOOP 2000

!       -----------------------------------------------------------------
    ELSE IF(tresid > 0.0_8)THEN
!       -----------------------------------------------------------------
!       INITIAL GUESS IS TOO HIGH 3000
        DO

!           SET INITIAL GUESS AS UPPER LIMIT AND TRY AGAIN
            tupper = tempor
            tempor = tempor - tininc
!           ---------------------------------------------------------------
!           INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
            tcoeff(0) = tuk2me
            DO icp = 1, nctmax
                tcoeff(icp) = 0.0_8
            END DO

!           FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
            DO ispec = 1, nspec

                itint = 1
                DO WHILE (tempor > tinthi(itint,ispec) .and. itint < ntint(ispec))
                    itint = itint + 1
                END DO

!               CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
                tcoeff(0) = tcoeff(0) + yrhs(OPS_ACC_MD7(ispec,0,0,0))*  &
                                        amascp(ncenth(itint,ispec),itint,ispec)
                tcoeff(1) = tcoeff(1) + yrhs(OPS_ACC_MD7(ispec,0,0,0))*amasct(1,itint,ispec)
                DO icp = 2, ncpoly(itint,ispec)
                    tcoeff(icp) = tcoeff(icp)  &
                                + yrhs(OPS_ACC_MD7(ispec,0,0,0))*amasct(icp,itint,ispec)
                END DO

            END DO
!           END OF RUN THROUGH ALL SPECIES

!           ---------------------------------------------------------------
!           EVALUATE TEMPERATURE RESIDUAL
            tresid = tcoeff(nctmax)
            DO icp = nctmm1,1,-1
                tresid = tcoeff(icp) + tresid*tempor
            END DO
            tresid = tcoeff(0) + tresid*tempor

!           ---------------------------------------------------------------
!           CHECK NEW GUESS FOR CONVERGENCE
            IF(ABS(tresid) < toltmp) THEN
!           -------------------------------------------------------------
!               NEW GUESS HAS CONVERGED
                fnconv = .false.
                EXIT
!           -------------------------------------------------------------
            ELSE IF(tresid > 0.0_8) THEN
!           -------------------------------------------------------------
!               NEW GUESS IS STILL TOO HIGH: GO ROUND AGAIN
                IF(tempor <= tlimlo) THEN
                    print *,'Fatal: TEMPIN: T lower bracket failed to converge'
                    print *,'processor:',iproc
                    print *,'at point:',ic,jc,kc
                    print *,'with values:',tempor,tresid
                    print *,drhs(OPS_ACC2(0,0,0))
                    print *,urhs(OPS_ACC3(0,0,0))
                    print *,vrhs(OPS_ACC4(0,0,0))
                    print *,wrhs(OPS_ACC5(0,0,0))
                    print *,erhs(OPS_ACC6(0,0,0))
                    DO ispec = 1, nspec
                        print *,yrhs(OPS_ACC_MD7(ispec,0,0,0))
                    END DO
                    STOP
                END IF
!           -------------------------------------------------------------
            ELSE IF(tresid < 0.0_8) THEN
!           -------------------------------------------------------------
!               ROOT IS BRACKETED
                tlower = tempor
                EXIT
!           -------------------------------------------------------------
            END IF
!       ---------------------------------------------------------------
        END DO
!       END OF LOOP 3000
!       -----------------------------------------------------------------
    END IF
!   END OF CHECK INITIAL GUESS FOR CONVERGENCE

!   ===================================================================
!   ROOT IS BRACKETED
!   =================
!   NOW USE BISECTION TO REFINE THE ROOT

    IF(fnconv) THEN

!       infinite DO loop
        DO

!           BISECT
            tempor = 0.5_8*(tlower+tupper)
!           ---------------------------------------------------------------

!           INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
            tcoeff(0) = tuk2me
            DO icp = 1, nctmax
                tcoeff(icp) = 0.0_8
            END DO

!           FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
            DO ispec = 1, nspec

                itint = 1
                DO WHILE (tempor > tinthi(itint,ispec) .and. itint < ntint(ispec))
                    itint = itint + 1
                END DO

!               CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
                tcoeff(0) = tcoeff(0) + yrhs(OPS_ACC_MD7(ispec,0,0,0))*  &
                                        amascp(ncenth(itint,ispec),itint,ispec)
                tcoeff(1) = tcoeff(1) + yrhs(OPS_ACC_MD7(ispec,0,0,0))*amasct(1,itint,ispec)

                DO icp = 2, ncpoly(itint,ispec)
                    tcoeff(icp) = tcoeff(icp)  &
                                + yrhs(OPS_ACC_MD7(ispec,0,0,0))*amasct(icp,itint,ispec)
                END DO

            END DO
!           END OF RUN THROUGH ALL SPECIES

!           ---------------------------------------------------------------

!           EVALUATE TEMPERATURE RESIDUAL
            tresid = tcoeff(nctmax)
            DO icp = nctmm1,1,-1
                tresid = tcoeff(icp) + tresid*tempor
            END DO
            tresid = tcoeff(0) + tresid*tempor

!           ---------------------------------------------------------------

            IF(ABS(tresid) < toltmp)THEN
!               CONVERGED
                trun(OPS_ACC1(0,0,0)) = tempor
!               Break inifinite DO loop
                EXIT
            ELSE IF(tresid < 0.0_8) THEN
                tlower = tempor
            ELSE IF(tresid > 0.0_8) THEN
                tupper = tempor
            END IF

!       -----------------------------------------------------------------

        END DO
!       CONTINUE LOOP UNTIL CONVERGENCE OR TERMINATION

    END IF
!	END OF BISECTION

!   ===================================================================

!   SET THE NEW TEMPERATURE
    trun(OPS_ACC1(0,0,0)) = tempor

!   ===================================================================

END SUBROUTINE tempin_kernel_main
