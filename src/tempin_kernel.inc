SUBROUTINE tempin_kernel_main(trun,store7,transp,prun,itndex,urhs,vrhs,wrhs,drhs,erhs,yrhs,amascp,amasct,ncpoly,ncpom1,ncenth,tinthi,rgspec,ntint,trin,nctmax,ncofmx,ntinmx,nspcmx,nintmx,nspec,nspimx,ntbase,nctmm1,iproc,idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: urhs,vrhs,wrhs,drhs,erhs,yrhs
    real(kind=8), dimension(1) :: trun,store7,transp,prun
    integer, dimension(1) :: itndex

    real(kind=8), intent(in) :: trin
    integer, dimension(3), intent(in) :: idx    

    integer, intent(in) :: nctmax,ncofmx,ntinmx,nspcmx,nintmx,nspec,nspimx,ntbase,nctmm1,iproc
    
    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amascp,amasct
    integer, dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth
    real(kind=8), dimension(nintmx,nspcmx), intent(in) :: tinthi
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer, dimension(nspcmx), intent(in) :: ntint

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
    integer :: ic,jc,kc,ispec,itint,icp
    integer :: iindex,ipower
    logical :: fnconv    

    ic = idx(1)
    jc = idx(2)
    kc = idx(3)

!   KINETIC ENERGY TERM
    ukuk = (urhs(OPS_ACC6(0,0,0))*urhs(OPS_ACC6(0,0,0))  &
          +  vrhs(OPS_ACC7(0,0,0))*vrhs(OPS_ACC7(0,0,0))  &
          +  wrhs(OPS_ACC8(0,0,0))*wrhs(OPS_ACC8(0,0,0)))/drhs(OPS_ACC9(0,0,0))

    tuk2me = 0.5_8*ukuk - erhs(OPS_ACC10(0,0,0))

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
        tcoeff(0) = tcoeff(0) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*  &
                                amascp(ncenth(itint,ispec),itint,ispec)
        
        tcoeff(1) = tcoeff(1) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*amasct(1,itint,ispec)
        DO icp = 2, ncpoly(itint,ispec)
            tcoeff(icp) = tcoeff(icp)  &
                        + yrhs(OPS_ACC_MD11(ispec,0,0,0))*amasct(icp,itint,ispec)
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
!       INITIAL GUESS IS TOO LOW
        2000          CONTINUE

!       SET INITIAL GUESS AS LOWER LIMIT AND TRY AGAIN
        tlower = tempor
        tempor = tempor + tininc

!       ---------------------------------------------------------------
!       INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
        tcoeff(0) = tuk2me
        DO icp = 1, nctmax
            tcoeff(icp) = 0.0_8
        END DO

!       FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
        DO ispec = 1, nspec

            itint = 1
            DO WHILE (tempor > tinthi(itint,ispec) .and. itint < ntint(ispec))
                itint = itint + 1
            END DO

!           CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
            tcoeff(0) = tcoeff(0) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*  &
                                    amascp(ncenth(itint,ispec),itint,ispec)
            tcoeff(1) = tcoeff(1) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*amasct(1,itint,ispec)

            DO icp = 2, ncpoly(itint,ispec)
                tcoeff(icp) = tcoeff(icp)  &
                            + yrhs(OPS_ACC_MD11(ispec,0,0,0))*amasct(icp,itint,ispec)
            END DO

        END DO
!       END OF RUN THROUGH ALL SPECIES

!       ---------------------------------------------------------------

!       EVALUATE TEMPERATURE RESIDUAL
        tresid = tcoeff(nctmax)
        DO icp = nctmm1,1,-1
            tresid = tcoeff(icp) + tresid*tempor
        END DO
        tresid = tcoeff(0) + tresid*tempor

!       ---------------------------------------------------------------
!       CHECK NEW GUESS FOR CONVERGENCE
        IF(ABS(tresid) < toltmp) THEN
!       -------------------------------------------------------------
!           NEW GUESS HAS CONVERGED
            fnconv = .false.
!       -------------------------------------------------------------
        ELSE IF(tresid < 0.0_8) THEN
!       -------------------------------------------------------------
!           NEW GUESS IS STILL TOO LOW: GO ROUND AGAIN
            IF(tempor < tlimhi) THEN
                GO TO 2000
            ELSE
                print *,'Fatal: TEMPIN: T upper bracket failed to converge'
                print *,'processor:',iproc
                print *,'at point:',ic,jc,kc
                print *,'with values:',tempor,tresid
                print *,drhs(OPS_ACC9(0,0,0))
                print *,urhs(OPS_ACC6(0,0,0))
                print *,vrhs(OPS_ACC7(0,0,0))
                print *,wrhs(OPS_ACC8(0,0,0))
                print *,erhs(OPS_ACC10(0,0,0))
                DO ispec = 1, nspec
                    print *,yrhs(OPS_ACC_MD11(ispec,0,0,0))
                END DO
                STOP
            END IF
!       -------------------------------------------------------------
        ELSE IF(tresid > 0.0_8) THEN
!       -------------------------------------------------------------
!           ROOT IS BRACKETED
            tupper = tempor
!       -------------------------------------------------------------
        END IF
!       ---------------------------------------------------------------
!       END OF LOOP 2000

!       -----------------------------------------------------------------
        ELSE IF(tresid > 0.0_8)THEN
!       -----------------------------------------------------------------
!       INITIAL GUESS IS TOO HIGH
            3000          CONTINUE
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
                tcoeff(0) = tcoeff(0) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*  &
                                        amascp(ncenth(itint,ispec),itint,ispec)
                tcoeff(1) = tcoeff(1) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*amasct(1,itint,ispec)
                DO icp = 2, ncpoly(itint,ispec)
                    tcoeff(icp) = tcoeff(icp)  &
                                + yrhs(OPS_ACC_MD11(ispec,0,0,0))*amasct(icp,itint,ispec)
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
!           -------------------------------------------------------------
            ELSE IF(tresid > 0.0_8) THEN
!           -------------------------------------------------------------
!               NEW GUESS IS STILL TOO HIGH: GO ROUND AGAIN
                IF(tempor > tlimlo) THEN
                    GO TO 3000
                ELSE
                    print *,'Fatal: TEMPIN: T lower bracket failed to converge'
                    print *,'processor:',iproc
                    print *,'at point:',ic,jc,kc
                    print *,'with values:',tempor,tresid
                    print *,drhs(OPS_ACC9(0,0,0))
                    print *,urhs(OPS_ACC6(0,0,0))
                    print *,vrhs(OPS_ACC7(0,0,0))
                    print *,wrhs(OPS_ACC8(0,0,0))
                    print *,erhs(OPS_ACC10(0,0,0))
                    DO ispec = 1, nspec
                        print *,yrhs(OPS_ACC_MD11(ispec,0,0,0))
                    END DO
                    STOP
                END IF
!           -------------------------------------------------------------
            ELSE IF(tresid < 0.0_8) THEN
!           -------------------------------------------------------------
!               ROOT IS BRACKETED
                tlower = tempor
!           -------------------------------------------------------------
        END IF
!       ---------------------------------------------------------------
!       END OF LOOP 3000
!       -----------------------------------------------------------------
    END IF
!   END OF CHECK INITIAL GUESS FOR CONVERGENCE

!   ===================================================================
!   ROOT IS BRACKETED
!   =================
!   NOW USE BISECTION TO REFINE THE ROOT

    IF(fnconv) THEN

        4000          CONTINUE

!       BISECT
        tempor = 0.5_8*(tlower+tupper)
!       ---------------------------------------------------------------

!       INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
        tcoeff(0) = tuk2me
        DO icp = 1, nctmax
            tcoeff(icp) = 0.0_8
        END DO

!       FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
        DO ispec = 1, nspec

            itint = 1
            DO WHILE (tempor > tinthi(itint,ispec) .and. itint < ntint(ispec))
                itint = itint + 1
            END DO

!           CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
            tcoeff(0) = tcoeff(0) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*  &
                                    amascp(ncenth(itint,ispec),itint,ispec)
            tcoeff(1) = tcoeff(1) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*amasct(1,itint,ispec)

            DO icp = 2, ncpoly(itint,ispec)
                tcoeff(icp) = tcoeff(icp)  &
                            + yrhs(OPS_ACC_MD11(ispec,0,0,0))*amasct(icp,itint,ispec)
            END DO

        END DO
!       END OF RUN THROUGH ALL SPECIES

!       ---------------------------------------------------------------

!       EVALUATE TEMPERATURE RESIDUAL
        tresid = tcoeff(nctmax)
        DO icp = nctmm1,1,-1
            tresid = tcoeff(icp) + tresid*tempor
        END DO
        tresid = tcoeff(0) + tresid*tempor

!       ---------------------------------------------------------------

        IF(ABS(tresid) < toltmp)THEN
!           CONVERGED
            trun(OPS_ACC1(0,0,0)) = tempor

        ELSE IF(tresid < 0.0_8) THEN
            tlower = tempor
            GO TO 4000

        ELSE IF(tresid > 0.0_8) THEN
            tupper = tempor
            GO TO 4000
        END IF

!       -----------------------------------------------------------------
!       END OF LOOP 4000

        END IF
!       END OF BISECTION

!   ===================================================================

!   SET THE NEW TEMPERATURE
    trun(OPS_ACC1(0,0,0)) = tempor

!   ===================================================================

!   CONSTRUCT THE TEMPERATURE INTERVAL INDEX
!   EVALUATE PRESSURE
!   EVALUATE MIXTURE SPECIFIC HEAT CP
    DO iindex = 1,nintmx
        itndex(OPS_ACC_MD5(iindex,0,0,0)) = 0
    END DO

    store7(OPS_ACC2(0,0,0)) = 0.0_8
    transp(OPS_ACC3(0,0,0)) = 0.0_8

    DO ispec = 1, nspec

        itint = 1
        DO WHILE (trun(OPS_ACC1(0,0,0)) > tinthi(itint,ispec) .and. itint < ntint(ispec))
            itint = itint + 1
        END DO

!       SET THE TEMPERATURE INDEX
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1
        itndex(OPS_ACC_MD5(iindex,0,0,0)) = itndex(OPS_ACC_MD5(iindex,0,0,0))  &
                                + (itint-1)*ntbase**ipower

!       =================================================================

!       EVALUATE MIXTURE SPECIFIC HEAT CP
        cpfory = amascp(ncpoly(itint,ispec),itint,ispec)
        DO icp = ncpom1(itint,ispec),1,-1
            cpfory = cpfory*trun(OPS_ACC1(0,0,0)) + amascp(icp,itint,ispec)
        END DO

        transp(OPS_ACC3(0,0,0)) = transp(OPS_ACC3(0,0,0)) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*cpfory

!       =================================================================

!       EVALUATE (DENSITY TIMES) MIXTURE GAS CONSTANT FOR PRESSURE
        store7(OPS_ACC2(0,0,0)) = store7(OPS_ACC2(0,0,0))  &
                         + yrhs(OPS_ACC_MD11(ispec,0,0,0))*rgspec(ispec)

!   ===================================================================

    END DO

    transp(OPS_ACC3(0,0,0)) = transp(OPS_ACC3(0,0,0))/drhs(OPS_ACC9(0,0,0))

!   ===================================================================

!   EVALUATE PRESSURE
    prun(OPS_ACC4(0,0,0)) = trun(OPS_ACC1(0,0,0))*store7(OPS_ACC2(0,0,0))

!   ===================================================================

END SUBROUTINE tempin_kernel_main
