SUBROUTINE temper_kernel_main(store7,trun,transp,prun,itndex,urhs,vrhs,wrhs,drhs,erhs,yrhs,amascp,amasct,ncenth,ncpom1,ncpoly,tinthi,rgspec,ntint,nctmax,ncofmx,ntinmx,nspcmx,nintmx,nspec,nspimx,ntbase,nctmm1,idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: urhs,vrhs,wrhs,drhs,erhs,yrhs
    real(kind=8), dimension(1) :: store7,trun,transp,prun
    integer, dimension(1) :: itndex
    integer, dimension(3), intent(in) :: idx

    integer, intent(in) :: nctmax,ncofmx,ntinmx,nspcmx,nintmx,nspec,nspimx,ntbase,nctmm1

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amascp,amasct 
    integer, dimension(ntinmx,nspec), intent(in) :: ncpoly,ncpom1,ncenth
    real(kind=8), dimension(nintmx,nspec), intent(in) :: tinthi
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer, dimension(nspec), intent(in) :: ntint


!   PARAMETERS
!   ==========
    real(kind=8), parameter :: toltmp=1.0E-10
    integer, parameter :: ntitrs=100

!   LOCAL DATA
!   ==========
    real(kind=8) :: ukuk
    real(kind=8) :: tempor,tfpoly,tdpoly,deltmp,cpfory
    real(kind=8) :: tcoeff(0:nctmax),tderiv(1:nctmax)

    integer :: ic,jc,kc,ispec,itint,icp,ititrs
    integer :: iindex,ipower,icoef1,icoef2

    ic = idx(1)
    jc = idx(2)
    kc = idx(3)

!   KINETIC ENERGY TERM
    ukuk = (urhs(OPS_ACC6(0,0,0))*urhs(OPS_ACC6(0,0,0))  &
          +  vrhs(OPS_ACC7(0,0,0))*vrhs(OPS_ACC7(0,0,0))  &
          +  wrhs(OPS_ACC8(0,0,0))*wrhs(OPS_ACC8(0,0,0)))/drhs(OPS_ACC9(0,0,0))

!   ===================================================================

!   INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
!   AND ITS DERIVATIVE
    tcoeff(0) = 0.5_8*ukuk - erhs(OPS_ACC10(0,0,0))
    DO icp = 1, nctmax
        tcoeff(icp) = 0.0_8
        tderiv(icp) = 0.0_8
    END DO

!   ===================================================================

!   USE STORE7 TO ACCUMULATE MIXTURE SPECIFIC GAS CONSTANT
!   INITIALISE STORE7
    store7(OPS_ACC1(0,0,0)) = 0.0_8

!   ===================================================================

!   RUN THROUGH ALL SPECIES
    DO ispec = 1,nspec

!       =================================================================

!       LOCATE TEMPERATURE IN AN INTERVAL
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1
        icoef2 = ntbase**ipower
        icoef1 = icoef2*ntbase
        itint = 1 + MOD(itndex(OPS_ACC_MD5(iindex,0,0,0)),icoef1)/icoef2

!       =================================================================

!       CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
        tcoeff(0) = tcoeff(0) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*  &
                    amascp(ncenth(itint,ispec),itint,ispec)
        tcoeff(1) = tcoeff(1) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*amasct(1,itint,ispec)
        tderiv(1) = tcoeff(1)
        DO icp = 2, ncpoly(itint,ispec)
            tcoeff(icp) = tcoeff(icp) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*  &
                          amasct(icp,itint,ispec)
            tderiv(icp) = tderiv(icp) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*  &
                          amascp(icp,itint,ispec)
        END DO

!       =================================================================

!       USE STORE7
!       TO ACCUMULATE (DENSITY TIMES) MIXTURE SPECIFIC GAS CONSTANT
        store7(OPS_ACC1(0,0,0)) = store7(OPS_ACC1(0,0,0))  &
            + yrhs(OPS_ACC_MD11(ispec,0,0,0))*rgspec(ispec)

!       =================================================================

    END DO

!   END OF RUN THROUGH ALL SPECIES

!    ===================================================================

!   SOLVE FOR TEMPERATURE
!   USING NEWTON-RAPHSON
    tempor = trun(OPS_ACC2(0,0,0))
    ititrs = 1

!   EVALUATE TEMPERATURE POLYNOMIAL AND ITS DERIVATIVE
    tfpoly = tcoeff(nctmax)
    tdpoly = tderiv(nctmax)

    DO icp = nctmm1,1,-1
        tfpoly = tcoeff(icp) + tfpoly*tempor
        tdpoly = tderiv(icp) + tdpoly*tempor
    END DO

    tfpoly = tcoeff(0) + tfpoly*tempor

!   EVALUATE TEMPERATURE CORRECTION
    deltmp = -tfpoly/tdpoly

!   CHECK FOR CONVERGENCE
    DO WHILE (ABS(deltmp) > toltmp)
        IF (ititrs < ntitrs) THEN
            tempor = tempor + deltmp
            ititrs = ititrs + 1

!           RE-EVALUATE TEMPERATURE POLYNOMIAL AND ITS DERIVATIVE
            tfpoly = tcoeff(nctmax)
            tdpoly = tderiv(nctmax)

            DO icp = nctmm1,1,-1
                tfpoly = tcoeff(icp) + tfpoly*tempor
                tdpoly = tderiv(icp) + tdpoly*tempor
            END DO

            tfpoly = tcoeff(0) + tfpoly*tempor

!           RE-EVALUATE TEMPERATURE CORRECTION
            deltmp = -tfpoly/tdpoly

        ELSE
            WRITE(6,*) 'Fatal: TEMPER: T iteration failed to converge'
!            WRITE(6,*)'processor:',iproc
            WRITE(6,*) 'at point:',ic,jc,kc
            WRITE(6,*) 'with values:',tempor,deltmp
            WRITE(6,*) drhs(OPS_ACC9(0,0,0))
            WRITE(6,*) urhs(OPS_ACC6(0,0,0))
            WRITE(6,*) vrhs(OPS_ACC7(0,0,0))
            WRITE(6,*) wrhs(OPS_ACC8(0,0,0))
            WRITE(6,*) erhs(OPS_ACC10(0,0,0))
            STOP
        END IF
    END DO

!   ===================================================================

!   SET THE NEW TEMPERATURE
    trun(OPS_ACC2(0,0,0)) = tempor

!   ===================================================================

!   FOR ALL SPECIES RELOCATE TEMPERATURE IN AN INTERVAL
!   EVALUATE MIXTURE SPECIFIC HEAT CP
    DO iindex = 1,nintmx
        itndex(OPS_ACC_MD5(iindex,0,0,0)) = 0
    END DO

    transp(OPS_ACC3(0,0,0)) = 0.0_8

    DO ispec = 1,nspec
        itint = 1

        DO WHILE (trun(OPS_ACC2(0,0,0)) > tinthi(itint,ispec) .and. itint < ntint(ispec))
            itint = itint + 1
        END DO

!       SET THE TEMPERATURE INTERVAL INDEX
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1
        itndex(OPS_ACC_MD5(iindex,0,0,0)) = itndex(OPS_ACC_MD5(iindex,0,0,0))  &
                                    +(itint-1)*ntbase**ipower

!       =================================================================

!       EVALUATE MIXTURE SPECIFIC HEAT CP
        cpfory = amascp(ncpoly(itint,ispec),itint,ispec)

        DO icp = ncpom1(itint,ispec),1,-1
          cpfory = cpfory*trun(OPS_ACC2(0,0,0)) + amascp(icp,itint,ispec)
        END DO

        transp(OPS_ACC3(0,0,0)) = transp(OPS_ACC3(0,0,0)) + yrhs(OPS_ACC_MD11(ispec,0,0,0))*cpfory

    END DO

    transp(OPS_ACC3(0,0,0)) = transp(OPS_ACC3(0,0,0))/drhs(OPS_ACC9(0,0,0))    

!   ===================================================================

!   EVALUATE MIXTURE PRESSURE
    prun(OPS_ACC4(0,0,0)) = trun(OPS_ACC2(0,0,0))*store7(OPS_ACC1(0,0,0))

!   ===================================================================

END SUBROUTINE temper_kernel_main
