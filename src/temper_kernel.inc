SUBROUTINE temper_kernel_eqA(tcoeff, drhs, urhs, vrhs, wrhs, erhs)

    real(kind=8), dimension(1), intent(in) :: drhs, urhs, vrhs, wrhs, erhs
    real(kind=8), dimension(1) :: tcoeff

!   LOCAL DATA
!   ==========
    real(kind=8) :: ukuk
    integer(kind=4) :: icp

!   KINETIC ENERGY TERM
    ukuk = (urhs(OPS_ACC3(0,0,0))*urhs(OPS_ACC3(0,0,0))  &
          + vrhs(OPS_ACC4(0,0,0))*vrhs(OPS_ACC4(0,0,0))  &
          + wrhs(OPS_ACC5(0,0,0))*wrhs(OPS_ACC5(0,0,0)))/drhs(OPS_ACC2(0,0,0))

    tcoeff(OPS_ACC_MD1(1,0,0,0)) = 0.5_8*ukuk - erhs(OPS_ACC6(0,0,0))

    DO icp = 1, nctmax
        tcoeff(OPS_ACC_MD1(icp+1,0,0,0)) = 0.0_8
    END DO

END SUBROUTINE temper_kernel_eqA

SUBROUTINE temper_kernel_eqB(tcoeff, tderiv, itndex, yrhs, amascp, amasct, ncpoly, ncenth, ipower, icoef1, icoef2, ispec)

    real(kind=8), dimension(1), intent(in) :: yrhs
    real(kind=8), dimension(1) :: tcoeff, tderiv
    integer(kind=4), dimension(1), intent(in) :: itndex
    integer(kind=4), intent(in) :: ipower, icoef1, icoef2, ispec

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amascp,amasct
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly, ncenth

    integer(kind=4) :: itint, icp

!   LOCATE TEMPERATURE IN AN INTERVAL
    itint = 1 + MOD(itndex(OPS_ACC3(0,0,0)),icoef1)/icoef2

!   CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
    tcoeff(OPS_ACC_MD1(1,0,0,0)) = tcoeff(OPS_ACC_MD1(1,0,0,0)) + yrhs(OPS_ACC4(0,0,0))* &
                                    amascp(ncenth(itint,ispec),itint,ispec)

    tcoeff(OPS_ACC_MD1(2,0,0,0)) = tcoeff(OPS_ACC_MD1(2,0,0,0)) + yrhs(OPS_ACC4(0,0,0))*amasct(1,itint,ispec)

    tderiv(OPS_ACC_MD2(1,0,0,0)) = tcoeff(OPS_ACC_MD1(2,0,0,0))

    DO icp = 2, ncpoly(itint,ispec)
        tcoeff(OPS_ACC_MD1(icp+1,0,0,0)) = tcoeff(OPS_ACC_MD1(icp+1,0,0,0)) + yrhs(OPS_ACC4(0,0,0))*  &
                                                                                amasct(icp,itint,ispec)
        tderiv(OPS_ACC_MD2(icp,0,0,0))    = tderiv(OPS_ACC_MD2(icp,0,0,0)) + yrhs(OPS_ACC4(0,0,0))*  &
                                                                             amascp(icp,itint,ispec)
    END DO

END SUBROUTINE temper_kernel_eqB

SUBROUTINE temper_kernel_eqC(store7, yrhs, rgspec, ispec)

    real(kind=8), dimension(1), intent(in) :: yrhs
    real(kind=8), dimension(1) :: store7

    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    store7(OPS_ACC1(0,0,0)) = store7(OPS_ACC1(0,0,0)) + yrhs(OPS_ACC2(0,0,0))*rgspec(ispec)

END SUBROUTINE temper_kernel_eqC

SUBROUTINE temper_kernel_eqD(trun, tcoeff, tderiv, drhs, urhs, vrhs, wrhs, erhs, idx)

    real(kind=8), dimension(1), intent(in) :: tcoeff, tderiv, drhs, urhs, vrhs, wrhs, erhs
    real(kind=8), dimension(1) :: trun
    integer(kind=4), dimension(3), intent(in) :: idx

!   PARAMETERS
!   ==========
    real(kind=8), parameter :: toltmp=1.0E-10
    integer(kind=4), parameter :: ntitrs=100

!   LOCAL DATA
!   ==========
    real(kind=8) :: tempor,tfpoly,tdpoly,deltmp
    integer(kind=4) :: ititrs,icp, ic,jc,kc

    ic = idx(1)
    jc = idx(2)
    kc = idx(3)

!   USING NEWTON-RAPHSON
    tempor = trun(OPS_ACC1(0,0,0))
    ititrs = 1

!   EVALUATE TEMPERATURE POLYNOMIAL AND ITS DERIVATIVE
    tfpoly = tcoeff(OPS_ACC_MD2(nctmax+1,0,0,0))
    tdpoly = tderiv(OPS_ACC_MD3(nctmax,0,0,0))

    DO icp = nctmm1,1,-1
        tfpoly = tcoeff(OPS_ACC_MD2(icp+1,0,0,0)) + tfpoly*tempor
        tdpoly = tderiv(OPS_ACC_MD3(icp,0,0,0))   + tdpoly*tempor
    END DO

    tfpoly = tcoeff(OPS_ACC_MD2(1,0,0,0)) + tfpoly*tempor

!   EVALUATE TEMPERATURE CORRECTION
    deltmp = -tfpoly/tdpoly

!   CHECK FOR CONVERGENCE
    DO WHILE (ABS(deltmp) > toltmp)
        IF (ititrs < ntitrs) THEN
            tempor = tempor + deltmp
            ititrs = ititrs + 1

!           RE-EVALUATE TEMPERATURE POLYNOMIAL AND ITS DERIVATIVE
            tfpoly = tcoeff(OPS_ACC_MD2(nctmax+1,0,0,0))
            tdpoly = tderiv(OPS_ACC_MD3(nctmax,0,0,0))

            DO icp = nctmm1,1,-1
                tfpoly = tcoeff(OPS_ACC_MD2(icp+1,0,0,0)) + tfpoly*tempor
                tdpoly = tderiv(OPS_ACC_MD3(icp,0,0,0))   + tdpoly*tempor
            END DO

            tfpoly = tcoeff(OPS_ACC_MD2(1,0,0,0)) + tfpoly*tempor

!           RE-EVALUATE TEMPERATURE CORRECTION
            deltmp = -tfpoly/tdpoly

        ELSE
            print *,'Fatal: TEMPER: T iteration failed to converge'
!            WRITE(6,*)'processor:',iproc
            print *,'at point:',ic,jc,kc
            print *,'with values:',tempor,deltmp
            print *,drhs(OPS_ACC4(0,0,0))
            print *,urhs(OPS_ACC5(0,0,0))
            print *,vrhs(OPS_ACC6(0,0,0))
            print *,wrhs(OPS_ACC7(0,0,0))
            print *,erhs(OPS_ACC8(0,0,0))
            STOP
        END IF
    END DO

!   ===================================================================

!   SET THE NEW TEMPERATURE
    trun(OPS_ACC1(0,0,0)) = tempor

END SUBROUTINE temper_kernel_eqD

SUBROUTINE temper_kernel_eqE(transp, itndex, yrhs, trun, amascp, ncpoly, ncpom1, tinthi, ntint, ispec, ipower)
    
    real(kind=8), dimension(1), intent(in) :: yrhs, trun
    real(kind=8), dimension(1) :: transp
    integer(kind=4), dimension(1) :: itndex
    integer(kind=4), intent(in) :: ispec, ipower

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amascp
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1
    real(kind=8), dimension(ntinmx,nspcmx), intent(in) :: tinthi
    integer(kind=4), dimension(nspcmx), intent(in) :: ntint

!   LOCAL DATA
!   ==========
    real(kind=8) :: cpfory
    integer(kind=4) :: itint,icp

    itint = 1

    DO WHILE (trun(OPS_ACC4(0,0,0)) > tinthi(itint,ispec) .and. itint < ntint(ispec))
        itint = itint + 1
    END DO

    itndex(OPS_ACC2(0,0,0)) = itndex(OPS_ACC2(0,0,0))  &
                                    +(itint-1)*ntbase**ipower

!   =================================================================

!   EVALUATE MIXTURE SPECIFIC HEAT CP
    cpfory = amascp(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        cpfory = cpfory*trun(OPS_ACC4(0,0,0)) + amascp(icp,itint,ispec)
    END DO

    transp(OPS_ACC1(0,0,0)) = transp(OPS_ACC1(0,0,0)) + yrhs(OPS_ACC3(0,0,0))*cpfory

END SUBROUTINE temper_kernel_eqE

SUBROUTINE temper_kernel_eqF(transp, prun, drhs, trun, store7)

    real(kind=8), dimension(1), intent(in) :: drhs, trun, store7
    real(kind=8), dimension(1) :: transp, prun

    transp(OPS_ACC1(0,0,0)) = transp(OPS_ACC1(0,0,0))/drhs(OPS_ACC3(0,0,0))

!   EVALUATE MIXTURE PRESSURE
    prun(OPS_ACC2(0,0,0)) = trun(OPS_ACC4(0,0,0))*store7(OPS_ACC5(0,0,0))

END SUBROUTINE temper_kernel_eqF
