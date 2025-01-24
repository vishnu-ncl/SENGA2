SUBROUTINE inflow_kernel_eqA(uinf2,vinf2,winf2,ufilt,vfilt,wfilt,uinf1,vinf1,winf1,phi)

    real(kind=8), dimension(1) :: uinf2,vinf2,winf2
    real(kind=8), dimension(1), intent(in) :: ufilt,vfilt,wfilt,uinf1,vinf1
    real(kind=8), intent(in) :: phi

    uinf2(OPS_ACC1(0,0,0)) = ufilt(OPS_ACC4(0,0,0)) + phi*uinf1(OPS_ACC7(0,0,0))
    vinf2(OPS_ACC2(0,0,0)) = vfilt(OPS_ACC5(0,0,0)) + phi*vinf1(OPS_ACC8(0,0,0))
    winf2(OPS_ACC3(0,0,0)) = wfilt(OPS_ACC6(0,0,0)) + phi*winf1(OPS_ACC9(0,0,0))

END SUBROUTINE inflow_kernel_eqA

SUBROUTINE inflow_kernel_eqB(yinf2,yfilt,yinf1,phi)

    real(kind=8), dimension(1) :: yinf2
    real(kind=8), dimension(1), intent(in) :: yfilt,yinf1
    real(kind=8), intent(in) :: phi

    yinf2(OPS_ACC1(0,0,0)) = yfilt(OPS_ACC2(0,0,0)) + phi*yinf1(OPS_ACC3(0,0,0))

END SUBROUTINE inflow_kernel_eqB

SUBROUTINE inflow_kernel_eqC(ustead)

    real(kind=8), dimension(1) :: ustead

    ustead(OPS_ACC1(0,0,0)) = 1.0_8

END SUBROUTINE inflow_kernel_eqC

SUBROUTINE inflow_kernel_eqD(uinf2,vinf2,winf2,ustead)

    real(kind=8), dimension(1) :: uinf2,vinf2,winf2
    real(kind=8), dimension(1), intent(in) :: ustead

    uinf2(OPS_ACC1(0,0,0)) = uinf2(OPS_ACC1(0,0,0)) * ustead(OPS_ACC4(0,0,0))
    vinf2(OPS_ACC2(0,0,0)) = vinf2(OPS_ACC2(0,0,0)) * ustead(OPS_ACC4(0,0,0))
    winf2(OPS_ACC3(0,0,0)) = 0.0_8      ! winf2(OPS_ACC3(0,0,0))*ustead(OPS_ACC4(0,0,0))

END SUBROUTINE inflow_kernel_eqD

SUBROUTINE inflow_kernel_eqE(yinf2,ustead)

    real(kind=8), dimension(1) :: yinf2
    real(kind=8), dimension(1), intent(in) :: ustead

    yinf2(OPS_ACC1(0,0,0)) = yinf2(OPS_ACC1(0,0,0)) * ustead(OPS_ACC2(0,0,0))

END SUBROUTINE inflow_kernel_eqE

SUBROUTINE inflow_kernel_eqF(arr,sum_val)

    real(kind=8), dimension(1), intent(in) :: arr
    real(kind=8) :: sum_val

    sum_val = sum_val + arr(OPS_ACC1(0,0,0))

END SUBROUTINE inflow_kernel_eqF

SUBROUTINE inflow_kernel_eqG(arr,sum_umean,sum_denom)

    real(kind=8), dimension(1) :: arr
    real(kind=8), intent(in) :: sum_umean,sum_denom

    arr(OPS_ACC1(0,0,0)) = arr(OPS_ACC1(0,0,0)) - sum_umean/(sum_denom)

END SUBROUTINE inflow_kernel_eqG

SUBROUTINE inflow_kernel_eqH(ufilt,vfilt,wfilt,urms,phi)

    real(kind=8), dimension(1) :: ufilt,vfilt,wfilt
    real(kind=8), intent(in) :: urms,phi

    ufilt(OPS_ACC1(0,0,0)) = ufilt(OPS_ACC1(0,0,0)) * SQRT(urms**2.0_8*(1.0_8-phi**2.0_8))
    vfilt(OPS_ACC2(0,0,0)) = vfilt(OPS_ACC2(0,0,0)) * SQRT(urms**2.0_8*(1.0_8-phi**2.0_8))
    wfilt(OPS_ACC3(0,0,0)) = wfilt(OPS_ACC3(0,0,0)) * SQRT(urms**2.0_8*(1.0_8-phi**2.0_8))

END SUBROUTINE inflow_kernel_eqH

SUBROUTINE inflow_kernel_eqI(yfilt,yrms,yref_ispec,phi)

    real(kind=8), dimension(1) :: yfilt
    real(kind=8), intent(in) :: yrms,yref_ispec,phi

    yfilt(OPS_ACC1(0,0,0)) = yfilt(OPS_ACC1(0,0,0)) &
                           * SQRT((yrms*yref_ispec)**2.0_8*(1.0_8-phi**2.0_8))

END SUBROUTINE inflow_kernel_eqI
