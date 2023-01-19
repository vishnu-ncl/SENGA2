SUBROUTINE bcyt_kernel_xdir(stryx, dydtx, ispec, yrin)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: stryx, dydtx
    integer, intent(in) :: ispec
    real(kind=8), intent(in) :: yrin

    stryx(OPS_ACC_MD1(ispec,0,0,0)) = yrin
    dydtx(OPS_ACC_MD2(ispec,0,0,0)) = 0.0_8

END SUBROUTINE bcyt_kernel_xdir

SUBROUTINE bcyt_kernel_ydir(stryy, dydty, ispec, yrin)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: stryy, dydty
    integer, intent(in) :: ispec
    real(kind=8), intent(in) :: yrin

    stryy(OPS_ACC_MD1(ispec,0,0,0)) = yrin
    dydty(OPS_ACC_MD2(ispec,0,0,0)) = 0.0_8

END SUBROUTINE bcyt_kernel_ydir

SUBROUTINE bcyt_kernel_zdir(stryz, dydtz, ispec, yrin)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: stryz, dydtz
    integer, intent(in) :: ispec
    real(kind=8), intent(in) :: yrin

    stryz(OPS_ACC_MD1(ispec,0,0,0)) = yrin
    dydtz(OPS_ACC_MD2(ispec,0,0,0)) = 0.0_8

END SUBROUTINE bcyt_kernel_zdir
