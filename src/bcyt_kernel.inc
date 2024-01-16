SUBROUTINE bcyt_kernel_xdir(stryx,dydtx,yrin,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: stryx,dydtx
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: yrin

    stryx(OPS_ACC1(0,0,0)) = yrin(ispec)
    dydtx(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcyt_kernel_xdir

SUBROUTINE bcyt_kernel_ydir(stryy,dydty,yrin,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: stryy,dydty
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: yrin

    stryy(OPS_ACC1(0,0,0)) = yrin(ispec)
    dydty(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcyt_kernel_ydir

SUBROUTINE bcyt_kernel_zdir(stryz,dydtz,yrin,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: stryz,dydtz
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: yrin

    stryz(OPS_ACC1(0,0,0)) = yrin(ispec)
    dydtz(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcyt_kernel_zdir
