SUBROUTINE bcut_kernel_xdir_const_uvel(strux, strvx, strwx, dudtx, dvdtx, dwdtx, init_val)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strux, strvx, strwx, dudtx, dvdtx, dwdtx
    real(kind=dp), intent(in) :: init_val

    strux(OPS_ACC1(0,0,0)) = init_val
    strvx(OPS_ACC2(0,0,0)) = 0.0_dp
    strwx(OPS_ACC3(0,0,0)) = 0.0_dp

    dudtx(OPS_ACC4(0,0,0)) = 0.0_dp
    dvdtx(OPS_ACC5(0,0,0)) = 0.0_dp
    dwdtx(OPS_ACC6(0,0,0)) = 0.0_dp

END SUBROUTINE bcut_kernel_xdir_const_uvel

SUBROUTINE bcut_kernel_xdir_sinusoidal_uvel(strux, strvx, strwx, dudtx, dvdtx, dwdtx, init_val1, init_val2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strux, strvx, strwx, dudtx, dvdtx, dwdtx
    real(kind=dp), intent(in) :: init_val1, init_val2

    strux(OPS_ACC1(0,0,0)) = init_val1
    strvx(OPS_ACC2(0,0,0)) = 0.0_dp
    strwx(OPS_ACC3(0,0,0)) = 0.0_dp

    dudtx(OPS_ACC4(0,0,0)) = init_val2
    dvdtx(OPS_ACC5(0,0,0)) = 0.0_dp
    dwdtx(OPS_ACC6(0,0,0)) = 0.0_dp

END SUBROUTINE bcut_kernel_xdir_sinusoidal_uvel

SUBROUTINE bcut_kernel_ydir(struy, strvy, strwy, dudty, dvdty, dwdty, init_val)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: struy, strvy, strwy, dudty, dvdty, dwdty
    real(kind=dp), intent(in) :: init_val

    struy(OPS_ACC1(0,0,0)) = 0.0_dp
    strvy(OPS_ACC2(0,0,0)) = init_val
    strwy(OPS_ACC3(0,0,0)) = 0.0_dp

    dudty(OPS_ACC4(0,0,0)) = 0.0_dp
    dvdty(OPS_ACC5(0,0,0)) = 0.0_dp
    dwdty(OPS_ACC6(0,0,0)) = 0.0_dp

END SUBROUTINE bcut_kernel_ydir

SUBROUTINE bcut_kernel_zdir(struz, strvz, strwz, dudtz, dvdtz, dwdtz, init_val)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: struz, strvz, strwz, dudtz, dvdtz, dwdtz
    real(kind=dp), intent(in) :: init_val

    struz(OPS_ACC1(0,0,0)) = 0.0_dp
    strvz(OPS_ACC2(0,0,0)) = 0.0_dp
    strwz(OPS_ACC3(0,0,0)) = init_val

    dudtz(OPS_ACC4(0,0,0)) = 0.0_dp
    dvdtz(OPS_ACC5(0,0,0)) = 0.0_dp
    dwdtz(OPS_ACC6(0,0,0)) = 0.0_dp

END SUBROUTINE bcut_kernel_zdir
