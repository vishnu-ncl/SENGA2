SUBROUTINE bcdt_kernel_xdir(strdx, dddtx, init_val)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strdx, dddtx
    real(kind=dp), intent(in) :: init_val

    strdx(OPS_ACC1(0,0,0)) = init_val
    dddtx(OPS_ACC2(0,0,0)) = 0.0_dp

END SUBROUTINE bcdt_kernel_xdir

SUBROUTINE bcdt_kernel_ydir(strdy, dddty, init_val)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strdy, dddty
    real(kind=dp), intent(in) :: init_val

    strdy(OPS_ACC1(0,0,0)) = init_val
    dddty(OPS_ACC2(0,0,0)) = 0.0_dp

END SUBROUTINE bcdt_kernel_ydir

SUBROUTINE bcdt_kernel_zdir(strdz, dddtz, init_val)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strdz, dddtz
    real(kind=dp), intent(in) :: init_val

    strdz(OPS_ACC1(0,0,0)) = init_val
    dddtz(OPS_ACC2(0,0,0)) = 0.0_dp

END SUBROUTINE bcdt_kernel_zdir
