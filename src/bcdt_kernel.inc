SUBROUTINE bcdt_kernel_xdir(strdx, dddtx, init_val)

    real(kind=8), dimension(1) :: strdx, dddtx
    real(kind=8), intent(in) :: init_val

    strdx(OPS_ACC1(0,0,0)) = init_val
    dddtx(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcdt_kernel_xdir

SUBROUTINE bcdt_kernel_ydir(strdy, dddty, init_val)

    real(kind=8), dimension(1) :: strdy, dddty
    real(kind=8), intent(in) :: init_val

    strdy(OPS_ACC1(0,0,0)) = init_val
    dddty(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcdt_kernel_ydir

SUBROUTINE bcdt_kernel_zdir(strdz, dddtz, init_val)

    real(kind=8), dimension(1) :: strdz, dddtz
    real(kind=8), intent(in) :: init_val

    strdz(OPS_ACC1(0,0,0)) = init_val
    dddtz(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcdt_kernel_zdir

SUBROUTINE bcdt_kernel_xdir_eqA(strdx, dddtx, init_val)

    real(kind=8), dimension(1) :: strdx, dddtx
    real(kind=8), dimension(nbcprr), intent(in) :: init_val

    strdx(OPS_ACC1(0,0,0)) = init_val(1)
    dddtx(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcdt_kernel_xdir_eqA

SUBROUTINE bcdt_kernel_ydir_eqA(strdy, dddty, init_val)

    real(kind=8), dimension(1) :: strdy, dddty
    real(kind=8), dimension(nbcprr), intent(in) :: init_val

    strdy(OPS_ACC1(0,0,0)) = init_val(1)
    dddty(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcdt_kernel_ydir_eqA

SUBROUTINE bcdt_kernel_zdir_eqA(strdz, dddtz, init_val)

    real(kind=8), dimension(1) :: strdz, dddtz
    real(kind=8), dimension(nbcprr), intent(in) :: init_val

    strdz(OPS_ACC1(0,0,0)) = init_val(1)
    dddtz(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcdt_kernel_zdir_eqA
