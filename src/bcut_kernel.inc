SUBROUTINE bcut_kernel_xdir_const_uvel(strux, strvx, strwx, dudtx, dvdtx, dwdtx, init_val)
    
    

    real(kind=8), dimension(1) :: strux, strvx, strwx, dudtx, dvdtx, dwdtx
    real(kind=8), dimension(nbcprr), intent(in) :: init_val

    strux(OPS_ACC1(0,0,0)) = init_val(1)
    strvx(OPS_ACC2(0,0,0)) = 0.0_8
    strwx(OPS_ACC3(0,0,0)) = 0.0_8

    dudtx(OPS_ACC4(0,0,0)) = 0.0_8
    dvdtx(OPS_ACC5(0,0,0)) = 0.0_8
    dwdtx(OPS_ACC6(0,0,0)) = 0.0_8

END SUBROUTINE bcut_kernel_xdir_const_uvel

SUBROUTINE bcut_kernel_xdir_sinusoidal_uvel(strux, strvx, strwx, dudtx, dvdtx, dwdtx, init_val1, init_val2)
    
    

    real(kind=8), dimension(1) :: strux, strvx, strwx, dudtx, dvdtx, dwdtx
    real(kind=8), intent(in) :: init_val1, init_val2

    strux(OPS_ACC1(0,0,0)) = init_val1
    strvx(OPS_ACC2(0,0,0)) = 0.0_8
    strwx(OPS_ACC3(0,0,0)) = 0.0_8

    dudtx(OPS_ACC4(0,0,0)) = init_val2
    dvdtx(OPS_ACC5(0,0,0)) = 0.0_8
    dwdtx(OPS_ACC6(0,0,0)) = 0.0_8

END SUBROUTINE bcut_kernel_xdir_sinusoidal_uvel

SUBROUTINE bcut_kernel_ydir(struy, strvy, strwy, dudty, dvdty, dwdty, init_val)
    
    

    real(kind=8), dimension(1) :: struy, strvy, strwy, dudty, dvdty, dwdty
    real(kind=8), dimension(nbcprr), intent(in) :: init_val

    struy(OPS_ACC1(0,0,0)) = 0.0_8
    strvy(OPS_ACC2(0,0,0)) = init_val(1)
    strwy(OPS_ACC3(0,0,0)) = 0.0_8

    dudty(OPS_ACC4(0,0,0)) = 0.0_8
    dvdty(OPS_ACC5(0,0,0)) = 0.0_8
    dwdty(OPS_ACC6(0,0,0)) = 0.0_8

END SUBROUTINE bcut_kernel_ydir

SUBROUTINE bcut_kernel_zdir(struz, strvz, strwz, dudtz, dvdtz, dwdtz, init_val)
    
    

    real(kind=8), dimension(1) :: struz, strvz, strwz, dudtz, dvdtz, dwdtz
    real(kind=8), dimension(nbcprr), intent(in) :: init_val

    struz(OPS_ACC1(0,0,0)) = 0.0_8
    strvz(OPS_ACC2(0,0,0)) = 0.0_8
    strwz(OPS_ACC3(0,0,0)) = init_val(1)

    dudtz(OPS_ACC4(0,0,0)) = 0.0_8
    dvdtz(OPS_ACC5(0,0,0)) = 0.0_8
    dwdtz(OPS_ACC6(0,0,0)) = 0.0_8

END SUBROUTINE bcut_kernel_zdir

!A = A + B*val1
SUBROUTINE bcut_kernel_xdir_eqA(out_arr, in_arr, val1)
    
    

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    real(kind=8), intent(in) :: val1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr(OPS_ACC2(0,0,0))*val1

END SUBROUTINE bcut_kernel_xdir_eqA

!A = A - val1*B*val2
SUBROUTINE bcut_kernel_xdir_eqB(out_arr, in_arr, val1, val2)
    
    

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    real(kind=8), intent(in) :: val1, val2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - val1*in_arr(OPS_ACC2(0,0,0))*val2

END SUBROUTINE bcut_kernel_xdir_eqB

!A = A + val1*B*val2
SUBROUTINE bcut_kernel_xdir_eqC(out_arr, in_arr, val1, val2)
    
    

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    real(kind=8), intent(in) :: val1, val2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + val1*in_arr(OPS_ACC2(0,0,0))*val2

END SUBROUTINE bcut_kernel_xdir_eqC

!A = A * val1
SUBROUTINE bcut_kernel_xdir_eqD(out_arr, val1)
    
    

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: val1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * val1

END SUBROUTINE bcut_kernel_xdir_eqD

!A = A + val1
SUBROUTINE bcut_kernel_xdir_eqE(out_arr, val1)
    
    

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: val1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + val1

END SUBROUTINE bcut_kernel_xdir_eqE

