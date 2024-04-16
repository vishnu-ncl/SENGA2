SUBROUTINE set_zero_kernel(out_arr)

    real(kind=8), dimension(1) :: out_arr

    out_arr(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE set_zero
