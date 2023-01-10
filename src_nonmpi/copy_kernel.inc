SUBROUTINE copy_kernel(out_arr, in_arr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE copy_kernel
