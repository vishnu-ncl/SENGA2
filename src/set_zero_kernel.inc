SUBROUTINE set_zero_kernel(farray)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: farray
    farray(OPS_ACC1(0,0,0)) = 0.0_dp

END SUBROUTINE set_zero_kernel
