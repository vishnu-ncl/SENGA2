SUBROUTINE set_zero_kernel(farray)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_dp

END SUBROUTINE set_zero_kernel

SUBROUTINE set_zero_kernel_multidim(farray, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: farray
    integer :: ispec

    farray(OPS_ACC_MD1(ispec,0,0,0)) = 0.0_dp

END SUBROUTINE set_zero_kernel_multidim
