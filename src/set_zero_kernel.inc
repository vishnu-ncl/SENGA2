!set zero for x,y,z dimension
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

!set zero for y,z dimension
SUBROUTINE set_zero_kernel_xdir(farray)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_dp

END SUBROUTINE set_zero_kernel_xdir

!set zero for x,z dimension
SUBROUTINE set_zero_kernel_ydir(farray)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_dp

END SUBROUTINE set_zero_kernel_ydir

!set zero for x,y dimension
SUBROUTINE set_zero_kernel_zdir(farray)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_dp

END SUBROUTINE set_zero_kernel_zdir
