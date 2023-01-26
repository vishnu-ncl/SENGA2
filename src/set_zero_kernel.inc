!set zero for x,y,z dimension
SUBROUTINE set_zero_kernel(farray)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel

SUBROUTINE set_zero_kernel_MD(farray, ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: farray
    integer :: ispec

    farray(OPS_ACC_MD1(ispec,0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_MD

SUBROUTINE set_zero_kernel_MD_int(farray, ispec)
    use data_types
    implicit none

    integer, dimension(1) :: farray
    integer :: ispec

    farray(OPS_ACC_MD1(ispec,0,0,0)) = 0

END SUBROUTINE set_zero_kernel_MD_int

SUBROUTINE set_zero_kernel_MD_xdir(farray, ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: farray
    integer :: ispec

    farray(OPS_ACC_MD1(ispec,0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_MD_xdir

SUBROUTINE set_zero_kernel_MD_ydir(farray, ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: farray
    integer :: ispec

    farray(OPS_ACC_MD1(ispec,0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_MD_ydir

SUBROUTINE set_zero_kernel_MD_zdir(farray, ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: farray
    integer :: ispec

    farray(OPS_ACC_MD1(ispec,0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_MD_zdir

!set zero for y,z dimension
SUBROUTINE set_zero_kernel_xdir(farray)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_xdir

!set zero for x,z dimension
SUBROUTINE set_zero_kernel_ydir(farray)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_ydir

!set zero for x,y dimension
SUBROUTINE set_zero_kernel_zdir(farray)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_zdir
