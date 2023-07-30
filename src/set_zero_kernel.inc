!set zero for x,y,z dimension
SUBROUTINE set_zero_kernel(farray)
    
    implicit none

    real(kind=8), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel

SUBROUTINE set_zero_kernel_MD5(farray, ispec)
    
    implicit none

    real(kind=8), dimension(1) :: farray
    integer(kind=4) :: ispec

    farray(OPS_ACC_MD1(ispec,0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_MD5

SUBROUTINE set_zero_kernel_MD6(farray, ispec)
    
    implicit none

    real(kind=8), dimension(1) :: farray
    integer(kind=4) :: ispec

    farray(OPS_ACC_MD1(ispec,0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_MD6

SUBROUTINE set_zero_kernel_int(farray)
    
    implicit none

    integer(kind=4), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0

END SUBROUTINE set_zero_kernel_int

!range x=1,1   y=1,nyglbl   z=1,nzglbl
!set zero for y,z dimension
SUBROUTINE set_zero_kernel_xdir(farray)
    
    implicit none

    real(kind=8), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_xdir

!range x=1,nxglbl   y=1,1   z=1,nzglbl
!set zero for x,z dimension
SUBROUTINE set_zero_kernel_ydir(farray)
    
    implicit none

    real(kind=8), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_ydir

!range x=1,nxglbl   y=1,nyglbl   z=1,1
!set zero for x,y dimension
SUBROUTINE set_zero_kernel_zdir(farray)
    
    implicit none

    real(kind=8), dimension(1) :: farray

    farray(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE set_zero_kernel_zdir
