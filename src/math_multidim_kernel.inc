!A = B_multidim
SUBROUTINE math_multidim_kernel_eqA(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_multidim_kernel_eqA

!A_multidim = B_multidim
SUBROUTINE math_multidim_kernel_eqB(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_multidim_kernel_eqB

!A = B*C_multidim
SUBROUTINE math_multidim_kernel_eqC(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    integer :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))

END SUBROUTINE math_multidim_kernel_eqC

!A = B_multidim*C
SUBROUTINE math_multidim_kernel_eqD(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    integer :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_multidim_kernel_eqD

!A_multidim = A_multidim/B
SUBROUTINE math_multidim_kernel_eqE(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) / in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_multidim_kernel_eqE

!A = B*C_multidim*D
SUBROUTINE math_multidim_kernel_eqF(out_arr,in_arr1,in_arr2,in_arr3,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))*in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_multidim_kernel_eqF

!A = B*C_multidim*(one-val1/D)
SUBROUTINE math_multidim_kernel_eqG(out_arr,in_arr1,in_arr2,in_arr3,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    real(kind=dp), intent(in) :: val1
    integer :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))*(1.0_dp - val1/in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_multidim_kernel_eqG

!A_multidim = A_multidim - half*(B+C+D)
SUBROUTINE math_multidim_kernel_eqH(out_arr,in_arr1,in_arr2,in_arr3,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_dp*(in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_multidim_kernel_eqH

!A_multidim = A_multidim - half*B_multidim*C
SUBROUTINE math_multidim_kernel_eqI(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    integer :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_dp*in_arr1(OPS_ACC_MD2(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_multidim_kernel_eqI

!A_multidim = A_multidim - half*(B*C+D*E+F*G)
SUBROUTINE math_multidim_kernel_eqJ(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6
    integer :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_dp*( in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                                                   + in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                                                   + in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) )

END SUBROUTINE math_multidim_kernel_eqJ

!A_multidim = B_multidim - A_multidim*C - D*E - F*G - H*I
SUBROUTINE math_multidim_kernel_eqK(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,in_arr7,in_arr8,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8
    integer :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0)) - out_arr(OPS_ACC_MD1(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                                       - in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                                       - in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) &
                                                                       - in_arr7(OPS_ACC8(0,0,0))*in_arr8(OPS_ACC9(0,0,0))

END SUBROUTINE math_multidim_kernel_eqK