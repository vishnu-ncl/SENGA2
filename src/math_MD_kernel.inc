!A = B_multidim
SUBROUTINE math_MD_kernel_eqA(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqA

!A_multidim = B_multidim
SUBROUTINE math_MD_kernel_eqB(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqB

!A = B*C_multidim
SUBROUTINE math_MD_kernel_eqC(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqC

!A = B_multidim*C
SUBROUTINE math_MD_kernel_eqD(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_MD_kernel_eqD

!A_multidim = A_multidim/B
SUBROUTINE math_MD_kernel_eqE(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) / in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_MD_kernel_eqE

!A = B*C_multidim*D
SUBROUTINE math_MD_kernel_eqF(out_arr,in_arr1,in_arr2,in_arr3,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))*in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_MD_kernel_eqF

!A = B*C_multidim*(one-val1/D)
SUBROUTINE math_MD_kernel_eqG(out_arr,in_arr1,in_arr2,in_arr3,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))*(1.0_dp - val1/in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_MD_kernel_eqG

!A_multidim = A_multidim - half*(B+C+D)
SUBROUTINE math_MD_kernel_eqH(out_arr,in_arr1,in_arr2,in_arr3,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_dp*(in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_MD_kernel_eqH

!A_multidim = A_multidim - half*B_multidim*C
SUBROUTINE math_MD_kernel_eqI(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_dp*in_arr1(OPS_ACC_MD2(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_MD_kernel_eqI

!A_multidim = A_multidim - half*(B*C+D*E+F*G)
SUBROUTINE math_MD_kernel_eqJ(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_dp*( in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                                                   + in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                                                   + in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) )

END SUBROUTINE math_MD_kernel_eqJ

!A_multidim = B_multidim - A_multidim*C - D*E - F*G - H*I
SUBROUTINE math_MD_kernel_eqK(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,in_arr7,in_arr8,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0)) - out_arr(OPS_ACC_MD1(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                                       - in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                                       - in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) &
                                                                       - in_arr7(OPS_ACC8(0,0,0))*in_arr8(OPS_ACC9(0,0,0))

END SUBROUTINE math_MD_kernel_eqK

SUBROUTINE math_MD_kernel_eqL(out_arr1,out_arr2,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1,out_arr2
    real(kind=dp), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = in_arr1(OPS_ACC3(0,0,0))*in_arr2(OPS_ACC4(0,0,0))  &
            + in_arr3(OPS_ACC5(0,0,0))*in_arr4(OPS_ACC6(0,0,0))  &
            + in_arr5(OPS_ACC7(0,0,0))*in_arr6(OPS_ACC8(0,0,0))

!   Y EQUATION
    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = out_arr1(OPS_ACC_MD1(ispec,0,0,0)) + fornow

!   DIFFUSION CORRECTION VELOCITY DIVERGENCE
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow

END SUBROUTINE math_MD_kernel_eqL

SUBROUTINE math_MD_kernel_eqM(out_arr1,out_arr2,in_arr1,in_arr2,in_arr3,in_arr4,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1,out_arr2
    real(kind=dp), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = ( in_arr1(OPS_ACC3(0,0,0)) + in_arr2(OPS_ACC4(0,0,0)) + in_arr3(OPS_ACC5(0,0,0)) ) &
            * in_arr4(OPS_ACC6(0,0,0))

!   Y EQUATION
    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = out_arr1(OPS_ACC_MD1(ispec,0,0,0)) + fornow

!   DIFFUSION CORRECTION VELOCITY DIVERGENCE
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow

END SUBROUTINE math_MD_kernel_eqM

!A_multidim_lastIdx = A_multidim_lastIdx + A_multidim_loopIdx
SUBROUTINE math_MD_kernel_eqN(out_arr1,ispec,nspec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    integer, intent(in) :: ispec,nspec

    out_arr1(OPS_ACC_MD1(nspec,0,0,0)) = out_arr1(OPS_ACC_MD1(nspec,0,0,0)) + out_arr1(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqN

!A_multidim = 1.0 - A_multidim
SUBROUTINE math_MD_kernel_eqO(out_arr1,nspec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    integer, intent(in) :: nspec

    out_arr1(OPS_ACC_MD1(nspec,0,0,0)) = 1.0_dp - out_arr1(OPS_ACC_MD1(nspec,0,0,0))

END SUBROUTINE math_MD_kernel_eqO

!A = A + var*B_multidim
SUBROUTINE math_MD_kernel_eqP(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec
    
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + val1*in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqP

!A_multidim = var1 + B*(var2-var1)
SUBROUTINE math_MD_kernel_eqQ(out_arr1,in_arr1,val1,val2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: val1,val2
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = val1 + in_arr1(OPS_ACC2(0,0,0))*(val2-val1)
    
END SUBROUTINE math_MD_kernel_eqQ

!A = A*max(B_multidim*var    zero)
SUBROUTINE math_MD_kernel_eqR(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = max(in_arr1(OPS_ACC_MD2(ispec,0,0,0))*val1, 0.0_dp)
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) * fornow

END SUBROUTINE math_MD_kernel_eqR
    
!A_multidim = A_multidim + B*var
SUBROUTINE math_MD_kernel_eqS(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = out_arr1(OPS_ACC_MD1(ispec,0,0,0)) + in_arr1(OPS_ACC2(0,0,0))*val1

END SUBROUTINE math_MD_kernel_eqS

!A_multidim = var
SUBROUTINE math_MD_kernel_eqT(out_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = val1

END SUBROUTINE math_MD_kernel_eqT

!A_multidim = B*A_multidim
SUBROUTINE math_MD_kernel_eqU(out_arr1,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*out_arr1(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqU

