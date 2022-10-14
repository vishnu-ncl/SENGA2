!A = LOG(B)
SUBROUTINE math_kernel_eqA(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = LOG(in_arr1(OPS_ACC2(0,0,0)))

END SUBROUTINE math_kernel_eqA

!A = -B
SUBROUTINE math_kernel_eqB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = -in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqB

!A = A*var
SUBROUTINE math_kernel_eqC(out_arr, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * var

END SUBROUTINE math_kernel_eqC

!A = B*var
SUBROUTINE math_kernel_eqD(out_arr, in_arr1, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * var

END SUBROUTINE math_kernel_eqD

!A = B*C*var
SUBROUTINE math_kernel_eqE(out_arr, in_arr1, in_arr2, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * var

END SUBROUTINE math_kernel_eqE

!A = A+var*BC*C
SUBROUTINE math_kernel_eqF(out_arr, in_arr1, in_arr2, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + var*in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqF

!A = B-var
SUBROUTINE math_kernel_eqG(out_arr, in_arr1, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) - var

END SUBROUTINE math_kernel_eqG

!A = A-var*B*C
SUBROUTINE math_kernel_eqH(out_arr, in_arr1, in_arr2, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - var*in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqH

!A = var1*B-var2*A
SUBROUTINE math_kernel_eqI(out_arr, in_arr1, var1, var2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: var1, var2

    out_arr(OPS_ACC1(0,0,0)) = var1*in_arr1(OPS_ACC2(0,0,0)) - var2*out_arr(OPS_ACC1(0,0,0))

END SUBROUTINE math_kernel_eqI

!A = A+B
SUBROUTINE math_kernel_eqJ(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqJ

!A = B+C
SUBROUTINE math_kernel_eqK(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqK

!A = B+C+D
SUBROUTINE math_kernel_eqL(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_kernel_eqL

!A = A+B+C+D
SUBROUTINE math_kernel_eqM(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_kernel_eqM

!A = A+B*C
SUBROUTINE math_kernel_eqN(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE math_kernel_eqN

!A = A+B*C*C
SUBROUTINE math_kernel_eqO(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE math_kernel_eqO

!A = A+B*C*D
SUBROUTINE math_kernel_eqP(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_kernel_eqP

!A = B+C*D
SUBROUTINE math_kernel_eqQ(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + (in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_kernel_eqQ

!A = A*B+C*D
SUBROUTINE math_kernel_eqR(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = (out_arr(OPS_ACC1(0,0,0))*in_arr1(OPS_ACC2(0,0,0))) + (in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_kernel_eqR

!A = A/B
SUBROUTINE math_kernel_eqS(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) / in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqS

!A = B/C
SUBROUTINE math_kernel_eqT(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) / in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqT

!A = A*B
SUBROUTINE math_kernel_eqU(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqU

!A = B*C
SUBROUTINE math_kernel_eqV(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqV

!A = A-B
SUBROUTINE math_kernel_eqW(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqW

!A = A-B*C
SUBROUTINE math_kernel_eqX(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE math_kernel_eqX

!A = A+(B+C+D)*E
SUBROUTINE math_kernel_eqY(out_arr, in_arr1, in_arr2, in_arr3, in_arr4)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) &
                        + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))) * in_arr4(OPS_ACC5(0,0,0))

END SUBROUTINE math_kernel_eqY

!A = A+(B+C+D)*E*F
SUBROUTINE math_kernel_eqZ(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                            + (in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0))) &
                            *in_arr4(OPS_ACC5(0,0,0))*in_arr5(OPS_ACC6(0,0,0))

END SUBROUTINE math_kernel_eqZ

!A = A+(B*C+D*E+F*G)
SUBROUTINE math_kernel_eqAA(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                            +(in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                            + in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                            + in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)))

END SUBROUTINE math_kernel_eqAA

!A = A+(B*C+D*E+F*G)*H
SUBROUTINE math_kernel_eqAB(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                  + (  in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) &
                    +  in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                    +  in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0))) &
                  * in_arr7(OPS_ACC8(0,0,0))

END SUBROUTINE math_kernel_eqAB

!A = -half*(B+C)
SUBROUTINE math_kernel_eqAC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = -0.5_dp*(in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE math_kernel_eqAC

!A = -half*(B+C*D+E*F+G*H)
SUBROUTINE math_kernel_eqAD(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7

    out_arr(OPS_ACC1(0,0,0)) = -0.5_dp * (in_arr1(OPS_ACC2(0,0,0)) &
                        + in_arr2(OPS_ACC3(0,0,0))*in_arr3(OPS_ACC4(0,0,0)) &
                        + in_arr4(OPS_ACC5(0,0,0))*in_arr5(OPS_ACC6(0,0,0)) &
                        + in_arr6(OPS_ACC7(0,0,0))*in_arr7(OPS_ACC8(0,0,0)))

END SUBROUTINE math_kernel_eqAD

!A = A-half*B*C*D
SUBROUTINE math_kernel_eqAE(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - 0.5_dp*in_arr1(OPS_ACC2(0,0,0))&
                                                            *in_arr2(OPS_ACC3(0,0,0))*in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_kernel_eqAE

!A = A-B*(C+D+E)
SUBROUTINE math_kernel_eqAF(out_arr, in_arr1, in_arr2, in_arr3, in_arr4)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0)) &
                *(in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0))+in_arr4(OPS_ACC5(0,0,0)))

END SUBROUTINE math_kernel_eqAF

!A = A-B*C-D*E-F*G
SUBROUTINE math_kernel_eqAG(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                             - in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                             - in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0))

END SUBROUTINE math_kernel_eqAG

!A = A-B*C-D*E-F*G-H*I
SUBROUTINE math_kernel_eqAH(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                             - in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                             - in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0)) &
                             - in_arr7(OPS_ACC8(0,0,0)) * in_arr8(OPS_ACC9(0,0,0))

END SUBROUTINE math_kernel_eqAH

SUBROUTINE math_kernel_eqtau(out_arr1, out_arr2, in_arr1, in_arr2, in_arr3, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1, out_arr2
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    real(kind=dp), intent(in) :: var
    real(kind=dp) :: fornow

    fornow = in_arr1(OPS_ACC3(0,0,0))*var*in_arr2(OPS_ACC4(0,0,0))
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + fornow
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow*in_arr3(OPS_ACC5(0,0,0))

END SUBROUTINE math_kernel_eqtau
