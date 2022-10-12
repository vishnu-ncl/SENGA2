!eqone:   A = -half*(B+C*D+E*F+G*H) 
!eqtwo:   A = A+(B*C+D*E+F*G)
!eqthree: A = A+(B+C+D)*E
!eqfour:  A = A+var*B*C
!eqfive:  A = A+(B*C+D*E+F*G)*H
!eqsix:   A = A+(B+C+D)*E*F
!eqseven: A = A-var*B*C
!eqeight: A = A-B*C-D*E-F*G-H*I
!eqnine:  A = -half*(B+C)
!eqten:   A = A-B*C-D*E-F*G
!eqeleven: A = A-B*(C+D+E)
!eqtweleve: A = var1*B-var2*A
!eqthirteen: A=A-half*B*C*D


!eqone:   A = -half*(B+C*D+E*F+G*H)
SUBROUTINE complexmath_kernel_eqone(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7

    out_arr(OPS_ACC1(0,0,0)) = -0.5_dp * (in_arr1(OPS_ACC2(0,0,0)) &
                        + in_arr2(OPS_ACC3(0,0,0))*in_arr3(OPS_ACC4(0,0,0)) &
                        + in_arr4(OPS_ACC5(0,0,0))*in_arr5(OPS_ACC6(0,0,0)) &
                        + in_arr6(OPS_ACC7(0,0,0))*in_arr7(OPS_ACC8(0,0,0)))

END SUBROUTINE complexmath_kernel_eqone

!eqtwo:   A = A+(B*C+D*E+F*G)
SUBROUTINE complexmath_kernel_eqtwo(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                            +(in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                            + in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                            + in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)))

END SUBROUTINE complexmath_kernel_eqtwo

!eqthree: A = A+(B+C+D)*E
SUBROUTINE complexmath_kernel_eqthree(out_arr, in_arr1, in_arr2, in_arr3, in_arr4)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) &
                        + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))) * in_arr4(OPS_ACC5(0,0,0))

END SUBROUTINE complexmath_kernel_eqthree

!eqfour:  A = A+var*B*C
SUBROUTINE complexmath_kernel_eqfour(out_arr, in_arr1, in_arr2, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + var*in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE complexmath_kernel_eqfour

!eqfive:  A = A+(B*C+D*E+F*G)*H
SUBROUTINE complexmath_kernel_eqfive(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                  + (  in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) &
                    +  in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                    +  in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0))) &
                  * in_arr7(OPS_ACC8(0,0,0))

END SUBROUTINE complexmath_kernel_eqfive

!eqsix:   A = A+(B+C+D)*E*F
SUBROUTINE complexmath_kernel_eqsix(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                            + (in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0))) &
                            *in_arr4(OPS_ACC5(0,0,0))*in_arr5(OPS_ACC6(0,0,0))

END SUBROUTINE complexmath_kernel_eqsix

!eqseven: A = A-var*B*C
SUBROUTINE complexmath_kernel_eqseven(out_arr, in_arr1, in_arr2, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - var*in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE complexmath_kernel_eqseven

!eqeight: A = A-B*C-D*E-F*G-H*I
SUBROUTINE complexmath_kernel_eqeight(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                             - in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                             - in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0)) &
                             - in_arr7(OPS_ACC8(0,0,0)) * in_arr8(OPS_ACC9(0,0,0))

END SUBROUTINE complexmath_kernel_eqeight

!eqnine:  A = -half*(B+C)
SUBROUTINE complexmath_kernel_eqnine(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = -0.5_dp*(in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE complexmath_kernel_eqnine

!eqten: A = A-B*C-D*E-F*G
SUBROUTINE complexmath_kernel_eqten(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                             - in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                             - in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0))

END SUBROUTINE complexmath_kernel_eqten

!eqeleven A = A-B*(C+D+E)
SUBROUTINE complexmath_kernel_eqeleven(out_arr, in_arr1, in_arr2, in_arr3, in_arr4)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4
    
    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0)) &
                *(in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0))+in_arr4(OPS_ACC5(0,0,0)))

END SUBROUTINE complexmath_kernel_eqeleven

!eqtweleve: A = var1*B-var2*A
SUBROUTINE complexmath_kernel_eqtwelve(out_arr, in_arr1, var1, var2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: var1, var2

    out_arr(OPS_ACC1(0,0,0)) = var1*in_arr1(OPS_ACC2(0,0,0)) - var2*out_arr(OPS_ACC1(0,0,0))

END SUBROUTINE complexmath_kernel_eqtwelve

!eqthirteen: A = A-half*B*C*D
SUBROUTINE complexmath_kernel_eqthriteen(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - 0.5_dp*in_arr1(OPS_ACC2(0,0,0))&
                                                            *in_arr2(OPS_ACC3(0,0,0))*in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE complexmath_kernel_eqthriteen

SUBROUTINE complexmath_kernel_eqtau(out_arr1, out_arr2, in_arr1, in_arr2, in_arr3, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1, out_arr2
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    real(kind=dp), intent(in) :: var
    real(kind=dp) :: fornow

    fornow = in_arr1(OPS_ACC3(0,0,0))*var*in_arr2(OPS_ACC4(0,0,0))
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + fornow
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow*in_arr3(OPS_ACC5(0,0,0))

END SUBROUTINE complexmath_kernel_eqtau
