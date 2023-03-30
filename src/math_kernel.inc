!A = LOG(B)
SUBROUTINE math_kernel_eqA(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = log(in_arr1(OPS_ACC2(0,0,0)))

END SUBROUTINE math_kernel_eqA

!A = -B
SUBROUTINE math_kernel_eqB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = -in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqB

!A = A*var
SUBROUTINE math_kernel_eqC(out_arr, var)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * var

END SUBROUTINE math_kernel_eqC

!A = B*var
SUBROUTINE math_kernel_eqD(out_arr, in_arr1, var, nspcmx, ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec,nspcmx
    real(kind=8), dimension(nspcmx), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * var(ispec)

END SUBROUTINE math_kernel_eqD

!A = B*C*var
SUBROUTINE math_kernel_eqE(out_arr, in_arr1, in_arr2, var)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=8), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = var * in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqE

!A = A+var*B*C
SUBROUTINE math_kernel_eqF(out_arr, in_arr1, in_arr2, rgspec, nspcmx, ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec, nspcmx
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + rgspec(ispec)*in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqF

!A = B-var
SUBROUTINE math_kernel_eqG(out_arr, in_arr1, var)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) - var

END SUBROUTINE math_kernel_eqG

!A = A-var*B*C
SUBROUTINE math_kernel_eqH(out_arr, in_arr1, in_arr2, rgspec, nspcmx, ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec, nspcmx
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - rgspec(ispec)*in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqH

!A = var1*B-var2*A
SUBROUTINE math_kernel_eqI(out_arr, in_arr1, var1, var2)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: var1, var2

    out_arr(OPS_ACC1(0,0,0)) = var1*in_arr1(OPS_ACC2(0,0,0)) - var2*out_arr(OPS_ACC1(0,0,0))

END SUBROUTINE math_kernel_eqI

!A = A+B
SUBROUTINE math_kernel_eqJ(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqJ

!A = B+C
SUBROUTINE math_kernel_eqK(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqK

!A = B+C+D
SUBROUTINE math_kernel_eqL(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_kernel_eqL

!A = A+B+C+D
SUBROUTINE math_kernel_eqM(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_kernel_eqM

!A = A+B*C
SUBROUTINE math_kernel_eqN(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE math_kernel_eqN

!A = A+B*C*C
SUBROUTINE math_kernel_eqO(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE math_kernel_eqO

!A = A+B*C*D
SUBROUTINE math_kernel_eqP(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_kernel_eqP

!A = B+C*D
SUBROUTINE math_kernel_eqQ(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + (in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_kernel_eqQ

!A = A*B+C*D
SUBROUTINE math_kernel_eqR(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = (out_arr(OPS_ACC1(0,0,0))*in_arr1(OPS_ACC2(0,0,0))) + (in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_kernel_eqR

!A = A/B
SUBROUTINE math_kernel_eqS(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) / in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqS

!A = B/C
SUBROUTINE math_kernel_eqT(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) / in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqT

!A = A*B
SUBROUTINE math_kernel_eqU(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqU

!A = B*C
SUBROUTINE math_kernel_eqV(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_kernel_eqV

!A = A-B
SUBROUTINE math_kernel_eqW(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_kernel_eqW

!A = A-B*C
SUBROUTINE math_kernel_eqX(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE math_kernel_eqX

!A = A+(B+C+D)*E
SUBROUTINE math_kernel_eqY(out_arr, in_arr1, in_arr2, in_arr3, in_arr4)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) &
                        + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))) * in_arr4(OPS_ACC5(0,0,0))

END SUBROUTINE math_kernel_eqY

!A = A+(B+C+D)*E*F
SUBROUTINE math_kernel_eqZ(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                            + (in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0))) &
                            *in_arr4(OPS_ACC5(0,0,0))*in_arr5(OPS_ACC6(0,0,0))

END SUBROUTINE math_kernel_eqZ

!A = A+(B*C+D*E+F*G)
SUBROUTINE math_kernel_eqAA(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                            +(in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                            + in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                            + in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)))

END SUBROUTINE math_kernel_eqAA

!A = A+(B*C+D*E+F*G)*H
SUBROUTINE math_kernel_eqAB(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7

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

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = -0.5_8*(in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE math_kernel_eqAC

!A = -half*(B+C*D+E*F+G*H)
SUBROUTINE math_kernel_eqAD(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7

    out_arr(OPS_ACC1(0,0,0)) = -0.5_8 * (in_arr1(OPS_ACC2(0,0,0)) &
                        + in_arr2(OPS_ACC3(0,0,0))*in_arr3(OPS_ACC4(0,0,0)) &
                        + in_arr4(OPS_ACC5(0,0,0))*in_arr5(OPS_ACC6(0,0,0)) &
                        + in_arr6(OPS_ACC7(0,0,0))*in_arr7(OPS_ACC8(0,0,0)))

END SUBROUTINE math_kernel_eqAD

!A = A-half*B*C*D
SUBROUTINE math_kernel_eqAE(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - 0.5_8*in_arr1(OPS_ACC2(0,0,0))&
                                                            *in_arr2(OPS_ACC3(0,0,0))*in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_kernel_eqAE

!A = A-B*(C+D+E)
SUBROUTINE math_kernel_eqAF(out_arr, in_arr1, in_arr2, in_arr3, in_arr4)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0)) &
                *(in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0))+in_arr4(OPS_ACC5(0,0,0)))

END SUBROUTINE math_kernel_eqAF

!A = A-B*C-D*E-F*G
SUBROUTINE math_kernel_eqAG(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                             - in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                             - in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0))

END SUBROUTINE math_kernel_eqAG

!A = A-B*C-D*E-F*G-H*I
SUBROUTINE math_kernel_eqAH(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                             - in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                             - in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0)) &
                             - in_arr7(OPS_ACC8(0,0,0)) * in_arr8(OPS_ACC9(0,0,0))

END SUBROUTINE math_kernel_eqAH

!A = exp(A)
SUBROUTINE math_kernel_eqAL(out_arr)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8) :: fornow

    fornow = out_arr(OPS_ACC1(0,0,0))
    out_arr(OPS_ACC1(0,0,0)) = EXP(fornow)

END SUBROUTINE math_kernel_eqAL

!A = A + var
SUBROUTINE math_kernel_eqAM(out_arr,var)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: var
    
    out_arr(OPS_ACC1(0,0,0)) = var + out_arr(OPS_ACC1(0,0,0))

END SUBROUTINE math_kernel_eqAM

!A = var
SUBROUTINE math_kernel_eqAN(out_arr,var)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = var

END SUBROUTINE math_kernel_eqAN

!A = A - B + half*(C*C+D*D+E*E)
SUBROUTINE math_kernel_eqAO(out_arr,in_arr1,in_arr2,in_arr3,in_arr4)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4
    real(kind=8) :: fornow

    fornow = in_arr2(OPS_ACC3(0,0,0))*in_arr2(OPS_ACC3(0,0,0))  & 
           + in_arr3(OPS_ACC4(0,0,0))*in_arr3(OPS_ACC4(0,0,0))  &
           + in_arr4(OPS_ACC5(0,0,0))*in_arr4(OPS_ACC5(0,0,0))
    
    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0)) + 0.5_8*fornow

END SUBROUTINE math_kernel_eqAO

SUBROUTINE math_kernel_eqAP(out_arr1,out_arr2,in_arr1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1, out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8) :: fornow
    
    fornow = alamda*EXP(rlamda*LOG(in_arr1(OPS_ACC3(0,0,0))))

    out_arr1(OPS_ACC1(0,0,0)) = fornow*out_arr2(OPS_ACC2(0,0,0))
    out_arr2(OPS_ACC2(0,0,0)) = fornow

END SUBROUTINE math_kernel_eqAP

SUBROUTINE math_kernel_eqAQ(out_arr1,out_arr2,in_arr1,racnst,rncnst,reovrr)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1, out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: racnst,rncnst,reovrr

    out_arr2(OPS_ACC2(0,0,0)) = racnst + rncnst*LOG(in_arr1(OPS_ACC3(0,0,0)))  &
            - reovrr/in_arr1(OPS_ACC3(0,0,0))
    out_arr1(OPS_ACC1(0,0,0)) = EXP(out_arr2(OPS_ACC2(0,0,0)))

END SUBROUTINE math_kernel_eqAQ

SUBROUTINE math_kernel_eqAR(out_arr1,out_arr2,out_arr3,in_arr1,in_arr2,racnst,rncnst,reovrr,flcnst)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1, out_arr2, out_arr3
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=8), intent(in) :: racnst,rncnst,reovrr,flcnst

    real(kind=8) :: preduc,fornow

!   EVALUATE K0
    out_arr3(OPS_ACC3(0,0,0)) = racnst + rncnst*LOG(in_arr1(OPS_ACC4(0,0,0)))  &
                              - reovrr/in_arr1(OPS_ACC4(0,0,0))
    out_arr3(OPS_ACC3(0,0,0)) = EXP(out_arr3(OPS_ACC3(0,0,0)))

!   EVALUATE REDUCED PRESURE
    preduc = in_arr2(OPS_ACC5(0,0,0))*out_arr3(OPS_ACC3(0,0,0)) /out_arr1(OPS_ACC1(0,0,0))

!   EVALUATE UPDATED FORWARD RATE CONSTANT
    fornow = flcnst*preduc/(1.0_8+preduc)
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0))*fornow

!   RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
    out_arr2(OPS_ACC2(0,0,0)) = LOG(out_arr1(OPS_ACC1(0,0,0)))

END SUBROUTINE math_kernel_eqAR

SUBROUTINE math_kernel_eqAS(out_arr1,out_arr2,out_arr3,in_arr1,in_arr2,racnst,rncnst,reovrr,talpha,ovtst1,tstar2,ovtst3,cfcst1,cfcst2,encst1,encst2,dtcnst,omalph,clnten)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1, out_arr2, out_arr3
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=8), intent(in) :: racnst,rncnst,reovrr,talpha,ovtst1,tstar2,ovtst3,cfcst1,cfcst2,encst1,encst2,dtcnst,omalph,clnten

    real(kind=8) :: preduc,fornow,trats1,trats2,trats3,ftcent,cfactr,enfact,fbroad

!   EVALUATE K0
    out_arr3(OPS_ACC3(0,0,0)) = racnst + rncnst*LOG(in_arr1(OPS_ACC4(0,0,0)))  &
                              - reovrr/in_arr1(OPS_ACC4(0,0,0))
    out_arr3(OPS_ACC3(0,0,0)) = EXP(out_arr3(OPS_ACC3(0,0,0)))

!   EVALUATE REDUCED PRESURE
    preduc = in_arr2(OPS_ACC5(0,0,0))*out_arr3(OPS_ACC3(0,0,0)) /out_arr1(OPS_ACC1(0,0,0))

!   EVALUATE FCENT
    trats1 = in_arr1(OPS_ACC4(0,0,0))*ovtst1
    trats2 = tstar2/in_arr1(OPS_ACC4(0,0,0))
    trats3 = in_arr1(OPS_ACC4(0,0,0))*ovtst3
    ftcent = omalph*EXP(trats3) + talpha*EXP(trats1)  &
           + EXP(trats2)
    ftcent = LOG10(ftcent)

!   EVALUATE BROADENING FACTOR
    cfactr = cfcst1 + cfcst2*ftcent
    enfact = encst1 + encst2*ftcent
    fornow = LOG10(preduc) + cfactr
    fornow = fornow/(enfact - dtcnst*fornow)
    fornow = 1.0_8 + fornow*fornow
    fbroad = ftcent/fornow
    fbroad = EXP(fbroad*clnten)

!   EVALUATE UPDATED FORWARD RATE CONSTANT
    fornow = fbroad*preduc/(1.0_8+preduc)
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0))*fornow

!   RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
    out_arr2(OPS_ACC2(0,0,0)) = LOG(out_arr1(OPS_ACC1(0,0,0)))

END SUBROUTINE math_kernel_eqAS

SUBROUTINE math_kernel_eqAT(out_arr1,out_arr2,out_arr3,in_arr1,in_arr2,racnst,rncnst,reovrr,acfsri,bcfsrm,ovcsrm,dcfsri,ecfsri)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1, out_arr2, out_arr3
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=8), intent(in) :: racnst,rncnst,reovrr,acfsri,bcfsrm,ovcsrm,dcfsri,ecfsri

    real(kind=8) :: preduc,fornow,fbroad

!   EVALUATE K0
    out_arr3(OPS_ACC3(0,0,0)) = racnst + rncnst*LOG(in_arr1(OPS_ACC4(0,0,0)))  &
                              - reovrr/in_arr1(OPS_ACC4(0,0,0))
    out_arr3(OPS_ACC3(0,0,0)) = EXP(out_arr3(OPS_ACC3(0,0,0)))

!   EVALUATE REDUCED PRESURE
    preduc = in_arr2(OPS_ACC5(0,0,0))*out_arr3(OPS_ACC3(0,0,0)) /out_arr1(OPS_ACC1(0,0,0))

!   EVALUATE FCENT
    fornow = bcfsrm/in_arr1(OPS_ACC4(0,0,0))
    fbroad = acfsri*EXP(fornow)
    fornow = in_arr1(OPS_ACC4(0,0,0))*ovcsrm
    fbroad = fbroad + EXP(fornow)
    fornow = LOG10(preduc)
    fornow = 1.0_8/(1.0_8+LOG10(fornow))
    fbroad = dcfsri*EXP(fornow*LOG(fbroad))
    fornow = fbroad*EXP(ecfsri*LOG(in_arr1(OPS_ACC4(0,0,0))))

!   EVALUATE UPDATED FORWARD RATE CONSTANT
    fornow = fbroad*preduc/(1.0_8+preduc)
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0))*fornow

!   RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
    out_arr2(OPS_ACC2(0,0,0)) = LOG(out_arr1(OPS_ACC1(0,0,0)))

END SUBROUTINE math_kernel_eqAT

SUBROUTINE math_kernel_eqtau(out_arr1,out_arr2,in_arr1,in_arr2,in_arr3,var)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1, out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    real(kind=8), intent(in) :: var
    real(kind=8) :: fornow

    fornow = in_arr1(OPS_ACC3(0,0,0))*var*in_arr2(OPS_ACC4(0,0,0))
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + fornow
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow*in_arr3(OPS_ACC5(0,0,0))

END SUBROUTINE math_kernel_eqtau

SUBROUTINE math_kernel_getval(prun,prefer)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: prun
    real(kind=8) :: prefer

    prefer = prun(OPS_ACC1(0,0,0))

END SUBROUTINE math_kernel_getval

SUBROUTINE math_kernel_print_drhs(out_arr, itime)
    implicit none

    real(kind=8), dimension(1), intent(in) :: out_arr
    integer, intent(in) :: itime

    WRITE(*,'(a,I7,a,F16.8)')  &
      'test_drhs: (step=',itime,') value: ',out_arr(OPS_ACC1(0,0,0))
!    print *, 'test_drhs: (step=',itime,') value: ',out_arr(OPS_ACC1(0,0,0))

END SUBROUTINE math_kernel_print_drhs

SUBROUTINE math_kernel_print_erhs(out_arr, itime)
    implicit none

    real(kind=8), dimension(1), intent(in) :: out_arr
    integer, intent(in) :: itime

    WRITE(*,'(a,I7,a,F16.8)')  &
      'test_erhs: (step=',itime,') value: ',out_arr(OPS_ACC1(0,0,0))
!    print *, 'test_erhs: (step=',itime,') value: ',out_arr(OPS_ACC1(0,0,0))

END SUBROUTINE math_kernel_print_erhs

SUBROUTINE math_kernel_print_urhs(out_arr, itime)
    implicit none

    real(kind=8), dimension(1), intent(in) :: out_arr
    integer, intent(in) :: itime

    WRITE(*,'(a,I7,a,F16.8)')  &
      'test_urhs: (step=',itime,') value: ',out_arr(OPS_ACC1(0,0,0))
!    print *, 'test_urhs: (step=',itime,') value: ',out_arr(OPS_ACC1(0,0,0))

END SUBROUTINE math_kernel_print_urhs
