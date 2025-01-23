!A = LOG(B)
SUBROUTINE maths_kernel_eqA(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = LOG(in_arr(OPS_ACC2(0,0,0)))

END SUBROUTINE maths_kernel_eqA

!A = exp(A)
SUBROUTINE maths_kernel_eqB(out_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8) :: fornow

    fornow = out_arr(OPS_ACC1(0,0,0))
    out_arr(OPS_ACC1(0,0,0)) = EXP(fornow)

END SUBROUTINE maths_kernel_eqB

!A = -B
SUBROUTINE maths_kernel_eqC(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = -in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE maths_kernel_eqC

!A = var
SUBROUTINE maths_kernel_eqD(out_arr,var)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = var

END SUBROUTINE maths_kernel_eqD

!A = var(indx)
SUBROUTINE maths_kernel_eqE(out_arr,var,ispec)

    real(kind=8), dimension(1) :: out_arr
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = var(ispec)

END SUBROUTINE maths_kernel_eqE

!A = A + var
SUBROUTINE maths_kernel_eqF(out_arr,var)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = var + out_arr(OPS_ACC1(0,0,0))

END SUBROUTINE maths_kernel_eqF

!A = A*var
SUBROUTINE maths_kernel_eqG(out_arr, var)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * var

END SUBROUTINE maths_kernel_eqG

!A = B*var(indx)
SUBROUTINE maths_kernel_eqH(out_arr, in_arr, var, ispec)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC2(0,0,0)) * var(ispec)

END SUBROUTINE maths_kernel_eqH

!A = var*B*C
SUBROUTINE maths_kernel_eqI(out_arr, in_arr1, in_arr2, var)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=8), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = var * in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE maths_kernel_eqI

!A = A + var*B
SUBROUTINE maths_kernel_eqJ(out_arr,in_arr,var)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    real(kind=8), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + var*in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE maths_kernel_eqJ

!A = A + var(indx)*B
SUBROUTINE maths_kernel_eqK(out_arr,in_arr,rgspec,ispec)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + rgspec(ispec)*in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE maths_kernel_eqK

!A = A+var(indx)*B*C
SUBROUTINE maths_kernel_eqL(out_arr, in_arr1, in_arr2, rgspec, ispec)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + rgspec(ispec)*in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE maths_kernel_eqL

!A = A-var(indx)*B*C
SUBROUTINE maths_kernel_eqM(out_arr, in_arr1, in_arr2, rgspec, ispec)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - rgspec(ispec)*in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE maths_kernel_eqM

!A = var1*B-var2*A
SUBROUTINE maths_kernel_eqN(out_arr, in_arr1, var1, var2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: var1, var2

    out_arr(OPS_ACC1(0,0,0)) = var1*in_arr1(OPS_ACC2(0,0,0)) - var2*out_arr(OPS_ACC1(0,0,0))

END SUBROUTINE maths_kernel_eqN

!A = B*C*(one-var(indx)/D)
SUBROUTINE maths_kernel_eqO(out_arr,in_arr1,in_arr2,in_arr3,var,ispec)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))*(1.0_8 - var(ispec)/in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE maths_kernel_eqO

!A = A*max(B*var(indx)    zero)
SUBROUTINE maths_kernel_eqP(out_arr,in_arr,var,ispec)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: var
    real(kind=8) :: fornow

    fornow = MAX(in_arr(OPS_ACC2(0,0,0))*var(ispec), 0.0_8)
    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * fornow

END SUBROUTINE maths_kernel_eqP

!A = A+B
SUBROUTINE maths_kernel_eqQ(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE maths_kernel_eqQ

!A = B+C
SUBROUTINE maths_kernel_eqR(out_arr, in_arr1, in_arr2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE maths_kernel_eqR

!A = A-B
SUBROUTINE maths_kernel_eqS(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE maths_kernel_eqS

!A = A/B
SUBROUTINE maths_kernel_eqT(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) / in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE maths_kernel_eqT

!A = B/C
SUBROUTINE maths_kernel_eqU(out_arr, in_arr1, in_arr2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) / in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE maths_kernel_eqU

SUBROUTINE maths_kernel_eqU_fused(out_arr1, out_arr2, out_arr3, in_arr1, in_arr2, in_arr3, in_arr4)

    real(kind=8), dimension(1) :: out_arr1, out_arr2, out_arr3
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4

    out_arr1(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC4(0,0,0)) / in_arr4(OPS_ACC7(0,0,0))
    out_arr2(OPS_ACC2(0,0,0)) = in_arr2(OPS_ACC5(0,0,0)) / in_arr4(OPS_ACC7(0,0,0))
    out_arr3(OPS_ACC3(0,0,0)) = in_arr3(OPS_ACC6(0,0,0)) / in_arr4(OPS_ACC7(0,0,0))

END SUBROUTINE maths_kernel_eqU_fused

!A = A*B
SUBROUTINE maths_kernel_eqV(out_arr, in_arr)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * in_arr(OPS_ACC2(0,0,0))

END SUBROUTINE maths_kernel_eqV

SUBROUTINE maths_kernel_eqV_fused(out_arr1, out_arr2, out_arr3, out_arr4, in_arr)

    real(kind=8), dimension(1) :: out_arr1, out_arr2, out_arr3, out_arr4
    real(kind=8), dimension(1), intent(in) :: in_arr

    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) * in_arr(OPS_ACC5(0,0,0))
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) * in_arr(OPS_ACC5(0,0,0))
    out_arr3(OPS_ACC3(0,0,0)) = out_arr3(OPS_ACC3(0,0,0)) * in_arr(OPS_ACC5(0,0,0))

    out_arr4(OPS_ACC4(0,0,0)) = out_arr4(OPS_ACC4(0,0,0)) * in_arr(OPS_ACC5(0,0,0))

END SUBROUTINE maths_kernel_eqV_fused

!A = B*C
SUBROUTINE maths_kernel_eqW(out_arr, in_arr1, in_arr2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE maths_kernel_eqW

!A = B+C+D
SUBROUTINE maths_kernel_eqX(out_arr, in_arr1, in_arr2, in_arr3)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE maths_kernel_eqX

!A = B*C*D
SUBROUTINE maths_kernel_eqY(out_arr,in_arr1,in_arr2,in_arr3)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))*in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE maths_kernel_eqY

!A = A+B+C+D
SUBROUTINE maths_kernel_eqZ(out_arr, in_arr1, in_arr2, in_arr3)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE maths_kernel_eqZ

!A = A+B*C
SUBROUTINE maths_kernel_eqAA(out_arr, in_arr1, in_arr2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE maths_kernel_eqAA

SUBROUTINE maths_kernel_eqAA_fused(out_arr1, out_arr2, out_arr3, in_arr1, in_arr2, in_arr3, in_arr4)

    real(kind=8), dimension(1) :: out_arr1, out_arr2, out_arr3
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4

    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC4(0,0,0)) * in_arr2(OPS_ACC5(0,0,0)))
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + (in_arr1(OPS_ACC4(0,0,0)) * in_arr3(OPS_ACC6(0,0,0)))
    out_arr3(OPS_ACC3(0,0,0)) = out_arr3(OPS_ACC3(0,0,0)) + (in_arr1(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC7(0,0,0)))

END SUBROUTINE maths_kernel_eqAA_fused

!A = B+C*D
SUBROUTINE maths_kernel_eqAB(out_arr, in_arr1, in_arr2, in_arr3)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + (in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE maths_kernel_eqAB

!A = A+B*C*C
SUBROUTINE maths_kernel_eqAC(out_arr, in_arr1, in_arr2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE maths_kernel_eqAC

!A = A+B*C*D
SUBROUTINE maths_kernel_eqAD(out_arr, in_arr1, in_arr2, in_arr3)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE maths_kernel_eqAD

!A = A*B+C*D
SUBROUTINE maths_kernel_eqAE(out_arr, in_arr1, in_arr2, in_arr3)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = (out_arr(OPS_ACC1(0,0,0))*in_arr1(OPS_ACC2(0,0,0))) + (in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE maths_kernel_eqAE

!A = A+(B+C+D)*E
SUBROUTINE maths_kernel_eqAF(out_arr, in_arr1, in_arr2, in_arr3, in_arr4)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) &
                        + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))) * in_arr4(OPS_ACC5(0,0,0))

END SUBROUTINE maths_kernel_eqAF

!A = A+(B+C+D)*E*F
SUBROUTINE maths_kernel_eqAG(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                            + (in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0))) &
                            * in_arr4(OPS_ACC5(0,0,0))*in_arr5(OPS_ACC6(0,0,0))

END SUBROUTINE maths_kernel_eqAG

!A = A+(B*C+D*E+F*G)
SUBROUTINE maths_kernel_eqAH(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                            +(in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                            + in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                            + in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)))

END SUBROUTINE maths_kernel_eqAH

!A = A+(B*C+D*E+F*G)*H
SUBROUTINE maths_kernel_eqAI(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) &
                  + (  in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) &
                    +  in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                    +  in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0))) &
                  * in_arr7(OPS_ACC8(0,0,0))

END SUBROUTINE maths_kernel_eqAI

!A = A - B*C
SUBROUTINE maths_kernel_eqAJ(out_arr, in_arr1, in_arr2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE maths_kernel_eqAJ

!A = A-B*(C+D+E)
SUBROUTINE maths_kernel_eqAK(out_arr, in_arr1, in_arr2, in_arr3, in_arr4)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0)) &
                *(in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0))+in_arr4(OPS_ACC5(0,0,0)))

END SUBROUTINE maths_kernel_eqAK

!A = A-B*C-D*E-F*G
SUBROUTINE maths_kernel_eqAL(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                             - in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                             - in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0))

END SUBROUTINE maths_kernel_eqAL

!A = A-B*C-D*E-F*G-H*I
SUBROUTINE maths_kernel_eqAM(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                             - in_arr3(OPS_ACC4(0,0,0)) * in_arr4(OPS_ACC5(0,0,0)) &
                             - in_arr5(OPS_ACC6(0,0,0)) * in_arr6(OPS_ACC7(0,0,0)) &
                             - in_arr7(OPS_ACC8(0,0,0)) * in_arr8(OPS_ACC9(0,0,0))

END SUBROUTINE maths_kernel_eqAM

!A = B - A*C - D*E - F*G - H*I
SUBROUTINE maths_kernel_eqAN(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,in_arr7,in_arr8)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) - out_arr(OPS_ACC1(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                        - in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                        - in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) &
                                                        - in_arr7(OPS_ACC8(0,0,0))*in_arr8(OPS_ACC9(0,0,0))

END SUBROUTINE maths_kernel_eqAN

!A = -half*(B+C)
SUBROUTINE maths_kernel_eqAO(out_arr, in_arr1, in_arr2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = -0.5_8*(in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE maths_kernel_eqAO

!A = A - half*B*C
SUBROUTINE maths_kernel_eqAP(out_arr,in_arr1,in_arr2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - 0.5_8*in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE maths_kernel_eqAP

!A = A - half*(B+C+D)
SUBROUTINE maths_kernel_eqAQ(out_arr,in_arr1,in_arr2,in_arr3)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - 0.5_8*(in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE maths_kernel_eqAQ

!A = A-half*B*C*D
SUBROUTINE maths_kernel_eqAR(out_arr, in_arr1, in_arr2, in_arr3)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - 0.5_8*in_arr1(OPS_ACC2(0,0,0))&
                                                            *in_arr2(OPS_ACC3(0,0,0))*in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE maths_kernel_eqAR

!A = -half*(B+C*D+E*F+G*H)
SUBROUTINE maths_kernel_eqAS(out_arr, in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7

    out_arr(OPS_ACC1(0,0,0)) = -0.5_8 * (in_arr1(OPS_ACC2(0,0,0)) &
                        + in_arr2(OPS_ACC3(0,0,0))*in_arr3(OPS_ACC4(0,0,0)) &
                        + in_arr4(OPS_ACC5(0,0,0))*in_arr5(OPS_ACC6(0,0,0)) &
                        + in_arr6(OPS_ACC7(0,0,0))*in_arr7(OPS_ACC8(0,0,0)))

END SUBROUTINE maths_kernel_eqAS

!A = A - half*(B*C+D*E+F*G)
SUBROUTINE maths_kernel_eqAT(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - 0.5_8*( in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                                + in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                                + in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) )

END SUBROUTINE maths_kernel_eqAT

!A = A - B + half*(C*C+D*D+E*E)
SUBROUTINE maths_kernel_eqAU(out_arr,in_arr1,in_arr2,in_arr3,in_arr4)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4
    real(kind=8) :: fornow

    fornow = in_arr2(OPS_ACC3(0,0,0))*in_arr2(OPS_ACC3(0,0,0))  &
           + in_arr3(OPS_ACC4(0,0,0))*in_arr3(OPS_ACC4(0,0,0))  &
           + in_arr4(OPS_ACC5(0,0,0))*in_arr4(OPS_ACC5(0,0,0))

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0)) + 0.5_8*fornow

END SUBROUTINE maths_kernel_eqAU

SUBROUTINE maths_kernel_eqAV(out_arr1,out_arr2,in_arr1)

    real(kind=8), dimension(1) :: out_arr1, out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8) :: fornow

    fornow = alamda*EXP(rlamda*LOG(in_arr1(OPS_ACC3(0,0,0))))

    out_arr1(OPS_ACC1(0,0,0)) = fornow*out_arr2(OPS_ACC2(0,0,0))
    out_arr2(OPS_ACC2(0,0,0)) = fornow

END SUBROUTINE maths_kernel_eqAV

SUBROUTINE maths_kernel_eqAW(out_arr1,out_arr2,in_arr1,racnst,rncnst,reovrr)

    real(kind=8), dimension(1) :: out_arr1, out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: racnst,rncnst,reovrr

    out_arr2(OPS_ACC2(0,0,0)) = racnst + rncnst*LOG(in_arr1(OPS_ACC3(0,0,0)))  &
            - reovrr/in_arr1(OPS_ACC3(0,0,0))
    out_arr1(OPS_ACC1(0,0,0)) = EXP(out_arr2(OPS_ACC2(0,0,0)))

END SUBROUTINE maths_kernel_eqAW

SUBROUTINE maths_kernel_eqAX(out_arr1,out_arr2,out_arr3,in_arr1,in_arr2,racnst,rncnst,reovrr,flcnst)

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

END SUBROUTINE maths_kernel_eqAX

SUBROUTINE maths_kernel_eqAY(out_arr1,out_arr2,out_arr3,in_arr1,in_arr2,racnst,rncnst,reovrr,talpha,ovtst1,tstar2,ovtst3,cfcst1,cfcst2,encst1,encst2,dtcnst,omalph,clnten)

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
    trats2 = -tstar2/in_arr1(OPS_ACC4(0,0,0))
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

END SUBROUTINE maths_kernel_eqAY

SUBROUTINE maths_kernel_eqAZ(out_arr1,out_arr2,out_arr3,in_arr1,in_arr2,racnst,rncnst,reovrr,acfsri,bcfsrm,ovcsrm,dcfsri,ecfsri)

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

END SUBROUTINE maths_kernel_eqAZ

SUBROUTINE maths_kernel_eqBA(out_arr1,out_arr2,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6)

    real(kind=8), dimension(1) :: out_arr1,out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6
    real(kind=8) :: fornow

    fornow = in_arr1(OPS_ACC3(0,0,0))*in_arr2(OPS_ACC4(0,0,0))  &
            + in_arr3(OPS_ACC5(0,0,0))*in_arr4(OPS_ACC6(0,0,0))  &
            + in_arr5(OPS_ACC7(0,0,0))*in_arr6(OPS_ACC8(0,0,0))

!   Y EQUATION
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + fornow

!   DIFFUSION CORRECTION VELOCITY DIVERGENCE
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow

END SUBROUTINE maths_kernel_eqBA

SUBROUTINE maths_kernel_eqBB(out_arr1,out_arr2,in_arr1,in_arr2,in_arr3,in_arr4)

    real(kind=8), dimension(1) :: out_arr1,out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4
    real(kind=8) :: fornow

    fornow = ( in_arr1(OPS_ACC3(0,0,0)) + in_arr2(OPS_ACC4(0,0,0)) + in_arr3(OPS_ACC5(0,0,0)) ) &
            * in_arr4(OPS_ACC6(0,0,0))

!   Y EQUATION
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + fornow

!   DIFFUSION CORRECTION VELOCITY DIVERGENCE
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow

END SUBROUTINE maths_kernel_eqBB

SUBROUTINE maths_kernel_eqBC(transp,trun,tdifgb)

    real(kind=8), dimension(1), intent(in) :: trun
    real(kind=8), dimension(1) :: transp

    real(kind=8), intent(in) :: tdifgb

    transp(OPS_ACC1(0,0,0)) = LOG(trun(OPS_ACC2(0,0,0))/tdifgb)

END SUBROUTINE maths_kernel_eqBC

SUBROUTINE maths_kernel_eqBD(combo1,combo2,combo3,transp,yrhs,condco,ovwmol,ncocon,ncocm1,ispec)

    real(kind=8), dimension(1), intent(in) :: transp,yrhs
    real(kind=8), dimension(1) :: combo1,combo2,combo3
    integer(kind=4), intent(in) :: ncocon,ncocm1,ispec

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), dimension(nccfmx,nspcmx), intent(in) :: condco

    real(kind=8) :: fornow,ctrans
    integer(kind=4) :: icp

!   CONDUCTIVITY FOR EACH SPECIES
    fornow = condco(ncocon,ispec)
    DO icp = ncocm1,1,-1
        fornow = fornow*transp(OPS_ACC4(0,0,0)) + condco(icp,ispec)
    END DO
    ctrans = EXP(fornow)

!   COMBINATION RULE FOR CONDUCTIVITY
    fornow = yrhs(OPS_ACC5(0,0,0))*ovwmol(ispec)
    combo1(OPS_ACC1(0,0,0)) = combo1(OPS_ACC1(0,0,0)) + fornow*ctrans
    combo2(OPS_ACC2(0,0,0)) = combo2(OPS_ACC2(0,0,0)) + fornow/ctrans
    combo3(OPS_ACC3(0,0,0)) = combo3(OPS_ACC3(0,0,0)) + fornow

END SUBROUTINE maths_kernel_eqBD

SUBROUTINE maths_kernel_eqBE(combo1,combo2,combo3,store7,wmomix,drhs)

    real(kind=8), dimension(1), intent(in) :: drhs
    real(kind=8), dimension(1) :: combo1,combo2,combo3,store7,wmomix

    combo3(OPS_ACC3(0,0,0)) = 1.0_8/combo3(OPS_ACC3(0,0,0))
    combo1(OPS_ACC1(0,0,0)) = combo1(OPS_ACC1(0,0,0))*combo3(OPS_ACC3(0,0,0))
    combo2(OPS_ACC2(0,0,0)) = combo2(OPS_ACC2(0,0,0))*combo3(OPS_ACC3(0,0,0))

    store7(OPS_ACC4(0,0,0)) = 0.5_8*(combo1(OPS_ACC1(0,0,0)) + 1.0_8/combo2(OPS_ACC2(0,0,0)))
    wmomix(OPS_ACC5(0,0,0)) = drhs(OPS_ACC6(0,0,0))*combo3(OPS_ACC3(0,0,0))

END SUBROUTINE maths_kernel_eqBE

SUBROUTINE maths_kernel_eqBCDE(transp,store7,wmomix,trun,drhs,yrhs_md,condco,ovwmol,tdifgb,ncocon,ncocm1)

    real(kind=8), dimension(1) :: transp,store7,wmomix
    real(kind=8), dimension(1), intent(in) :: trun,drhs
    real(kind=8), dimension(1), intent(in) :: yrhs_md

    real(kind=8), dimension(nccfmx,nspcmx), intent(in) :: condco
    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), intent(in) :: tdifgb
    integer(kind=4), intent(in) :: ncocon,ncocm1

!   Local data
    real(kind=8) :: fornow, combo1, combo2, combo3
    real(kind=8) :: ctrans(nspcmx)
    integer(kind=4) :: ispec,icp

!   CONDUCTIVITY FOR EACH SPECIES
    transp(OPS_ACC1(0,0,0)) = LOG(trun(OPS_ACC4(0,0,0))/tdifgb)
    DO ispec = 1, nspcmx
        fornow = condco(ncocon,ispec)
        DO icp = ncocm1,1,-1
            fornow = fornow*transp(OPS_ACC1(0,0,0)) + condco(icp,ispec)
        END DO
        ctrans(ispec) = EXP(fornow)
    END DO

!   COMBINATION RULE FOR CONDUCTIVITY
    combo1 = 0.0_8
    combo2 = 0.0_8
    combo3 = 0.0_8
    DO ispec = 1, nspcmx
        fornow = yrhs_md(OPS_ACC_MD6(ispec,0,0,0))*ovwmol(ispec)
        combo1 = combo1 + fornow*ctrans(ispec)
        combo2 = combo2 + fornow/ctrans(ispec)
        combo3 = combo3 + fornow
    END DO

    combo3 = 1.0_8/combo3
    combo1 = combo1*combo3
    combo2 = combo2*combo3
    store7(OPS_ACC2(0,0,0)) = 0.5_8*(combo1 + 1.0_8/combo2)
    wmomix(OPS_ACC3(0,0,0)) = drhs(OPS_ACC5(0,0,0))*combo3

END SUBROUTINE maths_kernel_eqBCDE

SUBROUTINE maths_kernel_eqBF(ctrans,combo1,combo2,transp,prun,yrhs,diffco,ovwmol,pdifgb,dfctol,ncodif,ncodm1,jspec,ispec)

    real(kind=8), dimension(1), intent(in) :: transp,prun,yrhs
    real(kind=8), dimension(1) :: ctrans,combo1,combo2
    integer(kind=4), intent(in) :: ncodif,ncodm1,jspec,ispec

    real(kind=8), dimension(ndcfmx,nspcmx,nspcmx), intent(in) :: diffco
    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), intent(in) :: pdifgb,dfctol

    real(kind=8) :: fornow
    integer(kind=4) :: icp

    fornow = diffco(ncodif,jspec,ispec)
    DO icp = ncodm1,1,-1
        fornow = fornow*transp(OPS_ACC4(0,0,0)) + diffco(icp,jspec,ispec)
    END DO

    ctrans(OPS_ACC1(0,0,0)) = EXP(fornow)*pdifgb/prun(OPS_ACC5(0,0,0))

    fornow = yrhs(OPS_ACC6(0,0,0)) + dfctol
    combo1(OPS_ACC2(0,0,0)) = combo1(OPS_ACC2(0,0,0)) + fornow
    combo2(OPS_ACC3(0,0,0)) = combo2(OPS_ACC3(0,0,0)) + fornow*ovwmol(jspec)/ctrans(OPS_ACC1(0,0,0))

END SUBROUTINE maths_kernel_eqBF

SUBROUTINE maths_kernel_eqBG(combo1,combo2,difmix,store7,ctrans,yrhs,wmomix,drhs,ovwmol,dfctol,ispec)

    real(kind=8), dimension(1), intent(in) :: ctrans,yrhs,wmomix,drhs
    real(kind=8), dimension(1) :: combo1,combo2,difmix,store7

    integer(kind=4), intent(in) :: ispec

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), intent(in) :: dfctol

    real(kind=8) :: fornow

    fornow = yrhs(OPS_ACC6(0,0,0)) + dfctol
    combo1(OPS_ACC1(0,0,0)) = combo1(OPS_ACC1(0,0,0)) - fornow
    combo2(OPS_ACC2(0,0,0)) = combo2(OPS_ACC2(0,0,0)) - fornow*ovwmol(ispec)/ctrans(OPS_ACC5(0,0,0))
    combo2(OPS_ACC2(0,0,0)) = combo2(OPS_ACC2(0,0,0))*wmomix(OPS_ACC7(0,0,0))

    difmix(OPS_ACC3(0,0,0)) = drhs(OPS_ACC8(0,0,0))*combo1(OPS_ACC1(0,0,0))/combo2(OPS_ACC2(0,0,0))
    store7(OPS_ACC4(0,0,0)) = difmix(OPS_ACC3(0,0,0))

END SUBROUTINE maths_kernel_eqBG

SUBROUTINE maths_kernel_eqBFG(difmix,store7,transp,prun,yrhs_md,wmomix,drhs,diffco,ovwmol,pdifgb,dfctol,ncodif,ncodm1,ispec)

!   Function Arguments
    real(kind=8), dimension(1) :: difmix, store7
    real(kind=8), dimension(1), intent(in) :: transp,prun
    real(kind=8), dimension(1), intent(in) :: yrhs_md
    real(kind=8), dimension(1), intent(in) :: wmomix,drhs

    real(kind=8), dimension(ndcfmx,nspcmx,nspcmx), intent(in) :: diffco
    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), intent(in) :: pdifgb,dfctol

    integer(kind=4), intent(in) :: ncodif,ncodm1,ispec

!   Local data
    real(kind=8) :: fornow, combo1, combo2
    real(kind=8) :: ctrans(nspcmx)
    integer(kind=4) :: jspec,icp

    DO jspec = 1, nspcmx
        fornow = diffco(ncodif,jspec,ispec)
        DO icp = ncodm1,1,-1
            fornow = fornow*transp(OPS_ACC3(0,0,0)) + diffco(icp,jspec,ispec)
        END DO
        ctrans(jspec) = EXP(fornow)*pdifgb/prun(OPS_ACC4(0,0,0))
    END DO

!   COMBINATION RULE FOR MASS DIFFUSIVITY
    combo1 = 0.0_8
    combo2 = 0.0_8
    DO jspec = 1, nspcmx
        fornow = yrhs_md(OPS_ACC_MD5(jspec,0,0,0)) + dfctol
        combo1 = combo1 + fornow
        combo2 = combo2 + fornow*ovwmol(jspec)/ctrans(jspec)
    END DO
    fornow = yrhs_md(OPS_ACC_MD5(ispec,0,0,0)) + dfctol
    combo1 = combo1 - fornow
    combo2 = combo2 - fornow*ovwmol(ispec)/ctrans(ispec)
    combo2 = combo2*wmomix(OPS_ACC6(0,0,0))
    difmix(OPS_ACC1(0,0,0)) = drhs(OPS_ACC7(0,0,0))*combo1/combo2
    store7(OPS_ACC2(0,0,0)) = difmix(OPS_ACC1(0,0,0))

END SUBROUTINE maths_kernel_eqBFG

SUBROUTINE maths_kernel_eqBH(tdrmix,trun,yrhs,wmomix,tdrcco,ovwmol,tdifgb,ncotdr,jspec,ispec,ncotm1)

    real(kind=8), dimension(1) :: tdrmix
    real(kind=8), dimension(1), intent(in) :: trun,yrhs,wmomix

    integer(kind=4), intent(in) :: ncotdr,jspec,ispec,ncotm1

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), dimension(ndcfmx,nspcmx,nspcmx), intent(in) :: tdrcco
    real(kind=8), intent(in) :: tdifgb

    real(kind=8) :: fornow, combo2
    real(kind=8) :: ctrans
    integer(kind=4) :: icp

!   THERMAL DIFFUSION RATIO FOR THIS SPECIES PAIR
    combo2 = trun(OPS_ACC2(0,0,0))/tdifgb
    fornow = tdrcco(ncotdr,jspec,ispec)
    DO icp = ncotm1,1,-1
        fornow = fornow*combo2 + tdrcco(icp,jspec,ispec)
    END DO
    ctrans = fornow

!   COMBINATION RULE FOR THERMAL DIFFUSIION RATIO
    fornow = yrhs(OPS_ACC3(0,0,0))*ovwmol(jspec)
    tdrmix(OPS_ACC1(0,0,0)) = tdrmix(OPS_ACC1(0,0,0))  &
                     + fornow*wmomix(OPS_ACC4(0,0,0))*ctrans

END SUBROUTINE maths_kernel_eqBH

SUBROUTINE maths_kernel_eqBI(ctrans,transp,viscco,ncovis,ncovm1,ispec)

    real(kind=8), dimension(1), intent(in) :: transp
    real(kind=8), dimension(1) :: ctrans
    integer(kind=4), intent(in) :: ncovis,ncovm1,ispec

    real(kind=8), dimension(nvcfmx,nspcmx), intent(in) :: viscco

    real(kind=8) :: fornow
    integer(kind=4) :: icp

    fornow = viscco(ncovis,ispec)
    DO icp = ncovm1,1,-1
        fornow = fornow*transp(OPS_ACC2(0,0,0)) + viscco(icp,ispec)
    END DO

    ctrans(OPS_ACC1(0,0,0)) = EXP(fornow)

END SUBROUTINE maths_kernel_eqBI

SUBROUTINE maths_kernel_eqBJ(combo2,ctrans_ispec,ctrans_jspec,yrhs,ovwmol,wilko1,wilko2,ispec,jspec)

    real(kind=8), dimension(1), intent(in) :: ctrans_ispec,ctrans_jspec,yrhs
    real(kind=8), dimension(1) :: combo2
    integer(kind=4), intent(in) :: ispec,jspec

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), dimension(nspcmx,nspcmx), intent(in) :: wilko1,wilko2

    real(kind=8) :: fornow

    fornow = SQRT(ctrans_ispec(OPS_ACC2(0,0,0))/ctrans_jspec(OPS_ACC3(0,0,0)))
    fornow = 1.0_8 + fornow*wilko2(jspec,ispec)
    fornow = wilko1(jspec,ispec)*fornow*fornow
    combo2(OPS_ACC1(0,0,0)) = combo2(OPS_ACC1(0,0,0)) + yrhs(OPS_ACC4(0,0,0))*ovwmol(jspec)*fornow

END SUBROUTINE maths_kernel_eqBJ

SUBROUTINE maths_kernel_eqBK(combo1,ctrans,combo2,yrhs,ovwmol,ispec)

    real(kind=8), dimension(1), intent(in) :: ctrans,combo2,yrhs
    real(kind=8), dimension(1) :: combo1
    integer(kind=4), intent(in) :: ispec

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol

    real(kind=8) :: fornow

    fornow = ctrans(OPS_ACC2(0,0,0))/combo2(OPS_ACC3(0,0,0))
    combo1(OPS_ACC1(0,0,0)) = combo1(OPS_ACC1(0,0,0)) + yrhs(OPS_ACC4(0,0,0))*ovwmol(ispec)*fornow

END SUBROUTINE maths_kernel_eqBK

SUBROUTINE maths_kernel_eqBIJK(difmix,transp,yrhs_md,viscco,wilko1,wilko2,ovwmol,ncovis,ncovm1)

!   Function Arguments
    real(kind=8), dimension(1) :: difmix
    real(kind=8), dimension(1), intent(in) :: transp
    real(kind=8), dimension(1), intent(in) :: yrhs_md
    real(kind=8), dimension(nvcfmx,nspcmx), intent(in) :: viscco
    real(kind=8), dimension(nspcmx,nspcmx), intent(in) :: wilko1,wilko2
    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    integer(kind=4), intent(in) :: ncovis,ncovm1

!   Local data
    real(kind=8) :: fornow, combo1, combo2
    real(kind=8) :: ctrans(nspcmx)
    integer(kind=4) :: ispec,jspec,icp

!   VISCOSITY FOR EACH SPECIES
    DO ispec = 1, nspcmx
        fornow = viscco(ncovis,ispec)
        DO icp = ncovm1,1,-1
            fornow = fornow*transp(OPS_ACC2(0,0,0)) + viscco(icp,ispec)
        END DO
        ctrans(ispec) = EXP(fornow)
    END DO

!   COMBINATION RULE FOR VISCOSITY
    combo1 = 0.0_8
    DO ispec = 1, nspcmx
        combo2 = 0.0_8
        DO jspec = 1, nspcmx
            fornow = SQRT(ctrans(ispec)/ctrans(jspec))
            fornow = 1.0_8 + fornow*wilko2(jspec,ispec)
            fornow = wilko1(jspec,ispec)*fornow*fornow
            combo2 = combo2 + yrhs_md(OPS_ACC_MD3(jspec,0,0,0))*ovwmol(jspec)*fornow
        END DO
        fornow = ctrans(ispec)/combo2
        combo1 = combo1 + yrhs_md(OPS_ACC_MD3(ispec,0,0,0))*ovwmol(ispec)*fornow
    END DO

    difmix(OPS_ACC1(0,0,0)) = combo1

END SUBROUTINE maths_kernel_eqBIJK

SUBROUTINE maths_kernel_eqBL(utmp,wtmp,trun,itndex,yrhs,amasch,ncpoly,ncpom1,ncenth,icoef1,icoef2,ispec)

    real(kind=8), dimension(1) :: utmp,wtmp
    real(kind=8), dimension(1), intent(in) :: trun,yrhs
    integer(kind=4), dimension(1), intent(in) :: itndex

    integer(kind=4), intent(in) ::  icoef1,icoef2,ispec

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly, ncpom1, ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 + MOD(itndex(OPS_ACC4(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*trun(OPS_ACC3(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    utmp(OPS_ACC1(0,0,0)) = amasch(ncenth(itint,ispec),itint,ispec)  &
                    + fornow*trun(OPS_ACC3(0,0,0))

!   MIXTURE H
    wtmp(OPS_ACC2(0,0,0)) = wtmp(OPS_ACC2(0,0,0)) + utmp(OPS_ACC1(0,0,0))*yrhs(OPS_ACC5(0,0,0))

END SUBROUTINE maths_kernel_eqBL

SUBROUTINE maths_kernel_eqBM(store1,yrhs,ovwmol,scoef,ysmall,ydenom,ispec)

    real(kind=8), dimension(1) :: store1
    real(kind=8), dimension(1), intent(in) :: yrhs
    integer(kind=4), intent(in) ::  ispec
    real(kind=8), intent(in) :: scoef,ysmall,ydenom
    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8) :: fornow

    fornow = yrhs(OPS_ACC2(0,0,0))*ovwmol(ispec)
    fornow = MAX(fornow,ysmall)
    fornow = EXP(scoef*LOG(fornow))
    fornow = fornow/(1.0_8+ydenom*fornow)
    store1(OPS_ACC1(0,0,0)) = store1(OPS_ACC1(0,0,0))*fornow

END SUBROUTINE maths_kernel_eqBM

SUBROUTINE maths_kernel_eqBN(store2,trun,itndex,amolgb,ncpoly,ncpom1,ncenth,ncenpy,diffmu,isspec,ispec,istep,ipower,icoef1,icoef2)

    real(kind=8), dimension(1) :: store2
    real(kind=8), dimension(1), intent(in) :: trun
    integer(kind=4), dimension(1), intent(in) :: itndex
    integer(kind=4), intent(in) ::  isspec,ispec,istep,ipower,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amolgb
    real(kind=8), dimension(nssmax,nstpmx), intent(in) :: diffmu
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth,ncenpy

    integer(kind=4) :: itint,icp
    real(kind=8) :: gibbsp

!   LOCATE THE TEMPERATURE IN AN INTERVAL
    itint = 1 + MOD(itndex(OPS_ACC3(0,0,0)),icoef1)/icoef2

!   CONSTRUCT GIBBS FUNCTION FROM ITS POLYNOMIAL COEFFICIENTS
    gibbsp = amolgb(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        gibbsp = amolgb(icp,itint,ispec) + gibbsp*trun(OPS_ACC2(0,0,0))
    END DO

    gibbsp = amolgb(ncenth(itint,ispec),itint,ispec)  &
            /trun(OPS_ACC2(0,0,0)) - amolgb(ncenpy(itint,ispec),itint,ispec)  &
            *LOG(trun(OPS_ACC2(0,0,0))) - gibbsp

!   ADD GIBBS FUNCTION CONTRIBUTION TO RATE COEFFICIENT
!   USING STEP SPECIES DELTA-LIST
!   TO GET BACKWARD RATE COEFFICIENT
    store2(OPS_ACC1(0,0,0)) = store2(OPS_ACC1(0,0,0)) + diffmu(isspec,istep)*gibbsp

END SUBROUTINE maths_kernel_eqBN

SUBROUTINE maths_kernel_eqBO(itndex,trun,tinthi,ntint,ipower,ispec)

    real(kind=8), dimension(1), intent(in) :: trun
    integer(kind=4), dimension(1) :: itndex

    integer(kind=4), intent(in) ::  ipower,ispec

    real(kind=8), dimension(ntinmx,nspcmx), intent(in) :: tinthi
    integer(kind=4), dimension(nspcmx), intent(in) :: ntint

    integer(kind=4) :: itint

    itint = 1
    DO WHILE (trun(OPS_ACC2(0,0,0)) > tinthi(itint,ispec) .and. itint < ntint(ispec))
        itint = itint + 1
    END DO

    itndex(OPS_ACC1(0,0,0)) = itndex(OPS_ACC1(0,0,0)) + (itint-1)*ntbase**ipower

END SUBROUTINE maths_kernel_eqBO

SUBROUTINE maths_kernel_eqBP(erun,trun,itndex,yrun,amasch,ncpoly,ncpom1,ncenth,ipower,icoef1,icoef2,ispec)

    real(kind=8), dimension(1), intent(in) :: trun, yrun
    real(kind=8), dimension(1) :: erun
    integer(kind=4), dimension(1), intent(in) :: itndex

    integer(kind=4), intent(in) ::  ipower,icoef1,icoef2,ispec

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: icp,itint

    itint = 1 + MOD(itndex(OPS_ACC3(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*trun(OPS_ACC2(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
           + fornow*trun(OPS_ACC2(0,0,0))

    erun(OPS_ACC1(0,0,0)) = erun(OPS_ACC1(0,0,0)) + fornow*yrun(OPS_ACC4(0,0,0))

END SUBROUTINE maths_kernel_eqBP

SUBROUTINE maths_kernel_eqBQ(utgv, vtgv, wtgv, drun, deltagx, deltagy, deltagz, tkes)

    real(kind=8), dimension(1), intent(in) :: utgv, vtgv, wtgv, drun
    real(kind=8), intent(in) :: deltagx, deltagy, deltagz
    real(kind=8) :: tkes

    real(kind=8) :: tkel

    tkel = (utgv(OPS_ACC1(0,0,0))*utgv(OPS_ACC1(0,0,0))) + &
           (vtgv(OPS_ACC2(0,0,0))*vtgv(OPS_ACC2(0,0,0))) + &
           (wtgv(OPS_ACC3(0,0,0))*wtgv(OPS_ACC3(0,0,0)))

    tkes = tkes + (0.5_8*drun(OPS_ACC4(0,0,0))*tkel &
                  *deltagx*deltagy*deltagz)

END SUBROUTINE maths_kernel_eqBQ

SUBROUTINE maths_kernel_eqBR_xdir(t1bx, store2, store3)

    real(kind=8), dimension(1), intent(in) :: store2, store3
    real(kind=8), dimension(1) :: t1bx

    t1bx(OPS_ACC1(0,0,0)) = -(store2(OPS_ACC2(0,0,0)) + store3(OPS_ACC3(0,0,0)))

END SUBROUTINE maths_kernel_eqBR_xdir

SUBROUTINE maths_kernel_eqBR_ydir(t1by, store1, store3)

    real(kind=8), dimension(1), intent(in) :: store1, store3
    real(kind=8), dimension(1) :: t1by

    t1by(OPS_ACC1(0,0,0)) = -(store1(OPS_ACC2(0,0,0)) + store3(OPS_ACC3(0,0,0)))

END SUBROUTINE maths_kernel_eqBR_ydir

SUBROUTINE maths_kernel_eqBR_zdir(t1bz, store1, store2)

    real(kind=8), dimension(1), intent(in) :: store1, store2
    real(kind=8), dimension(1) :: t1bz

    t1bz(OPS_ACC1(0,0,0)) = -(store1(OPS_ACC2(0,0,0)) + store2(OPS_ACC3(0,0,0)))

END SUBROUTINE maths_kernel_eqBR_zdir

SUBROUTINE maths_kernel_eqBS(tdrmix, yrhs, wmomix, ovwmol, ispec)

    real(kind=8), dimension(1), intent(in) :: yrhs, wmomix
    real(kind=8), dimension(1) :: tdrmix

    integer(kind=4), intent(in) :: ispec

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol

    tdrmix(OPS_ACC1(0,0,0)) = tdrmix(OPS_ACC1(0,0,0))*yrhs(OPS_ACC2(0,0,0))  &
                            * wmomix(OPS_ACC3(0,0,0))*ovwmol(ispec)

END SUBROUTINE maths_kernel_eqBS

SUBROUTINE maths_kernel_eqCQ(dvtgvdx, dwtgvdx, dutgvdy, dwtgvdy, dutgvdz, dvtgvdz, drun, deltagx, deltagy, deltagz, enstro)

    real(kind=8), dimension(1), intent(in) :: dvtgvdx, dwtgvdx, dutgvdy, dwtgvdy, dutgvdz, dvtgvdz, drun
    real(kind=8), intent(in) :: deltagx, deltagy, deltagz
    real(kind=8) :: enstro

    real(kind=8) :: vort1, vort2, vort3, enstrs

    vort1 = dwtgvdy(OPS_ACC4(0,0,0))-dvtgvdz(OPS_ACC6(0,0,0))
    vort2 = dutgvdz(OPS_ACC5(0,0,0))-dwtgvdx(OPS_ACC2(0,0,0))
    vort3 = dvtgvdx(OPS_ACC1(0,0,0))-dutgvdy(OPS_ACC3(0,0,0))
    enstrs = vort1*vort1+vort2*vort2+vort3*vort3
    enstro = enstro+0.5_8*drun(OPS_ACC7(0,0,0))*enstrs*deltagx*deltagy*deltagz

END SUBROUTINE maths_kernel_eqCQ

SUBROUTINE maths_kernel_eqtau(out_arr1,out_arr2,in_arr1,in_arr2,in_arr3,var)

    real(kind=8), dimension(1) :: out_arr1, out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    real(kind=8), intent(in) :: var
    real(kind=8) :: fornow

    fornow = in_arr1(OPS_ACC3(0,0,0))*var*in_arr2(OPS_ACC4(0,0,0))
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + fornow
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow*in_arr3(OPS_ACC5(0,0,0))

END SUBROUTINE maths_kernel_eqtau

SUBROUTINE maths_kernel_print_drhs(out_arr, itime)

    real(kind=8), dimension(1), intent(in) :: out_arr
    integer(kind=4), intent(in) :: itime

!    WRITE(*,'(a,I7,a,F16.8)')  &
!      'test_drhs: (step=',itime,') value: ',out_arr(OPS_ACC1(0,0,0))
#if !defined(OPS_WITH_CUDAFOR) && !defined(OPS_WITH_OMPOFFLOADFOR)
    print '(a,I7,a,F16.8)', 'test_drhs: (step=', itime, ') value: ', out_arr(OPS_ACC1(0,0,0))
#endif

END SUBROUTINE maths_kernel_print_drhs

SUBROUTINE maths_kernel_print_erhs(out_arr, itime)

    real(kind=8), dimension(1), intent(in) :: out_arr
    integer(kind=4), intent(in) :: itime

!    WRITE(*,'(a,I7,a,F16.8)')  &
!      'test_erhs: (step=',itime,') value: ',out_arr(OPS_ACC1(0,0,0))
#if !defined(OPS_WITH_CUDAFOR) && !defined(OPS_WITH_OMPOFFLOADFOR)
    print '(a,I7,a,F16.8)', 'test_erhs: (step=', itime, ') value: ', out_arr(OPS_ACC1(0,0,0))
#endif

END SUBROUTINE maths_kernel_print_erhs

SUBROUTINE maths_kernel_print_urhs(out_arr, itime)

    real(kind=8), dimension(1), intent(in) :: out_arr
    integer(kind=4), intent(in) :: itime

!    WRITE(*,'(a,I7,a,F16.8)')  &
!      'test_urhs: (step=',itime,') value: ',out_arr(OPS_ACC1(0,0,0))
#if !defined(OPS_WITH_CUDAFOR) && !defined(OPS_WITH_OMPOFFLOADFOR)
    print '(a,I7,a,F16.8)', 'test_urhs: (step=', itime, ') value: ', out_arr(OPS_ACC1(0,0,0))
#endif

END SUBROUTINE maths_kernel_print_urhs
