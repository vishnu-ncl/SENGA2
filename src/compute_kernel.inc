!equation1 A = LOG(B)
!A = -B
!equation2 A = B-var
!A = A*var
! A = var*B*C
!equation3 A = B*var
!equation4 A = A+B
!equation5 A = B+C
!equation6 A = B+C+D
!equation7 A = A+B+C+D
!equation8 A = A+B*C
!A = A+B*C*C
!A = A+B*C*D
!A = A*B+C*D
!equation9 A = B+C*D
!equation10 A = A/B
!equation11 A = B/C
!equation12 A = A*B
!A = A*B+C*D
!equation13 A = B*C
!equation14 A = A-B
!A = A-B*C

!A = LOG(B)
SUBROUTINE compute_kernel_AequalslogB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    
    out_arr(OPS_ACC1(0,0,0)) = LOG(in_arr1(OPS_ACC2(0,0,0)))

END SUBROUTINE compute_kernel_AequalslogB

SUBROUTINE compute_kernel_AequalsnegateB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = -in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE compute_kernel_AequalsnegateB

!A = B-var
SUBROUTINE compute_kernel_AequalsBminusvar(out_arr, in_arr1, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) - var

END SUBROUTINE compute_kernel_AequalsBminusvar

!A = A*var
SUBROUTINE compute_kernel_AequalsAmulvar(out_arr, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * var

END SUBROUTINE compute_kernel_AequalsAmulvar

!A = B*var
SUBROUTINE compute_kernel_AequalsBmulvar(out_arr, in_arr1, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * var

END SUBROUTINE compute_kernel_AequalsBmulvar

!A = B*C*var
SUBROUTINE compute_kernel_AequalsBmulCmulvar(out_arr, in_arr1, in_arr2, var)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=dp), intent(in) :: var

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * var    

END SUBROUTINE compute_kernel_AequalsBmulCmulvar

!A = A+B
SUBROUTINE compute_kernel_AequalsAplusB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE compute_kernel_AequalsAplusB

!A = B+C 
SUBROUTINE compute_kernel_AequalsBplusC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE compute_kernel_AequalsBplusC

!A = B+C+D
SUBROUTINE compute_kernel_AequalsBCDplus(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE compute_kernel_AequalsBCDplus

!A = A+B+C+D
SUBROUTINE compute_kernel_AequalsABCDplus(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE compute_kernel_AequalsABCDplus

!A = A+B*C
SUBROUTINE compute_kernel_AequalsAplusBmulC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE compute_kernel_AequalsAplusBmulC

!A = A+B*C*C
SUBROUTINE compute_kernel_AequalsAplusBmulCmulC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE compute_kernel_AequalsAplusBmulCmulC

!A = A+B*C*D
SUBROUTINE compute_kernel_AequalsAplusBmulCmulD(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE compute_kernel_AequalsAplusBmulCmulD

!A = A*B+C*D
SUBROUTINE compute_kernel_AequalsAmulBplusCmulD(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = (out_arr(OPS_ACC1(0,0,0))*in_arr1(OPS_ACC2(0,0,0))) + (in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE compute_kernel_AequalsAmulBplusCmulD

!A = B+C*D
SUBROUTINE compute_kernel_AequalsBplusCmulD(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + (in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE compute_kernel_AequalsBplusCmulD

!A = A/B
SUBROUTINE compute_kernel_AequalsAdivB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) / in_arr1(OPS_ACC2(0,0,0)) 

END SUBROUTINE compute_kernel_AequalsAdivB

!A = B/C
SUBROUTINE compute_kernel_AequalsBdivC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) / in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE compute_kernel_AequalsBdivC

!A = A*B
SUBROUTINE compute_kernel_AequalsAmulB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE compute_kernel_AequalsAmulB

!A = B*C
SUBROUTINE compute_kernel_AequalsBmulC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE compute_kernel_AequalsBmulC

!A = A-B
SUBROUTINE compute_kernel_AequalsAminusB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE compute_kernel_AequalsAminusB

!A = A-B*C
SUBROUTINE compute_kernel_AequalsAminusBmulC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE compute_kernel_AequalsAminusBmulC


