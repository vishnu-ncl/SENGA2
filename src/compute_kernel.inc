SUBROUTINE compute_kernel_AequalslogB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    
    out_arr(OPS_ACC1(0,0,0)) = LOG(in_arr1(OPS_ACC2(0,0,0)))

END SUBROUTINE compute_kernel_AequalslogB

SUBROUTINE compute_kernel_AequalsAplusB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE compute_kernel_AequalsAplusB

SUBROUTINE compute_kernel_AequalsBplusC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE compute_kernel_AequalsBplusC

SUBROUTINE compute_kernel_AequalsBCDplus(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE compute_kernel_AequalsBCDplus

SUBROUTINE compute_kernel_AequalsABCDplus(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr1(OPS_ACC2(0,0,0)) + in_arr2(OPS_ACC3(0,0,0)) + in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE compute_kernel_AequalsABCDplus

SUBROUTINE compute_kernel_AequalsAplusBmulC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + (in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE compute_kernel_AequalsAplusBmulC

SUBROUTINE compute_kernel_AequalsBplusCmulD(out_arr, in_arr1, in_arr2, in_arr3)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) + (in_arr2(OPS_ACC3(0,0,0)) * in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE compute_kernel_AequalsBplusCmulD

SUBROUTINE compute_kernel_AequalsAdivB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) / in_arr1(OPS_ACC2(0,0,0)) 

END SUBROUTINE compute_kernel_AequalsAdivB

SUBROUTINE compute_kernel_AequalsBdivC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) / in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE compute_kernel_AequalsBdivC

SUBROUTINE compute_kernel_AequalsAmulB(out_arr, in_arr1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE compute_kernel_AequalsAmulB

SUBROUTINE compute_kernel_AequalsBmulC(out_arr, in_arr1, in_arr2)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0)) * in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE compute_kernel_AequalsBmulC

