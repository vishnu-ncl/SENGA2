SUBROUTINE fincom_kernel_main(err_arr,run_arr,rhs_arr,rkerr_nrkstp,rklhs_nrkstp)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: err_arr, run_arr, rhs_arr
    real(kind=dp), intent(in) :: rkerr_nrkstp, rklhs_nrkstp

    err_arr(OPS_ACC1(0,0,0)) = err_arr(OPS_ACC1(0,0,0)) + rkerr_nrkstp * rhs_arr(OPS_ACC3(0,0,0))

    run_arr(OPS_ACC2(0,0,0)) = run_arr(OPS_ACC2(0,0,0)) + rklhs_nrkstp * rhs_arr(OPS_ACC3(0,0,0))
    rhs_arr(OPS_ACC3(0,0,0)) = run_arr(OPS_ACC2(0,0,0))

END SUBROUTINE fincom_kernel_main

SUBROUTINE fincom_kernel_MD(err_arr,run_arr,rhs_arr,rkerr_irkstp,rklhs_irkstp,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: err_arr, run_arr, rhs_arr
    real(kind=dp), intent(in) :: rkerr_irkstp, rklhs_irkstp
    integer, intent(in) :: ispec

    err_arr(OPS_ACC_MD1(ispec,0,0,0)) = err_arr(OPS_ACC_MD1(ispec,0,0,0)) + rkerr_irkstp * rhs_arr(OPS_ACC_MD3(ispec,0,0,0))

    run_arr(OPS_ACC_MD2(ispec,0,0,0)) = run_arr(OPS_ACC_MD2(ispec,0,0,0)) + rklhs_irkstp * rhs_arr(OPS_ACC_MD3(ispec,0,0,0))
    rhs_arr(OPS_ACC_MD3(ispec,0,0,0)) = run_arr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE fincom_kernel_MD
