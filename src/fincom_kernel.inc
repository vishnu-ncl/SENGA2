SUBROUTINE fincom_kernel_main(err_arr,run_arr,rhs_arr,rkerr_nrkstp,rklhs_nrkstp)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: err_arr, run_arr, rhs_arr
    real(kind=dp), intent(in) :: rkerr_nrkstp, rklhs_nrkstp

    err_arr(OPS_ACC1(0,0,0)) = err_arr(OPS_ACC1(0,0,0)) + rkerr_nrkstp * rhs_arr(OPS_ACC3(0,0,0))

    run_arr(OPS_ACC2(0,0,0)) = run_arr(OPS_ACC2(0,0,0)) + rklhs_nrkstp * rhs_arr(OPS_ACC3(0,0,0))
    rhs_arr(OPS_ACC3(0,0,0)) = run_arr(OPS_ACC2(0,0,0))

END SUBROUTINE fincom_kernel_main
