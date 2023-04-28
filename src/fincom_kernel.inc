SUBROUTINE fincom_kernel_main(err_arr,run_arr,rhs_arr,rkerr,rklhs,nrkstp)
    use data_types
    implicit none

    real(8), dimension(1) :: err_arr, run_arr, rhs_arr
    integer(4), intent(in) :: nrkstp
    real(8), dimension(nrkmax), intent(in) :: rkerr, rklhs

    err_arr(OPS_ACC1(0,0,0)) = err_arr(OPS_ACC1(0,0,0)) + rkerr(nrkstp) * rhs_arr(OPS_ACC3(0,0,0))

    run_arr(OPS_ACC2(0,0,0)) = run_arr(OPS_ACC2(0,0,0)) + rklhs(nrkstp) * rhs_arr(OPS_ACC3(0,0,0))
    rhs_arr(OPS_ACC3(0,0,0)) = run_arr(OPS_ACC2(0,0,0))

END SUBROUTINE fincom_kernel_main
