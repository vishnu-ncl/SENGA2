SUBROUTINE lincom_kernel_main(err_arr,run_arr,rhs_arr,rkerr,rklhs,rkrhs,irkstp)
    use data_types
    implicit none

    real(8), dimension(1) :: err_arr, run_arr, rhs_arr
    integer(4), intent(in) :: irkstp
    real(8), dimension(nrkmax), intent(in) :: rkerr, rklhs, rkrhs
    real(8) :: fornow

    err_arr(OPS_ACC1(0,0,0)) = err_arr(OPS_ACC1(0,0,0)) + rkerr(irkstp) * rhs_arr(OPS_ACC3(0,0,0))

    fornow = run_arr(OPS_ACC2(0,0,0))
    run_arr(OPS_ACC2(0,0,0)) = fornow + rklhs(irkstp) * rhs_arr(OPS_ACC3(0,0,0))
    rhs_arr(OPS_ACC3(0,0,0)) = fornow + rkrhs(irkstp) * rhs_arr(OPS_ACC3(0,0,0))

END SUBROUTINE lincom_kernel_main
