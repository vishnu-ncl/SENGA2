SUBROUTINE fincom_kernel_main(err_arr,run_arr,rhs_arr,rkerr,rklhs,nrkmax,nrkstp)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: err_arr, run_arr, rhs_arr
    integer, intent(in) :: nrkmax,nrkstp
    real(kind=8), dimension(nrkmax), intent(in) :: rkerr, rklhs

    err_arr(OPS_ACC1(0,0,0)) = err_arr(OPS_ACC1(0,0,0)) + rkerr(nrkstp) * rhs_arr(OPS_ACC3(0,0,0))

    run_arr(OPS_ACC2(0,0,0)) = run_arr(OPS_ACC2(0,0,0)) + rklhs(nrkstp) * rhs_arr(OPS_ACC3(0,0,0))
    rhs_arr(OPS_ACC3(0,0,0)) = run_arr(OPS_ACC2(0,0,0))

END SUBROUTINE fincom_kernel_main

SUBROUTINE fincom_kernel_MD(err_arr,run_arr,rhs_arr,rkerr,rklhs,nrkmax,nrkstp,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: err_arr, run_arr, rhs_arr
    integer, intent(in) :: ispec,nrkmax,nrkstp
    real(kind=8), dimension(nrkmax), intent(in) :: rkerr, rklhs

    err_arr(OPS_ACC_MD1(ispec,0,0,0)) = err_arr(OPS_ACC_MD1(ispec,0,0,0)) + rkerr(nrkstp) * rhs_arr(OPS_ACC_MD3(ispec,0,0,0))

    run_arr(OPS_ACC_MD2(ispec,0,0,0)) = run_arr(OPS_ACC_MD2(ispec,0,0,0)) + rklhs(nrkstp) * rhs_arr(OPS_ACC_MD3(ispec,0,0,0))
    rhs_arr(OPS_ACC_MD3(ispec,0,0,0)) = run_arr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE fincom_kernel_MD
