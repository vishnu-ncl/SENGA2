SUBROUTINE lincom_kernel_main(err_arr,run_arr,rhs_arr,rkerr_irkstp,rklhs_irkstp,rkrhs_irkstp)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: err_arr, run_arr, rhs_arr
    real(kind=dp), intent(in) :: rkerr_irkstp, rklhs_irkstp, rkrhs_irkstp
    real(kind=dp) :: fornow

    err_arr(OPS_ACC1(0,0,0)) = err_arr(OPS_ACC1(0,0,0)) + rkerr_irkstp * rhs_arr(OPS_ACC3(0,0,0))

    fornow = run_arr(OPS_ACC2(0,0,0))
    run_arr(OPS_ACC2(0,0,0)) = fornow + rklhs_irkstp * rhs_arr(OPS_ACC3(0,0,0))
    rhs_arr(OPS_ACC3(0,0,0)) = fornow + rkrhs_irkstp * rhs_arr(OPS_ACC3(0,0,0))

END SUBROUTINE lincom_kernel_main

SUBROUTINE lincom_kernel_multidim(err_arr,run_arr,rhs_arr,rkerr_irkstp,rklhs_irkstp,rkrhs_irkstp,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: err_arr, run_arr, rhs_arr
    real(kind=dp), intent(in) :: rkerr_irkstp, rklhs_irkstp, rkrhs_irkstp
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    err_arr(OPS_ACC_MD1(ispec,0,0,0)) = err_arr(OPS_ACC_MD1(ispec,0,0,0)) + rkerr_irkstp * rhs_arr(OPS_ACC_MD3(ispec,0,0,0))

    fornow = run_arr(OPS_ACC_MD2(ispec,0,0,0))
    run_arr(OPS_ACC_MD2(ispec,0,0,0)) = fornow + rklhs_irkstp * rhs_arr(OPS_ACC_MD3(ispec,0,0,0))
    rhs_arr(OPS_ACC_MD3(ispec,0,0,0)) = fornow + rkrhs_irkstp * rhs_arr(OPS_ACC_MD3(ispec,0,0,0))

END SUBROUTINE lincom_kernel_multidim
