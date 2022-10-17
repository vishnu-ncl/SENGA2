SUBROUTINE adaptt_kernel_err_eval(err_arr, run_arr, ernrm, ertot)
    use data_types
    implicit none
 
    real(kind=dp), dimension(1), intent(in) :: err_arr, run_arr
    real(kind=dp), intent(in) :: ernrm
    real(kind=8) :: ertot
    real(kind=dp) :: fornow

    fornow = abs(err_arr(OPS_ACC1(0,0,0)))/(abs(run_arr(OPS_ACC2(0,0,0)))+ernrm)
    ertot = max( ertot, fornow )

END SUBROUTINE adaptt_kernel_err_eval

SUBROUTINE adaptt_kernel_err_eval_multidim(err_arr, run_arr, ernrm, ispec, ertot)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: err_arr, run_arr
    real(kind=dp), intent(in) :: ernrm
    real(kind=8) :: ertot
    integer :: ispec
    real(kind=dp) :: fornow

    fornow = abs(err_arr(OPS_ACC_MD1(ispec,0,0,0)))/(abs(run_arr(OPS_ACC_MD2(ispec,0,0,0)))+ernrm)
    ertot = max( ertot, fornow )

END SUBROUTINE adaptt_kernel_err_eval_multidim