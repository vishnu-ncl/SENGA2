SUBROUTINE adaptt_kernel_err_eval(err_arr, run_arr, ernrm, ertot)
    use data_types
    implicit none
 
    real(8), dimension(1), intent(in) :: err_arr, run_arr
    real(8), intent(in) :: ernrm
    real(8) :: ertot
    real(8) :: fornow

    fornow = abs(err_arr(OPS_ACC1(0,0,0)))/(abs(run_arr(OPS_ACC2(0,0,0)))+ernrm)
    ertot = max( ertot, fornow )

END SUBROUTINE adaptt_kernel_err_eval

SUBROUTINE adaptt_kernel_err_eval_MD(err_arr, run_arr, ernrm, ispec, ertot)
    use data_types
    implicit none

    real(8), dimension(1), intent(in) :: err_arr, run_arr
    real(8) :: ertot
    integer(4) :: ispec
    real(8), dimension(nspcmx), intent(in) :: ernrm
    real(8) :: fornow

    fornow = abs(err_arr(OPS_ACC1(0,0,0)))/(abs(run_arr(OPS_ACC2(0,0,0)))+ernrm(ispec))
    ertot = max( ertot, fornow )

END SUBROUTINE adaptt_kernel_err_eval_MD
