SUBROUTINE adaptt_kernel_err_eval(err_arr, run_arr, ernrm, ertot)

    real(kind=8), dimension(1), intent(in) :: err_arr, run_arr
    real(kind=8), intent(in) :: ernrm
    real(kind=8) :: ertot
    real(kind=8) :: fornow

    fornow = abs(err_arr(OPS_ACC1(0,0,0)))/(abs(run_arr(OPS_ACC2(0,0,0)))+ernrm)
    ertot = max( ertot, fornow )

END SUBROUTINE adaptt_kernel_err_eval

SUBROUTINE adaptt_kernel_err_eval_MD(err_arr, run_arr, ernrm, ispec, ertot)

    real(kind=8), dimension(1), intent(in) :: err_arr, run_arr
    real(kind=8) :: ertot
    integer(kind=4) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: ernrm
    real(kind=8) :: fornow

    fornow = abs(err_arr(OPS_ACC1(0,0,0)))/(abs(run_arr(OPS_ACC2(0,0,0)))+ernrm(ispec))
    ertot = max( ertot, fornow )

END SUBROUTINE adaptt_kernel_err_eval_MD
