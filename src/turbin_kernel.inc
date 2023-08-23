SUBROUTINE turbin_kernel_eqA(urun, utmp, vrun, vtmp, wrun, wtmp, tket)

    real(kind=8), dimension(1), intent(in) :: urun, utmp, vrun, vtmp, wrun, wtmp
    real(kind=8) :: tket

    tket = tket + urun(OPS_ACC1(0,0,0))*urun(OPS_ACC1(0,0,0))  &
          + utmp(OPS_ACC2(0,0,0))*utmp(OPS_ACC2(0,0,0)) + vrun(OPS_ACC3(0,0,0))*vrun(OPS_ACC3(0,0,0))  &
          + vtmp(OPS_ACC4(0,0,0))*vtmp(OPS_ACC4(0,0,0)) + wrun(OPS_ACC5(0,0,0))*wrun(OPS_ACC5(0,0,0))  &
          + wtmp(OPS_ACC6(0,0,0))*wtmp(OPS_ACC6(0,0,0))

END SUBROUTINE turbin_kernel_eqA

SUBROUTINE turbin_kernel_eqB(in_arr, reduction_var)

    real(kind=8), dimension(1), intent(in) :: in_arr
    real(kind=8) :: reduction_var

    reduction_var = reduction_var + in_arr(OPS_ACC1(0,0,0))

END SUBROUTINE turbin_kernel_eqB

SUBROUTINE turbin_kernel_eqC(in_arr, val1, reduction_var)

    real(kind=8), dimension(1), intent(in) :: in_arr
    real(kind=8), intent(in) :: val1
    real(kind=8) :: reduction_var
    real(kind=8) :: fornow

    fornow = in_arr(OPS_ACC1(0,0,0)) - val1

    reduction_var = reduction_var + fornow*fornow

END SUBROUTINE turbin_kernel_eqC
