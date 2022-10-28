SUBROUTINE bounds_kernel_reduced_enthalpy_xdir(strhx,strgx,strtx,strrx,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strhx
    real(kind=dp), dimension(1), intent(in) :: strgx,strtx,strrx
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

    strhx(OPS_ACC_MD1(ispec,0,0,0)) = strhx(OPS_ACC_MD1(ispec,0,0,0))  &
            - strgx(OPS_ACC2(0,0,0))*strtx(OPS_ACC3(0,0,0))*rgspec/strrx(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_reduced_enthalpy_xdir
