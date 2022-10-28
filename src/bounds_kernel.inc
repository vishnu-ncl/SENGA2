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

SUBROUTINE bounds_kernel_reduced_enthalpy_ydir(strhy,strgy,strty,strry,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strhy
    real(kind=dp), dimension(1), intent(in) :: strgy,strty,strry
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

    strhy(OPS_ACC_MD1(ispec,0,0,0)) = strhy(OPS_ACC_MD1(ispec,0,0,0))  &
            - strgy(OPS_ACC2(0,0,0))*strty(OPS_ACC3(0,0,0))*rgspec/strry(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_reduced_enthalpy_ydir

SUBROUTINE bounds_kernel_reduced_enthalpy_zdir(strhz,strgz,strtz,strrz,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strhz
    real(kind=dp), dimension(1), intent(in) :: strgz,strtz,strrz
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

    strhz(OPS_ACC_MD1(ispec,0,0,0)) = strhz(OPS_ACC_MD1(ispec,0,0,0))  &
            - strgz(OPS_ACC2(0,0,0))*strtz(OPS_ACC3(0,0,0))*rgspec/strrz(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_reduced_enthalpy_zdir

SUBROUTINE bounds_kernel_reduced_energy_xdir(gam1x,strex,ovgmx,strgx,strrx,strtx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: gam1x,strex,ovgmx
    real(kind=dp), dimension(1), intent(in) :: strgx,strrx,strtx

    gam1x(OPS_ACC1(0,0,0)) = strgx(OPS_ACC4(0,0,0)) - strrx(OPS_ACC5(0,0,0))
    strex(OPS_ACC2(0,0,0)) = strex(OPS_ACC2(0,0,0)) - gam1x(OPS_ACC1(0,0,0))*strtx(OPS_ACC6(0,0,0))

    gam1x(OPS_ACC1(0,0,0)) = strrx(OPS_ACC5(0,0,0))/gam1x(OPS_ACC1(0,0,0))
    ovgmx(OPS_ACC3(0,0,0)) = 1.0_dp/gam1x(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_reduced_energy_xdir

SUBROUTINE bounds_kernel_reduced_energy_ydir(gam1y,strey,ovgmy,strgy,strry,strty)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: gam1y,strey,ovgmy
    real(kind=dp), dimension(1), intent(in) :: strgy,strry,strty

    gam1y(OPS_ACC1(0,0,0)) = strgy(OPS_ACC4(0,0,0)) - strry(OPS_ACC5(0,0,0))
    strey(OPS_ACC2(0,0,0)) = strey(OPS_ACC2(0,0,0)) - gam1y(OPS_ACC1(0,0,0))*strty(OPS_ACC6(0,0,0))

    gam1y(OPS_ACC1(0,0,0)) = strry(OPS_ACC5(0,0,0))/gam1y(OPS_ACC1(0,0,0))
    ovgmy(OPS_ACC3(0,0,0)) = 1.0_dp/gam1y(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_reduced_energy_ydir

SUBROUTINE bounds_kernel_reduced_energy_zdir(gam1z,strez,ovgmz,strgz,strrz,strtz)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: gam1z,strez,ovgmz
    real(kind=dp), dimension(1), intent(in) :: strgz,strrz,strtz

    gam1z(OPS_ACC1(0,0,0)) = strgz(OPS_ACC4(0,0,0)) - strrz(OPS_ACC5(0,0,0))
    strez(OPS_ACC2(0,0,0)) = strez(OPS_ACC2(0,0,0)) - gam1z(OPS_ACC1(0,0,0))*strtz(OPS_ACC6(0,0,0))

    gam1z(OPS_ACC1(0,0,0)) = strrz(OPS_ACC5(0,0,0))/gam1z(OPS_ACC1(0,0,0))
    ovgmz(OPS_ACC3(0,0,0)) = 1.0_dp/gam1z(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_reduced_energy_zdir

