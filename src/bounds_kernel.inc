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

!A_yz = A_yz + B_mulditim_yz*C_multidim_yz
SUBROUTINE bounds_kernel_eqA_xdir(sorpx,strhx,ratex,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sorpx
    real(kind=dp), dimension(1), intent(in) :: strhx,ratex
    integer, intent(in) :: ispec

    sorpx(OPS_ACC1(0,0,0)) = sorpx(OPS_ACC1(0,0,0)) + strhx(OPS_ACC_MD2(ispec,0,0,0))*ratex(OPS_ACC_MD3(ispec,0,0,0))

END SUBROUTINE bounds_kernel_eqA_xdir

!A_xz = A_xz + B_mulditim_xz*C_multidim_xz
SUBROUTINE bounds_kernel_eqA_ydir(sorpy,strhy,ratey,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sorpy
    real(kind=dp), dimension(1), intent(in) :: strhy,ratey
    integer, intent(in) :: ispec

    sorpy(OPS_ACC1(0,0,0)) = sorpy(OPS_ACC1(0,0,0)) + strhy(OPS_ACC_MD2(ispec,0,0,0))*ratey(OPS_ACC_MD3(ispec,0,0,0))

END SUBROUTINE bounds_kernel_eqA_ydir

!A_xy = A_xy + B_mulditim_xy*C_multidim_xy
SUBROUTINE bounds_kernel_eqA_zdir(sorpz,strhz,ratez,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sorpz
    real(kind=dp), dimension(1), intent(in) :: strhz,ratez
    integer, intent(in) :: ispec

    sorpz(OPS_ACC1(0,0,0)) = sorpz(OPS_ACC1(0,0,0)) + strhz(OPS_ACC_MD2(ispec,0,0,0))*ratez(OPS_ACC_MD3(ispec,0,0,0))

END SUBROUTINE bounds_kernel_eqA_zdir

!A_yz = A_yz + B_mulditim_yz*val1
SUBROUTINE bounds_kernel_eqB_xdir(sydtx,dydtx,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sydtx
    real(kind=dp), dimension(1), intent(in) :: dydtx
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

    sydtx(OPS_ACC1(0,0,0)) = sydtx(OPS_ACC1(0,0,0)) + dydtx(OPS_ACC_MD2(ispec,0,0,0))*rgspec

END SUBROUTINE bounds_kernel_eqB_xdir

!A_xz = A_xz + B_mulditim_xz*val1
SUBROUTINE bounds_kernel_eqB_ydir(sydty,dydty,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sydty
    real(kind=dp), dimension(1), intent(in) :: dydty
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

    sydty(OPS_ACC1(0,0,0)) = sydty(OPS_ACC1(0,0,0)) + dydty(OPS_ACC_MD2(ispec,0,0,0))*rgspec

END SUBROUTINE bounds_kernel_eqB_ydir

!A_xy = A_xy + B_mulditim_xy*val1
SUBROUTINE bounds_kernel_eqB_zdir(sydtz,dydtz,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sydtz
    real(kind=dp), dimension(1), intent(in) :: dydtz
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

    sydtz(OPS_ACC1(0,0,0)) = sydtz(OPS_ACC1(0,0,0)) + dydtz(OPS_ACC_MD2(ispec,0,0,0))*rgspec

END SUBROUTINE bounds_kernel_eqB_zdir

!A_xy = -A_xy*B_xy
SUBROUTINE bounds_kernel_eqC_xdir(sorpx,gam1x)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sorpx
    real(kind=dp), dimension(1), intent(in) :: gam1x

    sorpx(OPS_ACC1(0,0,0)) = -sorpx(OPS_ACC1(0,0,0))*gam1x(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqC_xdir

!A_xz = -A_xz*B_xz
SUBROUTINE bounds_kernel_eqC_ydir(sorpy,gam1y)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sorpy
    real(kind=dp), dimension(1), intent(in) :: gam1y

    sorpy(OPS_ACC1(0,0,0)) = -sorpy(OPS_ACC1(0,0,0))*gam1y(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqC_ydir

!A_xy = -A_xy*B_xy
SUBROUTINE bounds_kernel_eqC_zdir(sorpz,gam1z)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sorpz
    real(kind=dp), dimension(1), intent(in) :: gam1z

    sorpz(OPS_ACC1(0,0,0)) = -sorpz(OPS_ACC1(0,0,0))*gam1z(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqC_zdir

!A_yz = A_yz/B_yz
SUBROUTINE bounds_kernel_eqD_xdir(sydtx,strrx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sydtx
    real(kind=dp), dimension(1), intent(in) :: strrx

    sydtx(OPS_ACC1(0,0,0)) = sydtx(OPS_ACC1(0,0,0))/strrx(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqD_xdir

!A_xz = A_xz/B_xz
SUBROUTINE bounds_kernel_eqD_ydir(sydty,strry)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sydty
    real(kind=dp), dimension(1), intent(in) :: strry

    sydty(OPS_ACC1(0,0,0)) = sydty(OPS_ACC1(0,0,0))/strry(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqD_ydir

!A_xy = A_xy/B_xy
SUBROUTINE bounds_kernel_eqD_zdir(sydtz,strrz)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: sydtz
    real(kind=dp), dimension(1), intent(in) :: strrz

    sydtz(OPS_ACC1(0,0,0)) = sydtz(OPS_ACC1(0,0,0))/strrz(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqD_zdir
