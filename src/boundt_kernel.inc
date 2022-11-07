SUBROUTINE boundt_kernel_eqA_xdir(urhs,vrhs,wrhs,drhs,strux,strvx,strwx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: urhs,vrhs,wrhs
    real(kind=dp), dimension(1), intent(in) :: drhs,strux,strvx,strwx

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strux(OPS_ACC5(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strvx(OPS_ACC6(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strwx(OPS_ACC7(0,0,0))

END SUBROUTINE boundt_kernel_eqA_xdir

SUBROUTINE boundt_kernel_eqA_ydir(urhs,vrhs,wrhs,drhs,struy,strvy,strwy)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: urhs,vrhs,wrhs
    real(kind=dp), dimension(1), intent(in) :: drhs,struy,strvy,strwy

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC4(0,0,0))*struy(OPS_ACC5(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strvy(OPS_ACC6(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strwy(OPS_ACC7(0,0,0))

END SUBROUTINE boundt_kernel_eqA_ydir

SUBROUTINE boundt_kernel_eqA_zdir(urhs,vrhs,wrhs,drhs,struz,strvz,strwz)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: urhs,vrhs,wrhs
    real(kind=dp), dimension(1), intent(in) :: drhs,struz,strvz,strwz

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC4(0,0,0))*struz(OPS_ACC5(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strvz(OPS_ACC6(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strwz(OPS_ACC7(0,0,0))

END SUBROUTINE boundt_kernel_eqA_zdir

SUBROUTINE boundt_kernel_eqB_xdir(erhs,drhs,strux,strvx,strwx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: drhs,strux,strvx,strwx

    erhs(OPS_ACC1(0,0,0)) = 0.5_dp*(strux(OPS_ACC3(0,0,0))*strux(OPS_ACC3(0,0,0))  &
            + strvx(OPS_ACC4(0,0,0))*strvx(OPS_ACC4(0,0,0)) + strwx(OPS_ACC5(0,0,0))*strwx(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundt_kernel_eqB_xdir

SUBROUTINE boundt_kernel_eqB_ydir(erhs,drhs,struy,strvy,strwy)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: drhs,struy,strvy,strwy

    erhs(OPS_ACC1(0,0,0)) = 0.5_dp*(struy(OPS_ACC3(0,0,0))*struy(OPS_ACC3(0,0,0))  &
            + strvy(OPS_ACC4(0,0,0))*strvy(OPS_ACC4(0,0,0)) + strwy(OPS_ACC5(0,0,0))*strwy(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundt_kernel_eqB_ydir

SUBROUTINE boundt_kernel_eqB_zdir(erhs,drhs,struz,strvz,strwz)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: drhs,struz,strvz,strwz

    erhs(OPS_ACC1(0,0,0)) = 0.5_dp*(struz(OPS_ACC3(0,0,0))*struz(OPS_ACC3(0,0,0))  &
            + strvz(OPS_ACC4(0,0,0))*strvz(OPS_ACC4(0,0,0)) + strwz(OPS_ACC5(0,0,0))*strwz(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundt_kernel_eqB_zdir

SUBROUTINE boundt_kernel_eqC_xdir(drhs,urhs,vrhs,wrhs,strdx,strux,strvx,strwx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs
    real(kind=dp), dimension(1), intent(in) :: strdx,strux,strvx,strwx

    drhs(OPS_ACC1(0,0,0)) = strdx(OPS_ACC5(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = strdx(OPS_ACC5(0,0,0))*strux(OPS_ACC6(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdx(OPS_ACC5(0,0,0))*strvx(OPS_ACC7(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdx(OPS_ACC5(0,0,0))*strwx(OPS_ACC8(0,0,0))

END SUBROUTINE boundt_kernel_eqC_xdir

SUBROUTINE boundt_kernel_eqC_ydir(drhs,urhs,vrhs,wrhs,strdy,struy,strvy,strwy)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs
    real(kind=dp), dimension(1), intent(in) :: strdy,struy,strvy,strwy

    drhs(OPS_ACC1(0,0,0)) = strdy(OPS_ACC5(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = strdy(OPS_ACC5(0,0,0))*struy(OPS_ACC6(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdy(OPS_ACC5(0,0,0))*strvy(OPS_ACC7(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdy(OPS_ACC5(0,0,0))*strwy(OPS_ACC8(0,0,0))

END SUBROUTINE boundt_kernel_eqC_ydir

SUBROUTINE boundt_kernel_eqC_zdir(drhs,urhs,vrhs,wrhs,strdz,struz,strvz,strwz)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs
    real(kind=dp), dimension(1), intent(in) :: strdz,struz,strvz,strwz

    drhs(OPS_ACC1(0,0,0)) = strdz(OPS_ACC5(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = strdz(OPS_ACC5(0,0,0))*struz(OPS_ACC6(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdz(OPS_ACC5(0,0,0))*strvz(OPS_ACC7(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdz(OPS_ACC5(0,0,0))*strwz(OPS_ACC8(0,0,0))

END SUBROUTINE boundt_kernel_eqC_zdir

SUBROUTINE boundt_kernel_eqD_xdir(yrhs,stryx,strdx,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryx,strdx
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = strdx(OPS_ACC3(0,0,0))*stryx(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE boundt_kernel_eqD_xdir

SUBROUTINE boundt_kernel_eqD_ydir(yrhs,stryy,strdy,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryy,strdy
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = strdy(OPS_ACC3(0,0,0))*stryy(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE boundt_kernel_eqD_ydir

SUBROUTINE boundt_kernel_eqD_zdir(yrhs,stryz,strdz,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryz,strdz
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = strdz(OPS_ACC3(0,0,0))*stryz(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE boundt_kernel_eqD_zdir
