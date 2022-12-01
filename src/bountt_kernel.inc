SUBROUTINE bountt_kernel_eqA_xdir(urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr,drhs,strux,strvx,strwx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr
    real(kind=dp), dimension(1), intent(in) :: drhs,strux,strvx,strwx

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strux(OPS_ACC11(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strvx(OPS_ACC12(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strwx(OPS_ACC13(0,0,0))

    urun(OPS_ACC4(0,0,0)) = urhs(OPS_ACC1(0,0,0))
    vrun(OPS_ACC5(0,0,0)) = vrhs(OPS_ACC2(0,0,0))
    wrun(OPS_ACC6(0,0,0)) = wrhs(OPS_ACC3(0,0,0))

    uerr(OPS_ACC7(0,0,0)) = 0.0_dp
    verr(OPS_ACC8(0,0,0)) = 0.0_dp
    werr(OPS_ACC9(0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqA_xdir

SUBROUTINE bountt_kernel_eqA_ydir(urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr,drhs,struy,strvy,strwy)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr
    real(kind=dp), dimension(1), intent(in) :: drhs,struy,strvy,strwy

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC10(0,0,0))*struy(OPS_ACC11(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strvy(OPS_ACC12(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strwy(OPS_ACC13(0,0,0))

    urun(OPS_ACC4(0,0,0)) = urhs(OPS_ACC1(0,0,0))
    vrun(OPS_ACC5(0,0,0)) = vrhs(OPS_ACC2(0,0,0))
    wrun(OPS_ACC6(0,0,0)) = wrhs(OPS_ACC3(0,0,0))

    uerr(OPS_ACC7(0,0,0)) = 0.0_dp
    verr(OPS_ACC8(0,0,0)) = 0.0_dp
    werr(OPS_ACC9(0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqA_ydir

SUBROUTINE bountt_kernel_eqA_zdir(urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr,drhs,struz,strvz,strwz)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr
    real(kind=dp), dimension(1), intent(in) :: drhs,struz,strvz,strwz

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC10(0,0,0))*struz(OPS_ACC11(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strvz(OPS_ACC12(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strwz(OPS_ACC13(0,0,0))

    urun(OPS_ACC4(0,0,0)) = urhs(OPS_ACC1(0,0,0))
    vrun(OPS_ACC5(0,0,0)) = vrhs(OPS_ACC2(0,0,0))
    wrun(OPS_ACC6(0,0,0)) = wrhs(OPS_ACC3(0,0,0))

    uerr(OPS_ACC7(0,0,0)) = 0.0_dp
    verr(OPS_ACC8(0,0,0)) = 0.0_dp
    werr(OPS_ACC9(0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqA_zdir

SUBROUTINE bountt_kernel_eqB_xdir(erhs,drhs,strux,strvx,strwx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: drhs,strux,strvx,strwx

    erhs(OPS_ACC1(0,0,0)) = 0.5_dp*(strux(OPS_ACC3(0,0,0))*strux(OPS_ACC3(0,0,0))  &
            + strvx(OPS_ACC4(0,0,0))*strvx(OPS_ACC4(0,0,0)) + strwx(OPS_ACC5(0,0,0))*strwx(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE bountt_kernel_eqB_xdir

SUBROUTINE bountt_kernel_eqB_ydir(erhs,drhs,struy,strvy,strwy)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: drhs,struy,strvy,strwy

    erhs(OPS_ACC1(0,0,0)) = 0.5_dp*(struy(OPS_ACC3(0,0,0))*struy(OPS_ACC3(0,0,0))  &
            + strvy(OPS_ACC4(0,0,0))*strvy(OPS_ACC4(0,0,0)) + strwy(OPS_ACC5(0,0,0))*strwy(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE bountt_kernel_eqB_ydir

SUBROUTINE bountt_kernel_eqB_zdir(erhs,drhs,struz,strvz,strwz)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: drhs,struz,strvz,strwz

    erhs(OPS_ACC1(0,0,0)) = 0.5_dp*(struz(OPS_ACC3(0,0,0))*struz(OPS_ACC3(0,0,0))  &
            + strvz(OPS_ACC4(0,0,0))*strvz(OPS_ACC4(0,0,0)) + strwz(OPS_ACC5(0,0,0))*strwz(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE bountt_kernel_eqB_zdir

SUBROUTINE bountt_kernel_eqC_xdir(drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr,strdx,strux,strvx,strwx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr
    real(kind=dp), dimension(1), intent(in) :: strdx,strux,strvx,strwx

    drhs(OPS_ACC1(0,0,0)) = strdx(OPS_ACC13(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = strdx(OPS_ACC13(0,0,0))*strux(OPS_ACC14(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdx(OPS_ACC13(0,0,0))*strvx(OPS_ACC15(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdx(OPS_ACC13(0,0,0))*strwx(OPS_ACC16(0,0,0))

    drun(OPS_ACC5(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    urun(OPS_ACC6(0,0,0)) = urhs(OPS_ACC2(0,0,0))
    vrun(OPS_ACC7(0,0,0)) = vrhs(OPS_ACC3(0,0,0))
    wrun(OPS_ACC8(0,0,0)) = wrhs(OPS_ACC4(0,0,0))

    derr(OPS_ACC9(0,0,0)) = 0.0_dp
    uerr(OPS_ACC10(0,0,0)) = 0.0_dp
    verr(OPS_ACC11(0,0,0)) = 0.0_dp
    werr(OPS_ACC12(0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqC_xdir

SUBROUTINE bountt_kernel_eqC_ydir(drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr,strdy,struy,strvy,strwy)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr
    real(kind=dp), dimension(1), intent(in) :: strdy,struy,strvy,strwy

    drhs(OPS_ACC1(0,0,0)) = strdy(OPS_ACC13(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = strdy(OPS_ACC13(0,0,0))*struy(OPS_ACC14(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdy(OPS_ACC13(0,0,0))*strvy(OPS_ACC15(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdy(OPS_ACC13(0,0,0))*strwy(OPS_ACC16(0,0,0))

    drun(OPS_ACC5(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    urun(OPS_ACC6(0,0,0)) = urhs(OPS_ACC2(0,0,0))
    vrun(OPS_ACC7(0,0,0)) = vrhs(OPS_ACC3(0,0,0))
    wrun(OPS_ACC8(0,0,0)) = wrhs(OPS_ACC4(0,0,0))

    derr(OPS_ACC9(0,0,0)) = 0.0_dp
    uerr(OPS_ACC10(0,0,0)) = 0.0_dp
    verr(OPS_ACC11(0,0,0)) = 0.0_dp
    werr(OPS_ACC12(0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqC_ydir

SUBROUTINE bountt_kernel_eqC_zdir(drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr,strdz,struz,strvz,strwz)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr
    real(kind=dp), dimension(1), intent(in) :: strdz,struz,strvz,strwz

    drhs(OPS_ACC1(0,0,0)) = strdz(OPS_ACC13(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = strdz(OPS_ACC13(0,0,0))*struz(OPS_ACC14(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdz(OPS_ACC13(0,0,0))*strvz(OPS_ACC15(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdz(OPS_ACC13(0,0,0))*strwz(OPS_ACC16(0,0,0))

    drun(OPS_ACC5(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    urun(OPS_ACC6(0,0,0)) = urhs(OPS_ACC2(0,0,0))
    vrun(OPS_ACC7(0,0,0)) = vrhs(OPS_ACC3(0,0,0))
    wrun(OPS_ACC8(0,0,0)) = wrhs(OPS_ACC4(0,0,0))

    derr(OPS_ACC9(0,0,0)) = 0.0_dp
    uerr(OPS_ACC10(0,0,0)) = 0.0_dp
    verr(OPS_ACC11(0,0,0)) = 0.0_dp
    werr(OPS_ACC12(0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqC_zdir

SUBROUTINE bountt_kernel_eqD(erun,eerr,erhs)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erun,eerr
    real(kind=dp), dimension(1), intent(in) :: erhs

    erun(OPS_ACC1(0,0,0)) = erhs(OPS_ACC3(0,0,0))

    eerr(OPS_ACC2(0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqD

SUBROUTINE bountt_kernel_eqE_xdir(yrhs,yrun,yerr,stryxl,strdxl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs,yrun,yerr
    real(kind=dp), dimension(1), intent(in) :: stryxl,strdxl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = strdxl(OPS_ACC5(0,0,0))*stryxl(OPS_ACC_MD4(ispec,0,0,0))

    yrun(OPS_ACC_MD2(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))

    yerr(OPS_ACC_MD3(ispec,0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqE_xdir

SUBROUTINE bountt_kernel_eqE_ydir(yrhs,yrun,yerr,stryyl,strdyl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs,yrun,yerr
    real(kind=dp), dimension(1), intent(in) :: stryyl,strdyl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = strdyl(OPS_ACC5(0,0,0))*stryyl(OPS_ACC_MD4(ispec,0,0,0))

    yrun(OPS_ACC_MD2(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))

    yerr(OPS_ACC_MD3(ispec,0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqE_ydir

SUBROUTINE bountt_kernel_eqE_zdir(yrhs,yrun,yerr,stryzl,strdzl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs,yrun,yerr
    real(kind=dp), dimension(1), intent(in) :: stryzl,strdzl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = strdzl(OPS_ACC5(0,0,0))*stryzl(OPS_ACC_MD4(ispec,0,0,0))

    yrun(OPS_ACC_MD2(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))

    yerr(OPS_ACC_MD3(ispec,0,0,0)) = 0.0_dp

END SUBROUTINE bountt_kernel_eqE_zdir
    