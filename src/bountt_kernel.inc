SUBROUTINE bountt_kernel_eqA_xdir(urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr,drhs,strux,strvx,strwx)

    real(kind=8), dimension(1) :: urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr
    real(kind=8), dimension(1), intent(in) :: drhs,strux,strvx,strwx

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strux(OPS_ACC11(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strvx(OPS_ACC12(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strwx(OPS_ACC13(0,0,0))

    urun(OPS_ACC4(0,0,0)) = urhs(OPS_ACC1(0,0,0))
    vrun(OPS_ACC5(0,0,0)) = vrhs(OPS_ACC2(0,0,0))
    wrun(OPS_ACC6(0,0,0)) = wrhs(OPS_ACC3(0,0,0))

    uerr(OPS_ACC7(0,0,0)) = 0.0_8
    verr(OPS_ACC8(0,0,0)) = 0.0_8
    werr(OPS_ACC9(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqA_xdir

SUBROUTINE bountt_kernel_eqA_ydir(urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr,drhs,struy,strvy,strwy)

    real(kind=8), dimension(1) :: urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr
    real(kind=8), dimension(1), intent(in) :: drhs,struy,strvy,strwy

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC10(0,0,0))*struy(OPS_ACC11(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strvy(OPS_ACC12(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strwy(OPS_ACC13(0,0,0))

    urun(OPS_ACC4(0,0,0)) = urhs(OPS_ACC1(0,0,0))
    vrun(OPS_ACC5(0,0,0)) = vrhs(OPS_ACC2(0,0,0))
    wrun(OPS_ACC6(0,0,0)) = wrhs(OPS_ACC3(0,0,0))

    uerr(OPS_ACC7(0,0,0)) = 0.0_8
    verr(OPS_ACC8(0,0,0)) = 0.0_8
    werr(OPS_ACC9(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqA_ydir

SUBROUTINE bountt_kernel_eqA_zdir(urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr,drhs,struz,strvz,strwz)

    real(kind=8), dimension(1) :: urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr
    real(kind=8), dimension(1), intent(in) :: drhs,struz,strvz,strwz

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC10(0,0,0))*struz(OPS_ACC11(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strvz(OPS_ACC12(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC10(0,0,0))*strwz(OPS_ACC13(0,0,0))

    urun(OPS_ACC4(0,0,0)) = urhs(OPS_ACC1(0,0,0))
    vrun(OPS_ACC5(0,0,0)) = vrhs(OPS_ACC2(0,0,0))
    wrun(OPS_ACC6(0,0,0)) = wrhs(OPS_ACC3(0,0,0))

    uerr(OPS_ACC7(0,0,0)) = 0.0_8
    verr(OPS_ACC8(0,0,0)) = 0.0_8
    werr(OPS_ACC9(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqA_zdir

SUBROUTINE bountt_kernel_eqB_xdir(erhs,drhs,strux,strvx,strwx)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: drhs,strux,strvx,strwx

    erhs(OPS_ACC1(0,0,0)) = 0.5_8*(strux(OPS_ACC3(0,0,0))*strux(OPS_ACC3(0,0,0))  &
            + strvx(OPS_ACC4(0,0,0))*strvx(OPS_ACC4(0,0,0)) + strwx(OPS_ACC5(0,0,0))*strwx(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE bountt_kernel_eqB_xdir

SUBROUTINE bountt_kernel_eqB_ydir(erhs,drhs,struy,strvy,strwy)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: drhs,struy,strvy,strwy

    erhs(OPS_ACC1(0,0,0)) = 0.5_8*(struy(OPS_ACC3(0,0,0))*struy(OPS_ACC3(0,0,0))  &
            + strvy(OPS_ACC4(0,0,0))*strvy(OPS_ACC4(0,0,0)) + strwy(OPS_ACC5(0,0,0))*strwy(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE bountt_kernel_eqB_ydir

SUBROUTINE bountt_kernel_eqB_zdir(erhs,drhs,struz,strvz,strwz)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: drhs,struz,strvz,strwz

    erhs(OPS_ACC1(0,0,0)) = 0.5_8*(struz(OPS_ACC3(0,0,0))*struz(OPS_ACC3(0,0,0))  &
            + strvz(OPS_ACC4(0,0,0))*strvz(OPS_ACC4(0,0,0)) + strwz(OPS_ACC5(0,0,0))*strwz(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE bountt_kernel_eqB_zdir

SUBROUTINE bountt_kernel_eqC_xdir(drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr,strdx,strux,strvx,strwx)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr
    real(kind=8), dimension(1), intent(in) :: strdx,strux,strvx,strwx

    drhs(OPS_ACC1(0,0,0)) = strdx(OPS_ACC13(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = strdx(OPS_ACC13(0,0,0))*strux(OPS_ACC14(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdx(OPS_ACC13(0,0,0))*strvx(OPS_ACC15(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdx(OPS_ACC13(0,0,0))*strwx(OPS_ACC16(0,0,0))

    drun(OPS_ACC5(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    urun(OPS_ACC6(0,0,0)) = urhs(OPS_ACC2(0,0,0))
    vrun(OPS_ACC7(0,0,0)) = vrhs(OPS_ACC3(0,0,0))
    wrun(OPS_ACC8(0,0,0)) = wrhs(OPS_ACC4(0,0,0))

    derr(OPS_ACC9(0,0,0)) = 0.0_8
    uerr(OPS_ACC10(0,0,0)) = 0.0_8
    verr(OPS_ACC11(0,0,0)) = 0.0_8
    werr(OPS_ACC12(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqC_xdir

SUBROUTINE bountt_kernel_eqC_ydir(drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr,strdy,struy,strvy,strwy)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr
    real(kind=8), dimension(1), intent(in) :: strdy,struy,strvy,strwy

    drhs(OPS_ACC1(0,0,0)) = strdy(OPS_ACC13(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = strdy(OPS_ACC13(0,0,0))*struy(OPS_ACC14(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdy(OPS_ACC13(0,0,0))*strvy(OPS_ACC15(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdy(OPS_ACC13(0,0,0))*strwy(OPS_ACC16(0,0,0))

    drun(OPS_ACC5(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    urun(OPS_ACC6(0,0,0)) = urhs(OPS_ACC2(0,0,0))
    vrun(OPS_ACC7(0,0,0)) = vrhs(OPS_ACC3(0,0,0))
    wrun(OPS_ACC8(0,0,0)) = wrhs(OPS_ACC4(0,0,0))

    derr(OPS_ACC9(0,0,0)) = 0.0_8
    uerr(OPS_ACC10(0,0,0)) = 0.0_8
    verr(OPS_ACC11(0,0,0)) = 0.0_8
    werr(OPS_ACC12(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqC_ydir

SUBROUTINE bountt_kernel_eqC_zdir(drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr,strdz,struz,strvz,strwz)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,drun,urun,vrun,wrun,derr,uerr,verr,werr
    real(kind=8), dimension(1), intent(in) :: strdz,struz,strvz,strwz

    drhs(OPS_ACC1(0,0,0)) = strdz(OPS_ACC13(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = strdz(OPS_ACC13(0,0,0))*struz(OPS_ACC14(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdz(OPS_ACC13(0,0,0))*strvz(OPS_ACC15(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdz(OPS_ACC13(0,0,0))*strwz(OPS_ACC16(0,0,0))

    drun(OPS_ACC5(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    urun(OPS_ACC6(0,0,0)) = urhs(OPS_ACC2(0,0,0))
    vrun(OPS_ACC7(0,0,0)) = vrhs(OPS_ACC3(0,0,0))
    wrun(OPS_ACC8(0,0,0)) = wrhs(OPS_ACC4(0,0,0))

    derr(OPS_ACC9(0,0,0)) = 0.0_8
    uerr(OPS_ACC10(0,0,0)) = 0.0_8
    verr(OPS_ACC11(0,0,0)) = 0.0_8
    werr(OPS_ACC12(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqC_zdir

SUBROUTINE bountt_kernel_eqD(erun,eerr,erhs)

    real(kind=8), dimension(1) :: erun,eerr
    real(kind=8), dimension(1), intent(in) :: erhs

    erun(OPS_ACC1(0,0,0)) = erhs(OPS_ACC3(0,0,0))

    eerr(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqD

SUBROUTINE bountt_kernel_eqE_xdir(yrhs,yrun,yerr,stryxl,strdxl)

    real(kind=8), dimension(1) :: yrhs,yrun,yerr
    real(kind=8), dimension(1), intent(in) :: stryxl,strdxl

    yrhs(OPS_ACC1(0,0,0)) = strdxl(OPS_ACC5(0,0,0))*stryxl(OPS_ACC4(0,0,0))

    yrun(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC1(0,0,0))

    yerr(OPS_ACC3(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqE_xdir

SUBROUTINE bountt_kernel_eqE_ydir(yrhs,yrun,yerr,stryyl,strdyl)

    real(kind=8), dimension(1) :: yrhs,yrun,yerr
    real(kind=8), dimension(1), intent(in) :: stryyl,strdyl

    yrhs(OPS_ACC1(0,0,0)) = strdyl(OPS_ACC5(0,0,0))*stryyl(OPS_ACC4(0,0,0))

    yrun(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC1(0,0,0))

    yerr(OPS_ACC3(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqE_ydir

SUBROUTINE bountt_kernel_eqE_zdir(yrhs,yrun,yerr,stryzl,strdzl)

    real(kind=8), dimension(1) :: yrhs,yrun,yerr
    real(kind=8), dimension(1), intent(in) :: stryzl,strdzl

    yrhs(OPS_ACC1(0,0,0)) = strdzl(OPS_ACC5(0,0,0))*stryzl(OPS_ACC4(0,0,0))

    yrun(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC1(0,0,0))

    yerr(OPS_ACC3(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqE_zdir

SUBROUTINE bountt_kernel_eqF_xdir(erhs,yrhs,yrun,yerr,itndex,drhs,strtx,stryx,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs,yrhs,yrun,yerr
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: drhs,strtx,stryx

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC5(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strtx(OPS_ACC7(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtx(OPS_ACC7(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC6(0,0,0))*stryx(OPS_ACC8(0,0,0))

    yrun(OPS_ACC3(0,0,0)) = yrhs(OPS_ACC2(0,0,0))

    yerr(OPS_ACC4(0,0,0)) = 0.0_8

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strtx(OPS_ACC7(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE bountt_kernel_eqF_xdir

SUBROUTINE bountt_kernel_eqF_ydir(erhs,yrhs,yrun,yerr,itndex,drhs,strty,stryy,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs,yrhs,yrun,yerr
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: drhs,strty,stryy

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC5(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strty(OPS_ACC7(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strty(OPS_ACC7(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC6(0,0,0))*stryy(OPS_ACC8(0,0,0))

    yrun(OPS_ACC3(0,0,0)) = yrhs(OPS_ACC2(0,0,0))

    yerr(OPS_ACC4(0,0,0)) = 0.0_8

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strty(OPS_ACC7(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE bountt_kernel_eqF_ydir

SUBROUTINE bountt_kernel_eqF_zdir(erhs,yrhs,yrun,yerr,itndex,drhs,strtz,stryz,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs,yrhs,yrun,yerr
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: drhs,strtz,stryz

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC5(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strtz(OPS_ACC7(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtz(OPS_ACC7(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC6(0,0,0))*stryz(OPS_ACC8(0,0,0))

    yrun(OPS_ACC3(0,0,0)) = yrhs(OPS_ACC2(0,0,0))

    yerr(OPS_ACC4(0,0,0)) = 0.0_8

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strtz(OPS_ACC7(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE bountt_kernel_eqF_zdir

SUBROUTINE bountt_kernel_eqG_xyz(urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr,erhs,drhs)

    real(kind=8), dimension(1) :: urhs,vrhs,wrhs,urun,vrun,wrun,uerr,verr,werr,erhs
    real(kind=8), dimension(1), intent(in) :: drhs

    urhs(OPS_ACC1(0,0,0)) = 0.0_8
    vrhs(OPS_ACC2(0,0,0)) = 0.0_8
    wrhs(OPS_ACC3(0,0,0)) = 0.0_8

    urun(OPS_ACC4(0,0,0)) = 0.0_8
    vrun(OPS_ACC5(0,0,0)) = 0.0_8
    wrun(OPS_ACC6(0,0,0)) = 0.0_8

    uerr(OPS_ACC7(0,0,0)) = 0.0_8
    verr(OPS_ACC8(0,0,0)) = 0.0_8
    werr(OPS_ACC9(0,0,0)) = 0.0_8

    erhs(OPS_ACC10(0,0,0)) = 0.0_8
    erhs(OPS_ACC10(0,0,0)) = drhs(OPS_ACC11(0,0,0))*erhs(OPS_ACC10(0,0,0))

END SUBROUTINE bountt_kernel_eqG_xyz

SUBROUTINE bountt_kernel_eqH_xdir(strtx,dtdtx,val)

    real(kind=8), dimension(1) :: strtx,dtdtx
    real(kind=8), intent(in) :: val

    strtx(OPS_ACC1(0,0,0)) = val
    dtdtx(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqH_xdir

SUBROUTINE bountt_kernel_eqH_ydir(strty,dtdty,val)

    real(kind=8), dimension(1) :: strty,dtdty
    real(kind=8), intent(in) :: val

    strty(OPS_ACC1(0,0,0)) = val
    dtdty(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqH_ydir

SUBROUTINE bountt_kernel_eqH_zdir(strtz,dtdtz,val)

    real(kind=8), dimension(1) :: strtz,dtdtz
    real(kind=8), intent(in) :: val

    strtz(OPS_ACC1(0,0,0)) = val
    dtdtz(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bountt_kernel_eqH_zdir

