SUBROUTINE boundt_kernel_eqA_xdir(urhs,vrhs,wrhs,drhs,strux,strvx,strwx)

    real(kind=8), dimension(1) :: urhs,vrhs,wrhs
    real(kind=8), dimension(1), intent(in) :: drhs,strux,strvx,strwx

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strux(OPS_ACC5(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strvx(OPS_ACC6(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strwx(OPS_ACC7(0,0,0))

END SUBROUTINE boundt_kernel_eqA_xdir

SUBROUTINE boundt_kernel_eqA_ydir(urhs,vrhs,wrhs,drhs,struy,strvy,strwy)

    real(kind=8), dimension(1) :: urhs,vrhs,wrhs
    real(kind=8), dimension(1), intent(in) :: drhs,struy,strvy,strwy

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC4(0,0,0))*struy(OPS_ACC5(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strvy(OPS_ACC6(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strwy(OPS_ACC7(0,0,0))

END SUBROUTINE boundt_kernel_eqA_ydir

SUBROUTINE boundt_kernel_eqA_zdir(urhs,vrhs,wrhs,drhs,struz,strvz,strwz)

    real(kind=8), dimension(1) :: urhs,vrhs,wrhs
    real(kind=8), dimension(1), intent(in) :: drhs,struz,strvz,strwz

    urhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC4(0,0,0))*struz(OPS_ACC5(0,0,0))
    vrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strvz(OPS_ACC6(0,0,0))
    wrhs(OPS_ACC3(0,0,0)) = drhs(OPS_ACC4(0,0,0))*strwz(OPS_ACC7(0,0,0))

END SUBROUTINE boundt_kernel_eqA_zdir

SUBROUTINE boundt_kernel_eqB_xdir(erhs,drhs,strux,strvx,strwx)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: drhs,strux,strvx,strwx

    erhs(OPS_ACC1(0,0,0)) = 0.5_8*(strux(OPS_ACC3(0,0,0))*strux(OPS_ACC3(0,0,0))  &
            + strvx(OPS_ACC4(0,0,0))*strvx(OPS_ACC4(0,0,0)) + strwx(OPS_ACC5(0,0,0))*strwx(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundt_kernel_eqB_xdir

SUBROUTINE boundt_kernel_eqB_ydir(erhs,drhs,struy,strvy,strwy)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: drhs,struy,strvy,strwy

    erhs(OPS_ACC1(0,0,0)) = 0.5_8*(struy(OPS_ACC3(0,0,0))*struy(OPS_ACC3(0,0,0))  &
            + strvy(OPS_ACC4(0,0,0))*strvy(OPS_ACC4(0,0,0)) + strwy(OPS_ACC5(0,0,0))*strwy(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundt_kernel_eqB_ydir

SUBROUTINE boundt_kernel_eqB_zdir(erhs,drhs,struz,strvz,strwz)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: drhs,struz,strvz,strwz

    erhs(OPS_ACC1(0,0,0)) = 0.5_8*(struz(OPS_ACC3(0,0,0))*struz(OPS_ACC3(0,0,0))  &
            + strvz(OPS_ACC4(0,0,0))*strvz(OPS_ACC4(0,0,0)) + strwz(OPS_ACC5(0,0,0))*strwz(OPS_ACC5(0,0,0)))

    erhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC2(0,0,0))*erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundt_kernel_eqB_zdir

SUBROUTINE boundt_kernel_eqC_xdir(drhs,urhs,vrhs,wrhs,strdx,strux,strvx,strwx)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs
    real(kind=8), dimension(1), intent(in) :: strdx,strux,strvx,strwx

    drhs(OPS_ACC1(0,0,0)) = strdx(OPS_ACC5(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = strdx(OPS_ACC5(0,0,0))*strux(OPS_ACC6(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdx(OPS_ACC5(0,0,0))*strvx(OPS_ACC7(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdx(OPS_ACC5(0,0,0))*strwx(OPS_ACC8(0,0,0))

END SUBROUTINE boundt_kernel_eqC_xdir

SUBROUTINE boundt_kernel_eqC_ydir(drhs,urhs,vrhs,wrhs,strdy,struy,strvy,strwy)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs
    real(kind=8), dimension(1), intent(in) :: strdy,struy,strvy,strwy

    drhs(OPS_ACC1(0,0,0)) = strdy(OPS_ACC5(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = strdy(OPS_ACC5(0,0,0))*struy(OPS_ACC6(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdy(OPS_ACC5(0,0,0))*strvy(OPS_ACC7(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdy(OPS_ACC5(0,0,0))*strwy(OPS_ACC8(0,0,0))

END SUBROUTINE boundt_kernel_eqC_ydir

SUBROUTINE boundt_kernel_eqC_zdir(drhs,urhs,vrhs,wrhs,strdz,struz,strvz,strwz)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs
    real(kind=8), dimension(1), intent(in) :: strdz,struz,strvz,strwz

    drhs(OPS_ACC1(0,0,0)) = strdz(OPS_ACC5(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = strdz(OPS_ACC5(0,0,0))*struz(OPS_ACC6(0,0,0))
    vrhs(OPS_ACC3(0,0,0)) = strdz(OPS_ACC5(0,0,0))*strvz(OPS_ACC7(0,0,0))
    wrhs(OPS_ACC4(0,0,0)) = strdz(OPS_ACC5(0,0,0))*strwz(OPS_ACC8(0,0,0))

END SUBROUTINE boundt_kernel_eqC_zdir

SUBROUTINE boundt_kernel_eqD_xdir(yrhs,stryx,strdx)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryx,strdx

    yrhs(OPS_ACC1(0,0,0)) = strdx(OPS_ACC3(0,0,0))*stryx(OPS_ACC2(0,0,0))

END SUBROUTINE boundt_kernel_eqD_xdir

SUBROUTINE boundt_kernel_eqD_ydir(yrhs,stryy,strdy)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryy,strdy

    yrhs(OPS_ACC1(0,0,0)) = strdy(OPS_ACC3(0,0,0))*stryy(OPS_ACC2(0,0,0))

END SUBROUTINE boundt_kernel_eqD_ydir

SUBROUTINE boundt_kernel_eqD_zdir(yrhs,stryz,strdz)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryz,strdz

    yrhs(OPS_ACC1(0,0,0)) = strdz(OPS_ACC3(0,0,0))*stryz(OPS_ACC2(0,0,0))

END SUBROUTINE boundt_kernel_eqD_zdir

SUBROUTINE boundt_kernel_eqE_xdir(strtx,itndex,tinthi,ntint,ipower,ispec)

    integer(kind=4), dimension(1) :: itndex
    real(kind=8), dimension(1), intent(in) :: strtx
    integer(kind=4), intent(in) :: ipower, ispec

    real(kind=8), dimension(ntinmx,nspcmx), intent(in) :: tinthi
    integer(kind=4), dimension(nspcmx), intent(in) :: ntint
    integer(kind=4) :: itint

    itint = 1
    DO WHILE (strtx(OPS_ACC1(0,0,0)) > tinthi(itint,ispec) .and. itint < ntint(ispec))
        itint = itint + 1
    END DO

    itndex(OPS_ACC2(0,0,0)) = itndex(OPS_ACC2(0,0,0)) + (itint-1)*ntbase**ipower

END SUBROUTINE boundt_kernel_eqE_xdir

SUBROUTINE boundt_kernel_eqE_ydir(strty,itndex,tinthi,ntint,ipower,ispec)

    integer(kind=4), dimension(1) :: itndex
    real(kind=8), dimension(1), intent(in) :: strty
    integer(kind=4), intent(in) :: ipower, ispec

    real(kind=8), dimension(ntinmx,nspcmx), intent(in) :: tinthi
    integer(kind=4), dimension(nspcmx), intent(in) :: ntint
    integer(kind=4) :: itint

    itint = 1
    DO WHILE (strty(OPS_ACC1(0,0,0)) > tinthi(itint,ispec) .and. itint < ntint(ispec))
        itint = itint + 1
    END DO

    itndex(OPS_ACC2(0,0,0)) = itndex(OPS_ACC2(0,0,0)) + (itint-1)*ntbase**ipower

END SUBROUTINE boundt_kernel_eqE_ydir

SUBROUTINE boundt_kernel_eqE_zdir(strtz,itndex,tinthi,ntint,ipower,ispec)

    integer(kind=4), dimension(1) :: itndex
    real(kind=8), dimension(1), intent(in) :: strtz
    integer(kind=4), intent(in) :: ipower, ispec

    real(kind=8), dimension(ntinmx,nspcmx), intent(in) :: tinthi
    integer(kind=4), dimension(nspcmx), intent(in) :: ntint
    integer(kind=4) :: itint

    itint = 1
    DO WHILE (strtz(OPS_ACC1(0,0,0)) > tinthi(itint,ispec) .and. itint < ntint(ispec))
        itint = itint + 1
    END DO

    itndex(OPS_ACC2(0,0,0)) = itndex(OPS_ACC2(0,0,0)) + (itint-1)*ntbase**ipower

END SUBROUTINE boundt_kernel_eqE_zdir

SUBROUTINE boundt_kernel_eqF_xdir(erhs,yrhs,itndex,drhs,strtx,stryx,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs,yrhs
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: drhs,strtx,stryx

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC3(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strtx(OPS_ACC5(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtx(OPS_ACC5(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC4(0,0,0))*stryx(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strtx(OPS_ACC5(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE boundt_kernel_eqF_xdir

SUBROUTINE boundt_kernel_eqF_ydir(erhs,yrhs,itndex,drhs,strty,stryy,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs,yrhs
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: drhs,strty,stryy

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC3(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strty(OPS_ACC5(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strty(OPS_ACC5(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC4(0,0,0))*stryy(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strty(OPS_ACC5(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE boundt_kernel_eqF_ydir

SUBROUTINE boundt_kernel_eqF_zdir(erhs,yrhs,itndex,drhs,strtz,stryz,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs,yrhs
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: drhs,strtz,stryz

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC3(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strtz(OPS_ACC5(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtz(OPS_ACC5(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC4(0,0,0))*stryz(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strtz(OPS_ACC5(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE boundt_kernel_eqF_zdir

SUBROUTINE boundt_kernel_eqG_xdir(erhs,yrhs,itndex,strtx,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: yrhs,strtx

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC3(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strtx(OPS_ACC4(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtx(OPS_ACC4(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strtx(OPS_ACC4(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE boundt_kernel_eqG_xdir

SUBROUTINE boundt_kernel_eqG_ydir(erhs,yrhs,itndex,strty,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: yrhs,strty

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC3(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strty(OPS_ACC4(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strty(OPS_ACC4(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strty(OPS_ACC4(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE boundt_kernel_eqG_ydir

SUBROUTINE boundt_kernel_eqG_zdir(erhs,yrhs,itndex,strtz,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: yrhs,strtz

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC3(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strtz(OPS_ACC4(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtz(OPS_ACC4(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strtz(OPS_ACC4(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE boundt_kernel_eqG_zdir

SUBROUTINE boundt_kernel_eqH_xyz(urhs,vrhs,wrhs,erhs,drhs)

    real(kind=8), dimension(1) :: urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: drhs

    urhs(OPS_ACC1(0,0,0)) = 0.0_8
    vrhs(OPS_ACC2(0,0,0)) = 0.0_8
    wrhs(OPS_ACC3(0,0,0)) = 0.0_8
    erhs(OPS_ACC4(0,0,0)) = 0.0_8
    erhs(OPS_ACC4(0,0,0)) = drhs(OPS_ACC5(0,0,0))*erhs(OPS_ACC4(0,0,0))

END SUBROUTINE boundt_kernel_eqH_xyz 
