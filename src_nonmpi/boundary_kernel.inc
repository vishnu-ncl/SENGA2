SUBROUTINE boundary_kernel_CPandGAS_xdir(transp,store7,drhs,strgx,strrx)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=dp), dimension(1) :: strgx,strrx

    strgx(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strrx(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_CPandGAS_xdir

SUBROUTINE boundary_kernel_CPandGAS_ydir(transp,store7,drhs,strgy,strry)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=dp), dimension(1) :: strgy,strry

    strgy(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strry(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_CPandGAS_ydir

SUBROUTINE boundary_kernel_CPandGAS_zdir(transp,store7,drhs,strgz,strrz)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=dp), dimension(1) :: strgz,strrz

    strgz(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strrz(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_CPandGAS_zdir

SUBROUTINE boundary_kernel_internalenergy_xdir(erhs,strex)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: erhs
    real(kind=dp), dimension(1) :: strex
    
    strex(OPS_ACC2(0,0,0)) = erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_internalenergy_xdir

SUBROUTINE boundary_kernel_internalenergy_ydir(erhs,strey)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: erhs
    real(kind=dp), dimension(1) :: strey

    strey(OPS_ACC2(0,0,0)) = erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_internalenergy_ydir

SUBROUTINE boundary_kernel_internalenergy_zdir(erhs,strez)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: erhs
    real(kind=dp), dimension(1) :: strez

    strez(OPS_ACC2(0,0,0)) = erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_internalenergy_zdir

SUBROUTINE boundary_kernel_temperature_xdir(trun,store1,strtx,bcltx)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: trun,store1
    real(kind=dp), dimension(1) :: strtx,bcltx

    strtx(OPS_ACC3(0,0,0)) = trun(OPS_ACC1(0,0,0))
    bcltx(OPS_ACC4(0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_temperature_xdir

SUBROUTINE boundary_kernel_temperature_ydir(trun,store1,strty,bclty)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: trun,store1
    real(kind=dp), dimension(1) :: strty,bclty

    strty(OPS_ACC3(0,0,0)) = trun(OPS_ACC1(0,0,0))
    bclty(OPS_ACC4(0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_temperature_ydir

SUBROUTINE boundary_kernel_temperature_zdir(trun,store1,strtz,bcltz)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: trun,store1
    real(kind=dp), dimension(1) :: strtz,bcltz

    strtz(OPS_ACC3(0,0,0)) = trun(OPS_ACC1(0,0,0))
    bcltz(OPS_ACC4(0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_temperature_zdir

SUBROUTINE boundary_kernel_density_xdir(drhs,store,strdx,bcl2x)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: drhs,store
    real(kind=dp), dimension(1) :: strdx,bcl2x

    strdx(OPS_ACC3(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    bcl2x(OPS_ACC4(0,0,0)) = store(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_density_xdir

SUBROUTINE boundary_kernel_density_ydir(drhs,store,strdy,bcl2y)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: drhs,store
    real(kind=dp), dimension(1) :: strdy,bcl2y

    strdy(OPS_ACC3(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    bcl2y(OPS_ACC4(0,0,0)) = store(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_density_ydir

SUBROUTINE boundary_kernel_density_zdir(drhs,store,strdz,bcl2z)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: drhs,store
    real(kind=dp), dimension(1) :: strdz,bcl2z

    strdz(OPS_ACC3(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    bcl2z(OPS_ACC4(0,0,0)) = store(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_density_zdir

SUBROUTINE boundary_kernel_velcomp_xdir(utmp,vtmp,wtmp,strux,strvx,strwx)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: utmp,vtmp,wtmp
    real(kind=dp), dimension(1) :: strux,strvx,strwx
    
    strux(OPS_ACC4(0,0,0)) = utmp(OPS_ACC1(0,0,0))
    strvx(OPS_ACC5(0,0,0)) = vtmp(OPS_ACC2(0,0,0))
    strwx(OPS_ACC6(0,0,0)) = wtmp(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_velcomp_xdir

SUBROUTINE boundary_kernel_velcomp_ydir(utmp,vtmp,wtmp,struy,strvy,strwy)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: utmp,vtmp,wtmp
    real(kind=dp), dimension(1) :: struy,strvy,strwy

    struy(OPS_ACC4(0,0,0)) = utmp(OPS_ACC1(0,0,0))
    strvy(OPS_ACC5(0,0,0)) = vtmp(OPS_ACC2(0,0,0))
    strwy(OPS_ACC6(0,0,0)) = wtmp(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_velcomp_ydir

SUBROUTINE boundary_kernel_velcomp_zdir(utmp,vtmp,wtmp,struz,strvz,strwz)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: utmp,vtmp,wtmp
    real(kind=dp), dimension(1) :: struz,strvz,strwz

    struz(OPS_ACC4(0,0,0)) = utmp(OPS_ACC1(0,0,0))
    strvz(OPS_ACC5(0,0,0)) = vtmp(OPS_ACC2(0,0,0))
    strwz(OPS_ACC6(0,0,0)) = wtmp(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_velcomp_zdir

SUBROUTINE boundary_kernel_pressure_xdir(prun,store4,strpx,bcl5x)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: prun,store4
    real(kind=dp), dimension(1) :: strpx,bcl5x

    strpx(OPS_ACC3(0,0,0)) = prun(OPS_ACC1(0,0,0))
    bcl5x(OPS_ACC4(0,0,0)) = store4(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_pressure_xdir

SUBROUTINE boundary_kernel_pressure_ydir(prun,store5,strpy,bcl5y)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: prun,store5
    real(kind=dp), dimension(1) :: strpy,bcl5y

    strpy(OPS_ACC3(0,0,0)) = prun(OPS_ACC1(0,0,0))
    bcl5y(OPS_ACC4(0,0,0)) = store5(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_pressure_ydir

SUBROUTINE boundary_kernel_pressure_zdir(prun,store6,strpz,bcl5z)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: prun,store6
    real(kind=dp), dimension(1) :: strpz,bcl5z

    strpz(OPS_ACC3(0,0,0)) = prun(OPS_ACC1(0,0,0))
    bcl5z(OPS_ACC4(0,0,0)) = store6(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_pressure_zdir

SUBROUTINE boundary_kernel_velderiv_xdir(store,bclx)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: store
    real(kind=dp), dimension(1) :: bclx

    bclx(OPS_ACC2(0,0,0)) = store(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_velderiv_xdir

SUBROUTINE boundary_kernel_velderiv_ydir(store,bcly)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: store
    real(kind=dp), dimension(1) :: bcly

    bcly(OPS_ACC2(0,0,0)) = store(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_velderiv_ydir

SUBROUTINE boundary_kernel_velderiv_zdir(store,bclz)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: store
    real(kind=dp), dimension(1) :: bclz

    bclz(OPS_ACC2(0,0,0)) = store(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_velderiv_zdir

SUBROUTINE boundary_kernel_reaction_xdir(rate, ratex, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: ratex
    real(kind=dp), dimension(1), intent(in) :: rate
    integer, intent(in) :: ispec

    ratex(OPS_ACC_MD2(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE boundary_kernel_reaction_xdir

SUBROUTINE boundary_kernel_reaction_ydir(rate, ratey, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: ratey
    real(kind=dp), dimension(1), intent(in) :: rate
    integer, intent(in) :: ispec

    ratey(OPS_ACC_MD2(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE boundary_kernel_reaction_ydir

SUBROUTINE boundary_kernel_reaction_zdir(rate, ratez, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: ratez
    real(kind=dp), dimension(1), intent(in) :: rate
    integer, intent(in) :: ispec

    ratez(OPS_ACC_MD2(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE boundary_kernel_reaction_zdir

SUBROUTINE boundary_kernel_mass_xdir(yrhs, store1, stryx, bclyx, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: stryx, bclyx
    real(kind=dp), dimension(1), intent(in) :: yrhs, store1
    integer, intent(in) :: ispec

    stryx(OPS_ACC_MD3(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))
    bclyx(OPS_ACC_MD4(ispec,0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_mass_xdir

SUBROUTINE boundary_kernel_mass_ydir(yrhs, store2, stryy, bclyy, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: stryy, bclyy
    real(kind=dp), dimension(1), intent(in) :: yrhs, store2
    integer, intent(in) :: ispec

    stryy(OPS_ACC_MD3(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))
    bclyy(OPS_ACC_MD4(ispec,0,0,0)) = store2(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_mass_ydir

SUBROUTINE boundary_kernel_mass_zdir(yrhs, store3, stryz, bclyz, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: stryz, bclyz
    real(kind=dp), dimension(1), intent(in) :: yrhs, store3
    integer, intent(in) :: ispec

    stryz(OPS_ACC_MD3(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))
    bclyz(OPS_ACC_MD4(ispec,0,0,0)) = store3(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_mass_zdir

SUBROUTINE boundary_kernel_speciesH_xdir(utmp, strhx, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strhx
    real(kind=dp), dimension(1), intent(in) :: utmp
    integer, intent(in) :: ispec

    strhx(OPS_ACC_MD2(ispec,0,0,0)) = utmp(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_speciesH_xdir

SUBROUTINE boundary_kernel_speciesH_ydir(utmp, strhy, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strhy
    real(kind=dp), dimension(1), intent(in) :: utmp
    integer, intent(in) :: ispec

    strhy(OPS_ACC_MD2(ispec,0,0,0)) = utmp(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_speciesH_ydir

SUBROUTINE boundary_kernel_speciesH_zdir(utmp, strhz, ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: strhz
    real(kind=dp), dimension(1), intent(in) :: utmp
    integer, intent(in) :: ispec

    strhz(OPS_ACC_MD2(ispec,0,0,0)) = utmp(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_speciesH_zdir

