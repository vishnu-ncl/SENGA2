SUBROUTINE boundary_kernel_CPandGAS_xdir(transp,store7,drhs,strgx,strrx)

    real(kind=8), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=8), dimension(1) :: strgx,strrx

    strgx(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strrx(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_CPandGAS_xdir

SUBROUTINE boundary_kernel_CPandGAS_ydir(transp,store7,drhs,strgy,strry)

    real(kind=8), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=8), dimension(1) :: strgy,strry

    strgy(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strry(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_CPandGAS_ydir

SUBROUTINE boundary_kernel_CPandGAS_zdir(transp,store7,drhs,strgz,strrz)

    real(kind=8), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=8), dimension(1) :: strgz,strrz

    strgz(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strrz(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_CPandGAS_zdir

SUBROUTINE boundary_kernel_internalenergy_xdir(erhs,strex)

    real(kind=8), dimension(1), intent(in) :: erhs
    real(kind=8), dimension(1) :: strex

    strex(OPS_ACC2(0,0,0)) = erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_internalenergy_xdir

SUBROUTINE boundary_kernel_internalenergy_ydir(erhs,strey)

    real(kind=8), dimension(1), intent(in) :: erhs
    real(kind=8), dimension(1) :: strey

    strey(OPS_ACC2(0,0,0)) = erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_internalenergy_ydir

SUBROUTINE boundary_kernel_internalenergy_zdir(erhs,strez)

    real(kind=8), dimension(1), intent(in) :: erhs
    real(kind=8), dimension(1) :: strez

    strez(OPS_ACC2(0,0,0)) = erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_internalenergy_zdir

SUBROUTINE boundary_kernel_temperature_xdir(trun,store1,strtx,bcltx)

    real(kind=8), dimension(1), intent(in) :: trun,store1
    real(kind=8), dimension(1) :: strtx,bcltx

    strtx(OPS_ACC3(0,0,0)) = trun(OPS_ACC1(0,0,0))
    bcltx(OPS_ACC4(0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_temperature_xdir

SUBROUTINE boundary_kernel_temperature_ydir(trun,store1,strty,bclty)

    real(kind=8), dimension(1), intent(in) :: trun,store1
    real(kind=8), dimension(1) :: strty,bclty

    strty(OPS_ACC3(0,0,0)) = trun(OPS_ACC1(0,0,0))
    bclty(OPS_ACC4(0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_temperature_ydir

SUBROUTINE boundary_kernel_temperature_zdir(trun,store1,strtz,bcltz)

    real(kind=8), dimension(1), intent(in) :: trun,store1
    real(kind=8), dimension(1) :: strtz,bcltz

    strtz(OPS_ACC3(0,0,0)) = trun(OPS_ACC1(0,0,0))
    bcltz(OPS_ACC4(0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_temperature_zdir

SUBROUTINE boundary_kernel_density_xdir(drhs,store,strdx,bcl2x)

    real(kind=8), dimension(1), intent(in) :: drhs,store
    real(kind=8), dimension(1) :: strdx,bcl2x

    strdx(OPS_ACC3(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    bcl2x(OPS_ACC4(0,0,0)) = store(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_density_xdir

SUBROUTINE boundary_kernel_density_ydir(drhs,store,strdy,bcl2y)

    real(kind=8), dimension(1), intent(in) :: drhs,store
    real(kind=8), dimension(1) :: strdy,bcl2y

    strdy(OPS_ACC3(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    bcl2y(OPS_ACC4(0,0,0)) = store(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_density_ydir

SUBROUTINE boundary_kernel_density_zdir(drhs,store,strdz,bcl2z)

    real(kind=8), dimension(1), intent(in) :: drhs,store
    real(kind=8), dimension(1) :: strdz,bcl2z

    strdz(OPS_ACC3(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    bcl2z(OPS_ACC4(0,0,0)) = store(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_density_zdir

SUBROUTINE boundary_kernel_velcomp_xdir(strux,strvx,strwx,t2bx,utmp,vtmp,wtmp,store1,store2)

    real(kind=8), dimension(1) :: strux,strvx,strwx,t2bx
    real(kind=8), dimension(1), intent(in) :: utmp,vtmp,wtmp,store1,store2

    strux(OPS_ACC1(0,0,0)) = utmp(OPS_ACC5(0,0,0))
    strvx(OPS_ACC2(0,0,0)) = vtmp(OPS_ACC6(0,0,0))
    strwx(OPS_ACC3(0,0,0)) = wtmp(OPS_ACC7(0,0,0))

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    t2bx(OPS_ACC4(0,0,0)) = -vtmp(OPS_ACC6(0,0,0))*store1(OPS_ACC8(0,0,0)) &
                            -wtmp(OPS_ACC7(0,0,0))*store2(OPS_ACC9(0,0,0))

END SUBROUTINE boundary_kernel_velcomp_xdir

SUBROUTINE boundary_kernel_velcomp_ydir(struy,strvy,strwy,t2by,utmp,vtmp,wtmp,store1,store2)

    real(kind=8), dimension(1) :: struy,strvy,strwy,t2by
    real(kind=8), dimension(1), intent(in) :: utmp,vtmp,wtmp,store1,store2

    struy(OPS_ACC1(0,0,0)) = utmp(OPS_ACC5(0,0,0))
    strvy(OPS_ACC2(0,0,0)) = vtmp(OPS_ACC6(0,0,0))
    strwy(OPS_ACC3(0,0,0)) = wtmp(OPS_ACC7(0,0,0))

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    t2by(OPS_ACC4(0,0,0)) = -utmp(OPS_ACC6(0,0,0))*store1(OPS_ACC8(0,0,0)) &
                            -wtmp(OPS_ACC7(0,0,0))*store2(OPS_ACC9(0,0,0))

END SUBROUTINE boundary_kernel_velcomp_ydir

SUBROUTINE boundary_kernel_velcomp_zdir(struz,strvz,strwz,t2bz,utmp,vtmp,wtmp,store1,store2)

    real(kind=8), dimension(1) :: struz,strvz,strwz,t2bz
    real(kind=8), dimension(1), intent(in) :: utmp,vtmp,wtmp,store1,store2

    struz(OPS_ACC1(0,0,0)) = utmp(OPS_ACC5(0,0,0))
    strvz(OPS_ACC2(0,0,0)) = vtmp(OPS_ACC6(0,0,0))
    strwz(OPS_ACC3(0,0,0)) = wtmp(OPS_ACC7(0,0,0))

!   EVALUATION OF TRANSVERSE CONDITIONS FOR B.C. (NC)
    t2bz(OPS_ACC4(0,0,0)) = -utmp(OPS_ACC6(0,0,0))*store1(OPS_ACC8(0,0,0)) &
                            -vtmp(OPS_ACC7(0,0,0))*store2(OPS_ACC9(0,0,0))

END SUBROUTINE boundary_kernel_velcomp_zdir

SUBROUTINE boundary_kernel_pressure_xdir(strpx,bcl5x,t3bx,t4bx,t51bx,t52bx,drhs,vtmp,wtmp,prun, &
                                         store2,store3,store4,store5,store6)

    real(kind=8), dimension(1) :: strpx,bcl5x,t3bx,t4bx,t51bx,t52bx
    real(kind=8), dimension(1), intent(in) :: drhs,vtmp,wtmp,prun,store2,store3,store4,store5,store6

    strpx(OPS_ACC1(0,0,0)) = prun(OPS_ACC10(0,0,0))
    bcl5x(OPS_ACC2(0,0,0)) = store4(OPS_ACC13(0,0,0))

!   BOUNDARY CONDITION DUE TO LODATO's FORMULATION
    t3bx(OPS_ACC3(0,0,0))  = t3bx(OPS_ACC3(0,0,0)) - store5(OPS_ACC14(0,0,0))/drhs(OPS_ACC7(0,0,0))
    t4bx(OPS_ACC4(0,0,0))  = t4bx(OPS_ACC4(0,0,0)) - store6(OPS_ACC15(0,0,0))/drhs(OPS_ACC7(0,0,0))
    t51bx(OPS_ACC5(0,0,0)) = -vtmp(OPS_ACC8(0,0,0))*store5(OPS_ACC14(0,0,0)) &
                             -wtmp(OPS_ACC9(0,0,0))*store6(OPS_ACC15(0,0,0))
    t52bx(OPS_ACC6(0,0,0)) = -prun(OPS_ACC10(0,0,0))*(store2(OPS_ACC11(0,0,0))+store3(OPS_ACC12(0,0,0)))

END SUBROUTINE boundary_kernel_pressure_xdir

SUBROUTINE boundary_kernel_pressure_ydir(strpy,bcl5y,t3by,t4by,t51by,t52by,drhs,utmp,wtmp,prun, &
                                        store1,store3,store4,store5,store6)

    real(kind=8), dimension(1) :: strpy,bcl5y,t3by,t4by,t51by,t52by
    real(kind=8), dimension(1), intent(in) :: drhs,utmp,wtmp,prun,store1,store3,store4,store5,store6

    strpy(OPS_ACC1(0,0,0)) = prun(OPS_ACC10(0,0,0))
    bcl5y(OPS_ACC2(0,0,0)) = store5(OPS_ACC14(0,0,0))

!   BOUNDARY CONDITION DUE TO LODATO's FORMULATION
!   NC and VM changed here: ycorrection
    t4by(OPS_ACC4(0,0,0))  = t4by(OPS_ACC4(0,0,0)) - store6(OPS_ACC15(0,0,0))/drhs(OPS_ACC7(0,0,0))
    t3by(OPS_ACC3(0,0,0))  = t3by(OPS_ACC3(0,0,0)) - store4(OPS_ACC13(0,0,0))/drhs(OPS_ACC7(0,0,0))
    t51by(OPS_ACC5(0,0,0)) = -utmp(OPS_ACC8(0,0,0))*store4(OPS_ACC13(0,0,0)) &
                             -wtmp(OPS_ACC9(0,0,0))*store6(OPS_ACC15(0,0,0))
    t52by(OPS_ACC6(0,0,0)) = -prun(OPS_ACC10(0,0,0))*(store1(OPS_ACC11(0,0,0))+store3(OPS_ACC12(0,0,0)))

END SUBROUTINE boundary_kernel_pressure_ydir

SUBROUTINE boundary_kernel_pressure_zdir(strpz,bcl5z,t3bz,t4bz,t51bz,t52bz,drhs,utmp,vtmp,prun, &
                                        store1,store2,store4,store5,store6)

    real(kind=8), dimension(1) :: strpz,bcl5z,t3bz,t4bz,t51bz,t52bz
    real(kind=8), dimension(1), intent(in) :: drhs,utmp,vtmp,prun,store1,store2,store4,store5,store6

    strpz(OPS_ACC1(0,0,0)) = prun(OPS_ACC10(0,0,0))
    bcl5z(OPS_ACC2(0,0,0)) = store6(OPS_ACC15(0,0,0))

!   BOUNDARY CONDITION DUE TO LODATO's FORMULATION
    t3bz(OPS_ACC3(0,0,0))  = t3bz(OPS_ACC3(0,0,0)) - store4(OPS_ACC13(0,0,0))/drhs(OPS_ACC7(0,0,0))
    t4bz(OPS_ACC4(0,0,0))  = t4bz(OPS_ACC4(0,0,0)) - store5(OPS_ACC14(0,0,0))/drhs(OPS_ACC7(0,0,0))
    t51bz(OPS_ACC5(0,0,0)) = -utmp(OPS_ACC8(0,0,0))*store4(OPS_ACC13(0,0,0)) &
                             -vtmp(OPS_ACC9(0,0,0))*store5(OPS_ACC14(0,0,0))
    t52bz(OPS_ACC6(0,0,0)) = -prun(OPS_ACC10(0,0,0))*(store1(OPS_ACC11(0,0,0))+store2(OPS_ACC12(0,0,0)))

END SUBROUTINE boundary_kernel_pressure_zdir

SUBROUTINE boundary_kernel_velderiv_xdir(store,bclx)

    real(kind=8), dimension(1), intent(in) :: store
    real(kind=8), dimension(1) :: bclx

    bclx(OPS_ACC2(0,0,0)) = store(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_velderiv_xdir

SUBROUTINE boundary_kernel_velderiv_ydir(store,bcly)

    real(kind=8), dimension(1), intent(in) :: store
    real(kind=8), dimension(1) :: bcly

    bcly(OPS_ACC2(0,0,0)) = store(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_velderiv_ydir

SUBROUTINE boundary_kernel_velderiv_zdir(store,bclz)

    real(kind=8), dimension(1), intent(in) :: store
    real(kind=8), dimension(1) :: bclz

    bclz(OPS_ACC2(0,0,0)) = store(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_velderiv_zdir

SUBROUTINE boundary_kernel_mass_xdir(yrhs, store1, store2, store3, drhs, vrhs, wrhs, stryx, bclyx, t6bx)

    real(kind=8), dimension(1) :: stryx, bclyx, t6bx
    real(kind=8), dimension(1), intent(in) :: yrhs, store1, store2, store3, drhs, vrhs, wrhs

    stryx(OPS_ACC8(0,0,0)) = yrhs(OPS_ACC1(0,0,0))
    bclyx(OPS_ACC9(0,0,0)) = store1(OPS_ACC2(0,0,0))

!   TRANSVERSE TERM FOR SPECIES TRANSPORT (see Eq. 3.74 of Lodato's thesis)
    t6bx(OPS_ACC10(0,0,0)) = - store2(OPS_ACC3(0,0,0))*vrhs(OPS_ACC6(0,0,0))/drhs(OPS_ACC5(0,0,0))  &
                             - store3(OPS_ACC4(0,0,0))*wrhs(OPS_ACC7(0,0,0))/drhs(OPS_ACC5(0,0,0))

END SUBROUTINE boundary_kernel_mass_xdir

SUBROUTINE boundary_kernel_mass_ydir(yrhs, store1, store2, store3, drhs, urhs, wrhs, stryy, bclyy, t6by)

    real(kind=8), dimension(1) :: stryy, bclyy, t6by
    real(kind=8), dimension(1), intent(in) :: yrhs, store1, store2, store3, drhs, urhs, wrhs

    stryy(OPS_ACC8(0,0,0)) = yrhs(OPS_ACC1(0,0,0))
    bclyy(OPS_ACC9(0,0,0)) = store2(OPS_ACC3(0,0,0))

!   TRANSVERSE TERM FOR SPECIES TRANSPORT (see Eq. 3.74 of Lodato's thesis)
    t6by(OPS_ACC10(0,0,0)) = - store1(OPS_ACC2(0,0,0))*urhs(OPS_ACC6(0,0,0))/drhs(OPS_ACC5(0,0,0))  &
                             - store3(OPS_ACC4(0,0,0))*wrhs(OPS_ACC7(0,0,0))/drhs(OPS_ACC5(0,0,0))

END SUBROUTINE boundary_kernel_mass_ydir

SUBROUTINE boundary_kernel_mass_zdir(yrhs, store1, store2, store3, drhs, urhs, vrhs, stryz, bclyz, t6bz)

    real(kind=8), dimension(1) :: stryz, bclyz, t6bz
    real(kind=8), dimension(1), intent(in) :: yrhs, store1, store2, store3, drhs, urhs, vrhs

    stryz(OPS_ACC8(0,0,0)) = yrhs(OPS_ACC1(0,0,0))
    bclyz(OPS_ACC9(0,0,0)) = store3(OPS_ACC4(0,0,0))

!   TRANSVERSE TERM FOR SPECIES TRANSPORT (see Eq. 3.74 of Lodato's thesis)
    t6bz(OPS_ACC10(0,0,0)) = - store1(OPS_ACC2(0,0,0))*urhs(OPS_ACC6(0,0,0))/drhs(OPS_ACC5(0,0,0))  &
                             - store2(OPS_ACC3(0,0,0))*vrhs(OPS_ACC7(0,0,0))/drhs(OPS_ACC5(0,0,0))

END SUBROUTINE boundary_kernel_mass_zdir

SUBROUTINE boundary_kernel_speciesH_xdir(utmp, strhx)

    real(kind=8), dimension(1) :: strhx
    real(kind=8), dimension(1), intent(in) :: utmp

    strhx(OPS_ACC2(0,0,0)) = utmp(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_speciesH_xdir

SUBROUTINE boundary_kernel_speciesH_ydir(utmp, strhy)

    real(kind=8), dimension(1) :: strhy
    real(kind=8), dimension(1), intent(in) :: utmp

    strhy(OPS_ACC2(0,0,0)) = utmp(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_speciesH_ydir

SUBROUTINE boundary_kernel_speciesH_zdir(utmp, strhz)

    real(kind=8), dimension(1) :: strhz
    real(kind=8), dimension(1), intent(in) :: utmp

    strhz(OPS_ACC2(0,0,0)) = utmp(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_speciesH_zdir

SUBROUTINE boundary_kernel_eqA_xdir(out_arr,vtmp,wtmp,store1,store2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: vtmp,wtmp,store1,store2

    out_arr(OPS_ACC1(0,0,0)) = -vtmp(OPS_ACC2(0,0,0))*store1(OPS_ACC4(0,0,0)) &
                               -wtmp(OPS_ACC3(0,0,0))*store2(OPS_ACC5(0,0,0))

END SUBROUTINE boundary_kernel_eqA_xdir

SUBROUTINE boundary_kernel_eqA_ydir(out_arr,utmp,wtmp,store1,store2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: utmp,wtmp,store1,store2

    out_arr(OPS_ACC1(0,0,0)) = -utmp(OPS_ACC2(0,0,0))*store1(OPS_ACC4(0,0,0)) &
                               -wtmp(OPS_ACC3(0,0,0))*store2(OPS_ACC5(0,0,0))

END SUBROUTINE boundary_kernel_eqA_ydir

SUBROUTINE boundary_kernel_eqA_zdir(out_arr,utmp,vtmp,store1,store2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: utmp,vtmp,store1,store2

    out_arr(OPS_ACC1(0,0,0)) = -utmp(OPS_ACC2(0,0,0))*store1(OPS_ACC4(0,0,0)) &
                               -vtmp(OPS_ACC3(0,0,0))*store2(OPS_ACC5(0,0,0))

END SUBROUTINE boundary_kernel_eqA_zdir
