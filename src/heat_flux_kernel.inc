SUBROUTINE heat_flux_kernel_thermal_fxlcnw(erhs,store1,store7,acbcxl)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store1, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8) :: fornow

    fornow = acbcxl(1)*store7(OPS_ACC3(1,0,0))*store1(OPS_ACC2(1,0,0)) &
           + acbcxl(2)*store7(OPS_ACC3(2,0,0))*store1(OPS_ACC2(2,0,0)) &
           + acbcxl(3)*store7(OPS_ACC3(3,0,0))*store1(OPS_ACC2(3,0,0)) &
           + acbcxl(4)*store7(OPS_ACC3(4,0,0))*store1(OPS_ACC2(4,0,0))
    
    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fxlcnw
    
SUBROUTINE heat_flux_kernel_thermal_fxrcnw(erhs,store1,store7,acbcxr)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store1, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8) :: fornow

    fornow = acbcxr(4)*store7(OPS_ACC3(-4,0,0))*store1(OPS_ACC2(-4,0,0)) &
           + acbcxr(3)*store7(OPS_ACC3(-3,0,0))*store1(OPS_ACC2(-3,0,0)) &
           + acbcxr(2)*store7(OPS_ACC3(-2,0,0))*store1(OPS_ACC2(-2,0,0)) &
           + acbcxr(1)*store7(OPS_ACC3(-1,0,0))*store1(OPS_ACC2(-1,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fxrcnw

SUBROUTINE heat_flux_kernel_thermal_fylcnw(erhs,store2,store7,acbcyl)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store2, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8) :: fornow

    fornow = acbcyl(1)*store7(OPS_ACC3(0,1,0))*store2(OPS_ACC2(0,1,0)) &
           + acbcyl(2)*store7(OPS_ACC3(0,2,0))*store2(OPS_ACC2(0,2,0)) &
           + acbcyl(3)*store7(OPS_ACC3(0,3,0))*store2(OPS_ACC2(0,3,0)) &
           + acbcyl(4)*store7(OPS_ACC3(0,4,0))*store2(OPS_ACC2(0,4,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fylcnw

SUBROUTINE heat_flux_kernel_thermal_fyrcnw(erhs,store2,store7,acbcyr)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store2, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8) :: fornow

    fornow = acbcyr(4)*store7(OPS_ACC3(0,-4,0))*store2(OPS_ACC2(0,-4,0)) &
           + acbcyr(3)*store7(OPS_ACC3(0,-3,0))*store2(OPS_ACC2(0,-3,0)) &
           + acbcyr(2)*store7(OPS_ACC3(0,-2,0))*store2(OPS_ACC2(0,-2,0)) &
           + acbcyr(1)*store7(OPS_ACC3(0,-1,0))*store2(OPS_ACC2(0,-1,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fyrcnw

SUBROUTINE heat_flux_kernel_thermal_fzlcnw(erhs,store3,store7,acbczl)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store3, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8) :: fornow

    fornow = acbczl(1)*store7(OPS_ACC3(0,0,1))*store3(OPS_ACC2(0,0,1)) &
           + acbczl(2)*store7(OPS_ACC3(0,0,2))*store3(OPS_ACC2(0,0,2)) &
           + acbczl(3)*store7(OPS_ACC3(0,0,3))*store3(OPS_ACC2(0,0,3)) &
           + acbczl(4)*store7(OPS_ACC3(0,0,4))*store3(OPS_ACC2(0,0,4))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fzlcnw

SUBROUTINE heat_flux_kernel_thermal_fzrcnw(erhs,store3,store7,acbczr)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store3, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8) :: fornow

    fornow = acbczr(4)*store7(OPS_ACC3(0,0,-4))*store3(OPS_ACC2(0,0,-4)) &
           + acbczr(3)*store7(OPS_ACC3(0,0,-3))*store3(OPS_ACC2(0,0,-3)) &
           + acbczr(2)*store7(OPS_ACC3(0,0,-2))*store3(OPS_ACC2(0,0,-2)) &
           + acbczr(1)*store7(OPS_ACC3(0,0,-1))*store3(OPS_ACC2(0,0,-1))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fzrcnw

SUBROUTINE heat_flux_kernel_soret_fxlcnw(erhs,trun,tdrmix,store7,td1x,acbcxl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(1,0,0))*tdrmix(OPS_ACC3(1,0,0))
    combo2 = trun(OPS_ACC2(2,0,0))*tdrmix(OPS_ACC3(2,0,0))
    combo3 = trun(OPS_ACC2(3,0,0))*tdrmix(OPS_ACC3(3,0,0))
    combo4 = trun(OPS_ACC2(4,0,0))*tdrmix(OPS_ACC3(4,0,0))

    fornow = acbcxl(1)*combo1*store7(OPS_ACC4(1,0,0))*td1x(OPS_ACC5(1,0,0)) &
           + acbcxl(2)*combo2*store7(OPS_ACC4(2,0,0))*td1x(OPS_ACC5(2,0,0)) &
           + acbcxl(3)*combo3*store7(OPS_ACC4(3,0,0))*td1x(OPS_ACC5(3,0,0)) &
           + acbcxl(4)*combo4*store7(OPS_ACC4(4,0,0))*td1x(OPS_ACC5(4,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_soret_fxlcnw

SUBROUTINE heat_flux_kernel_soret_fxrcnw(erhs,trun,tdrmix,store7,td1x,acbcxr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(-4,0,0))*tdrmix(OPS_ACC3(-4,0,0))
    combo2 = trun(OPS_ACC2(-3,0,0))*tdrmix(OPS_ACC3(-3,0,0))
    combo3 = trun(OPS_ACC2(-2,0,0))*tdrmix(OPS_ACC3(-2,0,0))
    combo4 = trun(OPS_ACC2(-1,0,0))*tdrmix(OPS_ACC3(-1,0,0))

    fornow = acbcxr(4)*combo1*store7(OPS_ACC4(-4,0,0))*td1x(OPS_ACC5(-4,0,0)) &
           + acbcxr(3)*combo2*store7(OPS_ACC4(-3,0,0))*td1x(OPS_ACC5(-3,0,0)) &
           + acbcxr(2)*combo3*store7(OPS_ACC4(-2,0,0))*td1x(OPS_ACC5(-2,0,0)) &
           + acbcxr(1)*combo4*store7(OPS_ACC4(-1,0,0))*td1x(OPS_ACC5(-1,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_soret_fxrcnw

SUBROUTINE heat_flux_kernel_soret_fylcnw(erhs,trun,tdrmix,store7,td1y,acbcyl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(0,1,0))*tdrmix(OPS_ACC3(0,1,0))
    combo2 = trun(OPS_ACC2(0,2,0))*tdrmix(OPS_ACC3(0,2,0))
    combo3 = trun(OPS_ACC2(0,3,0))*tdrmix(OPS_ACC3(0,3,0))
    combo4 = trun(OPS_ACC2(0,4,0))*tdrmix(OPS_ACC3(0,4,0))

    fornow = acbcyl(1)*combo1*store7(OPS_ACC4(0,1,0))*td1y(OPS_ACC5(0,1,0)) &
           + acbcyl(2)*combo2*store7(OPS_ACC4(0,2,0))*td1y(OPS_ACC5(0,2,0)) &
           + acbcyl(3)*combo3*store7(OPS_ACC4(0,3,0))*td1y(OPS_ACC5(0,3,0)) &
           + acbcyl(4)*combo4*store7(OPS_ACC4(0,4,0))*td1y(OPS_ACC5(0,4,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_soret_fylcnw

SUBROUTINE heat_flux_kernel_soret_fyrcnw(erhs,trun,tdrmix,store7,td1y,acbcyr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(0,-4,0))*tdrmix(OPS_ACC3(0,-4,0))
    combo2 = trun(OPS_ACC2(0,-3,0))*tdrmix(OPS_ACC3(0,-3,0))
    combo3 = trun(OPS_ACC2(0,-2,0))*tdrmix(OPS_ACC3(0,-2,0))
    combo4 = trun(OPS_ACC2(0,-1,0))*tdrmix(OPS_ACC3(0,-1,0))

    fornow = acbcyr(4)*combo1*store7(OPS_ACC4(0,-4,0))*td1y(OPS_ACC5(0,-4,0)) &
           + acbcyr(3)*combo2*store7(OPS_ACC4(0,-3,0))*td1y(OPS_ACC5(0,-3,0)) &
           + acbcyr(2)*combo3*store7(OPS_ACC4(0,-2,0))*td1y(OPS_ACC5(0,-2,0)) &
           + acbcyr(1)*combo4*store7(OPS_ACC4(0,-1,0))*td1y(OPS_ACC5(0,-1,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_soret_fyrcnw

SUBROUTINE heat_flux_kernel_soret_fzlcnw(erhs,trun,tdrmix,store7,td1z,acbczl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(0,0,1))*tdrmix(OPS_ACC3(0,0,1))
    combo2 = trun(OPS_ACC2(0,0,2))*tdrmix(OPS_ACC3(0,0,2))
    combo3 = trun(OPS_ACC2(0,0,3))*tdrmix(OPS_ACC3(0,0,3))
    combo4 = trun(OPS_ACC2(0,0,4))*tdrmix(OPS_ACC3(0,0,4))

    fornow = acbczl(1)*combo1*store7(OPS_ACC4(0,0,1))*td1z(OPS_ACC5(0,0,1)) &
           + acbczl(2)*combo2*store7(OPS_ACC4(0,0,2))*td1z(OPS_ACC5(0,0,2)) &
           + acbczl(3)*combo3*store7(OPS_ACC4(0,0,3))*td1z(OPS_ACC5(0,0,3)) &
           + acbczl(4)*combo4*store7(OPS_ACC4(0,0,4))*td1z(OPS_ACC5(0,0,4))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_soret_fzlcnw

SUBROUTINE heat_flux_kernel_soret_fzrcnw(erhs,trun,tdrmix,store7,td1z,acbczr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(0,0,-4))*tdrmix(OPS_ACC3(0,0,-4))
    combo2 = trun(OPS_ACC2(0,0,-3))*tdrmix(OPS_ACC3(0,0,-3))
    combo3 = trun(OPS_ACC2(0,0,-2))*tdrmix(OPS_ACC3(0,0,-2))
    combo4 = trun(OPS_ACC2(0,0,-1))*tdrmix(OPS_ACC3(0,0,-1))

    fornow = acbczr(4)*combo1*store7(OPS_ACC4(0,0,-4))*td1z(OPS_ACC5(0,0,-4)) &
           + acbczr(3)*combo2*store7(OPS_ACC4(0,0,-3))*td1z(OPS_ACC5(0,0,-3)) &
           + acbczr(2)*combo3*store7(OPS_ACC4(0,0,-2))*td1z(OPS_ACC5(0,0,-2)) &
           + acbczr(1)*combo4*store7(OPS_ACC4(0,0,-1))*td1z(OPS_ACC5(0,0,-1))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_soret_fzrcnw

SUBROUTINE heat_flux_kernel_enthalpy_fxldfw(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1x,acbcxl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1x
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbcxl(1)*store7(OPS_ACC7(1,0,0))*td1x(OPS_ACC8(1,0,0)) &
           + acbcxl(2)*store7(OPS_ACC7(2,0,0))*td1x(OPS_ACC8(2,0,0)) &
           + acbcxl(3)*store7(OPS_ACC7(3,0,0))*td1x(OPS_ACC8(3,0,0)) &
           + acbcxl(4)*store7(OPS_ACC7(4,0,0))*td1x(OPS_ACC8(4,0,0))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE heat_flux_kernel_enthalpy_fxldfw

SUBROUTINE heat_flux_kernel_enthalpy_fxrdfw(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1x,acbcxr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1x
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbcxr(4)*store7(OPS_ACC7(-4,0,0))*td1x(OPS_ACC8(-4,0,0)) &
           + acbcxr(3)*store7(OPS_ACC7(-3,0,0))*td1x(OPS_ACC8(-3,0,0)) &
           + acbcxr(2)*store7(OPS_ACC7(-2,0,0))*td1x(OPS_ACC8(-2,0,0)) &
           + acbcxr(1)*store7(OPS_ACC7(-1,0,0))*td1x(OPS_ACC8(-1,0,0))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE heat_flux_kernel_enthalpy_fxrdfw

SUBROUTINE heat_flux_kernel_enthalpy_fyldfw(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1y,acbcyl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1y
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbcyl(1)*store7(OPS_ACC7(0,1,0))*td1y(OPS_ACC8(0,1,0)) &
           + acbcyl(2)*store7(OPS_ACC7(0,2,0))*td1y(OPS_ACC8(0,2,0)) &
           + acbcyl(3)*store7(OPS_ACC7(0,3,0))*td1y(OPS_ACC8(0,3,0)) &
           + acbcyl(4)*store7(OPS_ACC7(0,4,0))*td1y(OPS_ACC8(0,4,0))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE heat_flux_kernel_enthalpy_fyldfw

SUBROUTINE heat_flux_kernel_enthalpy_fyrdfw(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1y,acbcyr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1y
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbcyr(4)*store7(OPS_ACC7(0,-4,0))*td1y(OPS_ACC8(0,-4,0)) &
           + acbcyr(3)*store7(OPS_ACC7(0,-3,0))*td1y(OPS_ACC8(0,-3,0)) &
           + acbcyr(2)*store7(OPS_ACC7(0,-2,0))*td1y(OPS_ACC8(0,-2,0)) &
           + acbcyr(1)*store7(OPS_ACC7(0,-1,0))*td1y(OPS_ACC8(0,-1,0))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE heat_flux_kernel_enthalpy_fyrdfw

SUBROUTINE heat_flux_kernel_enthalpy_fzldfw(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1z,acbczl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1z
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbczl(1)*store7(OPS_ACC7(0,0,1))*td1z(OPS_ACC8(0,0,1)) &
           + acbczl(2)*store7(OPS_ACC7(0,0,2))*td1z(OPS_ACC8(0,0,2)) &
           + acbczl(3)*store7(OPS_ACC7(0,0,3))*td1z(OPS_ACC8(0,0,3)) &
           + acbczl(4)*store7(OPS_ACC7(0,0,4))*td1z(OPS_ACC8(0,0,4))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE heat_flux_kernel_enthalpy_fzldfw

SUBROUTINE heat_flux_kernel_enthalpy_fzrdfw(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1z,acbczr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1z
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbczr(4)*store7(OPS_ACC7(0,0,-4))*td1z(OPS_ACC8(0,0,-4)) &
           + acbczr(3)*store7(OPS_ACC7(0,0,-3))*td1z(OPS_ACC8(0,0,-3)) &
           + acbczr(2)*store7(OPS_ACC7(0,0,-2))*td1z(OPS_ACC8(0,0,-2)) &
           + acbczr(1)*store7(OPS_ACC7(0,0,-1))*td1z(OPS_ACC8(0,0,-1))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE heat_flux_kernel_enthalpy_fzrdfw

SUBROUTINE heat_flux_kernel_enthalpy2_fxldfw(rate,erhs,vtmp,utmp,store7,arr1,acbcxl,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8) :: fornow

    fornow = acbcxl(1)*store7(OPS_ACC5(1,0,0))*arr1(OPS_ACC6(1,0,0)) &
           + acbcxl(2)*store7(OPS_ACC5(2,0,0))*arr1(OPS_ACC6(2,0,0)) &
           + acbcxl(3)*store7(OPS_ACC5(3,0,0))*arr1(OPS_ACC6(3,0,0)) &
           + acbcxl(4)*store7(OPS_ACC5(4,0,0))*arr1(OPS_ACC6(4,0,0))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE heat_flux_kernel_enthalpy2_fxldfw

SUBROUTINE heat_flux_kernel_enthalpy2_fxrdfw(rate,erhs,vtmp,utmp,store7,arr1,acbcxr,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8) :: fornow

    fornow = acbcxr(4)*store7(OPS_ACC5(-4,0,0))*arr1(OPS_ACC6(-4,0,0)) &
           + acbcxr(3)*store7(OPS_ACC5(-3,0,0))*arr1(OPS_ACC6(-3,0,0)) &
           + acbcxr(2)*store7(OPS_ACC5(-2,0,0))*arr1(OPS_ACC6(-2,0,0)) &
           + acbcxr(1)*store7(OPS_ACC5(-1,0,0))*arr1(OPS_ACC6(-1,0,0))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE heat_flux_kernel_enthalpy2_fxrdfw

SUBROUTINE heat_flux_kernel_enthalpy2_fyldfw(rate,erhs,vtmp,utmp,store7,arr1,acbcyl,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8) :: fornow

    fornow = acbcyl(1)*store7(OPS_ACC5(0,1,0))*arr1(OPS_ACC6(0,1,0)) &
           + acbcyl(2)*store7(OPS_ACC5(0,2,0))*arr1(OPS_ACC6(0,2,0)) &
           + acbcyl(3)*store7(OPS_ACC5(0,3,0))*arr1(OPS_ACC6(0,3,0)) &
           + acbcyl(4)*store7(OPS_ACC5(0,4,0))*arr1(OPS_ACC6(0,4,0))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE heat_flux_kernel_enthalpy2_fyldfw

SUBROUTINE heat_flux_kernel_enthalpy2_fyrdfw(rate,erhs,vtmp,utmp,store7,arr1,acbcyr,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8) :: fornow

    fornow = acbcyr(4)*store7(OPS_ACC5(0,-4,0))*arr1(OPS_ACC6(0,-4,0)) &
           + acbcyr(3)*store7(OPS_ACC5(0,-3,0))*arr1(OPS_ACC6(0,-3,0)) &
           + acbcyr(2)*store7(OPS_ACC5(0,-2,0))*arr1(OPS_ACC6(0,-2,0)) &
           + acbcyr(1)*store7(OPS_ACC5(0,-1,0))*arr1(OPS_ACC6(0,-1,0))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE heat_flux_kernel_enthalpy2_fyrdfw

SUBROUTINE heat_flux_kernel_enthalpy2_fzldfw(rate,erhs,vtmp,utmp,store7,arr1,acbczl,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8) :: fornow

    fornow = acbczl(1)*store7(OPS_ACC5(0,0,1))*arr1(OPS_ACC6(0,0,1)) &
           + acbczl(2)*store7(OPS_ACC5(0,0,2))*arr1(OPS_ACC6(0,0,2)) &
           + acbczl(3)*store7(OPS_ACC5(0,0,3))*arr1(OPS_ACC6(0,0,3)) &
           + acbczl(4)*store7(OPS_ACC5(0,0,4))*arr1(OPS_ACC6(0,0,4))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE heat_flux_kernel_enthalpy2_fzldfw

SUBROUTINE heat_flux_kernel_enthalpy2_fzrdfw(rate,erhs,vtmp,utmp,store7,arr1,acbczr,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8) :: fornow

    fornow = acbczr(4)*store7(OPS_ACC5(0,0,-4))*arr1(OPS_ACC6(0,0,-4)) &
           + acbczr(3)*store7(OPS_ACC5(0,0,-3))*arr1(OPS_ACC6(0,0,-3)) &
           + acbczr(2)*store7(OPS_ACC5(0,0,-2))*arr1(OPS_ACC6(0,0,-2)) &
           + acbczr(1)*store7(OPS_ACC5(0,0,-1))*arr1(OPS_ACC6(0,0,-1))

    rate(OPS_ACC_MD1(ispec,0,0,0)) = rate(OPS_ACC_MD1(ispec,0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE heat_flux_kernel_enthalpy2_fzrdfw

SUBROUTINE heat_flux_kernel_isothermal_fxladw(erhs,trun,tdrmix,store7,td1x,acbcxl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    integer , intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1x(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(1,0,0))*tdrmix(OPS_ACC3(1,0,0))*store7(OPS_ACC4(1,0,0))*td1x(OPS_ACC5(1,0,0))
    combo2 = trun(OPS_ACC2(2,0,0))*tdrmix(OPS_ACC3(2,0,0))*store7(OPS_ACC4(2,0,0))*td1x(OPS_ACC5(2,0,0))
    combo3 = trun(OPS_ACC2(3,0,0))*tdrmix(OPS_ACC3(3,0,0))*store7(OPS_ACC4(3,0,0))*td1x(OPS_ACC5(3,0,0))
    combo4 = trun(OPS_ACC2(4,0,0))*tdrmix(OPS_ACC3(4,0,0))*store7(OPS_ACC4(4,0,0))*td1x(OPS_ACC5(4,0,0))

    fornow = acbcxl(1)*rgspec(ispec)*(combo1-combo0) &
           + acbcxl(2)*rgspec(ispec)*(combo2-combo0) &
           + acbcxl(3)*rgspec(ispec)*(combo3-combo0) &
           + acbcxl(4)*rgspec(ispec)*(combo4-combo0)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow    

END SUBROUTINE heat_flux_kernel_isothermal_fxladw

SUBROUTINE heat_flux_kernel_isothermal_fxradw(erhs,trun,tdrmix,store7,td1x,acbcxr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1x(OPS_ACC5(0,0,0))

    combo4 = trun(OPS_ACC2(-4,0,0))*tdrmix(OPS_ACC3(-4,0,0))*store7(OPS_ACC4(-4,0,0))*td1x(OPS_ACC5(-4,0,0))
    combo3 = trun(OPS_ACC2(-3,0,0))*tdrmix(OPS_ACC3(-3,0,0))*store7(OPS_ACC4(-3,0,0))*td1x(OPS_ACC5(-3,0,0))
    combo2 = trun(OPS_ACC2(-2,0,0))*tdrmix(OPS_ACC3(-2,0,0))*store7(OPS_ACC4(-2,0,0))*td1x(OPS_ACC5(-2,0,0))
    combo1 = trun(OPS_ACC2(-1,0,0))*tdrmix(OPS_ACC3(-1,0,0))*store7(OPS_ACC4(-1,0,0))*td1x(OPS_ACC5(-1,0,0))

    fornow = acbcxr(4)*rgspec(ispec)*(combo0-combo4) &
           + acbcxr(3)*rgspec(ispec)*(combo0-combo3) &
           + acbcxr(2)*rgspec(ispec)*(combo0-combo2) &
           + acbcxr(1)*rgspec(ispec)*(combo0-combo1)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_isothermal_fxradw

SUBROUTINE heat_flux_kernel_isothermal_fyladw(erhs,trun,tdrmix,store7,td1y,acbcyl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec    
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1y(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(0,1,0))*tdrmix(OPS_ACC3(0,1,0))*store7(OPS_ACC4(0,1,0))*td1y(OPS_ACC5(0,1,0))
    combo2 = trun(OPS_ACC2(0,2,0))*tdrmix(OPS_ACC3(0,2,0))*store7(OPS_ACC4(0,2,0))*td1y(OPS_ACC5(0,2,0))
    combo3 = trun(OPS_ACC2(0,3,0))*tdrmix(OPS_ACC3(0,3,0))*store7(OPS_ACC4(0,3,0))*td1y(OPS_ACC5(0,3,0))
    combo4 = trun(OPS_ACC2(0,4,0))*tdrmix(OPS_ACC3(0,4,0))*store7(OPS_ACC4(0,4,0))*td1y(OPS_ACC5(0,4,0))

    fornow = acbcyl(1)*rgspec(ispec)*(combo1-combo0) &
           + acbcyl(2)*rgspec(ispec)*(combo2-combo0) &
           + acbcyl(3)*rgspec(ispec)*(combo3-combo0) &
           + acbcyl(4)*rgspec(ispec)*(combo4-combo0)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_isothermal_fyladw

SUBROUTINE heat_flux_kernel_isothermal_fyradw(erhs,trun,tdrmix,store7,td1y,acbcyr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1y(OPS_ACC5(0,0,0))

    combo4 = trun(OPS_ACC2(0,-4,0))*tdrmix(OPS_ACC3(0,-4,0))*store7(OPS_ACC4(0,-4,0))*td1y(OPS_ACC5(0,-4,0))
    combo3 = trun(OPS_ACC2(0,-3,0))*tdrmix(OPS_ACC3(0,-3,0))*store7(OPS_ACC4(0,-3,0))*td1y(OPS_ACC5(0,-3,0))
    combo2 = trun(OPS_ACC2(0,-2,0))*tdrmix(OPS_ACC3(0,-2,0))*store7(OPS_ACC4(0,-2,0))*td1y(OPS_ACC5(0,-2,0))
    combo1 = trun(OPS_ACC2(0,-1,0))*tdrmix(OPS_ACC3(0,-1,0))*store7(OPS_ACC4(0,-1,0))*td1y(OPS_ACC5(0,-1,0))

    fornow = acbcyr(4)*rgspec(ispec)*(combo0-combo4) &
           + acbcyr(3)*rgspec(ispec)*(combo0-combo3) &
           + acbcyr(2)*rgspec(ispec)*(combo0-combo2) &
           + acbcyr(1)*rgspec(ispec)*(combo0-combo1)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_isothermal_fyradw

SUBROUTINE heat_flux_kernel_isothermal_fzladw(erhs,trun,tdrmix,store7,td1z,acbczl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1z(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(0,0,1))*tdrmix(OPS_ACC3(0,0,1))*store7(OPS_ACC4(0,0,1))*td1z(OPS_ACC5(0,0,1))
    combo2 = trun(OPS_ACC2(0,0,2))*tdrmix(OPS_ACC3(0,0,2))*store7(OPS_ACC4(0,0,2))*td1z(OPS_ACC5(0,0,2))
    combo3 = trun(OPS_ACC2(0,0,3))*tdrmix(OPS_ACC3(0,0,3))*store7(OPS_ACC4(0,0,3))*td1z(OPS_ACC5(0,0,3))
    combo4 = trun(OPS_ACC2(0,0,4))*tdrmix(OPS_ACC3(0,0,4))*store7(OPS_ACC4(0,0,4))*td1z(OPS_ACC5(0,0,4))

    fornow = acbczl(1)*rgspec(ispec)*(combo1-combo0) &
           + acbczl(2)*rgspec(ispec)*(combo2-combo0) &
           + acbczl(3)*rgspec(ispec)*(combo3-combo0) &
           + acbczl(4)*rgspec(ispec)*(combo4-combo0)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_isothermal_fzladw

SUBROUTINE heat_flux_kernel_isothermal_fzradw(erhs,trun,tdrmix,store7,td1z,acbczr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    integer, intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1z(OPS_ACC5(0,0,0))

    combo4 = trun(OPS_ACC2(0,0,-4))*tdrmix(OPS_ACC3(0,0,-4))*store7(OPS_ACC4(0,0,-4))*td1z(OPS_ACC5(0,0,-4))
    combo3 = trun(OPS_ACC2(0,0,-3))*tdrmix(OPS_ACC3(0,0,-3))*store7(OPS_ACC4(0,0,-3))*td1z(OPS_ACC5(0,0,-3))
    combo2 = trun(OPS_ACC2(0,0,-2))*tdrmix(OPS_ACC3(0,0,-2))*store7(OPS_ACC4(0,0,-2))*td1z(OPS_ACC5(0,0,-2))
    combo1 = trun(OPS_ACC2(0,0,-1))*tdrmix(OPS_ACC3(0,0,-1))*store7(OPS_ACC4(0,0,-1))*td1z(OPS_ACC5(0,0,-1))

    fornow = acbczr(4)*rgspec(ispec)*(combo0-combo4) &
           + acbczr(3)*rgspec(ispec)*(combo0-combo3) &
           + acbczr(2)*rgspec(ispec)*(combo0-combo2) &
           + acbczr(1)*rgspec(ispec)*(combo0-combo1)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE heat_flux_kernel_isothermal_fzradw
