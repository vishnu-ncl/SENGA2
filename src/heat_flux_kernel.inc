SUBROUTINE heat_flux_kernel_thermal_fxlcnw(erhs,store1,store7,acbcxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: store1, store7
    real(kind=dp), dimension(1), intent(in) :: acbcxl
    real(kind=dp) :: fornow

    fornow = acbcxl(1)*store7(OPS_ACC3(1,0,0))*store1(OPS_ACC2(1,0,0)) &
           + acbcxl(2)*store7(OPS_ACC3(2,0,0))*store1(OPS_ACC2(2,0,0)) &
           + acbcxl(3)*store7(OPS_ACC3(3,0,0))*store1(OPS_ACC2(3,0,0)) &
           + acbcxl(4)*store7(OPS_ACC3(4,0,0))*store1(OPS_ACC2(4,0,0))
    
    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fxlcnw
    
SUBROUTINE heat_flux_kernel_thermal_fxrcnw(erhs,store1,store7,acbcxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: store1, store7
    real(kind=dp), dimension(1), intent(in) :: acbcxr
    real(kind=dp) :: fornow

    fornow = acbcxr(4)*store7(OPS_ACC3(-4,0,0))*store1(OPS_ACC2(-4,0,0)) &
           + acbcxr(3)*store7(OPS_ACC3(-3,0,0))*store1(OPS_ACC2(-3,0,0)) &
           + acbcxr(2)*store7(OPS_ACC3(-2,0,0))*store1(OPS_ACC2(-2,0,0)) &
           + acbcxr(1)*store7(OPS_ACC3(-1,0,0))*store1(OPS_ACC2(-1,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fxrcnw

SUBROUTINE heat_flux_kernel_thermal_fylcnw(erhs,store2,store7,acbcyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: store2, store7
    real(kind=dp), dimension(1), intent(in) :: acbcyl
    real(kind=dp) :: fornow

    fornow = acbcyl(1)*store7(OPS_ACC3(0,1,0))*store2(OPS_ACC2(0,1,0)) &
           + acbcyl(2)*store7(OPS_ACC3(0,2,0))*store2(OPS_ACC2(0,2,0)) &
           + acbcyl(3)*store7(OPS_ACC3(0,3,0))*store2(OPS_ACC2(0,3,0)) &
           + acbcyl(4)*store7(OPS_ACC3(0,4,0))*store2(OPS_ACC2(0,4,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fylcnw

SUBROUTINE heat_flux_kernel_thermal_fyrcnw(erhs,store2,store7,acbcyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: store2, store7
    real(kind=dp), dimension(1), intent(in) :: acbcyr
    real(kind=dp) :: fornow

    fornow = acbcyr(4)*store7(OPS_ACC3(0,-4,0))*store2(OPS_ACC2(0,-4,0)) &
           + acbcyr(3)*store7(OPS_ACC3(0,-3,0))*store2(OPS_ACC2(0,-3,0)) &
           + acbcyr(2)*store7(OPS_ACC3(0,-2,0))*store2(OPS_ACC2(0,-2,0)) &
           + acbcyr(1)*store7(OPS_ACC3(0,-1,0))*store2(OPS_ACC2(0,-1,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fyrcnw

SUBROUTINE heat_flux_kernel_thermal_fzlcnw(erhs,store3,store7,acbczl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: store3, store7
    real(kind=dp), dimension(1), intent(in) :: acbczl
    real(kind=dp) :: fornow

    fornow = acbczl(1)*store7(OPS_ACC3(0,0,1))*store3(OPS_ACC2(0,0,1)) &
           + acbczl(2)*store7(OPS_ACC3(0,0,2))*store3(OPS_ACC2(0,0,2)) &
           + acbczl(3)*store7(OPS_ACC3(0,0,3))*store3(OPS_ACC2(0,0,3)) &
           + acbczl(4)*store7(OPS_ACC3(0,0,4))*store3(OPS_ACC2(0,0,4))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fzlcnw

SUBROUTINE heat_flux_kernel_thermal_fzrcnw(erhs,store3,store7,acbczr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: store3, store7
    real(kind=dp), dimension(1), intent(in) :: acbczr
    real(kind=dp) :: fornow

    fornow = acbczr(4)*store7(OPS_ACC3(0,0,-4))*store3(OPS_ACC2(0,0,-4)) &
           + acbczr(3)*store7(OPS_ACC3(0,0,-3))*store3(OPS_ACC2(0,0,-3)) &
           + acbczr(2)*store7(OPS_ACC3(0,0,-2))*store3(OPS_ACC2(0,0,-2)) &
           + acbczr(1)*store7(OPS_ACC3(0,0,-1))*store3(OPS_ACC2(0,0,-1))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_thermal_fzrcnw

SUBROUTINE heat_flux_kernel_soret_fxlcnw(erhs,trun,tdrmix,store7,td1x,acbcxl,rgspec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    real(kind=dp), dimension(1), intent(in) :: acbcxl
    real(kind=dp), intent(in) :: rgspec
    real(kind=dp) :: fornow
    real(kind=dp) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(1,0,0))*tdrmix(OPS_ACC3(1,0,0))
    combo2 = trun(OPS_ACC2(2,0,0))*tdrmix(OPS_ACC3(2,0,0))
    combo3 = trun(OPS_ACC2(3,0,0))*tdrmix(OPS_ACC3(3,0,0))
    combo4 = trun(OPS_ACC2(4,0,0))*tdrmix(OPS_ACC3(4,0,0))

    fornow = acbcxl(1)*combo1*store7(OPS_ACC4(1,0,0))*td1x(OPS_ACC5(1,0,0)) &
           + acbcxl(2)*combo2*store7(OPS_ACC4(2,0,0))*td1x(OPS_ACC5(2,0,0)) &
           + acbcxl(3)*combo3*store7(OPS_ACC4(3,0,0))*td1x(OPS_ACC5(3,0,0)) &
           + acbcxl(4)*combo4*store7(OPS_ACC4(4,0,0))*td1x(OPS_ACC5(4,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE heat_flux_kernel_soret_fxlcnw
