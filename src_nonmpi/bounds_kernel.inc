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

SUBROUTINE bounds_kernel_sound_speed_xdir(acoux,ova2x,strgx,gam1x,strtx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: acoux,ova2x
    real(kind=dp), dimension(1), intent(in) :: strgx,gam1x,strtx
    real(kind=dp) :: fornow
    
    fornow = strgx(OPS_ACC3(0,0,0))*gam1x(OPS_ACC4(0,0,0))*strtx(OPS_ACC5(0,0,0))
    acoux(OPS_ACC1(0,0,0)) = sqrt(fornow)
    ova2x(OPS_ACC2(0,0,0)) = 1.0_dp/fornow

END SUBROUTINE bounds_kernel_sound_speed_xdir

SUBROUTINE bounds_kernel_sound_speed_ydir(acouy,ova2y,strgy,gam1y,strty)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: acouy,ova2y
    real(kind=dp), dimension(1), intent(in) :: strgy,gam1y,strty
    real(kind=dp) :: fornow

    fornow = strgy(OPS_ACC3(0,0,0))*gam1y(OPS_ACC4(0,0,0))*strty(OPS_ACC5(0,0,0))
    acouy(OPS_ACC1(0,0,0)) = sqrt(fornow)
    ova2y(OPS_ACC2(0,0,0)) = 1.0_dp/fornow

END SUBROUTINE bounds_kernel_sound_speed_ydir

SUBROUTINE bounds_kernel_sound_speed_zdir(acouz,ova2z,strgz,gam1z,strtz)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: acouz,ova2z
    real(kind=dp), dimension(1), intent(in) :: strgz,gam1z,strtz
    real(kind=dp) :: fornow

    fornow = strgz(OPS_ACC3(0,0,0))*gam1z(OPS_ACC4(0,0,0))*strtz(OPS_ACC5(0,0,0))
    acouz(OPS_ACC1(0,0,0)) = sqrt(fornow)
    ova2z(OPS_ACC2(0,0,0)) = 1.0_dp/fornow

END SUBROUTINE bounds_kernel_sound_speed_zdir

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

SUBROUTINE bounds_kernel_outflowBC1_computeL_xl(bcl5xl,struxl,acouxl,strdxl,bcl1xl,sorpxl,strpxl,cobcxl,pinfxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl5xl
    real(kind=dp), dimension(1), intent(in) :: struxl,acouxl,strdxl,bcl1xl,sorpxl,strpxl
    real(kind=dp), intent(in) :: cobcxl,pinfxl

!   OLD VALUE OF L5X
    bcl5xl(OPS_ACC1(0,0,0)) = 0.5_dp*(struxl(OPS_ACC2(0,0,0))+acouxl(OPS_ACC3(0,0,0)))  &
            *(bcl5xl(OPS_ACC1(0,0,0))+strdxl(OPS_ACC4(0,0,0))*acouxl(OPS_ACC3(0,0,0))*bcl1xl(OPS_ACC5(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L5X
    bcl5xl(OPS_ACC1(0,0,0))= 0.5_dp*sorpxl(OPS_ACC6(0,0,0))  &
            + cobcxl*acouxl(OPS_ACC3(0,0,0))*(strpxl(OPS_ACC7(0,0,0))-pinfxl) - bcl5xl(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_computeL_xl

SUBROUTINE bounds_kernel_outflowBC1_computeL_yl(bcl5yl,strvyl,acouyl,strdyl,bcl1yl,sorpyl,strpyl,cobcyl,pinfyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl5yl
    real(kind=dp), dimension(1), intent(in) :: strvyl,acouyl,strdyl,bcl1yl,sorpyl,strpyl
    real(kind=dp), intent(in) :: cobcyl,pinfyl

!   OLD VALUE OF L5Y
    bcl5yl(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC2(0,0,0))+acouyl(OPS_ACC3(0,0,0)))  &
            *(bcl5yl(OPS_ACC1(0,0,0))+strdyl(OPS_ACC4(0,0,0))*acouyl(OPS_ACC3(0,0,0))*bcl1yl(OPS_ACC5(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L5Y
    bcl5yl(OPS_ACC1(0,0,0))= 0.5_dp*sorpyl(OPS_ACC6(0,0,0))  &
            + cobcyl*acouyl(OPS_ACC3(0,0,0))*(strpyl(OPS_ACC7(0,0,0))-pinfyl) - bcl5yl(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_computeL_yl

SUBROUTINE bounds_kernel_outflowBC1_computeL_zl(bcl5zl,strwzl,acouzl,strdzl,bcl1zl,sorpzl,strpzl,cobczl,pinfzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl5zl
    real(kind=dp), dimension(1), intent(in) :: strwzl,acouzl,strdzl,bcl1zl,sorpzl,strpzl
    real(kind=dp), intent(in) :: cobczl,pinfzl

!   OLD VALUE OF L5Z
    bcl5zl(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC2(0,0,0))+acouzl(OPS_ACC3(0,0,0)))  &
            *(bcl5zl(OPS_ACC1(0,0,0))+strdzl(OPS_ACC4(0,0,0))*acouzl(OPS_ACC3(0,0,0))*bcl1zl(OPS_ACC5(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L5Z
    bcl5zl(OPS_ACC1(0,0,0))= 0.5_dp*sorpzl(OPS_ACC6(0,0,0))  &
            + cobczl*acouzl(OPS_ACC3(0,0,0))*(strpzl(OPS_ACC7(0,0,0))-pinfzl) - bcl5zl(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_computeL_zl

SUBROUTINE bounds_kernel_inflowBC1_computeL_xl(bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,bcl1xl,struxl,ova2xl,sorpxl,strpxl,cobcxl,pinfxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl2xl,bcl3xl,bcl4xl,bcl5xl
    real(kind=dp), dimension(1), intent(in) :: strdxl,acouxl,bcl1xl,struxl,ova2xl,sorpxl,strpxl
    real(kind=dp), intent(in) :: cobcxl,pinfxl
    real(kind=dp) :: fornow
    
!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC5(0,0,0))*acouxl(OPS_ACC6(0,0,0))*bcl1xl(OPS_ACC7(0,0,0))
    bcl2xl(OPS_ACC1(0,0,0)) = struxl(OPS_ACC8(0,0,0))  &
            *(bcl2xl(OPS_ACC1(0,0,0))-bcl5xl(OPS_ACC4(0,0,0))*ova2xl(OPS_ACC9(0,0,0)))
    bcl3xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl3xl(OPS_ACC2(0,0,0))
    bcl4xl(OPS_ACC3(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl4xl(OPS_ACC3(0,0,0))
    bcl5xl(OPS_ACC4(0,0,0)) = 0.5_dp*(struxl(OPS_ACC8(0,0,0))+acouxl(OPS_ACC6(0,0,0)))  &
            *(bcl5xl(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
!   L1X UNCHANGED
    bcl2xl(OPS_ACC1(0,0,0)) = -bcl2xl(OPS_ACC1(0,0,0))
    bcl3xl(OPS_ACC2(0,0,0)) = -bcl3xl(OPS_ACC2(0,0,0))
    bcl4xl(OPS_ACC3(0,0,0)) = -bcl4xl(OPS_ACC3(0,0,0))
    bcl5xl(OPS_ACC4(0,0,0)) = 0.5_dp*sorpxl(OPS_ACC10(0,0,0))  &
            + cobcxl*acouxl(OPS_ACC6(0,0,0))*(strpxl(OPS_ACC11(0,0,0))-pinfxl) - bcl5xl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_computeL_xl

SUBROUTINE bounds_kernel_inflowBC1_computeL_yl(bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,bcl1yl,strvyl,ova2yl,sorpyl,strpyl,cobcyl,pinfyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl2yl,bcl3yl,bcl4yl,bcl5yl
    real(kind=dp), dimension(1), intent(in) :: strdyl,acouyl,bcl1yl,strvyl,ova2yl,sorpyl,strpyl
    real(kind=dp), intent(in) :: cobcyl,pinfyl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC5(0,0,0))*acouyl(OPS_ACC6(0,0,0))*bcl1yl(OPS_ACC7(0,0,0))
    bcl2yl(OPS_ACC1(0,0,0)) = strvyl(OPS_ACC8(0,0,0))  &
            *(bcl2yl(OPS_ACC1(0,0,0))-bcl5yl(OPS_ACC4(0,0,0))*ova2yl(OPS_ACC9(0,0,0)))
    bcl3yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl3yl(OPS_ACC2(0,0,0))
    bcl4yl(OPS_ACC3(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl4yl(OPS_ACC3(0,0,0))
    bcl5yl(OPS_ACC4(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC8(0,0,0))+acouyl(OPS_ACC6(0,0,0)))  &
            *(bcl5yl(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
!   L1Y UNCHANGED
    bcl2yl(OPS_ACC1(0,0,0)) = -bcl2yl(OPS_ACC1(0,0,0))
    bcl3yl(OPS_ACC2(0,0,0)) = -bcl3yl(OPS_ACC2(0,0,0))
    bcl4yl(OPS_ACC3(0,0,0)) = -bcl4yl(OPS_ACC3(0,0,0))
    bcl5yl(OPS_ACC4(0,0,0)) = 0.5_dp*sorpyl(OPS_ACC10(0,0,0))  &
            + cobcyl*acouyl(OPS_ACC6(0,0,0))*(strpyl(OPS_ACC11(0,0,0))-pinfyl) - bcl5yl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_computeL_yl

SUBROUTINE bounds_kernel_inflowBC1_computeL_zl(bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,bcl1zl,strwzl,ova2zl,sorpzl,strpzl,cobczl,pinfzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl2zl,bcl3zl,bcl4zl,bcl5zl
    real(kind=dp), dimension(1), intent(in) :: strdzl,acouzl,bcl1zl,strwzl,ova2zl,sorpzl,strpzl
    real(kind=dp), intent(in) :: cobczl,pinfzl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC5(0,0,0))*acouzl(OPS_ACC6(0,0,0))*bcl1zl(OPS_ACC7(0,0,0))
    bcl2zl(OPS_ACC1(0,0,0)) = strwzl(OPS_ACC8(0,0,0))  &
            *(bcl2zl(OPS_ACC1(0,0,0))-bcl5zl(OPS_ACC4(0,0,0))*ova2zl(OPS_ACC9(0,0,0)))
    bcl3zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl3zl(OPS_ACC2(0,0,0))
    bcl4zl(OPS_ACC3(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl4zl(OPS_ACC3(0,0,0))
    bcl5zl(OPS_ACC4(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC8(0,0,0))+acouzl(OPS_ACC6(0,0,0)))  &
            *(bcl5zl(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
!   L1Z UNCHANGED
    bcl2zl(OPS_ACC1(0,0,0)) = -bcl2zl(OPS_ACC1(0,0,0))
    bcl3zl(OPS_ACC2(0,0,0)) = -bcl3zl(OPS_ACC2(0,0,0))
    bcl4zl(OPS_ACC3(0,0,0)) = -bcl4zl(OPS_ACC3(0,0,0))
    bcl5zl(OPS_ACC4(0,0,0)) = 0.5_dp*sorpzl(OPS_ACC10(0,0,0))  &
            + cobczl*acouzl(OPS_ACC6(0,0,0))*(strpzl(OPS_ACC11(0,0,0))-pinfzl) - bcl5zl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_computeL_zl

SUBROUTINE bounds_kernel_inflowBC1_LYX_xl(bclyxl,ratexl,struxl,strdxl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bclyxl
    real(kind=dp), dimension(1), intent(in) :: ratexl,struxl,strdxl
    integer, intent(in) :: ispec    

!   OLD VALUE OF L's
    bclyxl(OPS_ACC_MD1(ispec,0,0,0)) = struxl(OPS_ACC3(0,0,0))*bclyxl(OPS_ACC_MD1(ispec,0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
    bclyxl(OPS_ACC_MD1(ispec,0,0,0)) = ratexl(OPS_ACC_MD2(ispec,0,0,0))/strdxl(OPS_ACC4(0,0,0))  &
              - bclyxl(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_LYX_xl

SUBROUTINE bounds_kernel_inflowBC1_LYY_yl(bclyyl,rateyl,strvyl,strdyl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bclyyl
    real(kind=dp), dimension(1), intent(in) :: rateyl,strvyl,strdyl
    integer, intent(in) :: ispec

!   OLD VALUE OF L's
    bclyyl(OPS_ACC_MD1(ispec,0,0,0)) = strvyl(OPS_ACC3(0,0,0))*bclyyl(OPS_ACC_MD1(ispec,0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
    bclyyl(OPS_ACC_MD1(ispec,0,0,0)) = rateyl(OPS_ACC_MD2(ispec,0,0,0))/strdyl(OPS_ACC4(0,0,0))  &
              - bclyyl(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_LYY_yl

SUBROUTINE bounds_kernel_inflowBC1_LYZ_zl(bclyzl,ratezl,strwzl,strdzl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bclyzl
    real(kind=dp), dimension(1), intent(in) :: ratezl,strwzl,strdzl
    integer, intent(in) :: ispec

!   OLD VALUE OF L's
    bclyzl(OPS_ACC_MD1(ispec,0,0,0)) = strwzl(OPS_ACC3(0,0,0))*bclyzl(OPS_ACC_MD1(ispec,0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
    bclyzl(OPS_ACC_MD1(ispec,0,0,0)) = ratezl(OPS_ACC_MD2(ispec,0,0,0))/strdzl(OPS_ACC4(0,0,0))  &
              - bclyzl(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_LYZ_zl

SUBROUTINE bounds_kernel_inflowBC2_computeL_xl(bcl1xl,bcl2xl,bcl5xl,strdxl,acouxl,struxl,ova2xl,dudtxl,gam1xl,dtdtxl,strtxl,sorpxl,strpxl,sydtxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xl,bcl2xl,bcl5xl
    real(kind=dp), dimension(1), intent(in) :: strdxl,acouxl,struxl,ova2xl,dudtxl,gam1xl,dtdtxl,strtxl,sorpxl,strpxl,sydtxl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC4(0,0,0))*acouxl(OPS_ACC5(0,0,0))*bcl1xl(OPS_ACC1(0,0,0))
    bcl1xl(OPS_ACC1(0,0,0)) = 0.5_dp*(struxl(OPS_ACC6(0,0,0))-acouxl(OPS_ACC5(0,0,0)))  &
            *(bcl5xl(OPS_ACC3(0,0,0))-fornow)
    bcl2xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC6(0,0,0))  &
            *(bcl2xl(OPS_ACC2(0,0,0))-bcl5xl(OPS_ACC3(0,0,0))*ova2xl(OPS_ACC7(0,0,0)))
    bcl5xl(OPS_ACC3(0,0,0)) = 0.5_dp*(struxl(OPS_ACC6(0,0,0))+acouxl(OPS_ACC5(0,0,0)))  &
            *(bcl5xl(OPS_ACC3(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1X UNCHANGED
    bcl5xl(OPS_ACC3(0,0,0)) = bcl1xl(OPS_ACC1(0,0,0))  &
            - strdxl(OPS_ACC4(0,0,0))*acouxl(OPS_ACC5(0,0,0))*dudtxl(OPS_ACC8(0,0,0)) - bcl5xl(OPS_ACC3(0,0,0))
    bcl2xl(OPS_ACC2(0,0,0)) = gam1xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC7(0,0,0))  &
            *(bcl1xl(OPS_ACC1(0,0,0))+bcl5xl(OPS_ACC3(0,0,0)))  &
            + strdxl(OPS_ACC4(0,0,0))*(dtdtxl(OPS_ACC10(0,0,0))/strtxl(OPS_ACC11(0,0,0))  &
            - sorpxl(OPS_ACC12(0,0,0))/strpxl(OPS_ACC13(0,0,0)) + sydtxl(OPS_ACC14(0,0,0)))  &
            - bcl2xl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_computeL_xl

SUBROUTINE bounds_kernel_inflowBC2_computeL_yl(bcl1yl,bcl2yl,bcl5yl,strdyl,acouyl,strvyl,ova2yl,dvdtyl,gam1yl,dtdtyl,strtyl,sorpyl,strpyl,sydtyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yl,bcl2yl,bcl5yl
    real(kind=dp), dimension(1), intent(in) :: strdyl,acouyl,strvyl,ova2yl,dvdtyl,gam1yl,dtdtyl,strtyl,sorpyl,strpyl,sydtyl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC4(0,0,0))*acouyl(OPS_ACC5(0,0,0))*bcl1yl(OPS_ACC1(0,0,0))
    bcl1yl(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC6(0,0,0))-acouyl(OPS_ACC5(0,0,0)))  &
            *(bcl5yl(OPS_ACC3(0,0,0))-fornow)
    bcl2yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC6(0,0,0))  &
            *(bcl2yl(OPS_ACC2(0,0,0))-bcl5yl(OPS_ACC3(0,0,0))*ova2yl(OPS_ACC7(0,0,0)))
    bcl5yl(OPS_ACC3(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC6(0,0,0))+acouyl(OPS_ACC5(0,0,0)))  &
            *(bcl5yl(OPS_ACC3(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Y UNCHANGED
    bcl5yl(OPS_ACC3(0,0,0)) = bcl1yl(OPS_ACC1(0,0,0))  &
            - strdyl(OPS_ACC4(0,0,0))*acouyl(OPS_ACC5(0,0,0))*dvdtyl(OPS_ACC8(0,0,0)) - bcl5yl(OPS_ACC3(0,0,0))
    bcl2yl(OPS_ACC2(0,0,0)) = gam1yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC7(0,0,0))  &
            *(bcl1yl(OPS_ACC1(0,0,0))+bcl5yl(OPS_ACC3(0,0,0)))  &
            + strdyl(OPS_ACC4(0,0,0))*(dtdtyl(OPS_ACC10(0,0,0))/strtyl(OPS_ACC11(0,0,0))  &
            - sorpyl(OPS_ACC12(0,0,0))/strpyl(OPS_ACC13(0,0,0)) + sydtyl(OPS_ACC14(0,0,0)))  &
            - bcl2yl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_computeL_yl

SUBROUTINE bounds_kernel_inflowBC2_computeL_zl(bcl1zl,bcl2zl,bcl5zl,strdzl,acouzl,strwzl,ova2zl,dwdtzl,gam1zl,dtdtzl,strtzl,sorpzl,strpzl,sydtzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zl,bcl2zl,bcl5zl
    real(kind=dp), dimension(1), intent(in) :: strdzl,acouzl,strwzl,ova2zl,dwdtzl,gam1zl,dtdtzl,strtzl,sorpzl,strpzl,sydtzl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC4(0,0,0))*acouzl(OPS_ACC5(0,0,0))*bcl1zl(OPS_ACC1(0,0,0))
    bcl1zl(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC6(0,0,0))-acouzl(OPS_ACC5(0,0,0)))  &
            *(bcl5zl(OPS_ACC3(0,0,0))-fornow)
    bcl2zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC6(0,0,0))  &
            *(bcl2zl(OPS_ACC2(0,0,0))-bcl5zl(OPS_ACC3(0,0,0))*ova2zl(OPS_ACC7(0,0,0)))
    bcl5zl(OPS_ACC3(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC6(0,0,0))+acouzl(OPS_ACC5(0,0,0)))  &
            *(bcl5zl(OPS_ACC3(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Z UNCHANGED
    bcl5zl(OPS_ACC3(0,0,0)) = bcl1zl(OPS_ACC1(0,0,0))  &
            - strdzl(OPS_ACC4(0,0,0))*acouzl(OPS_ACC5(0,0,0))*dwdtzl(OPS_ACC8(0,0,0)) - bcl5zl(OPS_ACC3(0,0,0))
    bcl2zl(OPS_ACC2(0,0,0)) = gam1zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC7(0,0,0))  &
            *(bcl1zl(OPS_ACC1(0,0,0))+bcl5zl(OPS_ACC3(0,0,0)))  &
            + strdzl(OPS_ACC4(0,0,0))*(dtdtzl(OPS_ACC10(0,0,0))/strtzl(OPS_ACC11(0,0,0))  &
            - sorpzl(OPS_ACC12(0,0,0))/strpzl(OPS_ACC13(0,0,0)) + sydtzl(OPS_ACC14(0,0,0)))  &
            - bcl2zl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_computeL_zl

SUBROUTINE bounds_kernel_inflowBC3_computeL_xl(bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,struxl,ova2xl,dudtxl,dddtxl,dvdtxl,dwdtxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl
    real(kind=dp), dimension(1), intent(in) :: strdxl,acouxl,struxl,ova2xl,dudtxl,dddtxl,dvdtxl,dwdtxl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC6(0,0,0))*acouxl(OPS_ACC7(0,0,0))*bcl1xl(OPS_ACC1(0,0,0))
    bcl1xl(OPS_ACC1(0,0,0)) = 0.5_dp*(struxl(OPS_ACC8(0,0,0))-acouxl(OPS_ACC7(0,0,0)))  &
            *(bcl5xl(OPS_ACC5(0,0,0))-fornow)
    bcl2xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC8(0,0,0))  &
            *(bcl2xl(OPS_ACC2(0,0,0))-bcl5xl(OPS_ACC5(0,0,0))*ova2xl(OPS_ACC9(0,0,0)))
    bcl3xl(OPS_ACC3(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl3xl(OPS_ACC3(0,0,0))
    bcl4xl(OPS_ACC4(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl4xl(OPS_ACC4(0,0,0))
    bcl5xl(OPS_ACC5(0,0,0)) = 0.5_dp*(struxl(OPS_ACC8(0,0,0))+acouxl(OPS_ACC7(0,0,0)))  &
            *(bcl5xl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1X UNCHANGED
    fornow = bcl1xl(OPS_ACC1(0,0,0)) - strdxl(OPS_ACC6(0,0,0))*acouxl(OPS_ACC7(0,0,0))*dudtxl(OPS_ACC10(0,0,0))
    bcl2xl(OPS_ACC2(0,0,0)) = -dddtxl(OPS_ACC11(0,0,0))  &
            - ova2xl(OPS_ACC9(0,0,0))*(bcl1xl(OPS_ACC1(0,0,0))+fornow) - bcl2xl(OPS_ACC2(0,0,0))
    bcl3xl(OPS_ACC3(0,0,0)) = -dvdtxl(OPS_ACC12(0,0,0)) - bcl3xl(OPS_ACC3(0,0,0))
    bcl4xl(OPS_ACC4(0,0,0)) = -dwdtxl(OPS_ACC13(0,0,0)) - bcl4xl(OPS_ACC4(0,0,0))
    bcl5xl(OPS_ACC5(0,0,0)) = fornow - bcl5xl(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_computeL_xl

SUBROUTINE bounds_kernel_inflowBC3_computeL_yl(bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,strvyl,ova2yl,dudtyl,dddtyl,dvdtyl,dwdtyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl
    real(kind=dp), dimension(1), intent(in) :: strdyl,acouyl,strvyl,ova2yl,dudtyl,dddtyl,dvdtyl,dwdtyl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC6(0,0,0))*acouyl(OPS_ACC7(0,0,0))*bcl1yl(OPS_ACC1(0,0,0))
    bcl1yl(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC8(0,0,0))-acouyl(OPS_ACC7(0,0,0)))  &
            *(bcl5yl(OPS_ACC5(0,0,0))-fornow)
    bcl2yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC8(0,0,0))  &
            *(bcl2yl(OPS_ACC2(0,0,0))-bcl5yl(OPS_ACC5(0,0,0))*ova2yl(OPS_ACC9(0,0,0)))
    bcl3yl(OPS_ACC3(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl3yl(OPS_ACC3(0,0,0))
    bcl4yl(OPS_ACC4(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl4yl(OPS_ACC4(0,0,0))
    bcl5yl(OPS_ACC5(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC8(0,0,0))+acouyl(OPS_ACC7(0,0,0)))  &
            *(bcl5yl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Y UNCHANGED
    fornow = bcl1yl(OPS_ACC1(0,0,0)) - strdyl(OPS_ACC6(0,0,0))*acouyl(OPS_ACC7(0,0,0))*dvdtyl(OPS_ACC10(0,0,0))
    bcl2yl(OPS_ACC2(0,0,0)) = -dddtyl(OPS_ACC11(0,0,0))  &
            - ova2yl(OPS_ACC9(0,0,0))*(bcl1yl(OPS_ACC1(0,0,0))+fornow) - bcl2yl(OPS_ACC2(0,0,0))
    bcl3yl(OPS_ACC3(0,0,0)) = -dudtyl(OPS_ACC12(0,0,0)) - bcl3yl(OPS_ACC3(0,0,0))
    bcl4yl(OPS_ACC4(0,0,0)) = -dwdtyl(OPS_ACC13(0,0,0)) - bcl4yl(OPS_ACC4(0,0,0))
    bcl5yl(OPS_ACC5(0,0,0)) = fornow - bcl5yl(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_computeL_yl

SUBROUTINE bounds_kernel_inflowBC3_computeL_zl(bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,strwzl,ova2zl,dudtzl,dddtzl,dvdtzl,dwdtzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl
    real(kind=dp), dimension(1), intent(in) :: strdzl,acouzl,strwzl,ova2zl,dudtzl,dddtzl,dvdtzl,dwdtzl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC6(0,0,0))*acouzl(OPS_ACC7(0,0,0))*bcl1zl(OPS_ACC1(0,0,0))
    bcl1zl(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC8(0,0,0))-acouzl(OPS_ACC7(0,0,0)))  &
            *(bcl5zl(OPS_ACC5(0,0,0))-fornow)
    bcl2zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC8(0,0,0))  &
            *(bcl2zl(OPS_ACC2(0,0,0))-bcl5zl(OPS_ACC5(0,0,0))*ova2zl(OPS_ACC9(0,0,0)))
    bcl3zl(OPS_ACC3(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl3zl(OPS_ACC3(0,0,0))
    bcl4zl(OPS_ACC4(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl4zl(OPS_ACC4(0,0,0))
    bcl5zl(OPS_ACC5(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC8(0,0,0))+acouzl(OPS_ACC7(0,0,0)))  &
            *(bcl5zl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Z UNCHANGED
    fornow = bcl1zl(OPS_ACC1(0,0,0)) - strdzl(OPS_ACC6(0,0,0))*acouzl(OPS_ACC7(0,0,0))*dwdtzl(OPS_ACC10(0,0,0))
    bcl2zl(OPS_ACC2(0,0,0)) = -dddtzl(OPS_ACC11(0,0,0))  &
            - ova2zl(OPS_ACC9(0,0,0))*(bcl1zl(OPS_ACC1(0,0,0))+fornow) - bcl2zl(OPS_ACC2(0,0,0))
    bcl3zl(OPS_ACC3(0,0,0)) = -dudtzl(OPS_ACC12(0,0,0)) - bcl3zl(OPS_ACC3(0,0,0))
    bcl4zl(OPS_ACC4(0,0,0)) = -dvdtzl(OPS_ACC13(0,0,0)) - bcl4zl(OPS_ACC4(0,0,0))
    bcl5zl(OPS_ACC5(0,0,0)) = fornow - bcl5zl(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_computeL_zl

SUBROUTINE bounds_kernel_wallBC1_computeL_xl(bcl1xl,bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,struxl,dudtxl,dvdtxl,dwdtxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xl,bcl3xl,bcl4xl,bcl5xl
    real(kind=dp), dimension(1), intent(in) :: strdxl,acouxl,struxl,dudtxl,dvdtxl,dwdtxl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC5(0,0,0))*acouxl(OPS_ACC6(0,0,0))*bcl1xl(OPS_ACC1(0,0,0))
    bcl1xl(OPS_ACC1(0,0,0)) = 0.5_dp*(struxl(OPS_ACC7(0,0,0))-acouxl(OPS_ACC6(0,0,0)))  &
            *(bcl5xl(OPS_ACC4(0,0,0))-fornow)
    bcl3xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC7(0,0,0))*bcl3xl(OPS_ACC2(0,0,0))
    bcl4xl(OPS_ACC3(0,0,0)) = struxl(OPS_ACC7(0,0,0))*bcl4xl(OPS_ACC3(0,0,0))
    bcl5xl(OPS_ACC4(0,0,0)) = 0.5_dp*(struxl(OPS_ACC7(0,0,0))+acouxl(OPS_ACC6(0,0,0)))  &
            *(bcl5xl(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1X,L2X UNCHANGED
    bcl3xl(OPS_ACC2(0,0,0)) = -dvdtxl(OPS_ACC9(0,0,0)) - bcl3xl(OPS_ACC2(0,0,0))
    bcl4xl(OPS_ACC3(0,0,0)) = -dwdtxl(OPS_ACC10(0,0,0)) - bcl4xl(OPS_ACC3(0,0,0))
    bcl5xl(OPS_ACC4(0,0,0)) = bcl1xl(OPS_ACC1(0,0,0))  &
            - strdxl(OPS_ACC5(0,0,0))*acouxl(OPS_ACC6(0,0,0))*dudtxl(OPS_ACC8(0,0,0)) - bcl5xl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_computeL_xl

SUBROUTINE bounds_kernel_wallBC1_computeL_yl(bcl1yl,bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,strvyl,dudtyl,dvdtyl,dwdtyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yl,bcl3yl,bcl4yl,bcl5yl
    real(kind=dp), dimension(1), intent(in) :: strdyl,acouyl,strvyl,dudtyl,dvdtyl,dwdtyl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC5(0,0,0))*acouyl(OPS_ACC6(0,0,0))*bcl1yl(OPS_ACC1(0,0,0))
    bcl1yl(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC7(0,0,0))-acouyl(OPS_ACC6(0,0,0)))  &
            *(bcl5yl(OPS_ACC4(0,0,0))-fornow)
    bcl3yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC7(0,0,0))*bcl3yl(OPS_ACC2(0,0,0))
    bcl4yl(OPS_ACC3(0,0,0)) = strvyl(OPS_ACC7(0,0,0))*bcl4yl(OPS_ACC3(0,0,0))
    bcl5yl(OPS_ACC4(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC7(0,0,0))+acouyl(OPS_ACC6(0,0,0)))  &
            *(bcl5yl(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Y,L2Y UNCHANGED
    bcl3yl(OPS_ACC2(0,0,0)) = -dudtyl(OPS_ACC9(0,0,0)) - bcl3yl(OPS_ACC2(0,0,0))
    bcl4yl(OPS_ACC3(0,0,0)) = -dwdtyl(OPS_ACC10(0,0,0)) - bcl4yl(OPS_ACC3(0,0,0))
    bcl5yl(OPS_ACC4(0,0,0)) = bcl1yl(OPS_ACC1(0,0,0))  &
            - strdyl(OPS_ACC5(0,0,0))*acouyl(OPS_ACC6(0,0,0))*dvdtyl(OPS_ACC8(0,0,0)) - bcl5yl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_computeL_yl

SUBROUTINE bounds_kernel_wallBC1_computeL_zl(bcl1zl,bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,strwzl,dudtzl,dvdtzl,dwdtzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zl,bcl3zl,bcl4zl,bcl5zl
    real(kind=dp), dimension(1), intent(in) :: strdzl,acouzl,strwzl,dudtzl,dvdtzl,dwdtzl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC5(0,0,0))*acouzl(OPS_ACC6(0,0,0))*bcl1zl(OPS_ACC1(0,0,0))
    bcl1zl(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC7(0,0,0))-acouzl(OPS_ACC6(0,0,0)))  &
            *(bcl5zl(OPS_ACC4(0,0,0))-fornow)
    bcl3zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC7(0,0,0))*bcl3zl(OPS_ACC2(0,0,0))
    bcl4zl(OPS_ACC3(0,0,0)) = strwzl(OPS_ACC7(0,0,0))*bcl4zl(OPS_ACC3(0,0,0))
    bcl5zl(OPS_ACC4(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC7(0,0,0))+acouzl(OPS_ACC6(0,0,0)))  &
            *(bcl5zl(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Z,L2Z UNCHANGED
    bcl3zl(OPS_ACC2(0,0,0)) = -dudtzl(OPS_ACC9(0,0,0)) - bcl3zl(OPS_ACC2(0,0,0))
    bcl4zl(OPS_ACC3(0,0,0)) = -dvdtzl(OPS_ACC10(0,0,0)) - bcl4zl(OPS_ACC3(0,0,0))
    bcl5zl(OPS_ACC4(0,0,0)) = bcl1zl(OPS_ACC1(0,0,0))  &
            - strdzl(OPS_ACC5(0,0,0))*acouzl(OPS_ACC6(0,0,0))*dwdtzl(OPS_ACC8(0,0,0)) - bcl5zl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_computeL_zl

SUBROUTINE bounds_kernel_wallBC2_computeL_xl(bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,struxl,gam1xl,ova2xl,dtdtxl,strtxl,sorpxl,strpxl,dudtxl,dvdtxl,dwdtxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl
    real(kind=dp), dimension(1), intent(in) :: strdxl,acouxl,struxl,gam1xl,ova2xl,dtdtxl,strtxl,sorpxl,strpxl,dudtxl,dvdtxl,dwdtxl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC6(0,0,0))*acouxl(OPS_ACC7(0,0,0))*bcl1xl(OPS_ACC1(0,0,0))
    bcl1xl(OPS_ACC1(0,0,0)) = 0.5_dp*(struxl(OPS_ACC8(0,0,0))-acouxl(OPS_ACC7(0,0,0)))  &
            *(bcl5xl(OPS_ACC5(0,0,0))-fornow)
    bcl2xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC8(0,0,0))  &
            *(bcl2xl(OPS_ACC2(0,0,0))-bcl5xl(OPS_ACC5(0,0,0))*ova2xl(OPS_ACC10(0,0,0)))
    bcl3xl(OPS_ACC3(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl3xl(OPS_ACC3(0,0,0))
    bcl4xl(OPS_ACC4(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl4xl(OPS_ACC4(0,0,0))
    bcl5xl(OPS_ACC5(0,0,0)) = 0.5_dp*(struxl(OPS_ACC8(0,0,0))+acouxl(OPS_ACC7(0,0,0)))  &
            *(bcl5xl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1X UNCHANGED
    bcl3xl(OPS_ACC3(0,0,0)) = -dvdtxl(OPS_ACC16(0,0,0)) - bcl3xl(OPS_ACC3(0,0,0))
    bcl4xl(OPS_ACC4(0,0,0)) = -dwdtxl(OPS_ACC17(0,0,0)) - bcl4xl(OPS_ACC4(0,0,0))
    bcl5xl(OPS_ACC5(0,0,0)) = bcl1xl(OPS_ACC1(0,0,0))  &
            - strdxl(OPS_ACC6(0,0,0))*acouxl(OPS_ACC7(0,0,0))*dudtxl(OPS_ACC15(0,0,0)) - bcl5xl(OPS_ACC5(0,0,0))
    bcl2xl(OPS_ACC2(0,0,0)) = gam1xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC10(0,0,0))  &
            *(bcl1xl(OPS_ACC1(0,0,0))+bcl5xl(OPS_ACC5(0,0,0)))  &
            + strdxl(OPS_ACC6(0,0,0))*(dtdtxl(OPS_ACC11(0,0,0))/strtxl(OPS_ACC12(0,0,0))  &
            - sorpxl(OPS_ACC13(0,0,0))/strpxl(OPS_ACC14(0,0,0))) - bcl2xl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_computeL_xl

SUBROUTINE bounds_kernel_wallBC2_computeL_yl(bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,strvyl,gam1yl,ova2yl,dtdtyl,strtyl,sorpyl,strpyl,dudtyl,dvdtyl,dwdtyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl
    real(kind=dp), dimension(1), intent(in) :: strdyl,acouyl,strvyl,gam1yl,ova2yl,dtdtyl,strtyl,sorpyl,strpyl,dudtyl,dvdtyl,dwdtyl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC6(0,0,0))*acouyl(OPS_ACC7(0,0,0))*bcl1yl(OPS_ACC1(0,0,0))
    bcl1yl(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC8(0,0,0))-acouyl(OPS_ACC7(0,0,0)))  &
            *(bcl5yl(OPS_ACC5(0,0,0))-fornow)
    bcl2yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC8(0,0,0))  &
            *(bcl2yl(OPS_ACC2(0,0,0))-bcl5yl(OPS_ACC5(0,0,0))*ova2yl(OPS_ACC10(0,0,0)))
    bcl3yl(OPS_ACC3(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl3yl(OPS_ACC3(0,0,0))
    bcl4yl(OPS_ACC4(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl4yl(OPS_ACC4(0,0,0))
    bcl5yl(OPS_ACC5(0,0,0)) = 0.5_dp*(strvyl(OPS_ACC8(0,0,0))+acouyl(OPS_ACC7(0,0,0)))  &
            *(bcl5yl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Y UNCHANGED
    bcl3yl(OPS_ACC3(0,0,0)) = -dudtyl(OPS_ACC16(0,0,0)) - bcl3yl(OPS_ACC3(0,0,0))
    bcl4yl(OPS_ACC4(0,0,0)) = -dwdtyl(OPS_ACC17(0,0,0)) - bcl4yl(OPS_ACC4(0,0,0))
    bcl5yl(OPS_ACC5(0,0,0)) = bcl1yl(OPS_ACC1(0,0,0))  &
            - strdyl(OPS_ACC6(0,0,0))*acouyl(OPS_ACC7(0,0,0))*dvdtyl(OPS_ACC15(0,0,0)) - bcl5yl(OPS_ACC5(0,0,0))
    bcl2yl(OPS_ACC2(0,0,0)) = gam1yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC10(0,0,0))  &
            *(bcl1yl(OPS_ACC1(0,0,0))+bcl5yl(OPS_ACC5(0,0,0)))  &
            + strdyl(OPS_ACC6(0,0,0))*(dtdtyl(OPS_ACC11(0,0,0))/strtyl(OPS_ACC12(0,0,0))  &
            - sorpyl(OPS_ACC13(0,0,0))/strpyl(OPS_ACC14(0,0,0))) - bcl2yl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_computeL_yl

SUBROUTINE bounds_kernel_wallBC2_computeL_zl(bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,strwzl,gam1zl,ova2zl,dtdtzl,strtzl,sorpzl,strpzl,dudtzl,dvdtzl,dwdtzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl
    real(kind=dp), dimension(1), intent(in) :: strdzl,acouzl,strwzl,gam1zl,ova2zl,dtdtzl,strtzl,sorpzl,strpzl,dudtzl,dvdtzl,dwdtzl
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC6(0,0,0))*acouzl(OPS_ACC7(0,0,0))*bcl1zl(OPS_ACC1(0,0,0))
    bcl1zl(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC8(0,0,0))-acouzl(OPS_ACC7(0,0,0)))  &
            *(bcl5zl(OPS_ACC5(0,0,0))-fornow)
    bcl2zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC8(0,0,0))  &
            *(bcl2zl(OPS_ACC2(0,0,0))-bcl5zl(OPS_ACC5(0,0,0))*ova2zl(OPS_ACC10(0,0,0)))
    bcl3zl(OPS_ACC3(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl3zl(OPS_ACC3(0,0,0))
    bcl4zl(OPS_ACC4(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl4zl(OPS_ACC4(0,0,0))
    bcl5zl(OPS_ACC5(0,0,0)) = 0.5_dp*(strwzl(OPS_ACC8(0,0,0))+acouzl(OPS_ACC7(0,0,0)))  &
            *(bcl5zl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Z UNCHANGED
    bcl3zl(OPS_ACC3(0,0,0)) = -dudtzl(OPS_ACC16(0,0,0)) - bcl3zl(OPS_ACC3(0,0,0))
    bcl4zl(OPS_ACC4(0,0,0)) = -dvdtzl(OPS_ACC17(0,0,0)) - bcl4zl(OPS_ACC4(0,0,0))
    bcl5zl(OPS_ACC5(0,0,0)) = bcl1zl(OPS_ACC1(0,0,0))  &
            - strdzl(OPS_ACC6(0,0,0))*acouzl(OPS_ACC7(0,0,0))*dwdtzl(OPS_ACC15(0,0,0)) - bcl5zl(OPS_ACC5(0,0,0))
    bcl2zl(OPS_ACC2(0,0,0)) = gam1zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC10(0,0,0))  &
            *(bcl1zl(OPS_ACC1(0,0,0))+bcl5zl(OPS_ACC5(0,0,0)))  &
            + strdzl(OPS_ACC6(0,0,0))*(dtdtzl(OPS_ACC11(0,0,0))/strtzl(OPS_ACC12(0,0,0))  &
            - sorpzl(OPS_ACC13(0,0,0))/strpzl(OPS_ACC14(0,0,0))) - bcl2zl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_computeL_zl

SUBROUTINE bounds_kernel_wallBC2_LYX_xl(bclyxl,bcl2xl,ratexl,struxl,strdxl,strrxl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) ::  bclyxl,bcl2xl
    real(kind=dp), dimension(1), intent(in) :: ratexl,struxl,strdxl,strrxl
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

!   OLD VALUE OF LYX
    bclyxl(OPS_ACC_MD1(ispec,0,0,0)) = struxl(OPS_ACC4(0,0,0))*bclyxl(OPS_ACC_MD1(ispec,0,0,0))

!   UPDATE L2X
    bcl2xl(OPS_ACC2(0,0,0)) = bcl2xl(OPS_ACC2(0,0,0)) + (ratexl(OPS_ACC_MD3(ispec,0,0,0))  &
              - strdxl(OPS_ACC5(0,0,0))*bclyxl(OPS_ACC_MD1(ispec,0,0,0))) *rgspec/strrxl(OPS_ACC6(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_LYX_xl

SUBROUTINE bounds_kernel_wallBC2_LYY_yl(bclyyl,bcl2yl,rateyl,strvyl,strdyl,strryl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) ::  bclyyl,bcl2yl
    real(kind=dp), dimension(1), intent(in) :: rateyl,strvyl,strdyl,strryl
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

!   OLD VALUE OF LYY
    bclyyl(OPS_ACC_MD1(ispec,0,0,0)) = strvyl(OPS_ACC4(0,0,0))*bclyyl(OPS_ACC_MD1(ispec,0,0,0))

!   UPDATE L2Y
    bcl2yl(OPS_ACC2(0,0,0)) = bcl2yl(OPS_ACC2(0,0,0)) + (rateyl(OPS_ACC_MD3(ispec,0,0,0))  &
              - strdyl(OPS_ACC5(0,0,0))*bclyyl(OPS_ACC_MD1(ispec,0,0,0))) *rgspec/strryl(OPS_ACC6(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_LYY_yl

SUBROUTINE bounds_kernel_wallBC2_LYZ_zl(bclyzl,bcl2zl,ratezl,strwzl,strdzl,strrzl,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) ::  bclyzl,bcl2zl
    real(kind=dp), dimension(1), intent(in) :: ratezl,strwzl,strdzl,strrzl
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

!   OLD VALUE OF LYZ
    bclyzl(OPS_ACC_MD1(ispec,0,0,0)) = strwzl(OPS_ACC4(0,0,0))*bclyzl(OPS_ACC_MD1(ispec,0,0,0))

!   UPDATE L2Z
    bcl2zl(OPS_ACC2(0,0,0)) = bcl2zl(OPS_ACC2(0,0,0)) + (ratezl(OPS_ACC_MD3(ispec,0,0,0))  &
              - strdzl(OPS_ACC5(0,0,0))*bclyzl(OPS_ACC_MD1(ispec,0,0,0))) *rgspec/strrzl(OPS_ACC6(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_LYZ_zl

SUBROUTINE bounds_kernel_outflowBC1_computeL_xr(bcl1xr,struxr,acouxr,bcl5xr,strdxr,sorpxr,strpxr,cobcxr,pinfxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xr
    real(kind=dp), dimension(1), intent(in) :: struxr,acouxr,bcl5xr,strdxr,sorpxr,strpxr
    real(kind=dp), intent(in) :: cobcxr,pinfxr

!   OLD VALUE OF L1X
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_dp*(struxr(OPS_ACC2(0,0,0))-acouxr(OPS_ACC3(0,0,0)))  &
            *(bcl5xr(OPS_ACC4(0,0,0))-strdxr(OPS_ACC5(0,0,0))*acouxr(OPS_ACC3(0,0,0))*bcl1xr(OPS_ACC1(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L1X
    bcl1xr(OPS_ACC1(0,0,0))= 0.5_dp*sorpxr(OPS_ACC6(0,0,0))  &
            + cobcxr*acouxr(OPS_ACC3(0,0,0))*(strpxr(OPS_ACC7(0,0,0))-pinfxr) - bcl1xr(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_computeL_xr

SUBROUTINE bounds_kernel_outflowBC1_computeL_yr(bcl1yr,strvyr,acouyr,bcl5yr,strdyr,sorpyr,strpyr,cobcyr,pinfyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yr
    real(kind=dp), dimension(1), intent(in) :: strvyr,acouyr,bcl5yr,strdyr,sorpyr,strpyr
    real(kind=dp), intent(in) :: cobcyr,pinfyr

!   OLD VALUE OF L1Y
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC2(0,0,0))-acouyr(OPS_ACC3(0,0,0)))  &
            *(bcl5yr(OPS_ACC4(0,0,0))-strdyr(OPS_ACC5(0,0,0))*acouyr(OPS_ACC3(0,0,0))*bcl1yr(OPS_ACC1(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L1Y
    bcl1yr(OPS_ACC1(0,0,0))= 0.5_dp*sorpyr(OPS_ACC6(0,0,0))  &
            + cobcyr*acouyr(OPS_ACC3(0,0,0))*(strpyr(OPS_ACC7(0,0,0))-pinfyr) - bcl1yr(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_computeL_yr

SUBROUTINE bounds_kernel_outflowBC1_computeL_zr(bcl1zr,strwzr,acouzr,bcl5zr,strdzr,sorpzr,strpzr,cobczr,pinfzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zr
    real(kind=dp), dimension(1), intent(in) :: strwzr,acouzr,bcl5zr,strdzr,sorpzr,strpzr
    real(kind=dp), intent(in) :: cobczr,pinfzr

!   OLD VALUE OF L1Z
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC2(0,0,0))-acouzr(OPS_ACC3(0,0,0)))  &
            *(bcl5zr(OPS_ACC4(0,0,0))-strdzr(OPS_ACC5(0,0,0))*acouzr(OPS_ACC3(0,0,0))*bcl1zr(OPS_ACC1(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L1Z
    bcl1zr(OPS_ACC1(0,0,0))= 0.5_dp*sorpzr(OPS_ACC6(0,0,0))  &
            + cobczr*acouzr(OPS_ACC3(0,0,0))*(strpzr(OPS_ACC7(0,0,0))-pinfzr) - bcl1zr(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_computeL_zr

SUBROUTINE bounds_kernel_inflowBC1_computeL_xr(bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,acouxr,struxr,bcl5xr,ova2xr,sorpxr,strpxr,cobcxr,pinfxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr
    real(kind=dp), dimension(1), intent(in) :: strdxr,acouxr,struxr,bcl5xr,ova2xr,sorpxr,strpxr
    real(kind=dp), intent(in) :: cobcxr,pinfxr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC5(0,0,0))*acouxr(OPS_ACC6(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_dp*(struxr(OPS_ACC7(0,0,0))-acouxr(OPS_ACC6(0,0,0)))  &
            *(bcl5xr(OPS_ACC8(0,0,0))-fornow)
    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC7(0,0,0))  &
            *(bcl2xr(OPS_ACC2(0,0,0))-bcl5xr(OPS_ACC8(0,0,0))*ova2xr(OPS_ACC9(0,0,0)))
    bcl3xr(OPS_ACC3(0,0,0)) = struxr(OPS_ACC7(0,0,0))*bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = struxr(OPS_ACC7(0,0,0))*bcl4xr(OPS_ACC4(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
!   L5X UNCHANGED
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_dp*sorpxr(OPS_ACC10(0,0,0))  &
            + cobcxr*acouxr(OPS_ACC6(0,0,0))*(strpxr(OPS_ACC11(0,0,0))-pinfxr) - bcl1xr(OPS_ACC1(0,0,0))
    bcl2xr(OPS_ACC2(0,0,0)) = -bcl2xr(OPS_ACC2(0,0,0))
    bcl3xr(OPS_ACC3(0,0,0)) = -bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = -bcl4xr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_computeL_xr

SUBROUTINE bounds_kernel_inflowBC1_computeL_yr(bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,acouyr,strvyr,bcl5yr,ova2yr,sorpyr,strpyr,cobcyr,pinfyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr
    real(kind=dp), dimension(1), intent(in) :: strdyr,acouyr,strvyr,bcl5yr,ova2yr,sorpyr,strpyr
    real(kind=dp), intent(in) :: cobcyr,pinfyr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC5(0,0,0))*acouyr(OPS_ACC6(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC7(0,0,0))-acouyr(OPS_ACC6(0,0,0)))  &
            *(bcl5yr(OPS_ACC8(0,0,0))-fornow)
    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC7(0,0,0))  &
            *(bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC8(0,0,0))*ova2yr(OPS_ACC9(0,0,0)))
    bcl3yr(OPS_ACC3(0,0,0)) = strvyr(OPS_ACC7(0,0,0))*bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = strvyr(OPS_ACC7(0,0,0))*bcl4yr(OPS_ACC4(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
!   L5Y UNCHANGED
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_dp*sorpyr(OPS_ACC10(0,0,0))  &
            + cobcyr*acouyr(OPS_ACC6(0,0,0))*(strpyr(OPS_ACC11(0,0,0))-pinfyr) - bcl1yr(OPS_ACC1(0,0,0))
    bcl2yr(OPS_ACC2(0,0,0)) = -bcl2yr(OPS_ACC2(0,0,0))
    bcl3yr(OPS_ACC3(0,0,0)) = -bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = -bcl4yr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_computeL_yr

SUBROUTINE bounds_kernel_inflowBC1_computeL_zr(bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,acouzr,strwzr,bcl5zr,ova2zr,sorpzr,strpzr,cobczr,pinfzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr
    real(kind=dp), dimension(1), intent(in) :: strdzr,acouzr,strwzr,bcl5zr,ova2zr,sorpzr,strpzr
    real(kind=dp), intent(in) :: cobczr,pinfzr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC5(0,0,0))*acouzr(OPS_ACC6(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC7(0,0,0))-acouzr(OPS_ACC6(0,0,0)))  &
            *(bcl5zr(OPS_ACC8(0,0,0))-fornow)
    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC7(0,0,0))  &
            *(bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC8(0,0,0))*ova2zr(OPS_ACC9(0,0,0)))
    bcl3zr(OPS_ACC3(0,0,0)) = strwzr(OPS_ACC7(0,0,0))*bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = strwzr(OPS_ACC7(0,0,0))*bcl4zr(OPS_ACC4(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
!   L5Z UNCHANGED
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_dp*sorpzr(OPS_ACC10(0,0,0))  &
            + cobczr*acouzr(OPS_ACC6(0,0,0))*(strpzr(OPS_ACC11(0,0,0))-pinfzr) - bcl1zr(OPS_ACC1(0,0,0))
    bcl2zr(OPS_ACC2(0,0,0)) = -bcl2zr(OPS_ACC2(0,0,0))
    bcl3zr(OPS_ACC3(0,0,0)) = -bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = -bcl4zr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_computeL_zr

SUBROUTINE bounds_kernel_inflowBC1_LYX_xr(bclyxr,ratexr,struxr,strdxr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bclyxr
    real(kind=dp), dimension(1), intent(in) :: ratexr,struxr,strdxr
    integer, intent(in) :: ispec

!   OLD VALUE OF L's
    bclyxr(OPS_ACC_MD1(ispec,0,0,0)) = struxr(OPS_ACC3(0,0,0))*bclyxr(OPS_ACC_MD1(ispec,0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
    bclyxr(OPS_ACC_MD1(ispec,0,0,0)) = ratexr(OPS_ACC_MD2(ispec,0,0,0))/strdxr(OPS_ACC4(0,0,0))  &
              - bclyxr(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_LYX_xr

SUBROUTINE bounds_kernel_inflowBC1_LYY_yr(bclyyr,rateyr,strvyr,strdyr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bclyyr
    real(kind=dp), dimension(1), intent(in) :: rateyr,strvyr,strdyr
    integer, intent(in) :: ispec

!   OLD VALUE OF L's
    bclyyr(OPS_ACC_MD1(ispec,0,0,0)) = strvyr(OPS_ACC3(0,0,0))*bclyyr(OPS_ACC_MD1(ispec,0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
    bclyyr(OPS_ACC_MD1(ispec,0,0,0)) = rateyr(OPS_ACC_MD2(ispec,0,0,0))/strdyr(OPS_ACC4(0,0,0))  &
              - bclyyr(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_LYY_yr

SUBROUTINE bounds_kernel_inflowBC1_LYZ_zr(bclyzr,ratezr,strwzr,strdzr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bclyzr
    real(kind=dp), dimension(1), intent(in) :: ratezr,strwzr,strdzr
    integer, intent(in) :: ispec

!   OLD VALUE OF L's
    bclyzr(OPS_ACC_MD1(ispec,0,0,0)) = strwzr(OPS_ACC3(0,0,0))*bclyzr(OPS_ACC_MD1(ispec,0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
    bclyzr(OPS_ACC_MD1(ispec,0,0,0)) = ratezr(OPS_ACC_MD2(ispec,0,0,0))/strdzr(OPS_ACC4(0,0,0))  &
              - bclyzr(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_LYZ_zr

SUBROUTINE bounds_kernel_inflowBC2_computeL_xr(bcl1xr,bcl2xr,bcl5xr,strdxr,acouxr,struxr,ova2xr,dudtxr,gam1xr,dtdtxr,strtxr,sorpxr,strpxr,sydtxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xr,bcl2xr,bcl5xr
    real(kind=dp), dimension(1), intent(in) :: strdxr,acouxr,struxr,ova2xr,dudtxr,gam1xr,dtdtxr,strtxr,sorpxr,strpxr,sydtxr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC4(0,0,0))*acouxr(OPS_ACC5(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_dp*(struxr(OPS_ACC6(0,0,0))-acouxr(OPS_ACC5(0,0,0)))  &
            *(bcl5xr(OPS_ACC3(0,0,0))-fornow)
    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC6(0,0,0))  &
            *(bcl2xr(OPS_ACC2(0,0,0))-bcl5xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC7(0,0,0)))
    bcl5xr(OPS_ACC3(0,0,0)) = 0.5_dp*(struxr(OPS_ACC6(0,0,0))+acouxr(OPS_ACC5(0,0,0)))  &
            *(bcl5xr(OPS_ACC3(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5X UNCHANGED
    bcl1xr(OPS_ACC1(0,0,0)) = bcl5xr(OPS_ACC3(0,0,0))  &
            + strdxr(OPS_ACC4(0,0,0))*acouxr(OPS_ACC5(0,0,0))*dudtxr(OPS_ACC8(0,0,0)) - bcl1xr(OPS_ACC1(0,0,0))
    bcl2xr(OPS_ACC2(0,0,0)) = gam1xr(OPS_ACC9(0,0,0))*ova2xr(OPS_ACC7(0,0,0))  &
            *(bcl1xr(OPS_ACC1(0,0,0))+bcl5xr(OPS_ACC3(0,0,0)))  &
            + strdxr(OPS_ACC4(0,0,0))*(dtdtxr(OPS_ACC10(0,0,0))/strtxr(OPS_ACC11(0,0,0))  &
            - sorpxr(OPS_ACC12(0,0,0))/strpxr(OPS_ACC13(0,0,0)) + sydtxr(OPS_ACC14(0,0,0)))  &
            - bcl2xr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_computeL_xr

SUBROUTINE bounds_kernel_inflowBC2_computeL_yr(bcl1yr,bcl2yr,bcl5yr,strdyr,acouyr,strvyr,ova2yr,dvdtyr,gam1yr,dtdtyr,strtyr,sorpyr,strpyr,sydtyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yr,bcl2yr,bcl5yr
    real(kind=dp), dimension(1), intent(in) :: strdyr,acouyr,strvyr,ova2yr,dvdtyr,gam1yr,dtdtyr,strtyr,sorpyr,strpyr,sydtyr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC4(0,0,0))*acouyr(OPS_ACC5(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC6(0,0,0))-acouyr(OPS_ACC5(0,0,0)))  &
            *(bcl5yr(OPS_ACC3(0,0,0))-fornow)
    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC6(0,0,0))  &
            *(bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC7(0,0,0)))
    bcl5yr(OPS_ACC3(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC6(0,0,0))+acouyr(OPS_ACC5(0,0,0)))  &
            *(bcl5yr(OPS_ACC3(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Y UNCHANGED
    bcl1yr(OPS_ACC1(0,0,0)) = bcl5yr(OPS_ACC3(0,0,0))  &
            + strdyr(OPS_ACC4(0,0,0))*acouyr(OPS_ACC5(0,0,0))*dvdtyr(OPS_ACC8(0,0,0)) - bcl1yr(OPS_ACC1(0,0,0))
    bcl2yr(OPS_ACC2(0,0,0)) = gam1yr(OPS_ACC9(0,0,0))*ova2yr(OPS_ACC7(0,0,0))  &
            *(bcl1yr(OPS_ACC1(0,0,0))+bcl5yr(OPS_ACC3(0,0,0)))  &
            + strdyr(OPS_ACC4(0,0,0))*(dtdtyr(OPS_ACC10(0,0,0))/strtyr(OPS_ACC11(0,0,0))  &
            - sorpyr(OPS_ACC12(0,0,0))/strpyr(OPS_ACC13(0,0,0)) + sydtyr(OPS_ACC14(0,0,0)))  &
            - bcl2yr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_computeL_yr

SUBROUTINE bounds_kernel_inflowBC2_computeL_zr(bcl1zr,bcl2zr,bcl5zr,strdzr,acouzr,strwzr,ova2zr,dwdtzr,gam1zr,dtdtzr,strtzr,sorpzr,strpzr,sydtzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zr,bcl2zr,bcl5zr
    real(kind=dp), dimension(1), intent(in) :: strdzr,acouzr,strwzr,ova2zr,dwdtzr,gam1zr,dtdtzr,strtzr,sorpzr,strpzr,sydtzr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC4(0,0,0))*acouzr(OPS_ACC5(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC6(0,0,0))-acouzr(OPS_ACC5(0,0,0)))  &
            *(bcl5zr(OPS_ACC3(0,0,0))-fornow)
    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC6(0,0,0))  &
            *(bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC7(0,0,0)))
    bcl5zr(OPS_ACC3(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC6(0,0,0))+acouzr(OPS_ACC5(0,0,0)))  &
            *(bcl5zr(OPS_ACC3(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Z UNCHANGED
    bcl1zr(OPS_ACC1(0,0,0)) = bcl5zr(OPS_ACC3(0,0,0))  &
            + strdzr(OPS_ACC4(0,0,0))*acouzr(OPS_ACC5(0,0,0))*dwdtzr(OPS_ACC8(0,0,0)) - bcl1zr(OPS_ACC1(0,0,0))
    bcl2zr(OPS_ACC2(0,0,0)) = gam1zr(OPS_ACC9(0,0,0))*ova2zr(OPS_ACC7(0,0,0))  &
            *(bcl1zr(OPS_ACC1(0,0,0))+bcl5zr(OPS_ACC3(0,0,0)))  &
            + strdzr(OPS_ACC4(0,0,0))*(dtdtzr(OPS_ACC10(0,0,0))/strtzr(OPS_ACC11(0,0,0))  &
            - sorpzr(OPS_ACC12(0,0,0))/strpzr(OPS_ACC13(0,0,0)) + sydtzr(OPS_ACC14(0,0,0)))  &
            - bcl2zr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_computeL_zr

SUBROUTINE bounds_kernel_inflowBC3_computeL_xr(bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr,strdxr,acouxr,struxr,ova2xr,dudtxr,dddtxr,dvdtxr,dwdtxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr
    real(kind=dp), dimension(1), intent(in) :: strdxr,acouxr,struxr,ova2xr,dudtxr,dddtxr,dvdtxr,dwdtxr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC6(0,0,0))*acouxr(OPS_ACC7(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_dp*(struxr(OPS_ACC8(0,0,0))-acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))-fornow)
    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC8(0,0,0))  &
            *(bcl2xr(OPS_ACC2(0,0,0))-bcl5xr(OPS_ACC5(0,0,0))*ova2xr(OPS_ACC9(0,0,0)))
    bcl3xr(OPS_ACC3(0,0,0)) = struxr(OPS_ACC8(0,0,0))*bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = struxr(OPS_ACC8(0,0,0))*bcl4xr(OPS_ACC4(0,0,0))
    bcl5xr(OPS_ACC5(0,0,0)) = 0.5_dp*(struxr(OPS_ACC8(0,0,0))+acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5X UNCHANGED
    fornow = bcl5xr(OPS_ACC5(0,0,0)) + strdxr(OPS_ACC6(0,0,0))*acouxr(OPS_ACC7(0,0,0))*dudtxr(OPS_ACC10(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = fornow - bcl1xr(OPS_ACC1(0,0,0))
    bcl2xr(OPS_ACC2(0,0,0)) = -dddtxr(OPS_ACC11(0,0,0))  &
            - ova2xr(OPS_ACC9(0,0,0))*(bcl5xr(OPS_ACC5(0,0,0))+fornow) - bcl2xr(OPS_ACC2(0,0,0))
    bcl3xr(OPS_ACC3(0,0,0)) = -dvdtxr(OPS_ACC12(0,0,0)) - bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = -dwdtxr(OPS_ACC13(0,0,0)) - bcl4xr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_computeL_xr

SUBROUTINE bounds_kernel_inflowBC3_computeL_yr(bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr,strdyr,acouyr,strvyr,ova2yr,dudtyr,dddtyr,dvdtyr,dwdtyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr
    real(kind=dp), dimension(1), intent(in) :: strdyr,acouyr,strvyr,ova2yr,dudtyr,dddtyr,dvdtyr,dwdtyr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC6(0,0,0))*acouyr(OPS_ACC7(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC8(0,0,0))-acouyr(OPS_ACC7(0,0,0)))  &
            *(bcl5yr(OPS_ACC5(0,0,0))-fornow)
    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC8(0,0,0))  &
            *(bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC5(0,0,0))*ova2yr(OPS_ACC9(0,0,0)))
    bcl3yr(OPS_ACC3(0,0,0)) = strvyr(OPS_ACC8(0,0,0))*bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = strvyr(OPS_ACC8(0,0,0))*bcl4yr(OPS_ACC4(0,0,0))
    bcl5yr(OPS_ACC5(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC8(0,0,0))+acouyr(OPS_ACC7(0,0,0)))  &
            *(bcl5yr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Y UNCHANGED
    fornow = bcl5yr(OPS_ACC5(0,0,0)) + strdyr(OPS_ACC6(0,0,0))*acouyr(OPS_ACC7(0,0,0))*dvdtyr(OPS_ACC10(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = fornow - bcl1yr(OPS_ACC1(0,0,0))
    bcl2yr(OPS_ACC2(0,0,0)) = -dddtyr(OPS_ACC11(0,0,0))  &
            - ova2yr(OPS_ACC9(0,0,0))*(bcl5yr(OPS_ACC5(0,0,0))+fornow) - bcl2yr(OPS_ACC2(0,0,0))
    bcl3yr(OPS_ACC3(0,0,0)) = -dudtyr(OPS_ACC12(0,0,0)) - bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = -dwdtyr(OPS_ACC13(0,0,0)) - bcl4yr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_computeL_yr

SUBROUTINE bounds_kernel_inflowBC3_computeL_zr(bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr,strdzr,acouzr,strwzr,ova2zr,dudtzr,dddtzr,dvdtzr,dwdtzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr
    real(kind=dp), dimension(1), intent(in) :: strdzr,acouzr,strwzr,ova2zr,dudtzr,dddtzr,dvdtzr,dwdtzr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC6(0,0,0))*acouzr(OPS_ACC7(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC8(0,0,0))-acouzr(OPS_ACC7(0,0,0)))  &
            *(bcl5zr(OPS_ACC5(0,0,0))-fornow)
    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC8(0,0,0))  &
            *(bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC5(0,0,0))*ova2zr(OPS_ACC9(0,0,0)))
    bcl3zr(OPS_ACC3(0,0,0)) = strwzr(OPS_ACC8(0,0,0))*bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = strwzr(OPS_ACC8(0,0,0))*bcl4zr(OPS_ACC4(0,0,0))
    bcl5zr(OPS_ACC5(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC8(0,0,0))+acouzr(OPS_ACC7(0,0,0)))  &
            *(bcl5zr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Z UNCHANGED
    fornow = bcl5zr(OPS_ACC5(0,0,0)) + strdzr(OPS_ACC6(0,0,0))*acouzr(OPS_ACC7(0,0,0))*dwdtzr(OPS_ACC10(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = fornow - bcl1zr(OPS_ACC1(0,0,0))
    bcl2zr(OPS_ACC2(0,0,0)) = -dddtzr(OPS_ACC11(0,0,0))  &
            - ova2zr(OPS_ACC9(0,0,0))*(bcl5zr(OPS_ACC5(0,0,0))+fornow) - bcl2zr(OPS_ACC2(0,0,0))
    bcl3zr(OPS_ACC3(0,0,0)) = -dudtzr(OPS_ACC12(0,0,0)) - bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = -dvdtzr(OPS_ACC13(0,0,0)) - bcl4zr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_computeL_zr

SUBROUTINE bounds_kernel_wallBC1_computeL_xr(bcl1xr,bcl3xr,bcl4xr,bcl5xr,strdxr,acouxr,struxr,dudtxr,dvdtxr,dwdtxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xr,bcl3xr,bcl4xr,bcl5xr
    real(kind=dp), dimension(1), intent(in) :: strdxr,acouxr,struxr,dudtxr,dvdtxr,dwdtxr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC5(0,0,0))*acouxr(OPS_ACC6(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_dp*(struxr(OPS_ACC7(0,0,0))-acouxr(OPS_ACC6(0,0,0)))  &
            *(bcl5xr(OPS_ACC4(0,0,0))-fornow)
    bcl3xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC7(0,0,0))*bcl3xr(OPS_ACC2(0,0,0))
    bcl4xr(OPS_ACC3(0,0,0)) = struxr(OPS_ACC7(0,0,0))*bcl4xr(OPS_ACC3(0,0,0))
    bcl5xr(OPS_ACC4(0,0,0)) = 0.5_dp*(struxr(OPS_ACC7(0,0,0))+acouxr(OPS_ACC6(0,0,0)))  &
            *(bcl5xr(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L2X,L5X UNCHANGED
    bcl1xr(OPS_ACC1(0,0,0)) = bcl5xr(OPS_ACC4(0,0,0))  &
            + strdxr(OPS_ACC5(0,0,0))*acouxr(OPS_ACC6(0,0,0))*dudtxr(OPS_ACC8(0,0,0)) - bcl1xr(OPS_ACC1(0,0,0))
    bcl3xr(OPS_ACC2(0,0,0)) = -dvdtxr(OPS_ACC9(0,0,0)) - bcl3xr(OPS_ACC2(0,0,0))
    bcl4xr(OPS_ACC3(0,0,0)) = -dwdtxr(OPS_ACC10(0,0,0)) - bcl4xr(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_computeL_xr

SUBROUTINE bounds_kernel_wallBC1_computeL_yr(bcl1yr,bcl3yr,bcl4yr,bcl5yr,strdyr,acouyr,strvyr,dudtyr,dvdtyr,dwdtyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yr,bcl3yr,bcl4yr,bcl5yr
    real(kind=dp), dimension(1), intent(in) :: strdyr,acouyr,strvyr,dudtyr,dvdtyr,dwdtyr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC5(0,0,0))*acouyr(OPS_ACC6(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC7(0,0,0))-acouyr(OPS_ACC6(0,0,0)))  &
            *(bcl5yr(OPS_ACC4(0,0,0))-fornow)
    bcl3yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC7(0,0,0))*bcl3yr(OPS_ACC2(0,0,0))
    bcl4yr(OPS_ACC3(0,0,0)) = strvyr(OPS_ACC7(0,0,0))*bcl4yr(OPS_ACC3(0,0,0))
    bcl5yr(OPS_ACC4(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC7(0,0,0))+acouyr(OPS_ACC6(0,0,0)))  &
            *(bcl5yr(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L2Y,L5Y UNCHANGED
    bcl1yr(OPS_ACC1(0,0,0)) = bcl5yr(OPS_ACC4(0,0,0))  &
            + strdyr(OPS_ACC5(0,0,0))*acouyr(OPS_ACC6(0,0,0))*dvdtyr(OPS_ACC8(0,0,0)) - bcl1yr(OPS_ACC1(0,0,0))
    bcl3yr(OPS_ACC2(0,0,0)) = -dudtyr(OPS_ACC9(0,0,0)) - bcl3yr(OPS_ACC2(0,0,0))
    bcl4yr(OPS_ACC3(0,0,0)) = -dwdtyr(OPS_ACC10(0,0,0)) - bcl4yr(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_computeL_yr

SUBROUTINE bounds_kernel_wallBC1_computeL_zr(bcl1zr,bcl3zr,bcl4zr,bcl5zr,strdzr,acouzr,strwzr,dudtzr,dvdtzr,dwdtzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zr,bcl3zr,bcl4zr,bcl5zr
    real(kind=dp), dimension(1), intent(in) :: strdzr,acouzr,strwzr,dudtzr,dvdtzr,dwdtzr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC5(0,0,0))*acouzr(OPS_ACC6(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC7(0,0,0))-acouzr(OPS_ACC6(0,0,0)))  &
            *(bcl5zr(OPS_ACC4(0,0,0))-fornow)
    bcl3zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC7(0,0,0))*bcl3zr(OPS_ACC2(0,0,0))
    bcl4zr(OPS_ACC3(0,0,0)) = strwzr(OPS_ACC7(0,0,0))*bcl4zr(OPS_ACC3(0,0,0))
    bcl5zr(OPS_ACC4(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC7(0,0,0))+acouzr(OPS_ACC6(0,0,0)))  &
            *(bcl5zr(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L2Z,L5Z UNCHANGED
    bcl1zr(OPS_ACC1(0,0,0)) = bcl5zr(OPS_ACC4(0,0,0))  &
            + strdzr(OPS_ACC5(0,0,0))*acouzr(OPS_ACC6(0,0,0))*dwdtzr(OPS_ACC8(0,0,0)) - bcl1zr(OPS_ACC1(0,0,0))
    bcl3zr(OPS_ACC2(0,0,0)) = -dudtzr(OPS_ACC9(0,0,0)) - bcl3zr(OPS_ACC2(0,0,0))
    bcl4zr(OPS_ACC3(0,0,0)) = -dvdtzr(OPS_ACC10(0,0,0)) - bcl4zr(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_computeL_zr

SUBROUTINE bounds_kernel_wallBC2_computeL_xr(bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr,strdxr,acouxr,struxr,ova2xr,dudtxr,dvdtxr,dwdtxr,gam1xr,dtdtxr,strtxr,sorpxr,strpxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr
    real(kind=dp), dimension(1), intent(in) :: strdxr,acouxr,struxr,ova2xr,dudtxr,dvdtxr,dwdtxr,gam1xr,dtdtxr,strtxr,sorpxr,strpxr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC6(0,0,0))*acouxr(OPS_ACC7(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_dp*(struxr(OPS_ACC8(0,0,0))-acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))-fornow)
    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC8(0,0,0))  &
            *(bcl2xr(OPS_ACC2(0,0,0))-bcl5xr(OPS_ACC5(0,0,0))*ova2xr(OPS_ACC9(0,0,0)))
    bcl3xr(OPS_ACC3(0,0,0)) = struxr(OPS_ACC8(0,0,0))*bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = struxr(OPS_ACC8(0,0,0))*bcl4xr(OPS_ACC4(0,0,0))
    bcl5xr(OPS_ACC5(0,0,0)) = 0.5_dp*(struxr(OPS_ACC8(0,0,0))+acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5X UNCHANGED
    bcl1xr(OPS_ACC1(0,0,0)) = bcl5xr(OPS_ACC5(0,0,0))  &
            + strdxr(OPS_ACC6(0,0,0))*acouxr(OPS_ACC7(0,0,0))*dudtxr(OPS_ACC10(0,0,0)) - bcl1xr(OPS_ACC1(0,0,0))
    bcl3xr(OPS_ACC3(0,0,0)) = -dvdtxr(OPS_ACC11(0,0,0)) - bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = -dwdtxr(OPS_ACC12(0,0,0)) - bcl4xr(OPS_ACC4(0,0,0))
    bcl2xr(OPS_ACC2(0,0,0)) = gam1xr(OPS_ACC13(0,0,0))*ova2xr(OPS_ACC9(0,0,0))  &
            *(bcl1xr(OPS_ACC1(0,0,0))+bcl5xr(OPS_ACC5(0,0,0)))  &
            + strdxr(OPS_ACC6(0,0,0))*(dtdtxr(OPS_ACC14(0,0,0))/strtxr(OPS_ACC15(0,0,0))  &
            - sorpxr(OPS_ACC16(0,0,0))/strpxr(OPS_ACC17(0,0,0))) - bcl2xr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_computeL_xr

SUBROUTINE bounds_kernel_wallBC2_computeL_yr(bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr,strdyr,acouyr,strvyr,ova2yr,dudtyr,dvdtyr,dwdtyr,gam1yr,dtdtyr,strtyr,sorpyr,strpyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr
    real(kind=dp), dimension(1), intent(in) :: strdyr,acouyr,strvyr,ova2yr,dudtyr,dvdtyr,dwdtyr,gam1yr,dtdtyr,strtyr,sorpyr,strpyr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC6(0,0,0))*acouyr(OPS_ACC7(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC8(0,0,0))-acouyr(OPS_ACC7(0,0,0)))  &
            *(bcl5yr(OPS_ACC5(0,0,0))-fornow)
    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC8(0,0,0))  &
            *(bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC5(0,0,0))*ova2yr(OPS_ACC9(0,0,0)))
    bcl3yr(OPS_ACC3(0,0,0)) = strvyr(OPS_ACC8(0,0,0))*bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = strvyr(OPS_ACC8(0,0,0))*bcl4yr(OPS_ACC4(0,0,0))
    bcl5yr(OPS_ACC5(0,0,0)) = 0.5_dp*(strvyr(OPS_ACC8(0,0,0))+acouyr(OPS_ACC7(0,0,0)))  &
            *(bcl5yr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Y UNCHANGED
    bcl1yr(OPS_ACC1(0,0,0)) = bcl5yr(OPS_ACC5(0,0,0))  &
            + strdyr(OPS_ACC6(0,0,0))*acouyr(OPS_ACC7(0,0,0))*dvdtyr(OPS_ACC10(0,0,0)) - bcl1yr(OPS_ACC1(0,0,0))
    bcl3yr(OPS_ACC3(0,0,0)) = -dudtyr(OPS_ACC11(0,0,0)) - bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = -dwdtyr(OPS_ACC12(0,0,0)) - bcl4yr(OPS_ACC4(0,0,0))
    bcl2yr(OPS_ACC2(0,0,0)) = gam1yr(OPS_ACC13(0,0,0))*ova2yr(OPS_ACC9(0,0,0))  &
            *(bcl1yr(OPS_ACC1(0,0,0))+bcl5yr(OPS_ACC5(0,0,0)))  &
            + strdyr(OPS_ACC6(0,0,0))*(dtdtyr(OPS_ACC14(0,0,0))/strtyr(OPS_ACC15(0,0,0))  &
            - sorpyr(OPS_ACC16(0,0,0))/strpyr(OPS_ACC17(0,0,0))) - bcl2yr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_computeL_yr

SUBROUTINE bounds_kernel_wallBC2_computeL_zr(bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr,strdzr,acouzr,strwzr,ova2zr,dudtzr,dvdtzr,dwdtzr,gam1zr,dtdtzr,strtzr,sorpzr,strpzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr
    real(kind=dp), dimension(1), intent(in) :: strdzr,acouzr,strwzr,ova2zr,dudtzr,dvdtzr,dwdtzr,gam1zr,dtdtzr,strtzr,sorpzr,strpzr
    real(kind=dp) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC6(0,0,0))*acouzr(OPS_ACC7(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC8(0,0,0))-acouzr(OPS_ACC7(0,0,0)))  &
            *(bcl5zr(OPS_ACC5(0,0,0))-fornow)
    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC8(0,0,0))  &
            *(bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC5(0,0,0))*ova2zr(OPS_ACC9(0,0,0)))
    bcl3zr(OPS_ACC3(0,0,0)) = strwzr(OPS_ACC8(0,0,0))*bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = strwzr(OPS_ACC8(0,0,0))*bcl4zr(OPS_ACC4(0,0,0))
    bcl5zr(OPS_ACC5(0,0,0)) = 0.5_dp*(strwzr(OPS_ACC8(0,0,0))+acouzr(OPS_ACC7(0,0,0)))  &
            *(bcl5zr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Z UNCHANGED
    bcl1zr(OPS_ACC1(0,0,0)) = bcl5zr(OPS_ACC5(0,0,0))  &
            + strdzr(OPS_ACC6(0,0,0))*acouzr(OPS_ACC7(0,0,0))*dwdtzr(OPS_ACC10(0,0,0)) - bcl1zr(OPS_ACC1(0,0,0))
    bcl3zr(OPS_ACC3(0,0,0)) = -dudtzr(OPS_ACC11(0,0,0)) - bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = -dvdtzr(OPS_ACC12(0,0,0)) - bcl4zr(OPS_ACC4(0,0,0))
    bcl2zr(OPS_ACC2(0,0,0)) = gam1zr(OPS_ACC13(0,0,0))*ova2zr(OPS_ACC9(0,0,0))  &
            *(bcl1zr(OPS_ACC1(0,0,0))+bcl5zr(OPS_ACC5(0,0,0)))  &
            + strdzr(OPS_ACC6(0,0,0))*(dtdtzr(OPS_ACC14(0,0,0))/strtzr(OPS_ACC15(0,0,0))  &
            - sorpzr(OPS_ACC16(0,0,0))/strpzr(OPS_ACC17(0,0,0))) - bcl2zr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_computeL_zr

SUBROUTINE bounds_kernel_wallBC2_LYX_xr(bclyxr,bcl2xr,ratexr,struxr,strdxr,strrxr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bclyxr,bcl2xr
    real(kind=dp), dimension(1), intent(in) :: ratexr,struxr,strdxr,strrxr
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

!   OLD VALUE OF LYX
    bclyxr(OPS_ACC_MD1(ispec,0,0,0)) = struxr(OPS_ACC4(0,0,0))*bclyxr(OPS_ACC_MD1(ispec,0,0,0))

!   UPDATE L2X
    bcl2xr(OPS_ACC2(0,0,0)) = bcl2xr(OPS_ACC2(0,0,0)) + (ratexr(OPS_ACC_MD3(ispec,0,0,0))  &
              - strdxr(OPS_ACC5(0,0,0))*bclyxr(OPS_ACC_MD1(ispec,0,0,0))) *rgspec/strrxr(OPS_ACC6(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_LYX_xr

SUBROUTINE bounds_kernel_wallBC2_LYY_yr(bclyyr,bcl2yr,rateyr,strvyr,strdyr,strryr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bclyyr,bcl2yr
    real(kind=dp), dimension(1), intent(in) :: rateyr,strvyr,strdyr,strryr
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

!   OLD VALUE OF LYY
    bclyyr(OPS_ACC_MD1(ispec,0,0,0)) = strvyr(OPS_ACC4(0,0,0))*bclyyr(OPS_ACC_MD1(ispec,0,0,0))

!   UPDATE L2Y
    bcl2yr(OPS_ACC2(0,0,0)) = bcl2yr(OPS_ACC2(0,0,0)) + (rateyr(OPS_ACC_MD3(ispec,0,0,0))  &
              - strdyr(OPS_ACC5(0,0,0))*bclyyr(OPS_ACC_MD1(ispec,0,0,0))) *rgspec/strryr(OPS_ACC6(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_LYY_yr

SUBROUTINE bounds_kernel_wallBC2_LYZ_zr(bclyzr,bcl2zr,ratezr,strwzr,strdzr,strrzr,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: bclyzr,bcl2zr
    real(kind=dp), dimension(1), intent(in) :: ratezr,strwzr,strdzr,strrzr
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec

!   OLD VALUE OF LYZ
    bclyzr(OPS_ACC_MD1(ispec,0,0,0)) = strwzr(OPS_ACC4(0,0,0))*bclyzr(OPS_ACC_MD1(ispec,0,0,0))

!   UPDATE L2Z
    bcl2zr(OPS_ACC2(0,0,0)) = bcl2zr(OPS_ACC2(0,0,0)) + (ratezr(OPS_ACC_MD3(ispec,0,0,0))  &
              - strdzr(OPS_ACC5(0,0,0))*bclyzr(OPS_ACC_MD1(ispec,0,0,0))) *rgspec/strrzr(OPS_ACC6(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_LYZ_zr

SUBROUTINE bounds_kernel_outflowBC1_addsource_xl(drhs,urhs,vrhs,wrhs,erhs,bcl5xl,ova2xl,acouxl,struxl,strvxl,strwxl,strexl,ovgmxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl5xl,ova2xl,acouxl,struxl,strvxl,strwxl,strexl,ovgmxl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl5xl(OPS_ACC6(0,0,0))*ova2xl(OPS_ACC7(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0))  &
            - bcl5xl(OPS_ACC6(0,0,0))*ova2xl(OPS_ACC7(0,0,0))*(struxl(OPS_ACC9(0,0,0))+acouxl(OPS_ACC8(0,0,0)))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0))  &
            - bcl5xl(OPS_ACC6(0,0,0))*ova2xl(OPS_ACC7(0,0,0))*strvxl(OPS_ACC10(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0))  &
            - bcl5xl(OPS_ACC6(0,0,0))*ova2xl(OPS_ACC7(0,0,0))*strwxl(OPS_ACC11(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0))  &
            - bcl5xl(OPS_ACC6(0,0,0))*(ova2xl(OPS_ACC7(0,0,0))*strexl(OPS_ACC12(0,0,0))  &
            + struxl(OPS_ACC9(0,0,0))/acouxl(OPS_ACC8(0,0,0)) + ovgmxl(OPS_ACC13(0,0,0)))

END SUBROUTINE bounds_kernel_outflowBC1_addsource_xl

SUBROUTINE bounds_kernel_outflowBC1_addsource_yl(drhs,urhs,vrhs,wrhs,erhs,bcl5yl,ova2yl,acouyl,struyl,strvyl,strwyl,streyl,ovgmyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl5yl,ova2yl,acouyl,struyl,strvyl,strwyl,streyl,ovgmyl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl5yl(OPS_ACC6(0,0,0))*ova2yl(OPS_ACC7(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0))  &
            - bcl5yl(OPS_ACC6(0,0,0))*ova2yl(OPS_ACC7(0,0,0))*struyl(OPS_ACC9(0,0,0))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0))  &
            - bcl5yl(OPS_ACC6(0,0,0))*ova2yl(OPS_ACC7(0,0,0))*(strvyl(OPS_ACC10(0,0,0))+acouyl(OPS_ACC8(0,0,0)))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0))  &
            - bcl5yl(OPS_ACC6(0,0,0))*ova2yl(OPS_ACC7(0,0,0))*strwyl(OPS_ACC11(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0))  &
            - bcl5yl(OPS_ACC6(0,0,0))*(ova2yl(OPS_ACC7(0,0,0))*streyl(OPS_ACC12(0,0,0))  &
            + strvyl(OPS_ACC9(0,0,0))/acouyl(OPS_ACC8(0,0,0)) + ovgmyl(OPS_ACC13(0,0,0)))

END SUBROUTINE bounds_kernel_outflowBC1_addsource_yl

SUBROUTINE bounds_kernel_outflowBC1_addsource_zl(drhs,urhs,vrhs,wrhs,erhs,bcl5zl,ova2zl,acouzl,struzl,strvzl,strwzl,strezl,ovgmzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl5zl,ova2zl,acouzl,struzl,strvzl,strwzl,strezl,ovgmzl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl5zl(OPS_ACC6(0,0,0))*ova2zl(OPS_ACC7(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0))  &
            - bcl5zl(OPS_ACC6(0,0,0))*ova2zl(OPS_ACC7(0,0,0))*struzl(OPS_ACC9(0,0,0))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0))  &
            - bcl5zl(OPS_ACC6(0,0,0))*ova2zl(OPS_ACC7(0,0,0))*strvzl(OPS_ACC10(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0))  &
            - bcl5zl(OPS_ACC6(0,0,0))*ova2zl(OPS_ACC7(0,0,0))*(strwzl(OPS_ACC11(0,0,0))+acouzl(OPS_ACC8(0,0,0)))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0))  &
            - bcl5zl(OPS_ACC6(0,0,0))*(ova2zl(OPS_ACC7(0,0,0))*strezl(OPS_ACC12(0,0,0))  &
            + strwzl(OPS_ACC9(0,0,0))/acouzl(OPS_ACC8(0,0,0)) + ovgmzl(OPS_ACC13(0,0,0)))

END SUBROUTINE bounds_kernel_outflowBC1_addsource_zl

SUBROUTINE bounds_kernel_outflowBC1_eval_xl(yrhs,stryxl,bcl5xl,ova2xl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryxl,bcl5xl,ova2xl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl5xl(OPS_ACC3(0,0,0))*ova2xl(OPS_ACC4(0,0,0))*stryxl(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_eval_xl

SUBROUTINE bounds_kernel_outflowBC1_eval_yl(yrhs,stryyl,bcl5yl,ova2yl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryyl,bcl5yl,ova2yl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl5yl(OPS_ACC3(0,0,0))*ova2yl(OPS_ACC4(0,0,0))*stryyl(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_eval_yl

SUBROUTINE bounds_kernel_outflowBC1_eval_zl(yrhs,stryzl,bcl5zl,ova2zl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryzl,bcl5zl,ova2zl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl5zl(OPS_ACC3(0,0,0))*ova2zl(OPS_ACC4(0,0,0))*stryzl(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_eval_zl

SUBROUTINE bounds_kernel_inflowBC1_addsource_xl(drhs,urhs,vrhs,wrhs,erhs,bcl2xl,bcl3xl,bcl4xl,bcl5xl,struxl,strvxl,strwxl,strdxl,strexl,acouxl,ova2xl,ovgmxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl2xl,bcl3xl,bcl4xl,bcl5xl,struxl,strvxl,strwxl,strdxl,strexl,acouxl,ova2xl,ovgmxl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))  &
            - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC16(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))*struxl(OPS_ACC10(0,0,0))  &
            - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC16(0,0,0))*(struxl(OPS_ACC10(0,0,0))+acouxl(OPS_ACC15(0,0,0)))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))*strvxl(OPS_ACC11(0,0,0))  &
            - bcl3xl(OPS_ACC7(0,0,0))*strdxl(OPS_ACC13(0,0,0))  &
            - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC16(0,0,0))*strvxl(OPS_ACC11(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))*strwxl(OPS_ACC12(0,0,0))  &
            - bcl4xl(OPS_ACC8(0,0,0))*strdxl(OPS_ACC13(0,0,0))  &
            - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC16(0,0,0))*strwxl(OPS_ACC12(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))*strexl(OPS_ACC14(0,0,0))  &
            - bcl3xl(OPS_ACC7(0,0,0))*strdxl(OPS_ACC13(0,0,0))*strvxl(OPS_ACC11(0,0,0))  &
            - bcl4xl(OPS_ACC8(0,0,0))*strdxl(OPS_ACC13(0,0,0))*strwxl(OPS_ACC12(0,0,0))  &
            - bcl5xl(OPS_ACC9(0,0,0))*(ova2xl(OPS_ACC16(0,0,0))*strexl(OPS_ACC14(0,0,0))  &
            + struxl(OPS_ACC10(0,0,0))/acouxl(OPS_ACC15(0,0,0)) + ovgmxl(OPS_ACC17(0,0,0)))

END SUBROUTINE bounds_kernel_inflowBC1_addsource_xl

SUBROUTINE bounds_kernel_inflowBC1_addsource_yl(drhs,urhs,vrhs,wrhs,erhs,bcl2yl,bcl3yl,bcl4yl,bcl5yl,struyl,strvyl,strwyl,strdyl,streyl,acouyl,ova2yl,ovgmyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl2yl,bcl3yl,bcl4yl,bcl5yl,struyl,strvyl,strwyl,strdyl,streyl,acouyl,ova2yl,ovgmyl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))  &
            - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC16(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))*struyl(OPS_ACC10(0,0,0))  &
            - bcl3yl(OPS_ACC7(0,0,0))*strdyl(OPS_ACC13(0,0,0))  &
            - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC16(0,0,0))*struyl(OPS_ACC10(0,0,0))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))*strvyl(OPS_ACC11(0,0,0))  &
            - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC16(0,0,0))*(strvyl(OPS_ACC11(0,0,0))+acouyl(OPS_ACC15(0,0,0)))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))*strwyl(OPS_ACC12(0,0,0))  &
            - bcl4yl(OPS_ACC8(0,0,0))*strdyl(OPS_ACC13(0,0,0))  &
            - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC16(0,0,0))*strwyl(OPS_ACC12(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))*streyl(OPS_ACC14(0,0,0))  &
            - bcl3yl(OPS_ACC7(0,0,0))*strdyl(OPS_ACC13(0,0,0))*struyl(OPS_ACC11(0,0,0))  &
            - bcl4yl(OPS_ACC8(0,0,0))*strdyl(OPS_ACC13(0,0,0))*strwyl(OPS_ACC12(0,0,0))  &
            - bcl5yl(OPS_ACC9(0,0,0))*(ova2yl(OPS_ACC16(0,0,0))*streyl(OPS_ACC14(0,0,0))  &
            + strvyl(OPS_ACC10(0,0,0))/acouyl(OPS_ACC15(0,0,0)) + ovgmyl(OPS_ACC17(0,0,0)))

END SUBROUTINE bounds_kernel_inflowBC1_addsource_yl

SUBROUTINE bounds_kernel_inflowBC1_addsource_zl(drhs,urhs,vrhs,wrhs,erhs,bcl2zl,bcl3zl,bcl4zl,bcl5zl,struzl,strvzl,strwzl,strdzl,strezl,acouzl,ova2zl,ovgmzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl2zl,bcl3zl,bcl4zl,bcl5zl,struzl,strvzl,strwzl,strdzl,strezl,acouzl,ova2zl,ovgmzl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))  &
            - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC16(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))*struzl(OPS_ACC10(0,0,0))  &
            - bcl3zl(OPS_ACC7(0,0,0))*strdzl(OPS_ACC13(0,0,0))  &
            - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC16(0,0,0))*struzl(OPS_ACC10(0,0,0))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))*strvzl(OPS_ACC11(0,0,0))  &
            - bcl4zl(OPS_ACC8(0,0,0))*strdzl(OPS_ACC13(0,0,0))  &
            - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC16(0,0,0))*strvzl(OPS_ACC11(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))*strwzl(OPS_ACC12(0,0,0))  &
            - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC16(0,0,0))*(strwzl(OPS_ACC12(0,0,0))+acouzl(OPS_ACC15(0,0,0)))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))*strezl(OPS_ACC14(0,0,0))  &
            - bcl3zl(OPS_ACC7(0,0,0))*strdzl(OPS_ACC13(0,0,0))*struzl(OPS_ACC11(0,0,0))  &
            - bcl4zl(OPS_ACC8(0,0,0))*strdzl(OPS_ACC13(0,0,0))*strvzl(OPS_ACC12(0,0,0))  &
            - bcl5zl(OPS_ACC9(0,0,0))*(ova2zl(OPS_ACC16(0,0,0))*strezl(OPS_ACC14(0,0,0))  &
            + strwzl(OPS_ACC10(0,0,0))/acouzl(OPS_ACC15(0,0,0)) + ovgmzl(OPS_ACC17(0,0,0)))

END SUBROUTINE bounds_kernel_inflowBC1_addsource_zl

SUBROUTINE bounds_kernel_inflowBC1_eval_xl(erhs,yrhs,bclyxl,strhxl,stryxl,strdxl,bcl2xl,bcl5xl,ova2xl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,yrhs
    real(kind=dp), dimension(1), intent(in) :: bclyxl,strhxl,stryxl,strdxl,bcl2xl,bcl5xl,ova2xl
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = bclyxl(OPS_ACC_MD3(ispec,0,0,0))*strdxl(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhxl(OPS_ACC_MD4(ispec,0,0,0))

    yrhs(OPS_ACC_MD2(ispec,0,0,0)) = yrhs(OPS_ACC_MD2(ispec,0,0,0))  &
              - (bcl2xl(OPS_ACC7(0,0,0))+bcl5xl(OPS_ACC8(0,0,0))*ova2xl(OPS_ACC9(0,0,0)))*stryxl(OPS_ACC_MD5(ispec,0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_inflowBC1_eval_xl

SUBROUTINE bounds_kernel_inflowBC1_eval_yl(erhs,yrhs,bclyyl,strhyl,stryyl,strdyl,bcl2yl,bcl5yl,ova2yl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,yrhs
    real(kind=dp), dimension(1), intent(in) :: bclyyl,strhyl,stryyl,strdyl,bcl2yl,bcl5yl,ova2yl
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = bclyyl(OPS_ACC_MD3(ispec,0,0,0))*strdyl(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhyl(OPS_ACC_MD4(ispec,0,0,0))

    yrhs(OPS_ACC_MD2(ispec,0,0,0)) = yrhs(OPS_ACC_MD2(ispec,0,0,0))  &
              - (bcl2yl(OPS_ACC7(0,0,0))+bcl5yl(OPS_ACC8(0,0,0))*ova2yl(OPS_ACC9(0,0,0)))*stryyl(OPS_ACC_MD5(ispec,0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_inflowBC1_eval_yl

SUBROUTINE bounds_kernel_inflowBC1_eval_zl(erhs,yrhs,bclyzl,strhzl,stryzl,strdzl,bcl2zl,bcl5zl,ova2zl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,yrhs
    real(kind=dp), dimension(1), intent(in) :: bclyzl,strhzl,stryzl,strdzl,bcl2zl,bcl5zl,ova2zl
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = bclyzl(OPS_ACC_MD3(ispec,0,0,0))*strdzl(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhzl(OPS_ACC_MD4(ispec,0,0,0))

    yrhs(OPS_ACC_MD2(ispec,0,0,0)) = yrhs(OPS_ACC_MD2(ispec,0,0,0))  &
              - (bcl2zl(OPS_ACC7(0,0,0))+bcl5zl(OPS_ACC8(0,0,0))*ova2zl(OPS_ACC9(0,0,0)))*stryzl(OPS_ACC_MD5(ispec,0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_inflowBC1_eval_zl

SUBROUTINE bounds_kernel_inflowBC2_addsource_xl(drhs,bcl2xl,bcl5xl,ova2xl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl2xl,bcl5xl,ova2xl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2xl(OPS_ACC2(0,0,0))  &
            - bcl5xl(OPS_ACC3(0,0,0))*ova2xl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_addsource_xl

SUBROUTINE bounds_kernel_inflowBC2_addsource_yl(drhs,bcl2yl,bcl5yl,ova2yl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl2yl,bcl5yl,ova2yl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2yl(OPS_ACC2(0,0,0))  &
            - bcl5yl(OPS_ACC3(0,0,0))*ova2yl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_addsource_yl

SUBROUTINE bounds_kernel_inflowBC2_addsource_zl(drhs,bcl2zl,bcl5zl,ova2zl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl2zl,bcl5zl,ova2zl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2zl(OPS_ACC2(0,0,0))  &
            - bcl5zl(OPS_ACC3(0,0,0))*ova2zl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_addsource_zl

SUBROUTINE bounds_kernel_inflowBC3_addsource_xl(erhs,bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,strexl,struxl,strvxl,strwxl,ova2xl,ovgmxl,acouxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,strexl,struxl,strvxl,strwxl,ova2xl,ovgmxl,acouxl

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - bcl2xl(OPS_ACC2(0,0,0))*strexl(OPS_ACC7(0,0,0))  &
            - bcl3xl(OPS_ACC3(0,0,0))*strdxl(OPS_ACC6(0,0,0))*strvxl(OPS_ACC9(0,0,0))  &
            - bcl4xl(OPS_ACC4(0,0,0))*strdxl(OPS_ACC6(0,0,0))*strwxl(OPS_ACC10(0,0,0))  &
            - bcl5xl(OPS_ACC5(0,0,0))*(ova2xl(OPS_ACC11(0,0,0))*strexl(OPS_ACC7(0,0,0))  &
            + struxl(OPS_ACC8(0,0,0))/acouxl(OPS_ACC13(0,0,0)) + ovgmxl(OPS_ACC12(0,0,0)))

END SUBROUTINE bounds_kernel_inflowBC3_addsource_xl

SUBROUTINE bounds_kernel_inflowBC3_addsource_yl(erhs,bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,streyl,struyl,strvyl,strwyl,ova2yl,ovgmyl,acouyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,streyl,struyl,strvyl,strwyl,ova2yl,ovgmyl,acouyl

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - bcl2yl(OPS_ACC2(0,0,0))*streyl(OPS_ACC7(0,0,0))  &
            - bcl3yl(OPS_ACC3(0,0,0))*strdyl(OPS_ACC6(0,0,0))*struyl(OPS_ACC9(0,0,0))  &
            - bcl4yl(OPS_ACC4(0,0,0))*strdyl(OPS_ACC6(0,0,0))*strwyl(OPS_ACC10(0,0,0))  &
            - bcl5yl(OPS_ACC5(0,0,0))*(ova2yl(OPS_ACC11(0,0,0))*streyl(OPS_ACC7(0,0,0))  &
            + strvyl(OPS_ACC8(0,0,0))/acouyl(OPS_ACC13(0,0,0)) + ovgmyl(OPS_ACC12(0,0,0)))

END SUBROUTINE bounds_kernel_inflowBC3_addsource_yl

SUBROUTINE bounds_kernel_inflowBC3_addsource_zl(erhs,bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,strezl,struzl,strvzl,strwzl,ova2zl,ovgmzl,acouzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,strezl,struzl,strvzl,strwzl,ova2zl,ovgmzl,acouzl

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - bcl2zl(OPS_ACC2(0,0,0))*strezl(OPS_ACC7(0,0,0))  &
            - bcl3zl(OPS_ACC3(0,0,0))*strdzl(OPS_ACC6(0,0,0))*struzl(OPS_ACC9(0,0,0))  &
            - bcl4zl(OPS_ACC4(0,0,0))*strdzl(OPS_ACC6(0,0,0))*strvzl(OPS_ACC10(0,0,0))  &
            - bcl5zl(OPS_ACC5(0,0,0))*(ova2zl(OPS_ACC11(0,0,0))*strezl(OPS_ACC7(0,0,0))  &
            + strwzl(OPS_ACC8(0,0,0))/acouzl(OPS_ACC13(0,0,0)) + ovgmzl(OPS_ACC12(0,0,0)))

END SUBROUTINE bounds_kernel_inflowBC3_addsource_zl

SUBROUTINE bounds_kernel_inflowBC3_eval_xl(erhs,bclyxl,ratexl,dydtxl,strhxl,strdxl,struxl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,bclyxl
    real(kind=dp), dimension(1), intent(in) :: ratexl,dydtxl,strhxl,strdxl,struxl
    integer, intent(in) :: ispec

    bclyxl(OPS_ACC_MD2(ispec,0,0,0)) = ratexl(OPS_ACC_MD3(ispec,0,0,0))/strdxl(OPS_ACC6(0,0,0))  &
              - dydtxl(OPS_ACC_MD4(ispec,0,0,0)) - struxl(OPS_ACC7(0,0,0))*bclyxl(OPS_ACC_MD2(ispec,0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyxl(OPS_ACC_MD2(ispec,0,0,0))*strdxl(OPS_ACC6(0,0,0))*strhxl(OPS_ACC_MD5(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_eval_xl

SUBROUTINE bounds_kernel_inflowBC3_eval_yl(erhs,bclyyl,rateyl,dydtyl,strhyl,strdyl,strvyl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,bclyyl
    real(kind=dp), dimension(1), intent(in) :: rateyl,dydtyl,strhyl,strdyl,strvyl
    integer, intent(in) :: ispec

    bclyyl(OPS_ACC_MD2(ispec,0,0,0)) = rateyl(OPS_ACC_MD3(ispec,0,0,0))/strdyl(OPS_ACC6(0,0,0))  &
              - dydtyl(OPS_ACC_MD4(ispec,0,0,0)) - strvyl(OPS_ACC7(0,0,0))*bclyyl(OPS_ACC_MD2(ispec,0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyyl(OPS_ACC_MD2(ispec,0,0,0))*strdyl(OPS_ACC6(0,0,0))*strhyl(OPS_ACC_MD5(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_eval_yl

SUBROUTINE bounds_kernel_inflowBC3_eval_zl(erhs,bclyzl,ratezl,dydtzl,strhzl,strdzl,strwzl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,bclyzl
    real(kind=dp), dimension(1), intent(in) :: ratezl,dydtzl,strhzl,strdzl,strwzl
    integer, intent(in) :: ispec

    bclyzl(OPS_ACC_MD2(ispec,0,0,0)) = ratezl(OPS_ACC_MD3(ispec,0,0,0))/strdzl(OPS_ACC6(0,0,0))  &
              - dydtzl(OPS_ACC_MD4(ispec,0,0,0)) - strwzl(OPS_ACC7(0,0,0))*bclyzl(OPS_ACC_MD2(ispec,0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyzl(OPS_ACC_MD2(ispec,0,0,0))*strdzl(OPS_ACC6(0,0,0))*strhzl(OPS_ACC_MD5(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_eval_zl

SUBROUTINE bounds_kernel_wallBC1_addsource_xl(drhs,erhs,bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,strexl,struxl,strvxl,strwxl,ova2xl,ovgmxl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,strexl,struxl,strvxl,strwxl,ova2xl,ovgmxl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl5xl(OPS_ACC5(0,0,0))*ova2xl(OPS_ACC12(0,0,0))

    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0))  &
            - bcl3xl(OPS_ACC3(0,0,0))*strdxl(OPS_ACC6(0,0,0))*strvxl(OPS_ACC10(0,0,0))  &
            - bcl4xl(OPS_ACC4(0,0,0))*strdxl(OPS_ACC6(0,0,0))*strwxl(OPS_ACC11(0,0,0))  &
            - bcl5xl(OPS_ACC5(0,0,0))*(ova2xl(OPS_ACC12(0,0,0))*strexl(OPS_ACC8(0,0,0))  &
            + struxl(OPS_ACC9(0,0,0))/acouxl(OPS_ACC7(0,0,0)) + ovgmxl(OPS_ACC13(0,0,0)))

END SUBROUTINE bounds_kernel_wallBC1_addsource_xl

SUBROUTINE bounds_kernel_wallBC1_addsource_yl(drhs,erhs,bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,streyl,struyl,strvyl,strwyl,ova2yl,ovgmyl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,streyl,struyl,strvyl,strwyl,ova2yl,ovgmyl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl5yl(OPS_ACC5(0,0,0))*ova2yl(OPS_ACC12(0,0,0))

    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0))  &
            - bcl3yl(OPS_ACC3(0,0,0))*strdyl(OPS_ACC6(0,0,0))*struyl(OPS_ACC10(0,0,0))  &
            - bcl4yl(OPS_ACC4(0,0,0))*strdyl(OPS_ACC6(0,0,0))*strwyl(OPS_ACC11(0,0,0))  &
            - bcl5yl(OPS_ACC5(0,0,0))*(ova2yl(OPS_ACC12(0,0,0))*streyl(OPS_ACC8(0,0,0))  &
            + strvyl(OPS_ACC9(0,0,0))/acouyl(OPS_ACC7(0,0,0)) + ovgmyl(OPS_ACC13(0,0,0)))

END SUBROUTINE bounds_kernel_wallBC1_addsource_yl

SUBROUTINE bounds_kernel_wallBC1_addsource_zl(drhs,erhs,bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,strezl,struzl,strvzl,strwzl,ova2zl,ovgmzl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,strezl,struzl,strvzl,strwzl,ova2zl,ovgmzl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl5zl(OPS_ACC5(0,0,0))*ova2zl(OPS_ACC12(0,0,0))

    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0))  &
            - bcl3zl(OPS_ACC3(0,0,0))*strdzl(OPS_ACC6(0,0,0))*struzl(OPS_ACC10(0,0,0))  &
            - bcl4zl(OPS_ACC4(0,0,0))*strdzl(OPS_ACC6(0,0,0))*strvzl(OPS_ACC11(0,0,0))  &
            - bcl5zl(OPS_ACC5(0,0,0))*(ova2zl(OPS_ACC12(0,0,0))*strezl(OPS_ACC8(0,0,0))  &
            + strwzl(OPS_ACC9(0,0,0))/acouzl(OPS_ACC7(0,0,0)) + ovgmzl(OPS_ACC13(0,0,0)))

END SUBROUTINE bounds_kernel_wallBC1_addsource_zl

SUBROUTINE bounds_kernel_wallBC1_eval_xl(yrhs,stryxl,bcl5xl,ova2xl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryxl,bcl5xl,ova2xl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl5xl(OPS_ACC3(0,0,0))*ova2xl(OPS_ACC4(0,0,0))*stryxl(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_eval_xl

SUBROUTINE bounds_kernel_wallBC1_eval_yl(yrhs,stryyl,bcl5yl,ova2yl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryyl,bcl5yl,ova2yl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl5yl(OPS_ACC3(0,0,0))*ova2yl(OPS_ACC4(0,0,0))*stryyl(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_eval_yl

SUBROUTINE bounds_kernel_wallBC1_eval_zl(yrhs,stryzl,bcl5zl,ova2zl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryzl,bcl5zl,ova2zl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl5zl(OPS_ACC3(0,0,0))*ova2zl(OPS_ACC4(0,0,0))*stryzl(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_eval_zl

SUBROUTINE bounds_kernel_wallBC2_addsource_xl(drhs,bcl2xl,bcl5xl,ova2xl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl2xl,bcl5xl,ova2xl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2xl(OPS_ACC2(0,0,0))  &
            - bcl5xl(OPS_ACC3(0,0,0))*ova2xl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_addsource_xl

SUBROUTINE bounds_kernel_wallBC2_addsource_yl(drhs,bcl2yl,bcl5yl,ova2yl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl2yl,bcl5yl,ova2yl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2yl(OPS_ACC2(0,0,0))  &
            - bcl5yl(OPS_ACC3(0,0,0))*ova2yl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_addsource_yl

SUBROUTINE bounds_kernel_wallBC2_addsource_zl(drhs,bcl2zl,bcl5zl,ova2zl)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl2zl,bcl5zl,ova2zl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2zl(OPS_ACC2(0,0,0))  &
            - bcl5zl(OPS_ACC3(0,0,0))*ova2zl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_addsource_zl

SUBROUTINE bounds_kernel_wallBC2_eval_xl(yrhs,stryxl,bcl2xl,bcl5xl,ova2xl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryxl,bcl2xl,bcl5xl,ova2xl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - (bcl2xl(OPS_ACC3(0,0,0))+bcl5xl(OPS_ACC4(0,0,0))*ova2xl(OPS_ACC5(0,0,0)))*stryxl(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_eval_xl

SUBROUTINE bounds_kernel_wallBC2_eval_yl(yrhs,stryyl,bcl2yl,bcl5yl,ova2yl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryyl,bcl2yl,bcl5yl,ova2yl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - (bcl2yl(OPS_ACC3(0,0,0))+bcl5yl(OPS_ACC4(0,0,0))*ova2yl(OPS_ACC5(0,0,0)))*stryyl(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_eval_yl

SUBROUTINE bounds_kernel_wallBC2_eval_zl(yrhs,stryzl,bcl2zl,bcl5zl,ova2zl,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryzl,bcl2zl,bcl5zl,ova2zl
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - (bcl2zl(OPS_ACC3(0,0,0))+bcl5zl(OPS_ACC4(0,0,0))*ova2zl(OPS_ACC5(0,0,0)))*stryzl(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_eval_zl

SUBROUTINE bounds_kernel_outflowBC1_addsource_xr(drhs,urhs,vrhs,wrhs,erhs,bcl1xr,ova2xr,acouxr,struxr,strvxr,strwxr,strexr,ovgmxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1xr,ova2xr,acouxr,struxr,strvxr,strwxr,strexr,ovgmxr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC7(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0))  &
            - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC7(0,0,0))*(struxr(OPS_ACC9(0,0,0))-acouxr(OPS_ACC8(0,0,0)))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0))  &
            - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC7(0,0,0))*strvxr(OPS_ACC10(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0))  &
            - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC7(0,0,0))*strwxr(OPS_ACC11(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0))  &
            - bcl1xr(OPS_ACC6(0,0,0))*(ova2xr(OPS_ACC7(0,0,0))*strexr(OPS_ACC12(0,0,0))  &
            - struxr(OPS_ACC9(0,0,0))/acouxr(OPS_ACC8(0,0,0)) + ovgmxr(OPS_ACC13(0,0,0)))

END SUBROUTINE bounds_kernel_outflowBC1_addsource_xr

SUBROUTINE bounds_kernel_outflowBC1_addsource_yr(drhs,urhs,vrhs,wrhs,erhs,bcl1yr,ova2yr,acouyr,struyr,strvyr,strwyr,streyr,ovgmyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1yr,ova2yr,acouyr,struyr,strvyr,strwyr,streyr,ovgmyr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC7(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0))  &
            - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC7(0,0,0))*struyr(OPS_ACC9(0,0,0))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0))  &
            - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC7(0,0,0))*(strvyr(OPS_ACC10(0,0,0))-acouyr(OPS_ACC8(0,0,0)))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0))  &
            - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC7(0,0,0))*strwyr(OPS_ACC11(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0))  &
            - bcl1yr(OPS_ACC6(0,0,0))*(ova2yr(OPS_ACC7(0,0,0))*streyr(OPS_ACC12(0,0,0))  &
            - strvyr(OPS_ACC9(0,0,0))/acouyr(OPS_ACC8(0,0,0)) + ovgmyr(OPS_ACC13(0,0,0)))

END SUBROUTINE bounds_kernel_outflowBC1_addsource_yr

SUBROUTINE bounds_kernel_outflowBC1_addsource_zr(drhs,urhs,vrhs,wrhs,erhs,bcl1zr,ova2zr,acouzr,struzr,strvzr,strwzr,strezr,ovgmzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1zr,ova2zr,acouzr,struzr,strvzr,strwzr,strezr,ovgmzr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC7(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0))  &
            - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC7(0,0,0))*struzr(OPS_ACC9(0,0,0))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0))  &
            - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC7(0,0,0))*strvzr(OPS_ACC10(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0))  &
            - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC7(0,0,0))*(strwzr(OPS_ACC11(0,0,0))-acouzr(OPS_ACC8(0,0,0)))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0))  &
            - bcl1zr(OPS_ACC6(0,0,0))*(ova2zr(OPS_ACC7(0,0,0))*strezr(OPS_ACC12(0,0,0))  &
            - strwzr(OPS_ACC9(0,0,0))/acouzr(OPS_ACC8(0,0,0)) + ovgmzr(OPS_ACC13(0,0,0)))

END SUBROUTINE bounds_kernel_outflowBC1_addsource_zr

SUBROUTINE bounds_kernel_outflowBC1_eval_xr(yrhs,stryxr,bcl1xr,ova2xr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryxr,bcl1xr,ova2xr
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl1xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC4(0,0,0))*stryxr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_eval_xr

SUBROUTINE bounds_kernel_outflowBC1_eval_yr(yrhs,stryyr,bcl1yr,ova2yr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryyr,bcl1yr,ova2yr
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl1yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC4(0,0,0))*stryyr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_eval_yr

SUBROUTINE bounds_kernel_outflowBC1_eval_zr(yrhs,stryzr,bcl1zr,ova2zr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryzr,bcl1zr,ova2zr
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl1zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC4(0,0,0))*stryzr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_outflowBC1_eval_zr

SUBROUTINE bounds_kernel_inflowBC1_addsource_xr(drhs,urhs,vrhs,wrhs,erhs,bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,strexr,acouxr,struxr,strvxr,strwxr,ova2xr,ovgmxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,strexr,acouxr,struxr,strvxr,strwxr,ova2xr,ovgmxr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC16(0,0,0))  &
            - bcl2xr(OPS_ACC7(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0))  &
            - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC16(0,0,0))*(struxr(OPS_ACC13(0,0,0))-acouxr(OPS_ACC12(0,0,0)))  &
            - bcl2xr(OPS_ACC7(0,0,0))*struxr(OPS_ACC13(0,0,0))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0))  &
            - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC16(0,0,0))*strvxr(OPS_ACC14(0,0,0))  &
            - bcl2xr(OPS_ACC7(0,0,0))*strvxr(OPS_ACC14(0,0,0)) - bcl3xr(OPS_ACC8(0,0,0))*strdxr(OPS_ACC10(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0))  &
            - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC16(0,0,0))*strwxr(OPS_ACC15(0,0,0))  &
            - bcl2xr(OPS_ACC7(0,0,0))*strwxr(OPS_ACC15(0,0,0)) - bcl4xr(OPS_ACC9(0,0,0))*strdxr(OPS_ACC10(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0))  &
            - bcl1xr(OPS_ACC6(0,0,0))*(ova2xr(OPS_ACC16(0,0,0))*strexr(OPS_ACC11(0,0,0))  &
            - struxr(OPS_ACC13(0,0,0))/acouxr(OPS_ACC12(0,0,0)) + ovgmxr(OPS_ACC17(0,0,0)))  &
            - bcl2xr(OPS_ACC7(0,0,0))*strexr(OPS_ACC11(0,0,0))  &
            - bcl3xr(OPS_ACC8(0,0,0))*strdxr(OPS_ACC10(0,0,0))*strvxr(OPS_ACC14(0,0,0))  &
            - bcl4xr(OPS_ACC9(0,0,0))*strdxr(OPS_ACC10(0,0,0))*strwxr(OPS_ACC15(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_addsource_xr

SUBROUTINE bounds_kernel_inflowBC1_addsource_yr(drhs,urhs,vrhs,wrhs,erhs,bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,streyr,acouyr,struyr,strvyr,strwyr,ova2yr,ovgmyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,streyr,acouyr,struyr,strvyr,strwyr,ova2yr,ovgmyr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC16(0,0,0))  &
            - bcl2yr(OPS_ACC7(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0))  &
            - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC16(0,0,0))*struyr(OPS_ACC13(0,0,0))  &
            - bcl2yr(OPS_ACC7(0,0,0))*struyr(OPS_ACC13(0,0,0)) - bcl3yr(OPS_ACC8(0,0,0))*strdyr(OPS_ACC10(0,0,0))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0))  &
            - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC16(0,0,0))*(strvyr(OPS_ACC14(0,0,0))-acouyr(OPS_ACC12(0,0,0)))  &
            - bcl2yr(OPS_ACC7(0,0,0))*strvyr(OPS_ACC14(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0))  &
            - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC16(0,0,0))*strwyr(OPS_ACC15(0,0,0))  &
            - bcl2yr(OPS_ACC7(0,0,0))*strwyr(OPS_ACC15(0,0,0)) - bcl4yr(OPS_ACC9(0,0,0))*strdyr(OPS_ACC10(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0))  &
            - bcl1yr(OPS_ACC6(0,0,0))*(ova2yr(OPS_ACC16(0,0,0))*streyr(OPS_ACC11(0,0,0))  &
            + strvyr(OPS_ACC13(0,0,0))/acouyr(OPS_ACC12(0,0,0)) + ovgmyr(OPS_ACC17(0,0,0)))  &
            - bcl2yr(OPS_ACC7(0,0,0))*streyr(OPS_ACC11(0,0,0))  &
            - bcl3yr(OPS_ACC8(0,0,0))*strdyr(OPS_ACC10(0,0,0))*struyr(OPS_ACC14(0,0,0))  &
            - bcl4yr(OPS_ACC9(0,0,0))*strdyr(OPS_ACC10(0,0,0))*strwyr(OPS_ACC15(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_addsource_yr

SUBROUTINE bounds_kernel_inflowBC1_addsource_zr(drhs,urhs,vrhs,wrhs,erhs,bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,strezr,acouzr,struzr,strvzr,strwzr,ova2zr,ovgmzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,strezr,acouzr,struzr,strvzr,strwzr,ova2zr,ovgmzr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC16(0,0,0))  &
            - bcl2zr(OPS_ACC7(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0))  &
            - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC16(0,0,0))*struzr(OPS_ACC13(0,0,0))  &
            - bcl2zr(OPS_ACC7(0,0,0))*struzr(OPS_ACC13(0,0,0)) - bcl3zr(OPS_ACC8(0,0,0))*strdzr(OPS_ACC10(0,0,0))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0))  &
            - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC16(0,0,0))*strvzr(OPS_ACC14(0,0,0))  &
            - bcl2zr(OPS_ACC7(0,0,0))*strvzr(OPS_ACC14(0,0,0)) - bcl4zr(OPS_ACC9(0,0,0))*strdzr(OPS_ACC10(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0))  &
            - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC16(0,0,0))*(strwzr(OPS_ACC15(0,0,0))-acouzr(OPS_ACC12(0,0,0)))  &
            - bcl2zr(OPS_ACC7(0,0,0))*strwzr(OPS_ACC15(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0))  &
            - bcl1zr(OPS_ACC6(0,0,0))*(ova2zr(OPS_ACC16(0,0,0))*strezr(OPS_ACC11(0,0,0))  &
            + strwzr(OPS_ACC13(0,0,0))/acouzr(OPS_ACC12(0,0,0)) + ovgmzr(OPS_ACC17(0,0,0)))  &
            - bcl2zr(OPS_ACC7(0,0,0))*strezr(OPS_ACC11(0,0,0))  &
            - bcl3zr(OPS_ACC8(0,0,0))*strdzr(OPS_ACC10(0,0,0))*struzr(OPS_ACC14(0,0,0))  &
            - bcl4zr(OPS_ACC9(0,0,0))*strdzr(OPS_ACC10(0,0,0))*strvzr(OPS_ACC15(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC1_addsource_zr

SUBROUTINE bounds_kernel_inflowBC1_eval_xr(erhs,yrhs,bclyxr,strhxr,stryxr,strdxr,bcl2xr,bcl1xr,ova2xr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,yrhs
    real(kind=dp), dimension(1), intent(in) :: bclyxr,strhxr,stryxr,strdxr,bcl2xr,bcl1xr,ova2xr
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = bclyxr(OPS_ACC_MD3(ispec,0,0,0))*strdxr(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhxr(OPS_ACC_MD4(ispec,0,0,0))

    yrhs(OPS_ACC_MD2(ispec,0,0,0)) = yrhs(OPS_ACC_MD2(ispec,0,0,0))  &
              - (bcl2xr(OPS_ACC7(0,0,0))+bcl1xr(OPS_ACC8(0,0,0))*ova2xr(OPS_ACC9(0,0,0)))*stryxr(OPS_ACC_MD5(ispec,0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_inflowBC1_eval_xr

SUBROUTINE bounds_kernel_inflowBC1_eval_yr(erhs,yrhs,bclyyr,strhyr,stryyr,strdyr,bcl2yr,bcl1yr,ova2yr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,yrhs
    real(kind=dp), dimension(1), intent(in) :: bclyyr,strhyr,stryyr,strdyr,bcl2yr,bcl1yr,ova2yr
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = bclyyr(OPS_ACC_MD3(ispec,0,0,0))*strdyr(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhyr(OPS_ACC_MD4(ispec,0,0,0))

    yrhs(OPS_ACC_MD2(ispec,0,0,0)) = yrhs(OPS_ACC_MD2(ispec,0,0,0))  &
              - (bcl2yr(OPS_ACC7(0,0,0))+bcl1yr(OPS_ACC8(0,0,0))*ova2yr(OPS_ACC9(0,0,0)))*stryyr(OPS_ACC_MD5(ispec,0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_inflowBC1_eval_yr

SUBROUTINE bounds_kernel_inflowBC1_eval_zr(erhs,yrhs,bclyzr,strhzr,stryzr,strdzr,bcl2zr,bcl5zr,ova2zr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,yrhs
    real(kind=dp), dimension(1), intent(in) :: bclyzr,strhzr,stryzr,strdzr,bcl2zr,bcl5zr,ova2zr
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = bclyzr(OPS_ACC_MD3(ispec,0,0,0))*strdzr(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhzr(OPS_ACC_MD4(ispec,0,0,0))

    yrhs(OPS_ACC_MD2(ispec,0,0,0)) = yrhs(OPS_ACC_MD2(ispec,0,0,0))  &
              - (bcl2zr(OPS_ACC7(0,0,0))+bcl5zr(OPS_ACC8(0,0,0))*ova2zr(OPS_ACC9(0,0,0)))*stryzr(OPS_ACC_MD5(ispec,0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_inflowBC1_eval_zr

SUBROUTINE bounds_kernel_inflowBC2_addsource_xr(drhs,bcl2xr,bcl1xr,ova2xr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl2xr,bcl1xr,ova2xr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2xr(OPS_ACC2(0,0,0))  &
            - bcl1xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_addsource_xr

SUBROUTINE bounds_kernel_inflowBC2_addsource_yr(drhs,bcl2yr,bcl1yr,ova2yr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl2yr,bcl1yr,ova2yr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2yr(OPS_ACC2(0,0,0))  &
            - bcl1yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_addsource_yr

SUBROUTINE bounds_kernel_inflowBC2_addsource_zr(drhs,bcl2zr,bcl1zr,ova2zr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl2zr,bcl1zr,ova2zr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2zr(OPS_ACC2(0,0,0))  &
            - bcl1zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC2_addsource_zr

SUBROUTINE bounds_kernel_inflowBC3_addsource_xr(erhs,bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,strexr,struxr,strvxr,strwxr,ova2xr,ovgmxr,acouxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,strexr,struxr,strvxr,strwxr,ova2xr,ovgmxr,acouxr

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
            - bcl1xr(OPS_ACC2(0,0,0))*(ova2xr(OPS_ACC11(0,0,0))*strexr(OPS_ACC7(0,0,0))  &
            - struxr(OPS_ACC8(0,0,0))/acouxr(OPS_ACC13(0,0,0)) + ovgmxr(OPS_ACC12(0,0,0)))  &
            - bcl2xr(OPS_ACC3(0,0,0))*strexr(OPS_ACC7(0,0,0))  &
            - bcl3xr(OPS_ACC4(0,0,0))*strdxr(OPS_ACC6(0,0,0))*strvxr(OPS_ACC9(0,0,0))  &
            - bcl4xr(OPS_ACC5(0,0,0))*strdxr(OPS_ACC6(0,0,0))*strwxr(OPS_ACC10(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_addsource_xr

SUBROUTINE bounds_kernel_inflowBC3_addsource_yr(erhs,bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,streyr,struyr,strvyr,strwyr,ova2yr,ovgmyr,acouyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,streyr,struyr,strvyr,strwyr,ova2yr,ovgmyr,acouyr

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
            - bcl1yr(OPS_ACC2(0,0,0))*(ova2yr(OPS_ACC11(0,0,0))*streyr(OPS_ACC7(0,0,0))  &
            + strvyr(OPS_ACC8(0,0,0))/acouyr(OPS_ACC13(0,0,0)) + ovgmyr(OPS_ACC12(0,0,0)))  &
            - bcl2yr(OPS_ACC3(0,0,0))*streyr(OPS_ACC7(0,0,0))  &
            - bcl3yr(OPS_ACC4(0,0,0))*strdyr(OPS_ACC6(0,0,0))*struyr(OPS_ACC9(0,0,0))  &
            - bcl4yr(OPS_ACC5(0,0,0))*strdyr(OPS_ACC6(0,0,0))*strwyr(OPS_ACC10(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_addsource_yr

SUBROUTINE bounds_kernel_inflowBC3_addsource_zr(erhs,bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,strezr,struzr,strvzr,strwzr,ova2zr,ovgmzr,acouzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,strezr,struzr,strvzr,strwzr,ova2zr,ovgmzr,acouzr

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
            - bcl1zr(OPS_ACC2(0,0,0))*(ova2zr(OPS_ACC11(0,0,0))*strezr(OPS_ACC7(0,0,0))  &
            + strwzr(OPS_ACC8(0,0,0))/acouzr(OPS_ACC13(0,0,0)) + ovgmzr(OPS_ACC12(0,0,0)))  &
            - bcl2zr(OPS_ACC3(0,0,0))*strezr(OPS_ACC7(0,0,0))  &
            - bcl3zr(OPS_ACC4(0,0,0))*strdzr(OPS_ACC6(0,0,0))*struzr(OPS_ACC9(0,0,0))  &
            - bcl4zr(OPS_ACC5(0,0,0))*strdzr(OPS_ACC6(0,0,0))*strvzr(OPS_ACC10(0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_addsource_zr

SUBROUTINE bounds_kernel_inflowBC3_eval_xr(erhs,bclyxr,ratexr,dydtxr,strhxr,strdxr,struxr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,bclyxr
    real(kind=dp), dimension(1), intent(in) :: ratexr,dydtxr,strhxr,strdxr,struxr
    integer, intent(in) :: ispec

    bclyxr(OPS_ACC_MD2(ispec,0,0,0)) = ratexr(OPS_ACC_MD3(ispec,0,0,0))/strdxr(OPS_ACC6(0,0,0))  &
              - dydtxr(OPS_ACC_MD4(ispec,0,0,0)) - struxr(OPS_ACC7(0,0,0))*bclyxr(OPS_ACC_MD2(ispec,0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyxr(OPS_ACC_MD2(ispec,0,0,0))*strdxr(OPS_ACC6(0,0,0))*strhxr(OPS_ACC_MD5(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_eval_xr

SUBROUTINE bounds_kernel_inflowBC3_eval_yr(erhs,bclyyr,rateyr,dydtyr,strhyr,strdyr,strvyr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,bclyyr
    real(kind=dp), dimension(1), intent(in) :: rateyr,dydtyr,strhyr,strdyr,strvyr
    integer, intent(in) :: ispec

    bclyyr(OPS_ACC_MD2(ispec,0,0,0)) = rateyr(OPS_ACC_MD3(ispec,0,0,0))/strdyr(OPS_ACC6(0,0,0))  &
              - dydtyr(OPS_ACC_MD4(ispec,0,0,0)) - strvyr(OPS_ACC7(0,0,0))*bclyyr(OPS_ACC_MD2(ispec,0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyyr(OPS_ACC_MD2(ispec,0,0,0))*strdyr(OPS_ACC6(0,0,0))*strhyr(OPS_ACC_MD5(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_eval_yr

SUBROUTINE bounds_kernel_inflowBC3_eval_zr(erhs,bclyzr,ratezr,dydtzr,strhzr,strdzr,strwzr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs,bclyzr
    real(kind=dp), dimension(1), intent(in) :: ratezr,dydtzr,strhzr,strdzr,strwzr
    integer, intent(in) :: ispec

    bclyzr(OPS_ACC_MD2(ispec,0,0,0)) = ratezr(OPS_ACC_MD3(ispec,0,0,0))/strdzr(OPS_ACC6(0,0,0))  &
              - dydtzr(OPS_ACC_MD4(ispec,0,0,0)) - strwzr(OPS_ACC7(0,0,0))*bclyzr(OPS_ACC_MD2(ispec,0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyzr(OPS_ACC_MD2(ispec,0,0,0))*strdzr(OPS_ACC6(0,0,0))*strhzr(OPS_ACC_MD5(ispec,0,0,0))

END SUBROUTINE bounds_kernel_inflowBC3_eval_zr

SUBROUTINE bounds_kernel_wallBC1_addsource_xr(drhs,erhs,bcl1xr,bcl3xr,bcl4xr,strdxr,acouxr,strexr,struxr,strvxr,strwxr,ova2xr,ovgmxr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1xr,bcl3xr,bcl4xr,strdxr,acouxr,strexr,struxr,strvxr,strwxr,ova2xr,ovgmxr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC12(0,0,0))

    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0))  &
            - bcl1xr(OPS_ACC3(0,0,0))*(ova2xr(OPS_ACC12(0,0,0))*strexr(OPS_ACC8(0,0,0))  &
            + struxr(OPS_ACC9(0,0,0))/acouxr(OPS_ACC7(0,0,0)) + ovgmxr(OPS_ACC13(0,0,0)))  &
            - bcl3xr(OPS_ACC4(0,0,0))*strdxr(OPS_ACC6(0,0,0))*strvxr(OPS_ACC10(0,0,0))  &
            - bcl4xr(OPS_ACC5(0,0,0))*strdxr(OPS_ACC6(0,0,0))*strwxr(OPS_ACC11(0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_addsource_xr

SUBROUTINE bounds_kernel_wallBC1_addsource_yr(drhs,erhs,bcl1yr,bcl3yr,bcl4yr,strdyr,acouyr,streyr,struyr,strvyr,strwyr,ova2yr,ovgmyr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1yr,bcl3yr,bcl4yr,strdyr,acouyr,streyr,struyr,strvyr,strwyr,ova2yr,ovgmyr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC12(0,0,0))

    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0))  &
            - bcl1yr(OPS_ACC3(0,0,0))*(ova2yr(OPS_ACC12(0,0,0))*streyr(OPS_ACC8(0,0,0))  &
            + strvyr(OPS_ACC9(0,0,0))/acouyr(OPS_ACC7(0,0,0)) + ovgmyr(OPS_ACC13(0,0,0)))  &
            - bcl3yr(OPS_ACC4(0,0,0))*strdyr(OPS_ACC6(0,0,0))*struyr(OPS_ACC10(0,0,0))  &
            - bcl4yr(OPS_ACC5(0,0,0))*strdyr(OPS_ACC6(0,0,0))*strwyr(OPS_ACC11(0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_addsource_yr

SUBROUTINE bounds_kernel_wallBC1_addsource_zr(drhs,erhs,bcl1zr,bcl3zr,bcl4zr,strdzr,acouzr,strezr,struzr,strvzr,strwzr,ova2zr,ovgmzr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs,erhs
    real(kind=dp), dimension(1), intent(in) :: bcl1zr,bcl3zr,bcl4zr,strdzr,acouzr,strezr,struzr,strvzr,strwzr,ova2zr,ovgmzr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC12(0,0,0))

    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0))  &
            - bcl1zr(OPS_ACC3(0,0,0))*(ova2zr(OPS_ACC12(0,0,0))*strezr(OPS_ACC8(0,0,0))  &
            + strwzr(OPS_ACC9(0,0,0))/acouzr(OPS_ACC7(0,0,0)) + ovgmzr(OPS_ACC13(0,0,0)))  &
            - bcl3zr(OPS_ACC4(0,0,0))*strdzr(OPS_ACC6(0,0,0))*struzr(OPS_ACC10(0,0,0))  &
            - bcl4zr(OPS_ACC5(0,0,0))*strdzr(OPS_ACC6(0,0,0))*strvzr(OPS_ACC11(0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_addsource_zr

SUBROUTINE bounds_kernel_wallBC1_eval_xr(yrhs,stryxr,bcl1xr,ova2xr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryxr,bcl1xr,ova2xr
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl1xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC4(0,0,0))*stryxr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_eval_xr

SUBROUTINE bounds_kernel_wallBC1_eval_yr(yrhs,stryyr,bcl1yr,ova2yr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryyr,bcl1yr,ova2yr
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl1yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC4(0,0,0))*stryyr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_eval_yr

SUBROUTINE bounds_kernel_wallBC1_eval_zr(yrhs,stryzr,bcl1zr,ova2zr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryzr,bcl1zr,ova2zr
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - bcl1zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC4(0,0,0))*stryzr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC1_eval_zr

SUBROUTINE bounds_kernel_wallBC2_addsource_xr(drhs,bcl1xr,bcl2xr,ova2xr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl1xr,bcl2xr,ova2xr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1xr(OPS_ACC2(0,0,0))*ova2xr(OPS_ACC4(0,0,0))  &
                    - bcl2xr(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_addsource_xr

SUBROUTINE bounds_kernel_wallBC2_addsource_yr(drhs,bcl1yr,bcl2yr,ova2yr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl1yr,bcl2yr,ova2yr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1yr(OPS_ACC2(0,0,0))*ova2yr(OPS_ACC4(0,0,0))  &
                    - bcl2yr(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_addsource_yr

SUBROUTINE bounds_kernel_wallBC2_addsource_zr(drhs,bcl1zr,bcl2zr,ova2zr)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: drhs
    real(kind=dp), dimension(1), intent(in) :: bcl1zr,bcl2zr,ova2zr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1zr(OPS_ACC2(0,0,0))*ova2zr(OPS_ACC4(0,0,0))  &
                    - bcl2zr(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_addsource_zr

SUBROUTINE bounds_kernel_wallBC2_eval_xr(yrhs,stryxr,bcl1xr,bcl2xr,ova2xr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryxr,bcl1xr,bcl2xr,ova2xr
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - (bcl2xr(OPS_ACC4(0,0,0))+bcl1xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC5(0,0,0)))*stryxr(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_eval_xr

SUBROUTINE bounds_kernel_wallBC2_eval_yr(yrhs,stryyr,bcl1yr,bcl2yr,ova2yr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryyr,bcl1yr,bcl2yr,ova2yr
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - (bcl2yr(OPS_ACC4(0,0,0))+bcl1yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC5(0,0,0)))*stryyr(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_eval_yr

SUBROUTINE bounds_kernel_wallBC2_eval_zr(yrhs,stryzr,bcl1zr,bcl2zr,ova2zr,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: yrhs
    real(kind=dp), dimension(1), intent(in) :: stryzr,bcl1zr,bcl2zr,ova2zr
    integer, intent(in) :: ispec

    yrhs(OPS_ACC_MD1(ispec,0,0,0)) = yrhs(OPS_ACC_MD1(ispec,0,0,0))  &
              - (bcl2zr(OPS_ACC4(0,0,0))+bcl1zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC5(0,0,0)))*stryzr(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE bounds_kernel_wallBC2_eval_zr




