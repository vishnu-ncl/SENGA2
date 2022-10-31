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
