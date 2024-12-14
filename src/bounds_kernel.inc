SUBROUTINE bounds_kernel_eqA_xdir(strhx,strgx,strtx,strrx,rgspec,ispec)

    real(kind=8), dimension(1) :: strhx
    real(kind=8), dimension(1), intent(in) :: strgx,strtx,strrx
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    strhx(OPS_ACC1(0,0,0)) = strhx(OPS_ACC1(0,0,0))  &
            - strgx(OPS_ACC2(0,0,0))*strtx(OPS_ACC3(0,0,0))*rgspec(ispec)/strrx(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqA_xdir

SUBROUTINE bounds_kernel_eqA_ydir(strhy,strgy,strty,strry,rgspec,ispec)

    real(kind=8), dimension(1) :: strhy
    real(kind=8), dimension(1), intent(in) :: strgy,strty,strry
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    strhy(OPS_ACC1(0,0,0)) = strhy(OPS_ACC1(0,0,0))  &
            - strgy(OPS_ACC2(0,0,0))*strty(OPS_ACC3(0,0,0))*rgspec(ispec)/strry(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqA_ydir

SUBROUTINE bounds_kernel_eqA_zdir(strhz,strgz,strtz,strrz,rgspec,ispec)

    real(kind=8), dimension(1) :: strhz
    real(kind=8), dimension(1), intent(in) :: strgz,strtz,strrz
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    strhz(OPS_ACC1(0,0,0)) = strhz(OPS_ACC1(0,0,0))  &
            - strgz(OPS_ACC2(0,0,0))*strtz(OPS_ACC3(0,0,0))*rgspec(ispec)/strrz(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqA_zdir

SUBROUTINE bounds_kernel_eqB_xdir(gam1x,strex,ovgmx,strgx,strrx,strtx)

    real(kind=8), dimension(1) :: gam1x,strex,ovgmx
    real(kind=8), dimension(1), intent(in) :: strgx,strrx,strtx

    gam1x(OPS_ACC1(0,0,0)) = strgx(OPS_ACC4(0,0,0)) - strrx(OPS_ACC5(0,0,0))
    strex(OPS_ACC2(0,0,0)) = strex(OPS_ACC2(0,0,0)) - gam1x(OPS_ACC1(0,0,0))*strtx(OPS_ACC6(0,0,0))

    gam1x(OPS_ACC1(0,0,0)) = strrx(OPS_ACC5(0,0,0))/gam1x(OPS_ACC1(0,0,0))
    ovgmx(OPS_ACC3(0,0,0)) = 1.0_8/gam1x(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqB_xdir

SUBROUTINE bounds_kernel_eqB_ydir(gam1y,strey,ovgmy,strgy,strry,strty)

    real(kind=8), dimension(1) :: gam1y,strey,ovgmy
    real(kind=8), dimension(1), intent(in) :: strgy,strry,strty

    gam1y(OPS_ACC1(0,0,0)) = strgy(OPS_ACC4(0,0,0)) - strry(OPS_ACC5(0,0,0))
    strey(OPS_ACC2(0,0,0)) = strey(OPS_ACC2(0,0,0)) - gam1y(OPS_ACC1(0,0,0))*strty(OPS_ACC6(0,0,0))

    gam1y(OPS_ACC1(0,0,0)) = strry(OPS_ACC5(0,0,0))/gam1y(OPS_ACC1(0,0,0))
    ovgmy(OPS_ACC3(0,0,0)) = 1.0_8/gam1y(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqB_ydir

SUBROUTINE bounds_kernel_eqB_zdir(gam1z,strez,ovgmz,strgz,strrz,strtz)

    real(kind=8), dimension(1) :: gam1z,strez,ovgmz
    real(kind=8), dimension(1), intent(in) :: strgz,strrz,strtz

    gam1z(OPS_ACC1(0,0,0)) = strgz(OPS_ACC4(0,0,0)) - strrz(OPS_ACC5(0,0,0))
    strez(OPS_ACC2(0,0,0)) = strez(OPS_ACC2(0,0,0)) - gam1z(OPS_ACC1(0,0,0))*strtz(OPS_ACC6(0,0,0))

    gam1z(OPS_ACC1(0,0,0)) = strrz(OPS_ACC5(0,0,0))/gam1z(OPS_ACC1(0,0,0))
    ovgmz(OPS_ACC3(0,0,0)) = 1.0_8/gam1z(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqB_zdir

SUBROUTINE bounds_kernel_eqC_xdir(acoux,ova2x,strgx,gam1x,strtx)

    real(kind=8), dimension(1) :: acoux,ova2x
    real(kind=8), dimension(1), intent(in) :: strgx,gam1x,strtx
    real(kind=8) :: fornow

    fornow = strgx(OPS_ACC3(0,0,0))*gam1x(OPS_ACC4(0,0,0))*strtx(OPS_ACC5(0,0,0))
    acoux(OPS_ACC1(0,0,0)) = sqrt(fornow)
    ova2x(OPS_ACC2(0,0,0)) = 1.0_8/fornow

END SUBROUTINE bounds_kernel_eqC_xdir

SUBROUTINE bounds_kernel_eqC_ydir(acouy,ova2y,strgy,gam1y,strty)

    real(kind=8), dimension(1) :: acouy,ova2y
    real(kind=8), dimension(1), intent(in) :: strgy,gam1y,strty
    real(kind=8) :: fornow

    fornow = strgy(OPS_ACC3(0,0,0))*gam1y(OPS_ACC4(0,0,0))*strty(OPS_ACC5(0,0,0))
    acouy(OPS_ACC1(0,0,0)) = sqrt(fornow)
    ova2y(OPS_ACC2(0,0,0)) = 1.0_8/fornow

END SUBROUTINE bounds_kernel_eqC_ydir

SUBROUTINE bounds_kernel_eqC_zdir(acouz,ova2z,strgz,gam1z,strtz)

    real(kind=8), dimension(1) :: acouz,ova2z
    real(kind=8), dimension(1), intent(in) :: strgz,gam1z,strtz
    real(kind=8) :: fornow

    fornow = strgz(OPS_ACC3(0,0,0))*gam1z(OPS_ACC4(0,0,0))*strtz(OPS_ACC5(0,0,0))
    acouz(OPS_ACC1(0,0,0)) = sqrt(fornow)
    ova2z(OPS_ACC2(0,0,0)) = 1.0_8/fornow

END SUBROUTINE bounds_kernel_eqC_zdir

!A_yz = A_yz + B_mulditim_yz*C_multidim_yz
SUBROUTINE bounds_kernel_eqD_xdir(sorpx,strhx,ratex)

    real(kind=8), dimension(1) :: sorpx
    real(kind=8), dimension(1), intent(in) :: strhx,ratex

    sorpx(OPS_ACC1(0,0,0)) = sorpx(OPS_ACC1(0,0,0)) + strhx(OPS_ACC2(0,0,0))*ratex(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_eqD_xdir

!A_xz = A_xz + B_mulditim_xz*C_multidim_xz
SUBROUTINE bounds_kernel_eqD_ydir(sorpy,strhy,ratey)

    real(kind=8), dimension(1) :: sorpy
    real(kind=8), dimension(1), intent(in) :: strhy,ratey

    sorpy(OPS_ACC1(0,0,0)) = sorpy(OPS_ACC1(0,0,0)) + strhy(OPS_ACC2(0,0,0))*ratey(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_eqD_ydir

!A_xy = A_xy + B_mulditim_xy*C_multidim_xy
SUBROUTINE bounds_kernel_eqD_zdir(sorpz,strhz,ratez)

    real(kind=8), dimension(1) :: sorpz
    real(kind=8), dimension(1), intent(in) :: strhz,ratez

    sorpz(OPS_ACC1(0,0,0)) = sorpz(OPS_ACC1(0,0,0)) + strhz(OPS_ACC2(0,0,0))*ratez(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_eqD_zdir

!A_yz = A_yz + B_mulditim_yz*val1
SUBROUTINE bounds_kernel_eqE_xdir(sydtx,dydtx,rgspec,ispec)

    real(kind=8), dimension(1) :: sydtx
    real(kind=8), dimension(1), intent(in) :: dydtx
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    sydtx(OPS_ACC1(0,0,0)) = sydtx(OPS_ACC1(0,0,0)) + dydtx(OPS_ACC2(0,0,0))*rgspec(ispec)

END SUBROUTINE bounds_kernel_eqE_xdir

!A_xz = A_xz + B_mulditim_xz*val1
SUBROUTINE bounds_kernel_eqE_ydir(sydty,dydty,rgspec,ispec)

    real(kind=8), dimension(1) :: sydty
    real(kind=8), dimension(1), intent(in) :: dydty
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    sydty(OPS_ACC1(0,0,0)) = sydty(OPS_ACC1(0,0,0)) + dydty(OPS_ACC2(0,0,0))*rgspec(ispec)

END SUBROUTINE bounds_kernel_eqE_ydir

!A_xy = A_xy + B_mulditim_xy*val1
SUBROUTINE bounds_kernel_eqE_zdir(sydtz,dydtz,rgspec,ispec)

    real(kind=8), dimension(1) :: sydtz
    real(kind=8), dimension(1), intent(in) :: dydtz
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    sydtz(OPS_ACC1(0,0,0)) = sydtz(OPS_ACC1(0,0,0)) + dydtz(OPS_ACC2(0,0,0))*rgspec(ispec)

END SUBROUTINE bounds_kernel_eqE_zdir

!A_xy = -A_xy*B_xy
SUBROUTINE bounds_kernel_eqF_xdir(sorpx,gam1x)

    real(kind=8), dimension(1) :: sorpx
    real(kind=8), dimension(1), intent(in) :: gam1x

    sorpx(OPS_ACC1(0,0,0)) = -sorpx(OPS_ACC1(0,0,0))*gam1x(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqF_xdir

!A_xz = -A_xz*B_xz
SUBROUTINE bounds_kernel_eqF_ydir(sorpy,gam1y)

    real(kind=8), dimension(1) :: sorpy
    real(kind=8), dimension(1), intent(in) :: gam1y

    sorpy(OPS_ACC1(0,0,0)) = -sorpy(OPS_ACC1(0,0,0))*gam1y(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqF_ydir

!A_xy = -A_xy*B_xy
SUBROUTINE bounds_kernel_eqF_zdir(sorpz,gam1z)

    real(kind=8), dimension(1) :: sorpz
    real(kind=8), dimension(1), intent(in) :: gam1z

    sorpz(OPS_ACC1(0,0,0)) = -sorpz(OPS_ACC1(0,0,0))*gam1z(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqF_zdir

!A_yz = A_yz/B_yz
SUBROUTINE bounds_kernel_eqG_xdir(sydtx,strrx)

    real(kind=8), dimension(1) :: sydtx
    real(kind=8), dimension(1), intent(in) :: strrx

    sydtx(OPS_ACC1(0,0,0)) = sydtx(OPS_ACC1(0,0,0))/strrx(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqG_xdir

!A_xz = A_xz/B_xz
SUBROUTINE bounds_kernel_eqG_ydir(sydty,strry)

    real(kind=8), dimension(1) :: sydty
    real(kind=8), dimension(1), intent(in) :: strry

    sydty(OPS_ACC1(0,0,0)) = sydty(OPS_ACC1(0,0,0))/strry(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqG_ydir

!A_xy = A_xy/B_xy
SUBROUTINE bounds_kernel_eqG_zdir(sydtz,strrz)

    real(kind=8), dimension(1) :: sydtz
    real(kind=8), dimension(1), intent(in) :: strrz

    sydtz(OPS_ACC1(0,0,0)) = sydtz(OPS_ACC1(0,0,0))/strrz(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqG_zdir

SUBROUTINE bounds_kernel_eqH_xl(bcl2xl,bcl3xl,bcl4xl,bcl5xl,bcl1xl,strdxl,struxl,strvxl,strwxl, &
                                strpxl,ova2xl,acouxl,sorpxl,tt5xl,cobcxl,pinfxl,xgdlen,bet,flag_bet_xl,flag_pio_xl)

    real(kind=8), dimension(1) :: bcl2xl,bcl3xl,bcl4xl,bcl5xl
    real(kind=8), dimension(1), intent(in) :: bcl1xl,strdxl,struxl,strvxl,strwxl
    real(kind=8), dimension(1), intent(in) :: strpxl,ova2xl,acouxl,sorpxl,tt5xl
    real(kind=8), intent(in) :: cobcxl,pinfxl,xgdlen,bet
    integer(kind=4), intent(in) :: flag_bet_xl,flag_pio_xl
    real(kind=8) :: bet_now

    bet_now = bet

    bcl2xl(OPS_ACC1(0,0,0)) = struxl(OPS_ACC7(0,0,0)) * (bcl2xl(OPS_ACC1(0,0,0))-bcl5xl(OPS_ACC4(0,0,0))*ova2xl(OPS_ACC11(0,0,0)))
    bcl3xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC7(0,0,0)) * bcl3xl(OPS_ACC2(0,0,0))
    bcl4xl(OPS_ACC3(0,0,0)) = struxl(OPS_ACC7(0,0,0)) * bcl4xl(OPS_ACC3(0,0,0))

!   OLD VALUE OF L5X
    bcl5xl(OPS_ACC4(0,0,0)) = 0.5_8*(struxl(OPS_ACC7(0,0,0))+acouxl(OPS_ACC12(0,0,0))) &
                            * (bcl5xl(OPS_ACC4(0,0,0)) + strdxl(OPS_ACC6(0,0,0))*acouxl(OPS_ACC12(0,0,0))*bcl1xl(OPS_ACC5(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L5X
!   TERM INTRODUCED FOR LODATO'S BC- NC
    IF ( flag_bet_xl==1 ) THEN
        bet_now = struxl(OPS_ACC7(0,0,0))*struxl(OPS_ACC7(0,0,0)) + strvxl(OPS_ACC8(0,0,0))*strvxl(OPS_ACC8(0,0,0)) &
                + strwxl(OPS_ACC9(0,0,0))*strwxl(OPS_ACC9(0,0,0))
        bet_now = SQRT(bet_now)/acouxl(OPS_ACC12(0,0,0))
    END IF

    IF ( (struxl(OPS_ACC7(0,0,0)) > 0.0_8) .AND. (flag_pio_xl==1) ) THEN

!       SUBTRACT FROM NEW VALUE OF L5X
        bcl2xl(OPS_ACC1(0,0,0)) = -bcl2xl(OPS_ACC1(0,0,0)) - ova2xl(OPS_ACC11(0,0,0))*sorpxl(OPS_ACC13(0,0,0))

        bcl3xl(OPS_ACC2(0,0,0)) = 0.1_8*(strvxl(OPS_ACC8(0,0,0))-0.0_8) - bcl3xl(OPS_ACC2(0,0,0))

        bcl4xl(OPS_ACC3(0,0,0)) = 0.1_8*(strwxl(OPS_ACC9(0,0,0))-0.0_8) - bcl4xl(OPS_ACC3(0,0,0))

        bcl5xl(OPS_ACC4(0,0,0)) = 0.5_8*sorpxl(OPS_ACC13(0,0,0)) &
                                + cobcxl*acouxl(OPS_ACC12(0,0,0))*(strpxl(OPS_ACC10(0,0,0))-pinfxl) &
                                + 0.5_8*(1.0_8-bet_now)*tt5xl(OPS_ACC14(0,0,0)) - bcl5xl(OPS_ACC4(0,0,0)) &
                                + 100.0_8*strdxl(OPS_ACC6(0,0,0))*(1.0_8/xgdlen) &
                                * (acouxl(OPS_ACC12(0,0,0))*acouxl(OPS_ACC12(0,0,0)) - struxl(OPS_ACC7(0,0,0))*struxl(OPS_ACC7(0,0,0))) &
                                * (struxl(OPS_ACC7(0,0,0))-0.0_8)
    ELSE

        bcl5xl(OPS_ACC4(0,0,0)) = 0.5_8*sorpxl(OPS_ACC13(0,0,0))  &
                                + cobcxl*acouxl(OPS_ACC12(0,0,0))*(strpxl(OPS_ACC10(0,0,0))-pinfxl)  &
                                + 0.5_8*(1.0_8-bet_now)*tt5xl(OPS_ACC14(0,0,0)) - bcl5xl(OPS_ACC4(0,0,0))
    END IF

END SUBROUTINE bounds_kernel_eqH_xl

SUBROUTINE bounds_kernel_eqH_yl(bcl2yl,bcl3yl,bcl4yl,bcl5yl,bcl1yl,strdyl,struyl,strvyl,strwyl, &
                                strpyl,ova2yl,acouyl,sorpyl,tt5yl,cobcyl,pinfyl,ygdlen,bet,flag_bet_yl,flag_pio_yl)

    real(kind=8), dimension(1) :: bcl2yl,bcl3yl,bcl4yl,bcl5yl
    real(kind=8), dimension(1), intent(in) :: bcl1yl,strdyl,struyl,strvyl,strwyl
    real(kind=8), dimension(1), intent(in) :: strpyl,ova2yl,acouyl,sorpyl,tt5yl
    real(kind=8), intent(in) :: cobcyl, pinfyl,ygdlen,bet
    integer(kind=4), intent(in) :: flag_bet_yl,flag_pio_yl
    real(kind=8) :: bet_now

    bet_now = bet

    bcl2yl(OPS_ACC1(0,0,0)) = strvyl(OPS_ACC8(0,0,0)) * (bcl2yl(OPS_ACC1(0,0,0))-bcl5yl(OPS_ACC4(0,0,0))*ova2yl(OPS_ACC11(0,0,0)))
    bcl3yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC8(0,0,0)) * bcl3yl(OPS_ACC2(0,0,0))
    bcl4yl(OPS_ACC3(0,0,0)) = strvyl(OPS_ACC8(0,0,0)) * bcl4yl(OPS_ACC3(0,0,0))

!   OLD VALUE OF L5Y
    bcl5yl(OPS_ACC4(0,0,0)) = 0.5_8*(strvyl(OPS_ACC8(0,0,0))+acouyl(OPS_ACC12(0,0,0))) &
                            * (bcl5yl(OPS_ACC4(0,0,0)) + strdyl(OPS_ACC6(0,0,0))*acouyl(OPS_ACC12(0,0,0))*bcl1yl(OPS_ACC5(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L5Y
!   TERMS INTRODUCED FOR LODATO'S BC- NC
    IF ( flag_bet_yl==1 ) THEN
        bet_now = struyl(OPS_ACC7(0,0,0))*struyl(OPS_ACC7(0,0,0)) + strvyl(OPS_ACC8(0,0,0))*strvyl(OPS_ACC8(0,0,0)) &
                + strwyl(OPS_ACC9(0,0,0))*strwyl(OPS_ACC9(0,0,0))
        bet_now = SQRT(bet_now)/acouyl(OPS_ACC12(0,0,0))
    END IF

    IF ( (strvyl(OPS_ACC8(0,0,0)) > 0.0_8) .AND. (flag_pio_yl==1) ) THEN

        bcl2yl(OPS_ACC1(0,0,0)) = -bcl2yl(OPS_ACC1(0,0,0)) - ova2yl(OPS_ACC11(0,0,0))*sorpyl(OPS_ACC13(0,0,0))

        bcl3yl(OPS_ACC2(0,0,0)) = 0.1_8*(struyl(OPS_ACC7(0,0,0))-0.0_8) - bcl3yl(OPS_ACC2(0,0,0))

        bcl4yl(OPS_ACC3(0,0,0)) = 0.1_8*(strwyl(OPS_ACC9(0,0,0))-0.0_8) - bcl4yl(OPS_ACC3(0,0,0))

        bcl5yl(OPS_ACC4(0,0,0)) = 0.5_8*sorpyl(OPS_ACC13(0,0,0)) &
                                + cobcyl*acouyl(OPS_ACC12(0,0,0))*(strpyl(OPS_ACC10(0,0,0))-pinfyl) &
                                + 0.5_8*(1.0-bet_now)*tt5yl(OPS_ACC14(0,0,0)) - bcl5yl(OPS_ACC4(0,0,0)) &
                                + 100.0_8*strdyl(OPS_ACC6(0,0,0))*(1.0_8/ygdlen) &
                                * (acouyl(OPS_ACC12(0,0,0))*acouyl(OPS_ACC12(0,0,0)) - strvyl(OPS_ACC8(0,0,0))*strvyl(OPS_ACC8(0,0,0))) &
                                * (strvyl(OPS_ACC8(0,0,0))-0.0_8)

    ELSE

        bcl5yl(OPS_ACC4(0,0,0)) = 0.5_8*sorpyl(OPS_ACC13(0,0,0)) &
                                + cobcyl*acouyl(OPS_ACC12(0,0,0))*(strpyl(OPS_ACC10(0,0,0))-pinfyl) &
                                + 0.5_8*(1.0_8-bet_now)*tt5yl(OPS_ACC14(0,0,0)) - bcl5yl(OPS_ACC4(0,0,0))
    END IF

END SUBROUTINE bounds_kernel_eqH_yl

SUBROUTINE bounds_kernel_eqH_zl(bcl2zl,bcl3zl,bcl4zl,bcl5zl,bcl1zl,strdzl,struzl,strvzl,strwzl, &
                                strpzl,ova2zl,acouzl,sorpzl,tt5zl,cobczl,pinfzl,zgdlen,bet,flag_bet_zl,flag_pio_zl)

    real(kind=8), dimension(1) :: bcl2zl,bcl3zl,bcl4zl,bcl5zl
    real(kind=8), dimension(1), intent(in) :: bcl1zl,strdzl,struzl,strvzl,strwzl
    real(kind=8), dimension(1), intent(in) :: strpzl,ova2zl,acouzl,sorpzl,tt5zl
    real(kind=8), intent(in) :: cobczl,pinfzl,zgdlen,bet
    integer(kind=4), intent(in) :: flag_bet_zl,flag_pio_zl
    real(kind=8) :: bet_now

    bet_now = bet

    bcl2zl(OPS_ACC1(0,0,0)) = strwzl(OPS_ACC9(0,0,0)) * (bcl2zl(OPS_ACC1(0,0,0))-bcl5zl(OPS_ACC4(0,0,0))*ova2zl(OPS_ACC11(0,0,0)))
    bcl3zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC9(0,0,0)) * bcl3zl(OPS_ACC2(0,0,0))
    bcl4zl(OPS_ACC3(0,0,0)) = strwzl(OPS_ACC9(0,0,0)) * bcl4zl(OPS_ACC3(0,0,0))

!   OLD VALUE OF L5Z
    bcl5zl(OPS_ACC4(0,0,0)) = 0.5_8*(strwzl(OPS_ACC9(0,0,0))+acouzl(OPS_ACC12(0,0,0))) &
                            * (bcl5zl(OPS_ACC4(0,0,0)) + strdzl(OPS_ACC6(0,0,0))*acouzl(OPS_ACC12(0,0,0))*bcl1zl(OPS_ACC5(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L5Z
!   TERMS INTRODUCED FOR LODATO'S BC- NC
    IF ( flag_bet_zl==1 ) THEN
        bet_now = struzl(OPS_ACC7(0,0,0))*struzl(OPS_ACC7(0,0,0)) + strvzl(OPS_ACC8(0,0,0))*strvzl(OPS_ACC8(0,0,0)) &
              + strwzl(OPS_ACC9(0,0,0))*strwzl(OPS_ACC9(0,0,0))
        bet_now = SQRT(bet_now)/acouzl(OPS_ACC12(0,0,0))
    END IF

    IF ( (strwzl(OPS_ACC9(0,0,0)) > 0.0_8) .AND. (flag_pio_zl==1) ) THEN

        bcl2zl(OPS_ACC1(0,0,0)) = -bcl2zl(OPS_ACC1(0,0,0)) - ova2zl(OPS_ACC11(0,0,0))*sorpzl(OPS_ACC13(0,0,0))

        bcl3zl(OPS_ACC2(0,0,0)) = 0.1_8*(struzl(OPS_ACC7(0,0,0))-0.0_8) - bcl3zl(OPS_ACC2(0,0,0))

        bcl4zl(OPS_ACC3(0,0,0)) = 0.1_8*(strvzl(OPS_ACC8(0,0,0))-0.0_8) - bcl4zl(OPS_ACC3(0,0,0))

        bcl5zl(OPS_ACC4(0,0,0)) = 0.5_8*sorpzl(OPS_ACC13(0,0,0)) &
                                + cobczl*acouzl(OPS_ACC12(0,0,0))*(strpzl(OPS_ACC10(0,0,0))-pinfzl) &
                                + 0.5_8*(1.0_8-bet_now)*tt5zl(OPS_ACC14(0,0,0)) - bcl5zl(OPS_ACC4(0,0,0)) &
                                + 100.0_8*strdzl(OPS_ACC6(0,0,0))*(1.0_8/zgdlen) &
                                * (acouzl(OPS_ACC12(0,0,0))*acouzl(OPS_ACC12(0,0,0)) - strwzl(OPS_ACC9(0,0,0))*strwzl(OPS_ACC9(0,0,0))) &
                                * (strwzl(OPS_ACC9(0,0,0))-0.0_8)
    ELSE

        bcl5zl(OPS_ACC4(0,0,0)) = 0.5_8*sorpzl(OPS_ACC13(0,0,0)) &
                                + cobczl*acouzl(OPS_ACC12(0,0,0))*(strpzl(OPS_ACC10(0,0,0))-pinfzl) &
                                + 0.5_8*(1.0_8-bet_now)*tt5zl(OPS_ACC14(0,0,0)) - bcl5zl(OPS_ACC4(0,0,0))
    END IF

END SUBROUTINE bounds_kernel_eqH_zl

SUBROUTINE bounds_kernel_eqI_xl(bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,bcl1xl,struxl,ova2xl,sorpxl,strpxl,cobcxl,pinfxl)

    real(kind=8), dimension(1) :: bcl2xl,bcl3xl,bcl4xl,bcl5xl
    real(kind=8), dimension(1), intent(in) :: strdxl,acouxl,bcl1xl,struxl,ova2xl,sorpxl,strpxl
    real(kind=8), intent(in) :: cobcxl,pinfxl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC5(0,0,0))*acouxl(OPS_ACC6(0,0,0))*bcl1xl(OPS_ACC7(0,0,0))
    bcl2xl(OPS_ACC1(0,0,0)) = struxl(OPS_ACC8(0,0,0))  &
            *(bcl2xl(OPS_ACC1(0,0,0))-bcl5xl(OPS_ACC4(0,0,0))*ova2xl(OPS_ACC9(0,0,0)))
    bcl3xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl3xl(OPS_ACC2(0,0,0))
    bcl4xl(OPS_ACC3(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl4xl(OPS_ACC3(0,0,0))
    bcl5xl(OPS_ACC4(0,0,0)) = 0.5_8*(struxl(OPS_ACC8(0,0,0))+acouxl(OPS_ACC6(0,0,0)))  &
            *(bcl5xl(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
!   L1X UNCHANGED
    bcl2xl(OPS_ACC1(0,0,0)) = -bcl2xl(OPS_ACC1(0,0,0))
    bcl3xl(OPS_ACC2(0,0,0)) = -bcl3xl(OPS_ACC2(0,0,0))
    bcl4xl(OPS_ACC3(0,0,0)) = -bcl4xl(OPS_ACC3(0,0,0))
    bcl5xl(OPS_ACC4(0,0,0)) = 0.5_8*sorpxl(OPS_ACC10(0,0,0))  &
            + cobcxl*acouxl(OPS_ACC6(0,0,0))*(strpxl(OPS_ACC11(0,0,0))-pinfxl) - bcl5xl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqI_xl

SUBROUTINE bounds_kernel_eqI_yl(bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,bcl1yl,strvyl,ova2yl,sorpyl,strpyl,cobcyl,pinfyl)

    real(kind=8), dimension(1) :: bcl2yl,bcl3yl,bcl4yl,bcl5yl
    real(kind=8), dimension(1), intent(in) :: strdyl,acouyl,bcl1yl,strvyl,ova2yl,sorpyl,strpyl
    real(kind=8), intent(in) :: cobcyl,pinfyl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC5(0,0,0))*acouyl(OPS_ACC6(0,0,0))*bcl1yl(OPS_ACC7(0,0,0))
    bcl2yl(OPS_ACC1(0,0,0)) = strvyl(OPS_ACC8(0,0,0))  &
            *(bcl2yl(OPS_ACC1(0,0,0))-bcl5yl(OPS_ACC4(0,0,0))*ova2yl(OPS_ACC9(0,0,0)))
    bcl3yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl3yl(OPS_ACC2(0,0,0))
    bcl4yl(OPS_ACC3(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl4yl(OPS_ACC3(0,0,0))
    bcl5yl(OPS_ACC4(0,0,0)) = 0.5_8*(strvyl(OPS_ACC8(0,0,0))+acouyl(OPS_ACC6(0,0,0)))  &
            *(bcl5yl(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
!   L1Y UNCHANGED
    bcl2yl(OPS_ACC1(0,0,0)) = -bcl2yl(OPS_ACC1(0,0,0))
    bcl3yl(OPS_ACC2(0,0,0)) = -bcl3yl(OPS_ACC2(0,0,0))
    bcl4yl(OPS_ACC3(0,0,0)) = -bcl4yl(OPS_ACC3(0,0,0))
    bcl5yl(OPS_ACC4(0,0,0)) = 0.5_8*sorpyl(OPS_ACC10(0,0,0))  &
            + cobcyl*acouyl(OPS_ACC6(0,0,0))*(strpyl(OPS_ACC11(0,0,0))-pinfyl) - bcl5yl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqI_yl

SUBROUTINE bounds_kernel_eqI_zl(bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,bcl1zl,strwzl,ova2zl,sorpzl,strpzl,cobczl,pinfzl)

    real(kind=8), dimension(1) :: bcl2zl,bcl3zl,bcl4zl,bcl5zl
    real(kind=8), dimension(1), intent(in) :: strdzl,acouzl,bcl1zl,strwzl,ova2zl,sorpzl,strpzl
    real(kind=8), intent(in) :: cobczl,pinfzl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC5(0,0,0))*acouzl(OPS_ACC6(0,0,0))*bcl1zl(OPS_ACC7(0,0,0))
    bcl2zl(OPS_ACC1(0,0,0)) = strwzl(OPS_ACC8(0,0,0))  &
            *(bcl2zl(OPS_ACC1(0,0,0))-bcl5zl(OPS_ACC4(0,0,0))*ova2zl(OPS_ACC9(0,0,0)))
    bcl3zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl3zl(OPS_ACC2(0,0,0))
    bcl4zl(OPS_ACC3(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl4zl(OPS_ACC3(0,0,0))
    bcl5zl(OPS_ACC4(0,0,0)) = 0.5_8*(strwzl(OPS_ACC8(0,0,0))+acouzl(OPS_ACC6(0,0,0)))  &
            *(bcl5zl(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
!   L1Z UNCHANGED
    bcl2zl(OPS_ACC1(0,0,0)) = -bcl2zl(OPS_ACC1(0,0,0))
    bcl3zl(OPS_ACC2(0,0,0)) = -bcl3zl(OPS_ACC2(0,0,0))
    bcl4zl(OPS_ACC3(0,0,0)) = -bcl4zl(OPS_ACC3(0,0,0))
    bcl5zl(OPS_ACC4(0,0,0)) = 0.5_8*sorpzl(OPS_ACC10(0,0,0))  &
            + cobczl*acouzl(OPS_ACC6(0,0,0))*(strpzl(OPS_ACC11(0,0,0))-pinfzl) - bcl5zl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqI_zl

SUBROUTINE bounds_kernel_eqJ_xl(bclyxl,ratexl,struxl,strdxl)

    real(kind=8), dimension(1) :: bclyxl
    real(kind=8), dimension(1), intent(in) :: ratexl,struxl,strdxl

!   OLD VALUE OF L's
    bclyxl(OPS_ACC1(0,0,0)) = struxl(OPS_ACC3(0,0,0))*bclyxl(OPS_ACC1(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
    bclyxl(OPS_ACC1(0,0,0)) = ratexl(OPS_ACC2(0,0,0))/strdxl(OPS_ACC4(0,0,0))  &
              - bclyxl(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqJ_xl

SUBROUTINE bounds_kernel_eqJ_yl(bclyyl,rateyl,strvyl,strdyl)

    real(kind=8), dimension(1) :: bclyyl
    real(kind=8), dimension(1), intent(in) :: rateyl,strvyl,strdyl

!   OLD VALUE OF L's
    bclyyl(OPS_ACC1(0,0,0)) = strvyl(OPS_ACC3(0,0,0))*bclyyl(OPS_ACC1(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
    bclyyl(OPS_ACC1(0,0,0)) = rateyl(OPS_ACC2(0,0,0))/strdyl(OPS_ACC4(0,0,0))  &
              - bclyyl(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqJ_yl

SUBROUTINE bounds_kernel_eqJ_zl(bclyzl,ratezl,strwzl,strdzl)

    real(kind=8), dimension(1) :: bclyzl
    real(kind=8), dimension(1), intent(in) :: ratezl,strwzl,strdzl

!   OLD VALUE OF L's
    bclyzl(OPS_ACC1(0,0,0)) = strwzl(OPS_ACC3(0,0,0))*bclyzl(OPS_ACC1(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
    bclyzl(OPS_ACC1(0,0,0)) = ratezl(OPS_ACC2(0,0,0))/strdzl(OPS_ACC4(0,0,0))  &
              - bclyzl(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqJ_zl

SUBROUTINE bounds_kernel_eqK_xl(bcl1xl,bcl2xl,bcl5xl,strdxl,acouxl,struxl,ova2xl,dudtxl,gam1xl,dtdtxl,strtxl,sorpxl,strpxl,sydtxl)

    real(kind=8), dimension(1) :: bcl1xl,bcl2xl,bcl5xl
    real(kind=8), dimension(1), intent(in) :: strdxl,acouxl,struxl,ova2xl,dudtxl,gam1xl,dtdtxl,strtxl,sorpxl,strpxl,sydtxl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC4(0,0,0))*acouxl(OPS_ACC5(0,0,0))*bcl1xl(OPS_ACC1(0,0,0))
    bcl1xl(OPS_ACC1(0,0,0)) = 0.5_8*(struxl(OPS_ACC6(0,0,0))-acouxl(OPS_ACC5(0,0,0)))  &
            *(bcl5xl(OPS_ACC3(0,0,0))-fornow)
    bcl2xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC6(0,0,0))  &
            *(bcl2xl(OPS_ACC2(0,0,0))-bcl5xl(OPS_ACC3(0,0,0))*ova2xl(OPS_ACC7(0,0,0)))
    bcl5xl(OPS_ACC3(0,0,0)) = 0.5_8*(struxl(OPS_ACC6(0,0,0))+acouxl(OPS_ACC5(0,0,0)))  &
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

END SUBROUTINE bounds_kernel_eqK_xl

SUBROUTINE bounds_kernel_eqK_yl(bcl1yl,bcl2yl,bcl5yl,strdyl,acouyl,strvyl,ova2yl,dvdtyl,gam1yl,dtdtyl,strtyl,sorpyl,strpyl,sydtyl)

    real(kind=8), dimension(1) :: bcl1yl,bcl2yl,bcl5yl
    real(kind=8), dimension(1), intent(in) :: strdyl,acouyl,strvyl,ova2yl,dvdtyl,gam1yl,dtdtyl,strtyl,sorpyl,strpyl,sydtyl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC4(0,0,0))*acouyl(OPS_ACC5(0,0,0))*bcl1yl(OPS_ACC1(0,0,0))
    bcl1yl(OPS_ACC1(0,0,0)) = 0.5_8*(strvyl(OPS_ACC6(0,0,0))-acouyl(OPS_ACC5(0,0,0)))  &
            *(bcl5yl(OPS_ACC3(0,0,0))-fornow)
    bcl2yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC6(0,0,0))  &
            *(bcl2yl(OPS_ACC2(0,0,0))-bcl5yl(OPS_ACC3(0,0,0))*ova2yl(OPS_ACC7(0,0,0)))
    bcl5yl(OPS_ACC3(0,0,0)) = 0.5_8*(strvyl(OPS_ACC6(0,0,0))+acouyl(OPS_ACC5(0,0,0)))  &
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

END SUBROUTINE bounds_kernel_eqK_yl

SUBROUTINE bounds_kernel_eqK_zl(bcl1zl,bcl2zl,bcl5zl,strdzl,acouzl,strwzl,ova2zl,dwdtzl,gam1zl,dtdtzl,strtzl,sorpzl,strpzl,sydtzl)

    real(kind=8), dimension(1) :: bcl1zl,bcl2zl,bcl5zl
    real(kind=8), dimension(1), intent(in) :: strdzl,acouzl,strwzl,ova2zl,dwdtzl,gam1zl,dtdtzl,strtzl,sorpzl,strpzl,sydtzl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC4(0,0,0))*acouzl(OPS_ACC5(0,0,0))*bcl1zl(OPS_ACC1(0,0,0))
    bcl1zl(OPS_ACC1(0,0,0)) = 0.5_8*(strwzl(OPS_ACC6(0,0,0))-acouzl(OPS_ACC5(0,0,0)))  &
            *(bcl5zl(OPS_ACC3(0,0,0))-fornow)
    bcl2zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC6(0,0,0))  &
            *(bcl2zl(OPS_ACC2(0,0,0))-bcl5zl(OPS_ACC3(0,0,0))*ova2zl(OPS_ACC7(0,0,0)))
    bcl5zl(OPS_ACC3(0,0,0)) = 0.5_8*(strwzl(OPS_ACC6(0,0,0))+acouzl(OPS_ACC5(0,0,0)))  &
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

END SUBROUTINE bounds_kernel_eqK_zl

SUBROUTINE bounds_kernel_eqL_xl(bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,struxl,ova2xl,dudtxl,dddtxl,dvdtxl,dwdtxl)

    real(kind=8), dimension(1) :: bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl
    real(kind=8), dimension(1), intent(in) :: strdxl,acouxl,struxl,ova2xl,dudtxl,dddtxl,dvdtxl,dwdtxl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC6(0,0,0))*acouxl(OPS_ACC7(0,0,0))*bcl1xl(OPS_ACC1(0,0,0))
    bcl1xl(OPS_ACC1(0,0,0)) = 0.5_8*(struxl(OPS_ACC8(0,0,0))-acouxl(OPS_ACC7(0,0,0)))  &
            *(bcl5xl(OPS_ACC5(0,0,0))-fornow)
    bcl2xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC8(0,0,0))  &
            *(bcl2xl(OPS_ACC2(0,0,0))-bcl5xl(OPS_ACC5(0,0,0))*ova2xl(OPS_ACC9(0,0,0)))
    bcl3xl(OPS_ACC3(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl3xl(OPS_ACC3(0,0,0))
    bcl4xl(OPS_ACC4(0,0,0)) = struxl(OPS_ACC8(0,0,0))*bcl4xl(OPS_ACC4(0,0,0))
    bcl5xl(OPS_ACC5(0,0,0)) = 0.5_8*(struxl(OPS_ACC8(0,0,0))+acouxl(OPS_ACC7(0,0,0)))  &
            *(bcl5xl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1X UNCHANGED
    fornow = bcl1xl(OPS_ACC1(0,0,0)) - strdxl(OPS_ACC6(0,0,0))*acouxl(OPS_ACC7(0,0,0))*dudtxl(OPS_ACC10(0,0,0))
    bcl2xl(OPS_ACC2(0,0,0)) = -dddtxl(OPS_ACC11(0,0,0))  &
            - ova2xl(OPS_ACC9(0,0,0))*(bcl1xl(OPS_ACC1(0,0,0))+fornow) - bcl2xl(OPS_ACC2(0,0,0))
    bcl3xl(OPS_ACC3(0,0,0)) = -dvdtxl(OPS_ACC12(0,0,0)) - bcl3xl(OPS_ACC3(0,0,0))
    bcl4xl(OPS_ACC4(0,0,0)) = -dwdtxl(OPS_ACC13(0,0,0)) - bcl4xl(OPS_ACC4(0,0,0))
    bcl5xl(OPS_ACC5(0,0,0)) = fornow - bcl5xl(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqL_xl

SUBROUTINE bounds_kernel_eqL_yl(bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,strvyl,ova2yl,dudtyl,dddtyl,dvdtyl,dwdtyl)

    real(kind=8), dimension(1) :: bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl
    real(kind=8), dimension(1), intent(in) :: strdyl,acouyl,strvyl,ova2yl,dudtyl,dddtyl,dvdtyl,dwdtyl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC6(0,0,0))*acouyl(OPS_ACC7(0,0,0))*bcl1yl(OPS_ACC1(0,0,0))
    bcl1yl(OPS_ACC1(0,0,0)) = 0.5_8*(strvyl(OPS_ACC8(0,0,0))-acouyl(OPS_ACC7(0,0,0)))  &
            *(bcl5yl(OPS_ACC5(0,0,0))-fornow)
    bcl2yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC8(0,0,0))  &
            *(bcl2yl(OPS_ACC2(0,0,0))-bcl5yl(OPS_ACC5(0,0,0))*ova2yl(OPS_ACC9(0,0,0)))
    bcl3yl(OPS_ACC3(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl3yl(OPS_ACC3(0,0,0))
    bcl4yl(OPS_ACC4(0,0,0)) = strvyl(OPS_ACC8(0,0,0))*bcl4yl(OPS_ACC4(0,0,0))
    bcl5yl(OPS_ACC5(0,0,0)) = 0.5_8*(strvyl(OPS_ACC8(0,0,0))+acouyl(OPS_ACC7(0,0,0)))  &
            *(bcl5yl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Y UNCHANGED
    fornow = bcl1yl(OPS_ACC1(0,0,0)) - strdyl(OPS_ACC6(0,0,0))*acouyl(OPS_ACC7(0,0,0))*dvdtyl(OPS_ACC10(0,0,0))
    bcl2yl(OPS_ACC2(0,0,0)) = -dddtyl(OPS_ACC11(0,0,0))  &
            - ova2yl(OPS_ACC9(0,0,0))*(bcl1yl(OPS_ACC1(0,0,0))+fornow) - bcl2yl(OPS_ACC2(0,0,0))
    bcl3yl(OPS_ACC3(0,0,0)) = -dudtyl(OPS_ACC12(0,0,0)) - bcl3yl(OPS_ACC3(0,0,0))
    bcl4yl(OPS_ACC4(0,0,0)) = -dwdtyl(OPS_ACC13(0,0,0)) - bcl4yl(OPS_ACC4(0,0,0))
    bcl5yl(OPS_ACC5(0,0,0)) = fornow - bcl5yl(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqL_yl

SUBROUTINE bounds_kernel_eqL_zl(bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,strwzl,ova2zl,dudtzl,dddtzl,dvdtzl,dwdtzl)

    real(kind=8), dimension(1) :: bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl
    real(kind=8), dimension(1), intent(in) :: strdzl,acouzl,strwzl,ova2zl,dudtzl,dddtzl,dvdtzl,dwdtzl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC6(0,0,0))*acouzl(OPS_ACC7(0,0,0))*bcl1zl(OPS_ACC1(0,0,0))
    bcl1zl(OPS_ACC1(0,0,0)) = 0.5_8*(strwzl(OPS_ACC8(0,0,0))-acouzl(OPS_ACC7(0,0,0)))  &
            *(bcl5zl(OPS_ACC5(0,0,0))-fornow)
    bcl2zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC8(0,0,0))  &
            *(bcl2zl(OPS_ACC2(0,0,0))-bcl5zl(OPS_ACC5(0,0,0))*ova2zl(OPS_ACC9(0,0,0)))
    bcl3zl(OPS_ACC3(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl3zl(OPS_ACC3(0,0,0))
    bcl4zl(OPS_ACC4(0,0,0)) = strwzl(OPS_ACC8(0,0,0))*bcl4zl(OPS_ACC4(0,0,0))
    bcl5zl(OPS_ACC5(0,0,0)) = 0.5_8*(strwzl(OPS_ACC8(0,0,0))+acouzl(OPS_ACC7(0,0,0)))  &
            *(bcl5zl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Z UNCHANGED
    fornow = bcl1zl(OPS_ACC1(0,0,0)) - strdzl(OPS_ACC6(0,0,0))*acouzl(OPS_ACC7(0,0,0))*dwdtzl(OPS_ACC10(0,0,0))
    bcl2zl(OPS_ACC2(0,0,0)) = -dddtzl(OPS_ACC11(0,0,0))  &
            - ova2zl(OPS_ACC9(0,0,0))*(bcl1zl(OPS_ACC1(0,0,0))+fornow) - bcl2zl(OPS_ACC2(0,0,0))
    bcl3zl(OPS_ACC3(0,0,0)) = -dudtzl(OPS_ACC12(0,0,0)) - bcl3zl(OPS_ACC3(0,0,0))
    bcl4zl(OPS_ACC4(0,0,0)) = -dvdtzl(OPS_ACC13(0,0,0)) - bcl4zl(OPS_ACC4(0,0,0))
    bcl5zl(OPS_ACC5(0,0,0)) = fornow - bcl5zl(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqL_zl

SUBROUTINE bounds_kernel_eqM_xl(bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,struxl,ova2xl)

    real(kind=8), dimension(1) :: bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl
    real(kind=8), dimension(1), intent(in) :: strdxl,acouxl,struxl,ova2xl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC6(0,0,0))*acouxl(OPS_ACC7(0,0,0))*bcl1xl(OPS_ACC1(0,0,0))
    bcl1xl(OPS_ACC1(0,0,0)) = 0.5_8*(struxl(OPS_ACC8(0,0,0))-acouxl(OPS_ACC7(0,0,0)))  &
                            * (bcl5xl(OPS_ACC5(0,0,0))-fornow)
    bcl2xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC8(0,0,0))  &
                            * (bcl2xl(OPS_ACC2(0,0,0))-bcl5xl(OPS_ACC5(0,0,0))*ova2xl(OPS_ACC9(0,0,0)))
    bcl3xl(OPS_ACC3(0,0,0)) = struxl(OPS_ACC8(0,0,0)) * bcl3xl(OPS_ACC3(0,0,0))
    bcl4xl(OPS_ACC4(0,0,0)) = struxl(OPS_ACC8(0,0,0)) * bcl4xl(OPS_ACC4(0,0,0))
    bcl5xl(OPS_ACC5(0,0,0)) = 0.5_8*(struxl(OPS_ACC8(0,0,0))+acouxl(OPS_ACC7(0,0,0)))  &
                            * (bcl5xl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1X UNCHANGED
    bcl3xl(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4xl(OPS_ACC4(0,0,0)) = 0.0_8
!   NC SUGGESTION: IF NOT WORKING REMOVE SECOND LINE AS WELL
    bcl5xl(OPS_ACC5(0,0,0)) = bcl1xl(OPS_ACC1(0,0,0)) - bcl5xl(OPS_ACC5(0,0,0))
    bcl2xl(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqM_xl

SUBROUTINE bounds_kernel_eqM_yl(bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,strvyl,ova2yl)

    real(kind=8), dimension(1) :: bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl
    real(kind=8), dimension(1), intent(in) :: strdyl,acouyl,strvyl,ova2yl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC6(0,0,0))*acouyl(OPS_ACC7(0,0,0))*bcl1yl(OPS_ACC1(0,0,0))
    bcl1yl(OPS_ACC1(0,0,0)) = 0.5_8*(strvyl(OPS_ACC8(0,0,0))-acouyl(OPS_ACC7(0,0,0)))  &
                            * (bcl5yl(OPS_ACC5(0,0,0))-fornow)
    bcl2yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC8(0,0,0))  &
                            * (bcl2yl(OPS_ACC2(0,0,0))-bcl5yl(OPS_ACC5(0,0,0))*ova2yl(OPS_ACC9(0,0,0)))
    bcl3yl(OPS_ACC3(0,0,0)) = strvyl(OPS_ACC8(0,0,0)) * bcl3yl(OPS_ACC3(0,0,0))
    bcl4yl(OPS_ACC4(0,0,0)) = strvyl(OPS_ACC8(0,0,0)) * bcl4yl(OPS_ACC4(0,0,0))
    bcl5yl(OPS_ACC5(0,0,0)) = 0.5_8*(strvyl(OPS_ACC8(0,0,0))+acouyl(OPS_ACC7(0,0,0)))  &
                            * (bcl5yl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Y UNCHANGED
    bcl3yl(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4yl(OPS_ACC4(0,0,0)) = 0.0_8
    bcl5yl(OPS_ACC5(0,0,0)) = bcl1yl(OPS_ACC1(0,0,0)) - bcl5yl(OPS_ACC5(0,0,0))
    bcl2yl(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqM_yl

SUBROUTINE bounds_kernel_eqM_zl(bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,strwzl,ova2zl)

    real(kind=8), dimension(1) :: bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl
    real(kind=8), dimension(1), intent(in) :: strdzl,acouzl,strwzl,ova2zl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC6(0,0,0))*acouzl(OPS_ACC7(0,0,0))*bcl1zl(OPS_ACC1(0,0,0))
    bcl1zl(OPS_ACC1(0,0,0)) = 0.5_8*(strwzl(OPS_ACC8(0,0,0))-acouzl(OPS_ACC7(0,0,0)))  &
                            * (bcl5zl(OPS_ACC5(0,0,0))-fornow)
    bcl2zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC8(0,0,0))  &
                            * (bcl2zl(OPS_ACC2(0,0,0))-bcl5zl(OPS_ACC5(0,0,0))*ova2zl(OPS_ACC9(0,0,0)))
    bcl3zl(OPS_ACC3(0,0,0)) = strwzl(OPS_ACC8(0,0,0)) * bcl3zl(OPS_ACC3(0,0,0))
    bcl4zl(OPS_ACC4(0,0,0)) = strwzl(OPS_ACC8(0,0,0)) * bcl4zl(OPS_ACC4(0,0,0))
    bcl5zl(OPS_ACC5(0,0,0)) = 0.5_8*(strwzl(OPS_ACC8(0,0,0))+acouzl(OPS_ACC7(0,0,0)))  &
                            * (bcl5zl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Y UNCHANGED
    bcl3zl(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4zl(OPS_ACC4(0,0,0)) = 0.0_8
    bcl5zl(OPS_ACC5(0,0,0)) = bcl1zl(OPS_ACC1(0,0,0)) - bcl5zl(OPS_ACC5(0,0,0))
    bcl2zl(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqM_zl

SUBROUTINE bounds_kernel_eqN_xl(bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,acouxl,struxl,ova2xl)

    real(kind=8), dimension(1) :: bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl
    real(kind=8), dimension(1), intent(in) :: strdxl,acouxl,struxl,ova2xl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxl(OPS_ACC6(0,0,0)) * acouxl(OPS_ACC7(0,0,0)) * bcl1xl(OPS_ACC1(0,0,0))
    bcl1xl(OPS_ACC1(0,0,0)) = 0.5_8 * (struxl(OPS_ACC8(0,0,0))-acouxl(OPS_ACC7(0,0,0)))  &
                            * (bcl5xl(OPS_ACC5(0,0,0))-fornow)
    bcl2xl(OPS_ACC2(0,0,0)) = struxl(OPS_ACC8(0,0,0))  &
                            * (bcl2xl(OPS_ACC2(0,0,0))-bcl5xl(OPS_ACC5(0,0,0))*ova2xl(OPS_ACC9(0,0,0)))
    bcl3xl(OPS_ACC3(0,0,0)) = struxl(OPS_ACC8(0,0,0)) * bcl3xl(OPS_ACC3(0,0,0))
    bcl4xl(OPS_ACC4(0,0,0)) = struxl(OPS_ACC8(0,0,0)) * bcl4xl(OPS_ACC4(0,0,0))
    bcl5xl(OPS_ACC5(0,0,0)) = 0.5_8 * (struxl(OPS_ACC8(0,0,0))+acouxl(OPS_ACC7(0,0,0)))  &
                            * (bcl5xl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1X UNCHANGED
    bcl3xl(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4xl(OPS_ACC4(0,0,0)) = 0.0_8
!   NC SUGGESTION: IF NOT WORKING REMOVE SECOND LINE AS WELL
    bcl5xl(OPS_ACC5(0,0,0)) = bcl1xl(OPS_ACC1(0,0,0)) - bcl5xl(OPS_ACC5(0,0,0))
    bcl2xl(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqN_xl

SUBROUTINE bounds_kernel_eqN_yl(bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,acouyl,strvyl,ova2yl)

    real(kind=8), dimension(1) :: bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl
    real(kind=8), dimension(1), intent(in) :: strdyl,acouyl,strvyl,ova2yl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyl(OPS_ACC6(0,0,0)) * acouyl(OPS_ACC7(0,0,0)) * bcl1yl(OPS_ACC1(0,0,0))
    bcl1yl(OPS_ACC1(0,0,0)) = 0.5_8 * (strvyl(OPS_ACC8(0,0,0))-acouyl(OPS_ACC7(0,0,0)))  &
                            * (bcl5yl(OPS_ACC5(0,0,0))-fornow)
    bcl2yl(OPS_ACC2(0,0,0)) = strvyl(OPS_ACC8(0,0,0))  &
                            * (bcl2yl(OPS_ACC2(0,0,0))-bcl5yl(OPS_ACC5(0,0,0))*ova2yl(OPS_ACC9(0,0,0)))
    bcl3yl(OPS_ACC3(0,0,0)) = strvyl(OPS_ACC8(0,0,0)) * bcl3yl(OPS_ACC3(0,0,0))
    bcl4yl(OPS_ACC4(0,0,0)) = strvyl(OPS_ACC8(0,0,0)) * bcl4yl(OPS_ACC4(0,0,0))
    bcl5yl(OPS_ACC5(0,0,0)) = 0.5_8 * (strvyl(OPS_ACC8(0,0,0))+acouyl(OPS_ACC7(0,0,0)))  &
                            * (bcl5yl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Y UNCHANGED
    bcl3yl(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4yl(OPS_ACC4(0,0,0)) = 0.0_8
    bcl5yl(OPS_ACC5(0,0,0)) = bcl1yl(OPS_ACC1(0,0,0)) - bcl5yl(OPS_ACC5(0,0,0))
    bcl2yl(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqN_yl

SUBROUTINE bounds_kernel_eqN_zl(bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,acouzl,strwzl,ova2zl)

    real(kind=8), dimension(1) :: bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl
    real(kind=8), dimension(1), intent(in) :: strdzl,acouzl,strwzl,ova2zl
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzl(OPS_ACC6(0,0,0)) * acouzl(OPS_ACC7(0,0,0)) * bcl1zl(OPS_ACC1(0,0,0))
    bcl1zl(OPS_ACC1(0,0,0)) = 0.5_8 * (strwzl(OPS_ACC8(0,0,0))-acouzl(OPS_ACC7(0,0,0)))  &
                            * (bcl5zl(OPS_ACC5(0,0,0))-fornow)
    bcl2zl(OPS_ACC2(0,0,0)) = strwzl(OPS_ACC8(0,0,0))  &
                            * (bcl2zl(OPS_ACC2(0,0,0))-bcl5zl(OPS_ACC5(0,0,0))*ova2zl(OPS_ACC9(0,0,0)))
    bcl3zl(OPS_ACC3(0,0,0)) = strwzl(OPS_ACC8(0,0,0)) * bcl3zl(OPS_ACC3(0,0,0))
    bcl4zl(OPS_ACC4(0,0,0)) = strwzl(OPS_ACC8(0,0,0)) * bcl4zl(OPS_ACC4(0,0,0))
    bcl5zl(OPS_ACC5(0,0,0)) = 0.5_8 * (strwzl(OPS_ACC8(0,0,0))+acouzl(OPS_ACC7(0,0,0)))  &
                            * (bcl5zl(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1Z UNCHANGED
    bcl3zl(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4zl(OPS_ACC4(0,0,0)) = 0.0_8
    bcl5zl(OPS_ACC5(0,0,0)) = bcl1zl(OPS_ACC1(0,0,0)) - bcl5zl(OPS_ACC5(0,0,0))
    bcl2zl(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqN_zl

SUBROUTINE bounds_kernel_eqO_xl(bclyxl)

    real(kind=8), dimension(1) ::  bclyxl

    bclyxl(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqO_xl

SUBROUTINE bounds_kernel_eqO_yl(bclyyl)

    real(kind=8), dimension(1) ::  bclyyl

    bclyyl(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqO_yl

SUBROUTINE bounds_kernel_eqO_zl(bclyzl)

    real(kind=8), dimension(1) ::  bclyzl

    bclyzl(OPS_ACC1(0,0,0)) = 0.0_8 

END SUBROUTINE bounds_kernel_eqO_zl

SUBROUTINE bounds_kernel_eqH_xr(bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr,strdxr,struxr,strvxr,strwxr, &
                                strpxr,ova2xr,acouxr,sorpxr,tt1xr,cobcxr,pinfxr,xgdlen,bet,flag_bet_xr,flag_pio_xr)

    real(kind=8), dimension(1) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr
    real(kind=8), dimension(1), intent(in) :: bcl5xr,strdxr,struxr,strvxr,strwxr
    real(kind=8), dimension(1), intent(in) :: strpxr,ova2xr,acouxr,sorpxr,tt1xr
    real(kind=8), intent(in) :: cobcxr,pinfxr,xgdlen,bet
    integer(kind=4), intent(in) :: flag_bet_xr,flag_pio_xr
    real(kind=8) :: bet_now

    bet_now = bet

    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC7(0,0,0)) * (bcl2xr(OPS_ACC2(0,0,0))-ova2xr(OPS_ACC11(0,0,0))*bcl5xr(OPS_ACC5(0,0,0)))
    bcl3xr(OPS_ACC3(0,0,0)) = struxr(OPS_ACC7(0,0,0)) * bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = struxr(OPS_ACC7(0,0,0)) * bcl4xr(OPS_ACC4(0,0,0))

!   OLD VALUE OF L1X
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_8*(struxr(OPS_ACC7(0,0,0))-acouxr(OPS_ACC12(0,0,0))) &
                            * (bcl5xr(OPS_ACC5(0,0,0)) - strdxr(OPS_ACC6(0,0,0))*acouxr(OPS_ACC12(0,0,0))*bcl1xr(OPS_ACC1(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L1X
!   TERMS INTRODUCED FOR LODATO'S BC- NC
    IF ( flag_bet_xr==1 ) THEN
        bet_now = struxr(OPS_ACC7(0,0,0))*struxr(OPS_ACC7(0,0,0)) + strvxr(OPS_ACC8(0,0,0))*strvxr(OPS_ACC8(0,0,0)) &
                + strwxr(OPS_ACC9(0,0,0))*strwxr(OPS_ACC9(0,0,0))
        bet_now = SQRT(bet_now)/acouxr(OPS_ACC12(0,0,0))
    END IF

    IF ( (struxr(OPS_ACC7(0,0,0)) < 0.0_8) .AND. (flag_pio_xr==1) ) THEN

        bcl2xr(OPS_ACC2(0,0,0)) = -ova2xr(OPS_ACC11(0,0,0))*sorpxr(OPS_ACC13(0,0,0)) - bcl2xr(OPS_ACC2(0,0,0))

        bcl3xr(OPS_ACC3(0,0,0)) = 0.1_8*(strvxr(OPS_ACC8(0,0,0))-0.0_8) - bcl3xr(OPS_ACC3(0,0,0))

        bcl4xr(OPS_ACC4(0,0,0)) = 0.1_8*(strwxr(OPS_ACC9(0,0,0))-0.0_8) - bcl4xr(OPS_ACC4(0,0,0))

        bcl1xr(OPS_ACC1(0,0,0)) = 0.5_8*sorpxr(OPS_ACC13(0,0,0)) &
                                + cobcxr*acouxr(OPS_ACC12(0,0,0))*(strpxr(OPS_ACC10(0,0,0))-pinfxr) &
                                + 0.5_8*(1.0_8-bet_now)*tt1xr(OPS_ACC14(0,0,0)) - bcl1xr(OPS_ACC1(0,0,0)) &
                                - 100.0_8*strdxr(OPS_ACC6(0,0,0))*(1.0_8/xgdlen) &
                                * (acouxr(OPS_ACC12(0,0,0))*acouxr(OPS_ACC12(0,0,0))-struxr(OPS_ACC7(0,0,0))*struxr(OPS_ACC7(0,0,0))) &
                                * (struxr(OPS_ACC7(0,0,0))-0.0_8)
    ELSE

        bcl1xr(OPS_ACC1(0,0,0)) = 0.5_8*sorpxr(OPS_ACC13(0,0,0)) &
                                + cobcxr*acouxr(OPS_ACC12(0,0,0))*(strpxr(OPS_ACC10(0,0,0))-pinfxr) &
                                + 0.5_8*(1.0_8-bet_now)*tt1xr(OPS_ACC14(0,0,0)) - bcl1xr(OPS_ACC1(0,0,0))
    END IF

END SUBROUTINE bounds_kernel_eqH_xr

SUBROUTINE bounds_kernel_eqH_yr(bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr,strdyr,struyr,strvyr,strwyr, &
                                strpyr,ova2yr,acouyr,sorpyr,tt1yr,cobcyr,pinfyr,ygdlen,bet,flag_bet_yr,flag_pio_yr)

    real(kind=8), dimension(1) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr
    real(kind=8), dimension(1), intent(in) :: bcl5yr,strdyr,struyr,strvyr,strwyr
    real(kind=8), dimension(1), intent(in) :: strpyr,ova2yr,acouyr,sorpyr,tt1yr
    real(kind=8), intent(in) :: cobcyr,pinfyr,ygdlen,bet
    integer(kind=4), intent(in) :: flag_bet_yr,flag_pio_yr
    real(kind=8) :: bet_now

    bet_now = bet

    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC8(0,0,0)) * (bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC5(0,0,0))*ova2yr(OPS_ACC11(0,0,0)))
    bcl3yr(OPS_ACC3(0,0,0)) = strvyr(OPS_ACC8(0,0,0)) * bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = strvyr(OPS_ACC8(0,0,0)) * bcl4yr(OPS_ACC4(0,0,0))

!   OLD VALUE OF L1Y
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_8*(strvyr(OPS_ACC8(0,0,0))-acouyr(OPS_ACC12(0,0,0))) &
                            * (bcl5yr(OPS_ACC5(0,0,0))-strdyr(OPS_ACC6(0,0,0))*acouyr(OPS_ACC12(0,0,0))*bcl1yr(OPS_ACC1(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L1Y
!   TERMS INTRODUCED FOR LODATO'S BC- NC
    IF( flag_bet_yr==1 ) THEN
        bet_now = struyr(OPS_ACC7(0,0,0))*struyr(OPS_ACC7(0,0,0)) + strvyr(OPS_ACC8(0,0,0))*strvyr(OPS_ACC8(0,0,0)) &
                + strwyr(OPS_ACC9(0,0,0))*strwyr(OPS_ACC9(0,0,0))
        bet_now = SQRT(bet_now)/acouyr(OPS_ACC12(0,0,0))
    END IF

    IF( (strvyr(OPS_ACC8(0,0,0)) < 0.0_8) .AND. (flag_pio_yr==1) ) THEN

        bcl2yr(OPS_ACC2(0,0,0)) = -ova2yr(OPS_ACC11(0,0,0))*sorpyr(OPS_ACC13(0,0,0)) - bcl2yr(OPS_ACC2(0,0,0))

        bcl3yr(OPS_ACC3(0,0,0)) = 0.1_8*(struyr(OPS_ACC7(0,0,0))-0.0_8) - bcl3yr(OPS_ACC3(0,0,0))

        bcl4yr(OPS_ACC4(0,0,0)) = 0.1_8*(strwyr(OPS_ACC9(0,0,0))-0.0_8) - bcl4yr(OPS_ACC4(0,0,0))

        bcl1yr(OPS_ACC1(0,0,0)) = 0.5_8*sorpyr(OPS_ACC13(0,0,0)) &
                                + cobcyr*acouyr(OPS_ACC12(0,0,0))*(strpyr(OPS_ACC10(0,0,0))-pinfyr) &
                                + 0.5_8*(1.0_8-bet_now)*tt1yr(OPS_ACC14(0,0,0)) - bcl1yr(OPS_ACC1(0,0,0)) &
                                - 100.0_8*strdyr(OPS_ACC6(0,0,0))*(1.0_8/ygdlen) &
                                * (acouyr(OPS_ACC12(0,0,0))*acouyr(OPS_ACC12(0,0,0))-strvyr(OPS_ACC8(0,0,0))*strvyr(OPS_ACC8(0,0,0))) &
                                * (strvyr(OPS_ACC8(0,0,0))-0.0_8)
    ELSE

        bcl1yr(OPS_ACC1(0,0,0)) = 0.5_8*sorpyr(OPS_ACC13(0,0,0)) &
                                + cobcyr*acouyr(OPS_ACC12(0,0,0))*(strpyr(OPS_ACC10(0,0,0))-pinfyr) &
                                + 0.5_8*(1.0_8-bet_now)*tt1yr(OPS_ACC14(0,0,0)) - bcl1yr(OPS_ACC1(0,0,0))
    END IF

END SUBROUTINE bounds_kernel_eqH_yr

SUBROUTINE bounds_kernel_eqH_zr(bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr,strdzr,struzr,strvzr,strwzr, &
                                strpzr,ova2zr,acouzr,sorpzr,tt1zr,cobczr,pinfzr,zgdlen,bet,flag_bet_zr,flag_pio_zr)

    real(kind=8), dimension(1) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr
    real(kind=8), dimension(1), intent(in) :: bcl5zr,strdzr,struzr,strvzr,strwzr
    real(kind=8), dimension(1), intent(in) :: strpzr,ova2zr,acouzr,sorpzr,tt1zr
    real(kind=8), intent(in) :: cobczr,pinfzr,zgdlen,bet
    integer(kind=4), intent(in) :: flag_bet_zr,flag_pio_zr
    real(kind=8) :: bet_now

    bet_now = bet

    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC9(0,0,0)) * (bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC5(0,0,0))*ova2zr(OPS_ACC11(0,0,0)))
    bcl3zr(OPS_ACC3(0,0,0)) = strwzr(OPS_ACC9(0,0,0)) * bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = strwzr(OPS_ACC9(0,0,0)) * bcl4zr(OPS_ACC4(0,0,0))

!   OLD VALUE OF L1Z
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_8*(strwzr(OPS_ACC9(0,0,0))-acouzr(OPS_ACC12(0,0,0))) &
                            * (bcl5zr(OPS_ACC5(0,0,0))-strdzr(OPS_ACC6(0,0,0))*acouzr(OPS_ACC12(0,0,0))*bcl1zr(OPS_ACC1(0,0,0)))

!   SUBTRACT FROM NEW VALUE OF L1Z
!   TERMS INTRODUCED FOR LODATO'S BC- NC
    IF ( flag_bet_zr==1 ) THEN
        bet_now = struzr(OPS_ACC7(0,0,0))*struzr(OPS_ACC7(0,0,0)) + strvzr(OPS_ACC8(0,0,0))*strvzr(OPS_ACC8(0,0,0)) &
                + strwzr(OPS_ACC9(0,0,0))*strwzr(OPS_ACC9(0,0,0))
        bet_now = SQRT(bet_now)/acouzr(OPS_ACC12(0,0,0))
    END IF

    IF ( (strwzr(OPS_ACC9(0,0,0)) < 0.0_8) .AND. (flag_pio_zr==1) ) THEN

        bcl2zr(OPS_ACC2(0,0,0)) = -ova2zr(OPS_ACC11(0,0,0))*sorpzr(OPS_ACC13(0,0,0)) - bcl2zr(OPS_ACC2(0,0,0))

        bcl3zr(OPS_ACC3(0,0,0)) = 0.1_8*(struzr(OPS_ACC7(0,0,0))-0.0_8) - bcl3zr(OPS_ACC3(0,0,0))

        bcl4zr(OPS_ACC4(0,0,0)) = 0.1_8*(strvzr(OPS_ACC8(0,0,0))-0.0_8) - bcl4zr(OPS_ACC4(0,0,0))

        bcl1zr(OPS_ACC1(0,0,0)) = 0.5_8*sorpzr(OPS_ACC13(0,0,0)) &
                                + cobczr*acouzr(OPS_ACC12(0,0,0))*(strpzr(OPS_ACC10(0,0,0))-pinfzr) &
                                + 0.5_8*(1.0_8-bet_now)*tt1zr(OPS_ACC14(0,0,0)) - bcl1zr(OPS_ACC1(0,0,0)) &
                                - 100.0_8*strdzr(OPS_ACC6(0,0,0))*(1.0_8/zgdlen) &
                                * (acouzr(OPS_ACC12(0,0,0))*acouzr(OPS_ACC12(0,0,0))-strwzr(OPS_ACC9(0,0,0))*strwzr(OPS_ACC9(0,0,0))) &
                                * (strwzr(OPS_ACC9(0,0,0))-0.0_8)
    ELSE

        bcl1zr(OPS_ACC1(0,0,0)) = 0.5_8*sorpzr(OPS_ACC13(0,0,0)) &
                                + cobczr*acouzr(OPS_ACC12(0,0,0))*(strpzr(OPS_ACC10(0,0,0))-pinfzr) &
                                + 0.5_8*(1.0_8-bet_now)*tt1zr(OPS_ACC14(0,0,0)) - bcl1zr(OPS_ACC1(0,0,0))
    END IF

END SUBROUTINE bounds_kernel_eqH_zr

SUBROUTINE bounds_kernel_eqI_xr(bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,acouxr,struxr,bcl5xr,ova2xr,sorpxr,strpxr,cobcxr,pinfxr)

    real(kind=8), dimension(1) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr
    real(kind=8), dimension(1), intent(in) :: strdxr,acouxr,struxr,bcl5xr,ova2xr,sorpxr,strpxr
    real(kind=8), intent(in) :: cobcxr,pinfxr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC5(0,0,0))*acouxr(OPS_ACC6(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_8*(struxr(OPS_ACC7(0,0,0))-acouxr(OPS_ACC6(0,0,0)))  &
            *(bcl5xr(OPS_ACC8(0,0,0))-fornow)
    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC7(0,0,0))  &
            *(bcl2xr(OPS_ACC2(0,0,0))-bcl5xr(OPS_ACC8(0,0,0))*ova2xr(OPS_ACC9(0,0,0)))
    bcl3xr(OPS_ACC3(0,0,0)) = struxr(OPS_ACC7(0,0,0))*bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = struxr(OPS_ACC7(0,0,0))*bcl4xr(OPS_ACC4(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2X-L4X)
!   L5X UNCHANGED
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_8*sorpxr(OPS_ACC10(0,0,0))  &
            + cobcxr*acouxr(OPS_ACC6(0,0,0))*(strpxr(OPS_ACC11(0,0,0))-pinfxr) - bcl1xr(OPS_ACC1(0,0,0))
    bcl2xr(OPS_ACC2(0,0,0)) = -bcl2xr(OPS_ACC2(0,0,0))
    bcl3xr(OPS_ACC3(0,0,0)) = -bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = -bcl4xr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqI_xr

SUBROUTINE bounds_kernel_eqI_yr(bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,acouyr,strvyr,bcl5yr,ova2yr,sorpyr,strpyr,cobcyr,pinfyr)

    real(kind=8), dimension(1) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr
    real(kind=8), dimension(1), intent(in) :: strdyr,acouyr,strvyr,bcl5yr,ova2yr,sorpyr,strpyr
    real(kind=8), intent(in) :: cobcyr,pinfyr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC5(0,0,0))*acouyr(OPS_ACC6(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_8*(strvyr(OPS_ACC7(0,0,0))-acouyr(OPS_ACC6(0,0,0)))  &
            *(bcl5yr(OPS_ACC8(0,0,0))-fornow)
    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC7(0,0,0))  &
            *(bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC8(0,0,0))*ova2yr(OPS_ACC9(0,0,0)))
    bcl3yr(OPS_ACC3(0,0,0)) = strvyr(OPS_ACC7(0,0,0))*bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = strvyr(OPS_ACC7(0,0,0))*bcl4yr(OPS_ACC4(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Y-L4Y)
!   L5Y UNCHANGED
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_8*sorpyr(OPS_ACC10(0,0,0))  &
            + cobcyr*acouyr(OPS_ACC6(0,0,0))*(strpyr(OPS_ACC11(0,0,0))-pinfyr) - bcl1yr(OPS_ACC1(0,0,0))
    bcl2yr(OPS_ACC2(0,0,0)) = -bcl2yr(OPS_ACC2(0,0,0))
    bcl3yr(OPS_ACC3(0,0,0)) = -bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = -bcl4yr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqI_yr

SUBROUTINE bounds_kernel_eqI_zr(bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,acouzr,strwzr,bcl5zr,ova2zr,sorpzr,strpzr,cobczr,pinfzr)

    real(kind=8), dimension(1) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr
    real(kind=8), dimension(1), intent(in) :: strdzr,acouzr,strwzr,bcl5zr,ova2zr,sorpzr,strpzr
    real(kind=8), intent(in) :: cobczr,pinfzr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC5(0,0,0))*acouzr(OPS_ACC6(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_8*(strwzr(OPS_ACC7(0,0,0))-acouzr(OPS_ACC6(0,0,0)))  &
            *(bcl5zr(OPS_ACC8(0,0,0))-fornow)
    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC7(0,0,0))  &
            *(bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC8(0,0,0))*ova2zr(OPS_ACC9(0,0,0)))
    bcl3zr(OPS_ACC3(0,0,0)) = strwzr(OPS_ACC7(0,0,0))*bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = strwzr(OPS_ACC7(0,0,0))*bcl4zr(OPS_ACC4(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR L2Z-L4Z)
!   L5Z UNCHANGED
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_8*sorpzr(OPS_ACC10(0,0,0))  &
            + cobczr*acouzr(OPS_ACC6(0,0,0))*(strpzr(OPS_ACC11(0,0,0))-pinfzr) - bcl1zr(OPS_ACC1(0,0,0))
    bcl2zr(OPS_ACC2(0,0,0)) = -bcl2zr(OPS_ACC2(0,0,0))
    bcl3zr(OPS_ACC3(0,0,0)) = -bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = -bcl4zr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqI_zr

SUBROUTINE bounds_kernel_eqJ_xr(bclyxr,ratexr,struxr,strdxr)

    real(kind=8), dimension(1) :: bclyxr
    real(kind=8), dimension(1), intent(in) :: ratexr,struxr,strdxr

!   OLD VALUE OF L's
    bclyxr(OPS_ACC1(0,0,0)) = struxr(OPS_ACC3(0,0,0))*bclyxr(OPS_ACC1(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYX)
    bclyxr(OPS_ACC1(0,0,0)) = ratexr(OPS_ACC2(0,0,0))/strdxr(OPS_ACC4(0,0,0))  &
              - bclyxr(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqJ_xr

SUBROUTINE bounds_kernel_eqJ_yr(bclyyr,rateyr,strvyr,strdyr)

    real(kind=8), dimension(1) :: bclyyr
    real(kind=8), dimension(1), intent(in) :: rateyr,strvyr,strdyr

!   OLD VALUE OF L's
    bclyyr(OPS_ACC1(0,0,0)) = strvyr(OPS_ACC3(0,0,0))*bclyyr(OPS_ACC1(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYY)
    bclyyr(OPS_ACC1(0,0,0)) = rateyr(OPS_ACC2(0,0,0))/strdyr(OPS_ACC4(0,0,0))  &
              - bclyyr(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqJ_yr

SUBROUTINE bounds_kernel_eqJ_zr(bclyzr,ratezr,strwzr,strdzr)

    real(kind=8), dimension(1) :: bclyzr
    real(kind=8), dimension(1), intent(in) :: ratezr,strwzr,strdzr

!   OLD VALUE OF L's
    bclyzr(OPS_ACC1(0,0,0)) = strwzr(OPS_ACC3(0,0,0))*bclyzr(OPS_ACC1(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's (=0 FOR LYZ)
    bclyzr(OPS_ACC1(0,0,0)) = ratezr(OPS_ACC2(0,0,0))/strdzr(OPS_ACC4(0,0,0))  &
              - bclyzr(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqJ_zr

SUBROUTINE bounds_kernel_eqK_xr(bcl1xr,bcl2xr,bcl5xr,strdxr,acouxr,struxr,ova2xr,dudtxr,gam1xr,dtdtxr,strtxr,sorpxr,strpxr,sydtxr)

    real(kind=8), dimension(1) :: bcl1xr,bcl2xr,bcl5xr
    real(kind=8), dimension(1), intent(in) :: strdxr,acouxr,struxr,ova2xr,dudtxr,gam1xr,dtdtxr,strtxr,sorpxr,strpxr,sydtxr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC4(0,0,0))*acouxr(OPS_ACC5(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_8*(struxr(OPS_ACC6(0,0,0))-acouxr(OPS_ACC5(0,0,0)))  &
            *(bcl5xr(OPS_ACC3(0,0,0))-fornow)
    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC6(0,0,0))  &
            *(bcl2xr(OPS_ACC2(0,0,0))-bcl5xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC7(0,0,0)))
    bcl5xr(OPS_ACC3(0,0,0)) = 0.5_8*(struxr(OPS_ACC6(0,0,0))+acouxr(OPS_ACC5(0,0,0)))  &
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

END SUBROUTINE bounds_kernel_eqK_xr

SUBROUTINE bounds_kernel_eqK_yr(bcl1yr,bcl2yr,bcl5yr,strdyr,acouyr,strvyr,ova2yr,dvdtyr,gam1yr,dtdtyr,strtyr,sorpyr,strpyr,sydtyr)

    real(kind=8), dimension(1) :: bcl1yr,bcl2yr,bcl5yr
    real(kind=8), dimension(1), intent(in) :: strdyr,acouyr,strvyr,ova2yr,dvdtyr,gam1yr,dtdtyr,strtyr,sorpyr,strpyr,sydtyr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC4(0,0,0))*acouyr(OPS_ACC5(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_8*(strvyr(OPS_ACC6(0,0,0))-acouyr(OPS_ACC5(0,0,0)))  &
            *(bcl5yr(OPS_ACC3(0,0,0))-fornow)
    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC6(0,0,0))  &
            *(bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC7(0,0,0)))
    bcl5yr(OPS_ACC3(0,0,0)) = 0.5_8*(strvyr(OPS_ACC6(0,0,0))+acouyr(OPS_ACC5(0,0,0)))  &
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

END SUBROUTINE bounds_kernel_eqK_yr

SUBROUTINE bounds_kernel_eqK_zr(bcl1zr,bcl2zr,bcl5zr,strdzr,acouzr,strwzr,ova2zr,dwdtzr,gam1zr,dtdtzr,strtzr,sorpzr,strpzr,sydtzr)

    real(kind=8), dimension(1) :: bcl1zr,bcl2zr,bcl5zr
    real(kind=8), dimension(1), intent(in) :: strdzr,acouzr,strwzr,ova2zr,dwdtzr,gam1zr,dtdtzr,strtzr,sorpzr,strpzr,sydtzr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC4(0,0,0))*acouzr(OPS_ACC5(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_8*(strwzr(OPS_ACC6(0,0,0))-acouzr(OPS_ACC5(0,0,0)))  &
            *(bcl5zr(OPS_ACC3(0,0,0))-fornow)
    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC6(0,0,0))  &
            *(bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC7(0,0,0)))
    bcl5zr(OPS_ACC3(0,0,0)) = 0.5_8*(strwzr(OPS_ACC6(0,0,0))+acouzr(OPS_ACC5(0,0,0)))  &
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

END SUBROUTINE bounds_kernel_eqK_zr

SUBROUTINE bounds_kernel_eqL_xr(bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr,strdxr,acouxr,struxr,ova2xr,dudtxr,dddtxr,dvdtxr,dwdtxr)

    real(kind=8), dimension(1) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr
    real(kind=8), dimension(1), intent(in) :: strdxr,acouxr,struxr,ova2xr,dudtxr,dddtxr,dvdtxr,dwdtxr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC6(0,0,0))*acouxr(OPS_ACC7(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_8*(struxr(OPS_ACC8(0,0,0))-acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))-fornow)
    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC8(0,0,0))  &
            *(bcl2xr(OPS_ACC2(0,0,0))-bcl5xr(OPS_ACC5(0,0,0))*ova2xr(OPS_ACC9(0,0,0)))
    bcl3xr(OPS_ACC3(0,0,0)) = struxr(OPS_ACC8(0,0,0))*bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = struxr(OPS_ACC8(0,0,0))*bcl4xr(OPS_ACC4(0,0,0))
    bcl5xr(OPS_ACC5(0,0,0)) = 0.5_8*(struxr(OPS_ACC8(0,0,0))+acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5X UNCHANGED
    fornow = bcl5xr(OPS_ACC5(0,0,0)) + strdxr(OPS_ACC6(0,0,0))*acouxr(OPS_ACC7(0,0,0))*dudtxr(OPS_ACC10(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = fornow - bcl1xr(OPS_ACC1(0,0,0))
    bcl2xr(OPS_ACC2(0,0,0)) = -dddtxr(OPS_ACC11(0,0,0))  &
            - ova2xr(OPS_ACC9(0,0,0))*(bcl5xr(OPS_ACC5(0,0,0))+fornow) - bcl2xr(OPS_ACC2(0,0,0))
    bcl3xr(OPS_ACC3(0,0,0)) = -dvdtxr(OPS_ACC12(0,0,0)) - bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = -dwdtxr(OPS_ACC13(0,0,0)) - bcl4xr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqL_xr

SUBROUTINE bounds_kernel_eqL_yr(bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr,strdyr,acouyr,strvyr,ova2yr,dudtyr,dddtyr,dvdtyr,dwdtyr)

    real(kind=8), dimension(1) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr
    real(kind=8), dimension(1), intent(in) :: strdyr,acouyr,strvyr,ova2yr,dudtyr,dddtyr,dvdtyr,dwdtyr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC6(0,0,0))*acouyr(OPS_ACC7(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_8*(strvyr(OPS_ACC8(0,0,0))-acouyr(OPS_ACC7(0,0,0)))  &
            *(bcl5yr(OPS_ACC5(0,0,0))-fornow)
    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC8(0,0,0))  &
            *(bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC5(0,0,0))*ova2yr(OPS_ACC9(0,0,0)))
    bcl3yr(OPS_ACC3(0,0,0)) = strvyr(OPS_ACC8(0,0,0))*bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = strvyr(OPS_ACC8(0,0,0))*bcl4yr(OPS_ACC4(0,0,0))
    bcl5yr(OPS_ACC5(0,0,0)) = 0.5_8*(strvyr(OPS_ACC8(0,0,0))+acouyr(OPS_ACC7(0,0,0)))  &
            *(bcl5yr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Y UNCHANGED
    fornow = bcl5yr(OPS_ACC5(0,0,0)) + strdyr(OPS_ACC6(0,0,0))*acouyr(OPS_ACC7(0,0,0))*dvdtyr(OPS_ACC10(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = fornow - bcl1yr(OPS_ACC1(0,0,0))
    bcl2yr(OPS_ACC2(0,0,0)) = -dddtyr(OPS_ACC11(0,0,0))  &
            - ova2yr(OPS_ACC9(0,0,0))*(bcl5yr(OPS_ACC5(0,0,0))+fornow) - bcl2yr(OPS_ACC2(0,0,0))
    bcl3yr(OPS_ACC3(0,0,0)) = -dudtyr(OPS_ACC12(0,0,0)) - bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = -dwdtyr(OPS_ACC13(0,0,0)) - bcl4yr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqL_yr

SUBROUTINE bounds_kernel_eqL_zr(bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr,strdzr,acouzr,strwzr,ova2zr,dudtzr,dddtzr,dvdtzr,dwdtzr)

    real(kind=8), dimension(1) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr
    real(kind=8), dimension(1), intent(in) :: strdzr,acouzr,strwzr,ova2zr,dudtzr,dddtzr,dvdtzr,dwdtzr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC6(0,0,0))*acouzr(OPS_ACC7(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_8*(strwzr(OPS_ACC8(0,0,0))-acouzr(OPS_ACC7(0,0,0)))  &
            *(bcl5zr(OPS_ACC5(0,0,0))-fornow)
    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC8(0,0,0))  &
            *(bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC5(0,0,0))*ova2zr(OPS_ACC9(0,0,0)))
    bcl3zr(OPS_ACC3(0,0,0)) = strwzr(OPS_ACC8(0,0,0))*bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = strwzr(OPS_ACC8(0,0,0))*bcl4zr(OPS_ACC4(0,0,0))
    bcl5zr(OPS_ACC5(0,0,0)) = 0.5_8*(strwzr(OPS_ACC8(0,0,0))+acouzr(OPS_ACC7(0,0,0)))  &
            *(bcl5zr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Z UNCHANGED
    fornow = bcl5zr(OPS_ACC5(0,0,0)) + strdzr(OPS_ACC6(0,0,0))*acouzr(OPS_ACC7(0,0,0))*dwdtzr(OPS_ACC10(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = fornow - bcl1zr(OPS_ACC1(0,0,0))
    bcl2zr(OPS_ACC2(0,0,0)) = -dddtzr(OPS_ACC11(0,0,0))  &
            - ova2zr(OPS_ACC9(0,0,0))*(bcl5zr(OPS_ACC5(0,0,0))+fornow) - bcl2zr(OPS_ACC2(0,0,0))
    bcl3zr(OPS_ACC3(0,0,0)) = -dudtzr(OPS_ACC12(0,0,0)) - bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = -dvdtzr(OPS_ACC13(0,0,0)) - bcl4zr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqL_zr

SUBROUTINE bounds_kernel_eqM_xr(bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr,strdxr,acouxr,struxr,ova2xr)

    real(kind=8), dimension(1) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr
    real(kind=8), dimension(1), intent(in) :: strdxr,acouxr,struxr,ova2xr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC6(0,0,0))*acouxr(OPS_ACC7(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_8*(struxr(OPS_ACC8(0,0,0))-acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))-fornow)
    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC8(0,0,0))  &
            *(bcl2xr(OPS_ACC2(0,0,0))-bcl5xr(OPS_ACC5(0,0,0))*ova2xr(OPS_ACC9(0,0,0)))
    bcl3xr(OPS_ACC3(0,0,0)) = struxr(OPS_ACC8(0,0,0)) * bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = struxr(OPS_ACC8(0,0,0)) * bcl4xr(OPS_ACC4(0,0,0))
    bcl5xr(OPS_ACC5(0,0,0)) = 0.5_8*(struxr(OPS_ACC8(0,0,0))+acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5X UNCHANGED
!   NC SUGGESTION: IF NOT WORKING REMOVE SECOND LINE AS WELL
    bcl1xr(OPS_ACC1(0,0,0)) = bcl5xr(OPS_ACC5(0,0,0)) - bcl1xr(OPS_ACC1(0,0,0))
    bcl3xr(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4xr(OPS_ACC4(0,0,0)) = 0.0_8
    bcl2xr(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqM_xr

SUBROUTINE bounds_kernel_eqM_yr(bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr,strdyr,acouyr,strvyr,ova2yr)

    real(kind=8), dimension(1) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr
    real(kind=8), dimension(1), intent(in) :: strdyr,acouyr,strvyr,ova2yr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC6(0,0,0))*acouyr(OPS_ACC7(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_8*(strvyr(OPS_ACC8(0,0,0))-acouyr(OPS_ACC7(0,0,0)))  &
                            * (bcl5yr(OPS_ACC5(0,0,0))-fornow)
    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC8(0,0,0))  &
                            * (bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC5(0,0,0))*ova2yr(OPS_ACC9(0,0,0)))
    bcl3yr(OPS_ACC3(0,0,0)) = strvyr(OPS_ACC8(0,0,0)) * bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = strvyr(OPS_ACC8(0,0,0)) * bcl4yr(OPS_ACC4(0,0,0))
    bcl5yr(OPS_ACC5(0,0,0)) = 0.5_8*(strvyr(OPS_ACC8(0,0,0))+acouyr(OPS_ACC7(0,0,0)))  &
                            * (bcl5yr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Y UNCHANGED
    bcl1yr(OPS_ACC1(0,0,0)) = bcl5yr(OPS_ACC5(0,0,0)) - bcl1yr(OPS_ACC1(0,0,0))
    bcl3yr(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4yr(OPS_ACC4(0,0,0)) = 0.0_8
    bcl2yr(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqM_yr

SUBROUTINE bounds_kernel_eqM_zr(bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr,strdzr,acouzr,strwzr,ova2zr)

    real(kind=8), dimension(1) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr
    real(kind=8), dimension(1), intent(in) :: strdzr,acouzr,strwzr,ova2zr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC6(0,0,0))*acouzr(OPS_ACC7(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_8*(strwzr(OPS_ACC8(0,0,0))-acouzr(OPS_ACC7(0,0,0)))  &
                            * (bcl5zr(OPS_ACC5(0,0,0))-fornow)
    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC8(0,0,0))  &
                            * (bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC5(0,0,0))*ova2zr(OPS_ACC9(0,0,0)))
    bcl3zr(OPS_ACC3(0,0,0)) = strwzr(OPS_ACC8(0,0,0)) * bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = strwzr(OPS_ACC8(0,0,0)) * bcl4zr(OPS_ACC4(0,0,0))
    bcl5zr(OPS_ACC5(0,0,0)) = 0.5_8*(strwzr(OPS_ACC8(0,0,0))+acouzr(OPS_ACC7(0,0,0)))  &
                            * (bcl5zr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Z UNCHANGED
    bcl1zr(OPS_ACC1(0,0,0)) = bcl5zr(OPS_ACC5(0,0,0)) - bcl1zr(OPS_ACC1(0,0,0))
    bcl3zr(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4zr(OPS_ACC4(0,0,0)) = 0.0_8
    bcl2zr(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqM_zr

SUBROUTINE bounds_kernel_eqN_xr(bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr,strdxr,acouxr,struxr,ova2xr)

    real(kind=8), dimension(1) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr
    real(kind=8), dimension(1), intent(in) :: strdxr,acouxr,struxr,ova2xr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdxr(OPS_ACC6(0,0,0))*acouxr(OPS_ACC7(0,0,0))*bcl1xr(OPS_ACC1(0,0,0))
    bcl1xr(OPS_ACC1(0,0,0)) = 0.5_8*(struxr(OPS_ACC8(0,0,0))-acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))-fornow)
    bcl2xr(OPS_ACC2(0,0,0)) = struxr(OPS_ACC8(0,0,0))  &
            *(bcl2xr(OPS_ACC2(0,0,0))-bcl5xr(OPS_ACC5(0,0,0))*ova2xr(OPS_ACC9(0,0,0)))
    bcl3xr(OPS_ACC3(0,0,0)) = struxr(OPS_ACC8(0,0,0))*bcl3xr(OPS_ACC3(0,0,0))
    bcl4xr(OPS_ACC4(0,0,0)) = struxr(OPS_ACC8(0,0,0))*bcl4xr(OPS_ACC4(0,0,0))
    bcl5xr(OPS_ACC5(0,0,0)) = 0.5_8*(struxr(OPS_ACC8(0,0,0))+acouxr(OPS_ACC7(0,0,0)))  &
            *(bcl5xr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5X UNCHANGED
!   NC SUGGESTION: IF NOT WORKING REMOVE SECOND LINE AS WELL
    bcl1xr(OPS_ACC1(0,0,0)) = bcl5xr(OPS_ACC5(0,0,0)) - bcl1xr(OPS_ACC1(0,0,0))
    bcl3xr(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4xr(OPS_ACC4(0,0,0)) = 0.0_8
    bcl2xr(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqN_xr

SUBROUTINE bounds_kernel_eqN_yr(bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr,strdyr,acouyr,strvyr,ova2yr)

    real(kind=8), dimension(1) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr
    real(kind=8), dimension(1), intent(in) :: strdyr,acouyr,strvyr,ova2yr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdyr(OPS_ACC6(0,0,0))*acouyr(OPS_ACC7(0,0,0))*bcl1yr(OPS_ACC1(0,0,0))
    bcl1yr(OPS_ACC1(0,0,0)) = 0.5_8*(strvyr(OPS_ACC8(0,0,0))-acouyr(OPS_ACC7(0,0,0)))  &
            *(bcl5yr(OPS_ACC5(0,0,0))-fornow)
    bcl2yr(OPS_ACC2(0,0,0)) = strvyr(OPS_ACC8(0,0,0))  &
            *(bcl2yr(OPS_ACC2(0,0,0))-bcl5yr(OPS_ACC5(0,0,0))*ova2yr(OPS_ACC9(0,0,0)))
    bcl3yr(OPS_ACC3(0,0,0)) = strvyr(OPS_ACC8(0,0,0))*bcl3yr(OPS_ACC3(0,0,0))
    bcl4yr(OPS_ACC4(0,0,0)) = strvyr(OPS_ACC8(0,0,0))*bcl4yr(OPS_ACC4(0,0,0))
    bcl5yr(OPS_ACC5(0,0,0)) = 0.5_8*(strvyr(OPS_ACC8(0,0,0))+acouyr(OPS_ACC7(0,0,0)))  &
            *(bcl5yr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Y UNCHANGED
    bcl1yr(OPS_ACC1(0,0,0)) = bcl5yr(OPS_ACC5(0,0,0)) - bcl1yr(OPS_ACC1(0,0,0))
    bcl3yr(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4yr(OPS_ACC4(0,0,0)) = 0.0_8
    bcl2yr(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqN_yr

SUBROUTINE bounds_kernel_eqN_zr(bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr,strdzr,acouzr,strwzr,ova2zr)

    real(kind=8), dimension(1) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr
    real(kind=8), dimension(1), intent(in) :: strdzr,acouzr,strwzr,ova2zr
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdzr(OPS_ACC6(0,0,0))*acouzr(OPS_ACC7(0,0,0))*bcl1zr(OPS_ACC1(0,0,0))
    bcl1zr(OPS_ACC1(0,0,0)) = 0.5_8*(strwzr(OPS_ACC8(0,0,0))-acouzr(OPS_ACC7(0,0,0)))  &
            *(bcl5zr(OPS_ACC5(0,0,0))-fornow)
    bcl2zr(OPS_ACC2(0,0,0)) = strwzr(OPS_ACC8(0,0,0))  &
            *(bcl2zr(OPS_ACC2(0,0,0))-bcl5zr(OPS_ACC5(0,0,0))*ova2zr(OPS_ACC9(0,0,0)))
    bcl3zr(OPS_ACC3(0,0,0)) = strwzr(OPS_ACC8(0,0,0))*bcl3zr(OPS_ACC3(0,0,0))
    bcl4zr(OPS_ACC4(0,0,0)) = strwzr(OPS_ACC8(0,0,0))*bcl4zr(OPS_ACC4(0,0,0))
    bcl5zr(OPS_ACC5(0,0,0)) = 0.5_8*(strwzr(OPS_ACC8(0,0,0))+acouzr(OPS_ACC7(0,0,0)))  &
            *(bcl5zr(OPS_ACC5(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L5Z UNCHANGED
    bcl1zr(OPS_ACC1(0,0,0)) = bcl5zr(OPS_ACC5(0,0,0)) - bcl1zr(OPS_ACC1(0,0,0))
    bcl3zr(OPS_ACC3(0,0,0)) = 0.0_8
    bcl4zr(OPS_ACC4(0,0,0)) = 0.0_8
    bcl2zr(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqN_zr

SUBROUTINE bounds_kernel_eqO_xr(bclyxr)

    real(kind=8), dimension(1) :: bclyxr

    bclyxr(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqO_xr

SUBROUTINE bounds_kernel_eqO_yr(bclyyr)

    real(kind=8), dimension(1) :: bclyyr

    bclyyr(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqO_yr

SUBROUTINE bounds_kernel_eqO_zr(bclyzr)

    real(kind=8), dimension(1) :: bclyzr

    bclyzr(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqO_zr

SUBROUTINE bounds_kernel_eqP_xl(drhs,urhs,vrhs,wrhs,erhs,bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,struxl, &
                                strvxl,strwxl,strexl,ova2xl,ovgmxl,acouxl,flag_pio_xl)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,struxl
    real(kind=8), dimension(1), intent(in) :: strvxl,strwxl,strexl,ova2xl,ovgmxl,acouxl
    integer(kind=4), intent(in) :: flag_pio_xl

    IF ( (struxl(OPS_ACC11(0,0,0)) > 0.0_8) .AND. (flag_pio_xl==1) ) THEN

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC15(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))*struxl(OPS_ACC11(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC15(0,0,0))*(struxl(OPS_ACC11(0,0,0))+acouxl(OPS_ACC17(0,0,0)))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))*strvxl(OPS_ACC12(0,0,0)) &
                              - bcl3xl(OPS_ACC7(0,0,0))*strdxl(OPS_ACC10(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC15(0,0,0))*strvxl(OPS_ACC12(0,0,0))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))*strwxl(OPS_ACC13(0,0,0)) &
                              - bcl4xl(OPS_ACC8(0,0,0))*strdxl(OPS_ACC10(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC15(0,0,0))*strwxl(OPS_ACC13(0,0,0))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))*strexl(OPS_ACC14(0,0,0)) &
                              - bcl3xl(OPS_ACC7(0,0,0))*strdxl(OPS_ACC10(0,0,0))*strvxl(OPS_ACC12(0,0,0)) &
                              - bcl4xl(OPS_ACC8(0,0,0))*strdxl(OPS_ACC10(0,0,0))*strwxl(OPS_ACC13(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*(ova2xl(OPS_ACC15(0,0,0))*strexl(OPS_ACC14(0,0,0)) &
                              + struxl(OPS_ACC11(0,0,0))/acouxl(OPS_ACC17(0,0,0)) + ovgmxl(OPS_ACC16(0,0,0)))

    ELSE

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC15(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC15(0,0,0))*(struxl(OPS_ACC11(0,0,0))+acouxl(OPS_ACC17(0,0,0)))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC15(0,0,0))*strvxl(OPS_ACC12(0,0,0))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC15(0,0,0))*strwxl(OPS_ACC13(0,0,0))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*(ova2xl(OPS_ACC15(0,0,0))*strexl(OPS_ACC14(0,0,0)) &
                              + struxl(OPS_ACC11(0,0,0))/acouxl(OPS_ACC17(0,0,0)) + ovgmxl(OPS_ACC16(0,0,0)))

    END IF

END SUBROUTINE bounds_kernel_eqP_xl

SUBROUTINE bounds_kernel_eqP_yl(drhs,urhs,vrhs,wrhs,erhs,bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,struyl, &
                                strvyl,strwyl,streyl,ova2yl,ovgmyl,acouyl,flag_pio_yl)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,struyl
    real(kind=8), dimension(1), intent(in) :: strvyl,strwyl,streyl,ova2yl,ovgmyl,acouyl
    integer(kind=4), intent(in) :: flag_pio_yl

    IF ( (strvyl(OPS_ACC12(0,0,0)) > 0.0_8) .AND. (flag_pio_yl==1) ) THEN

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC15(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))*struyl(OPS_ACC11(0,0,0)) &
                              - bcl3yl(OPS_ACC7(0,0,0))*strdyl(OPS_ACC10(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC15(0,0,0))*struyl(OPS_ACC11(0,0,0))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))*strvyl(OPS_ACC12(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC15(0,0,0))*(strvyl(OPS_ACC12(0,0,0))+acouyl(OPS_ACC17(0,0,0)))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))*strwyl(OPS_ACC13(0,0,0)) &
                              - bcl4yl(OPS_ACC8(0,0,0))*strdyl(OPS_ACC10(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC15(0,0,0))*strwyl(OPS_ACC13(0,0,0))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))*streyl(OPS_ACC14(0,0,0)) &
                              - bcl3yl(OPS_ACC7(0,0,0))*strdyl(OPS_ACC10(0,0,0))*struyl(OPS_ACC11(0,0,0)) &
                              - bcl4yl(OPS_ACC8(0,0,0))*strdyl(OPS_ACC10(0,0,0))*strwyl(OPS_ACC13(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*(ova2yl(OPS_ACC15(0,0,0))*streyl(OPS_ACC14(0,0,0)) &
                              + strvyl(OPS_ACC12(0,0,0))/acouyl(OPS_ACC17(0,0,0)) + ovgmyl(OPS_ACC16(0,0,0)))

    ELSE

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC15(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC15(0,0,0))*struyl(OPS_ACC11(0,0,0))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC15(0,0,0))*(strvyl(OPS_ACC12(0,0,0))+acouyl(OPS_ACC17(0,0,0)))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC15(0,0,0))*strwyl(OPS_ACC13(0,0,0))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*(ova2yl(OPS_ACC15(0,0,0))*streyl(OPS_ACC14(0,0,0)) &
                              + strvyl(OPS_ACC12(0,0,0))/acouyl(OPS_ACC17(0,0,0)) + ovgmyl(OPS_ACC16(0,0,0)))

    END IF

END SUBROUTINE bounds_kernel_eqP_yl

SUBROUTINE bounds_kernel_eqP_zl(drhs,urhs,vrhs,wrhs,erhs,bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,struzl, &
                                strvzl,strwzl,strezl,ova2zl,ovgmzl,acouzl,flag_pio_zl)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,struzl
    real(kind=8), dimension(1), intent(in) :: strvzl,strwzl,strezl,ova2zl,ovgmzl,acouzl
    integer(kind=4), intent(in) :: flag_pio_zl

    IF( (strwzl(OPS_ACC13(0,0,0)) > 0.0_8) .AND. (flag_pio_zl==1) ) THEN

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))  &
                              - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC15(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))*struzl(OPS_ACC11(0,0,0)) &
                              - bcl3zl(OPS_ACC7(0,0,0))*strdzl(OPS_ACC10(0,0,0)) &
                              - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC15(0,0,0))*struzl(OPS_ACC11(0,0,0))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))*strvzl(OPS_ACC12(0,0,0)) &
                              - bcl4zl(OPS_ACC8(0,0,0))*strdzl(OPS_ACC10(0,0,0)) &
                              - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC15(0,0,0))*strvzl(OPS_ACC12(0,0,0))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))*strwzl(OPS_ACC13(0,0,0)) &
                              - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC15(0,0,0))*(strwzl(OPS_ACC13(0,0,0))+acouzl(OPS_ACC17(0,0,0)))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))*strezl(OPS_ACC14(0,0,0)) &
                              - bcl3zl(OPS_ACC7(0,0,0))*strdzl(OPS_ACC10(0,0,0))*struzl(OPS_ACC11(0,0,0)) &
                              - bcl4zl(OPS_ACC8(0,0,0))*strdzl(OPS_ACC10(0,0,0))*strvzl(OPS_ACC12(0,0,0)) &
                              - bcl5zl(OPS_ACC9(0,0,0))*(ova2zl(OPS_ACC15(0,0,0))*strezl(OPS_ACC14(0,0,0)) &
                              + strwzl(OPS_ACC13(0,0,0))/acouzl(OPS_ACC17(0,0,0)) + ovgmzl(OPS_ACC16(0,0,0)))

    ELSE

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC15(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) &
                              - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC15(0,0,0))*struzl(OPS_ACC11(0,0,0))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) &
                              - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC15(0,0,0))*strvzl(OPS_ACC12(0,0,0))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) &
                              - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC15(0,0,0))*(strwzl(OPS_ACC13(0,0,0))+acouzl(OPS_ACC17(0,0,0)))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) &
                              - bcl5zl(OPS_ACC9(0,0,0))*(ova2zl(OPS_ACC15(0,0,0))*strezl(OPS_ACC14(0,0,0)) &
                              + strwzl(OPS_ACC13(0,0,0))/acouzl(OPS_ACC17(0,0,0)) + ovgmzl(OPS_ACC16(0,0,0)))

    END IF

END SUBROUTINE bounds_kernel_eqP_zl

SUBROUTINE bounds_kernel_eqQ_xl(erhs,yrhs,strdxl,struxl,strhxl,stryxl,bclyxl,bcl2xl,bcl5xl,ova2xl,flag_pio_xl)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: strdxl,struxl,strhxl,stryxl,bclyxl,bcl2xl,bcl5xl,ova2xl
    integer(kind=4), intent(in) :: flag_pio_xl

    real(kind=8) :: fornow

    IF ( (struxl(OPS_ACC4(0,0,0)) > 0.0_8) .AND. (flag_pio_xl==1) ) THEN

        fornow = bclyxl(OPS_ACC7(0,0,0))*strdxl(OPS_ACC3(0,0,0))

        erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhxl(OPS_ACC5(0,0,0))

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - (bcl2xl(OPS_ACC8(0,0,0))+bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC10(0,0,0)))*stryxl(OPS_ACC6(0,0,0)) &
                              - fornow

    ELSE

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - bcl5xl(OPS_ACC9(0,0,0))*ova2xl(OPS_ACC10(0,0,0))*stryxl(OPS_ACC6(0,0,0))

    END IF

END SUBROUTINE bounds_kernel_eqQ_xl

SUBROUTINE bounds_kernel_eqQ_yl(erhs,yrhs,strdyl,strvyl,strhyl,stryyl,bclyyl,bcl2yl,bcl5yl,ova2yl,flag_pio_yl)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: strdyl,strvyl,strhyl,stryyl,bclyyl,bcl2yl,bcl5yl,ova2yl
    integer(kind=4), intent(in) :: flag_pio_yl

    real(kind=8) :: fornow

    IF ( (strvyl(OPS_ACC4(0,0,0)) > 0.0_8) .AND. (flag_pio_yl==1) ) THEN

        fornow = bclyyl(OPS_ACC7(0,0,0))*strdyl(OPS_ACC3(0,0,0))

        erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhyl(OPS_ACC5(0,0,0))

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - (bcl2yl(OPS_ACC8(0,0,0))+bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC10(0,0,0)))*stryyl(OPS_ACC6(0,0,0)) &
                              - fornow

    ELSE

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - bcl5yl(OPS_ACC9(0,0,0))*ova2yl(OPS_ACC10(0,0,0))*stryyl(OPS_ACC6(0,0,0))

    END IF

END SUBROUTINE bounds_kernel_eqQ_yl

SUBROUTINE bounds_kernel_eqQ_zl(erhs,yrhs,strdzl,strwzl,strhzl,stryzl,bclyzl,bcl2zl,bcl5zl,ova2zl,flag_pio_zl)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: strdzl,strwzl,strhzl,stryzl,bclyzl,bcl2zl,bcl5zl,ova2zl
    integer(kind=4), intent(in) :: flag_pio_zl

    real(kind=8) :: fornow

    IF ( (strwzl(OPS_ACC4(0,0,0)) > 0.0_8) .AND. (flag_pio_zl==1) ) THEN

        fornow = bclyzl(OPS_ACC7(0,0,0))*strdzl(OPS_ACC3(0,0,0))

        erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhzl(OPS_ACC5(0,0,0))

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - (bcl2zl(OPS_ACC8(0,0,0))+bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC10(0,0,0)))*stryzl(OPS_ACC6(0,0,0)) &
                              - fornow

    ELSE

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - bcl5zl(OPS_ACC9(0,0,0))*ova2zl(OPS_ACC10(0,0,0))*stryzl(OPS_ACC6(0,0,0))

    END IF

END SUBROUTINE bounds_kernel_eqQ_zl

SUBROUTINE bounds_kernel_eqR_xl(drhs,urhs,vrhs,wrhs,erhs,bcl2xl,bcl3xl,bcl4xl,bcl5xl,struxl,strvxl,strwxl,strdxl,strexl,acouxl,ova2xl,ovgmxl)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2xl,bcl3xl,bcl4xl,bcl5xl,struxl,strvxl,strwxl,strdxl,strexl,acouxl,ova2xl,ovgmxl

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

END SUBROUTINE bounds_kernel_eqR_xl

SUBROUTINE bounds_kernel_eqR_yl(drhs,urhs,vrhs,wrhs,erhs,bcl2yl,bcl3yl,bcl4yl,bcl5yl,struyl,strvyl,strwyl,strdyl,streyl,acouyl,ova2yl,ovgmyl)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2yl,bcl3yl,bcl4yl,bcl5yl,struyl,strvyl,strwyl,strdyl,streyl,acouyl,ova2yl,ovgmyl

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

END SUBROUTINE bounds_kernel_eqR_yl

SUBROUTINE bounds_kernel_eqR_zl(drhs,urhs,vrhs,wrhs,erhs,bcl2zl,bcl3zl,bcl4zl,bcl5zl,struzl,strvzl,strwzl,strdzl,strezl,acouzl,ova2zl,ovgmzl)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2zl,bcl3zl,bcl4zl,bcl5zl,struzl,strvzl,strwzl,strdzl,strezl,acouzl,ova2zl,ovgmzl

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

END SUBROUTINE bounds_kernel_eqR_zl

SUBROUTINE bounds_kernel_eqS_xl(erhs,yrhs,bclyxl,strhxl,stryxl,strdxl,bcl2xl,bcl5xl,ova2xl)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: bclyxl,strhxl,stryxl,strdxl,bcl2xl,bcl5xl,ova2xl
    real(kind=8) :: fornow

    fornow = bclyxl(OPS_ACC3(0,0,0))*strdxl(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhxl(OPS_ACC4(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0))  &
              - (bcl2xl(OPS_ACC7(0,0,0))+bcl5xl(OPS_ACC8(0,0,0))*ova2xl(OPS_ACC9(0,0,0)))*stryxl(OPS_ACC5(0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_eqS_xl

SUBROUTINE bounds_kernel_eqS_yl(erhs,yrhs,bclyyl,strhyl,stryyl,strdyl,bcl2yl,bcl5yl,ova2yl)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: bclyyl,strhyl,stryyl,strdyl,bcl2yl,bcl5yl,ova2yl
    real(kind=8) :: fornow

    fornow = bclyyl(OPS_ACC3(0,0,0))*strdyl(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhyl(OPS_ACC4(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0))  &
              - (bcl2yl(OPS_ACC7(0,0,0))+bcl5yl(OPS_ACC8(0,0,0))*ova2yl(OPS_ACC9(0,0,0)))*stryyl(OPS_ACC5(0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_eqS_yl

SUBROUTINE bounds_kernel_eqS_zl(erhs,yrhs,bclyzl,strhzl,stryzl,strdzl,bcl2zl,bcl5zl,ova2zl)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: bclyzl,strhzl,stryzl,strdzl,bcl2zl,bcl5zl,ova2zl
    real(kind=8) :: fornow

    fornow = bclyzl(OPS_ACC3(0,0,0))*strdzl(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhzl(OPS_ACC4(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0))  &
              - (bcl2zl(OPS_ACC7(0,0,0))+bcl5zl(OPS_ACC8(0,0,0))*ova2zl(OPS_ACC9(0,0,0)))*stryzl(OPS_ACC5(0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_eqS_zl

SUBROUTINE bounds_kernel_eqT_xl(drhs,bcl2xl,bcl5xl,ova2xl)

    real(kind=8), dimension(1) :: drhs
    real(kind=8), dimension(1), intent(in) :: bcl2xl,bcl5xl,ova2xl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2xl(OPS_ACC2(0,0,0))  &
            - bcl5xl(OPS_ACC3(0,0,0))*ova2xl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqT_xl

SUBROUTINE bounds_kernel_eqT_yl(drhs,bcl2yl,bcl5yl,ova2yl)

    real(kind=8), dimension(1) :: drhs
    real(kind=8), dimension(1), intent(in) :: bcl2yl,bcl5yl,ova2yl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2yl(OPS_ACC2(0,0,0))  &
            - bcl5yl(OPS_ACC3(0,0,0))*ova2yl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqT_yl

SUBROUTINE bounds_kernel_eqT_zl(drhs,bcl2zl,bcl5zl,ova2zl)

    real(kind=8), dimension(1) :: drhs
    real(kind=8), dimension(1), intent(in) :: bcl2zl,bcl5zl,ova2zl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2zl(OPS_ACC2(0,0,0))  &
            - bcl5zl(OPS_ACC3(0,0,0))*ova2zl(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqT_zl

SUBROUTINE bounds_kernel_eqU_xl(erhs,bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,strexl,struxl,strvxl,strwxl,ova2xl,ovgmxl,acouxl)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: bcl2xl,bcl3xl,bcl4xl,bcl5xl,strdxl,strexl,struxl,strvxl,strwxl,ova2xl,ovgmxl,acouxl

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - bcl2xl(OPS_ACC2(0,0,0))*strexl(OPS_ACC7(0,0,0))  &
            - bcl3xl(OPS_ACC3(0,0,0))*strdxl(OPS_ACC6(0,0,0))*strvxl(OPS_ACC9(0,0,0))  &
            - bcl4xl(OPS_ACC4(0,0,0))*strdxl(OPS_ACC6(0,0,0))*strwxl(OPS_ACC10(0,0,0))  &
            - bcl5xl(OPS_ACC5(0,0,0))*(ova2xl(OPS_ACC11(0,0,0))*strexl(OPS_ACC7(0,0,0))  &
            + struxl(OPS_ACC8(0,0,0))/acouxl(OPS_ACC13(0,0,0)) + ovgmxl(OPS_ACC12(0,0,0)))

END SUBROUTINE bounds_kernel_eqU_xl

SUBROUTINE bounds_kernel_eqU_yl(erhs,bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,streyl,struyl,strvyl,strwyl,ova2yl,ovgmyl,acouyl)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: bcl2yl,bcl3yl,bcl4yl,bcl5yl,strdyl,streyl,struyl,strvyl,strwyl,ova2yl,ovgmyl,acouyl

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - bcl2yl(OPS_ACC2(0,0,0))*streyl(OPS_ACC7(0,0,0))  &
            - bcl3yl(OPS_ACC3(0,0,0))*strdyl(OPS_ACC6(0,0,0))*struyl(OPS_ACC9(0,0,0))  &
            - bcl4yl(OPS_ACC4(0,0,0))*strdyl(OPS_ACC6(0,0,0))*strwyl(OPS_ACC10(0,0,0))  &
            - bcl5yl(OPS_ACC5(0,0,0))*(ova2yl(OPS_ACC11(0,0,0))*streyl(OPS_ACC7(0,0,0))  &
            + strvyl(OPS_ACC8(0,0,0))/acouyl(OPS_ACC13(0,0,0)) + ovgmyl(OPS_ACC12(0,0,0)))

END SUBROUTINE bounds_kernel_eqU_yl

SUBROUTINE bounds_kernel_eqU_zl(erhs,bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,strezl,struzl,strvzl,strwzl,ova2zl,ovgmzl,acouzl)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: bcl2zl,bcl3zl,bcl4zl,bcl5zl,strdzl,strezl,struzl,strvzl,strwzl,ova2zl,ovgmzl,acouzl

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - bcl2zl(OPS_ACC2(0,0,0))*strezl(OPS_ACC7(0,0,0))  &
            - bcl3zl(OPS_ACC3(0,0,0))*strdzl(OPS_ACC6(0,0,0))*struzl(OPS_ACC9(0,0,0))  &
            - bcl4zl(OPS_ACC4(0,0,0))*strdzl(OPS_ACC6(0,0,0))*strvzl(OPS_ACC10(0,0,0))  &
            - bcl5zl(OPS_ACC5(0,0,0))*(ova2zl(OPS_ACC11(0,0,0))*strezl(OPS_ACC7(0,0,0))  &
            + strwzl(OPS_ACC8(0,0,0))/acouzl(OPS_ACC13(0,0,0)) + ovgmzl(OPS_ACC12(0,0,0)))

END SUBROUTINE bounds_kernel_eqU_zl

SUBROUTINE bounds_kernel_eqV_xl(erhs,bclyxl,ratexl,dydtxl,strhxl,strdxl,struxl)

    real(kind=8), dimension(1) :: erhs,bclyxl
    real(kind=8), dimension(1), intent(in) :: ratexl,dydtxl,strhxl,strdxl,struxl

    bclyxl(OPS_ACC2(0,0,0)) = ratexl(OPS_ACC3(0,0,0))/strdxl(OPS_ACC6(0,0,0))  &
              - dydtxl(OPS_ACC4(0,0,0)) - struxl(OPS_ACC7(0,0,0))*bclyxl(OPS_ACC2(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyxl(OPS_ACC2(0,0,0))*strdxl(OPS_ACC6(0,0,0))*strhxl(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqV_xl

SUBROUTINE bounds_kernel_eqV_yl(erhs,bclyyl,rateyl,dydtyl,strhyl,strdyl,strvyl)

    real(kind=8), dimension(1) :: erhs,bclyyl
    real(kind=8), dimension(1), intent(in) :: rateyl,dydtyl,strhyl,strdyl,strvyl

    bclyyl(OPS_ACC2(0,0,0)) = rateyl(OPS_ACC3(0,0,0))/strdyl(OPS_ACC6(0,0,0))  &
              - dydtyl(OPS_ACC4(0,0,0)) - strvyl(OPS_ACC7(0,0,0))*bclyyl(OPS_ACC2(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyyl(OPS_ACC2(0,0,0))*strdyl(OPS_ACC6(0,0,0))*strhyl(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqV_yl

SUBROUTINE bounds_kernel_eqV_zl(erhs,bclyzl,ratezl,dydtzl,strhzl,strdzl,strwzl)

    real(kind=8), dimension(1) :: erhs,bclyzl
    real(kind=8), dimension(1), intent(in) :: ratezl,dydtzl,strhzl,strdzl,strwzl

    bclyzl(OPS_ACC2(0,0,0)) = ratezl(OPS_ACC3(0,0,0))/strdzl(OPS_ACC6(0,0,0))  &
              - dydtzl(OPS_ACC4(0,0,0)) - strwzl(OPS_ACC7(0,0,0))*bclyzl(OPS_ACC2(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyzl(OPS_ACC2(0,0,0))*strdzl(OPS_ACC6(0,0,0))*strhzl(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqV_zl

SUBROUTINE bounds_kernel_eqW_xl(drhs,urhs,vrhs,wrhs,erhs,bcl2xl,bcl5xl,ova2xl)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2xl,bcl5xl,ova2xl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2xl(OPS_ACC6(0,0,0))  &
                          - bcl5xl(OPS_ACC7(0,0,0))*ova2xl(OPS_ACC8(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = 0.0_8
    vrhs(OPS_ACC3(0,0,0)) = 0.0_8
    wrhs(OPS_ACC4(0,0,0)) = 0.0_8
    erhs(OPS_ACC5(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqW_xl

SUBROUTINE bounds_kernel_eqW_yl(drhs,urhs,vrhs,wrhs,erhs,bcl2yl,bcl5yl,ova2yl)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2yl,bcl5yl,ova2yl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2yl(OPS_ACC6(0,0,0))  &
                          - bcl5yl(OPS_ACC7(0,0,0))*ova2yl(OPS_ACC8(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = 0.0_8
    vrhs(OPS_ACC3(0,0,0)) = 0.0_8
    wrhs(OPS_ACC4(0,0,0)) = 0.0_8
    erhs(OPS_ACC5(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqW_yl

SUBROUTINE bounds_kernel_eqW_zl(drhs,urhs,vrhs,wrhs,erhs,bcl2zl,bcl5zl,ova2zl)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2zl,bcl5zl,ova2zl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2zl(OPS_ACC6(0,0,0))  &
                          - bcl5zl(OPS_ACC7(0,0,0))*ova2zl(OPS_ACC8(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = 0.0_8
    vrhs(OPS_ACC3(0,0,0)) = 0.0_8
    wrhs(OPS_ACC4(0,0,0)) = 0.0_8
    erhs(OPS_ACC5(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqW_zl

SUBROUTINE bounds_kernel_eqX_xl(yrhs,stryxl,bcl2xl,bcl5xl,ova2xl)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryxl,bcl2xl,bcl5xl,ova2xl

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
                          - (bcl2xl(OPS_ACC3(0,0,0)) + bcl5xl(OPS_ACC4(0,0,0))*ova2xl(OPS_ACC5(0,0,0))) &
                            * stryxl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqX_xl

SUBROUTINE bounds_kernel_eqX_yl(yrhs,stryyl,bcl2yl,bcl5yl,ova2yl)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryyl,bcl2yl,bcl5yl,ova2yl

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
                          - (bcl2yl(OPS_ACC3(0,0,0)) + bcl5yl(OPS_ACC4(0,0,0))*ova2yl(OPS_ACC5(0,0,0))) &
                            * stryyl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqX_yl

SUBROUTINE bounds_kernel_eqX_zl(yrhs,stryzl,bcl2zl,bcl5zl,ova2zl)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryzl,bcl2zl,bcl5zl,ova2zl

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
                          - (bcl2zl(OPS_ACC3(0,0,0)) + bcl5zl(OPS_ACC4(0,0,0))*ova2zl(OPS_ACC5(0,0,0))) &
                            * stryzl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqX_zl

SUBROUTINE bounds_kernel_eqY_xl(drhs,bcl2xl,bcl5xl,ova2xl,urhs,vrhs,wrhs,erhs)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2xl,bcl5xl,ova2xl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2xl(OPS_ACC2(0,0,0))  &
            - bcl5xl(OPS_ACC3(0,0,0))*ova2xl(OPS_ACC4(0,0,0))
    urhs(OPS_ACC5(0,0,0)) = 0.0_8
    vrhs(OPS_ACC6(0,0,0)) = 0.0_8
    wrhs(OPS_ACC7(0,0,0)) = 0.0_8
    erhs(OPS_ACC8(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqY_xl

SUBROUTINE bounds_kernel_eqY_yl(drhs,bcl2yl,bcl5yl,ova2yl,urhs,vrhs,wrhs,erhs)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2yl,bcl5yl,ova2yl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2yl(OPS_ACC2(0,0,0))  &
            - bcl5yl(OPS_ACC3(0,0,0))*ova2yl(OPS_ACC4(0,0,0))
    urhs(OPS_ACC5(0,0,0)) = 0.0_8
    vrhs(OPS_ACC6(0,0,0)) = 0.0_8
    wrhs(OPS_ACC7(0,0,0)) = 0.0_8
    erhs(OPS_ACC8(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqY_yl

SUBROUTINE bounds_kernel_eqY_zl(drhs,bcl2zl,bcl5zl,ova2zl,urhs,vrhs,wrhs,erhs)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2zl,bcl5zl,ova2zl

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2zl(OPS_ACC2(0,0,0))  &
            - bcl5zl(OPS_ACC3(0,0,0))*ova2zl(OPS_ACC4(0,0,0))
    urhs(OPS_ACC5(0,0,0)) = 0.0_8
    vrhs(OPS_ACC6(0,0,0)) = 0.0_8
    wrhs(OPS_ACC7(0,0,0)) = 0.0_8
    erhs(OPS_ACC8(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqY_zl

SUBROUTINE bounds_kernel_eqZ_xl(yrhs,stryxl,bcl2xl,bcl5xl,ova2xl)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryxl,bcl2xl,bcl5xl,ova2xl

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
              - (bcl2xl(OPS_ACC3(0,0,0))+bcl5xl(OPS_ACC4(0,0,0))*ova2xl(OPS_ACC5(0,0,0)))*stryxl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqZ_xl

SUBROUTINE bounds_kernel_eqZ_yl(yrhs,stryyl,bcl2yl,bcl5yl,ova2yl)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryyl,bcl2yl,bcl5yl,ova2yl

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
              - (bcl2yl(OPS_ACC3(0,0,0))+bcl5yl(OPS_ACC4(0,0,0))*ova2yl(OPS_ACC5(0,0,0)))*stryyl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqZ_yl

SUBROUTINE bounds_kernel_eqZ_zl(yrhs,stryzl,bcl2zl,bcl5zl,ova2zl)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryzl,bcl2zl,bcl5zl,ova2zl

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
              - (bcl2zl(OPS_ACC3(0,0,0))+bcl5zl(OPS_ACC4(0,0,0))*ova2zl(OPS_ACC5(0,0,0)))*stryzl(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqZ_zl

SUBROUTINE bounds_kernel_eqP_xr(drhs,urhs,vrhs,wrhs,erhs,bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,struxr, &
                                strvxr,strwxr,strexr,ova2xr,ovgmxr,acouxr,flag_pio_xr)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,struxr
    real(kind=8), dimension(1), intent(in) :: strvxr,strwxr,strexr,ova2xr,ovgmxr,acouxr
    integer(kind=4), intent(in) :: flag_pio_xr

    IF ( (struxr(OPS_ACC11(0,0,0)) < 0.0_8) .AND. (flag_pio_xr==1) ) THEN

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC15(0,0,0)) &
                              - bcl2xr(OPS_ACC7(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) &
                              - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC15(0,0,0))*(struxr(OPS_ACC11(0,0,0))-acouxr(OPS_ACC17(0,0,0))) &
                              - bcl2xr(OPS_ACC7(0,0,0))*struxr(OPS_ACC11(0,0,0))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) &
                              - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC15(0,0,0))*strvxr(OPS_ACC12(0,0,0)) &
                              - bcl2xr(OPS_ACC7(0,0,0))*strvxr(OPS_ACC12(0,0,0)) &
                              - bcl3xr(OPS_ACC8(0,0,0))*strdxr(OPS_ACC10(0,0,0))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) &
                              - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC15(0,0,0))*strwxr(OPS_ACC13(0,0,0)) &
                              - bcl2xr(OPS_ACC7(0,0,0))*strwxr(OPS_ACC13(0,0,0)) &
                              - bcl4xr(OPS_ACC9(0,0,0))*strdxr(OPS_ACC10(0,0,0))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) &
                              - bcl1xr(OPS_ACC6(0,0,0))*(ova2xr(OPS_ACC15(0,0,0))*strexr(OPS_ACC14(0,0,0)) &
                              - struxr(OPS_ACC11(0,0,0))/acouxr(OPS_ACC17(0,0,0)) + ovgmxr(OPS_ACC16(0,0,0))) &
                              - bcl2xr(OPS_ACC7(0,0,0))*strexr(OPS_ACC14(0,0,0)) &
                              - bcl3xr(OPS_ACC8(0,0,0))*strdxr(OPS_ACC10(0,0,0))*strvxr(OPS_ACC12(0,0,0)) &
                              - bcl4xr(OPS_ACC9(0,0,0))*strdxr(OPS_ACC10(0,0,0))*strwxr(OPS_ACC13(0,0,0))

    ELSE

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC15(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) &
                              - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC15(0,0,0))*(struxr(OPS_ACC11(0,0,0))-acouxr(OPS_ACC17(0,0,0)))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) &
                              - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC15(0,0,0))*strvxr(OPS_ACC12(0,0,0))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) &
                              - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC15(0,0,0))*strwxr(OPS_ACC13(0,0,0))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) &
                              - bcl1xr(OPS_ACC6(0,0,0))*(ova2xr(OPS_ACC15(0,0,0))*strexr(OPS_ACC14(0,0,0)) &
                              - struxr(OPS_ACC11(0,0,0))/acouxr(OPS_ACC17(0,0,0)) + ovgmxr(OPS_ACC16(0,0,0)))

    END IF

END SUBROUTINE bounds_kernel_eqP_xr

SUBROUTINE bounds_kernel_eqP_yr(drhs,urhs,vrhs,wrhs,erhs,bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,struyr, &
                                strvyr,strwyr,streyr,ova2yr,ovgmyr,acouyr,flag_pio_yr)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,struyr
    real(kind=8), dimension(1), intent(in) :: strvyr,strwyr,streyr,ova2yr,ovgmyr,acouyr
    integer(kind=4), intent(in) :: flag_pio_yr

    IF( (strvyr(OPS_ACC12(0,0,0)) < 0.0_8) .AND. (flag_pio_yr==1) ) THEN

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC15(0,0,0)) &
                              - bcl2yr(OPS_ACC7(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) &
                              - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC15(0,0,0))*struyr(OPS_ACC11(0,0,0)) &
                              - bcl2yr(OPS_ACC7(0,0,0))*struyr(OPS_ACC11(0,0,0)) &
                              - bcl3yr(OPS_ACC8(0,0,0))*strdyr(OPS_ACC10(0,0,0))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) &
                              - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC15(0,0,0))*(strvyr(OPS_ACC12(0,0,0))-acouyr(OPS_ACC17(0,0,0))) &
                              - bcl2yr(OPS_ACC7(0,0,0))*strvyr(OPS_ACC12(0,0,0))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) &
                              - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC15(0,0,0))*strwyr(OPS_ACC13(0,0,0)) &
                              - bcl2yr(OPS_ACC7(0,0,0))*strwyr(OPS_ACC13(0,0,0)) &
                              - bcl4yr(OPS_ACC9(0,0,0))*strdyr(OPS_ACC10(0,0,0))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) &
                              - bcl1yr(OPS_ACC6(0,0,0))*(ova2yr(OPS_ACC15(0,0,0))*streyr(OPS_ACC14(0,0,0)) &
                              - strvyr(OPS_ACC12(0,0,0))/acouyr(OPS_ACC17(0,0,0)) + ovgmyr(OPS_ACC16(0,0,0))) &
                              - bcl2yr(OPS_ACC7(0,0,0))*streyr(OPS_ACC14(0,0,0)) &
                              - bcl3yr(OPS_ACC8(0,0,0))*strdyr(OPS_ACC10(0,0,0))*struyr(OPS_ACC11(0,0,0)) &
                              - bcl4yr(OPS_ACC9(0,0,0))*strdyr(OPS_ACC10(0,0,0))*strwyr(OPS_ACC13(0,0,0))

    ELSE

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC15(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) &
                              - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC15(0,0,0))*struyr(OPS_ACC11(0,0,0))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) &
                              - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC15(0,0,0))*(strvyr(OPS_ACC12(0,0,0))-acouyr(OPS_ACC17(0,0,0)))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) &
                              - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC15(0,0,0))*strwyr(OPS_ACC13(0,0,0))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) &
                              - bcl1yr(OPS_ACC6(0,0,0))*(ova2yr(OPS_ACC15(0,0,0))*streyr(OPS_ACC14(0,0,0)) &
                              - strvyr(OPS_ACC12(0,0,0))/acouyr(OPS_ACC17(0,0,0)) + ovgmyr(OPS_ACC16(0,0,0)))

    END IF

END SUBROUTINE bounds_kernel_eqP_yr

SUBROUTINE bounds_kernel_eqP_zr(drhs,urhs,vrhs,wrhs,erhs,bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,struzr, &
                                strvzr,strwzr,strezr,ova2zr,ovgmzr,acouzr,flag_pio_zr)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,struzr
    real(kind=8), dimension(1), intent(in) :: strvzr,strwzr,strezr,ova2zr,ovgmzr,acouzr
    integer(kind=4), intent(in) :: flag_pio_zr

    IF( (strwzr(OPS_ACC13(0,0,0)) < 0.0_8) .AND. (flag_pio_zr==1) ) THEN

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC15(0,0,0)) &
                              - bcl2zr(OPS_ACC7(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) &
                              - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC15(0,0,0))*struzr(OPS_ACC11(0,0,0)) &
                              - bcl2zr(OPS_ACC7(0,0,0))*struzr(OPS_ACC11(0,0,0)) &
                              - bcl3zr(OPS_ACC8(0,0,0))*strdzr(OPS_ACC10(0,0,0))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) &
                              - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC15(0,0,0))*strvzr(OPS_ACC12(0,0,0)) &
                              - bcl2zr(OPS_ACC7(0,0,0))*strvzr(OPS_ACC12(0,0,0)) &
                              - bcl4zr(OPS_ACC9(0,0,0))*strdzr(OPS_ACC10(0,0,0))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) &
                              - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC15(0,0,0))*(strwzr(OPS_ACC13(0,0,0))-acouzr(OPS_ACC17(0,0,0))) &
                              - bcl2zr(OPS_ACC7(0,0,0))*strwzr(OPS_ACC13(0,0,0))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) &
                              - bcl1zr(OPS_ACC6(0,0,0))*(ova2zr(OPS_ACC15(0,0,0))*strezr(OPS_ACC14(0,0,0)) &
                              - strwzr(OPS_ACC13(0,0,0))/acouzr(OPS_ACC17(0,0,0)) + ovgmzr(OPS_ACC16(0,0,0))) &
                              - bcl2zr(OPS_ACC7(0,0,0))*strezr(OPS_ACC14(0,0,0)) &
                              - bcl3zr(OPS_ACC8(0,0,0))*strdzr(OPS_ACC10(0,0,0))*struzr(OPS_ACC11(0,0,0)) &
                              - bcl4zr(OPS_ACC9(0,0,0))*strdzr(OPS_ACC10(0,0,0))*strvzr(OPS_ACC12(0,0,0))

    ELSE

        drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC15(0,0,0))

        urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) &
                              - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC15(0,0,0))*struzr(OPS_ACC11(0,0,0))

        vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) &
                              - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC15(0,0,0))*strvzr(OPS_ACC12(0,0,0))

        wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) &
                              - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC15(0,0,0))*(strwzr(OPS_ACC13(0,0,0))-acouzr(OPS_ACC17(0,0,0)))

        erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) &
                              - bcl1zr(OPS_ACC6(0,0,0))*(ova2zr(OPS_ACC15(0,0,0))*strezr(OPS_ACC14(0,0,0))  &
                              - strwzr(OPS_ACC13(0,0,0))/acouzr(OPS_ACC17(0,0,0)) + ovgmzr(OPS_ACC16(0,0,0)))

    END IF

END SUBROUTINE bounds_kernel_eqP_zr

SUBROUTINE bounds_kernel_eqQ_xr(erhs,yrhs,strdxr,struxr,strhxr,stryxr,bclyxr,bcl1xr,bcl2xr,ova2xr,flag_pio_xr)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: strdxr,struxr,strhxr,stryxr,bclyxr,bcl1xr,bcl2xr,ova2xr
    integer(kind=4), intent(in) :: flag_pio_xr

    real(kind=8) :: fornow

    IF ( (struxr(OPS_ACC4(0,0,0)) < 0.0_8) .AND. (flag_pio_xr==1) ) THEN

        fornow = bclyxr(OPS_ACC7(0,0,0))*strdxr(OPS_ACC3(0,0,0))

        erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhxr(OPS_ACC5(0,0,0))

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - (bcl2xr(OPS_ACC9(0,0,0))+bcl1xr(OPS_ACC8(0,0,0))*ova2xr(OPS_ACC10(0,0,0)))*stryxr(OPS_ACC6(0,0,0)) &
                              - fornow

    ELSE

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - bcl1xr(OPS_ACC8(0,0,0))*ova2xr(OPS_ACC10(0,0,0))*stryxr(OPS_ACC6(0,0,0))

    END IF

END SUBROUTINE bounds_kernel_eqQ_xr

SUBROUTINE bounds_kernel_eqQ_yr(erhs,yrhs,strdyr,strvyr,strhyr,stryyr,bclyyr,bcl1yr,bcl2yr,ova2yr,flag_pio_yr)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: strdyr,strvyr,strhyr,stryyr,bclyyr,bcl1yr,bcl2yr,ova2yr
    integer(kind=4), intent(in) :: flag_pio_yr

    real(kind=8) :: fornow

    IF ( (strvyr(OPS_ACC4(0,0,0)) < 0.0_8) .AND. (flag_pio_yr==1) ) THEN

        fornow = bclyyr(OPS_ACC7(0,0,0))*strdyr(OPS_ACC3(0,0,0))

        erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhyr(OPS_ACC5(0,0,0))

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - (bcl2yr(OPS_ACC9(0,0,0))+bcl1yr(OPS_ACC8(0,0,0))*ova2yr(OPS_ACC10(0,0,0)))*stryyr(OPS_ACC6(0,0,0)) &
                              - fornow

    ELSE

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - bcl1yr(OPS_ACC8(0,0,0))*ova2yr(OPS_ACC10(0,0,0))*stryyr(OPS_ACC6(0,0,0))

    END IF

END SUBROUTINE bounds_kernel_eqQ_yr

SUBROUTINE bounds_kernel_eqQ_zr(erhs,yrhs,strdzr,strwzr,strhzr,stryzr,bclyzr,bcl1zr,bcl2zr,ova2zr,flag_pio_zr)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: strdzr,strwzr,strhzr,stryzr,bclyzr,bcl1zr,bcl2zr,ova2zr
    integer(kind=4), intent(in) :: flag_pio_zr

    real(kind=8) :: fornow

    IF ( (strwzr(OPS_ACC4(0,0,0)) < 0.0_8) .AND. (flag_pio_zr==1) ) THEN

        fornow = bclyzr(OPS_ACC7(0,0,0))*strdzr(OPS_ACC3(0,0,0))

        erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhzr(OPS_ACC5(0,0,0))

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - (bcl2zr(OPS_ACC9(0,0,0))+bcl1zr(OPS_ACC8(0,0,0))*ova2zr(OPS_ACC10(0,0,0)))*stryzr(OPS_ACC6(0,0,0)) &
                              - fornow

    ELSE

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - bcl1zr(OPS_ACC8(0,0,0))*ova2zr(OPS_ACC10(0,0,0))*stryzr(OPS_ACC6(0,0,0))

    END IF

END SUBROUTINE bounds_kernel_eqQ_zr

SUBROUTINE bounds_kernel_eqR_xr(drhs,urhs,vrhs,wrhs,erhs,bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,strexr,acouxr,struxr,strvxr,strwxr,ova2xr,ovgmxr)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,strexr,acouxr,struxr,strvxr,strwxr,ova2xr,ovgmxr

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

END SUBROUTINE bounds_kernel_eqR_xr

SUBROUTINE bounds_kernel_eqR_yr(drhs,urhs,vrhs,wrhs,erhs,bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,streyr,acouyr,struyr,strvyr,strwyr,ova2yr,ovgmyr)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,streyr,acouyr,struyr,strvyr,strwyr,ova2yr,ovgmyr

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

END SUBROUTINE bounds_kernel_eqR_yr

SUBROUTINE bounds_kernel_eqR_zr(drhs,urhs,vrhs,wrhs,erhs,bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,strezr,acouzr,struzr,strvzr,strwzr,ova2zr,ovgmzr)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,strezr,acouzr,struzr,strvzr,strwzr,ova2zr,ovgmzr

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

END SUBROUTINE bounds_kernel_eqR_zr

SUBROUTINE bounds_kernel_eqS_xr(erhs,yrhs,bclyxr,strhxr,stryxr,strdxr,bcl2xr,bcl1xr,ova2xr)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: bclyxr,strhxr,stryxr,strdxr,bcl2xr,bcl1xr,ova2xr
    real(kind=8) :: fornow

    fornow = bclyxr(OPS_ACC3(0,0,0))*strdxr(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhxr(OPS_ACC4(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0))  &
              - (bcl2xr(OPS_ACC7(0,0,0))+bcl1xr(OPS_ACC8(0,0,0))*ova2xr(OPS_ACC9(0,0,0)))*stryxr(OPS_ACC5(0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_eqS_xr

SUBROUTINE bounds_kernel_eqS_yr(erhs,yrhs,bclyyr,strhyr,stryyr,strdyr,bcl2yr,bcl1yr,ova2yr)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: bclyyr,strhyr,stryyr,strdyr,bcl2yr,bcl1yr,ova2yr
    real(kind=8) :: fornow

    fornow = bclyyr(OPS_ACC3(0,0,0))*strdyr(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhyr(OPS_ACC4(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0))  &
              - (bcl2yr(OPS_ACC7(0,0,0))+bcl1yr(OPS_ACC8(0,0,0))*ova2yr(OPS_ACC9(0,0,0)))*stryyr(OPS_ACC5(0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_eqS_yr

SUBROUTINE bounds_kernel_eqS_zr(erhs,yrhs,bclyzr,strhzr,stryzr,strdzr,bcl2zr,bcl5zr,ova2zr)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: bclyzr,strhzr,stryzr,strdzr,bcl2zr,bcl5zr,ova2zr
    real(kind=8) :: fornow

    fornow = bclyzr(OPS_ACC3(0,0,0))*strdzr(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhzr(OPS_ACC4(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0))  &
              - (bcl2zr(OPS_ACC7(0,0,0))+bcl5zr(OPS_ACC8(0,0,0))*ova2zr(OPS_ACC9(0,0,0)))*stryzr(OPS_ACC5(0,0,0))  &
              - fornow

END SUBROUTINE bounds_kernel_eqS_zr

SUBROUTINE bounds_kernel_eqT_xr(drhs,bcl2xr,bcl1xr,ova2xr)

    real(kind=8), dimension(1) :: drhs
    real(kind=8), dimension(1), intent(in) :: bcl2xr,bcl1xr,ova2xr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2xr(OPS_ACC2(0,0,0))  &
            - bcl1xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqT_xr

SUBROUTINE bounds_kernel_eqT_yr(drhs,bcl2yr,bcl1yr,ova2yr)

    real(kind=8), dimension(1) :: drhs
    real(kind=8), dimension(1), intent(in) :: bcl2yr,bcl1yr,ova2yr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2yr(OPS_ACC2(0,0,0))  &
            - bcl1yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqT_yr

SUBROUTINE bounds_kernel_eqT_zr(drhs,bcl2zr,bcl1zr,ova2zr)

    real(kind=8), dimension(1) :: drhs
    real(kind=8), dimension(1), intent(in) :: bcl2zr,bcl1zr,ova2zr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2zr(OPS_ACC2(0,0,0))  &
            - bcl1zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqT_zr

SUBROUTINE bounds_kernel_eqU_xr(erhs,bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,strexr,struxr,strvxr,strwxr,ova2xr,ovgmxr,acouxr)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: bcl1xr,bcl2xr,bcl3xr,bcl4xr,strdxr,strexr,struxr,strvxr,strwxr,ova2xr,ovgmxr,acouxr

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
            - bcl1xr(OPS_ACC2(0,0,0))*(ova2xr(OPS_ACC11(0,0,0))*strexr(OPS_ACC7(0,0,0))  &
            - struxr(OPS_ACC8(0,0,0))/acouxr(OPS_ACC13(0,0,0)) + ovgmxr(OPS_ACC12(0,0,0)))  &
            - bcl2xr(OPS_ACC3(0,0,0))*strexr(OPS_ACC7(0,0,0))  &
            - bcl3xr(OPS_ACC4(0,0,0))*strdxr(OPS_ACC6(0,0,0))*strvxr(OPS_ACC9(0,0,0))  &
            - bcl4xr(OPS_ACC5(0,0,0))*strdxr(OPS_ACC6(0,0,0))*strwxr(OPS_ACC10(0,0,0))

END SUBROUTINE bounds_kernel_eqU_xr

SUBROUTINE bounds_kernel_eqU_yr(erhs,bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,streyr,struyr,strvyr,strwyr,ova2yr,ovgmyr,acouyr)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: bcl1yr,bcl2yr,bcl3yr,bcl4yr,strdyr,streyr,struyr,strvyr,strwyr,ova2yr,ovgmyr,acouyr

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
            - bcl1yr(OPS_ACC2(0,0,0))*(ova2yr(OPS_ACC11(0,0,0))*streyr(OPS_ACC7(0,0,0))  &
            + strvyr(OPS_ACC8(0,0,0))/acouyr(OPS_ACC13(0,0,0)) + ovgmyr(OPS_ACC12(0,0,0)))  &
            - bcl2yr(OPS_ACC3(0,0,0))*streyr(OPS_ACC7(0,0,0))  &
            - bcl3yr(OPS_ACC4(0,0,0))*strdyr(OPS_ACC6(0,0,0))*struyr(OPS_ACC9(0,0,0))  &
            - bcl4yr(OPS_ACC5(0,0,0))*strdyr(OPS_ACC6(0,0,0))*strwyr(OPS_ACC10(0,0,0))

END SUBROUTINE bounds_kernel_eqU_yr

SUBROUTINE bounds_kernel_eqU_zr(erhs,bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,strezr,struzr,strvzr,strwzr,ova2zr,ovgmzr,acouzr)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: bcl1zr,bcl2zr,bcl3zr,bcl4zr,strdzr,strezr,struzr,strvzr,strwzr,ova2zr,ovgmzr,acouzr

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
            - bcl1zr(OPS_ACC2(0,0,0))*(ova2zr(OPS_ACC11(0,0,0))*strezr(OPS_ACC7(0,0,0))  &
            + strwzr(OPS_ACC8(0,0,0))/acouzr(OPS_ACC13(0,0,0)) + ovgmzr(OPS_ACC12(0,0,0)))  &
            - bcl2zr(OPS_ACC3(0,0,0))*strezr(OPS_ACC7(0,0,0))  &
            - bcl3zr(OPS_ACC4(0,0,0))*strdzr(OPS_ACC6(0,0,0))*struzr(OPS_ACC9(0,0,0))  &
            - bcl4zr(OPS_ACC5(0,0,0))*strdzr(OPS_ACC6(0,0,0))*strvzr(OPS_ACC10(0,0,0))

END SUBROUTINE bounds_kernel_eqU_zr

SUBROUTINE bounds_kernel_eqV_xr(erhs,bclyxr,ratexr,dydtxr,strhxr,strdxr,struxr)

    real(kind=8), dimension(1) :: erhs,bclyxr
    real(kind=8), dimension(1), intent(in) :: ratexr,dydtxr,strhxr,strdxr,struxr

    bclyxr(OPS_ACC2(0,0,0)) = ratexr(OPS_ACC3(0,0,0))/strdxr(OPS_ACC6(0,0,0))  &
              - dydtxr(OPS_ACC4(0,0,0)) - struxr(OPS_ACC7(0,0,0))*bclyxr(OPS_ACC2(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyxr(OPS_ACC2(0,0,0))*strdxr(OPS_ACC6(0,0,0))*strhxr(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqV_xr

SUBROUTINE bounds_kernel_eqV_yr(erhs,bclyyr,rateyr,dydtyr,strhyr,strdyr,strvyr)

    real(kind=8), dimension(1) :: erhs,bclyyr
    real(kind=8), dimension(1), intent(in) :: rateyr,dydtyr,strhyr,strdyr,strvyr

    bclyyr(OPS_ACC2(0,0,0)) = rateyr(OPS_ACC3(0,0,0))/strdyr(OPS_ACC6(0,0,0))  &
              - dydtyr(OPS_ACC4(0,0,0)) - strvyr(OPS_ACC7(0,0,0))*bclyyr(OPS_ACC2(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyyr(OPS_ACC2(0,0,0))*strdyr(OPS_ACC6(0,0,0))*strhyr(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqV_yr

SUBROUTINE bounds_kernel_eqV_zr(erhs,bclyzr,ratezr,dydtzr,strhzr,strdzr,strwzr)

    real(kind=8), dimension(1) :: erhs,bclyzr
    real(kind=8), dimension(1), intent(in) :: ratezr,dydtzr,strhzr,strdzr,strwzr

    bclyzr(OPS_ACC2(0,0,0)) = ratezr(OPS_ACC3(0,0,0))/strdzr(OPS_ACC6(0,0,0))  &
              - dydtzr(OPS_ACC4(0,0,0)) - strwzr(OPS_ACC7(0,0,0))*bclyzr(OPS_ACC2(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
              - bclyzr(OPS_ACC2(0,0,0))*strdzr(OPS_ACC6(0,0,0))*strhzr(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqV_zr

SUBROUTINE bounds_kernel_eqW_xr(drhs,urhs,vrhs,wrhs,erhs,bcl1xr,bcl2xr,ova2xr)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1xr,bcl2xr,ova2xr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1xr(OPS_ACC6(0,0,0))*ova2xr(OPS_ACC8(0,0,0))  &
                          - bcl2xr(OPS_ACC7(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = 0.0_8
    vrhs(OPS_ACC3(0,0,0)) = 0.0_8
    wrhs(OPS_ACC4(0,0,0)) = 0.0_8
    erhs(OPS_ACC5(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqW_xr

SUBROUTINE bounds_kernel_eqW_yr(drhs,urhs,vrhs,wrhs,erhs,bcl1yr,bcl2yr,ova2yr)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1yr,bcl2yr,ova2yr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1yr(OPS_ACC6(0,0,0))*ova2yr(OPS_ACC8(0,0,0))  &
                          - bcl2yr(OPS_ACC7(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = 0.0_8
    vrhs(OPS_ACC3(0,0,0)) = 0.0_8
    wrhs(OPS_ACC4(0,0,0)) = 0.0_8
    erhs(OPS_ACC5(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqW_yr

SUBROUTINE bounds_kernel_eqW_zr(drhs,urhs,vrhs,wrhs,erhs,bcl1zr,bcl2zr,ova2zr)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1zr,bcl2zr,ova2zr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1zr(OPS_ACC6(0,0,0))*ova2zr(OPS_ACC8(0,0,0))  &
                          - bcl2zr(OPS_ACC7(0,0,0))
    urhs(OPS_ACC2(0,0,0)) = 0.0_8
    vrhs(OPS_ACC3(0,0,0)) = 0.0_8
    wrhs(OPS_ACC4(0,0,0)) = 0.0_8
    erhs(OPS_ACC5(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqW_zr

SUBROUTINE bounds_kernel_eqX_xr(yrhs,stryxr,bcl1xr,bcl2xr,ova2xr)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryxr,bcl1xr,bcl2xr,ova2xr

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
                          - (bcl2xr(OPS_ACC4(0,0,0)) + bcl1xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC5(0,0,0))) &
                             * stryxr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqX_xr

SUBROUTINE bounds_kernel_eqX_yr(yrhs,stryyr,bcl1yr,bcl2yr,ova2yr)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryyr,bcl1yr,bcl2yr,ova2yr

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
                          - (bcl2yr(OPS_ACC4(0,0,0)) + bcl1yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC5(0,0,0))) &
                             * stryyr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqX_yr

SUBROUTINE bounds_kernel_eqX_zr(yrhs,stryzr,bcl1zr,bcl2zr,ova2zr)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryzr,bcl1zr,bcl2zr,ova2zr

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
                          - (bcl2zr(OPS_ACC4(0,0,0)) + bcl1zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC5(0,0,0))) &
                             * stryzr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqX_zr

SUBROUTINE bounds_kernel_eqY_xr(drhs,bcl1xr,bcl2xr,ova2xr,urhs,vrhs,wrhs,erhs)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1xr,bcl2xr,ova2xr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1xr(OPS_ACC2(0,0,0))*ova2xr(OPS_ACC4(0,0,0))  &
                    - bcl2xr(OPS_ACC3(0,0,0))
    urhs(OPS_ACC5(0,0,0)) = 0.0_8
    vrhs(OPS_ACC6(0,0,0)) = 0.0_8
    wrhs(OPS_ACC7(0,0,0)) = 0.0_8
    erhs(OPS_ACC8(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqY_xr

SUBROUTINE bounds_kernel_eqY_yr(drhs,bcl1yr,bcl2yr,ova2yr,urhs,vrhs,wrhs,erhs)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1yr,bcl2yr,ova2yr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1yr(OPS_ACC2(0,0,0))*ova2yr(OPS_ACC4(0,0,0))  &
                    - bcl2yr(OPS_ACC3(0,0,0))
    urhs(OPS_ACC5(0,0,0)) = 0.0_8
    vrhs(OPS_ACC6(0,0,0)) = 0.0_8
    wrhs(OPS_ACC7(0,0,0)) = 0.0_8
    erhs(OPS_ACC8(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqY_yr

SUBROUTINE bounds_kernel_eqY_zr(drhs,bcl1zr,bcl2zr,ova2zr,urhs,vrhs,wrhs,erhs)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl1zr,bcl2zr,ova2zr

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl1zr(OPS_ACC2(0,0,0))*ova2zr(OPS_ACC4(0,0,0))  &
                    - bcl2zr(OPS_ACC3(0,0,0))
    urhs(OPS_ACC5(0,0,0)) = 0.0_8
    vrhs(OPS_ACC6(0,0,0)) = 0.0_8
    wrhs(OPS_ACC7(0,0,0)) = 0.0_8
    erhs(OPS_ACC8(0,0,0)) = 0.0_8

END SUBROUTINE bounds_kernel_eqY_zr

SUBROUTINE bounds_kernel_eqZ_xr(yrhs,stryxr,bcl1xr,bcl2xr,ova2xr)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryxr,bcl1xr,bcl2xr,ova2xr

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
              - (bcl2xr(OPS_ACC4(0,0,0))+bcl1xr(OPS_ACC3(0,0,0))*ova2xr(OPS_ACC5(0,0,0)))*stryxr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqZ_xr

SUBROUTINE bounds_kernel_eqZ_yr(yrhs,stryyr,bcl1yr,bcl2yr,ova2yr)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryyr,bcl1yr,bcl2yr,ova2yr

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
              - (bcl2yr(OPS_ACC4(0,0,0))+bcl1yr(OPS_ACC3(0,0,0))*ova2yr(OPS_ACC5(0,0,0)))*stryyr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqZ_yr

SUBROUTINE bounds_kernel_eqZ_zr(yrhs,stryzr,bcl1zr,bcl2zr,ova2zr)

    real(kind=8), dimension(1) :: yrhs
    real(kind=8), dimension(1), intent(in) :: stryzr,bcl1zr,bcl2zr,ova2zr

    yrhs(OPS_ACC1(0,0,0)) = yrhs(OPS_ACC1(0,0,0))  &
              - (bcl2zr(OPS_ACC4(0,0,0))+bcl1zr(OPS_ACC3(0,0,0))*ova2zr(OPS_ACC5(0,0,0)))*stryzr(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqZ_zr

SUBROUTINE bounds_kernel_eqAA_xdir(tt1x, tt2x, tt5x, t1bx, t2bx, t51bx, t52bx, gam1x, strdx, acoux)

    real(kind=8), dimension(1) :: tt1x, tt2x, tt5x
    real(kind=8), dimension(1), intent(in) :: t1bx, t2bx, t51bx, t52bx, gam1x, strdx, acoux

    tt1x(OPS_ACC1(0,0,0)) = t51bx(OPS_ACC6(0,0,0)) + t52bx(OPS_ACC7(0,0,0))*(gam1x(OPS_ACC8(0,0,0))+1.0_8) &
                          - strdx(OPS_ACC9(0,0,0))*acoux(OPS_ACC10(0,0,0))*t2bx(OPS_ACC5(0,0,0))

    tt2x(OPS_ACC2(0,0,0)) = acoux(OPS_ACC10(0,0,0))*acoux(OPS_ACC10(0,0,0))*t1bx(OPS_ACC4(0,0,0)) &
                          - t51bx(OPS_ACC6(0,0,0)) - (gam1x(OPS_ACC8(0,0,0))+1.0_8)*t52bx(OPS_ACC7(0,0,0))

    tt5x(OPS_ACC3(0,0,0)) = t51bx(OPS_ACC6(0,0,0)) + t52bx(OPS_ACC7(0,0,0))*(gam1x(OPS_ACC8(0,0,0))+1.0_8) &
                          + strdx(OPS_ACC9(0,0,0))*acoux(OPS_ACC10(0,0,0))*t2bx(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqAA_xdir

SUBROUTINE bounds_kernel_eqAA_ydir(tt1y, tt2y, tt5y, t1by, t2by, t51by, t52by, gam1y, strdy, acouy)

    real(kind=8), dimension(1) :: tt1y, tt2y, tt5y
    real(kind=8), dimension(1), intent(in) :: t1by, t2by, t51by, t52by, gam1y, strdy, acouy

    tt1y(OPS_ACC1(0,0,0)) = t51by(OPS_ACC6(0,0,0)) + t52by(OPS_ACC7(0,0,0))*(gam1y(OPS_ACC8(0,0,0))+1.0_8) &
                          - strdy(OPS_ACC9(0,0,0))*acouy(OPS_ACC10(0,0,0))*t2by(OPS_ACC5(0,0,0))

    tt2y(OPS_ACC2(0,0,0)) = acouy(OPS_ACC10(0,0,0))*acouy(OPS_ACC10(0,0,0))*t1by(OPS_ACC4(0,0,0)) &
                          - t51by(OPS_ACC6(0,0,0)) - (gam1y(OPS_ACC8(0,0,0))+1.0_8)*t52by(OPS_ACC7(0,0,0))

    tt5y(OPS_ACC3(0,0,0)) = t51by(OPS_ACC6(0,0,0)) + t52by(OPS_ACC7(0,0,0))*(gam1y(OPS_ACC8(0,0,0))+1.0_8) &
                          + strdy(OPS_ACC9(0,0,0))*acouy(OPS_ACC10(0,0,0))*t2by(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqAA_ydir

SUBROUTINE bounds_kernel_eqAA_zdir(tt1z, tt2z, tt5z, t1bz, t2bz, t51bz, t52bz, gam1z, strdz, acouz)

    real(kind=8), dimension(1) :: tt1z, tt2z, tt5z
    real(kind=8), dimension(1), intent(in) :: t1bz, t2bz, t51bz, t52bz, gam1z, strdz, acouz

    tt1z(OPS_ACC1(0,0,0)) = t51bz(OPS_ACC6(0,0,0)) + t52bz(OPS_ACC7(0,0,0))*(gam1z(OPS_ACC8(0,0,0))+1.0_8) &
                          - strdz(OPS_ACC9(0,0,0))*acouz(OPS_ACC10(0,0,0))*t2bz(OPS_ACC5(0,0,0))

    tt2z(OPS_ACC2(0,0,0)) = acouz(OPS_ACC10(0,0,0))*acouz(OPS_ACC10(0,0,0))*t1bz(OPS_ACC4(0,0,0)) &
                          - t51bz(OPS_ACC6(0,0,0)) - (gam1z(OPS_ACC8(0,0,0))+1.0_8)*t52bz(OPS_ACC7(0,0,0))

    tt5z(OPS_ACC3(0,0,0)) = t51bz(OPS_ACC6(0,0,0)) + t52bz(OPS_ACC7(0,0,0))*(gam1z(OPS_ACC8(0,0,0))+1.0_8) &
                          + strdz(OPS_ACC9(0,0,0))*acouz(OPS_ACC10(0,0,0))*t2bz(OPS_ACC5(0,0,0))

END SUBROUTINE bounds_kernel_eqAA_zdir

SUBROUTINE bounds_kernel_eqAB_xdir(bcl2x,bcl3x,bcl4x,bcl5x,bcl1x,strdx,acoux,ova2x,tt2x,tt3x,tt4x,tt5x,strrx,stlux,strux,stlvx,strvx,stlwx,strwx,stltx,strtx,xgdlen,nrieta2,nrieta3,nrieta4,nrieta5,m2max)

    real(kind=8), dimension(1) :: bcl2x,bcl3x,bcl4x,bcl5x
    real(kind=8), dimension(1), intent(in) :: bcl1x,strdx,acoux,ova2x,tt2x,tt3x,tt4x,tt5x
    real(kind=8), dimension(1), intent(in) :: strrx,stlux,strux,stlvx,strvx,stlwx,strwx,stltx,strtx
    real(kind=8), intent(in) :: xgdlen,nrieta2,nrieta3,nrieta4,nrieta5,m2max
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdx(OPS_ACC6(0,0,0))*acoux(OPS_ACC7(0,0,0))*bcl1x(OPS_ACC5(0,0,0))

    bcl2x(OPS_ACC1(0,0,0)) = stlux(OPS_ACC14(0,0,0))  &
                  * (bcl2x(OPS_ACC1(0,0,0))-bcl5x(OPS_ACC4(0,0,0))*ova2x(OPS_ACC8(0,0,0)))

    bcl3x(OPS_ACC2(0,0,0)) = stlux(OPS_ACC14(0,0,0))*bcl3x(OPS_ACC2(0,0,0))

    bcl4x(OPS_ACC3(0,0,0)) = stlux(OPS_ACC14(0,0,0))*bcl4x(OPS_ACC3(0,0,0))

    bcl5x(OPS_ACC4(0,0,0)) = 0.5_8*(stlux(OPS_ACC14(0,0,0))+acoux(OPS_ACC7(0,0,0)))  &
                  * (bcl5x(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1X UNCHANGED
    bcl2x(OPS_ACC1(0,0,0)) = nrieta2 * strdx(OPS_ACC6(0,0,0))  &
                  * strrx(OPS_ACC13(0,0,0))/(xgdlen*acoux(OPS_ACC7(0,0,0)))  &
                  * (stltx(OPS_ACC20(0,0,0))-strtx(OPS_ACC21(0,0,0)))  &
                  + tt2x(OPS_ACC9(0,0,0))*ova2x(OPS_ACC8(0,0,0))-bcl2x(OPS_ACC1(0,0,0))

    bcl3x(OPS_ACC2(0,0,0)) = nrieta3 * acoux(OPS_ACC7(0,0,0))  &
                  / xgdlen*(stlvx(OPS_ACC16(0,0,0))-strvx(OPS_ACC17(0,0,0))) +tt3x(OPS_ACC10(0,0,0))-bcl3x(OPS_ACC2(0,0,0))

    bcl4x(OPS_ACC3(0,0,0)) = nrieta4 * acoux(OPS_ACC7(0,0,0))  &
                  / xgdlen*(stlwx(OPS_ACC18(0,0,0))-strwx(OPS_ACC19(0,0,0))) +tt4x(OPS_ACC11(0,0,0))-bcl4x(OPS_ACC3(0,0,0))

    bcl5x(OPS_ACC4(0,0,0)) = nrieta5 *strdx(OPS_ACC6(0,0,0))*acoux(OPS_ACC7(0,0,0))  &
                  * acoux(OPS_ACC7(0,0,0))*(1.0_8-m2max)/(2.0_8*xgdlen)  &
                  * (stlux(OPS_ACC14(0,0,0))-strux(OPS_ACC15(0,0,0))) +0.5_8*tt5x(OPS_ACC12(0,0,0))  &
                  - bcl5x(OPS_ACC4(0,0,0))

END SUBROUTINE bounds_kernel_eqAB_xdir

SUBROUTINE bounds_kernel_eqAC_xdir(bclyx,stlyx,stryx,tt6x,stlux,acoux,xgdlen,nrieta6)

    real(kind=8), dimension(1) :: bclyx
    real(kind=8), dimension(1), intent(in) :: stlyx,stryx,tt6x,stlux,acoux
    real(kind=8), intent(in) :: xgdlen,nrieta6

!   OLD VALUE OF L's
    bclyx(OPS_ACC1(0,0,0)) = stlux(OPS_ACC5(0,0,0)) *bclyx(OPS_ACC1(0,0,0))

!   SUBTRACT FROM NEW VALUE OF L's
    bclyx(OPS_ACC1(0,0,0)) = nrieta6 *acoux(OPS_ACC6(0,0,0))/xgdlen  &
                        * (stlyx(OPS_ACC2(0,0,0)) - stryx(OPS_ACC3(0,0,0)))  &
                        + tt6x(OPS_ACC4(0,0,0)) - bclyx(OPS_ACC1(0,0,0))


END SUBROUTINE bounds_kernel_eqAC_xdir

SUBROUTINE bounds_kernel_eqAD_xdir(drhs,urhs,vrhs,wrhs,erhs,bcl2x,bcl3x,bcl4x,bcl5x,ova2x,acoux,ovgmx,stlux,stlvx,stlwx,strdx,strex)

    real(kind=8), dimension(1) :: drhs,urhs,vrhs,wrhs,erhs
    real(kind=8), dimension(1), intent(in) :: bcl2x,bcl3x,bcl4x,bcl5x
    real(kind=8), dimension(1), intent(in) :: ova2x,acoux,ovgmx,stlux,stlvx,stlwx,strdx,strex

    drhs(OPS_ACC1(0,0,0)) = drhs(OPS_ACC1(0,0,0)) - bcl2x(OPS_ACC6(0,0,0))  &
                      - bcl5x(OPS_ACC9(0,0,0))*ova2x(OPS_ACC10(0,0,0))

    urhs(OPS_ACC2(0,0,0)) = urhs(OPS_ACC2(0,0,0)) - bcl2x(OPS_ACC6(0,0,0))*stlux(OPS_ACC13(0,0,0))  &
                      - bcl5x(OPS_ACC9(0,0,0))*ova2x(OPS_ACC10(0,0,0))*(stlux(OPS_ACC13(0,0,0)) + acoux(OPS_ACC11(0,0,0)))

    vrhs(OPS_ACC3(0,0,0)) = vrhs(OPS_ACC3(0,0,0)) - bcl2x(OPS_ACC6(0,0,0))*stlvx(OPS_ACC14(0,0,0))  &
                      - bcl3x(OPS_ACC7(0,0,0))*strdx(OPS_ACC16(0,0,0))  &
                      - bcl5x(OPS_ACC9(0,0,0))*ova2x(OPS_ACC10(0,0,0))*stlvx(OPS_ACC14(0,0,0))

    wrhs(OPS_ACC4(0,0,0)) = wrhs(OPS_ACC4(0,0,0)) - bcl2x(OPS_ACC6(0,0,0))*stlwx(OPS_ACC15(0,0,0))  &
                      - bcl4x(OPS_ACC8(0,0,0))*strdx(OPS_ACC16(0,0,0))  &
                      - bcl5x(OPS_ACC9(0,0,0))*ova2x(OPS_ACC10(0,0,0))*stlwx(OPS_ACC15(0,0,0))

    erhs(OPS_ACC5(0,0,0)) = erhs(OPS_ACC5(0,0,0)) - bcl2x(OPS_ACC6(0,0,0))*strex(OPS_ACC17(0,0,0))  &
                      - bcl3x(OPS_ACC7(0,0,0))*strdx(OPS_ACC16(0,0,0))*stlvx(OPS_ACC14(0,0,0))  &
                      - bcl4x(OPS_ACC8(0,0,0))*strdx(OPS_ACC16(0,0,0))*stlwx(OPS_ACC15(0,0,0))  &
                      - bcl5x(OPS_ACC9(0,0,0))*(ova2x(OPS_ACC10(0,0,0))*strex(OPS_ACC17(0,0,0))  &
                                     + stlux(OPS_ACC13(0,0,0))/acoux(OPS_ACC11(0,0,0)) + ovgmx(OPS_ACC12(0,0,0)))

END SUBROUTINE bounds_kernel_eqAD_xdir

SUBROUTINE bounds_kernel_eqAE_xdir(erhs,yrhs,bclyx,stlyx,strhx,strdx,bcl2x,bcl5x,ova2x)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: bclyx,stlyx,strhx,strdx,bcl2x,bcl5x,ova2x
    real(kind=8) :: fornow

    fornow = bclyx(OPS_ACC3(0,0,0))*strdx(OPS_ACC6(0,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhx(OPS_ACC5(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0))  &
                            - (bcl2x(OPS_ACC7(0,0,0)) + bcl5x(OPS_ACC8(0,0,0))*ova2x(OPS_ACC9(0,0,0)))*stlyx(OPS_ACC4(0,0,0))  &
                            - fornow

END SUBROUTINE bounds_kernel_eqAE_xdir

SUBROUTINE bounds_kernel_eqAF_xl(strtxl,trun)

    real(kind=8), dimension(1) :: strtxl, trun

    strtxl(OPS_ACC1(0,0,0)) = (48.0_8*trun(OPS_ACC2(1,0,0)) - 36.0_8*trun(OPS_ACC2(2,0,0))  &
                             + 16.0_8*trun(OPS_ACC2(3,0,0)) - 3.0_8*trun(OPS_ACC2(4,0,0)))/25.0_8

    trun(OPS_ACC2(0,0,0)) = strtxl(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqAF_xl

SUBROUTINE bounds_kernel_eqAF_xr(strtxr,trun)

    real(kind=8), dimension(1) :: strtxr, trun

    strtxr(OPS_ACC1(0,0,0)) = (48.0_8*trun(OPS_ACC2(-1,0,0)) - 36.0_8*trun(OPS_ACC2(-2,0,0))  &
                             + 16.0_8*trun(OPS_ACC2(-3,0,0)) - 3.0_8*trun(OPS_ACC2(-4,0,0)))/25.0_8

    trun(OPS_ACC2(0,0,0)) = strtxr(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqAF_xr

SUBROUTINE bounds_kernel_eqAF_yl(strtyl,trun)

    real(kind=8), dimension(1) :: strtyl, trun

    strtyl(OPS_ACC1(0,0,0)) = (48.0_8*trun(OPS_ACC2(0,1,0)) - 36.0_8*trun(OPS_ACC2(0,2,0))  &
                             + 16.0_8*trun(OPS_ACC2(0,3,0)) - 3.0_8*trun(OPS_ACC2(0,4,0)))/25.0_8

    trun(OPS_ACC2(0,0,0)) = strtyl(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqAF_yl

SUBROUTINE bounds_kernel_eqAF_yr(strtyr,trun)

    real(kind=8), dimension(1) :: strtyr, trun

    strtyr(OPS_ACC1(0,0,0)) = (48.0_8*trun(OPS_ACC2(0,-1,0)) - 36.0_8*trun(OPS_ACC2(0,-2,0))  &
                             + 16.0_8*trun(OPS_ACC2(0,-3,0)) - 3.0_8*trun(OPS_ACC2(0,-4,0)))/25.0_8

    trun(OPS_ACC2(0,0,0)) = strtyr(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqAF_yr

SUBROUTINE bounds_kernel_eqAF_zl(strtzl,trun)

    real(kind=8), dimension(1) :: strtzl, trun

    strtzl(OPS_ACC1(0,0,0)) = (48.0_8*trun(OPS_ACC2(0,0,1)) - 36.0_8*trun(OPS_ACC2(0,0,2))  &
                             + 16.0_8*trun(OPS_ACC2(0,0,3)) - 3.0_8*trun(OPS_ACC2(0,0,4)))/25.0_8

    trun(OPS_ACC2(0,0,0)) = strtzl(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqAF_zl

SUBROUTINE bounds_kernel_eqAF_zr(strtzr,trun)

    real(kind=8), dimension(1) :: strtzr, trun

    strtzr(OPS_ACC1(0,0,0)) = (48.0_8*trun(OPS_ACC2(0,0,-1)) - 36.0_8*trun(OPS_ACC2(0,0,-2))  &
                             + 16.0_8*trun(OPS_ACC2(0,0,-3)) - 3.0_8*trun(OPS_ACC2(0,0,-4)))/25.0_8

    trun(OPS_ACC2(0,0,0)) = strtzr(OPS_ACC1(0,0,0))

END SUBROUTINE bounds_kernel_eqAF_zr

SUBROUTINE bounds_kernel_eqAG_xdir(sorpx,strhx,ratex)

    real(kind=8), dimension(1) :: sorpx
    real(kind=8), dimension(1), intent(in) :: strhx,ratex

    sorpx(OPS_ACC1(0,0,0)) = sorpx(OPS_ACC1(0,0,0))  &
                           + strhx(OPS_ACC2(0,0,0))*ratex(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_eqAG_xdir

SUBROUTINE bounds_kernel_eqAG_ydir(sorpy,strhy,ratey)

    real(kind=8), dimension(1) :: sorpy
    real(kind=8), dimension(1), intent(in) :: strhy,ratey

    sorpy(OPS_ACC1(0,0,0)) = sorpy(OPS_ACC1(0,0,0))  &
                           + strhy(OPS_ACC2(0,0,0))*ratey(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_eqAG_ydi

SUBROUTINE bounds_kernel_eqAG_zdir(sorpz,strhz,ratez)

    real(kind=8), dimension(1) :: sorpz
    real(kind=8), dimension(1), intent(in) :: strhz,ratez

    sorpz(OPS_ACC1(0,0,0)) = sorpz(OPS_ACC1(0,0,0))  &
                           + strhz(OPS_ACC2(0,0,0))*ratez(OPS_ACC3(0,0,0))

END SUBROUTINE bounds_kernel_eqAG_zdir

SUBROUTINE bounds_kernel_eqAH_xdir(sorpx,gam1x)

    real(kind=8), dimension(1) :: sorpx
    real(kind=8), dimension(1), intent(in) :: gam1x

    sorpx(OPS_ACC1(0,0,0)) = -sorpx(OPS_ACC1(0,0,0))*gam1x(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqAH_xdir

SUBROUTINE bounds_kernel_eqAH_ydir(sorpy,gam1y)

    real(kind=8), dimension(1) :: sorpy
    real(kind=8), dimension(1), intent(in) :: gam1y

    sorpy(OPS_ACC1(0,0,0)) = -sorpy(OPS_ACC1(0,0,0))*gam1y(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqAH_ydir

SUBROUTINE bounds_kernel_eqAH_zdir(sorpz,gam1z)

    real(kind=8), dimension(1) :: sorpz
    real(kind=8), dimension(1), intent(in) :: gam1z

    sorpz(OPS_ACC1(0,0,0)) = -sorpz(OPS_ACC1(0,0,0))*gam1z(OPS_ACC2(0,0,0))

END SUBROUTINE bounds_kernel_eqAH_zdir
