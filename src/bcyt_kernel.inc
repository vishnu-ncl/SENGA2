SUBROUTINE bcyt_kernel_xdir_eqA(stryx,dydtx,yrin,ispec)

    real(kind=8), dimension(1) :: stryx,dydtx
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: yrin

!   SET MASS FRACTIONS TO CONSTANT (INITIAL) VALUES
    stryx(OPS_ACC1(0,0,0)) = yrin(ispec)

!   SET MASS FRACTION TIME DERIVATIVES TO ZERO
    dydtx(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcyt_kernel_xdir_eqA

SUBROUTINE bcyt_kernel_xdir_eqB(stryx,dydtx,yinf2x,yinf1x,yrin,tstep,ispec)

    real(kind=8), dimension(1) :: stryx, dydtx, yinf2x
    real(kind=8), dimension(1), intent(in) :: yinf1x
    real(kind=8), dimension(nspcmx), intent(in) :: yrin
    real(kind=8), intent(in) :: tstep
    integer(kind=4), intent(in) :: ispec

    stryx(OPS_ACC1(0,0,0)) = yrin(ispec) + yinf2x(OPS_ACC3(0,0,0))

    IF(stryx(OPS_ACC1(0,0,0)) > 1.0_8) THEN
          yinf2x(OPS_ACC3(0,0,0)) = 1.0_8 -yrin(ispec)
          stryx(OPS_ACC1(0,0,0)) = 1.0_8
    END IF

    IF(stryx(OPS_ACC1(0,0,0)) < 0.0_8) THEN
        yinf2x(OPS_ACC3(0,0,0)) = yrin(ispec)-0.0_8
        stryx(OPS_ACC1(0,0,0)) = 0.0_8
    END IF

    dydtx(OPS_ACC2(0,0,0)) = (yinf2x(OPS_ACC3(0,0,0))- yinf1x(OPS_ACC4(0,0,0)))/tstep

END SUBROUTINE bcyt_kernel_xdir_eqB

SUBROUTINE bcyt_kernel_xdir_eqC(totyx,stryx)

    real(kind=8), dimension(1) :: totyx
    real(kind=8), dimension(1), intent(in) :: stryx

     totyx(OPS_ACC1(0,0,0)) = totyx(OPS_ACC1(0,0,0)) + stryx(OPS_ACC2(0,0,0))

END SUBROUTINE bcyt_kernel_xdir_eqC

SUBROUTINE bcyt_kernel_xdir_eqD(stryx,dydtx,yinf2x,yinf1x,totyx,tstep)

    real(kind=8), dimension(1) :: stryx, dydtx
    real(kind=8), dimension(1), intent(in) :: yinf2x, yinf1x, totyx
    real(kind=8), intent(in) :: tstep

    stryx(OPS_ACC1(0,0,0)) = 1.0_8 - totyx(OPS_ACC5(0,0,0))

    stryx(OPS_ACC1(0,0,0)) = MAX(0.0_8, stryx(OPS_ACC1(0,0,0)))

    stryx(OPS_ACC1(0,0,0)) = MIN(1.0_8, stryx(OPS_ACC1(0,0,0)))

    dydtx(OPS_ACC2(0,0,0)) = (yinf2x(OPS_ACC3(0,0,0)) -yinf1x(OPS_ACC4(0,0,0)))/tstep

END SUBROUTINE bcyt_kernel_xdir_eqD

SUBROUTINE bcyt_kernel_ydir(stryy,dydty,yrin,ispec)

    real(kind=8), dimension(1) :: stryy,dydty
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: yrin

    stryy(OPS_ACC1(0,0,0)) = yrin(ispec)
    dydty(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcyt_kernel_ydir

SUBROUTINE bcyt_kernel_zdir(stryz,dydtz,yrin,ispec)

    real(kind=8), dimension(1) :: stryz,dydtz
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: yrin

    stryz(OPS_ACC1(0,0,0)) = yrin(ispec)
    dydtz(OPS_ACC2(0,0,0)) = 0.0_8

END SUBROUTINE bcyt_kernel_zdir
