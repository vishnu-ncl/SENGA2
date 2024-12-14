SUBROUTINE bcut_kernel_ydir(struy, strvy, strwy, dudty, dvdty, dwdty, init_val)

    real(kind=8), dimension(1) :: struy, strvy, strwy, dudty, dvdty, dwdty
    real(kind=8), dimension(nbcprr), intent(in) :: init_val

    struy(OPS_ACC1(0,0,0)) = 0.0_8
    strvy(OPS_ACC2(0,0,0)) = init_val(1)
    strwy(OPS_ACC3(0,0,0)) = 0.0_8

    dudty(OPS_ACC4(0,0,0)) = 0.0_8
    dvdty(OPS_ACC5(0,0,0)) = 0.0_8
    dwdty(OPS_ACC6(0,0,0)) = 0.0_8

END SUBROUTINE bcut_kernel_ydir

SUBROUTINE bcut_kernel_zdir(struz, strvz, strwz, dudtz, dvdtz, dwdtz, init_val)

    real(kind=8), dimension(1) :: struz, strvz, strwz, dudtz, dvdtz, dwdtz
    real(kind=8), dimension(nbcprr), intent(in) :: init_val

    struz(OPS_ACC1(0,0,0)) = 0.0_8
    strvz(OPS_ACC2(0,0,0)) = 0.0_8
    strwz(OPS_ACC3(0,0,0)) = init_val(1)

    dudtz(OPS_ACC4(0,0,0)) = 0.0_8
    dvdtz(OPS_ACC5(0,0,0)) = 0.0_8
    dwdtz(OPS_ACC6(0,0,0)) = 0.0_8

END SUBROUTINE bcut_kernel_zdir

!A = A + B*val1
SUBROUTINE bcut_kernel_xdir_eqA(out_arr, in_arr, val1)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    real(kind=8), intent(in) :: val1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + in_arr(OPS_ACC2(0,0,0))*val1

END SUBROUTINE bcut_kernel_xdir_eqA

!A = A - val1*B*val2
SUBROUTINE bcut_kernel_xdir_eqB(out_arr, in_arr, val1, val2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    real(kind=8), intent(in) :: val1, val2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) - val1*in_arr(OPS_ACC2(0,0,0))*val2

END SUBROUTINE bcut_kernel_xdir_eqB

!A = A + val1*B*val2
SUBROUTINE bcut_kernel_xdir_eqC(out_arr, in_arr, val1, val2)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    real(kind=8), intent(in) :: val1, val2

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + val1*in_arr(OPS_ACC2(0,0,0))*val2

END SUBROUTINE bcut_kernel_xdir_eqC

!A = A * val1
SUBROUTINE bcut_kernel_xdir_eqD(out_arr, val1)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: val1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) * val1

END SUBROUTINE bcut_kernel_xdir_eqD

!A = A + val1
SUBROUTINE bcut_kernel_xdir_eqE(out_arr, val1)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), intent(in) :: val1

    out_arr(OPS_ACC1(0,0,0)) = out_arr(OPS_ACC1(0,0,0)) + val1

END SUBROUTINE bcut_kernel_xdir_eqE

SUBROUTINE bcut_kernel_xdir_eqF(strux, strvx, strwx, dudtx, dvdtx, dwdtx, init_val)

    real(kind=8), dimension(1) :: strux, strvx, strwx, dudtx, dvdtx, dwdtx
    real(kind=8), intent(in) :: init_val

    strux(OPS_ACC1(0,0,0)) = init_val
    strvx(OPS_ACC2(0,0,0)) = 0.0_8
    strwx(OPS_ACC3(0,0,0)) = 0.0_8

    dudtx(OPS_ACC4(0,0,0)) = 0.0_8
    dvdtx(OPS_ACC5(0,0,0)) = 0.0_8
    dwdtx(OPS_ACC6(0,0,0)) = 0.0_8

END SUBROUTINE bcut_kernel_xdir_eqF

SUBROUTINE bcut_kernel_xdir_eqG(strux, strvx, strwx, dudtx, dvdtx, dwdtx, init_val1, init_val2)

    real(kind=8), dimension(1) :: strux, strvx, strwx, dudtx, dvdtx, dwdtx
    real(kind=8), intent(in) :: init_val1, init_val2

    strux(OPS_ACC1(0,0,0)) = init_val1
    strvx(OPS_ACC2(0,0,0)) = 0.0_8
    strwx(OPS_ACC3(0,0,0)) = 0.0_8

    dudtx(OPS_ACC4(0,0,0)) = init_val2
    dvdtx(OPS_ACC5(0,0,0)) = 0.0_8
    dwdtx(OPS_ACC6(0,0,0)) = 0.0_8

END SUBROUTINE bcut_kernel_xdir_eqG

SUBROUTINE bcut_kernel_xdir_eqH(strux,strvx,strwx,dudtx,dvdtx,dwdtx,slope,widthp,ptly,rxlprm_1,rxlprc_7,fornow,argmnt,deltagy,ygdlen,idx)

    real(kind=8), dimension(1) :: strux, strvx, strwx, dudtx, dvdtx, dwdtx
    real(kind=8), intent(in) :: slope,widthp,ptly,rxlprm_1,rxlprc_7,fornow,argmnt,deltagy,ygdlen
    integer(kind=4), dimension(3), intent(in) :: idx

    integer(kind=4) :: jc
    real(kind=8) :: hyptan  ! FOR HYPERBOLIC TANGENT PROFILE

    jc = idx(2)

!   DEFINE HYPERBOLIC TANGENT PROFILE
    hyptan = ( 0.5_8*TANH(slope * (REAL(jc-1,kind=8)*deltagy-(ygdlen/2.0_8-widthp*ptly) )) + 0.5_8) &
           * (-0.5_8*TANH(slope * (REAL(jc-1,kind=8)*deltagy-(ygdlen/2.0_8+widthp*ptly) )) + 0.5_8)

    strux(OPS_ACC1(0,0,0)) = rxlprm_1+hyptan*rxlprc_7*SIN(argmnt)
    strvx(OPS_ACC2(0,0,0)) = 0.0_8
    strwx(OPS_ACC3(0,0,0)) = 0.0_8

    dudtx(OPS_ACC4(0,0,0)) = hyptan*fornow*rxlprc_7*COS(argmnt)
    dvdtx(OPS_ACC5(0,0,0)) = 0.0_8
    dwdtx(OPS_ACC6(0,0,0)) = 0.0_8

END SUBROUTINE bcut_kernel_xdir_eqH

SUBROUTINE bcut_kernel_xdir_eqI(strux,strvx,strwx,dudtx,dvdtx,dwdtx,ustead,uinf1,vinf1,winf1,uinf2,vinf2,winf2,rxlprm_1,vfac,coflow,tstep,dvfdt)

    real(kind=8), dimension(1) :: strux,strvx,strwx,dudtx,dvdtx,dwdtx
    real(kind=8), dimension(1), intent(in) :: ustead,uinf1,vinf1,winf1,uinf2,vinf2,winf2
    real(kind=8), intent(in) :: rxlprm_1,vfac,coflow,tstep,dvfdt

    strux(OPS_ACC1(0,0,0)) = (rxlprm_1*ustead(OPS_ACC7(0,0,0)) + uinf2(OPS_ACC11(0,0,0))) * vfac + coflow
    strvx(OPS_ACC2(0,0,0)) = vinf2(OPS_ACC12(0,0,0)) * vfac
    strwx(OPS_ACC3(0,0,0)) = winf2(OPS_ACC13(0,0,0)) * vfac


    dudtx(OPS_ACC4(0,0,0)) = (uinf2(OPS_ACC11(0,0,0))-uinf1(OPS_ACC8(0,0,0)))/tstep  &
                            * vfac + (rxlprm_1*ustead(OPS_ACC7(0,0,0))+uinf2(OPS_ACC11(0,0,0))) * dvfdt
    dvdtx(OPS_ACC5(0,0,0)) = (vinf2(OPS_ACC12(0,0,0))-vinf1(OPS_ACC9(0,0,0)))/tstep  &
                            * vfac + vinf2(OPS_ACC12(0,0,0)) * dvfdt
    dwdtx(OPS_ACC6(0,0,0)) = (winf2(OPS_ACC13(0,0,0))-winf1(OPS_ACC10(0,0,0)))/tstep  &
                            * vfac + winf2(OPS_ACC13(0,0,0)) * dvfdt

END SUBROUTINE bcut_kernel_xdir_eqI
