SUBROUTINE flamin_kernel_set_velocity_tgv(urun,vrun,wrun,prun,trun,psi,prin,drin,u0,deltagx,deltagy,deltagz,angfrx,angfry,angfrz,xgdlen,cpsi,rpsi,idx)

    real(kind=8), dimension(1) :: urun, vrun, wrun, prun, trun, psi
    real(kind=8), intent(in) :: prin,drin,u0, deltagx, deltagy, deltagz
    real(kind=8), intent(in) :: angfrx, angfry, angfrz
    real(kind=8), intent(in) :: xgdlen, cpsi, rpsi
    integer(kind=4), dimension(3), intent(in) :: idx

    real(kind=8) :: fornow,rd
    real(kind=8) :: xx,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11
    real(kind=8) :: xcoord,ycoord,zcoord
    real(kind=8) :: xrgmnt,yrgmnt,zrgmnt
    integer(kind=4) :: ix, iy, iz

    ix = idx(1)
    iy = idx(2)
    iz = idx(3)

    xcoord = REAL(ix-1,kind=8)*deltagx
    ycoord = REAL(iy-1,kind=8)*deltagy
    zcoord = REAL(iz-1,kind=8)*deltagz

    xrgmnt = angfrx*xcoord
    yrgmnt = angfry*ycoord
    zrgmnt = angfrz*zcoord

!   set taylor-green vortex velocity field
    urun(OPS_ACC1(0,0,0)) = u0 * SIN(xrgmnt) * COS(yrgmnt) * COS(zrgmnt)
    vrun(OPS_ACC2(0,0,0)) = -u0 * COS(xrgmnt) * SIN(yrgmnt) * COS(zrgmnt)
    wrun(OPS_ACC3(0,0,0)) = 0.0_8

!   set pressure profile assuming constant density
    xrgmnt = 2.0_8*xrgmnt
    yrgmnt = 2.0_8*yrgmnt
    zrgmnt = 2.0_8*zrgmnt

    prun(OPS_ACC4(0,0,0)) = prin + ((drin*u0*u0)/16.0_8) * (COS(xrgmnt)+COS(yrgmnt)) &
                                   * (COS(zrgmnt)+2.0_8)

    rd = ABS(xcoord-0.5_8*xgdlen)
    fornow = cpsi*(rd-rpsi)/rpsi
    fornow = (EXP(fornow)-EXP(-fornow)) / (EXP(fornow)+EXP(-fornow))
    psi(OPS_ACC6(0,0,0)) = 0.5_8 * (1.0_8+fornow)

    xx=ABS(xcoord-0.5*xgdlen)/xgdlen

    p1 = 307.0597_8
    p2 = -1200.3_8
    p3 = 246980.0_8
    p4 = -14803000.0_8
    p5 = 458540000.0_8
    p6 = -7629300000.0_8
    p7 = 70736000000.0_8
    p8 = 336430000000.0_8
    p9 = 637340000000.0_8

    IF((xx >= 0.0_8).AND.(xx <= 0.129_8)) THEN

        trun(OPS_ACC5(0,0,0)) = p1 + p2*xx &
                       + p3*xx**2.0_8 &
                       + p4*xx**3.0_8 &
                       + p5*xx**4.0_8 &
                       + p6*xx**5.0_8 &
                       + p7*xx**6.0_8 &
                       + p8*xx**7.0_8 &
                       + p9*xx**8.0_8

    END IF

    p1 = -15833000.0_8
    p2 = 733200000.0_8
    p3 = -15154000000.0_8
    p4 = 184100000000.0_8
    p5 = -1456500000000.0_8
    p6 = 7842400000000.0_8
    p7 = -29107000000000.0_8
    p8 = 73543000000000.0_8
    p9 = -121070000000000.0_8
    p10 = 117300000000000.0_8
    p11 = -50789000000000.0_8

    IF((xx >= 0.129_8).AND.(xx <= 0.3057_8)) THEN

        trun(OPS_ACC5(0,0,0)) = p1 + p2*xx + p3*xx**2.0_8 &
                       + p4*xx**3.0_8 &
                       + p5*xx**4.0_8 &
                       + p6*xx**5.0_8 &
                       + p7*xx**6.0_8 &
                       + p8*xx**7.0_8 &
                       + p9*xx**8.0_8 &
                       + p10*xx**9.0_8 &
                       + p11*xx**10.0_8

    END IF

    p1 = -50428000.0
    p2 = 1273000000.0
    p3 = -14416000000.0
    p4 = 96442000000.0
    p5 = -422100000000.0
    p6 = 1262900000000.0
    p7 = -2616200000000.0
    p8 = 3705300000000.0
    p9 = -3433700000000.0
    p10 = 1880200000000.0
    p11 = -461970000000.0

    IF(xx > 0.3057_8) THEN

        fornow = p1 + p2*xx + p3*xx**2.0_8 &
               + p4*xx**3.0_8 &
               + p5*xx**4.0_8 &
               + p6*xx**5.0_8 &
               + p7*xx**6.0_8 &
               + p8*xx**7.0_8 &
               + p9*xx**8.0_8 &
               + p10*xx**9.0_8 &
               + p11*xx**10.0_8

        trun(OPS_ACC5(0,0,0)) = MAX(300.0_8,fornow)

    END IF

END SUBROUTINE flamin_kernel_set_velocity_tgv

SUBROUTINE flamin_kernel_eqA(yrun_1,yrun_2,yrun_9,psi,cfo)

    real(kind=8), dimension(1) :: yrun_1,yrun_2,yrun_9
    real(kind=8), dimension(1), intent(in) :: psi 
    real(kind=8), intent(in) :: cfo

    yrun_1(OPS_ACC1(0,0,0)) = cfo * (1.0_8 - psi(OPS_ACC4(0,0,0)))
    yrun_2(OPS_ACC2(0,0,0)) = 0.233_8 * psi(OPS_ACC4(0,0,0))
    yrun_9(OPS_ACC3(0,0,0)) = 1.0_8 - yrun_1(OPS_ACC1(0,0,0)) - yrun_2(OPS_ACC2(0,0,0)) 

END SUBROUTINE flamin_kernel_eqA

SUBROUTINE flamin_kernel_eqB(store1,yrun,rgspec_ispec)

    real(kind=8), dimension(1) :: store1
    real(kind=8), dimension(1), intent(in) :: yrun
    real(kind=8), intent(in) :: rgspec_ispec

    store1(OPS_ACC1(0,0,0)) = store1(OPS_ACC1(0,0,0)) + rgspec_ispec*yrun(OPS_ACC2(0,0,0))

END SUBROUTINE flamin_kernel_eqB

SUBROUTINE flamin_kernel_eqC(drun,prun,trun,store1)

    real(kind=8), dimension(1) :: drun
    real(kind=8), dimension(1), intent(in) :: prun,trun,store1

    drun(OPS_ACC1(0,0,0)) = prun(OPS_ACC2(0,0,0))/(store1(OPS_ACC4(0,0,0))*trun(OPS_ACC3(0,0,0)))

END SUBROUTINE flamin_kernel_eqC
