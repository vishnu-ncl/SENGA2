SUBROUTINE flamin_kernel_set_velocity_tgv(urun,vrun,wrun,psi,prin,drin,u0,deltagx,deltagy,deltagz,angfrx,angfry,angfrz,xgdlen,cpsi,rpsi,idx)

    real(kind=8), dimension(1) :: urun, vrun, wrun, psi
    real(kind=8), intent(in) :: prin,drin,u0, deltagx, deltagy, deltagz
    real(kind=8), intent(in) :: angfrx, angfry, angfrz
    real(kind=8), intent(in) :: xgdlen, cpsi, rpsi
    integer(kind=4), dimension(3), intent(in) :: idx

    real(kind=8) :: fornow,rd
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

    rd = ABS(xcoord - 0.5_8*xgdlen)

    fornow = cpsi*(rd-rpsi)/rpsi
    fornow = (EXP(fornow)-EXP(-fornow)) / (EXP(fornow)+EXP(-fornow))
    psi(OPS_ACC4(0,0,0)) = 0.5_8 * (1.0_8+fornow)

END SUBROUTINE flamin_kernel_set_velocity_tgv

SUBROUTINE flamin_kernel_eqA(yrun,yrunl,xl,deltagx,nspec,idx)

    real(kind=8), dimension(1) :: yrun

    real(kind=8), intent(in), dimension(nl,nspcmx) :: yrunl
    real(kind=8), intent(in), dimension(nl) :: xl
    real(kind=8), intent(in) :: deltagx
    integer(kind=4), intent(in) :: nspec
    integer(kind=4), dimension(3), intent(in) :: idx

    real(kind=8) :: fornow
    real(kind=8) :: xcoord
    integer(kind=4) :: ix,il,ispec

    ix = idx(1)

    xcoord = REAL(ix-1,kind=8)*deltagx

    DO il = 1,nl-1
        IF( (xl(il) <= xcoord) .AND. (xl(il+1) > xcoord) )THEN
            yrun(OPS_ACC_MD1(nspec,0,0,0)) = 1.0_8

            DO ispec = 1,nspec-1
                fornow = (xcoord-xl(il)) / (xl(il+1)-xl(il))

                yrun(OPS_ACC_MD1(ispec,0,0,0)) = yrunl(il,ispec) &
                                               + fornow*(yrunl(il+1,ispec)-yrunl(il,ispec))

                yrun(OPS_ACC_MD1(nspec,0,0,0)) = yrun(OPS_ACC_MD1(nspec,0,0,0))-yrun(OPS_ACC_MD1(ispec,0,0,0))
            END DO
        END IF
    END DO

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

SUBROUTINE flamin_kernel_eqD(trun,xl,theta,deltagx,idx)

    real(kind=8), dimension(1) :: trun

    real(kind=8), intent(in), dimension(nl) :: xl,theta
    real(kind=8), intent(in) :: deltagx
    integer(kind=4), dimension(3), intent(in) :: idx

    real(kind=8) :: fornow,rd
    real(kind=8) :: xcoord
    integer(kind=4) :: ix,il

    ix = idx(1)

    xcoord = REAL(ix-1,kind=8)*deltagx

!   NEW ADDITION FOR REACTING CASE
    DO il = 1,nl-1
        IF( (xl(il) <= xcoord) .AND. (xl(il+1) > xcoord) )THEN

            fornow = (xcoord-xl(il)) / (xl(il+1)-xl(il))

            trun(OPS_ACC1(0,0,0)) = theta(il) + fornow*(theta(il+1)-theta(il))

        END IF
    END DO

END SUBROUTINE flamin_kernel_eqD

SUBROUTINE flamin_kernel_eqE(drun,store1,prin,trin)

    real(kind=8), dimension(1) :: drun
    real(kind=8), dimension(1), intent(in) :: store1

    real(kind=8), intent(in) :: prin,trin

    drun(OPS_ACC1(0,0,0)) = prin/(store1(OPS_ACC2(0,0,0))*trin)

END SUBROUTINE flamin_kernel_eqE

SUBROUTINE flamin_kernel_eqF(prun,drun,prin,u0,deltagx,deltagy,deltagz,angfrx,angfry,angfrz,idx)

    real(kind=8), dimension(1) :: prun
    real(kind=8), dimension(1), intent(in) :: drun
    real(kind=8), intent(in) :: prin,u0,deltagx,deltagy,deltagz
    real(kind=8), intent(in) :: angfrx, angfry, angfrz
    integer(kind=4), dimension(3), intent(in) :: idx

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

!   set pressure profile assuming constant density
    xrgmnt = 2.0_8*xrgmnt
    yrgmnt = 2.0_8*yrgmnt
    zrgmnt = 2.0_8*zrgmnt

    prun(OPS_ACC1(0,0,0)) = prin + ((drun(OPS_ACC2(0,0,0))*u0*u0)/16.0_8) * (COS(xrgmnt)+COS(yrgmnt)) &
                                   * (COS(zrgmnt)+2.0_8)

END SUBROUTINE flamin_kernel_eqF
