SUBROUTINE flamin_kernel_set_velocity_tgv(urun, vrun, wrun, prun, prin, drin, u0, deltagx, deltagy, deltagz, angfrx, angfry, angfrz, idx)

    real(kind=8), dimension(1) :: urun, vrun, wrun, prun
    real(kind=8), intent(in) :: prin,drin,u0, deltagx, deltagy, deltagz
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

!  set taylor-green vortex velocity field
   urun(OPS_ACC1(0,0,0)) = u0 * dsin(xrgmnt) * dcos(yrgmnt) * dcos(zrgmnt)
   vrun(OPS_ACC2(0,0,0)) = -u0 * dcos(xrgmnt) * dsin(yrgmnt) * dcos(zrgmnt)
   wrun(OPS_ACC3(0,0,0)) = 0.0_8

!  set pressure profile assuming constant density
   xrgmnt = 2.0_8*xrgmnt
   yrgmnt = 2.0_8*yrgmnt
   zrgmnt = 2.0_8*zrgmnt

   prun(OPS_ACC4(0,0,0)) = prin+((drin*u0*u0)/16.0_8)*(dcos(xrgmnt)+dcos(yrgmnt))*(dcos(zrgmnt)+2.0_8)

END SUBROUTINE flamin_kernel_set_velocity_tgv

!A = B/var
SUBROUTINE flamin_kernel_eqA(trun, prun, rglocl)

    real(kind=8), dimension(1) :: trun
    real(kind=8), dimension(1), intent(in) :: prun
    real(kind=8), intent(in) :: rglocl

    trun(OPS_ACC1(0,0,0)) = prun(OPS_ACC2(0,0,0))/rglocl

END SUBROUTINE flamin_kernel_eqA
