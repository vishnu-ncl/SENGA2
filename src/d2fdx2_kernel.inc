SUBROUTINE d2fdx2_kernel_main(functn, fderiv, nxglbl, nendxl, nendxr, nbound, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm,fdifem

    integer(kind=4), intent(in) :: nxglbl, nendxl, nendxr, nbound
    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4) :: ic, istart, ifinis

    istart = 1
    ifinis = nxglbl
    IF(nendxl == nbound) istart = 6
    IF(nendxr == nbound) ifinis = nxglbl-5

    ic = idx(1)

    IF (ic >= istart .and. ic <= ifinis) THEN
!       interior scheme
        fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
        fdifdp = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-4,0,0))
        fdifep = functn(OPS_ACC1(5,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifem = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-5,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acofsx*(fdifap-fdifam) + bcofsx*(fdifbp-fdifbm)  &
                                + ccofsx*(fdifcp-fdifcm) + dcofsx*(fdifdp-fdifdm)  &
                                + ecofsx*(fdifep-fdifem)

    ELSE IF (ic == 1) THEN
!       lh point: 4th order one-sided
        fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(5,0,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs1x*fdifap + bcfs1x*fdifbp  &
                                + ccfs1x*fdifcp + dcfs1x*fdifdp  &
                                + ecfs1x*fdifep

    ELSE IF (ic == 2) THEN
!       lh point plus 1: 4th order mixed
        fdifap = functn(OPS_ACC1(-1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs2x*fdifap + bcfs2x*fdifbp  &
                                + ccfs2x*fdifcp + dcfs2x*fdifdp  &
                                + ecfs2x*fdifep

    ELSE IF (ic == 3) THEN
!       lh point plus 2: 4th order centered
        fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs3x*(fdifap-fdifam) + bcfs3x*(fdifbp-fdifbm)

    ELSE IF (ic == 4) THEN
!       lh point plus 3: 6th order centered
        fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs4x*(fdifap-fdifam)  &
                                + bcfs4x*(fdifbp-fdifbm) + ccfs4x*(fdifcp-fdifcm)

    ELSE IF (ic == 5) THEN
!       lh point plus 4: 8th order centered
        fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
        fdifdp = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-4,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs5x*(fdifap-fdifam)  &
                                + bcfs5x*(fdifbp-fdifbm) + ccfs5x*(fdifcp-fdifcm)  &
                                + dcfs5x*(fdifdp-fdifdm)

    ELSE IF (ic == nxglbl-4) THEN
!       rh point minus 4: 8th order centered
        fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
        fdifdp = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-4,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs5x*(fdifap-fdifam)  &
                                + bcfs5x*(fdifbp-fdifbm) + ccfs5x*(fdifcp-fdifcm)  &
                                + dcfs5x*(fdifdp-fdifdm)

    ELSE IF (ic == nxglbl-3) THEN
!       rh point minus 3: 6th order centered
        fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdifcp = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs4x*(fdifap-fdifam)  &
                                + bcfs4x*(fdifbp-fdifbm) + ccfs4x*(fdifcp-fdifcm)

    ELSE IF (ic == nxglbl-2) THEN
!       rh point minus 2: 4th order centered
        fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdifbp = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs3x*(fdifap-fdifam) + bcfs3x*(fdifbp-fdifbm)

    ELSE IF (ic == nxglbl-1) THEN
!       rh point minus 1: 4th order mixed
        fdifap = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(-1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(-2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(-3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(-4,0,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs2x*fdifap + bcfs2x*fdifbp  &
                                + ccfs2x*fdifcp + dcfs2x*fdifdp  &
                                + ecfs2x*fdifep

    ELSE IF (ic == nxglbl) THEN
!       rh point: 4th order one-sided
        fdifap = functn(OPS_ACC1(-1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(-2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(-3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(-4,0,0)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(-5,0,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs1x*fdifap + bcfs1x*fdifbp  &
                                + ccfs1x*fdifcp + dcfs1x*fdifdp  &
                                + ecfs1x*fdifep

    END IF

!   scaling
    fderiv(OPS_ACC2(0,0,0)) = fderiv(OPS_ACC2(0,0,0))*ovdlx2

END SUBROUTINE d2fdx2_kernel_main
