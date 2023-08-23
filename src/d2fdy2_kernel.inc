SUBROUTINE d2fdy2_kernel_null(fderiv)

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdy2_kernel_null

SUBROUTINE d2fdy2_kernel_main(functn, fderiv, nyglbl, nendyl, nendyr, nbound, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm,fdifem

    integer(kind=4), intent(in) :: nyglbl, nendyl, nendyr, nbound
    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4) :: jc, jstart, jfinis

    jstart = 1
    jfinis = nyglbl
    IF(nendyl == nbound) jstart = 6
    IF(nendyr == nbound) jfinis = nyglbl-5

    jc = idx(2)

    IF (jc >= jstart .and. jc <= jfinis) THEN
!       interior scheme
        fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
        fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
        fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))
        fdifdp = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-4,0))
        fdifep = functn(OPS_ACC1(0,5,0)) - functn(OPS_ACC1(0,0,0))
        fdifem = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-5,0))

        fderiv(OPS_ACC2(0,0,0)) = acofsy*(fdifap-fdifam) + bcofsy*(fdifbp-fdifbm)  &
                                + ccofsy*(fdifcp-fdifcm) + dcofsy*(fdifdp-fdifdm)  &
                                + ecofsy*(fdifep-fdifem)

    ELSE IF (jc == 1) THEN
!       lh point: 4th order one-sided
        fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,5,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs1y*fdifap + bcfs1y*fdifbp  &
                                + ccfs1y*fdifcp + dcfs1y*fdifdp  &
                                + ecfs1y*fdifep

    ELSE IF (jc == 2) THEN
!       lh point plus 1: 4th order mixed
        fdifap = functn(OPS_ACC1(0,-1,0)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs2y*fdifap + bcfs2y*fdifbp  &
                                + ccfs2y*fdifcp + dcfs2y*fdifdp  &
                                + ecfs2y*fdifep

    ELSE IF (jc == 3) THEN
!       lh point plus 2: 4th order centered
        fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
        fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs3y*(fdifap-fdifam) + bcfs3y*(fdifbp-fdifbm)

    ELSE IF (jc == 4) THEN
!       lh point plus 3: 6th order centered
        fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
        fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
        fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs4y*(fdifap-fdifam)  &
                                + bcfs4y*(fdifbp-fdifbm) + ccfs4y*(fdifcp-fdifcm)

    ELSE IF (jc == 5) THEN
!       lh point plus 4: 8th order centered
        fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
        fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
        fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))
        fdifdp = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-4,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs5y*(fdifap-fdifam)  &
                                + bcfs5y*(fdifbp-fdifbm) + ccfs5y*(fdifcp-fdifcm)  &
                                + dcfs5y*(fdifdp-fdifdm)

    ELSE IF (jc == nyglbl-4) THEN
!       rh point minus 4: 8th order centered
        fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
        fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
        fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))
        fdifdp = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-4,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs5y*(fdifap-fdifam)  &
                                + bcfs5y*(fdifbp-fdifbm) + ccfs5y*(fdifcp-fdifcm)  &
                                + dcfs5y*(fdifdp-fdifdm)

    ELSE IF (jc == nyglbl-3) THEN
!       rh point minus 3: 6th order centered
        fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
        fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
        fdifcp = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs4y*(fdifap-fdifam)  &
                                + bcfs4y*(fdifbp-fdifbm) + ccfs4y*(fdifcp-fdifcm)

    ELSE IF (jc == nyglbl-2) THEN
!       rh point minus 2: 4th order centered
        fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
        fdifbp = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs3y*(fdifap-fdifam) + bcfs3y*(fdifbp-fdifbm)

    ELSE IF (jc == nyglbl-1) THEN
!       rh point minus 1: 4th order mixed
        fdifap = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,-1,0)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,-2,0)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,-3,0)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,-4,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs2y*fdifap + bcfs2y*fdifbp  &
                                + ccfs2y*fdifcp + dcfs2y*fdifdp  &
                                + ecfs2y*fdifep

    ELSE IF (jc == nyglbl) THEN
!       rh point: 4th order one-sided
        fdifap = functn(OPS_ACC1(0,-1,0)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,-2,0)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,-3,0)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,-4,0)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,-5,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs1y*fdifap + bcfs1y*fdifbp  &
                                + ccfs1y*fdifcp + dcfs1y*fdifdp  &
                                + ecfs1y*fdifep

    END IF

!   scaling
    fderiv(OPS_ACC2(0,0,0)) = fderiv(OPS_ACC2(0,0,0))*ovdly2

END SUBROUTINE d2fdy2_kernel_main
