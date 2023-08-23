SUBROUTINE d2fdz2_kernel_null(fderiv)

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdz2_kernel_null

SUBROUTINE d2fdz2_kernel_main(functn, fderiv, nzglbl, nendzl, nendzr, nbound, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm,fdifem

    integer(kind=4), intent(in) :: nzglbl, nendzl, nendzr, nbound
    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4) :: kc, kstart, kfinis

    kstart = 1
    kfinis = nzglbl
    IF(nendzl == nbound) kstart = 6
    IF(nendzr == nbound) kfinis = nzglbl-5

    kc = idx(3)

    IF (kc >= kstart .and. kc <= kfinis) THEN
!       interior scheme
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
        fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))
        fdifep = functn(OPS_ACC1(0,0,5)) - functn(OPS_ACC1(0,0,0))
        fdifem = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-5))

        fderiv(OPS_ACC2(0,0,0)) = acofsz*(fdifap-fdifam) + bcofsz*(fdifbp-fdifbm)  &
                                + ccofsz*(fdifcp-fdifcm) + dcofsz*(fdifdp-fdifdm)  &
                                + ecofsz*(fdifep-fdifem)

    ELSE IF (kc == 1) THEN
!       lh point: 4th order one-sided
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,0,5)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs1z*fdifap + bcfs1z*fdifbp  &
                                + ccfs1z*fdifcp + dcfs1z*fdifdp  &
                                + ecfs1z*fdifep

    ELSE IF (kc == 2) THEN
!       lh point plus 1: 4th order mixed
        fdifap = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs2z*fdifap + bcfs2z*fdifbp  &
                                + ccfs2z*fdifcp + dcfs2z*fdifdp  &
                                + ecfs2z*fdifep

    ELSE IF (kc == 3) THEN
!       lh point plus 2: 4th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acfs3z*(fdifap-fdifam) + bcfs3z*(fdifbp-fdifbm)

    ELSE IF (kc == 4) THEN
!       lh point plus 3: 6th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acfs4z*(fdifap-fdifam)  &
                                + bcfs4z*(fdifbp-fdifbm) + ccfs4z*(fdifcp-fdifcm)

    ELSE IF (kc == 5) THEN
!       lh point plus 4: 8th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
        fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acfs5z*(fdifap-fdifam)  &
                                + bcfs5z*(fdifbp-fdifbm) + ccfs5z*(fdifcp-fdifcm)  &
                                + dcfs5z*(fdifdp-fdifdm)

    ELSE IF (kc == nzglbl-4) THEN
!       rh point minus 4: 8th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
        fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acfs5z*(fdifap-fdifam)  &
                                + bcfs5z*(fdifbp-fdifbm) + ccfs5z*(fdifcp-fdifcm)  &
                                + dcfs5z*(fdifdp-fdifdm)

    ELSE IF (kc == nzglbl-3) THEN
!       rh point minus 3: 6th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acfs4z*(fdifap-fdifam)  &
                                + bcfs4z*(fdifbp-fdifbm) + ccfs4z*(fdifcp-fdifcm)

    ELSE IF (kc == nzglbl-2) THEN
!       rh point minus 2: 4th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acfs3z*(fdifap-fdifam) + bcfs3z*(fdifbp-fdifbm)

    ELSE IF (kc == nzglbl-1) THEN
!       rh point minus 1: 4th order mixed
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,0,-2)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,0,-3)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,0,-4)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs2z*fdifap + bcfs2z*fdifbp  &
                                + ccfs2z*fdifcp + dcfs2z*fdifdp  &
                                + ecfs2z*fdifep

    ELSE IF (kc == nzglbl) THEN
!       rh point: 4th order one-sided
        fdifap = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,0,-2)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,0,-3)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,0,-4)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,0,-5)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs1z*fdifap + bcfs1z*fdifbp  &
                                + ccfs1z*fdifcp + dcfs1z*fdifdp  &
                                + ecfs1z*fdifep

    END IF

!   scaling
    fderiv(OPS_ACC2(0,0,0)) = fderiv(OPS_ACC2(0,0,0))*ovdlz2

END SUBROUTINE d2fdz2_kernel_main
