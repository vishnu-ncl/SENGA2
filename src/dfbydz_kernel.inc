SUBROUTINE dfbydz_kernel_null(fderiv)

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE dfbydz_kernel_null

SUBROUTINE dfbydz_kernel_main(functn, fderiv, nzglbl, nendzl, nendzr, nbound, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd, fdiffe

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
        fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
        fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
        fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))
        fdiffd = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,-4))
        fdiffe = functn(OPS_ACC1(0,0,5)) - functn(OPS_ACC1(0,0,-5))

        fderiv(OPS_ACC2(0,0,0)) = acoffz*fdiffa + bcoffz*fdiffb  &
                                + ccoffz*fdiffc + dcoffz*fdiffd  &
                                + ecoffz*fdiffe

    ELSE IF (kc == 1) THEN
!       lh point: 4th order one-sided
        fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdiffd = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof1z*fdiffa + bcof1z*fdiffb  &
                                + ccof1z*fdiffc + dcof1z*fdiffd

    ELSE IF (kc == 2) THEN
!       lh point plus 1: 4th order mixed
        fdiffa = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
        fdiffb = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdiffc = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdiffd = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof2z*fdiffa + bcof2z*fdiffb  &
                                + ccof2z*fdiffc + dcof2z*fdiffd

    ELSE IF (kc == 3) THEN
!       lh point plus 2: 4th order centered
        fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
        fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acof3z*fdiffa + bcof3z*fdiffb

    ELSE IF (kc == 4) THEN
!       lh point plus 3: 6th order centered
        fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
        fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
        fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acof4z*fdiffa + bcof4z*fdiffb  &
                                + ccof4z*fdiffc

    ELSE IF (kc == 5) THEN
!       lh point plus 4: 8th order centered
        fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
        fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
        fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))
        fdiffd = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acof5z*fdiffa + bcof5z*fdiffb  &
                                + ccof5z*fdiffc + dcof5z*fdiffd

    ELSE IF (kc == nzglbl-4) THEN
!       rh point minus 4: 8th order centered
        fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
        fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
        fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))
        fdiffd = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acof5z*fdiffa + bcof5z*fdiffb  &
                                + ccof5z*fdiffc + dcof5z*fdiffd

    ELSE IF (kc == nzglbl-3) THEN
!       rh point minus 3: 6th order centered
        fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
        fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
        fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acof4z*fdiffa + bcof4z*fdiffb  &
                                + ccof4z*fdiffc

    ELSE IF (kc == nzglbl-2) THEN
!       rh point minus 2: 4th order centered
        fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
        fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acof3z*fdiffa + bcof3z*fdiffb

    ELSE IF (kc == nzglbl-1) THEN
!       rh point minus 1: 4th order mixed
        fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,1))
        fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acof2z*fdiffa + bcof2z*fdiffb  &
                                + ccof2z*fdiffc + dcof2z*fdiffd

    ELSE IF (kc == nzglbl) THEN
!       rh point: 4th order one-sided
        fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
        fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acof1z*fdiffa + bcof1z*fdiffb  &
                                + ccof1z*fdiffc + dcof1z*fdiffd

    END IF

!   scaling
    fderiv(OPS_ACC2(0,0,0)) = fderiv(OPS_ACC2(0,0,0))*ovdelz

END SUBROUTINE dfbydz_kernel_main
