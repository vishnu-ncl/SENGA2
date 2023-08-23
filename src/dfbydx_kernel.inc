SUBROUTINE dfbydx_kernel_main(functn, fderiv, nxglbl, nendxl, nendxr, nbound, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd, fdiffe

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
        fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))
        fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(-4,0,0))
        fdiffe = functn(OPS_ACC1(5,0,0)) - functn(OPS_ACC1(-5,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acoffx*fdiffa + bcoffx*fdiffb  &
                                + ccoffx*fdiffc + dcoffx*fdiffd  &
                                + ecoffx*fdiffe

    ELSE IF (ic == 1) THEN
!       lh point: 4th order one-sided
        fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
        fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) =  acof1x*fdiffa + bcof1x*fdiffb  &
                                 + ccof1x*fdiffc + dcof1x*fdiffd

    ELSE IF (ic == 2) THEN
!       lh point plus 1: 4th order mixed
        fdiffa = functn(OPS_ACC1(-1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdiffb = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
        fdiffc = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
        fdiffd = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof2x*fdiffa + bcof2x*fdiffb  &
                                + ccof2x*fdiffc + dcof2x*fdiffd

    ELSE IF (ic == 3) THEN
!       lh point plus 2: 4th order centered
        fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof3x*fdiffa + bcof3x*fdiffb

    ELSE IF (ic == 4) THEN
!       lh point plus 3: 6th order centered
        fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof4x*fdiffa + bcof4x*fdiffb  &
                                + ccof4x*fdiffc

    ELSE IF (ic == 5) THEN
!       lh point plus 4: 8th order centered
        fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))
        fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(-4,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof5x*fdiffa + bcof5x*fdiffb  &
                                + ccof5x*fdiffc + dcof5x*fdiffd

    ELSE IF (ic == nxglbl-4) THEN
!       rh point minus 4: 8th order centered
        fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))
        fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(-4,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof5x*fdiffa + bcof5x*fdiffb  &
                                + ccof5x*fdiffc + dcof5x*fdiffd

    ELSE IF (ic == nxglbl-3) THEN
!       rh point minus 3: 6th order centered
        fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof4x*fdiffa + bcof4x*fdiffb  &
                                + ccof4x*fdiffc

    ELSE IF (ic == nxglbl-2) THEN
!       rh point minus 2: 4th order centered
        fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof3x*fdiffa + bcof3x*fdiffb

    ELSE IF (ic == nxglbl-1) THEN
!       rh point minus 1: 4th order mixed
        fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(1,0,0))
        fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))

        fderiv(OPS_ACC2(0,0,0)) =  acof2x*fdiffa + bcof2x*fdiffb  &
                                 + ccof2x*fdiffc + dcof2x*fdiffd

    ELSE IF (ic == nxglbl) THEN
!       rh point: 4th order one-sided
        fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
        fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
        fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
        fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-4,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof1x*fdiffa + bcof1x*fdiffb  &
                                + ccof1x*fdiffc + dcof1x*fdiffd

    END IF

!   scaling
    fderiv(OPS_ACC2(0,0,0)) = fderiv(OPS_ACC2(0,0,0))*ovdelx

END SUBROUTINE dfbydx_kernel_main
