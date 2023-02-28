SUBROUTINE dfbydy_kernel_null(fderiv)
    use data_types
    implicit none
    
    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE dfbydy_kernel_null

SUBROUTINE dfbydy_kernel_main(functn, fderiv, nyglbl, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd, fdiffe

    integer, intent(in) :: nyglbl
    integer, dimension(3), intent(in) :: idx
    integer :: jc

    jc = idx(2)

    IF (jc >= 6 .and. jc <= nyglbl-5) THEN
!       interior scheme
        fdiffa = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))
        fdiffb = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))
        fdiffc = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,-3,0))
        fdiffd = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,-4,0))
        fdiffe = functn(OPS_ACC1(0,5,0)) - functn(OPS_ACC1(0,-5,0))

        fderiv(OPS_ACC2(0,0,0)) = acoffy*fdiffa + bcoffy*fdiffb  &
                                + ccoffy*fdiffc + dcoffy*fdiffd  &
                                + ecoffy*fdiffe

    ELSE IF (jc == 1) THEN
!       lh point: 4th order one-sided
        fdiffa = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdiffb = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdiffc = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))
        fdiffd = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof1y*fdiffa + bcof1y*fdiffb  &
                                + ccof1y*fdiffc + dcof1y*fdiffd

    ELSE IF (jc == 2) THEN
!       lh point plus 1: 4th order mixed
        fdiffa = functn(OPS_ACC1(0,-1,0)) - functn(OPS_ACC1(0,0,0))
        fdiffb = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))
        fdiffc = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))
        fdiffd = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acof2y*fdiffa + bcof2y*fdiffb  &
                                + ccof2y*fdiffc + dcof2y*fdiffd

    ELSE IF (jc == 3) THEN
!       lh point plus 2: 4th order centered
        fdiffa = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))
        fdiffb = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))

        fderiv(OPS_ACC2(0,0,0)) = acof3y*fdiffa + bcof3y*fdiffb

    ELSE IF (jc == 4) THEN
!       lh point plus 3: 6th order centered
        fdiffa = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))
        fdiffb = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))
        fdiffc = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,-3,0))

        fderiv(OPS_ACC2(0,0,0)) = acof4y*fdiffa + bcof4y*fdiffb  &
                                + ccof4y*fdiffc

    ELSE IF (jc == 5) THEN
!       lh point plus 4: 8th order centered
        fdiffa = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))
        fdiffb = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))
        fdiffc = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,-3,0))
        fdiffd = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,-4,0))

        fderiv(OPS_ACC2(0,0,0)) = acof5y*fdiffa + bcof5y*fdiffb  &
                                + ccof5y*fdiffc + dcof5y*fdiffd

    ELSE IF (jc == nyglbl-4) THEN
!       rh point minus 4: 8th order centered
        fdiffa = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))
        fdiffb = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))
        fdiffc = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,-3,0))
        fdiffd = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,-4,0))

        fderiv(OPS_ACC2(0,0,0)) = acof5y*fdiffa + bcof5y*fdiffb  &
                                + ccof5y*fdiffc + dcof5y*fdiffd

    ELSE IF (jc == nyglbl-3) THEN
!       rh point minus 3: 6th order centered
        fdiffa = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))
        fdiffb = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))
        fdiffc = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,-3,0))

        fderiv(OPS_ACC2(0,0,0)) = acof4y*fdiffa + bcof4y*fdiffb  &
                                + ccof4y*fdiffc

    ELSE IF (jc == nyglbl-2) THEN
!       rh point minus 2: 4th order centered
        fdiffa = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))
        fdiffb = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))

        fderiv(OPS_ACC2(0,0,0)) = acof3y*fdiffa + bcof3y*fdiffb

    ELSE IF (jc == nyglbl-1) THEN
!       rh point minus 1: 4th order mixed
        fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,1,0))
        fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
        fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
        fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))

        fderiv(OPS_ACC2(0,0,0)) = acof2y*fdiffa + bcof2y*fdiffb  &
                                + ccof2y*fdiffc + dcof2y*fdiffd

    ELSE IF (jc == nyglbl) THEN
!       rh point: 4th order one-sided
        fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-1,0))
        fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-2,0))
        fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-3,0))
        fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,-4,0))

        fderiv(OPS_ACC2(0,0,0)) = acof1y*fdiffa + bcof1y*fdiffb  &
                                + ccof1y*fdiffc + dcof1y*fdiffd

    END IF

!   scaling
    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdely

END SUBROUTINE dfbydy_kernel_main