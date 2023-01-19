SUBROUTINE dfbydx_kernel_interior(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd, fdiffe

    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))
    fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(-4,0,0))
    fdiffe = functn(OPS_ACC1(5,0,0)) - functn(OPS_ACC1(-5,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acoffx*fdiffa + bcoffx*fdiffb  &
          + ccoffx*fdiffc + dcoffx*fdiffd  &
          + ecoffx*fdiffe

END SUBROUTINE dfbydx_kernel_interior

SUBROUTINE dfbydx_kernel_lhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   lh point: 4th order one-sided
    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(0,0,0))
    fderiv(OPS_ACC2(0,0,0)) =  acof1x*fdiffa + bcof1x*fdiffb  &
          + ccof1x*fdiffc + dcof1x*fdiffd

END SUBROUTINE dfbydx_kernel_lhpoint_4th_onesided

SUBROUTINE dfbydx_kernel_lhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   lh point plus 1: 4th order mixed
    fdiffa = functn(OPS_ACC1(-1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(0,0,0))
    fderiv(OPS_ACC2(0,0,0)) = acof2x*fdiffa + bcof2x*fdiffb  &
          + ccof2x*fdiffc + dcof2x*fdiffd

END SUBROUTINE dfbydx_kernel_lhpoint_4th_mixed

SUBROUTINE dfbydx_kernel_lhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb

!   lh point plus 2: 4th order centered
    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
    fderiv(OPS_ACC2(0,0,0)) = acof3x*fdiffa + bcof3x*fdiffb

END SUBROUTINE dfbydx_kernel_lhpoint_4th_centered

SUBROUTINE dfbydx_kernel_lhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc

!   lh point plus 3: 6th order centered
    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))
    fderiv(OPS_ACC2(0,0,0)) = acof4x*fdiffa + bcof4x*fdiffb  &
          + ccof4x*fdiffc

END SUBROUTINE dfbydx_kernel_lhpoint_6th_centered

SUBROUTINE dfbydx_kernel_lhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   lh point plus 4: 8th order centered
    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))
    fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(-4,0,0))
    fderiv(OPS_ACC2(0,0,0)) = acof5x*fdiffa + bcof5x*fdiffb  &
          + ccof5x*fdiffc + dcof5x*fdiffd

END SUBROUTINE dfbydx_kernel_lhpoint_8th_centered

SUBROUTINE dfbydx_kernel_rhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   rh point minus 4: 8th order centered
    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))
    fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(-4,0,0))
    fderiv(OPS_ACC2(0,0,0)) = acof5x*fdiffa + bcof5x*fdiffb  &
          + ccof5x*fdiffc + dcof5x*fdiffd

END SUBROUTINE dfbydx_kernel_rhpoint_8th_centered

SUBROUTINE dfbydx_kernel_rhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc

!   rh point minus 3: 6th order centered
    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(-3,0,0))
    fderiv(OPS_ACC2(0,0,0)) = acof4x*fdiffa + bcof4x*fdiffb  &
          + ccof4x*fdiffc

END SUBROUTINE dfbydx_kernel_rhpoint_6th_centered

SUBROUTINE dfbydx_kernel_rhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb

!   rh point minus 2: 4th order centered
    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))
    fderiv(OPS_ACC2(0,0,0)) = acof3x*fdiffa + bcof3x*fdiffb

END SUBROUTINE dfbydx_kernel_rhpoint_4th_centered

SUBROUTINE dfbydx_kernel_rhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   rh point minus 1: 4th order mixed
    fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(1,0,0))
    fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
    fderiv(OPS_ACC2(0,0,0)) =  acof2x*fdiffa + bcof2x*fdiffb  &
          + ccof2x*fdiffc + dcof2x*fdiffd

END SUBROUTINE dfbydx_kernel_rhpoint_4th_mixed

SUBROUTINE dfbydx_kernel_rhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   rh point: 4th order one-sided
    fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-2,0,0))
    fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-3,0,0))
    fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(-4,0,0))
    fderiv(OPS_ACC2(0,0,0)) = acof1x*fdiffa + bcof1x*fdiffb  &
          + ccof1x*fdiffc + dcof1x*fdiffd

END SUBROUTINE dfbydx_kernel_rhpoint_4th_onesided

SUBROUTINE dfbydx_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdelx

END SUBROUTINE dfbydx_kernel_scaling
