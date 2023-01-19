SUBROUTINE dfbydz_kernel_null(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE dfbydz_kernel_null

SUBROUTINE dfbydz_kernel_interior(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd, fdiffe

    fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
    fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
    fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))
    fdiffd = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,-4))
    fdiffe = functn(OPS_ACC1(0,0,5)) - functn(OPS_ACC1(0,0,-5))

    fderiv(OPS_ACC2(0,0,0)) = acoffz*fdiffa + bcoffz*fdiffb  &
          + ccoffz*fdiffc + dcoffz*fdiffd  &
          + ecoffz*fdiffe

END SUBROUTINE dfbydz_kernel_interior

SUBROUTINE dfbydz_kernel_lhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   lh point: 4th order one-sided
    fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
    fderiv(OPS_ACC2(0,0,0)) =  acof1z*fdiffa + bcof1z*fdiffb  &
          + ccof1z*fdiffc + dcof1z*fdiffd

END SUBROUTINE dfbydz_kernel_lhpoint_4th_onesided

SUBROUTINE dfbydz_kernel_lhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   lh point plus 1: 4th order mixed
    fdiffa = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
    fderiv(OPS_ACC2(0,0,0)) = acof2z*fdiffa + bcof2z*fdiffb  &
          + ccof2z*fdiffc + dcof2z*fdiffd

END SUBROUTINE dfbydz_kernel_lhpoint_4th_mixed

SUBROUTINE dfbydz_kernel_lhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb

!   lh point plus 2: 4th order centered
    fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
    fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
    fderiv(OPS_ACC2(0,0,0)) = acof3z*fdiffa + bcof3z*fdiffb

END SUBROUTINE dfbydz_kernel_lhpoint_4th_centered

SUBROUTINE dfbydz_kernel_lhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc

!   lh point plus 3: 6th order centered
    fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
    fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
    fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))
    fderiv(OPS_ACC2(0,0,0)) = acof4z*fdiffa + bcof4z*fdiffb  &
          + ccof4z*fdiffc

END SUBROUTINE dfbydz_kernel_lhpoint_6th_centered

SUBROUTINE dfbydz_kernel_lhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   lh point plus 4: 8th order centered
    fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
    fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
    fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))
    fdiffd = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,-4))
    fderiv(OPS_ACC2(0,0,0)) = acof5z*fdiffa + bcof5z*fdiffb  &
          + ccof5z*fdiffc + dcof5z*fdiffd

END SUBROUTINE dfbydz_kernel_lhpoint_8th_centered

SUBROUTINE dfbydz_kernel_rhpoint_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   rh point minus 4: 8th order centered
    fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
    fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
    fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))
    fdiffd = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,-4))
    fderiv(OPS_ACC2(0,0,0)) = acof5z*fdiffa + bcof5z*fdiffb  &
          + ccof5z*fdiffc + dcof5z*fdiffd

END SUBROUTINE dfbydz_kernel_rhpoint_8th_centered

SUBROUTINE dfbydz_kernel_rhpoint_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc

!   rh point minus 3: 6th order centered
    fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
    fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
    fdiffc = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,-3))
    fderiv(OPS_ACC2(0,0,0)) = acof4z*fdiffa + bcof4z*fdiffb  &
          + ccof4z*fdiffc

END SUBROUTINE dfbydz_kernel_rhpoint_6th_centered

SUBROUTINE dfbydz_kernel_rhpoint_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb

!   rh point minus 2: 4th order centered
    fdiffa = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))
    fdiffb = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))
    fderiv(OPS_ACC2(0,0,0)) = acof3z*fdiffa + bcof3z*fdiffb

END SUBROUTINE dfbydz_kernel_rhpoint_4th_centered

SUBROUTINE dfbydz_kernel_rhpoint_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   rh point minus 1: 4th order mixed
    fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,1))
    fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
    fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
    fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
    fderiv(OPS_ACC2(0,0,0)) =  acof2z*fdiffa + bcof2z*fdiffb  &
          + ccof2z*fdiffc + dcof2z*fdiffd

END SUBROUTINE dfbydz_kernel_rhpoint_4th_mixed

SUBROUTINE dfbydz_kernel_rhpoint_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd

!   rh point: 4th order one-sided
    fdiffa = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
    fdiffb = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
    fdiffc = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
    fdiffd = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))
    fderiv(OPS_ACC2(0,0,0)) = acof1z*fdiffa + bcof1z*fdiffb  &
          + ccof1z*fdiffc + dcof1z*fdiffd

END SUBROUTINE dfbydz_kernel_rhpoint_4th_onesided

SUBROUTINE dfbydz_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdelz

END SUBROUTINE dfbydz_kernel_scaling
