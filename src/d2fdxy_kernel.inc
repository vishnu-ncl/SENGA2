SUBROUTINE d2fdxy_kernel_null(fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_dp

END SUBROUTINE d2fdxy_kernel_null

SUBROUTINE d2fdxy_kernel_interior(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd,fdiffe
    
    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0)) &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0)) &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0)) &
           - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))
    fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0)) &
           - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,-4,0))
    fdiffe = functn(OPS_ACC1(5,5,0)) - functn(OPS_ACC1(5,-5,0)) &
           - functn(OPS_ACC1(-5,5,0)) + functn(OPS_ACC1(-5,-5,0))

    fderiv(OPS_ACC2(0,0,0)) = acofxy*fdiffa + bcofxy*fdiffb  &
          + ccofxy*fdiffc + dcofxy*fdiffd  &
          + ecofxy*fdiffe

END SUBROUTINE d2fdxy_kernel_interior

SUBROUTINE d2fdxy_kernel_lh_xdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0)) &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(1,-2,0)) &
                    - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(2,-1,0)) &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0)) &
                    - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(3,1,0)) - functn(OPS_ACC1(3,-1,0)) &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(3,2,0)) - functn(OPS_ACC1(3,-2,0)) &
                   - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(4,1,0)) - functn(OPS_ACC1(4,-1,0)) &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(4,2,0)) - functn(OPS_ACC1(4,-2,0)) &
                    - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_xdir_4th_onesided

SUBROUTINE d2fdxy_kernel_lh_xdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(-1,1,0)) - functn(OPS_ACC1(-1,-1,0)) &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(-1,2,0)) - functn(OPS_ACC1(-1,-2,0)) &
                    - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0)) &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(1,-2,0)) &
                    - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(2,-1,0)) &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0)) &
                    - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(3,1,0)) - functn(OPS_ACC1(3,-1,0)) &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(3,2,0)) - functn(OPS_ACC1(3,-2,0)) &
                    - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
          + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_xdir_4th_mixed

SUBROUTINE d2fdxy_kernel_lh_xdir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
          - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
          - functn(OPS_ACC1(-2,2,0))  + functn(OPS_ACC1(-2,-2,0))
    fderiv(OPS_ACC2(0,0,0))  = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_lh_xdir_4th_centered

SUBROUTINE d2fdxy_kernel_lh_xdir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0))  + functn(OPS_ACC1(-3,-3,0))
    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc

END SUBROUTINE d2fdxy_kernel_lh_xdir_6th_centered

SUBROUTINE d2fdxy_kernel_lh_xdir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))
    fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
           - functn(OPS_ACC1(-4,4,0))  + functn(OPS_ACC1(-4,-4,0))
    fderiv(OPS_ACC2(0,0,0)) = acf5xy*fdiffa + bcf5xy*fdiffb  &
          + ccf5xy*fdiffc + dcf5xy*fdiffd
    
END SUBROUTINE d2fdxy_kernel_lh_xdir_8th_centered

SUBROUTINE d2fdxy_kernel_rh_xdir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))
    fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
           - functn(OPS_ACC1(-4,4,0))  + functn(OPS_ACC1(-4,-4,0))
    fderiv(OPS_ACC2(0,0,0)) = acf5xy*fdiffa + bcf5xy*fdiffb  &
          + ccf5xy*fdiffc + dcf5xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_xdir_8th_centered

SUBROUTINE d2fdxy_kernel_rh_xdir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0))  + functn(OPS_ACC1(-3,-3,0))
    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc

END SUBROUTINE d2fdxy_kernel_rh_xdir_6th_centered

SUBROUTINE d2fdxy_kernel_rh_xdir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
          - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
          - functn(OPS_ACC1(-2,2,0))  + functn(OPS_ACC1(-2,-2,0))
    fderiv(OPS_ACC2(0,0,0))  = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_rh_xdir_4th_centered

SUBROUTINE d2fdxy_kernel_rh_xdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(1,1,0)) + functn(OPS_ACC1(1,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(1,2,0)) + functn(OPS_ACC1(1,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(-1,2,0)) + functn(OPS_ACC1(-1,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-2,1,0)) + functn(OPS_ACC1(-2,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-3,1,0)) + functn(OPS_ACC1(-3,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(-3,2,0)) + functn(OPS_ACC1(-3,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
          + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_xdir_4th_mixed

SUBROUTINE d2fdxy_kernel_rh_xdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(-1,2,0)) + functn(OPS_ACC1(-1,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-2,1,0)) + functn(OPS_ACC1(-2,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-3,1,0)) + functn(OPS_ACC1(-3,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                   - functn(OPS_ACC1(-3,2,0)) + functn(OPS_ACC1(-3,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-4,1,0)) + functn(OPS_ACC1(-4,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(-4,2,0)) + functn(OPS_ACC1(-4,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_xdir_4th_onesided



