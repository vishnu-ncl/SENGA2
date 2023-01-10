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

SUBROUTINE d2fdxy_kernel_lh_ydir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(-2,1,0))  &
                  - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(-1,2,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
                  - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,3,0)) - functn(OPS_ACC1(-1,3,0))  &
                   - functn(OPS_ACC1(1,0,0)  + functn(OPS_ACC1(-1,0,0))  &
           + bcofx1*(functn(OPS_ACC1(2,3,0)) - functn(OPS_ACC1(-2,3,0))  &
                  - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,4,0)) - functn(OPS_ACC1(-1,4,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,4,0)) - functn(OPS_ACC1(-2,4,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_4th_onesided

SUBROUTINE d2fdxy_kernel_lh_ydir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,-1,0))  - functn(OPS_ACC1(-1,-1,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,-1,0))  - functn(OPS_ACC1(-2,-1,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(-2,1,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(-1,2,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,3,0)) - functn(OPS_ACC1(-1,3,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,3,0)) - functn(OPS_ACC1(-2,3,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
          + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_4th_mixed

SUBROUTINE d2fdxy_kernel_lh_ydir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0))  + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_lh_ydir_4th_centered

SUBROUTINE d2fdxy_kernel_lh_ydir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(-3,3,0))  &
           - functn(OPS_ACC1(3,-3,0))  + functn(OPS_ACC1(-3,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc

END SUBROUTINE d2fdxy_kernel_lh_ydir_6th_centered

SUBROUTINE d2fdxy_kernel_lh_ydir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(-3,3,0))  &
           - functn(OPS_ACC1(3,-3,0)) + functn(OPS_ACC1(-3,-3,0))
    fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(-4,4,0))  &
           - functn(OPS_ACC1(4,-4,0))  + functn(OPS_ACC1(-4,-4,0))

    fderiv(OPS_ACC2(0,0,0)) = acf5xy*fdiffa + bcf5xy*fdiffb  &
          + ccf5xy*fdiffc + dcf5xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_8th_centered

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqA(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,0,0))  &
           - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,0,0))  &
           - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,0,0))  &
           - functn(OPS_ACC1(0,3,0))  + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,0,0))  &
           - functn(OPS_ACC1(0,4,0))  + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1xy*fdiffa + bcc1xy*fdiffb  &
          + ccc1xy*fdiffc + dcc1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqA

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqB(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(-1,-1,0))   - functn(OPS_ACC1(-1,0,0))  &
           - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,0,0))  &
           - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,0,0))  &
           - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,0,0))  &
           - functn(OPS_ACC1(0,3,0)) + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc2xy*fdiffa + bcc2xy*fdiffb  &
          + ccc2xy*fdiffc + dcc2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqB

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqC(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xy*(functn(OPS_ACC1(1,-1,0))  - functn(OPS_ACC1(1,0,0))  &
                   - functn(OPS_ACC1(0,-1,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,0,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(1,0,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(1,3,0)) - functn(OPS_ACC1(1,0,0))  &
                   - functn(OPS_ACC1(0,3,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2xy*(functn(OPS_ACC1(2,-1,0))  - functn(OPS_ACC1(2,0,0))  &
                   - functn(OPS_ACC1(0,-1,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(2,0,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,0,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(2,3,0)) - functn(OPS_ACC1(2,0,0))  &
                   - functn(OPS_ACC1(0,3,0))  + functn(OPS_ACC1(0,0,0)))

    fdiffc = acf2xy*(functn(OPS_ACC1(3,-1,0))  - functn(OPS_ACC1(3,0,0))  &
                   - functn(OPS_ACC1(0,-1,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(3,1,0)) - functn(OPS_ACC1(3,0,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(3,2,0)) - functn(OPS_ACC1(3,0,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,0,0))  &
                   - functn(OPS_ACC1(0,3,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2xy*(functn(OPS_ACC1(4,-1,0))  - functn(OPS_ACC1(4,0,0))  &
                   - functn(OPS_ACC1(0,-1,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(4,1,0)) - functn(OPS_ACC1(4,0,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(4,2,0)) - functn(OPS_ACC1(4,0,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(4,3,0)) - functn(OPS_ACC1(4,0,0))  &
                   - functn(OPS_ACC1(0,3,0))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqC

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqD(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xy*(functn(OPS_ACC1(-1,1,0))  - functn(OPS_ACC1(0,1,0))  &
                   - functn(OPS_ACC1(-1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(0,1,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(0,1,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,1,0)) - functn(OPS_ACC1(0,1,0))  &
                   - functn(OPS_ACC1(3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2xy*(functn(OPS_ACC1(-1,2,0))  - functn(OPS_ACC1(0,2,0))  &
                   - functn(OPS_ACC1(-1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(0,2,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(0,2,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,2,0)) - functn(OPS_ACC1(0,2,0))  &
                   - functn(OPS_ACC1(3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffc = acf2xy*(functn(OPS_ACC1(-1,3,0))  - functn(OPS_ACC1(0,3,0))  &
                   - functn(OPS_ACC1(-1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,3,0)) - functn(OPS_ACC1(0,3,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,3,0)) - functn(OPS_ACC1(0,3,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(0,3,0))  &
                   - functn(OPS_ACC1(3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2xy*(functn(OPS_ACC1(-1,4,0))  - functn(OPS_ACC1(0,4,0))  &
                   - functn(OPS_ACC1(-1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,4,0)) - functn(OPS_ACC1(0,4,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,4,0)) - functn(OPS_ACC1(0,4,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,4,0)) - functn(OPS_ACC1(0,4,0))  &
                   - functn(OPS_ACC1(3,0,0))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqD

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqE(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(-2,1,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(-1,2,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,3,0)) - functn(OPS_ACC1(-1,3,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,3,0)) - functn(OPS_ACC1(-2,3,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,4,0)) - functn(OPS_ACC1(-1,4,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,4,0)) - functn(OPS_ACC1(-2,4,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqE

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqF(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,-1,0))  - functn(OPS_ACC1(-1,-1,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,-1,0))  - functn(OPS_ACC1(-2,-1,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(-2,1,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(-1,2,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,3,0)) - functn(OPS_ACC1(-1,3,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,3,0)) - functn(OPS_ACC1(-2,3,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqF

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqG(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(1,-2,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(2,-1,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(3,1,0)) - functn(OPS_ACC1(3,-1,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(3,2,0)) - functn(OPS_ACC1(3,-2,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(4,1,0)) - functn(OPS_ACC1(4,-1,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(4,2,0)) - functn(OPS_ACC1(4,-2,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqG

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqH(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(-1,1,0))  - functn(OPS_ACC1(-1,-1,0))  &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(-1,2,0))  - functn(OPS_ACC1(-1,-2,0))  &
                   - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(1,-2,0))  &
                   - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(2,-1,0))  &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
                   - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(3,1,0)) - functn(OPS_ACC1(3,-1,0))  &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(3,2,0)) - functn(OPS_ACC1(3,-2,0))  &
                   - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqH

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqL(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(OPS_ACC1(-2,0,0))
    fdiffc = functn(OPS_ACC1(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(OPS_ACC1(-3,0,0))
    fdiffd = functn(OPS_ACC1(OPS_ACC1(0,4,0))  - functn(OPS_ACC1(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(OPS_ACC1(-4,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1xy*fdiffa + bcc1xy*fdiffb  &
          + ccc1xy*fdiffc + dcc1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqL

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqM(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,-1,0))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(1,-1,0))   + functn(OPS_ACC1(1,0,0))
    fdiffb = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,0,0))
    fdiffc = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,0,0))
    fdiffd = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc2xy*fdiffa + bcc2xy*fdiffb  &
          + ccc2xy*fdiffc + dcc2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqM

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqN(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xy*(functn(OPS_ACC1(0,-1,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,-1,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,2,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,3,0)) + functn(OPS_ACC1(-1,0,0)))
    fdiffb = acf2xy*(functn(OPS_ACC1(0,-1,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-2,-1,0))  + functn(OPS_ACC1(-2,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-2,1,0)) + functn(OPS_ACC1(-2,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-2,3,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acf2xy*(functn(OPS_ACC1(0,-1,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-3,-1,0))  + functn(OPS_ACC1(-3,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-3,1,0)) + functn(OPS_ACC1(-3,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-3,2,0)) + functn(OPS_ACC1(-3,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,0,0)))
    fdiffd = acf2xy*(functn(OPS_ACC1(0,-1,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-4,-1,0))  + functn(OPS_ACC1(-4,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-4,1,0)) + functn(OPS_ACC1(-4,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-4,2,0)) + functn(OPS_ACC1(-4,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-4,3,0)) + functn(OPS_ACC1(-4,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqN

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqO(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xy*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(1,1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(1,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(-1,1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(-2,1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(-3,1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-3,0,0)))
    fdiffb = acf2xy*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(1,2,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(1,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(-1,2,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(-2,2,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(-3,2,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-3,0,0)))
    fdiffc = acf2xy*(functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(1,3,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(1,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(-1,3,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(-2,3,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(-3,3,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-3,0,0)))
    fdiffd = acf2xy*(functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(1,4,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(1,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(-1,4,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(-2,4,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(-3,4,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-3,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqO

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqP(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(-2,1,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(-1,2,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,3,0)) - functn(OPS_ACC1(-1,3,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,3,0)) - functn(OPS_ACC1(-2,3,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,4,0)) - functn(OPS_ACC1(-1,4,0))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,4,0)) - functn(OPS_ACC1(-2,4,0))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqP

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqQ(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,-1,0))  - functn(OPS_ACC1(-1,-1,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,-1,0))  - functn(OPS_ACC1(-2,-1,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(-2,1,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(-1,2,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,3,0)) - functn(OPS_ACC1(-1,3,0))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,3,0)) - functn(OPS_ACC1(-2,3,0))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqQ

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqR(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-1,2,0)) + functn(OPS_ACC1(-1,-2,0)))

    fdiffb = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-2,1,0)) + functn(OPS_ACC1(-2,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0)))

    fdiffc = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-3,1,0)) + functn(OPS_ACC1(-3,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-3,2,0)) + functn(OPS_ACC1(-3,-2,0)))

    fdiffd = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-4,1,0)) + functn(OPS_ACC1(-4,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-4,2,0)) + functn(OPS_ACC1(-4,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqR

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqS(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(1,1,0))  + functn(OPS_ACC1(1,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(1,2,0))  + functn(OPS_ACC1(1,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-1,2,0)) + functn(OPS_ACC1(-1,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-2,1,0)) + functn(OPS_ACC1(-2,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-3,1,0)) + functn(OPS_ACC1(-3,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-3,2,0)) + functn(OPS_ACC1(-3,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqS

SUBROUTINE d2fdxy_kernel_rh_ydir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(-3,3,0))  &
           - functn(OPS_ACC1(3,-3,0)) + functn(OPS_ACC1(-3,-3,0))
    fdiffd = functn(OPS_ACC1(4,4,0))  - functn(OPS_ACC1(-4,4,0))  &
           - functn(OPS_ACC1(4,-4,0)) + functn(OPS_ACC1(-4,-4,0))

    fderiv(OPS_ACC2(0,0,0)) = acf5xy*fdiffa + bcf5xy*fdiffb  &
          + ccf5xy*fdiffc + dcf5xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_8th_centered

SUBROUTINE d2fdxy_kernel_rh_ydir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0))  - functn(OPS_ACC1(-3,3,0))  &
           - functn(OPS_ACC1(3,-3,0)) + functn(OPS_ACC1(-3,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc

END SUBROUTINE d2fdxy_kernel_rh_ydir_6th_centered

SUBROUTINE d2fdxy_kernel_rh_ydir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0))  - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_rh_ydir_4th_centered

SUBROUTINE d2fdxy_kernel_rh_ydir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,1,0))  + functn(OPS_ACC1(-1,1,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,1,0))  + functn(OPS_ACC1(-2,1,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-1,0)) + functn(OPS_ACC1(-2,-1,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-2,0)) + functn(OPS_ACC1(-1,-2,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-3,0)) + functn(OPS_ACC1(-1,-3,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-3,0)) + functn(OPS_ACC1(-2,-3,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
          + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_4th_mixed

SUBROUTINE d2fdxy_kernel_rh_ydir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-1,0)) + functn(OPS_ACC1(-2,-1,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-2,0)) + functn(OPS_ACC1(-1,-2,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-3,0)) + functn(OPS_ACC1(-1,-3,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-3,0)) + functn(OPS_ACC1(-2,-3,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-4,0)) + functn(OPS_ACC1(-1,-4,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-4,0)) + functn(OPS_ACC1(-2,-4,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_4th_onesided

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqA(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-1,0))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-2,0))
    fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-3,0))
    fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,-4,0))  &
           - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-4,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1xy*fdiffa + bcc1xy*fdiffb  &
          + ccc1xy*fdiffc + dcc1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqA

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqB(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(-1,0,0))  - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,1,0))
    fdiffb = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,-1,0))
    fdiffc = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,-2,0))
    fdiffd = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acc2xy*fdiffa + bcc2xy*fdiffb  &
          + ccc2xy*fdiffc + dcc2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqB

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqC(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xy*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,1,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,-1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + ccf2xy*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,-2,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-2,0)))  &
           + dcf2xy*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,-3,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-3,0)))
    fdiffb = acf2xy*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,1,0)))  &
           + bcf2xy*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,-1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,-2,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-2,0)))  &
           + dcf2xy*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,-3,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-3,0)))
    fdiffc = acf2xy*(functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,1,0)))  &
           + bcf2xy*(functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,-1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + ccf2xy*(functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,-2,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-2,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,-3,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-3,0)))
    fdiffd = acf2xy*(functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,1,0)))  &
           + bcf2xy*(functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,-1,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + ccf2xy*(functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,-2,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-2,0)))  &
           + dcf2xy*(functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,-3,0))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-3,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqC

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqD(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xy*(functn(OPS_ACC1(-1,0,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,-1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(2,-1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(3,-1,0)) + functn(OPS_ACC1(0,-1,0)))
    fdiffb = acf2xy*(functn(OPS_ACC1(-1,0,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,-2,0))  + functn(OPS_ACC1(0,-2,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(1,-2,0)) + functn(OPS_ACC1(0,-2,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(0,-2,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(3,-2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acf2xy*(functn(OPS_ACC1(-1,0,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,-3,0))  + functn(OPS_ACC1(0,-3,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(1,-3,0)) + functn(OPS_ACC1(0,-3,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(2,-3,0)) + functn(OPS_ACC1(0,-3,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(3,-3,0)) + functn(OPS_ACC1(0,-3,0)))
    fdiffd = acf2xy*(functn(OPS_ACC1(-1,0,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,-4,0))  + functn(OPS_ACC1(0,-4,0)))  &
           + bcf2xy*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(1,-4,0)) + functn(OPS_ACC1(0,-4,0)))  &
           + ccf2xy*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(2,-4,0)) + functn(OPS_ACC1(0,-4,0)))  &
           + dcf2xy*(functn(OPS_ACC1(3,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(3,-4,0)) + functn(OPS_ACC1(0,-4,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqD

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqE(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,1,0))  + functn(OPS_ACC1(-1,1,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,1,0))  + functn(OPS_ACC1(-2,1,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-1,0)) + functn(OPS_ACC1(-2,-1,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-2,0)) + functn(OPS_ACC1(-1,-2,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-3,0)) + functn(OPS_ACC1(-1,-3,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-3,0)) + functn(OPS_ACC1(-2,-3,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqE

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqF(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-1,0)) + functn(OPS_ACC1(-2,-1,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-2,0)) + functn(OPS_ACC1(-1,-2,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-3,0)) + functn(OPS_ACC1(-1,-3,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-3,0)) + functn(OPS_ACC1(-2,-3,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-4,0)) + functn(OPS_ACC1(-1,-4,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-4,0)) + functn(OPS_ACC1(-2,-4,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqF

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqG(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(1,-2,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(2,-1,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(3,1,0)) - functn(OPS_ACC1(3,-1,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(3,2,0)) - functn(OPS_ACC1(3,-2,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(4,1,0)) - functn(OPS_ACC1(4,-1,0))  &
                   - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(4,2,0)) - functn(OPS_ACC1(4,-2,0))  &
                   - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqG

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqH(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(-1,1,0))  - functn(OPS_ACC1(-1,-1,0))  &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(-1,2,0))  - functn(OPS_ACC1(-1,-2,0))  &
                   - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(1,2,0)) - functn(OPS_ACC1(1,-2,0))  &
                   - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(2,1,0)) - functn(OPS_ACC1(2,-1,0))  &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
                   - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(3,1,0)) - functn(OPS_ACC1(3,-1,0))  &
                   - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(3,2,0)) - functn(OPS_ACC1(3,-2,0))  &
                   - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqH

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqL(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(-1,-1,0)) - functn(OPS_ACC1(-1,0,0))  &
           - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(-2,-2,0)) - functn(OPS_ACC1(-2,0,0))  &
           - functn(OPS_ACC1(0,-2,0))  + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(-3,-3,0)) - functn(OPS_ACC1(-3,0,0))  &
           - functn(OPS_ACC1(0,-3,0))  + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(-4,-4,0)) - functn(OPS_ACC1(-4,0,0))  &
           - functn(OPS_ACC1(0,-4,0))  + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1xy*fdiffa + bcc1xy*fdiffb  &
          + ccc1xy*fdiffc + dcc1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqL

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqM(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,1,0))   - functn(OPS_ACC1(1,0,0))  &
           - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(-1,-1,0)) - functn(OPS_ACC1(-1,0,0))  &
           - functn(OPS_ACC1(0,-1,0)) + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(-2,-2,0)) - functn(OPS_ACC1(-2,0,0))  &
           - functn(OPS_ACC1(0,-2,0)) + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(-3,-3,0)) - functn(OPS_ACC1(-3,0,0))  &
           - functn(OPS_ACC1(0,-3,0)) + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc2xy*fdiffa + bcc2xy*fdiffb  &
          + ccc2xy*fdiffc + dcc2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqM

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqN(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xy*(functn(OPS_ACC1(-1,1,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(0,1,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(-1,-1,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(-1,-2,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(0,-2,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(-1,-3,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(0,-3,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2xy*(functn(OPS_ACC1(-2,1,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(0,1,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(-2,-1,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(-2,-2,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(0,-2,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(-2,-3,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(0,-3,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffc = acf2xy*(functn(OPS_ACC1(-3,1,0))  - functn(OPS_ACC1(-3,0,0))  &
                   - functn(OPS_ACC1(0,1,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(-3,-1,0)) - functn(OPS_ACC1(-3,0,0))  &
                   - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(-3,-2,0)) - functn(OPS_ACC1(-3,0,0))  &
                   - functn(OPS_ACC1(0,-2,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(-3,-3,0)) - functn(OPS_ACC1(-3,0,0))  &
                   - functn(OPS_ACC1(0,-3,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2xy*(functn(OPS_ACC1(-4,1,0))  - functn(OPS_ACC1(-4,0,0))  &
                   - functn(OPS_ACC1(0,1,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(-4,-1,0)) - functn(OPS_ACC1(-4,0,0))  &
                   - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(-4,-2,0)) - functn(OPS_ACC1(-4,0,0))  &
                   - functn(OPS_ACC1(0,-2,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(-4,-3,0)) - functn(OPS_ACC1(-4,0,0))  &
                   - functn(OPS_ACC1(0,-3,0))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqN

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqO(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xy*(functn(OPS_ACC1(1,-1,0))  - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(-1,-1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(-2,-1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(-3,-1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2xy*(functn(OPS_ACC1(1,-2,0))  - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(-1,-2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(-2,-2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(-3,-2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffc = acf2xy*(functn(OPS_ACC1(1,-3,0))  - functn(OPS_ACC1(0,-3,0))  &
                   - functn(OPS_ACC1(1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(-1,-3,0)) - functn(OPS_ACC1(0,-3,0))  &
                   - functn(OPS_ACC1(-1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(-2,-3,0)) - functn(OPS_ACC1(0,-3,0))  &
                   - functn(OPS_ACC1(-2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(-3,-3,0)) - functn(OPS_ACC1(0,-3,0))  &
                   - functn(OPS_ACC1(-3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2xy*(functn(OPS_ACC1(1,-4,0))  - functn(OPS_ACC1(0,-4,0))  &
                   - functn(OPS_ACC1(1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xy*(functn(OPS_ACC1(-1,-4,0)) - functn(OPS_ACC1(0,-4,0))  &
                   - functn(OPS_ACC1(-1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xy*(functn(OPS_ACC1(-2,-4,0)) - functn(OPS_ACC1(0,-4,0))  &
                   - functn(OPS_ACC1(-2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xy*(functn(OPS_ACC1(-3,-4,0)) - functn(OPS_ACC1(0,-4,0))  &
                   - functn(OPS_ACC1(-3,0,0))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqO

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqP(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,1,0))  + functn(OPS_ACC1(-1,1,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,1,0))  + functn(OPS_ACC1(-2,1,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-1,0)) + functn(OPS_ACC1(-2,-1,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-2,0)) + functn(OPS_ACC1(-1,-2,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-3,0)) + functn(OPS_ACC1(-1,-3,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-3,0)) + functn(OPS_ACC1(-2,-3,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqP

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqQ(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-1,0)) + functn(OPS_ACC1(-2,-1,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-2,0)) + functn(OPS_ACC1(-1,-2,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-3,0)) + functn(OPS_ACC1(-1,-3,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-3,0)) + functn(OPS_ACC1(-2,-3,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,-4,0)) + functn(OPS_ACC1(-1,-4,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,-4,0)) + functn(OPS_ACC1(-2,-4,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqQ

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqR(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-1,2,0)) + functn(OPS_ACC1(-1,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-2,1,0)) + functn(OPS_ACC1(-2,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-3,1,0)) + functn(OPS_ACC1(-3,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-3,2,0)) + functn(OPS_ACC1(-3,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-4,1,0)) + functn(OPS_ACC1(-4,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-4,2,0)) + functn(OPS_ACC1(-4,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqR

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqS(functn, fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: functn
    real(kind=dp), dimension(1) :: fderiv

    real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(1,1,0))  + functn(OPS_ACC1(1,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(1,2,0))  + functn(OPS_ACC1(1,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-1,2,0)) + functn(OPS_ACC1(-1,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-2,1,0)) + functn(OPS_ACC1(-2,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
                   - functn(OPS_ACC1(-3,1,0)) + functn(OPS_ACC1(-3,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
                   - functn(OPS_ACC1(-3,2,0)) + functn(OPS_ACC1(-3,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqS

SUBROUTINE d2fdxy_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdelx*ovdely

END SUBROUTINE d2fdxy_kernel_scaling
