SUBROUTINE d2fdxy_kernel_null(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdxy_kernel_null

!a3d_p550_m550_mixed_xy
SUBROUTINE d2fdxy_kernel_interior(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd,fdiffe
    
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

!a3d_p420_m020_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_xdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p320_m120_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_xdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p220_m220_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_xdir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
          - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
          - functn(OPS_ACC1(-2,2,0))  + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0))  = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_lh_xdir_4th_centered

!a3d_p330_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_xdir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0))  + functn(OPS_ACC1(-3,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc

END SUBROUTINE d2fdxy_kernel_lh_xdir_6th_centered

!a3d_p440_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_xdir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p440_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_xdir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p330_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_xdir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0))  + functn(OPS_ACC1(-3,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc

END SUBROUTINE d2fdxy_kernel_rh_xdir_6th_centered

!a3d_p220_m220_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_xdir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
          - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
          - functn(OPS_ACC1(-2,2,0))  + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0))  = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_rh_xdir_4th_centered

!a3d_p120_m320_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_xdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p020_m420_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_xdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p240_m200_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

END SUBROUTINE d2fdxy_kernel_lh_ydir_4th_onesided

!a3d_p230_m210_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p220_m220_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0))  + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_lh_ydir_4th_centered

!a3d_p330_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(-3,3,0))  &
           - functn(OPS_ACC1(3,-3,0))  + functn(OPS_ACC1(-3,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc

END SUBROUTINE d2fdxy_kernel_lh_ydir_6th_centered

!a3d_p440_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p440_p000_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqA(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p330_m110_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqB(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p430_m010_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqC(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p340_m100_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqD(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p240_m200_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqE(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p230_m210_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqF(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p420_m020_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqG(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p320_m120_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqH(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqI(functn, fderiv, fstora, fstorb, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv, fstora, fstorb
    real(kind=8) :: fdiffa, fdiffb

    integer, dimension(3), intent(in) :: idx
    integer :: is, js, ijs
    is = idx(1) -3
    js = idx(2) -3
    ijs = is*3+js +1

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

    fstora(OPS_ACC_MD3(ijs,0,0,0)) = fdiffa
    fstorb(OPS_ACC_MD4(ijs,0,0,0)) = fdiffb
    
END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqI

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqJ(functn, fderiv, fstora, fstorb, fstorc, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb
    real(kind=8), dimension(1) :: fderiv, fstorc
    real(kind=8) :: fdiffc

    integer, dimension(3), intent(in) :: idx
    integer :: is, js, ijs
    is = idx(1) -3
    js = idx(2) -3
    ijs = is*3+js +1

    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))
    
    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fstora(OPS_ACC_MD3(ijs,0,0,0)) + bcf4xy*fstorb(OPS_ACC_MD4(ijs,0,0,0))  &
              + ccf4xy*fdiffc

    is = is-1
    js = js-1
    ijs = is*2+js +1
    fstorc(OPS_ACC_MD5(ijs,0,0,0)) = fdiffc

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqJ

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqK(functn, fderiv, fstora, fstorb, fstorc, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb, fstorc
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffd

    integer, dimension(3), intent(in) :: idx
    integer :: is, js, ijs1, ijs2
    is = idx(1) -3
    js = idx(2) -3
    ijs1 = is*3+js +1

    is = is-1
    js = js-1
    ijs2 = is*2+js +1

    fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
           - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,-4,0))

    fderiv(OPS_ACC2(0,0,0)) = acf5xy*fstora(OPS_ACC_MD3(ijs1,0,0,0)) + bcf5xy*fstorb(OPS_ACC_MD4(ijs1,0,0,0))  &
                            + ccf5xy*fstorc(OPS_ACC_MD5(ijs2,0,0,0)) + dcf5xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqK

!a3d_p040_m400_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqL(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,0,0))
    fdiffc = functn(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,0,0))
    fdiffd = functn(OPS_ACC1(0,4,0))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1xy*fdiffa + bcc1xy*fdiffb  &
          + ccc1xy*fdiffc + dcc1xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqL

!a3d_p100_m300_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqM(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p030_m410_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqN(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p140_m300_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqO(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p240_m200_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqP(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p230_m210_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqQ(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p020_m420_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqR(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p120_m320_mixed_xy
SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqS(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqT(functn, fderiv, fstora, fstorb, nxglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv, fstora, fstorb
    real(kind=8) :: fdiffa, fdiffb

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nxglblm4

    integer :: is, js, ijs
    is = idx(1) -nxglblm4
    js = idx(2) -3
    ijs = is*3+js +1

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

    fstora(OPS_ACC_MD3(ijs,0,0,0)) = fdiffa
    fstorb(OPS_ACC_MD4(ijs,0,0,0)) = fdiffb

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqT

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqU(functn, fderiv, fstora, fstorb, fstorc, nxglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb
    real(kind=8), dimension(1) :: fderiv, fstorc
    real(kind=8) :: fdiffc

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nxglblm4

    integer :: is, js, ijs
    is = idx(1) -nxglblm4
    js = idx(2) -3
    ijs = is*3+js +1

    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fstora(OPS_ACC_MD3(ijs,0,0,0)) + bcf4xy*fstorb(OPS_ACC_MD4(ijs,0,0,0))  &
              + ccf4xy*fdiffc

    js = js-1
    ijs = is*2+js +1
    fstorc(OPS_ACC_MD5(ijs,0,0,0)) = fdiffc

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqU

SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqV(functn, fderiv, fstora, fstorb, fstorc, nxglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb, fstorc
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffd

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nxglblm4

    integer :: is, js, ijs1, ijs2
    is = idx(1) -nxglblm4
    js = idx(2) -3
    ijs1 = is*3+js +1

    js = js-1
    ijs2 = is*2+js +1

    fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
           - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,-4,0))

    fderiv(OPS_ACC2(0,0,0)) = acf5xy*fstora(OPS_ACC_MD3(ijs1,0,0,0)) + bcf5xy*fstorb(OPS_ACC_MD4(ijs1,0,0,0))  &
                            + ccf5xy*fstorc(OPS_ACC_MD5(ijs2,0,0,0)) + dcf5xy*fdiffd

END SUBROUTINE d2fdxy_kernel_lh_ydir_corner_eqV

!a3d_p440_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p330_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0))
    fdiffc = functn(OPS_ACC1(3,3,0))  - functn(OPS_ACC1(-3,3,0))  &
           - functn(OPS_ACC1(3,-3,0)) + functn(OPS_ACC1(-3,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc

END SUBROUTINE d2fdxy_kernel_rh_ydir_6th_centered

!a3d_p220_m220_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0))  - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_rh_ydir_4th_centered

!a3d_p210_m230_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p200_m240_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p400_p040_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqA(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p300_m100_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqB(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p410_p030_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqC(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p300_m140_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqD(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p210_m230_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqE(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p200_m240_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqF(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p420_m020_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqG(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p320_m120_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqH(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqI(functn, fderiv, fstora, fstorb, nyglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv, fstora, fstorb
    real(kind=8) :: fdiffa, fdiffb

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nyglblm4

    integer :: is, js, ijs
    is = idx(1) -3
    js = idx(2) -nyglblm4
    ijs = is*3+js +1

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

    fstora(OPS_ACC_MD3(ijs,0,0,0)) = fdiffa
    fstorb(OPS_ACC_MD4(ijs,0,0,0)) = fdiffb

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqI

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqJ(functn, fderiv, fstora, fstorb, fstorc, nyglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb
    real(kind=8), dimension(1) :: fderiv, fstorc
    real(kind=8) :: fdiffc

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nyglblm4

    integer :: is, js, ijs
    is = idx(1) -3
    js = idx(2) -nyglblm4
    ijs = is*3+js +1

    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fstora(OPS_ACC_MD3(ijs,0,0,0)) + bcf4xy*fstorb(OPS_ACC_MD4(ijs,0,0,0))  &
              + ccf4xy*fdiffc

    is = is-1
    ijs = is*2+js +1
    fstorc(OPS_ACC_MD5(ijs,0,0,0)) = fdiffc

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqJ

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqK(functn, fderiv, fstora, fstorb, fstorc, nyglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb, fstorc
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffd

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nyglblm4

    integer :: is, js, ijs1, ijs2
    is = idx(1) -3
    js = idx(2) -nyglblm4
    ijs1 = is*3+js +1

    is = is-1
    ijs2 = is*2+js +1

    fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
           - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,-4,0))

    fderiv(OPS_ACC2(0,0,0)) = acf5xy*fstora(OPS_ACC_MD3(ijs1,0,0,0)) + bcf5xy*fstorb(OPS_ACC_MD4(ijs1,0,0,0))  &
                            + ccf5xy*fstorc(OPS_ACC_MD5(ijs2,0,0,0)) + dcf5xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqK

!a3d_p000_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqL(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p110_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqM(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p010_m430_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqN(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p100_m340_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqO(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p210_m230_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqP(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p200_m240_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqQ(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p020_m420_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqR(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

!a3d_p120_m320_mixed_xy
SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqS(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

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

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqT(functn, fderiv, fstora, fstorb, nxglblm4, nyglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv, fstora, fstorb
    real(kind=8) :: fdiffa, fdiffb

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nxglblm4, nyglblm4

    integer :: is, js, ijs
    is = idx(1) -nxglblm4
    js = idx(2) -nyglblm4
    ijs = is*3+js +1

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
           - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
           - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

    fstora(OPS_ACC_MD3(ijs,0,0,0)) = fdiffa
    fstorb(OPS_ACC_MD4(ijs,0,0,0)) = fdiffb

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqT

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqU(functn, fderiv, fstora, fstorb, fstorc, nxglblm4, nyglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb
    real(kind=8), dimension(1) :: fderiv, fstorc
    real(kind=8) :: fdiffc

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nxglblm4, nyglblm4

    integer :: is, js, ijs
    is = idx(1) -nxglblm4
    js = idx(2) -nyglblm4
    ijs = is*3+js +1

    fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
           - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acf4xy*fstora(OPS_ACC_MD3(ijs,0,0,0)) + bcf4xy*fstorb(OPS_ACC_MD4(ijs,0,0,0))  &
              + ccf4xy*fdiffc

    ijs = is*2+js +1
    fstorc(OPS_ACC_MD5(ijs,0,0,0)) = fdiffc

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqU

SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqV(functn, fderiv, fstora, fstorb, fstorc, nxglblm4, nyglblm4, idx)
     use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb, fstorc
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffd

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nxglblm4, nyglblm4

    integer :: is, js, ijs1, ijs2
    is = idx(1) -nxglblm4
    js = idx(2) -nyglblm4

    ijs1 = is*3+js +1
    ijs2 = is*2+js +1

    fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
           - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,-4,0))

    fderiv(OPS_ACC2(0,0,0)) = acf5xy*fstora(OPS_ACC_MD3(ijs1,0,0,0)) + bcf5xy*fstorb(OPS_ACC_MD4(ijs1,0,0,0))  &
                            + ccf5xy*fstorc(OPS_ACC_MD5(ijs2,0,0,0)) + dcf5xy*fdiffd

END SUBROUTINE d2fdxy_kernel_rh_ydir_corner_eqV

SUBROUTINE d2fdxy_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdelx*ovdely

END SUBROUTINE d2fdxy_kernel_scaling
