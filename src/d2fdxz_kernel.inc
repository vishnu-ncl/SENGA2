SUBROUTINE d2fdxz_kernel_null(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdxz_kernel_null

!a3d_p505_m505_mixed_xz
SUBROUTINE d2fdxz_kernel_interior(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd,fdiffe

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,-3))  &
           - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,-3))
    fdiffd = functn(OPS_ACC1(4,0,4)) - functn(OPS_ACC1(4,0,-4))  &
           - functn(OPS_ACC1(-4,0,4)) + functn(OPS_ACC1(-4,0,-4))
    fdiffe = functn(OPS_ACC1(5,0,5)) - functn(OPS_ACC1(5,0,-5))  &
           - functn(OPS_ACC1(-5,0,5)) + functn(OPS_ACC1(-5,0,-5))

    fderiv(OPS_ACC2(0,0,0)) = acofxz*fdiffa + bcofxz*fdiffb  &
          + ccofxz*fdiffc + dcofxz*fdiffd  &
          + ecofxz*fdiffe

END SUBROUTINE d2fdxz_kernel_interior

!a3d_p402_m002_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_xdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(1,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(2,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(3,0,1)) - functn(OPS_ACC1(3,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(3,0,2)) - functn(OPS_ACC1(3,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(4,0,1)) - functn(OPS_ACC1(4,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(4,0,2)) - functn(OPS_ACC1(4,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_xdir_4th_onesided

!a3d_p302_m102_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_xdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(-1,0,1))  - functn(OPS_ACC1(-1,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(-1,0,2))  - functn(OPS_ACC1(-1,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(1,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(2,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(3,0,1)) - functn(OPS_ACC1(3,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(3,0,2)) - functn(OPS_ACC1(3,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
          + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_xdir_4th_mixed

!a3d_p202_m202_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_xdir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(-2,0,2))  + functn(OPS_ACC1(-2,0,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

END SUBROUTINE d2fdxz_kernel_lh_xdir_4th_centered

!a3d_p303_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_xdir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,-3))  &
           - functn(OPS_ACC1(-3,0,3))  + functn(OPS_ACC1(-3,0,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4xz*fdiffa + bcf4xz*fdiffb  &
          + ccf4xz*fdiffc

END SUBROUTINE d2fdxz_kernel_lh_xdir_6th_centered

!a3d_p404_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_xdir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,-3))  &
           - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,-3))
    fdiffd = functn(OPS_ACC1(4,0,4)) - functn(OPS_ACC1(4,0,-4))  &
           - functn(OPS_ACC1(-4,0,4))  + functn(OPS_ACC1(-4,0,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5xz*fdiffa + bcf5xz*fdiffb  &
          + ccf5xz*fdiffc + dcf5xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_xdir_8th_centered

!a3d_p404_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_xdir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,-3))  &
           - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,-3))
    fdiffd = functn(OPS_ACC1(4,0,4))  - functn(OPS_ACC1(4,0,-4))  &
           - functn(OPS_ACC1(-4,0,4)) + functn(OPS_ACC1(-4,0,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5xz*fdiffa + bcf5xz*fdiffb  &
          + ccf5xz*fdiffc + dcf5xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_xdir_8th_centered

!a3d_p303_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_xdir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,3))  - functn(OPS_ACC1(3,0,-3))  &
           - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4xz*fdiffa + bcf4xz*fdiffb  &
          + ccf4xz*fdiffc

END SUBROUTINE d2fdxz_kernel_rh_xdir_6th_centered

!a3d_p202_m202_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_xdir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2))  - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

END SUBROUTINE d2fdxz_kernel_rh_xdir_4th_centered

!a3d_p102_m302_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_xdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(1,0,1))  + functn(OPS_ACC1(1,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(1,0,2))  + functn(OPS_ACC1(1,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-1,0,2)) + functn(OPS_ACC1(-1,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-2,0,1)) + functn(OPS_ACC1(-2,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-3,0,1)) + functn(OPS_ACC1(-3,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-3,0,2)) + functn(OPS_ACC1(-3,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
          + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_xdir_4th_mixed

!a3d_p002_m402_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_xdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-1,0,2)) + functn(OPS_ACC1(-1,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-2,0,1)) + functn(OPS_ACC1(-2,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-3,0,1)) + functn(OPS_ACC1(-3,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-3,0,2)) + functn(OPS_ACC1(-3,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-4,0,1)) + functn(OPS_ACC1(-4,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-4,0,2)) + functn(OPS_ACC1(-4,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_xdir_4th_onesided

!a3d_p204_m200_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(-2,0,1))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(-1,0,2))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,3)) - functn(OPS_ACC1(-1,0,3))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,3)) - functn(OPS_ACC1(-2,0,3))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,4)) - functn(OPS_ACC1(-1,0,4))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,4)) - functn(OPS_ACC1(-2,0,4))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_4th_onesided

!a3d_p203_m201_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,-1))  - functn(OPS_ACC1(-1,0,-1))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,-1))  - functn(OPS_ACC1(-2,0,-1))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(-2,0,1))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(-1,0,2))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,3)) - functn(OPS_ACC1(-1,0,3))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,3)) - functn(OPS_ACC1(-2,0,3))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
          + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_4th_mixed

!a3d_p202_m202_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
           - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
           - functn(OPS_ACC1(2,0,-2))  + functn(OPS_ACC1(-2,0,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

END SUBROUTINE d2fdxz_kernel_lh_zdir_4th_centered

!a3d_p303_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
           - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
           - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(-3,0,3))  &
           - functn(OPS_ACC1(3,0,-3))  + functn(OPS_ACC1(-3,0,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4xz*fdiffa + bcf4xz*fdiffb  &
          + ccf4xz*fdiffc

END SUBROUTINE d2fdxz_kernel_lh_zdir_6th_centered

!a3d_p404_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
           - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
           - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(-3,0,3))  &
           - functn(OPS_ACC1(3,0,-3)) + functn(OPS_ACC1(-3,0,-3))
    fdiffd = functn(OPS_ACC1(4,0,4)) - functn(OPS_ACC1(-4,0,4))  &
           - functn(OPS_ACC1(4,0,-4))  + functn(OPS_ACC1(-4,0,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5xz*fdiffa + bcf5xz*fdiffb  &
          + ccf5xz*fdiffc + dcf5xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_8th_centered

!a3d_p404_p000_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqA(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,0))  &
           - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,0))  &
           - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,0))  &
           - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(4,0,4)) - functn(OPS_ACC1(4,0,0))  &
           - functn(OPS_ACC1(0,0,4))  + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1xz*fdiffa + bcc1xz*fdiffb  &
          + ccc1xz*fdiffc + dcc1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqA

!a3d_p303_m101_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqB(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(-1,0,-1))   - functn(OPS_ACC1(-1,0,0))  &
           - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,0))  &
           - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,0))  &
           - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,0))  &
           - functn(OPS_ACC1(0,0,3)) + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc2xz*fdiffa + bcc2xz*fdiffb  &
          + ccc2xz*fdiffc + dcc2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqB

!a3d_p403_m001_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqC(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xz*(functn(OPS_ACC1(1,0,-1))  - functn(OPS_ACC1(1,0,0))  &
                   - functn(OPS_ACC1(0,0,-1))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,0))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(1,0,0))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(1,0,3)) - functn(OPS_ACC1(1,0,0))  &
                   - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2xz*(functn(OPS_ACC1(2,0,-1))  - functn(OPS_ACC1(2,0,0))  &
                     - functn(OPS_ACC1(0,0,-1))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(2,0,0))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,0))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(2,0,3)) - functn(OPS_ACC1(2,0,0))  &
                   - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0)))
    fdiffc = acf2xz*(functn(OPS_ACC1(3,0,-1))  - functn(OPS_ACC1(3,0,0))  &
                     - functn(OPS_ACC1(0,0,-1))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(3,0,1)) - functn(OPS_ACC1(3,0,0))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(3,0,2)) - functn(OPS_ACC1(3,0,0))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,0))  &
                   - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2xz*(functn(OPS_ACC1(4,0,-1))  - functn(OPS_ACC1(4,0,0))  &
                   - functn(OPS_ACC1(0,0,-1))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(4,0,1)) - functn(OPS_ACC1(4,0,0))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(4,0,2)) - functn(OPS_ACC1(4,0,0))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(4,0,3)) - functn(OPS_ACC1(4,0,0))  &
                   - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqC

!a3d_p304_m100_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqD(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xz*(functn(OPS_ACC1(-1,0,1))  - functn(OPS_ACC1(0,0,1))  &
                   - functn(OPS_ACC1(-1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(0,0,1))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(0,0,1))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,1)) - functn(OPS_ACC1(0,0,1))  &
                   - functn(OPS_ACC1(3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2xz*(functn(OPS_ACC1(-1,0,2))  - functn(OPS_ACC1(0,0,2))  &
                   - functn(OPS_ACC1(-1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(0,0,2))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(0,0,2))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,2)) - functn(OPS_ACC1(0,0,2))  &
                   - functn(OPS_ACC1(3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffc = acf2xz*(functn(OPS_ACC1(-1,0,3))  - functn(OPS_ACC1(0,0,3))  &
                   - functn(OPS_ACC1(-1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,3)) - functn(OPS_ACC1(0,0,3))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,3)) - functn(OPS_ACC1(0,0,3))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(0,0,3))  &
                   - functn(OPS_ACC1(3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2xz*(functn(OPS_ACC1(-1,0,4))  - functn(OPS_ACC1(0,0,4))  &
                   - functn(OPS_ACC1(-1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,4)) - functn(OPS_ACC1(0,0,4))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,4)) - functn(OPS_ACC1(0,0,4))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,4)) - functn(OPS_ACC1(0,0,4))  &
                   - functn(OPS_ACC1(3,0,0))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqD

!a3d_p204_m200_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqE(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(-2,0,1))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(-1,0,2))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,3)) - functn(OPS_ACC1(-1,0,3))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,3)) - functn(OPS_ACC1(-2,0,3))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,4)) - functn(OPS_ACC1(-1,0,4))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,4)) - functn(OPS_ACC1(-2,0,4))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqE

!a3d_p203_m201_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqF(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,-1))  - functn(OPS_ACC1(-1,0,-1))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,-1))  - functn(OPS_ACC1(-2,0,-1))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(-2,0,1))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(-1,0,2))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,3)) - functn(OPS_ACC1(-1,0,3))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,3)) - functn(OPS_ACC1(-2,0,3))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqF

!a3d_p402_m002_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqG(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(1,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(2,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(3,0,1)) - functn(OPS_ACC1(3,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(3,0,2)) - functn(OPS_ACC1(3,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(4,0,1)) - functn(OPS_ACC1(4,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(4,0,2)) - functn(OPS_ACC1(4,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqG

!a3d_p302_m102_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqH(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(-1,0,1))  - functn(OPS_ACC1(-1,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(-1,0,2))  - functn(OPS_ACC1(-1,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(1,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(2,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(3,0,1)) - functn(OPS_ACC1(3,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(3,0,2)) - functn(OPS_ACC1(3,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqH

!a3d_p004_m400_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqL(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,0))
    fdiffb = functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,0))
    fdiffc = functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,0))
    fdiffd = functn(OPS_ACC1(0,0,4))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-4,0,4)) + functn(OPS_ACC1(-4,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1xz*fdiffa + bcc1xz*fdiffb  &
          + ccc1xz*fdiffc + dcc1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqL

!a3d_p100_m300_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqM(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,0,-1))  - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(1,0,-1))   + functn(OPS_ACC1(1,0,0))
    fdiffb = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,0))
    fdiffc = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,0))
    fdiffd = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))  &
           - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc2xz*fdiffa + bcc2xz*fdiffb  &
          + ccc2xz*fdiffc + dcc2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqM

!a3d_p003_m401_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqN(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xz*(functn(OPS_ACC1(0,0,-1))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,0,-1))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,0,2)) + functn(OPS_ACC1(-1,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,0,3)) + functn(OPS_ACC1(-1,0,0)))
    fdiffb = acf2xz*(functn(OPS_ACC1(0,0,-1))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-2,0,-1))  + functn(OPS_ACC1(-2,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-2,0,1)) + functn(OPS_ACC1(-2,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-2,0,3)) + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acf2xz*(functn(OPS_ACC1(0,0,-1))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-3,0,-1))  + functn(OPS_ACC1(-3,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-3,0,1)) + functn(OPS_ACC1(-3,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-3,0,2)) + functn(OPS_ACC1(-3,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,0)))
    fdiffd = acf2xz*(functn(OPS_ACC1(0,0,-1))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-4,0,-1))  + functn(OPS_ACC1(-4,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-4,0,1)) + functn(OPS_ACC1(-4,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-4,0,2)) + functn(OPS_ACC1(-4,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-4,0,3)) + functn(OPS_ACC1(-4,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqN

!a3d_p104_m300_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqO(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xz*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(1,0,1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(1,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(-1,0,1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(-2,0,1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(-3,0,1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-3,0,0)))
    fdiffb = acf2xz*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(1,0,2))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(1,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(-1,0,2))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(-2,0,2))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(-3,0,2))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-3,0,0)))
    fdiffc = acf2xz*(functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(1,0,3))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(1,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(-1,0,3))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(-2,0,3))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(-3,0,3))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-3,0,0)))
    fdiffd = acf2xz*(functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(1,0,4))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(1,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(-1,0,4))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(-2,0,4))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-2,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(-3,0,4))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(-3,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqO

!a3d_p204_m200_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqP(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(-2,0,1))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(-1,0,2))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,3)) - functn(OPS_ACC1(-1,0,3))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,3)) - functn(OPS_ACC1(-2,0,3))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,4)) - functn(OPS_ACC1(-1,0,4))  &
                   - functn(OPS_ACC1(1,0,0))  + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,4)) - functn(OPS_ACC1(-2,0,4))  &
                   - functn(OPS_ACC1(2,0,0))  + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqP

!a3d_p203_m201_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqQ(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,-1))  - functn(OPS_ACC1(-1,0,-1))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,-1))  - functn(OPS_ACC1(-2,0,-1))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(-2,0,1))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(-1,0,2))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,3)) - functn(OPS_ACC1(-1,0,3))  &
                   - functn(OPS_ACC1(1,0,0)) + functn(OPS_ACC1(-1,0,0)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,3)) - functn(OPS_ACC1(-2,0,3))  &
                   - functn(OPS_ACC1(2,0,0)) + functn(OPS_ACC1(-2,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqQ

!a3d_p002_m402_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqR(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-1,0,2)) + functn(OPS_ACC1(-1,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-2,0,1)) + functn(OPS_ACC1(-2,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-3,0,1)) + functn(OPS_ACC1(-3,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-3,0,2)) + functn(OPS_ACC1(-3,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-4,0,1)) + functn(OPS_ACC1(-4,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-4,0,2)) + functn(OPS_ACC1(-4,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqR

!a3d_p102_m302_mixed_xz
SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqS(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(1,0,1))  + functn(OPS_ACC1(1,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(1,0,2))  + functn(OPS_ACC1(1,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-1,0,2)) + functn(OPS_ACC1(-1,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-2,0,1)) + functn(OPS_ACC1(-2,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-3,0,1)) + functn(OPS_ACC1(-3,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-3,0,2)) + functn(OPS_ACC1(-3,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_lh_zdir_corner_eqS

!a3d_p404_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_8th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
           - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
           - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(-3,0,3))  &
           - functn(OPS_ACC1(3,0,-3)) + functn(OPS_ACC1(-3,0,-3))
    fdiffd = functn(OPS_ACC1(4,0,4))  - functn(OPS_ACC1(-4,0,4))  &
           - functn(OPS_ACC1(4,0,-4)) + functn(OPS_ACC1(-4,0,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5xz*fdiffa + bcf5xz*fdiffb  &
          + ccf5xz*fdiffc + dcf5xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_8th_centered

!a3d_p303_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_6th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
           - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
           - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,3))  - functn(OPS_ACC1(-3,0,3))  &
           - functn(OPS_ACC1(3,0,-3)) + functn(OPS_ACC1(-3,0,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4xz*fdiffa + bcf4xz*fdiffb  &
          + ccf4xz*fdiffc

END SUBROUTINE d2fdxz_kernel_rh_zdir_6th_centered

!a3d_p202_m202_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_4th_centered(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
           - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2))  - functn(OPS_ACC1(-2,0,2))  &
           - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

END SUBROUTINE d2fdxz_kernel_rh_zdir_4th_centered

!a3d_p201_m203_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,1))  + functn(OPS_ACC1(-1,0,1)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,1))  + functn(OPS_ACC1(-2,0,1)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-1)) + functn(OPS_ACC1(-2,0,-1)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-2)) + functn(OPS_ACC1(-1,0,-2)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-3)) + functn(OPS_ACC1(-1,0,-3)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-3)) + functn(OPS_ACC1(-2,0,-3)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
          + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_4th_mixed

!a3d_p200_m204_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-1)) + functn(OPS_ACC1(-2,0,-1)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-2)) + functn(OPS_ACC1(-1,0,-2)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-3)) + functn(OPS_ACC1(-1,0,-3)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-3)) + functn(OPS_ACC1(-2,0,-3)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-4)) + functn(OPS_ACC1(-1,0,-4)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-4)) + functn(OPS_ACC1(-2,0,-4)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_4th_onesided

!a3d_p400_p004_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqA(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2))
    fdiffc = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,0,-3))  &
           - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3))
    fdiffd = functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,0,-4))  &
           - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-4))

    fderiv(OPS_ACC2(0,0,0)) = acc1xz*fdiffa + bcc1xz*fdiffb  &
          + ccc1xz*fdiffc + dcc1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqA

!a3d_p300_m100_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqB(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(-1,0,0))  - functn(OPS_ACC1(-1,0,1))  &
           - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,0,1))
    fdiffb = functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,0,-1))
    fdiffc = functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,0,-2))
    fdiffd = functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,0,-3))  &
           - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,0,-3))

    fderiv(OPS_ACC2(0,0,0)) = acc2xz*fdiffa + bcc2xz*fdiffb  &
          + ccc2xz*fdiffc + dcc2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqB

!a3d_p401_p003_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqC(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xz*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,0,1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,1)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,0,-1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1)))  &
           + ccf2xz*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,0,-2))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2)))  &
           + dcf2xz*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(1,0,-3))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3)))
    fdiffb = acf2xz*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,0,1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,1)))  &
           + bcf2xz*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,0,-1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,0,-2))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2)))  &
           + dcf2xz*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(2,0,-3))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3)))
    fdiffc = acf2xz*(functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,0,1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,1)))  &
           + bcf2xz*(functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,0,-1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1)))  &
           + ccf2xz*(functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,0,-2))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,0)) - functn(OPS_ACC1(3,0,-3))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3)))
    fdiffd = acf2xz*(functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,0,1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,1)))  &
           + bcf2xz*(functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,0,-1))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1)))  &
           + ccf2xz*(functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,0,-2))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2)))  &
           + dcf2xz*(functn(OPS_ACC1(4,0,0)) - functn(OPS_ACC1(4,0,-3))  &
                   - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqC

!a3d_p300_m104_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqD(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xz*(functn(OPS_ACC1(-1,0,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,0,-1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(0,0,-1)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(2,0,-1)) + functn(OPS_ACC1(0,0,-1)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(3,0,-1)) + functn(OPS_ACC1(0,0,-1)))
    fdiffb = acf2xz*(functn(OPS_ACC1(-1,0,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,0,-2))  + functn(OPS_ACC1(0,0,-2)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(1,0,-2)) + functn(OPS_ACC1(0,0,-2)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(0,0,-2)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(3,0,-2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acf2xz*(functn(OPS_ACC1(-1,0,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,0,-3))  + functn(OPS_ACC1(0,0,-3)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(1,0,-3)) + functn(OPS_ACC1(0,0,-3)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(2,0,-3)) + functn(OPS_ACC1(0,0,-3)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(3,0,-3)) + functn(OPS_ACC1(0,0,-3)))
    fdiffd = acf2xz*(functn(OPS_ACC1(-1,0,0))   - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(-1,0,-4))  + functn(OPS_ACC1(0,0,-4)))  &
           + bcf2xz*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(1,0,-4)) + functn(OPS_ACC1(0,0,-4)))  &
           + ccf2xz*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(2,0,-4)) + functn(OPS_ACC1(0,0,-4)))  &
           + dcf2xz*(functn(OPS_ACC1(3,0,0))  - functn(OPS_ACC1(0,0,0))  &
                   - functn(OPS_ACC1(3,0,-4)) + functn(OPS_ACC1(0,0,-4)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqD

!a3d_p201_m203_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqE(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,1))  + functn(OPS_ACC1(-1,0,1)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,1))  + functn(OPS_ACC1(-2,0,1)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-1)) + functn(OPS_ACC1(-2,0,-1)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-2)) + functn(OPS_ACC1(-1,0,-2)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-3)) + functn(OPS_ACC1(-1,0,-3)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-3)) + functn(OPS_ACC1(-2,0,-3)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqE

!a3d_p200_m204_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqF(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-1)) + functn(OPS_ACC1(-2,0,-1)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-2)) + functn(OPS_ACC1(-1,0,-2)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-3)) + functn(OPS_ACC1(-1,0,-3)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-3)) + functn(OPS_ACC1(-2,0,-3)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-4)) + functn(OPS_ACC1(-1,0,-4)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-4)) + functn(OPS_ACC1(-2,0,-4)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqF

!a3d_p402_m002_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqG(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(1,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(2,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(3,0,1)) - functn(OPS_ACC1(3,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(3,0,2)) - functn(OPS_ACC1(3,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(4,0,1)) - functn(OPS_ACC1(4,0,-1))  &
                   - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(4,0,2)) - functn(OPS_ACC1(4,0,-2))  &
                   - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqG

!a3d_p302_m102_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqH(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(-1,0,1))  - functn(OPS_ACC1(-1,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(-1,0,2))  - functn(OPS_ACC1(-1,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(1,0,2)) - functn(OPS_ACC1(1,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(2,0,1)) - functn(OPS_ACC1(2,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(3,0,1)) - functn(OPS_ACC1(3,0,-1))  &
                   - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(3,0,2)) - functn(OPS_ACC1(3,0,-2))  &
                   - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqH

!a3d_p000_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqL(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(-1,0,-1)) - functn(OPS_ACC1(-1,0,0))  &
           - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(-2,0,-2)) - functn(OPS_ACC1(-2,0,0))  &
           - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(-3,0,-3)) - functn(OPS_ACC1(-3,0,0))  &
           - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(-4,0,-4)) - functn(OPS_ACC1(-4,0,0))  &
           - functn(OPS_ACC1(0,0,-4))  + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1xz*fdiffa + bcc1xz*fdiffb  &
          + ccc1xz*fdiffc + dcc1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqL

!a3d_p101_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqM(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(1,0,1))   - functn(OPS_ACC1(1,0,0))  &
           - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(-1,0,-1)) - functn(OPS_ACC1(-1,0,0))  &
           - functn(OPS_ACC1(0,0,-1)) + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(-2,0,-2)) - functn(OPS_ACC1(-2,0,0))  &
           - functn(OPS_ACC1(0,0,-2)) + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(-3,0,-3)) - functn(OPS_ACC1(-3,0,0))  &
           - functn(OPS_ACC1(0,0,-3)) + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc2xz*fdiffa + bcc2xz*fdiffb  &
          + ccc2xz*fdiffc + dcc2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqM

!a3d_p001_m403_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqN(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xz*(functn(OPS_ACC1(-1,0,1))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(0,0,1))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(-1,0,-1)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(-1,0,-2)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(-1,0,-3)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2xz*(functn(OPS_ACC1(-2,0,1))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(0,0,1))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(-2,0,-1)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(-2,0,-2)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(-2,0,-3)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0)))
    fdiffc = acf2xz*(functn(OPS_ACC1(-3,0,1))  - functn(OPS_ACC1(-3,0,0))  &
                   - functn(OPS_ACC1(0,0,1))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(-3,0,-1)) - functn(OPS_ACC1(-3,0,0))  &
                   - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(-3,0,-2)) - functn(OPS_ACC1(-3,0,0))  &
                   - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(-3,0,-3)) - functn(OPS_ACC1(-3,0,0))  &
                   - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2xz*(functn(OPS_ACC1(-4,0,1))  - functn(OPS_ACC1(-4,0,0))  &
                   - functn(OPS_ACC1(0,0,1))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(-4,0,-1)) - functn(OPS_ACC1(-4,0,0))  &
                   - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(-4,0,-2)) - functn(OPS_ACC1(-4,0,0))  &
                   - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(-4,0,-3)) - functn(OPS_ACC1(-4,0,0))  &
                   - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqN

!a3d_p100_m304_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqO(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2xz*(functn(OPS_ACC1(1,0,-1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(-1,0,-1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(-2,0,-1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(-3,0,-1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2xz*(functn(OPS_ACC1(1,0,-2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(-1,0,-2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(-2,0,-2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(-3,0,-2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffc = acf2xz*(functn(OPS_ACC1(1,0,-3))  - functn(OPS_ACC1(0,0,-3))  &
                   - functn(OPS_ACC1(1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(-1,0,-3)) - functn(OPS_ACC1(0,0,-3))  &
                   - functn(OPS_ACC1(-1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(-2,0,-3)) - functn(OPS_ACC1(0,0,-3))  &
                   - functn(OPS_ACC1(-2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(-3,0,-3)) - functn(OPS_ACC1(0,0,-3))  &
                   - functn(OPS_ACC1(-3,0,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2xz*(functn(OPS_ACC1(1,0,-4))  - functn(OPS_ACC1(0,0,-4))  &
                   - functn(OPS_ACC1(1,0,0))   + functn(OPS_ACC1(0,0,0)))  &
           + bcf2xz*(functn(OPS_ACC1(-1,0,-4)) - functn(OPS_ACC1(0,0,-4))  &
                   - functn(OPS_ACC1(-1,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + ccf2xz*(functn(OPS_ACC1(-2,0,-4)) - functn(OPS_ACC1(0,0,-4))  &
                   - functn(OPS_ACC1(-2,0,0))  + functn(OPS_ACC1(0,0,0)))  &
           + dcf2xz*(functn(OPS_ACC1(-3,0,-4)) - functn(OPS_ACC1(0,0,-4))  &
                   - functn(OPS_ACC1(-3,0,0))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqO

!a3d_p201_m203_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqP(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,1))  + functn(OPS_ACC1(-1,0,1)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,1))  + functn(OPS_ACC1(-2,0,1)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-1)) + functn(OPS_ACC1(-2,0,-1)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-2)) + functn(OPS_ACC1(-1,0,-2)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0)) - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-3)) + functn(OPS_ACC1(-1,0,-3)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0)) - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-3)) + functn(OPS_ACC1(-2,0,-3)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqP

!a3d_p200_m204_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqQ(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-1)) + functn(OPS_ACC1(-2,0,-1)))
    fdiffb = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-2)) + functn(OPS_ACC1(-1,0,-2)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffc = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-3)) + functn(OPS_ACC1(-1,0,-3)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-3)) + functn(OPS_ACC1(-2,0,-3)))
    fdiffd = acofx1*(functn(OPS_ACC1(1,0,0))  - functn(OPS_ACC1(-1,0,0))  &
                   - functn(OPS_ACC1(1,0,-4)) + functn(OPS_ACC1(-1,0,-4)))  &
           + bcofx1*(functn(OPS_ACC1(2,0,0))  - functn(OPS_ACC1(-2,0,0))  &
                   - functn(OPS_ACC1(2,0,-4)) + functn(OPS_ACC1(-2,0,-4)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqQ

!a3d_p002_m402_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqR(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-1,0,2)) + functn(OPS_ACC1(-1,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-2,0,1)) + functn(OPS_ACC1(-2,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-3,0,1)) + functn(OPS_ACC1(-3,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-3,0,2)) + functn(OPS_ACC1(-3,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-4,0,1)) + functn(OPS_ACC1(-4,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-4,0,2)) + functn(OPS_ACC1(-4,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqR

!a3d_p102_m302_mixed_xz
SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqS(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(1,0,1))  + functn(OPS_ACC1(1,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(1,0,2))  + functn(OPS_ACC1(1,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-1,0,2)) + functn(OPS_ACC1(-1,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-2,0,1)) + functn(OPS_ACC1(-2,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
                   - functn(OPS_ACC1(-3,0,1)) + functn(OPS_ACC1(-3,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
                   - functn(OPS_ACC1(-3,0,2)) + functn(OPS_ACC1(-3,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd

END SUBROUTINE d2fdxz_kernel_rh_zdir_corner_eqS

SUBROUTINE d2fdxz_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC2(0,0,0)) = fderiv(OPS_ACC2(0,0,0))*ovdelx*ovdelz

END SUBROUTINE d2fdxz_kernel_scaling

