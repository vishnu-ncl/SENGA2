SUBROUTINE d2fdyz_kernel_null(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdyz_kernel_null

!a3d_p055_m055_mixed_yz
SUBROUTINE d2fdyz_kernel_interior(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd,fdiffe

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))
    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
           - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))
    fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
           - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))
    fdiffe = functn(OPS_ACC1(0,5,5)) - functn(OPS_ACC1(0,5,-5))  &
           - functn(OPS_ACC1(0,-5,5)) + functn(OPS_ACC1(0,-5,-5))

    fderiv(OPS_ACC2(0,0,0)) = acofyz*fdiffa + bcofyz*fdiffb  &
                            + ccofyz*fdiffc + dcofyz*fdiffd  &
                            + ecofyz*fdiffe

END SUBROUTINE d2fdyz_kernel_interior

!a3d_p042_m002_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_ydir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,1,-2))  &
           - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,2,-1))  &
           - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,3,1)) - functn(OPS_ACC1(0,3,-1))  &
           - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,3,2)) - functn(OPS_ACC1(0,3,-2))  &
           - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,4,1)) - functn(OPS_ACC1(0,4,-1))  &
           - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,4,2)) - functn(OPS_ACC1(0,4,-2))  &
           - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
                                + ccf1yz*fdiffc + dcf1yz*fdiffd


END SUBROUTINE d2fdyz_kernel_lh_ydir_4th_onesided

!a3d_p032_m012_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_ydir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,-1,1))  - functn(OPS_ACC1(0,-1,-1))  &
           - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,-1,2))  - functn(OPS_ACC1(0,-1,-2))  &
           - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,1,-2))  &
           - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,2,-1))  &
           - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,3,1)) - functn(OPS_ACC1(0,3,-1))  &
           - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,3,2)) - functn(OPS_ACC1(0,3,-2))  &
           - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
                            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_ydir_4th_mixed

!a3d_p022_m022_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_ydir_4th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2))  + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

END SUBROUTINE d2fdyz_kernel_lh_ydir_4th_centred

!a3d_p033_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_ydir_6th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))
    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
           - functn(OPS_ACC1(0,-3,3))  + functn(OPS_ACC1(0,-3,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4yz*fdiffa + bcf4yz*fdiffb  &
          + ccf4yz*fdiffc

END SUBROUTINE d2fdyz_kernel_lh_ydir_6th_centred

!a3d_p044_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_ydir_8th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))
    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
           - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))
    fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
           - functn(OPS_ACC1(0,-4,4))  + functn(OPS_ACC1(0,-4,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5yz*fdiffa + bcf5yz*fdiffb  &
                            + ccf5yz*fdiffc + dcf5yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_ydir_8th_centred

!a3d_p044_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_ydir_8th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))
    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
           - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))
    fdiffd = functn(OPS_ACC1(0,4,4))  - functn(OPS_ACC1(0,4,-4))  &
           - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5yz*fdiffa + bcf5yz*fdiffb  &
                        + ccf5yz*fdiffc + dcf5yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_ydir_8th_centred

!a3d_p033_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_ydir_6th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))
    fdiffc = functn(OPS_ACC1(0,3,3))  - functn(OPS_ACC1(0,3,-3))  &
           - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4yz*fdiffa + bcf4yz*fdiffb  &
                            + ccf4yz*fdiffc

END SUBROUTINE d2fdyz_kernel_rh_ydir_6th_centred

!a3d_p022_m022_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_ydir_4th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2))  - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

END SUBROUTINE d2fdyz_kernel_rh_ydir_4th_centred

!a3d_p012_m032_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_ydir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
           - functn(OPS_ACC1(0,1,1))  + functn(OPS_ACC1(0,1,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
           - functn(OPS_ACC1(0,1,2))  + functn(OPS_ACC1(0,1,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
           - functn(OPS_ACC1(0,-1,2)) + functn(OPS_ACC1(0,-1,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
           - functn(OPS_ACC1(0,-2,1)) + functn(OPS_ACC1(0,-2,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
           - functn(OPS_ACC1(0,-3,1)) + functn(OPS_ACC1(0,-3,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
           - functn(OPS_ACC1(0,-3,2)) + functn(OPS_ACC1(0,-3,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
                            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_ydir_4th_mixed

!a3d_p002_m042_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_ydir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
           - functn(OPS_ACC1(0,-1,2)) + functn(OPS_ACC1(0,-1,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
           - functn(OPS_ACC1(0,-2,1)) + functn(OPS_ACC1(0,-2,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
           - functn(OPS_ACC1(0,-3,1)) + functn(OPS_ACC1(0,-3,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
           - functn(OPS_ACC1(0,-3,2)) + functn(OPS_ACC1(0,-3,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
           - functn(OPS_ACC1(0,-4,1)) + functn(OPS_ACC1(0,-4,-1)))  &
           + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
           - functn(OPS_ACC1(0,-4,2)) + functn(OPS_ACC1(0,-4,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
                            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_ydir_4th_onesided

!a3d_p024_m020_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
           - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,-2,1))  &
           - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,-1,2))  &
           - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
           - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,3)) - functn(OPS_ACC1(0,-1,3))  &
           - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,3)) - functn(OPS_ACC1(0,-2,3))  &
           - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,4)) - functn(OPS_ACC1(0,-1,4))  &
           - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,4)) - functn(OPS_ACC1(0,-2,4))  &
           - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
                            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_4th_onesided

!a3d_p023_m021_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,-1))  - functn(OPS_ACC1(0,-1,-1))  &
           - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,-1))  - functn(OPS_ACC1(0,-2,-1))  &
           - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
           - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,-2,1))  &
           - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,-1,2))  &
           - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
           - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,3)) - functn(OPS_ACC1(0,-1,3))  &
           - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
           + bcofy1*(functn(OPS_ACC1(0,2,3)) - functn(OPS_ACC1(0,-2,3))  &
           - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
                            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_4th_mixed

!a3d_p022_m022_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_4th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
           - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
           - functn(OPS_ACC1(0,2,-2))  + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

END SUBROUTINE d2fdyz_kernel_lh_zdir_4th_centred

!a3d_p033_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_6th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
           - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
           - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2))
    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,-3,3))  &
           - functn(OPS_ACC1(0,3,-3))  + functn(OPS_ACC1(0,-3,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4yz*fdiffa + bcf4yz*fdiffb  &
                            + ccf4yz*fdiffc

END SUBROUTINE d2fdyz_kernel_lh_zdir_6th_centred

!a3d_p044_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_8th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
           - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
           - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2))
    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,-3,3))  &
           - functn(OPS_ACC1(0,3,-3)) + functn(OPS_ACC1(0,-3,-3))
    fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,-4,4))  &
           - functn(OPS_ACC1(0,4,-4))  + functn(OPS_ACC1(0,-4,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5yz*fdiffa + bcf5yz*fdiffb  &
          + ccf5yz*fdiffc + dcf5yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_8th_centred

!a3d_p044_p000_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqA(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,0))  &
           - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,0))  &
           - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,0))  &
           - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,0))  &
           - functn(OPS_ACC1(0,0,4))  + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1yz*fdiffa + bcc1yz*fdiffb  &
                            + ccc1yz*fdiffc + dcc1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqA

!a3d_p033_m011_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqB(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,-1,-1)) - functn(OPS_ACC1(0,-1,0))  &
           - functn(OPS_ACC1(0,0,-1)) + functn(OPS_ACC1(0,0,0))
    fdiffb = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,0))  &
           - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,0))
    fdiffc = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,0))  &
           - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,0))
    fdiffd = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,0))  &
           - functn(OPS_ACC1(0,0,3)) + functn(OPS_ACC1(0,0,0))

    fderiv(OPS_ACC2(0,0,0)) = acc2yz*fdiffa + bcc2yz*fdiffb  &
                            + ccc2yz*fdiffc + dcc2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqB

!a3d_p043_m001_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqC(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2yz*(functn(OPS_ACC1(0,1,-1))  - functn(OPS_ACC1(0,1,0))  &
                    - functn(OPS_ACC1(0,0,-1))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,0))  &
          - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,1,0))  &
          - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,1,3)) - functn(OPS_ACC1(0,1,0))  &
          - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2yz*(functn(OPS_ACC1(0,2,-1))  - functn(OPS_ACC1(0,2,0))  &
          - functn(OPS_ACC1(0,0,-1))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,2,0))  &
          - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,0))  &
          - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,2,3)) - functn(OPS_ACC1(0,2,0))  &
          - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0)))
    fdiffc = acf2yz*(functn(OPS_ACC1(0,3,-1))  - functn(OPS_ACC1(0,3,0))  &
          - functn(OPS_ACC1(0,0,-1))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,3,1)) - functn(OPS_ACC1(0,3,0))  &
          - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,3,2)) - functn(OPS_ACC1(0,3,0))  &
          - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,0))  &
          - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2yz*(functn(OPS_ACC1(0,4,-1))  - functn(OPS_ACC1(0,4,0))  &
          - functn(OPS_ACC1(0,0,-1))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,4,1)) - functn(OPS_ACC1(0,4,0))  &
          - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,4,2)) - functn(OPS_ACC1(0,4,0))  &
          - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,4,3)) - functn(OPS_ACC1(0,4,0))  &
          - functn(OPS_ACC1(0,0,3))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqC

!a3d_p034_m010_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqD(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2yz*(functn(OPS_ACC1(0,-1,1))  - functn(OPS_ACC1(0,0,1))  &
          - functn(OPS_ACC1(0,-1,0))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,0,1))  &
          - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,0,1))  &
          - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,1)) - functn(OPS_ACC1(0,0,1))  &
          - functn(OPS_ACC1(0,3,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffb = acf2yz*(functn(OPS_ACC1(0,-1,2))  - functn(OPS_ACC1(0,0,2))  &
          - functn(OPS_ACC1(0,-1,0))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,0,2))  &
          - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,0,2))  &
          - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,2)) - functn(OPS_ACC1(0,0,2))  &
          - functn(OPS_ACC1(0,3,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffc = acf2yz*(functn(OPS_ACC1(0,-1,3))  - functn(OPS_ACC1(0,0,3))  &
          - functn(OPS_ACC1(0,-1,0))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,3)) - functn(OPS_ACC1(0,0,3))  &
          - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,3)) - functn(OPS_ACC1(0,0,3))  &
          - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,0,3))  &
          - functn(OPS_ACC1(0,3,0))  + functn(OPS_ACC1(0,0,0)))
    fdiffd = acf2yz*(functn(OPS_ACC1(0,-1,4))  - functn(OPS_ACC1(0,0,4))  &
          - functn(OPS_ACC1(0,-1,0))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,4)) - functn(OPS_ACC1(0,0,4))  &
          - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,4)) - functn(OPS_ACC1(0,0,4))  &
          - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,4)) - functn(OPS_ACC1(0,0,4))  &
          - functn(OPS_ACC1(0,3,0))  + functn(OPS_ACC1(0,0,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqD

!a3d_p024_m020_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqE(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
            - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,-2,1))  &
            - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,-1,2))  &
            - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
            - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,3)) - functn(OPS_ACC1(0,-1,3))  &
            - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,3)) - functn(OPS_ACC1(0,-2,3))  &
            - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,4)) - functn(OPS_ACC1(0,-1,4))  &
            - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,4)) - functn(OPS_ACC1(0,-2,4))  &
            - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqE

!a3d_p023_m021_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqF(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,-1))  - functn(OPS_ACC1(0,-1,-1))  &
            - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,-1))  - functn(OPS_ACC1(0,-2,-1))  &
            - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
            - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,-2,1))  &
            - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,-1,2))  &
            - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
            - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,3)) - functn(OPS_ACC1(0,-1,3))  &
            - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,3)) - functn(OPS_ACC1(0,-2,3))  &
            - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqF

!a3d_p042_m002_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqG(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
            - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,1,-2))  &
            - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,2,-1))  &
            - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
            - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,3,1)) - functn(OPS_ACC1(0,3,-1))  &
            - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,3,2)) - functn(OPS_ACC1(0,3,-2))  &
            - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,4,1)) - functn(OPS_ACC1(0,4,-1))  &
            - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,4,2)) - functn(OPS_ACC1(0,4,-2))  &
            - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqG

!a3d_p032_m012_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqH(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,-1,1))  - functn(OPS_ACC1(0,-1,-1))  &
            - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,-1,2))  - functn(OPS_ACC1(0,-1,-2))  &
            - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffb = acofz1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
            - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,1,-2))  &
            - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffc = acofz1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,2,-1))  &
            - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
            - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
    fdiffd = acofz1*(functn(OPS_ACC1(0,3,1)) - functn(OPS_ACC1(0,3,-1))  &
            - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,3,2)) - functn(OPS_ACC1(0,3,-2))  &
            - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))

    fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqH

SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqI(functn, fderiv, fstora, fstorb, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv, fstora, fstorb
    real(kind=8) :: fdiffa, fdiffb

    integer, dimension(3), intent(in) :: idx
    integer :: js, ks, jks
    js = idx(2) -3
    ks = idx(3) -3
    jks = js*3+ks +1

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

    fstora(OPS_ACC_MD3(jks,0,0,0)) = fdiffa
    fstorb(OPS_ACC_MD4(jks,0,0,0)) = fdiffb

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqI

SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqJ(functn, fderiv, fstora, fstorb, fstorc, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb
    real(kind=8), dimension(1) :: fderiv, fstorc
    real(kind=8) :: fdiffc

    integer, dimension(3), intent(in) :: idx
    integer :: js, ks, jks
    js = idx(2) -3
    ks = idx(3) -3
    jks = js*3+ks +1

    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
           - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4yz*fstora(OPS_ACC_MD3(jks,0,0,0)) + bcf4yz*fstorb(OPS_ACC_MD4(jks,0,0,0))  &
              + ccf4yz*fdiffc

    js = js-1
    ks = ks-1
    jks = js*2+ks +1
    fstorc(OPS_ACC_MD5(jks,0,0,0)) = fdiffc

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqJ

SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqK(functn, fderiv, fstora, fstorb, fstorc, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb, fstorc
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffd

    integer, dimension(3), intent(in) :: idx
    integer :: js, ks, jks1, jks2
    js = idx(2) -3
    ks = idx(3) -3
    jks1 = js*3+ks +1

    js = js-1
    ks = ks-1
    jks2 = js*2+ks +1

    fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
           - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5yz*fstora(OPS_ACC_MD3(jks1,0,0,0)) + bcf5yz*fstorb(OPS_ACC_MD4(jks1,0,0,0))  &
                            + ccf5yz*fstorc(OPS_ACC_MD5(jks2,0,0,0)) + dcf5yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqK

!a3d_p004_m040_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqL(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,0))
    fdiffb = functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,0))
    fdiffc = functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,0))
    fdiffd = functn(OPS_ACC1(0,0,4))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,0))

    fderiv(OPS_ACC2(0,0,0)) = acc1yz*fdiffa + bcc1yz*fdiffb  &
          + ccc1yz*fdiffc + dcc1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqL

!a3d_p010_m030_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqM(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,0,-1))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,1,-1))   + functn(OPS_ACC1(0,1,0))
    fdiffb = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,0))
    fdiffc = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,0))
    fdiffd = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,0))

    fderiv(OPS_ACC2(0,0,0)) = acc2yz*fdiffa + bcc2yz*fdiffb  &
          + ccc2yz*fdiffc + dcc2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqM

!a3d_p003_m041_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqN(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2yz*(functn(OPS_ACC1(0,0,-1))   - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,-1))  + functn(OPS_ACC1(0,-1,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,2)) + functn(OPS_ACC1(0,-1,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,3)) + functn(OPS_ACC1(0,-1,0)))
      fdiffb = acf2yz*(functn(OPS_ACC1(0,0,-1))   - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-2,-1))  + functn(OPS_ACC1(0,-2,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-2,1)) + functn(OPS_ACC1(0,-2,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-2,3)) + functn(OPS_ACC1(0,-2,0)))
      fdiffc = acf2yz*(functn(OPS_ACC1(0,0,-1))   - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-3,-1))  + functn(OPS_ACC1(0,-3,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-3,1)) + functn(OPS_ACC1(0,-3,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-3,2)) + functn(OPS_ACC1(0,-3,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,0)))
      fdiffd = acf2yz*(functn(OPS_ACC1(0,0,-1))   - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-4,-1))  + functn(OPS_ACC1(0,-4,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-4,1)) + functn(OPS_ACC1(0,-4,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-4,2)) + functn(OPS_ACC1(0,-4,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,0,3))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-4,3)) + functn(OPS_ACC1(0,-4,0)))

      fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqN

!a3d_p014_m030_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqO(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2yz*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,1,1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,1,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,-1,1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-1,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,-2,1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-2,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,-3,1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-3,0)))
      fdiffb = acf2yz*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,1,2))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,1,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,-1,2))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-1,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,-2,2))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-2,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,-3,2))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-3,0)))
      fdiffc = acf2yz*(functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,1,3))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,1,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,-1,3))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-1,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,-2,3))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-2,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,-3,3))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-3,0)))
      fdiffd = acf2yz*(functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,1,4))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,1,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,-1,4))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-1,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,-2,4))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-2,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,-3,4))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,-3,0)))

      fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqO

!a3d_p024_m020_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqP(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
            - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,-2,1))  &
            - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
        fdiffb = acofy1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,-1,2))  &
            - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
            - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
        fdiffc = acofy1*(functn(OPS_ACC1(0,1,3)) - functn(OPS_ACC1(0,-1,3))  &
            - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,3)) - functn(OPS_ACC1(0,-2,3))  &
            - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))
        fdiffd = acofy1*(functn(OPS_ACC1(0,1,4)) - functn(OPS_ACC1(0,-1,4))  &
            - functn(OPS_ACC1(0,1,0))  + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,4)) - functn(OPS_ACC1(0,-2,4))  &
            - functn(OPS_ACC1(0,2,0))  + functn(OPS_ACC1(0,-2,0)))

        fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqP

!a3d_p023_m021_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqQ(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,-1))  - functn(OPS_ACC1(0,-1,-1))  &
            - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,-1))  - functn(OPS_ACC1(0,-2,-1))  &
            - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
        fdiffb = acofy1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
            - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,-2,1))  &
            - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
        fdiffc = acofy1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,-1,2))  &
            - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
            - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))
        fdiffd = acofy1*(functn(OPS_ACC1(0,1,3)) - functn(OPS_ACC1(0,-1,3))  &
            - functn(OPS_ACC1(0,1,0)) + functn(OPS_ACC1(0,-1,0)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,3)) - functn(OPS_ACC1(0,-2,3))  &
            - functn(OPS_ACC1(0,2,0)) + functn(OPS_ACC1(0,-2,0)))

        fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqQ

!a3d_p002_m042_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqR(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-1,2)) + functn(OPS_ACC1(0,-1,-2)))
        fdiffb = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-2,1)) + functn(OPS_ACC1(0,-2,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2)))
        fdiffc = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-3,1)) + functn(OPS_ACC1(0,-3,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-3,2)) + functn(OPS_ACC1(0,-3,-2)))
        fdiffd = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-4,1)) + functn(OPS_ACC1(0,-4,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-4,2)) + functn(OPS_ACC1(0,-4,-2)))

        fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqR

!a3d_p012_m032_mixed_yz
SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqS(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,1,1))  + functn(OPS_ACC1(0,1,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,1,2))  + functn(OPS_ACC1(0,1,-2)))
        fdiffb = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-1,2)) + functn(OPS_ACC1(0,-1,-2)))
        fdiffc = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-2,1)) + functn(OPS_ACC1(0,-2,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2)))
        fdiffd = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-3,1)) + functn(OPS_ACC1(0,-3,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-3,2)) + functn(OPS_ACC1(0,-3,-2)))

        fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqS

SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqT(functn, fderiv, fstora, fstorb, nyglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv, fstora, fstorb
    real(kind=8) :: fdiffa, fdiffb

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nyglblm4

    integer :: js, ks, jks
    js = idx(2) -nyglblm4
    ks = idx(3) -3
    jks = js*3+ks +1

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

    fstora(OPS_ACC_MD3(jks,0,0,0)) = fdiffa
    fstorb(OPS_ACC_MD4(jks,0,0,0)) = fdiffb

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqT

SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqU(functn, fderiv, fstora, fstorb, fstorc, nyglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb
    real(kind=8), dimension(1) :: fderiv, fstorc
    real(kind=8) :: fdiffc

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nyglblm4

    integer :: js, ks, jks
    js = idx(2) -nyglblm4
    ks = idx(3) -3
    jks = js*3+ks +1

    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
           - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4yz*fstora(OPS_ACC_MD3(jks,0,0,0)) + bcf4yz*fstorb(OPS_ACC_MD4(jks,0,0,0))  &
              + ccf4yz*fdiffc

    ks = ks-1
    jks = js*2+ks +1
    fstorc(OPS_ACC_MD5(jks,0,0,0)) = fdiffc

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqU

SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqV(functn, fderiv, fstora, fstorb, fstorc, nyglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb, fstorc
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffd

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nyglblm4

    integer :: js, ks, jks1, jks2
    js = idx(2) -nyglblm4
    ks = idx(3) -3
    jks1 = js*3+ks +1

    ks = ks-1
    jks2 = js*2+ks +1

    fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
           - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5yz*fstora(OPS_ACC_MD3(jks1,0,0,0)) + bcf5yz*fstorb(OPS_ACC_MD4(jks1,0,0,0))  &
                            + ccf5yz*fstorc(OPS_ACC_MD5(jks2,0,0,0)) + dcf5yz*fdiffd

END SUBROUTINE d2fdyz_kernel_lh_zdir_corner_eqV

!a3d_p044_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_8th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
          - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1))
      fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
          - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2))
      fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,-3,3))  &
          - functn(OPS_ACC1(0,3,-3)) + functn(OPS_ACC1(0,-3,-3))
      fdiffd = functn(OPS_ACC1(0,4,4))  - functn(OPS_ACC1(0,-4,4))  &
          - functn(OPS_ACC1(0,4,-4)) + functn(OPS_ACC1(0,-4,-4))

      fderiv(OPS_ACC2(0,0,0)) = acf5yz*fdiffa + bcf5yz*fdiffb  &
          + ccf5yz*fdiffc + dcf5yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_8th_centred

!a3d_p033_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_6th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
          - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1))
      fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
          - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2))
      fdiffc = functn(OPS_ACC1(0,3,3))  - functn(OPS_ACC1(0,-3,3))  &
          - functn(OPS_ACC1(0,3,-3)) + functn(OPS_ACC1(0,-3,-3))

      fderiv(OPS_ACC2(0,0,0)) = acf4yz*fdiffa + bcf4yz*fdiffb  &
          + ccf4yz*fdiffc

END SUBROUTINE d2fdyz_kernel_rh_zdir_6th_centred

!a3d_p022_m022_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_4th_centred(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
          - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1))
      fdiffb = functn(OPS_ACC1(0,2,2))  - functn(OPS_ACC1(0,-2,2))  &
          - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2))

      fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

END SUBROUTINE d2fdyz_kernel_rh_zdir_4th_centred

!a3d_p021_m023_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_4th_mixed(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,1,1))  + functn(OPS_ACC1(0,-1,1)))  &
          + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,2,1))  + functn(OPS_ACC1(0,-2,1)))
      fdiffb = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1)))  &
          + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,2,-1)) + functn(OPS_ACC1(0,-2,-1)))
      fdiffc = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,1,-2)) + functn(OPS_ACC1(0,-1,-2)))  &
          + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2)))
      fdiffd = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,1,-3)) + functn(OPS_ACC1(0,-1,-3)))  &
          + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,2,-3)) + functn(OPS_ACC1(0,-2,-3)))

    fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
          + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_4th_mixed

!a3d_p020_m024_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_4th_onesided(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1)))  &
          + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,2,-1)) + functn(OPS_ACC1(0,-2,-1)))
      fdiffb = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,1,-2)) + functn(OPS_ACC1(0,-1,-2)))  &
          + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2)))
      fdiffc = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,1,-3)) + functn(OPS_ACC1(0,-1,-3)))  &
          + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,2,-3)) + functn(OPS_ACC1(0,-2,-3)))
      fdiffd = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,1,-4)) + functn(OPS_ACC1(0,-1,-4)))  &
          + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,2,-4)) + functn(OPS_ACC1(0,-2,-4)))

      fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_4th_onesided

!a3d_p040_p004_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqA(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,1,-1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1))
      fdiffb = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,2,-2))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2))
      fdiffc = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,3,-3))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3))
      fdiffd = functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,4,-4))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-4))

      fderiv(OPS_ACC2(0,0,0)) = acc1yz*fdiffa + bcc1yz*fdiffb  &
          + ccc1yz*fdiffc + dcc1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqA

!a3d_p030_m010_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqB(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,-1,0))  - functn(OPS_ACC1(0,-1,1))  &
          - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,0,1))
      fdiffb = functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,1,-1))  &
          - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,0,-1))
      fdiffc = functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,2,-2))  &
          - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,0,-2))
      fdiffd = functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,3,-3))  &
          - functn(OPS_ACC1(0,0,0)) + functn(OPS_ACC1(0,0,-3))

      fderiv(OPS_ACC2(0,0,0)) = acc2yz*fdiffa + bcc2yz*fdiffb  &
          + ccc2yz*fdiffc + dcc2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqB

!a3d_p041_p003_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqC(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2yz*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,1,1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,1)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,1,-1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1)))  &
          + ccf2yz*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,1,-2))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2)))  &
          + dcf2yz*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,1,-3))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3)))
      fdiffb = acf2yz*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,2,1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,1)))  &
          + bcf2yz*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,2,-1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,2,-2))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2)))  &
          + dcf2yz*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,2,-3))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3)))
      fdiffc = acf2yz*(functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,3,1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,1)))  &
          + bcf2yz*(functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,3,-1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1)))  &
          + ccf2yz*(functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,3,-2))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,0)) - functn(OPS_ACC1(0,3,-3))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3)))
      fdiffd = acf2yz*(functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,4,1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,1)))  &
          + bcf2yz*(functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,4,-1))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-1)))  &
          + ccf2yz*(functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,4,-2))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-2)))  &
          + dcf2yz*(functn(OPS_ACC1(0,4,0)) - functn(OPS_ACC1(0,4,-3))  &
          - functn(OPS_ACC1(0,0,0))  + functn(OPS_ACC1(0,0,-3)))

      fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqC

!a3d_p030_m014_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqD(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2yz*(functn(OPS_ACC1(0,-1,0))   - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,-1))  + functn(OPS_ACC1(0,0,-1)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,0,-1)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,2,-1)) + functn(OPS_ACC1(0,0,-1)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,3,-1)) + functn(OPS_ACC1(0,0,-1)))
      fdiffb = acf2yz*(functn(OPS_ACC1(0,-1,0))   - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,-2))  + functn(OPS_ACC1(0,0,-2)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,1,-2)) + functn(OPS_ACC1(0,0,-2)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,0,-2)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,3,-2)) + functn(OPS_ACC1(0,0,-2)))
      fdiffc = acf2yz*(functn(OPS_ACC1(0,-1,0))   - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,-3))  + functn(OPS_ACC1(0,0,-3)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,1,-3)) + functn(OPS_ACC1(0,0,-3)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,2,-3)) + functn(OPS_ACC1(0,0,-3)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,3,-3)) + functn(OPS_ACC1(0,0,-3)))
      fdiffd = acf2yz*(functn(OPS_ACC1(0,-1,0))   - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,-1,-4))  + functn(OPS_ACC1(0,0,-4)))  &
          + bcf2yz*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,1,-4)) + functn(OPS_ACC1(0,0,-4)))  &
          + ccf2yz*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,2,-4)) + functn(OPS_ACC1(0,0,-4)))  &
          + dcf2yz*(functn(OPS_ACC1(0,3,0))  - functn(OPS_ACC1(0,0,0))  &
          - functn(OPS_ACC1(0,3,-4)) + functn(OPS_ACC1(0,0,-4)))

      fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqD

!a3d_p021_m023_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqE(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,1))  + functn(OPS_ACC1(0,-1,1)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,1))  + functn(OPS_ACC1(0,-2,1)))
        fdiffb = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-1)) + functn(OPS_ACC1(0,-2,-1)))
        fdiffc = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-2)) + functn(OPS_ACC1(0,-1,-2)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2)))
        fdiffd = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-3)) + functn(OPS_ACC1(0,-1,-3)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-3)) + functn(OPS_ACC1(0,-2,-3)))

        fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqE

!a3d_p020_m024_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqF(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-1)) + functn(OPS_ACC1(0,-2,-1)))
        fdiffb = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-2)) + functn(OPS_ACC1(0,-1,-2)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2)))
        fdiffc = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-3)) + functn(OPS_ACC1(0,-1,-3)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-3)) + functn(OPS_ACC1(0,-2,-3)))
        fdiffd = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-4)) + functn(OPS_ACC1(0,-1,-4)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-4)) + functn(OPS_ACC1(0,-2,-4)))

        fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqF

!a3d_p042_m002_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqG(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
            - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,1,-2))  &
            - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
        fdiffb = acofz1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,2,-1))  &
            - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
            - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
        fdiffc = acofz1*(functn(OPS_ACC1(0,3,1)) - functn(OPS_ACC1(0,3,-1))  &
            - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,3,2)) - functn(OPS_ACC1(0,3,-2))  &
            - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))
        fdiffd = acofz1*(functn(OPS_ACC1(0,4,1)) - functn(OPS_ACC1(0,4,-1))  &
            - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,4,2)) - functn(OPS_ACC1(0,4,-2))  &
            - functn(OPS_ACC1(0,0,2))  + functn(OPS_ACC1(0,0,-2)))

        fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqG

!a3d_p032_m012_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqH(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,-1,1))  - functn(OPS_ACC1(0,-1,-1))  &
            - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,-1,2))  - functn(OPS_ACC1(0,-1,-2))  &
            - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
        fdiffb = acofz1*(functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
            - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,1,2)) - functn(OPS_ACC1(0,1,-2))  &
            - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
        fdiffc = acofz1*(functn(OPS_ACC1(0,2,1)) - functn(OPS_ACC1(0,2,-1))  &
            - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
            - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))
        fdiffd = acofz1*(functn(OPS_ACC1(0,3,1)) - functn(OPS_ACC1(0,3,-1))  &
            - functn(OPS_ACC1(0,0,1)) + functn(OPS_ACC1(0,0,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,3,2)) - functn(OPS_ACC1(0,3,-2))  &
            - functn(OPS_ACC1(0,0,2)) + functn(OPS_ACC1(0,0,-2)))

        fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqH

SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqI(functn, fderiv, fstora, fstorb, nzglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv, fstora, fstorb
    real(kind=8) :: fdiffa, fdiffb

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nzglblm4

    integer :: js, ks, jks
    js = idx(2) -3
    ks = idx(3) -nzglblm4
    jks = js*3+ks +1

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

    fstora(OPS_ACC_MD3(jks,0,0,0)) = fdiffa
    fstorb(OPS_ACC_MD4(jks,0,0,0)) = fdiffb

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqI

SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqJ(functn, fderiv, fstora, fstorb, fstorc, nzglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb
    real(kind=8), dimension(1) :: fderiv, fstorc
    real(kind=8) :: fdiffc

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nzglblm4

    integer :: js, ks, jks
    js = idx(2) -3
    ks = idx(3) -nzglblm4
    jks = js*3+ks +1

    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
           - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4yz*fstora(OPS_ACC_MD3(jks,0,0,0)) + bcf4yz*fstorb(OPS_ACC_MD4(jks,0,0,0))  &
              + ccf4yz*fdiffc

    js = js-1
    jks = js*2+ks +1
    fstorc(OPS_ACC_MD5(jks,0,0,0)) = fdiffc

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqJ

SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqK(functn, fderiv, fstora, fstorb, fstorc, nzglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb, fstorc
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffd

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nzglblm4

    integer :: js, ks, jks1, jks2
    js = idx(2) -3
    ks = idx(3) -nzglblm4
    jks1 = js*3+ks +1

    js = js-1
    jks2 = js*2+ks +1

    fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
           - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5yz*fstora(OPS_ACC_MD3(jks1,0,0,0)) + bcf5yz*fstorb(OPS_ACC_MD4(jks1,0,0,0))  &
                            + ccf5yz*fstorc(OPS_ACC_MD5(jks2,0,0,0)) + dcf5yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqK

!a3d_p000_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqL(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,-1,-1)) - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0))
      fdiffb = functn(OPS_ACC1(0,-2,-2)) - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0))
      fdiffc = functn(OPS_ACC1(0,-3,-3)) - functn(OPS_ACC1(0,-3,0))  &
          - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0))
      fdiffd = functn(OPS_ACC1(0,-4,-4)) - functn(OPS_ACC1(0,-4,0))  &
          - functn(OPS_ACC1(0,0,-4))  + functn(OPS_ACC1(0,0,0))

      fderiv(OPS_ACC2(0,0,0)) = acc1yz*fdiffa + bcc1yz*fdiffb  &
          + ccc1yz*fdiffc + dcc1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqL

!a3d_p011_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqM(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = functn(OPS_ACC1(0,1,1))   - functn(OPS_ACC1(0,1,0))  &
          - functn(OPS_ACC1(0,0,1))  + functn(OPS_ACC1(0,0,0))
      fdiffb = functn(OPS_ACC1(0,-1,-1)) - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,0,-1)) + functn(OPS_ACC1(0,0,0))
      fdiffc = functn(OPS_ACC1(0,-2,-2)) - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,0,-2)) + functn(OPS_ACC1(0,0,0))
      fdiffd = functn(OPS_ACC1(0,-3,-3)) - functn(OPS_ACC1(0,-3,0))  &
          - functn(OPS_ACC1(0,0,-3)) + functn(OPS_ACC1(0,0,0))

      fderiv(OPS_ACC2(0,0,0)) = acc2yz*fdiffa + bcc2yz*fdiffb  &
          + ccc2yz*fdiffc + dcc2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqM

!a3d_p001_m043_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqN(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2yz*(functn(OPS_ACC1(0,-1,1))  - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,0,1))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,-1,-1)) - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,-1,-2)) - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,-1,-3)) - functn(OPS_ACC1(0,-1,0))  &
          - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0)))
      fdiffb = acf2yz*(functn(OPS_ACC1(0,-2,1))  - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,0,1))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,-2,-1)) - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,-2,-2)) - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,-2,-3)) - functn(OPS_ACC1(0,-2,0))  &
          - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0)))
      fdiffc = acf2yz*(functn(OPS_ACC1(0,-3,1))  - functn(OPS_ACC1(0,-3,0))  &
          - functn(OPS_ACC1(0,0,1))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,-3,-1)) - functn(OPS_ACC1(0,-3,0))  &
          - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,-3,-2)) - functn(OPS_ACC1(0,-3,0))  &
          - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,-3,-3)) - functn(OPS_ACC1(0,-3,0))  &
          - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0)))
      fdiffd = acf2yz*(functn(OPS_ACC1(0,-4,1))  - functn(OPS_ACC1(0,-4,0))  &
          - functn(OPS_ACC1(0,0,1))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,-4,-1)) - functn(OPS_ACC1(0,-4,0))  &
          - functn(OPS_ACC1(0,0,-1))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,-4,-2)) - functn(OPS_ACC1(0,-4,0))  &
          - functn(OPS_ACC1(0,0,-2))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,-4,-3)) - functn(OPS_ACC1(0,-4,0))  &
          - functn(OPS_ACC1(0,0,-3))  + functn(OPS_ACC1(0,0,0)))

      fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqN

!a3d_p010_m034_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqO(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acf2yz*(functn(OPS_ACC1(0,1,-1))  - functn(OPS_ACC1(0,0,-1))  &
          - functn(OPS_ACC1(0,1,0))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,-1,-1)) - functn(OPS_ACC1(0,0,-1))  &
          - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,-2,-1)) - functn(OPS_ACC1(0,0,-1))  &
          - functn(OPS_ACC1(0,-2,0))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,-3,-1)) - functn(OPS_ACC1(0,0,-1))  &
          - functn(OPS_ACC1(0,-3,0))  + functn(OPS_ACC1(0,0,0)))
      fdiffb = acf2yz*(functn(OPS_ACC1(0,1,-2))  - functn(OPS_ACC1(0,0,-2))  &
          - functn(OPS_ACC1(0,1,0))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,-1,-2)) - functn(OPS_ACC1(0,0,-2))  &
          - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,-2,-2)) - functn(OPS_ACC1(0,0,-2))  &
          - functn(OPS_ACC1(0,-2,0))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,-3,-2)) - functn(OPS_ACC1(0,0,-2))  &
          - functn(OPS_ACC1(0,-3,0))  + functn(OPS_ACC1(0,0,0)))
      fdiffc = acf2yz*(functn(OPS_ACC1(0,1,-3))  - functn(OPS_ACC1(0,0,-3))  &
          - functn(OPS_ACC1(0,1,0))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,-1,-3)) - functn(OPS_ACC1(0,0,-3))  &
          - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,-2,-3)) - functn(OPS_ACC1(0,0,-3))  &
          - functn(OPS_ACC1(0,-2,0))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,-3,-3)) - functn(OPS_ACC1(0,0,-3))  &
          - functn(OPS_ACC1(0,-3,0))  + functn(OPS_ACC1(0,0,0)))
      fdiffd = acf2yz*(functn(OPS_ACC1(0,1,-4))  - functn(OPS_ACC1(0,0,-4))  &
          - functn(OPS_ACC1(0,1,0))   + functn(OPS_ACC1(0,0,0)))  &
          + bcf2yz*(functn(OPS_ACC1(0,-1,-4)) - functn(OPS_ACC1(0,0,-4))  &
          - functn(OPS_ACC1(0,-1,0))  + functn(OPS_ACC1(0,0,0)))  &
          + ccf2yz*(functn(OPS_ACC1(0,-2,-4)) - functn(OPS_ACC1(0,0,-4))  &
          - functn(OPS_ACC1(0,-2,0))  + functn(OPS_ACC1(0,0,0)))  &
          + dcf2yz*(functn(OPS_ACC1(0,-3,-4)) - functn(OPS_ACC1(0,0,-4))  &
          - functn(OPS_ACC1(0,-3,0))  + functn(OPS_ACC1(0,0,0)))

      fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqO

!a3d_p021_m023_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqP(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,1))  + functn(OPS_ACC1(0,-1,1)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,1))  + functn(OPS_ACC1(0,-2,1)))
        fdiffb = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-1)) + functn(OPS_ACC1(0,-2,-1)))
        fdiffc = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-2)) + functn(OPS_ACC1(0,-1,-2)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2)))
        fdiffd = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-3)) + functn(OPS_ACC1(0,-1,-3)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-3)) + functn(OPS_ACC1(0,-2,-3)))

        fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqP

!a3d_p020_m024_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqQ(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-1)) + functn(OPS_ACC1(0,-2,-1)))
        fdiffb = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-2)) + functn(OPS_ACC1(0,-1,-2)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2)))
        fdiffc = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-3)) + functn(OPS_ACC1(0,-1,-3)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-3)) + functn(OPS_ACC1(0,-2,-3)))
        fdiffd = acofy1*(functn(OPS_ACC1(0,1,0))  - functn(OPS_ACC1(0,-1,0))  &
            - functn(OPS_ACC1(0,1,-4)) + functn(OPS_ACC1(0,-1,-4)))  &
            + bcofy1*(functn(OPS_ACC1(0,2,0))  - functn(OPS_ACC1(0,-2,0))  &
            - functn(OPS_ACC1(0,2,-4)) + functn(OPS_ACC1(0,-2,-4)))

        fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqQ

!a3d_p002_m042_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqR(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-1,2)) + functn(OPS_ACC1(0,-1,-2)))
        fdiffb = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-2,1)) + functn(OPS_ACC1(0,-2,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2)))
        fdiffc = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-3,1)) + functn(OPS_ACC1(0,-3,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-3,2)) + functn(OPS_ACC1(0,-3,-2)))
        fdiffd = acofz1*(functn(OPS_ACC1(0,0,1))  - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-4,1)) + functn(OPS_ACC1(0,-4,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2))  - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-4,2)) + functn(OPS_ACC1(0,-4,-2)))

        fderiv(OPS_ACC2(0,0,0)) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqR

!a3d_p012_m032_mixed_yz
SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqS(functn, fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,1,1))  + functn(OPS_ACC1(0,1,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,1,2))  + functn(OPS_ACC1(0,1,-2)))
        fdiffb = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-1,2)) + functn(OPS_ACC1(0,-1,-2)))
        fdiffc = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-2,1)) + functn(OPS_ACC1(0,-2,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2)))
        fdiffd = acofz1*(functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,-1))  &
            - functn(OPS_ACC1(0,-3,1)) + functn(OPS_ACC1(0,-3,-1)))  &
            + bcofz1*(functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,-2))  &
            - functn(OPS_ACC1(0,-3,2)) + functn(OPS_ACC1(0,-3,-2)))

        fderiv(OPS_ACC2(0,0,0)) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqS

SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqT(functn, fderiv, fstora, fstorb, nyglblm4, nzglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv, fstora, fstorb
    real(kind=8) :: fdiffa, fdiffb

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nyglblm4, nzglblm4

    integer :: js, ks, jks
    js = idx(2) -nyglblm4
    ks = idx(3) -nzglblm4
    jks = js*3+ks +1

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

    fstora(OPS_ACC_MD3(jks,0,0,0)) = fdiffa
    fstorb(OPS_ACC_MD4(jks,0,0,0)) = fdiffb

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqT

SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqU(functn, fderiv, fstora, fstorb, fstorc, nyglblm4, nzglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb
    real(kind=8), dimension(1) :: fderiv, fstorc
    real(kind=8) :: fdiffc

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nyglblm4, nzglblm4

    integer :: js, ks, jks
    js = idx(2) -nyglblm4
    ks = idx(3) -nzglblm4
    jks = js*3+ks +1

    fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
           - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))

    fderiv(OPS_ACC2(0,0,0)) = acf4yz*fstora(OPS_ACC_MD3(jks,0,0,0)) + bcf4yz*fstorb(OPS_ACC_MD4(jks,0,0,0))  &
              + ccf4yz*fdiffc

    jks = js*2+ks +1
    fstorc(OPS_ACC_MD5(jks,0,0,0)) = fdiffc

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqU

SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqV(functn, fderiv, fstora, fstorb, fstorc, nyglblm4, nzglblm4, idx)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: functn, fstora, fstorb, fstorc
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdiffd

    integer, dimension(3), intent(in) :: idx
    integer, intent(in) :: nyglblm4, nzglblm4

    integer :: js, ks, jks1, jks2
    js = idx(2) -nyglblm4
    ks = idx(3) -nzglblm4

    jks1 = js*3+ks +1
    jks2 = js*2+ks +1

    fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
           - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))

    fderiv(OPS_ACC2(0,0,0)) = acf5yz*fstora(OPS_ACC_MD3(jks1,0,0,0)) + bcf5yz*fstorb(OPS_ACC_MD4(jks1,0,0,0))  &
                            + ccf5yz*fstorc(OPS_ACC_MD5(jks2,0,0,0)) + dcf5yz*fdiffd

END SUBROUTINE d2fdyz_kernel_rh_zdir_corner_eqV

SUBROUTINE d2fdyz_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdely*ovdelz

END SUBROUTINE d2fdyz_kernel_scaling

