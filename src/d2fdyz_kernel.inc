SUBROUTINE d2fdyz_kernel_null(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdyz_kernel_null

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

SUBROUTINE d2fdyz_kernel_scaling(fderiv)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdely*ovdelz

END SUBROUTINE d2fdyz_kernel_scaling

