SUBROUTINE d2fdyz_kernel_null(fderiv)

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdyz_kernel_null

!a3d_p055_m055_mixed_yz
SUBROUTINE d2fdyz_kernel_interior(functn, fderiv)

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
SUBROUTINE d2fdyz_kernel_eqA(functn, fderiv)

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


END SUBROUTINE d2fdyz_kernel_eqA

!a3d_p032_m012_mixed_yz
SUBROUTINE d2fdyz_kernel_eqB(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqB

!a3d_p022_m022_mixed_yz
SUBROUTINE d2fdyz_kernel_eqC(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2))  + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

END SUBROUTINE d2fdyz_kernel_eqC

!a3d_p033_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_eqD(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqD

!a3d_p044_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_eqE(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqE

!a3d_p044_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_eqF(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqF

!a3d_p033_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_eqG(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqG

!a3d_p022_m022_mixed_yz
SUBROUTINE d2fdyz_kernel_eqH(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
           - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2))  - functn(OPS_ACC1(0,2,-2))  &
           - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

END SUBROUTINE d2fdyz_kernel_eqH

!a3d_p012_m032_mixed_yz
SUBROUTINE d2fdyz_kernel_eqI(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqI

!a3d_p002_m042_mixed_yz
SUBROUTINE d2fdyz_kernel_eqJ(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqJ

!a3d_p024_m020_mixed_yz
SUBROUTINE d2fdyz_kernel_eqK(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqK

!a3d_p023_m021_mixed_yz
SUBROUTINE d2fdyz_kernel_eqL(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqL

!a3d_p022_m022_mixed_yz
SUBROUTINE d2fdyz_kernel_eqM(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
           - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1))
    fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,-2,2))  &
           - functn(OPS_ACC1(0,2,-2))  + functn(OPS_ACC1(0,-2,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

END SUBROUTINE d2fdyz_kernel_eqM

!a3d_p033_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_eqN(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqN

!a3d_p044_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_eqO(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqO

!a3d_p044_p000_mixed_yz
SUBROUTINE d2fdyz_kernel_eqP(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqP

!a3d_p033_m011_mixed_yz
SUBROUTINE d2fdyz_kernel_eqQ(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqQ

!a3d_p043_m001_mixed_yz
SUBROUTINE d2fdyz_kernel_eqR(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqR

!a3d_p034_m010_mixed_yz
SUBROUTINE d2fdyz_kernel_eqS(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqS

!a3d_p024_m020_mixed_yz
SUBROUTINE d2fdyz_kernel_eqT(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqT

!a3d_p023_m021_mixed_yz
SUBROUTINE d2fdyz_kernel_eqU(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqU

!a3d_p042_m002_mixed_yz
SUBROUTINE d2fdyz_kernel_eqV(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqV

!a3d_p032_m012_mixed_yz
SUBROUTINE d2fdyz_kernel_eqW(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqW

SUBROUTINE d2fdyz_kernel_eqX(functn, fderiv, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: jc, kc

    jc = idx(2)
    kc = idx(3)

    IF (jc>=3 .and. jc<=5 .and. kc>=3 .and. kc<=5) THEN
        fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
               - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
        fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
               - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

        fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (jc>=4 .and. jc<=5 .and. kc>=4 .and. kc<=5) THEN
        fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
               - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))

        fderiv(OPS_ACC2(0,0,0)) = acf4yz*fstora + bcf4yz*fstorb + ccf4yz*fdiffc

        fstorc = fdiffc
    END IF

    IF (jc==5 .and. kc==5) THEN
        fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
               - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))

        fderiv(OPS_ACC2(0,0,0)) = acf5yz*fstora + bcf5yz*fstorb  &
                                + ccf5yz*fstorc + dcf5yz*fdiffd
    END IF

END SUBROUTINE d2fdyz_kernel_eqX

!a3d_p004_m040_mixed_yz
SUBROUTINE d2fdyz_kernel_eqY(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqY

!a3d_p010_m030_mixed_yz
SUBROUTINE d2fdyz_kernel_eqZ(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqZ

!a3d_p003_m041_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAA(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAA

!a3d_p014_m030_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAB(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAB

!a3d_p024_m020_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAC(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAC

!a3d_p023_m021_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAD(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAD

!a3d_p002_m042_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAE(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAE

!a3d_p012_m032_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAF(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAF

SUBROUTINE d2fdyz_kernel_eqAG(functn, fderiv, nyglbl, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4), intent(in) :: nyglbl

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: jc, kc

    jc = idx(2)
    kc = idx(3)

    IF (jc>=nyglbl-4 .and. jc<=nyglbl-2 .and. kc>=3 .and. kc<=5) THEN
        fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
               - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
        fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
               - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

        fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (jc>=nyglbl-4 .and. jc<=nyglbl-3 .and. kc>=4 .and. kc<=5) THEN
        fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
               - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))

        fderiv(OPS_ACC2(0,0,0)) = acf4yz*fstora + bcf4yz*fstorb + ccf4yz*fdiffc

        fstorc = fdiffc
    END IF

    IF (jc==nyglbl-4 .and. kc==5) THEN
        fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
               - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))

        fderiv(OPS_ACC2(0,0,0)) = acf5yz*fstora + bcf5yz*fstorb  &
                                + ccf5yz*fstorc + dcf5yz*fdiffd
    END IF

END SUBROUTINE d2fdyz_kernel_eqAG

!a3d_p044_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAH(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAH

!a3d_p033_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAI(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAI

!a3d_p022_m022_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAJ(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,-1,1))  &
          - functn(OPS_ACC1(0,1,-1)) + functn(OPS_ACC1(0,-1,-1))
      fdiffb = functn(OPS_ACC1(0,2,2))  - functn(OPS_ACC1(0,-2,2))  &
          - functn(OPS_ACC1(0,2,-2)) + functn(OPS_ACC1(0,-2,-2))

      fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

END SUBROUTINE d2fdyz_kernel_eqAJ

!a3d_p021_m023_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAK(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAK

!a3d_p020_m024_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAL(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAL

!a3d_p040_p004_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAM(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAM

!a3d_p030_m010_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAN(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAN

!a3d_p041_p003_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAO(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAO

!a3d_p030_m014_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAP(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAP

!a3d_p021_m023_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAQ(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAQ

!a3d_p020_m024_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAR(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAR

!a3d_p042_m002_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAS(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAS

!a3d_p032_m012_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAT(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAT

SUBROUTINE d2fdyz_kernel_eqAU(functn, fderiv, nzglbl, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4), intent(in) :: nzglbl

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: jc, kc

    jc = idx(2)
    kc = idx(3)

    IF (jc>=3 .and. jc<=5 .and. kc>=nzglbl-4 .and. kc<=nzglbl-2) THEN
        fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
               - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
        fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
               - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

        fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (jc>=4 .and. jc<=5 .and. kc>=nzglbl-4 .and. kc<=nzglbl-3) THEN
        fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
               - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))

        fderiv(OPS_ACC2(0,0,0)) = acf4yz*fstora + bcf4yz*fstorb + ccf4yz*fdiffc

        fstorc = fdiffc
    END IF

    IF (jc==5 .and. kc==nzglbl-4) THEN
        fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
               - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))

        fderiv(OPS_ACC2(0,0,0)) = acf5yz*fstora + bcf5yz*fstorb  &
                                + ccf5yz*fstorc + dcf5yz*fdiffd
    END IF

END SUBROUTINE d2fdyz_kernel_eqAU

!a3d_p000_m044_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAV(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAV

!a3d_p011_m033_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAW(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAW

!a3d_p001_m043_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAX(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAX

!a3d_p010_m034_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAY(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAY

!a3d_p021_m023_mixed_yz
SUBROUTINE d2fdyz_kernel_eqAZ(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqAZ

!a3d_p020_m024_mixed_yz
SUBROUTINE d2fdyz_kernel_eqBA(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqBA

!a3d_p002_m042_mixed_yz
SUBROUTINE d2fdyz_kernel_eqBB(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqBB

!a3d_p012_m032_mixed_yz
SUBROUTINE d2fdyz_kernel_eqBC(functn, fderiv)

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

END SUBROUTINE d2fdyz_kernel_eqBC

SUBROUTINE d2fdyz_kernel_eqBD(functn, fderiv, nyglbl, nzglbl, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4), intent(in) :: nyglbl, nzglbl

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: jc, kc

    jc = idx(2)
    kc = idx(3)

    IF (jc>=nyglbl-4 .and. jc<=nyglbl-2 .and. kc>=nzglbl-4 .and. kc<=nzglbl-2) THEN
        fdiffa = functn(OPS_ACC1(0,1,1)) - functn(OPS_ACC1(0,1,-1))  &
               - functn(OPS_ACC1(0,-1,1)) + functn(OPS_ACC1(0,-1,-1))
        fdiffb = functn(OPS_ACC1(0,2,2)) - functn(OPS_ACC1(0,2,-2))  &
               - functn(OPS_ACC1(0,-2,2)) + functn(OPS_ACC1(0,-2,-2))

        fderiv(OPS_ACC2(0,0,0)) = acf3yz*fdiffa + bcf3yz*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (jc>=nyglbl-4 .and. jc<=nyglbl-3 .and. kc>=nzglbl-4 .and. kc<=nzglbl-3) THEN
        fdiffc = functn(OPS_ACC1(0,3,3)) - functn(OPS_ACC1(0,3,-3))  &
               - functn(OPS_ACC1(0,-3,3)) + functn(OPS_ACC1(0,-3,-3))

        fderiv(OPS_ACC2(0,0,0)) = acf4yz*fstora + bcf4yz*fstorb + ccf4yz*fdiffc

        fstorc = fdiffc
    END IF

    IF (jc==nyglbl-4 .and. kc==nzglbl-4) THEN
        fdiffd = functn(OPS_ACC1(0,4,4)) - functn(OPS_ACC1(0,4,-4))  &
               - functn(OPS_ACC1(0,-4,4)) + functn(OPS_ACC1(0,-4,-4))

        fderiv(OPS_ACC2(0,0,0)) = acf5yz*fstora + bcf5yz*fstorb  &
                                + ccf5yz*fstorc + dcf5yz*fdiffd
    END IF

END SUBROUTINE d2fdyz_kernel_eqBD

SUBROUTINE d2fdyz_kernel_scaling(fderiv)

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdely*ovdelz

END SUBROUTINE d2fdyz_kernel_scaling
