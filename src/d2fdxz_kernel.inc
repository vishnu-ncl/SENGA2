SUBROUTINE d2fdxz_kernel_null(fderiv)

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdxz_kernel_null

!a3d_p505_m505_mixed_xz
SUBROUTINE d2fdxz_kernel_interior(functn, fderiv)

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
SUBROUTINE d2fdxz_kernel_eqA(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqA

!a3d_p302_m102_mixed_xz
SUBROUTINE d2fdxz_kernel_eqB(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqB

!a3d_p202_m202_mixed_xz
SUBROUTINE d2fdxz_kernel_eqC(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(-2,0,2))  + functn(OPS_ACC1(-2,0,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

END SUBROUTINE d2fdxz_kernel_eqC

!a3d_p303_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_eqD(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqD

!a3d_p404_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_eqE(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqE

!a3d_p404_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_eqF(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqF

!a3d_p303_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_eqG(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqG

!a3d_p202_m202_mixed_xz
SUBROUTINE d2fdxz_kernel_eqH(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
           - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2))  - functn(OPS_ACC1(2,0,-2))  &
           - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

END SUBROUTINE d2fdxz_kernel_eqH

!a3d_p102_m302_mixed_xz
SUBROUTINE d2fdxz_kernel_eqI(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqI

!a3d_p002_m402_mixed_xz
SUBROUTINE d2fdxz_kernel_eqJ(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqJ

!a3d_p204_m200_mixed_xz
SUBROUTINE d2fdxz_kernel_eqK(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqK

!a3d_p203_m201_mixed_xz
SUBROUTINE d2fdxz_kernel_eqL(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqL

!a3d_p202_m202_mixed_xz
SUBROUTINE d2fdxz_kernel_eqM(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
           - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(-2,0,2))  &
           - functn(OPS_ACC1(2,0,-2))  + functn(OPS_ACC1(-2,0,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

END SUBROUTINE d2fdxz_kernel_eqM

!a3d_p303_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_eqN(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqN

!a3d_p404_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_eqO(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqO

!a3d_p404_p000_mixed_xz
SUBROUTINE d2fdxz_kernel_eqP(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqP

!a3d_p303_m101_mixed_xz
SUBROUTINE d2fdxz_kernel_eqQ(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqQ

!a3d_p403_m001_mixed_xz
SUBROUTINE d2fdxz_kernel_eqR(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqR

!a3d_p304_m100_mixed_xz
SUBROUTINE d2fdxz_kernel_eqS(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqS

!a3d_p204_m200_mixed_xz
SUBROUTINE d2fdxz_kernel_eqT(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqT

!a3d_p203_m201_mixed_xz
SUBROUTINE d2fdxz_kernel_eqU(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqU

!a3d_p402_m002_mixed_xz
SUBROUTINE d2fdxz_kernel_eqV(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqV

!a3d_p302_m102_mixed_xz
SUBROUTINE d2fdxz_kernel_eqW(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqW

SUBROUTINE d2fdxz_kernel_eqX(functn, fderiv, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3) :: idx

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: ic, kc

    ic = idx(1)
    kc = idx(3)

    IF (ic>=3 .and. ic<=5 .and. kc>=3 .and. kc<=5) THEN
        fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
               - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
        fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
               - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (ic>=4 .and. ic<=5 .and. kc>=4 .and. kc<=5) THEN
        fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,-3))  &
               - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acf4xz*fstora + bcf4xz*fstorb + ccf4xz*fdiffc

        fstorc = fdiffc
    END IF

    IF (ic==5 .and. kc==5) THEN
        fdiffd = functn(OPS_ACC1(4,0,4)) - functn(OPS_ACC1(4,0,-4))  &
               - functn(OPS_ACC1(-4,0,4)) + functn(OPS_ACC1(-4,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acf5xz*fstora + bcf5xz*fstorb  &
                                + ccf5xz*fstorc + dcf5xz*fdiffd
    END IF

END SUBROUTINE d2fdxz_kernel_eqX

!a3d_p004_m400_mixed_xz
SUBROUTINE d2fdxz_kernel_eqY(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqY

!a3d_p100_m300_mixed_xz
SUBROUTINE d2fdxz_kernel_eqZ(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqZ

!a3d_p003_m401_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAA(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAA

!a3d_p104_m300_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAB(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAB

!a3d_p204_m200_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAC(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAC

!a3d_p203_m201_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAD(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAD

!a3d_p002_m402_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAE(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAE

!a3d_p102_m302_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAF(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAF

SUBROUTINE d2fdxz_kernel_eqAG(functn, fderiv, nxglbl, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4), intent(in) :: nxglbl

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: ic, kc

    ic = idx(1)
    kc = idx(3)

    IF (ic>=nxglbl-4 .and. ic<=nxglbl-2 .and. kc>=3 .and. kc<=5) THEN
        fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
               - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
        fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
               - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (ic>=nxglbl-4 .and. ic<=nxglbl-3 .and. kc>=4 .and. kc<=5) THEN
        fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,-3))  &
               - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acf4xz*fstora + bcf4xz*fstorb + ccf4xz*fdiffc

        fstorc = fdiffc
    END IF

    IF (ic==nxglbl-4 .and. kc==5) THEN
        fdiffd = functn(OPS_ACC1(4,0,4)) - functn(OPS_ACC1(4,0,-4))  &
               - functn(OPS_ACC1(-4,0,4)) + functn(OPS_ACC1(-4,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acf5xz*fstora + bcf5xz*fstorb  &
                                + ccf5xz*fstorc + dcf5xz*fdiffd
    END IF

END SUBROUTINE d2fdxz_kernel_eqAG

!a3d_p404_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAH(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAH

!a3d_p303_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAI(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAI

!a3d_p202_m202_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAJ(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(-1,0,1))  &
           - functn(OPS_ACC1(1,0,-1)) + functn(OPS_ACC1(-1,0,-1))
    fdiffb = functn(OPS_ACC1(2,0,2))  - functn(OPS_ACC1(-2,0,2))  &
           - functn(OPS_ACC1(2,0,-2)) + functn(OPS_ACC1(-2,0,-2))

    fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

END SUBROUTINE d2fdxz_kernel_eqAJ

!a3d_p201_m203_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAK(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAK

!a3d_p200_m204_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAL(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAL

!a3d_p400_p004_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAM(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAM

!a3d_p300_m100_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAN(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAN

!a3d_p401_p003_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAO(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAO

!a3d_p300_m104_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAP(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAP

!a3d_p201_m203_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAQ(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAQ

!a3d_p200_m204_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAR(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAR

!a3d_p402_m002_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAS(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAS

!a3d_p302_m102_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAT(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAT

SUBROUTINE d2fdxz_kernel_eqAU(functn, fderiv, nzglbl, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4), intent(in) :: nzglbl

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: ic, kc

    ic = idx(1)
    kc = idx(3)

    IF (ic>=3 .and. ic<=5 .and. kc>=nzglbl-4 .and. kc<=nzglbl-2) THEN
        fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
               - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
        fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
               - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (ic>=4 .and. ic<=5 .and. kc>=nzglbl-4 .and. kc<=nzglbl-3) THEN
        fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,-3))  &
               - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acf4xz*fstora + bcf4xz*fstorb + ccf4xz*fdiffc

        fstorc = fdiffc
    END IF

    IF (ic==5 .and. kc==nzglbl-4) THEN
        fdiffd = functn(OPS_ACC1(4,0,4)) - functn(OPS_ACC1(4,0,-4))  &
               - functn(OPS_ACC1(-4,0,4)) + functn(OPS_ACC1(-4,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acf5xz*fstora + bcf5xz*fstorb  &
                                + ccf5xz*fstorc + dcf5xz*fdiffd
    END IF

END SUBROUTINE d2fdxz_kernel_eqAU

!a3d_p000_m404_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAV(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAV

!a3d_p101_m303_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAW(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAW

!a3d_p001_m403_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAX(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAX

!a3d_p100_m304_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAY(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAY

!a3d_p201_m203_mixed_xz
SUBROUTINE d2fdxz_kernel_eqAZ(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqAZ

!a3d_p200_m204_mixed_xz
SUBROUTINE d2fdxz_kernel_eqBA(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqBA

!a3d_p002_m402_mixed_xz
SUBROUTINE d2fdxz_kernel_eqBB(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqBB

!a3d_p102_m302_mixed_xz
SUBROUTINE d2fdxz_kernel_eqBC(functn, fderiv)

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

END SUBROUTINE d2fdxz_kernel_eqBC

SUBROUTINE d2fdxz_kernel_eqBD(functn, fderiv, nxglbl, nzglbl, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4), intent(in) :: nxglbl, nzglbl

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: ic, kc

    ic = idx(1)
    kc = idx(3)

    IF (ic>=nxglbl-4 .and. ic<=nxglbl-2 .and. kc>=nzglbl-4 .and. kc<=nzglbl-2) THEN
        fdiffa = functn(OPS_ACC1(1,0,1)) - functn(OPS_ACC1(1,0,-1))  &
               - functn(OPS_ACC1(-1,0,1)) + functn(OPS_ACC1(-1,0,-1))
        fdiffb = functn(OPS_ACC1(2,0,2)) - functn(OPS_ACC1(2,0,-2))  &
               - functn(OPS_ACC1(-2,0,2)) + functn(OPS_ACC1(-2,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acf3xz*fdiffa + bcf3xz*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (ic>=nxglbl-4 .and. ic<=nxglbl-3 .and. kc>=nzglbl-4 .and. kc<=nzglbl-3) THEN
        fdiffc = functn(OPS_ACC1(3,0,3)) - functn(OPS_ACC1(3,0,-3))  &
               - functn(OPS_ACC1(-3,0,3)) + functn(OPS_ACC1(-3,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acf4xz*fstora + bcf4xz*fstorb + ccf4xz*fdiffc

        fstorc = fdiffc
    END IF

    IF (ic==nxglbl-4 .and. kc==nzglbl-4) THEN
        fdiffd = functn(OPS_ACC1(4,0,4)) - functn(OPS_ACC1(4,0,-4))  &
               - functn(OPS_ACC1(-4,0,4)) + functn(OPS_ACC1(-4,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acf5xz*fstora + bcf5xz*fstorb  &
                                + ccf5xz*fstorc + dcf5xz*fdiffd
    END IF

END SUBROUTINE d2fdxz_kernel_eqBD

SUBROUTINE d2fdxz_kernel_scaling(fderiv)

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdelx*ovdelz

END SUBROUTINE d2fdxz_kernel_scaling

