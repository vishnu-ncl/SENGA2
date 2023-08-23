SUBROUTINE d2fdxy_kernel_null(fderiv)

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_8

END SUBROUTINE d2fdxy_kernel_null

!a3d_p550_m550_mixed_xy
SUBROUTINE d2fdxy_kernel_interior(functn, fderiv)

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
SUBROUTINE d2fdxy_kernel_eqA(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqA

!a3d_p320_m120_mixed_xy
SUBROUTINE d2fdxy_kernel_eqB(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqB

!a3d_p220_m220_mixed_xy
SUBROUTINE d2fdxy_kernel_eqC(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
          - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
          - functn(OPS_ACC1(-2,2,0))  + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0))  = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_eqC

!a3d_p330_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_eqD(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqD

!a3d_p440_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_eqE(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqE

!a3d_p440_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_eqF(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqF

!a3d_p330_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_eqG(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqG

!a3d_p220_m220_mixed_xy
SUBROUTINE d2fdxy_kernel_eqH(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
          - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
          - functn(OPS_ACC1(-2,2,0))  + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0))  = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_eqH

!a3d_p120_m320_mixed_xy
SUBROUTINE d2fdxy_kernel_eqI(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqI

!a3d_p020_m420_mixed_xy
SUBROUTINE d2fdxy_kernel_eqJ(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqJ

!a3d_p240_m200_mixed_xy
SUBROUTINE d2fdxy_kernel_eqK(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqK

!a3d_p230_m210_mixed_xy
SUBROUTINE d2fdxy_kernel_eqL(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqL

!a3d_p220_m220_mixed_xy
SUBROUTINE d2fdxy_kernel_eqM(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0))  + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_eqM

!a3d_p330_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_eqN(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqN

!a3d_p440_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_eqO(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqO

!a3d_p440_p000_mixed_xy
SUBROUTINE d2fdxy_kernel_eqP(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqP

!a3d_p330_m110_mixed_xy
SUBROUTINE d2fdxy_kernel_eqQ(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqQ

!a3d_p430_m010_mixed_xy
SUBROUTINE d2fdxy_kernel_eqR(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqR

!a3d_p340_m100_mixed_xy
SUBROUTINE d2fdxy_kernel_eqS(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqS

!a3d_p240_m200_mixed_xy
SUBROUTINE d2fdxy_kernel_eqT(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqT

!a3d_p230_m210_mixed_xy
SUBROUTINE d2fdxy_kernel_eqU(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqU

!a3d_p420_m020_mixed_xy
SUBROUTINE d2fdxy_kernel_eqV(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqV

!a3d_p320_m120_mixed_xy
SUBROUTINE d2fdxy_kernel_eqW(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqW

SUBROUTINE d2fdxy_kernel_eqX(functn, fderiv, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: ic, jc

    ic = idx(1)
    jc = idx(2)

    IF (ic>=3 .and. ic<=5 .and. jc>=3 .and. jc<=5) THEN
        fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
               - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
        fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
               - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))

        fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (ic>=4 .and. ic<=5 .and. jc>=4 .and. jc<=5) THEN
        fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
               - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))

        fderiv(OPS_ACC2(0,0,0)) = acf4xy*fstora + bcf4xy*fstorb + ccf4xy*fdiffc

        fstorc = fdiffc
    END IF

    IF (ic==5 .and. jc==5) THEN
        fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
               - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,-4,0))

        fderiv(OPS_ACC2(0,0,0)) = acf5xy*fstora + bcf5xy*fstorb  &
                                + ccf5xy*fstorc + dcf5xy*fdiffd
    END IF

END SUBROUTINE d2fdxy_kernel_eqX

!a3d_p040_m400_mixed_xy
SUBROUTINE d2fdxy_kernel_eqY(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqY

!a3d_p100_m300_mixed_xy
SUBROUTINE d2fdxy_kernel_eqZ(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqZ

!a3d_p030_m410_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAA(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAA

!a3d_p140_m300_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAB(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAB

!a3d_p240_m200_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAC(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAC

!a3d_p230_m210_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAD(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAD

!a3d_p020_m420_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAE(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAE

!a3d_p120_m320_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAF(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAF

SUBROUTINE d2fdxy_kernel_eqAG(functn, fderiv, nxglbl, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4), intent(in) :: nxglbl

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: ic, jc

    ic = idx(1)
    jc = idx(2)

    IF (ic>=nxglbl-4 .and. ic<=nxglbl-2 .and. jc>=3 .and. jc<=5) THEN
        fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
               - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
        fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
               - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))

        fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (ic>=nxglbl-4 .and. ic<=nxglbl-3 .and. jc>=4 .and. jc<=5) THEN
        fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
               - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))

        fderiv(OPS_ACC2(0,0,0)) = acf4xy*fstora + bcf4xy*fstorb + ccf4xy*fdiffc

        fstorc = fdiffc
    END IF

    IF (ic==nxglbl-4 .and. jc==5) THEN
        fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
               - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,-4,0))

        fderiv(OPS_ACC2(0,0,0)) = acf5xy*fstora + bcf5xy*fstorb  &
                                + ccf5xy*fstorc + dcf5xy*fdiffd
    END IF

END SUBROUTINE d2fdxy_kernel_eqAG

!a3d_p440_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAH(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAH

!a3d_p330_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAI(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAI

!a3d_p220_m220_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAJ(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb

    fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(-1,1,0))  &
           - functn(OPS_ACC1(1,-1,0)) + functn(OPS_ACC1(-1,-1,0))
    fdiffb = functn(OPS_ACC1(2,2,0))  - functn(OPS_ACC1(-2,2,0))  &
           - functn(OPS_ACC1(2,-2,0)) + functn(OPS_ACC1(-2,-2,0))

    fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

END SUBROUTINE d2fdxy_kernel_eqAJ

!a3d_p210_m230_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAK(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAK

!a3d_p200_m240_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAL(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAL

!a3d_p400_p040_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAM(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAM

!a3d_p300_m100_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAN(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAN

!a3d_p410_p030_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAO(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAO

!a3d_p300_m140_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAP(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAP

!a3d_p210_m230_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAQ(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAQ

!a3d_p200_m240_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAR(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAR

!a3d_p420_m020_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAS(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAS

!a3d_p320_m120_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAT(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAT

SUBROUTINE d2fdxy_kernel_eqAU(functn, fderiv, nyglbl, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4), intent(in) :: nyglbl

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: ic, jc

    ic = idx(1)
    jc = idx(2)

    IF (ic>=3 .and. ic<=5 .and. jc>=nyglbl-4 .and. jc<=nyglbl-2) THEN
        fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
               - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
        fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
               - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))

        fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (ic>=4 .and. ic<=5 .and. jc>=nyglbl-4 .and. jc<=nyglbl-3) THEN
        fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
               - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))

        fderiv(OPS_ACC2(0,0,0)) = acf4xy*fstora + bcf4xy*fstorb + ccf4xy*fdiffc

        fstorc = fdiffc
    END IF
    IF (ic==5 .and. jc==nyglbl-4) THEN
        fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
               - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,-4,0))

        fderiv(OPS_ACC2(0,0,0)) = acf5xy*fstora + bcf5xy*fstorb  &
                                + ccf5xy*fstorc + dcf5xy*fdiffd
    END IF

END SUBROUTINE d2fdxy_kernel_eqAU

!a3d_p000_m440_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAV(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAV

!a3d_p110_m330_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAW(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAW

!a3d_p010_m430_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAX(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAX

!a3d_p100_m340_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAY(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAY

!a3d_p210_m230_mixed_xy
SUBROUTINE d2fdxy_kernel_eqAZ(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqAZ

!a3d_p200_m240_mixed_xy
SUBROUTINE d2fdxy_kernel_eqBA(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqBA

!a3d_p020_m420_mixed_xy
SUBROUTINE d2fdxy_kernel_eqBB(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqBB

!a3d_p120_m320_mixed_xy
SUBROUTINE d2fdxy_kernel_eqBC(functn, fderiv)

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

END SUBROUTINE d2fdxy_kernel_eqBC

SUBROUTINE d2fdxy_kernel_eqBD(functn, fderiv, nxglbl, nyglbl, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4), intent(in) :: nxglbl, nyglbl

    real(kind=8) :: fdiffa, fdiffb, fdiffc, fdiffd
    real(kind=8) :: fstora, fstorb, fstorc
    integer(kind=4) :: ic, jc

    ic = idx(1)
    jc = idx(2)

    IF (ic>=nxglbl-4 .and. ic<=nxglbl-2 .and. jc>=nyglbl-4 .and. jc<=nyglbl-2) THEN
        fdiffa = functn(OPS_ACC1(1,1,0)) - functn(OPS_ACC1(1,-1,0))  &
               - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))
        fdiffb = functn(OPS_ACC1(2,2,0)) - functn(OPS_ACC1(2,-2,0))  &
               - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0))

        fderiv(OPS_ACC2(0,0,0)) = acf3xy*fdiffa + bcf3xy*fdiffb

        fstora = fdiffa
        fstorb = fdiffb
    END IF

    IF (ic>=nxglbl-4 .and. ic<=nxglbl-3 .and. jc>=nyglbl-4 .and. jc<=nyglbl-3) THEN
        fdiffc = functn(OPS_ACC1(3,3,0)) - functn(OPS_ACC1(3,-3,0))  &
               - functn(OPS_ACC1(-3,3,0)) + functn(OPS_ACC1(-3,-3,0))

        fderiv(OPS_ACC2(0,0,0)) = acf4xy*fstora + bcf4xy*fstorb + ccf4xy*fdiffc

        fstorc = fdiffc
    END IF

    IF (ic==nxglbl-4 .and. jc==nyglbl-4) THEN
        fdiffd = functn(OPS_ACC1(4,4,0)) - functn(OPS_ACC1(4,-4,0))  &
               - functn(OPS_ACC1(-4,4,0)) + functn(OPS_ACC1(-4,-4,0))

        fderiv(OPS_ACC2(0,0,0)) = acf5xy*fstora + bcf5xy*fstorb  &
                            + ccf5xy*fstorc + dcf5xy*fdiffd
    END IF

END SUBROUTINE d2fdxy_kernel_eqBD

SUBROUTINE d2fdxy_kernel_scaling(fderiv)

    real(kind=8), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = fderiv(OPS_ACC1(0,0,0))*ovdelx*ovdely

END SUBROUTINE d2fdxy_kernel_scaling
