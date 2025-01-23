!HEAT FLUX: Thermal fxlcnw
SUBROUTINE hf_kernel_eqA(erhs,store1,store7,deltax)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store1, store7
    real(kind=8), intent(in) :: deltax
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC3(0,0,0))*store1(OPS_ACC2(0,0,0))
    fornow1 = store7(OPS_ACC3(1,0,0))*store1(OPS_ACC2(1,0,0))
    fornow2 = store7(OPS_ACC3(2,0,0))*store1(OPS_ACC2(2,0,0))
    fornow3 = store7(OPS_ACC3(3,0,0))*store1(OPS_ACC2(3,0,0))
    fornow4 = store7(OPS_ACC3(4,0,0))*store1(OPS_ACC2(4,0,0))

    fornow = 0.0_8
    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1 &
            - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3 &
            - (1.0_8/4.0_8)*fornow4) / deltax

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqA

!HEAT FLUX: Thermal fxrcnw
SUBROUTINE hf_kernel_eqB(erhs,store1,store7,deltax)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store1, store7
    real(kind=8), intent(in) :: deltax
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC3(0,0,0))*store1(OPS_ACC2(0,0,0))
    fornow1 = store7(OPS_ACC3(-1,0,0))*store1(OPS_ACC2(-1,0,0))
    fornow2 = store7(OPS_ACC3(-2,0,0))*store1(OPS_ACC2(-2,0,0))
    fornow3 = store7(OPS_ACC3(-3,0,0))*store1(OPS_ACC2(-3,0,0))
    fornow4 = store7(OPS_ACC3(-4,0,0))*store1(OPS_ACC2(-4,0,0))

    fornow = 0.0_8
    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1 &
            + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3 &
            + (1.0_8/4.0_8)*fornow4) / deltax

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqB

!HEAT FLUX: Thermal fylcnw
SUBROUTINE hf_kernel_eqC(erhs,store2,store7,deltay)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store2, store7
    real(kind=8), intent(in) :: deltay
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC3(0,0,0))*store2(OPS_ACC2(0,0,0))
    fornow1 = store7(OPS_ACC3(0,1,0))*store2(OPS_ACC2(0,1,0))
    fornow2 = store7(OPS_ACC3(0,2,0))*store2(OPS_ACC2(0,2,0))
    fornow3 = store7(OPS_ACC3(0,3,0))*store2(OPS_ACC2(0,3,0))
    fornow4 = store7(OPS_ACC3(0,4,0))*store2(OPS_ACC2(0,4,0))

    fornow = 0.0_8
    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1 &
            - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3 &
            - (1.0_8/4.0_8)*fornow4) / deltay

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqC

!HEAT FLUX: Thermal fyrcnw
SUBROUTINE hf_kernel_eqD(erhs,store2,store7,deltay)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store2, store7
    real(kind=8), intent(in) :: deltay
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC3(0,0,0))*store2(OPS_ACC2(0,0,0))
    fornow1 = store7(OPS_ACC3(0,-1,0))*store2(OPS_ACC2(0,-1,0))
    fornow2 = store7(OPS_ACC3(0,-2,0))*store2(OPS_ACC2(0,-2,0))
    fornow3 = store7(OPS_ACC3(0,-3,0))*store2(OPS_ACC2(0,-3,0))
    fornow4 = store7(OPS_ACC3(0,-4,0))*store2(OPS_ACC2(0,-4,0))

    fornow = 0.0_8
    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1 &
            + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3 &
            + (1.0_8/4.0_8)*fornow4) / deltay

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqD

!HEAT FLUX: Thermal fzlcnw
SUBROUTINE hf_kernel_eqE(erhs,store3,store7,deltaz)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store3, store7
    real(kind=8), intent(in) :: deltaz
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC3(0,0,0))*store3(OPS_ACC2(0,0,0))
    fornow1 = store7(OPS_ACC3(0,0,1))*store3(OPS_ACC2(0,0,1))
    fornow2 = store7(OPS_ACC3(0,0,2))*store3(OPS_ACC2(0,0,2))
    fornow3 = store7(OPS_ACC3(0,0,3))*store3(OPS_ACC2(0,0,3))
    fornow4 = store7(OPS_ACC3(0,0,4))*store3(OPS_ACC2(0,0,4))

    fornow = 0.0_8
    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1 &
            - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3 &
            - (1.0_8/4.0_8)*fornow4) / deltaz

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqE

!HEAT FLUX: Thermal fzrcnw
SUBROUTINE hf_kernel_eqF(erhs,store3,store7,deltaz)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store3, store7
    real(kind=8), intent(in) :: deltaz
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC3(0,0,0))*store3(OPS_ACC2(0,0,0))
    fornow1 = store7(OPS_ACC3(0,0,-1))*store3(OPS_ACC2(0,0,-1))
    fornow2 = store7(OPS_ACC3(0,0,-2))*store3(OPS_ACC2(0,0,-2))
    fornow3 = store7(OPS_ACC3(0,0,-3))*store3(OPS_ACC2(0,0,-3))
    fornow4 = store7(OPS_ACC3(0,0,-4))*store3(OPS_ACC2(0,0,-4))

    fornow = 0.0_8
    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1 &
            + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3 &
            + (1.0_8/4.0_8)*fornow4) / deltaz

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqF

SUBROUTINE hf_kernel_eqG(erhs,trun,tdrmix,store7,td1x,rgspec_ispec,deltax)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    real(kind=8), intent(in) :: rgspec_ispec,deltax

    real(kind=8) :: combo1
    real(kind=8) :: fornow
    real(kind=8) :: fornow0, fornow1, fornow2, fornow3, fornow4

    fornow = 0.0_8

    combo1 = trun(OPS_ACC2(0,0,0)) * tdrmix(OPS_ACC3(0,0,0))
    fornow0 = combo1 * store7(OPS_ACC4(0,0,0))*td1x(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(1,0,0)) * tdrmix(OPS_ACC3(1,0,0))
    fornow1 = combo1 * store7(OPS_ACC4(1,0,0)) * td1x(OPS_ACC5(1,0,0))

    combo1 = trun(OPS_ACC2(2,0,0)) * tdrmix(OPS_ACC3(2,0,0))
    fornow2 = combo1 * store7(OPS_ACC4(2,0,0)) * td1x(OPS_ACC5(2,0,0))

    combo1 = trun(OPS_ACC2(3,0,0)) * tdrmix(OPS_ACC3(3,0,0))
    fornow3 = combo1 * store7(OPS_ACC4(3,0,0)) * td1x(OPS_ACC5(3,0,0))

    combo1 = trun(OPS_ACC2(4,0,0)) * tdrmix(OPS_ACC3(4,0,0))
    fornow4 = combo1 * store7(OPS_ACC4(4,0,0)) * td1x(OPS_ACC5(4,0,0))

    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1  &
               - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3  &
               - (1.0_8/4.0_8)*fornow4)/deltax

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec_ispec*fornow

END SUBROUTINE hf_kernel_eqG

SUBROUTINE hf_kernel_eqH(erhs,trun,tdrmix,store7,td1x,rgspec_ispec,deltax)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    real(kind=8), intent(in) :: rgspec_ispec,deltax

    real(kind=8) :: combo1
    real(kind=8) :: fornow
    real(kind=8) :: fornow0, fornow1, fornow2, fornow3, fornow4

    fornow = 0.0_8

    combo1 = trun(OPS_ACC2(0,0,0)) * tdrmix(OPS_ACC3(0,0,0))
    fornow0 = combo1 * store7(OPS_ACC4(0,0,0)) * td1x(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(-1,0,0)) * tdrmix(OPS_ACC3(-1,0,0))
    fornow1 = combo1 * store7(OPS_ACC4(-1,0,0)) * td1x(OPS_ACC5(-1,0,0))

    combo1 = trun(OPS_ACC2(-2,0,0)) * tdrmix(OPS_ACC3(-2,0,0))
    fornow2 = combo1 * store7(OPS_ACC4(-2,0,0)) * td1x(OPS_ACC5(-2,0,0))

    combo1 = trun(OPS_ACC2(-3,0,0)) * tdrmix(OPS_ACC3(-3,0,0))
    fornow3 = combo1 * store7(OPS_ACC4(-3,0,0)) * td1x(OPS_ACC5(-3,0,0))

    combo1 = trun(OPS_ACC2(-4,0,0)) * tdrmix(OPS_ACC3(-4,0,0))
    fornow4 = combo1 * store7(OPS_ACC4(-4,0,0)) * td1x(OPS_ACC5(-4,0,0))

    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1  &
              + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3  &
              + (1.0_8/4.0_8)*fornow4)/deltax

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec_ispec*fornow

END SUBROUTINE hf_kernel_eqH

SUBROUTINE hf_kernel_eqI(erhs,trun,tdrmix,store7,td1y,rgspec_ispec,deltay)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    real(kind=8), intent(in) :: rgspec_ispec,deltay

    real(kind=8) :: combo1
    real(kind=8) :: fornow
    real(kind=8) :: fornow0, fornow1, fornow2, fornow3, fornow4

    fornow = 0.0_8

    combo1 = trun(OPS_ACC2(0,0,0)) * tdrmix(OPS_ACC3(0,0,0))
    fornow0 = combo1 * store7(OPS_ACC4(0,0,0)) * td1y(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(0,1,0)) * tdrmix(OPS_ACC3(0,1,0))
    fornow1 = combo1 * store7(OPS_ACC4(0,1,0)) * td1y(OPS_ACC5(0,1,0))

    combo1 = trun(OPS_ACC2(0,2,0)) * tdrmix(OPS_ACC3(0,2,0))
    fornow2 = combo1 * store7(OPS_ACC4(0,2,0)) * td1y(OPS_ACC5(0,2,0))

    combo1 = trun(OPS_ACC2(0,3,0)) * tdrmix(OPS_ACC3(0,3,0))
    fornow3 = combo1 * store7(OPS_ACC4(0,3,0)) * td1y(OPS_ACC5(0,3,0))

    combo1 = trun(OPS_ACC2(0,4,0)) * tdrmix(OPS_ACC3(0,4,0))
    fornow4 = combo1 * store7(OPS_ACC4(0,4,0)) * td1y(OPS_ACC5(0,4,0))

    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1  &
               - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3  &
               - (1.0_8/4.0_8)*fornow4)/deltay

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec_ispec*fornow

END SUBROUTINE hf_kernel_eqI

SUBROUTINE hf_kernel_eqJ(erhs,trun,tdrmix,store7,td1y,rgspec_ispec,deltay)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    real(kind=8), intent(in) :: rgspec_ispec,deltay

    real(kind=8) :: combo1
    real(kind=8) :: fornow
    real(kind=8) :: fornow0, fornow1, fornow2, fornow3, fornow4

    fornow = 0.0_8

    combo1 = trun(OPS_ACC2(0,0,0)) * tdrmix(OPS_ACC3(0,0,0))
    fornow0 = combo1 * store7(OPS_ACC4(0,0,0)) * td1y(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(0,-1,0)) * tdrmix(OPS_ACC3(0,-1,0))
    fornow1 = combo1 * store7(OPS_ACC4(0,-1,0)) * td1y(OPS_ACC5(0,-1,0))

    combo1 = trun(OPS_ACC2(0,-2,0)) * tdrmix(OPS_ACC3(0,-2,0))
    fornow2 = combo1 * store7(OPS_ACC4(0,-2,0)) * td1y(OPS_ACC5(0,-2,0))

    combo1 = trun(OPS_ACC2(0,-3,0)) * tdrmix(OPS_ACC3(0,-3,0))
    fornow3 = combo1 * store7(OPS_ACC4(0,-3,0)) * td1y(OPS_ACC5(0,-3,0))

    combo1 = trun(OPS_ACC2(0,-4,0)) * tdrmix(OPS_ACC3(0,-4,0))
    fornow4 = combo1 * store7(OPS_ACC4(0,-4,0)) * td1y(OPS_ACC5(0,-4,0))

    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1  &
              + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3  &
              + (1.0_8/4.0_8)*fornow4)/deltay

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec_ispec*fornow

END SUBROUTINE hf_kernel_eqJ

SUBROUTINE hf_kernel_eqK(erhs,trun,tdrmix,store7,td1z,rgspec_ispec,deltaz)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    real(kind=8), intent(in) :: rgspec_ispec,deltaz

    real(kind=8) :: combo1
    real(kind=8) :: fornow
    real(kind=8) :: fornow0, fornow1, fornow2, fornow3, fornow4

    fornow = 0.0_8

    combo1 = trun(OPS_ACC2(0,0,0)) * tdrmix(OPS_ACC3(0,0,0))
    fornow0 = combo1 * store7(OPS_ACC4(0,0,0)) * td1z(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(0,0,1)) * tdrmix(OPS_ACC3(0,0,1))
    fornow1 = combo1 * store7(OPS_ACC4(0,0,1)) * td1z(OPS_ACC5(0,0,1))

    combo1 = trun(OPS_ACC2(0,0,2)) * tdrmix(OPS_ACC3(0,0,2))
    fornow2 = combo1 * store7(OPS_ACC4(0,0,2)) * td1z(OPS_ACC5(0,0,2))

    combo1 = trun(OPS_ACC2(0,0,3)) * tdrmix(OPS_ACC3(0,0,3))
    fornow3 = combo1 * store7(OPS_ACC4(0,0,3)) * td1z(OPS_ACC5(0,0,3))

    combo1 = trun(OPS_ACC2(0,0,4)) * tdrmix(OPS_ACC3(0,0,4))
    fornow4 = combo1 * store7(OPS_ACC4(0,0,4)) * td1z(OPS_ACC5(0,0,4))

    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1  &
               - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3  &
               - (1.0_8/4.0_8)*fornow4)/deltaz

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec_ispec*fornow

END SUBROUTINE hf_kernel_eqK

SUBROUTINE hf_kernel_eqL(erhs,trun,tdrmix,store7,td1z,rgspec_ispec,deltaz)

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    real(kind=8), intent(in) :: rgspec_ispec,deltaz

    real(kind=8) :: combo1
    real(kind=8) :: fornow
    real(kind=8) :: fornow0, fornow1, fornow2, fornow3, fornow4

    fornow = 0.0_8

    combo1 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))
    fornow0 = combo1*store7(OPS_ACC4(0,0,0))*td1z(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(0,0,-1))*tdrmix(OPS_ACC3(0,0,-1))
    fornow1 = combo1*store7(OPS_ACC4(0,0,-1)) *td1z(OPS_ACC5(0,0,-1))

    combo1 = trun(OPS_ACC2(0,0,-2))*tdrmix(OPS_ACC3(0,0,-2))
    fornow2 = combo1*store7(OPS_ACC4(0,0,-2)) *td1z(OPS_ACC5(0,0,-2))

    combo1 = trun(OPS_ACC2(0,0,-3))*tdrmix(OPS_ACC3(0,0,-3))
    fornow3 = combo1*store7(OPS_ACC4(0,0,-3)) *td1z(OPS_ACC5(0,0,-3))

    combo1 = trun(OPS_ACC2(0,0,-4))*tdrmix(OPS_ACC3(0,0,-4))
    fornow4 = combo1*store7(OPS_ACC4(0,0,-4)) *td1z(OPS_ACC5(0,0,-4))

    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1  &
              + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3  &
              + (1.0_8/4.0_8)*fornow4)/deltaz

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec_ispec*fornow

END SUBROUTINE hf_kernel_eqL

!HEAT FLUX: Enthalpy fxldfw 
SUBROUTINE hf_kernel_eqM(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1x,deltax,rgspec_ispec)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1x
    real(kind=8), intent(in) :: deltax
    real(kind=8), intent(in) :: rgspec_ispec
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC7(0,0,0))*td1x(OPS_ACC8(0,0,0))
    fornow1 = store7(OPS_ACC7(1,0,0))*td1x(OPS_ACC8(1,0,0))
    fornow2 = store7(OPS_ACC7(2,0,0))*td1x(OPS_ACC8(2,0,0))
    fornow3 = store7(OPS_ACC7(3,0,0))*td1x(OPS_ACC8(3,0,0))
    fornow4 = store7(OPS_ACC7(4,0,0))*td1x(OPS_ACC8(4,0,0))

    fornow = 0.0_8
    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1 &
             - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3 &
             - (1.0_8/4.0_8)*fornow4) / deltax

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec_ispec*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqM

!HEAT FLUX: Enthalpy fxrdfw
SUBROUTINE hf_kernel_eqN(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1x,deltax,rgspec_ispec)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1x
    real(kind=8), intent(in) :: deltax
    real(kind=8), intent(in) :: rgspec_ispec
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC7(0,0,0))*td1x(OPS_ACC8(0,0,0))
    fornow1 = store7(OPS_ACC7(-1,0,0))*td1x(OPS_ACC8(-1,0,0))
    fornow2 = store7(OPS_ACC7(-2,0,0))*td1x(OPS_ACC8(-2,0,0))
    fornow3 = store7(OPS_ACC7(-3,0,0))*td1x(OPS_ACC8(-3,0,0))
    fornow4 = store7(OPS_ACC7(-4,0,0))*td1x(OPS_ACC8(-4,0,0))

    fornow = 0.0_8
    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1 &
             + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3 &
             + (1.0_8/4.0_8)*fornow4) / deltax

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec_ispec*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqN

!HEAT FLUX: Enthalpy fyldfw
SUBROUTINE hf_kernel_eqO(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1y,deltay,rgspec_ispec)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1y
    real(kind=8), intent(in) :: deltay
    real(kind=8), intent(in) :: rgspec_ispec
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC7(0,0,0))*td1y(OPS_ACC8(0,0,0))
    fornow1 = store7(OPS_ACC7(0,1,0))*td1y(OPS_ACC8(0,1,0))
    fornow2 = store7(OPS_ACC7(0,2,0))*td1y(OPS_ACC8(0,2,0))
    fornow3 = store7(OPS_ACC7(0,3,0))*td1y(OPS_ACC8(0,3,0))
    fornow4 = store7(OPS_ACC7(0,4,0))*td1y(OPS_ACC8(0,4,0))

    fornow = 0.0_8
    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1 &
            - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3 &
            - (1.0_8/4.0_8)*fornow4) / deltay

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec_ispec*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqO

!HEAT FLUX: Enthalpy fyrdfw
SUBROUTINE hf_kernel_eqP(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1y,deltay,rgspec_ispec)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1y
    real(kind=8), intent(in) :: deltay
    real(kind=8), intent(in) :: rgspec_ispec
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC7(0,0,0))*td1y(OPS_ACC8(0,0,0))
    fornow1 = store7(OPS_ACC7(0,-1,0))*td1y(OPS_ACC8(0,-1,0))
    fornow2 = store7(OPS_ACC7(0,-2,0))*td1y(OPS_ACC8(0,-2,0))
    fornow3 = store7(OPS_ACC7(0,-3,0))*td1y(OPS_ACC8(0,-3,0))
    fornow4 = store7(OPS_ACC7(0,-4,0))*td1y(OPS_ACC8(0,-4,0))

    fornow = 0.0_8
    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1 &
             + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3 &
             + (1.0_8/4.0_8)*fornow4) / deltay

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec_ispec*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqP

!HEAT FLUX: Enthalpy fzldfw
SUBROUTINE hf_kernel_eqQ(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1z,deltaz,rgspec_ispec)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1z
    real(kind=8), intent(in) :: deltaz
    real(kind=8), intent(in) :: rgspec_ispec
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC7(0,0,0))*td1z(OPS_ACC8(0,0,0))
    fornow1 = store7(OPS_ACC7(0,0,1))*td1z(OPS_ACC8(0,0,1))
    fornow2 = store7(OPS_ACC7(0,0,2))*td1z(OPS_ACC8(0,0,2))
    fornow3 = store7(OPS_ACC7(0,0,3))*td1z(OPS_ACC8(0,0,3))
    fornow4 = store7(OPS_ACC7(0,0,4))*td1z(OPS_ACC8(0,0,4))

    fornow = 0.0_8
    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1 &
             - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3 &
             - (1.0_8/4.0_8)*fornow4) / deltaz

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec_ispec*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqQ

!HEAT FLUX: Enthalpy fzrdfw
SUBROUTINE hf_kernel_eqR(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1z,deltaz,rgspec_ispec)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1z
    real(kind=8), intent(in) :: deltaz
    real(kind=8), intent(in) :: rgspec_ispec
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC7(0,0,0))*td1z(OPS_ACC8(0,0,0))
    fornow1 = store7(OPS_ACC7(0,0,-1))*td1z(OPS_ACC8(0,0,-1))
    fornow2 = store7(OPS_ACC7(0,0,-2))*td1z(OPS_ACC8(0,0,-2))
    fornow3 = store7(OPS_ACC7(0,0,-3))*td1z(OPS_ACC8(0,0,-3))
    fornow4 = store7(OPS_ACC7(0,0,-4))*td1z(OPS_ACC8(0,0,-4))

    fornow = 0.0_8
    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1 &
             + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3 &
             + (1.0_8/4.0_8)*fornow4) / deltaz

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec_ispec*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqR

!HEAT FLUX: Enthalpy fxldfw
SUBROUTINE hf_kernel_eqS(rate,erhs,vtmp,utmp,store7,arr1,deltax)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), intent(in) :: deltax
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC5(0,0,0))*arr1(OPS_ACC6(0,0,0))
    fornow1 = store7(OPS_ACC5(1,0,0))*arr1(OPS_ACC6(1,0,0))
    fornow2 = store7(OPS_ACC5(2,0,0))*arr1(OPS_ACC6(2,0,0))
    fornow3 = store7(OPS_ACC5(3,0,0))*arr1(OPS_ACC6(3,0,0))
    fornow4 = store7(OPS_ACC5(4,0,0))*arr1(OPS_ACC6(4,0,0))

    fornow = 0.0_8
    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1 &
            - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3 &
            - (1.0_8/4.0_8)*fornow4) / deltax

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqS

!HEAT FLUX: Enthalpy fxrdfw
SUBROUTINE hf_kernel_eqT(rate,erhs,vtmp,utmp,store7,arr1,deltax)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), intent(in) :: deltax
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC5(0,0,0))*arr1(OPS_ACC6(0,0,0))
    fornow1 = store7(OPS_ACC5(-1,0,0))*arr1(OPS_ACC6(-1,0,0))
    fornow2 = store7(OPS_ACC5(-2,0,0))*arr1(OPS_ACC6(-2,0,0))
    fornow3 = store7(OPS_ACC5(-3,0,0))*arr1(OPS_ACC6(-3,0,0))
    fornow4 = store7(OPS_ACC5(-4,0,0))*arr1(OPS_ACC6(-4,0,0))

    fornow = 0.0_8
    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1 &
            + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3 &
            + (1.0_8/4.0_8)*fornow4) / deltax

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqT

!HEAT FLUX: Enthalpy fxrdfw
SUBROUTINE hf_kernel_eqU(rate,erhs,vtmp,utmp,store7,arr1,deltay)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), intent(in) :: deltay
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC5(0,0,0))*arr1(OPS_ACC6(0,0,0))
    fornow1 = store7(OPS_ACC5(0,1,0))*arr1(OPS_ACC6(0,1,0))
    fornow2 = store7(OPS_ACC5(0,2,0))*arr1(OPS_ACC6(0,2,0))
    fornow3 = store7(OPS_ACC5(0,3,0))*arr1(OPS_ACC6(0,3,0))
    fornow4 = store7(OPS_ACC5(0,4,0))*arr1(OPS_ACC6(0,4,0))

    fornow = 0.0_8
    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1 &
            - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3 &
            - (1.0_8/4.0_8)*fornow4) / deltay

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqU

!HEAT FLUX: Enthalpy fyrdfw
SUBROUTINE hf_kernel_eqV(rate,erhs,vtmp,utmp,store7,arr1,deltay)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), intent(in) :: deltay
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC5(0,0,0))*arr1(OPS_ACC6(0,0,0))
    fornow1 = store7(OPS_ACC5(0,-1,0))*arr1(OPS_ACC6(0,-1,0))
    fornow2 = store7(OPS_ACC5(0,-2,0))*arr1(OPS_ACC6(0,-2,0))
    fornow3 = store7(OPS_ACC5(0,-3,0))*arr1(OPS_ACC6(0,-3,0))
    fornow4 = store7(OPS_ACC5(0,-4,0))*arr1(OPS_ACC6(0,-4,0))

    fornow = 0.0_8
    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1 &
            + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3 &
            + (1.0_8/4.0_8)*fornow4) / deltay

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqV

!HEAT FLUX: Enthalpy fzldfw
SUBROUTINE hf_kernel_eqW(rate,erhs,vtmp,utmp,store7,arr1,deltaz)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), intent(in) :: deltaz
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC5(0,0,0))*arr1(OPS_ACC6(0,0,0))
    fornow1 = store7(OPS_ACC5(0,0,1))*arr1(OPS_ACC6(0,0,1))
    fornow2 = store7(OPS_ACC5(0,0,2))*arr1(OPS_ACC6(0,0,2))
    fornow3 = store7(OPS_ACC5(0,0,3))*arr1(OPS_ACC6(0,0,3))
    fornow4 = store7(OPS_ACC5(0,0,4))*arr1(OPS_ACC6(0,0,4))

    fornow = 0.0_8
    fornow = ((-25.0_8/12.0_8)*fornow0 + (4.0_8)*fornow1 &
            - (3.0_8)*fornow2 + (4.0_8/3.0_8)*fornow3 &
            - (1.0_8/4.0_8)*fornow4) / deltaz

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqW

!HEAT FLUX: Enthalpy fzrdfw
SUBROUTINE hf_kernel_eqX(rate,erhs,vtmp,utmp,store7,arr1,deltaz)

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), intent(in) :: deltaz
    real(kind=8) :: fornow
    real(kind=8) :: fornow0,fornow1,fornow2,fornow3,fornow4

    fornow0 = store7(OPS_ACC5(0,0,0))*arr1(OPS_ACC6(0,0,0))
    fornow1 = store7(OPS_ACC5(0,0,-1))*arr1(OPS_ACC6(0,0,-1))
    fornow2 = store7(OPS_ACC5(0,0,-2))*arr1(OPS_ACC6(0,0,-2))
    fornow3 = store7(OPS_ACC5(0,0,-3))*arr1(OPS_ACC6(0,0,-3))
    fornow4 = store7(OPS_ACC5(0,0,-4))*arr1(OPS_ACC6(0,0,-4))

    fornow = 0.0_8
    fornow = ((25.0_8/12.0_8)*fornow0 - (4.0_8)*fornow1 &
            + (3.0_8)*fornow2 - (4.0_8/3.0_8)*fornow3 &
            + (1.0_8/4.0_8)*fornow4) / deltaz

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqX
