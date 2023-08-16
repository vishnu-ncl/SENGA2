!HEAT FLUX: Thermal fxlcnw
SUBROUTINE hf_kernel_eqA(erhs,store1,store7,acbcxl)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store1, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8) :: fornow

    fornow = acbcxl(1)*store7(OPS_ACC3(1,0,0))*store1(OPS_ACC2(1,0,0)) &
           + acbcxl(2)*store7(OPS_ACC3(2,0,0))*store1(OPS_ACC2(2,0,0)) &
           + acbcxl(3)*store7(OPS_ACC3(3,0,0))*store1(OPS_ACC2(3,0,0)) &
           + acbcxl(4)*store7(OPS_ACC3(4,0,0))*store1(OPS_ACC2(4,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqA

!HEAT FLUX: Thermal fxrcnw
SUBROUTINE hf_kernel_eqB(erhs,store1,store7,acbcxr)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store1, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8) :: fornow

    fornow = acbcxr(4)*store7(OPS_ACC3(-4,0,0))*store1(OPS_ACC2(-4,0,0)) &
           + acbcxr(3)*store7(OPS_ACC3(-3,0,0))*store1(OPS_ACC2(-3,0,0)) &
           + acbcxr(2)*store7(OPS_ACC3(-2,0,0))*store1(OPS_ACC2(-2,0,0)) &
           + acbcxr(1)*store7(OPS_ACC3(-1,0,0))*store1(OPS_ACC2(-1,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqB

!HEAT FLUX: Thermal fylcnw
SUBROUTINE hf_kernel_eqC(erhs,store2,store7,acbcyl)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store2, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8) :: fornow

    fornow = acbcyl(1)*store7(OPS_ACC3(0,1,0))*store2(OPS_ACC2(0,1,0)) &
           + acbcyl(2)*store7(OPS_ACC3(0,2,0))*store2(OPS_ACC2(0,2,0)) &
           + acbcyl(3)*store7(OPS_ACC3(0,3,0))*store2(OPS_ACC2(0,3,0)) &
           + acbcyl(4)*store7(OPS_ACC3(0,4,0))*store2(OPS_ACC2(0,4,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqC

!HEAT FLUX: Thermal fyrcnw
SUBROUTINE hf_kernel_eqD(erhs,store2,store7,acbcyr)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store2, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8) :: fornow

    fornow = acbcyr(4)*store7(OPS_ACC3(0,-4,0))*store2(OPS_ACC2(0,-4,0)) &
           + acbcyr(3)*store7(OPS_ACC3(0,-3,0))*store2(OPS_ACC2(0,-3,0)) &
           + acbcyr(2)*store7(OPS_ACC3(0,-2,0))*store2(OPS_ACC2(0,-2,0)) &
           + acbcyr(1)*store7(OPS_ACC3(0,-1,0))*store2(OPS_ACC2(0,-1,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqD

!HEAT FLUX: Thermal fzlcnw
SUBROUTINE hf_kernel_eqE(erhs,store3,store7,acbczl)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store3, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8) :: fornow

    fornow = acbczl(1)*store7(OPS_ACC3(0,0,1))*store3(OPS_ACC2(0,0,1)) &
           + acbczl(2)*store7(OPS_ACC3(0,0,2))*store3(OPS_ACC2(0,0,2)) &
           + acbczl(3)*store7(OPS_ACC3(0,0,3))*store3(OPS_ACC2(0,0,3)) &
           + acbczl(4)*store7(OPS_ACC3(0,0,4))*store3(OPS_ACC2(0,0,4))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqE

!HEAT FLUX: Thermal fzrcnw
SUBROUTINE hf_kernel_eqF(erhs,store3,store7,acbczr)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: store3, store7
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8) :: fornow

    fornow = acbczr(4)*store7(OPS_ACC3(0,0,-4))*store3(OPS_ACC2(0,0,-4)) &
           + acbczr(3)*store7(OPS_ACC3(0,0,-3))*store3(OPS_ACC2(0,0,-3)) &
           + acbczr(2)*store7(OPS_ACC3(0,0,-2))*store3(OPS_ACC2(0,0,-2)) &
           + acbczr(1)*store7(OPS_ACC3(0,0,-1))*store3(OPS_ACC2(0,0,-1))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + fornow

END SUBROUTINE hf_kernel_eqF

!HEAT FLUX: Soret and DUFOUR fxlcnw
SUBROUTINE hf_kernel_eqG(erhs,trun,tdrmix,store7,td1x,acbcxl,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(1,0,0))*tdrmix(OPS_ACC3(1,0,0))
    combo2 = trun(OPS_ACC2(2,0,0))*tdrmix(OPS_ACC3(2,0,0))
    combo3 = trun(OPS_ACC2(3,0,0))*tdrmix(OPS_ACC3(3,0,0))
    combo4 = trun(OPS_ACC2(4,0,0))*tdrmix(OPS_ACC3(4,0,0))

    fornow = acbcxl(1)*combo1*store7(OPS_ACC4(1,0,0))*td1x(OPS_ACC5(1,0,0)) &
           + acbcxl(2)*combo2*store7(OPS_ACC4(2,0,0))*td1x(OPS_ACC5(2,0,0)) &
           + acbcxl(3)*combo3*store7(OPS_ACC4(3,0,0))*td1x(OPS_ACC5(3,0,0)) &
           + acbcxl(4)*combo4*store7(OPS_ACC4(4,0,0))*td1x(OPS_ACC5(4,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqG

!HEAT FLUX: Soret and DUFOUR fxrcnw
SUBROUTINE hf_kernel_eqH(erhs,trun,tdrmix,store7,td1x,acbcxr,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(-4,0,0))*tdrmix(OPS_ACC3(-4,0,0))
    combo2 = trun(OPS_ACC2(-3,0,0))*tdrmix(OPS_ACC3(-3,0,0))
    combo3 = trun(OPS_ACC2(-2,0,0))*tdrmix(OPS_ACC3(-2,0,0))
    combo4 = trun(OPS_ACC2(-1,0,0))*tdrmix(OPS_ACC3(-1,0,0))

    fornow = acbcxr(4)*combo1*store7(OPS_ACC4(-4,0,0))*td1x(OPS_ACC5(-4,0,0)) &
           + acbcxr(3)*combo2*store7(OPS_ACC4(-3,0,0))*td1x(OPS_ACC5(-3,0,0)) &
           + acbcxr(2)*combo3*store7(OPS_ACC4(-2,0,0))*td1x(OPS_ACC5(-2,0,0)) &
           + acbcxr(1)*combo4*store7(OPS_ACC4(-1,0,0))*td1x(OPS_ACC5(-1,0,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqH

!HEAT FLUX: Soret and DUFOUR fylcnw
SUBROUTINE hf_kernel_eqI(erhs,trun,tdrmix,store7,td1y,acbcyl,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(0,1,0))*tdrmix(OPS_ACC3(0,1,0))
    combo2 = trun(OPS_ACC2(0,2,0))*tdrmix(OPS_ACC3(0,2,0))
    combo3 = trun(OPS_ACC2(0,3,0))*tdrmix(OPS_ACC3(0,3,0))
    combo4 = trun(OPS_ACC2(0,4,0))*tdrmix(OPS_ACC3(0,4,0))

    fornow = acbcyl(1)*combo1*store7(OPS_ACC4(0,1,0))*td1y(OPS_ACC5(0,1,0)) &
           + acbcyl(2)*combo2*store7(OPS_ACC4(0,2,0))*td1y(OPS_ACC5(0,2,0)) &
           + acbcyl(3)*combo3*store7(OPS_ACC4(0,3,0))*td1y(OPS_ACC5(0,3,0)) &
           + acbcyl(4)*combo4*store7(OPS_ACC4(0,4,0))*td1y(OPS_ACC5(0,4,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqI

!HEAT FLUX: Soret and DUFOUR fyrcnw
SUBROUTINE hf_kernel_eqJ(erhs,trun,tdrmix,store7,td1y,acbcyr,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(0,-4,0))*tdrmix(OPS_ACC3(0,-4,0))
    combo2 = trun(OPS_ACC2(0,-3,0))*tdrmix(OPS_ACC3(0,-3,0))
    combo3 = trun(OPS_ACC2(0,-2,0))*tdrmix(OPS_ACC3(0,-2,0))
    combo4 = trun(OPS_ACC2(0,-1,0))*tdrmix(OPS_ACC3(0,-1,0))

    fornow = acbcyr(4)*combo1*store7(OPS_ACC4(0,-4,0))*td1y(OPS_ACC5(0,-4,0)) &
           + acbcyr(3)*combo2*store7(OPS_ACC4(0,-3,0))*td1y(OPS_ACC5(0,-3,0)) &
           + acbcyr(2)*combo3*store7(OPS_ACC4(0,-2,0))*td1y(OPS_ACC5(0,-2,0)) &
           + acbcyr(1)*combo4*store7(OPS_ACC4(0,-1,0))*td1y(OPS_ACC5(0,-1,0))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqJ

!HEAT FLUX: Soret and DUFOUR fzlcnw
SUBROUTINE hf_kernel_eqK(erhs,trun,tdrmix,store7,td1z,acbczl,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(0,0,1))*tdrmix(OPS_ACC3(0,0,1))
    combo2 = trun(OPS_ACC2(0,0,2))*tdrmix(OPS_ACC3(0,0,2))
    combo3 = trun(OPS_ACC2(0,0,3))*tdrmix(OPS_ACC3(0,0,3))
    combo4 = trun(OPS_ACC2(0,0,4))*tdrmix(OPS_ACC3(0,0,4))

    fornow = acbczl(1)*combo1*store7(OPS_ACC4(0,0,1))*td1z(OPS_ACC5(0,0,1)) &
           + acbczl(2)*combo2*store7(OPS_ACC4(0,0,2))*td1z(OPS_ACC5(0,0,2)) &
           + acbczl(3)*combo3*store7(OPS_ACC4(0,0,3))*td1z(OPS_ACC5(0,0,3)) &
           + acbczl(4)*combo4*store7(OPS_ACC4(0,0,4))*td1z(OPS_ACC5(0,0,4))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqK

!HEAT FLUX: Soret and DUFOUR fzrcnw
SUBROUTINE hf_kernel_eqL(erhs,trun,tdrmix,store7,td1z,acbczr,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow
    real(kind=8) :: combo1, combo2, combo3, combo4

    combo1 = trun(OPS_ACC2(0,0,-4))*tdrmix(OPS_ACC3(0,0,-4))
    combo2 = trun(OPS_ACC2(0,0,-3))*tdrmix(OPS_ACC3(0,0,-3))
    combo3 = trun(OPS_ACC2(0,0,-2))*tdrmix(OPS_ACC3(0,0,-2))
    combo4 = trun(OPS_ACC2(0,0,-1))*tdrmix(OPS_ACC3(0,0,-1))

    fornow = acbczr(4)*combo1*store7(OPS_ACC4(0,0,-4))*td1z(OPS_ACC5(0,0,-4)) &
           + acbczr(3)*combo2*store7(OPS_ACC4(0,0,-3))*td1z(OPS_ACC5(0,0,-3)) &
           + acbczr(2)*combo3*store7(OPS_ACC4(0,0,-2))*td1z(OPS_ACC5(0,0,-2)) &
           + acbczr(1)*combo4*store7(OPS_ACC4(0,0,-1))*td1z(OPS_ACC5(0,0,-1))

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqL

!HEAT FLUX: Enthalpy fxldfw 
SUBROUTINE hf_kernel_eqM(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1x,acbcxl,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1x
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbcxl(1)*store7(OPS_ACC7(1,0,0))*td1x(OPS_ACC8(1,0,0)) &
           + acbcxl(2)*store7(OPS_ACC7(2,0,0))*td1x(OPS_ACC8(2,0,0)) &
           + acbcxl(3)*store7(OPS_ACC7(3,0,0))*td1x(OPS_ACC8(3,0,0)) &
           + acbcxl(4)*store7(OPS_ACC7(4,0,0))*td1x(OPS_ACC8(4,0,0))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqM

!HEAT FLUX: Enthalpy fxrdfw
SUBROUTINE hf_kernel_eqN(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1x,acbcxr,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1x
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbcxr(4)*store7(OPS_ACC7(-4,0,0))*td1x(OPS_ACC8(-4,0,0)) &
           + acbcxr(3)*store7(OPS_ACC7(-3,0,0))*td1x(OPS_ACC8(-3,0,0)) &
           + acbcxr(2)*store7(OPS_ACC7(-2,0,0))*td1x(OPS_ACC8(-2,0,0)) &
           + acbcxr(1)*store7(OPS_ACC7(-1,0,0))*td1x(OPS_ACC8(-1,0,0))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqN

!HEAT FLUX: Enthalpy fyldfw
SUBROUTINE hf_kernel_eqO(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1y,acbcyl,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1y
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbcyl(1)*store7(OPS_ACC7(0,1,0))*td1y(OPS_ACC8(0,1,0)) &
           + acbcyl(2)*store7(OPS_ACC7(0,2,0))*td1y(OPS_ACC8(0,2,0)) &
           + acbcyl(3)*store7(OPS_ACC7(0,3,0))*td1y(OPS_ACC8(0,3,0)) &
           + acbcyl(4)*store7(OPS_ACC7(0,4,0))*td1y(OPS_ACC8(0,4,0))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqO

!HEAT FLUX: Enthalpy fyrdfw
SUBROUTINE hf_kernel_eqP(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1y,acbcyr,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1y
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbcyr(4)*store7(OPS_ACC7(0,-4,0))*td1y(OPS_ACC8(0,-4,0)) &
           + acbcyr(3)*store7(OPS_ACC7(0,-3,0))*td1y(OPS_ACC8(0,-3,0)) &
           + acbcyr(2)*store7(OPS_ACC7(0,-2,0))*td1y(OPS_ACC8(0,-2,0)) &
           + acbcyr(1)*store7(OPS_ACC7(0,-1,0))*td1y(OPS_ACC8(0,-1,0))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqP

!HEAT FLUX: Enthalpy fzldfw
SUBROUTINE hf_kernel_eqQ(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1z,acbczl,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1z
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbczl(1)*store7(OPS_ACC7(0,0,1))*td1z(OPS_ACC8(0,0,1)) &
           + acbczl(2)*store7(OPS_ACC7(0,0,2))*td1z(OPS_ACC8(0,0,2)) &
           + acbczl(3)*store7(OPS_ACC7(0,0,3))*td1z(OPS_ACC8(0,0,3)) &
           + acbczl(4)*store7(OPS_ACC7(0,0,4))*td1z(OPS_ACC8(0,0,4))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqQ

!HEAT FLUX: Enthalpy fzrdfw
SUBROUTINE hf_kernel_eqR(rate,erhs,vtmp,trun,tdrmix,utmp,store7,td1z,acbczr,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,utmp,store7,td1z
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow

    fornow = acbczr(4)*store7(OPS_ACC7(0,0,-4))*td1z(OPS_ACC8(0,0,-4)) &
           + acbczr(3)*store7(OPS_ACC7(0,0,-3))*td1z(OPS_ACC8(0,0,-3)) &
           + acbczr(2)*store7(OPS_ACC7(0,0,-2))*td1z(OPS_ACC8(0,0,-2)) &
           + acbczr(1)*store7(OPS_ACC7(0,0,-1))*td1z(OPS_ACC8(0,0,-1))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*(utmp(OPS_ACC6(0,0,0))  &
                         + rgspec(ispec)*trun(OPS_ACC4(0,0,0))*tdrmix(OPS_ACC5(0,0,0)))

END SUBROUTINE hf_kernel_eqR

!HEAT FLUX: Enthalpy fxldfw
SUBROUTINE hf_kernel_eqS(rate,erhs,vtmp,utmp,store7,arr1,acbcxl)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8) :: fornow

    fornow = acbcxl(1)*store7(OPS_ACC5(1,0,0))*arr1(OPS_ACC6(1,0,0)) &
           + acbcxl(2)*store7(OPS_ACC5(2,0,0))*arr1(OPS_ACC6(2,0,0)) &
           + acbcxl(3)*store7(OPS_ACC5(3,0,0))*arr1(OPS_ACC6(3,0,0)) &
           + acbcxl(4)*store7(OPS_ACC5(4,0,0))*arr1(OPS_ACC6(4,0,0))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqS

!HEAT FLUX: Enthalpy fxrdfw
SUBROUTINE hf_kernel_eqT(rate,erhs,vtmp,utmp,store7,arr1,acbcxr)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8) :: fornow

    fornow = acbcxr(4)*store7(OPS_ACC5(-4,0,0))*arr1(OPS_ACC6(-4,0,0)) &
           + acbcxr(3)*store7(OPS_ACC5(-3,0,0))*arr1(OPS_ACC6(-3,0,0)) &
           + acbcxr(2)*store7(OPS_ACC5(-2,0,0))*arr1(OPS_ACC6(-2,0,0)) &
           + acbcxr(1)*store7(OPS_ACC5(-1,0,0))*arr1(OPS_ACC6(-1,0,0))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqT

!HEAT FLUX: Enthalpy fxrdfw
SUBROUTINE hf_kernel_eqU(rate,erhs,vtmp,utmp,store7,arr1,acbcyl)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8) :: fornow

    fornow = acbcyl(1)*store7(OPS_ACC5(0,1,0))*arr1(OPS_ACC6(0,1,0)) &
           + acbcyl(2)*store7(OPS_ACC5(0,2,0))*arr1(OPS_ACC6(0,2,0)) &
           + acbcyl(3)*store7(OPS_ACC5(0,3,0))*arr1(OPS_ACC6(0,3,0)) &
           + acbcyl(4)*store7(OPS_ACC5(0,4,0))*arr1(OPS_ACC6(0,4,0))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqU

!HEAT FLUX: Enthalpy fyrdfw
SUBROUTINE hf_kernel_eqV(rate,erhs,vtmp,utmp,store7,arr1,acbcyr)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8) :: fornow

    fornow = acbcyr(4)*store7(OPS_ACC5(0,-4,0))*arr1(OPS_ACC6(0,-4,0)) &
           + acbcyr(3)*store7(OPS_ACC5(0,-3,0))*arr1(OPS_ACC6(0,-3,0)) &
           + acbcyr(2)*store7(OPS_ACC5(0,-2,0))*arr1(OPS_ACC6(0,-2,0)) &
           + acbcyr(1)*store7(OPS_ACC5(0,-1,0))*arr1(OPS_ACC6(0,-1,0))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqV

!HEAT FLUX: Enthalpy fzldfw
SUBROUTINE hf_kernel_eqW(rate,erhs,vtmp,utmp,store7,arr1,acbczl)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8) :: fornow

    fornow = acbczl(1)*store7(OPS_ACC5(0,0,1))*arr1(OPS_ACC6(0,0,1)) &
           + acbczl(2)*store7(OPS_ACC5(0,0,2))*arr1(OPS_ACC6(0,0,2)) &
           + acbczl(3)*store7(OPS_ACC5(0,0,3))*arr1(OPS_ACC6(0,0,3)) &
           + acbczl(4)*store7(OPS_ACC5(0,0,4))*arr1(OPS_ACC6(0,0,4))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqW

!HEAT FLUX: Enthalpy fzrdfw
SUBROUTINE hf_kernel_eqX(rate,erhs,vtmp,utmp,store7,arr1,acbczr)
    
    implicit none

    real(kind=8), dimension(1) :: rate,erhs,vtmp
    real(kind=8), dimension(1), intent(in) :: utmp,store7,arr1
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8) :: fornow

    fornow = acbczr(4)*store7(OPS_ACC5(0,0,-4))*arr1(OPS_ACC6(0,0,-4)) &
           + acbczr(3)*store7(OPS_ACC5(0,0,-3))*arr1(OPS_ACC6(0,0,-3)) &
           + acbczr(2)*store7(OPS_ACC5(0,0,-2))*arr1(OPS_ACC6(0,0,-2)) &
           + acbczr(1)*store7(OPS_ACC5(0,0,-1))*arr1(OPS_ACC6(0,0,-1))

    rate(OPS_ACC1(0,0,0)) = rate(OPS_ACC1(0,0,0)) + fornow
    vtmp(OPS_ACC3(0,0,0)) = vtmp(OPS_ACC3(0,0,0)) + fornow
    erhs(OPS_ACC2(0,0,0)) = erhs(OPS_ACC2(0,0,0)) + fornow*utmp(OPS_ACC4(0,0,0))

END SUBROUTINE hf_kernel_eqX

!HEAT FLUX: Isothermal fxladw
SUBROUTINE hf_kernel_eqY(erhs,trun,tdrmix,store7,td1x,acbcxl,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1x(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(1,0,0))*tdrmix(OPS_ACC3(1,0,0))*store7(OPS_ACC4(1,0,0))*td1x(OPS_ACC5(1,0,0))
    combo2 = trun(OPS_ACC2(2,0,0))*tdrmix(OPS_ACC3(2,0,0))*store7(OPS_ACC4(2,0,0))*td1x(OPS_ACC5(2,0,0))
    combo3 = trun(OPS_ACC2(3,0,0))*tdrmix(OPS_ACC3(3,0,0))*store7(OPS_ACC4(3,0,0))*td1x(OPS_ACC5(3,0,0))
    combo4 = trun(OPS_ACC2(4,0,0))*tdrmix(OPS_ACC3(4,0,0))*store7(OPS_ACC4(4,0,0))*td1x(OPS_ACC5(4,0,0))

    fornow = acbcxl(1)*rgspec(ispec)*(combo1-combo0) &
           + acbcxl(2)*rgspec(ispec)*(combo2-combo0) &
           + acbcxl(3)*rgspec(ispec)*(combo3-combo0) &
           + acbcxl(4)*rgspec(ispec)*(combo4-combo0)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqY

!HEAT FLUX: Isothermal fxradw
SUBROUTINE hf_kernel_eqZ(erhs,trun,tdrmix,store7,td1x,acbcxr,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1x
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcxr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1x(OPS_ACC5(0,0,0))

    combo4 = trun(OPS_ACC2(-4,0,0))*tdrmix(OPS_ACC3(-4,0,0))*store7(OPS_ACC4(-4,0,0))*td1x(OPS_ACC5(-4,0,0))
    combo3 = trun(OPS_ACC2(-3,0,0))*tdrmix(OPS_ACC3(-3,0,0))*store7(OPS_ACC4(-3,0,0))*td1x(OPS_ACC5(-3,0,0))
    combo2 = trun(OPS_ACC2(-2,0,0))*tdrmix(OPS_ACC3(-2,0,0))*store7(OPS_ACC4(-2,0,0))*td1x(OPS_ACC5(-2,0,0))
    combo1 = trun(OPS_ACC2(-1,0,0))*tdrmix(OPS_ACC3(-1,0,0))*store7(OPS_ACC4(-1,0,0))*td1x(OPS_ACC5(-1,0,0))

    fornow = acbcxr(4)*rgspec(ispec)*(combo0-combo4) &
           + acbcxr(3)*rgspec(ispec)*(combo0-combo3) &
           + acbcxr(2)*rgspec(ispec)*(combo0-combo2) &
           + acbcxr(1)*rgspec(ispec)*(combo0-combo1)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqZ

!HEAT FLUX: Isothermal fyladw
SUBROUTINE hf_kernel_eqAA(erhs,trun,tdrmix,store7,td1y,acbcyl,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1y(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(0,1,0))*tdrmix(OPS_ACC3(0,1,0))*store7(OPS_ACC4(0,1,0))*td1y(OPS_ACC5(0,1,0))
    combo2 = trun(OPS_ACC2(0,2,0))*tdrmix(OPS_ACC3(0,2,0))*store7(OPS_ACC4(0,2,0))*td1y(OPS_ACC5(0,2,0))
    combo3 = trun(OPS_ACC2(0,3,0))*tdrmix(OPS_ACC3(0,3,0))*store7(OPS_ACC4(0,3,0))*td1y(OPS_ACC5(0,3,0))
    combo4 = trun(OPS_ACC2(0,4,0))*tdrmix(OPS_ACC3(0,4,0))*store7(OPS_ACC4(0,4,0))*td1y(OPS_ACC5(0,4,0))

    fornow = acbcyl(1)*rgspec(ispec)*(combo1-combo0) &
           + acbcyl(2)*rgspec(ispec)*(combo2-combo0) &
           + acbcyl(3)*rgspec(ispec)*(combo3-combo0) &
           + acbcyl(4)*rgspec(ispec)*(combo4-combo0)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqAA

!HEAT FLUX: Isothermal fyradw
SUBROUTINE hf_kernel_eqAB(erhs,trun,tdrmix,store7,td1y,acbcyr,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1y
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbcyr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1y(OPS_ACC5(0,0,0))

    combo4 = trun(OPS_ACC2(0,-4,0))*tdrmix(OPS_ACC3(0,-4,0))*store7(OPS_ACC4(0,-4,0))*td1y(OPS_ACC5(0,-4,0))
    combo3 = trun(OPS_ACC2(0,-3,0))*tdrmix(OPS_ACC3(0,-3,0))*store7(OPS_ACC4(0,-3,0))*td1y(OPS_ACC5(0,-3,0))
    combo2 = trun(OPS_ACC2(0,-2,0))*tdrmix(OPS_ACC3(0,-2,0))*store7(OPS_ACC4(0,-2,0))*td1y(OPS_ACC5(0,-2,0))
    combo1 = trun(OPS_ACC2(0,-1,0))*tdrmix(OPS_ACC3(0,-1,0))*store7(OPS_ACC4(0,-1,0))*td1y(OPS_ACC5(0,-1,0))

    fornow = acbcyr(4)*rgspec(ispec)*(combo0-combo4) &
           + acbcyr(3)*rgspec(ispec)*(combo0-combo3) &
           + acbcyr(2)*rgspec(ispec)*(combo0-combo2) &
           + acbcyr(1)*rgspec(ispec)*(combo0-combo1)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqAB

!HEAT FLUX: Isothermal fzladw
SUBROUTINE hf_kernel_eqAC(erhs,trun,tdrmix,store7,td1z,acbczl,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczl
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1z(OPS_ACC5(0,0,0))

    combo1 = trun(OPS_ACC2(0,0,1))*tdrmix(OPS_ACC3(0,0,1))*store7(OPS_ACC4(0,0,1))*td1z(OPS_ACC5(0,0,1))
    combo2 = trun(OPS_ACC2(0,0,2))*tdrmix(OPS_ACC3(0,0,2))*store7(OPS_ACC4(0,0,2))*td1z(OPS_ACC5(0,0,2))
    combo3 = trun(OPS_ACC2(0,0,3))*tdrmix(OPS_ACC3(0,0,3))*store7(OPS_ACC4(0,0,3))*td1z(OPS_ACC5(0,0,3))
    combo4 = trun(OPS_ACC2(0,0,4))*tdrmix(OPS_ACC3(0,0,4))*store7(OPS_ACC4(0,0,4))*td1z(OPS_ACC5(0,0,4))

    fornow = acbczl(1)*rgspec(ispec)*(combo1-combo0) &
           + acbczl(2)*rgspec(ispec)*(combo2-combo0) &
           + acbczl(3)*rgspec(ispec)*(combo3-combo0) &
           + acbczl(4)*rgspec(ispec)*(combo4-combo0)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqAC

!HEAT FLUX: Isothermal fzradw
SUBROUTINE hf_kernel_eqAD(erhs,trun,tdrmix,store7,td1z,acbczr,rgspec,ispec)
    
    implicit none

    real(kind=8), dimension(1) :: erhs
    real(kind=8), dimension(1), intent(in) :: trun,tdrmix,store7,td1z
    integer(kind=4), intent(in) :: ispec
    real(kind=8), dimension(ncbcsz), intent(in) :: acbczr
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    real(kind=8) :: fornow, combo0, combo1, combo2, combo3, combo4

    combo0 = trun(OPS_ACC2(0,0,0))*tdrmix(OPS_ACC3(0,0,0))*store7(OPS_ACC4(0,0,0))*td1z(OPS_ACC5(0,0,0))

    combo4 = trun(OPS_ACC2(0,0,-4))*tdrmix(OPS_ACC3(0,0,-4))*store7(OPS_ACC4(0,0,-4))*td1z(OPS_ACC5(0,0,-4))
    combo3 = trun(OPS_ACC2(0,0,-3))*tdrmix(OPS_ACC3(0,0,-3))*store7(OPS_ACC4(0,0,-3))*td1z(OPS_ACC5(0,0,-3))
    combo2 = trun(OPS_ACC2(0,0,-2))*tdrmix(OPS_ACC3(0,0,-2))*store7(OPS_ACC4(0,0,-2))*td1z(OPS_ACC5(0,0,-2))
    combo1 = trun(OPS_ACC2(0,0,-1))*tdrmix(OPS_ACC3(0,0,-1))*store7(OPS_ACC4(0,0,-1))*td1z(OPS_ACC5(0,0,-1))

    fornow = acbczr(4)*rgspec(ispec)*(combo0-combo4) &
           + acbczr(3)*rgspec(ispec)*(combo0-combo3) &
           + acbczr(2)*rgspec(ispec)*(combo0-combo2) &
           + acbczr(1)*rgspec(ispec)*(combo0-combo1)

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) + rgspec(ispec)*fornow

END SUBROUTINE hf_kernel_eqAD
