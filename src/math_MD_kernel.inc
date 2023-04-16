!A = B_multidim
SUBROUTINE math_MD_kernel_eqA(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqA

!A_multidim = B_multidim
SUBROUTINE math_MD_kernel_eqB(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqB

!A = B*C_multidim
SUBROUTINE math_MD_kernel_eqC(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqC

!A = B_multidim*C
SUBROUTINE math_MD_kernel_eqD(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_MD_kernel_eqD

!A_multidim = A_multidim/B
SUBROUTINE math_MD_kernel_eqE(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) / in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_MD_kernel_eqE

!A = B*C_multidim*D
SUBROUTINE math_MD_kernel_eqF(out_arr,in_arr1,in_arr2,in_arr3,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))*in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_MD_kernel_eqF

!A = B*C_multidim*(one-val1/D)
SUBROUTINE math_MD_kernel_eqG(out_arr,in_arr1,in_arr2,in_arr3,val1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer, intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: val1

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))*(1.0_8 - val1(ispec)/in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_MD_kernel_eqG

!A_multidim = A_multidim - half*(B+C+D)
SUBROUTINE math_MD_kernel_eqH(out_arr,in_arr1,in_arr2,in_arr3,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_8*(in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_MD_kernel_eqH

!A_multidim = A_multidim - half*B_multidim*C
SUBROUTINE math_MD_kernel_eqI(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_8*in_arr1(OPS_ACC_MD2(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_MD_kernel_eqI

!A_multidim = A_multidim - half*(B*C+D*E+F*G)
SUBROUTINE math_MD_kernel_eqJ(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_8*( in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                                                   + in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                                                   + in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) )

END SUBROUTINE math_MD_kernel_eqJ

!A_multidim = B_multidim - A_multidim*C - D*E - F*G - H*I
SUBROUTINE math_MD_kernel_eqK(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,in_arr7,in_arr8,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0)) - out_arr(OPS_ACC_MD1(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                                       - in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                                       - in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) &
                                                                       - in_arr7(OPS_ACC8(0,0,0))*in_arr8(OPS_ACC9(0,0,0))

END SUBROUTINE math_MD_kernel_eqK

SUBROUTINE math_MD_kernel_eqL(out_arr1,out_arr2,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1,out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6
    integer, intent(in) :: ispec
    real(kind=8) :: fornow

    fornow = in_arr1(OPS_ACC3(0,0,0))*in_arr2(OPS_ACC4(0,0,0))  &
            + in_arr3(OPS_ACC5(0,0,0))*in_arr4(OPS_ACC6(0,0,0))  &
            + in_arr5(OPS_ACC7(0,0,0))*in_arr6(OPS_ACC8(0,0,0))

!   Y EQUATION
    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = out_arr1(OPS_ACC_MD1(ispec,0,0,0)) + fornow

!   DIFFUSION CORRECTION VELOCITY DIVERGENCE
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow

END SUBROUTINE math_MD_kernel_eqL

SUBROUTINE math_MD_kernel_eqM(out_arr1,out_arr2,in_arr1,in_arr2,in_arr3,in_arr4,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1,out_arr2
    real(kind=8), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4
    integer, intent(in) :: ispec
    real(kind=8) :: fornow

    fornow = ( in_arr1(OPS_ACC3(0,0,0)) + in_arr2(OPS_ACC4(0,0,0)) + in_arr3(OPS_ACC5(0,0,0)) ) &
            * in_arr4(OPS_ACC6(0,0,0))

!   Y EQUATION
    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = out_arr1(OPS_ACC_MD1(ispec,0,0,0)) + fornow

!   DIFFUSION CORRECTION VELOCITY DIVERGENCE
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow

END SUBROUTINE math_MD_kernel_eqM

!A = A + var*B_multidim
SUBROUTINE math_MD_kernel_eqP(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: val1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + val1*in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqP

!A = A + var*B_multidim
SUBROUTINE math_MD_kernel_eqPP(out_arr1,in_arr1,rgspec,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    real(kind=8), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec

    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + rgspec(ispec)*in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqPP

!A = A*max(B_multidim*var    zero)
SUBROUTINE math_MD_kernel_eqR(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    real(kind=8), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: val1
    real(kind=8) :: fornow

    fornow = max(in_arr1(OPS_ACC_MD2(ispec,0,0,0))*val1(ispec), 0.0_8)
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) * fornow

END SUBROUTINE math_MD_kernel_eqR
    
!A_multidim = A_multidim + B*var
SUBROUTINE math_MD_kernel_eqS(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: val1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = out_arr1(OPS_ACC_MD1(ispec,0,0,0)) + in_arr1(OPS_ACC2(0,0,0))*val1

END SUBROUTINE math_MD_kernel_eqS

!A_multidim = var
SUBROUTINE math_MD_kernel_eqT(out_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    integer, intent(in) :: ispec
    real(kind=8), dimension(nspcmx), intent(in) :: val1

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = val1(ispec)

END SUBROUTINE math_MD_kernel_eqT

!A_multidim = B*A_multidim
SUBROUTINE math_MD_kernel_eqU(out_arr1,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    real(kind=8), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*out_arr1(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqU

SUBROUTINE math_MD_kernel_eqV1(transp,trun,tdifgb)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: trun
    real(kind=8), dimension(1) :: transp

    real(kind=8), intent(in) :: tdifgb

    transp(OPS_ACC1(0,0,0)) = LOG(trun(OPS_ACC2(0,0,0))/tdifgb)

END SUBROUTINE math_MD_kernel_eqV1

SUBROUTINE math_MD_kernel_eqV2(combo1,combo2,combo3,transp,yrhs,condco,ovwmol,ncocon,ncocm1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: transp,yrhs
    real(kind=8), dimension(1) :: combo1,combo2,combo3
    integer, intent(in) :: ncocon,ncocm1,ispec

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), dimension(nccfmx,nspcmx), intent(in) :: condco

    real(kind=8) :: fornow,ctrans
    integer :: icp

!   CONDUCTIVITY FOR EACH SPECIES
    fornow = condco(ncocon,ispec)
    DO icp = ncocm1,1,-1
        fornow = fornow*transp(OPS_ACC4(0,0,0)) + condco(icp,ispec)
    END DO
    ctrans = EXP(fornow)

!   COMBINATION RULE FOR CONDUCTIVITY
    fornow = yrhs(OPS_ACC_MD5(ispec,0,0,0))*ovwmol(ispec)
    combo1(OPS_ACC1(0,0,0)) = combo1(OPS_ACC1(0,0,0)) + fornow*ctrans
    combo2(OPS_ACC2(0,0,0)) = combo2(OPS_ACC2(0,0,0)) + fornow/ctrans
    combo3(OPS_ACC3(0,0,0)) = combo3(OPS_ACC3(0,0,0)) + fornow

END SUBROUTINE math_MD_kernel_eqV2

SUBROUTINE math_MD_kernel_eqV3(combo1,combo2,combo3,store7,wmomix,drhs)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: drhs
    real(kind=8), dimension(1) :: combo1,combo2,combo3,store7,wmomix

    combo3(OPS_ACC3(0,0,0)) = 1.0_8/combo3(OPS_ACC3(0,0,0))
    combo1(OPS_ACC1(0,0,0)) = combo1(OPS_ACC1(0,0,0))*combo3(OPS_ACC3(0,0,0))
    combo2(OPS_ACC2(0,0,0)) = combo2(OPS_ACC2(0,0,0))*combo3(OPS_ACC3(0,0,0))

    store7(OPS_ACC4(0,0,0)) = 0.5_8*(combo1(OPS_ACC1(0,0,0)) + 1.0_8/combo2(OPS_ACC2(0,0,0)))
    wmomix(OPS_ACC5(0,0,0)) = drhs(OPS_ACC6(0,0,0))*combo3(OPS_ACC3(0,0,0))

END SUBROUTINE math_MD_kernel_eqV3

SUBROUTINE math_MD_kernel_eqW(difmix,store7,transp,prun,yrhs,wmomix,drhs,diffco,ovwmol,pdifgb,dfctol,ispec,ncodif,ncodm1,nspec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: difmix,store7
    real(kind=8), dimension(1), intent(in) :: transp,prun,yrhs,wmomix,drhs

    integer, intent(in) :: ispec,ncodif,ncodm1,nspec

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), dimension(ndcfmx,nspcmx,nspcmx), intent(in) :: diffco
    real(kind=8), intent(in) :: pdifgb,dfctol

    real(kind=8) :: fornow, combo1, combo2
    real(kind=8), dimension(nspcmx) :: ctrans
    integer :: jspec,icp

    DO jspec = 1, nspec
        fornow = diffco(ncodif,jspec,ispec)
        DO icp = ncodm1,1,-1
            fornow = fornow*transp(OPS_ACC3(0,0,0)) + diffco(icp,jspec,ispec)
        END DO
        ctrans(jspec) = EXP(fornow)*pdifgb/prun(OPS_ACC4(0,0,0))
    END DO

!   COMBINATION RULE FOR MASS DIFFUSIVITY
    combo1 = 0.0_8
    combo2 = 0.0_8
    DO jspec = 1, nspec
        fornow = yrhs(OPS_ACC_MD5(jspec,0,0,0)) + dfctol
        combo1 = combo1 + fornow
        combo2 = combo2 + fornow*ovwmol(jspec)/ctrans(jspec)
    END DO

    fornow = yrhs(OPS_ACC_MD5(ispec,0,0,0)) + dfctol
    combo1 = combo1 - fornow
    combo2 = combo2 - fornow*ovwmol(ispec)/ctrans(ispec)
    combo2 = combo2*wmomix(OPS_ACC6(0,0,0))

    difmix(OPS_ACC1(0,0,0)) = drhs(OPS_ACC7(0,0,0))*combo1/combo2
    store7(OPS_ACC2(0,0,0)) = difmix(OPS_ACC1(0,0,0))

END SUBROUTINE math_MD_kernel_eqW

SUBROUTINE math_MD_kernel_eqX(tdrmix,trun,yrhs,wmomix,tdrcco,ovwmol,tdifgb,ncotdr,jspec,ispec,ncotm1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: tdrmix
    real(kind=8), dimension(1), intent(in) :: trun,yrhs,wmomix

    integer, intent(in) :: ncotdr,jspec,ispec,ncotm1

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), dimension(ndcfmx,nspcmx,nspcmx), intent(in) :: tdrcco
    real(kind=8), intent(in) :: tdifgb

    real(kind=8) :: fornow, combo2
    real(kind=8), dimension(nspcmx) :: ctrans
    integer :: icp

!   THERMAL DIFFUSION RATIO FOR THIS SPECIES PAIR
    combo2 = trun(OPS_ACC2(0,0,0))/tdifgb
    fornow = tdrcco(ncotdr,jspec,ispec)
    DO icp = ncotm1,1,-1
        fornow = fornow*combo2 + tdrcco(icp,jspec,ispec)
    END DO
    ctrans(jspec) = fornow

!   COMBINATION RULE FOR THERMAL DIFFUSIION RATIO
    fornow = yrhs(OPS_ACC_MD3(jspec,0,0,0))*ovwmol(jspec)
    tdrmix(OPS_ACC1(0,0,0)) = tdrmix(OPS_ACC1(0,0,0))  &
                     + fornow*wmomix(OPS_ACC4(0,0,0))*ctrans(jspec)

END SUBROUTINE math_MD_kernel_eqX

SUBROUTINE math_MD_kernel_eqY1(ctrans,transp,viscco,ncovis,ncovm1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: transp
    real(kind=8), dimension(1) :: ctrans
    integer, intent(in) :: ncovis,ncovm1,ispec

    real(kind=8), dimension(nvcfmx,nspcmx), intent(in) :: viscco

    real(kind=8) :: fornow
    integer :: icp

    fornow = viscco(ncovis,ispec)
    DO icp = ncovm1,1,-1
        fornow = fornow*transp(OPS_ACC2(0,0,0)) + viscco(icp,ispec)
    END DO

    ctrans(OPS_ACC_MD1(ispec,0,0,0)) = EXP(fornow)

END SUBROUTINE math_MD_kernel_eqY1

SUBROUTINE math_MD_kernel_eqY2(combo2,ctrans,yrhs,ovwmol,wilko1,wilko2,ispec,jspec)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: ctrans,yrhs
    real(kind=8), dimension(1) :: combo2
    integer, intent(in) :: ispec,jspec

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8), dimension(nspcmx,nspcmx), intent(in) :: wilko1,wilko2

    real(kind=8) :: fornow

    fornow = SQRT(ctrans(OPS_ACC_MD2(ispec,0,0,0))/ctrans(OPS_ACC_MD2(jspec,0,0,0)))
    fornow = 1.0_8 + fornow*wilko2(jspec,ispec)
    fornow = wilko1(jspec,ispec)*fornow*fornow
    combo2(OPS_ACC1(0,0,0)) = combo2(OPS_ACC1(0,0,0)) + yrhs(OPS_ACC_MD3(jspec,0,0,0))*ovwmol(jspec)*fornow

END SUBROUTINE math_MD_kernel_eqY2

SUBROUTINE math_MD_kernel_eqY3(combo1,ctrans,combo2,yrhs,ovwmol,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: ctrans,combo2,yrhs
    real(kind=8), dimension(1) :: combo1
    integer, intent(in) :: ispec

    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol

    real(kind=8) :: fornow

    fornow = ctrans(OPS_ACC_MD2(ispec,0,0,0))/combo2(OPS_ACC3(0,0,0))
    combo1(OPS_ACC1(0,0,0)) = combo1(OPS_ACC1(0,0,0)) + yrhs(OPS_ACC_MD4(ispec,0,0,0))*ovwmol(ispec)*fornow

END SUBROUTINE math_MD_kernel_eqY3

SUBROUTINE math_MD_kernel_eqZ(utmp,wtmp,trun,itndex,yrhs,amasch,ncpoly,ncpom1,ncenth,icoef1,icoef2,iindex,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: utmp,wtmp
    real(kind=8), dimension(1), intent(in) :: trun,yrhs
    integer, dimension(1), intent(in) :: itndex

    integer, intent(in) ::  iindex,icoef1,icoef2,ispec

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    integer, dimension(ntinmx,nspcmx), intent(in) :: ncpoly, ncpom1, ncenth

    real(kind=8) :: fornow
    integer :: itint,icp

    itint = 1 + MOD(itndex(OPS_ACC_MD4(iindex,0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)
    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*trun(OPS_ACC3(0,0,0)) + amasch(icp,itint,ispec)
    END DO
    utmp(OPS_ACC1(0,0,0)) = amasch(ncenth(itint,ispec),itint,ispec)  &
                    + fornow*trun(OPS_ACC3(0,0,0))

!   MIXTURE H
    wtmp(OPS_ACC2(0,0,0)) = wtmp(OPS_ACC2(0,0,0)) + utmp(OPS_ACC1(0,0,0))*yrhs(OPS_ACC_MD5(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqZ

SUBROUTINE math_MD_kernel_eqAA(store1,yrhs,ovwmol,scoef,ysmall,ydenom,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: store1
    real(kind=8), dimension(1), intent(in) :: yrhs
    integer, intent(in) ::  ispec
    real(kind=8), intent(in) :: scoef,ysmall,ydenom
    real(kind=8), dimension(nspcmx), intent(in) :: ovwmol
    real(kind=8) :: fornow

    fornow = yrhs(OPS_ACC_MD2(ispec,0,0,0))*ovwmol(ispec)
    fornow = MAX(fornow,ysmall)
    fornow = EXP(scoef*LOG(fornow))
    fornow = fornow/(1.0_8+ydenom*fornow)
    store1(OPS_ACC1(0,0,0)) = store1(OPS_ACC1(0,0,0))*fornow

END SUBROUTINE math_MD_kernel_eqAA

SUBROUTINE math_MD_kernel_eqAB(store2,trun,itndex,amolgb,ncpoly,ncpom1,ncenth,ncenpy,diffmu,isspec,ispec,istep)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: store2
    real(kind=8), dimension(1), intent(in) :: trun
    integer, dimension(1), intent(in) :: itndex
    integer, intent(in) ::  isspec,ispec,istep

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amolgb
    real(kind=8), dimension(nssmax,nstpmx), intent(in) :: diffmu
    integer, dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth,ncenpy

    integer :: iindex,ipower,icoef2,icoef1,itint,icp
    real(kind=8) :: gibbsp

!   LOCATE THE TEMPERATURE IN AN INTERVAL
    iindex = 1 + (ispec-1)/nspimx
    ipower = ispec - (iindex-1)*nspimx - 1
    icoef2 = ntbase**ipower
    icoef1 = icoef2*ntbase
    itint = 1 + MOD(itndex(OPS_ACC_MD3(iindex,0,0,0)),icoef1)/icoef2

!   CONSTRUCT GIBBS FUNCTION FROM ITS POLYNOMIAL COEFFICIENTS
    gibbsp = amolgb(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        gibbsp = amolgb(icp,itint,ispec) + gibbsp*trun(OPS_ACC2(0,0,0))
    END DO

    gibbsp = amolgb(ncenth(itint,ispec),itint,ispec)  &
            /trun(OPS_ACC2(0,0,0)) - amolgb(ncenpy(itint,ispec),itint,ispec)  &
            *LOG(trun(OPS_ACC2(0,0,0))) - gibbsp

!   ADD GIBBS FUNCTION CONTRIBUTION TO RATE COEFFICIENT
!   USING STEP SPECIES DELTA-LIST
!   TO GET BACKWARD RATE COEFFICIENT
    store2(OPS_ACC1(0,0,0)) = store2(OPS_ACC1(0,0,0)) + diffmu(isspec,istep)*gibbsp

END SUBROUTINE math_MD_kernel_eqAB

SUBROUTINE math_MD_kernel_eqAC(itndex,trun,tinthi,ntint,nspec)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: trun
    integer, dimension(1) :: itndex

    integer, intent(in) ::  nspec

    real(kind=8), dimension(ntinmx,nspcmx), intent(in) :: tinthi
    integer, dimension(nspcmx), intent(in) :: ntint

    integer :: iindex,ispec,itint,ipower

    DO iindex = 1,nintmx
        itndex(OPS_ACC_MD1(iindex,0,0,0)) = 0
    END DO

    DO ispec = 1, nspec

        itint = 1

        DO WHILE (trun(OPS_ACC2(0,0,0)) > tinthi(itint,ispec) .and. itint < ntint(ispec))
            itint = itint + 1
        END DO

!       SET THE TEMPERATURE INDEX
        iindex = 1 + (ispec-1)/nspimx
        ipower = ispec - (iindex-1)*nspimx - 1

        itndex(OPS_ACC_MD1(iindex,0,0,0)) = itndex(OPS_ACC_MD1(iindex,0,0,0))  &
            +(itint-1)*ntbase**ipower

    END DO

END SUBROUTINE math_MD_kernel_eqAC

SUBROUTINE math_MD_kernel_eqAD(erun,trun,itndex,yrun,amasch,ncpoly,ncpom1,ncenth,ispec,iindex,ipower,icoef1,icoef2)
    use data_types
    implicit none

    real(kind=8), dimension(1), intent(in) :: trun, yrun
    real(kind=8), dimension(1) :: erun
    integer, dimension(1), intent(in) :: itndex
    
    integer, intent(in) ::  ispec,iindex,ipower,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    integer, dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer :: icp,itint

    itint = 1 + MOD(itndex(OPS_ACC_MD3(iindex,0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*trun(OPS_ACC2(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
           + fornow*trun(OPS_ACC2(0,0,0))

    erun(OPS_ACC1(0,0,0)) = erun(OPS_ACC1(0,0,0)) + fornow*yrun(OPS_ACC_MD4(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqAD
