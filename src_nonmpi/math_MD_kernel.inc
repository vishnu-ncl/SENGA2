!A = B_multidim
SUBROUTINE math_MD_kernel_eqA(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqA

!A_multidim = B_multidim
SUBROUTINE math_MD_kernel_eqB(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqB

!A = B*C_multidim
SUBROUTINE math_MD_kernel_eqC(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqC

!A = B_multidim*C
SUBROUTINE math_MD_kernel_eqD(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_MD_kernel_eqD

!A_multidim = A_multidim/B
SUBROUTINE math_MD_kernel_eqE(out_arr,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) / in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE math_MD_kernel_eqE

!A = B*C_multidim*D
SUBROUTINE math_MD_kernel_eqF(out_arr,in_arr1,in_arr2,in_arr3,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))*in_arr3(OPS_ACC4(0,0,0))

END SUBROUTINE math_MD_kernel_eqF

!A = B*C_multidim*(one-val1/D)
SUBROUTINE math_MD_kernel_eqG(out_arr,in_arr1,in_arr2,in_arr3,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC_MD3(ispec,0,0,0))*(1.0_dp - val1/in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_MD_kernel_eqG

!A_multidim = A_multidim - half*(B+C+D)
SUBROUTINE math_MD_kernel_eqH(out_arr,in_arr1,in_arr2,in_arr3,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_dp*(in_arr1(OPS_ACC2(0,0,0))+in_arr2(OPS_ACC3(0,0,0))+in_arr3(OPS_ACC4(0,0,0)))

END SUBROUTINE math_MD_kernel_eqH

!A_multidim = A_multidim - half*B_multidim*C
SUBROUTINE math_MD_kernel_eqI(out_arr,in_arr1,in_arr2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_dp*in_arr1(OPS_ACC_MD2(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0))

END SUBROUTINE math_MD_kernel_eqI

!A_multidim = A_multidim - half*(B*C+D*E+F*G)
SUBROUTINE math_MD_kernel_eqJ(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = out_arr(OPS_ACC_MD1(ispec,0,0,0)) - 0.5_dp*( in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                                                   + in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                                                   + in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) )

END SUBROUTINE math_MD_kernel_eqJ

!A_multidim = B_multidim - A_multidim*C - D*E - F*G - H*I
SUBROUTINE math_MD_kernel_eqK(out_arr,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,in_arr7,in_arr8,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr
    real(kind=dp), dimension(1), intent(in) :: in_arr1, in_arr2, in_arr3, in_arr4, in_arr5, in_arr6, in_arr7, in_arr8
    integer, intent(in) :: ispec

    out_arr(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC_MD2(ispec,0,0,0)) - out_arr(OPS_ACC_MD1(ispec,0,0,0))*in_arr2(OPS_ACC3(0,0,0)) &
                                                                       - in_arr3(OPS_ACC4(0,0,0))*in_arr4(OPS_ACC5(0,0,0)) &
                                                                       - in_arr5(OPS_ACC6(0,0,0))*in_arr6(OPS_ACC7(0,0,0)) &
                                                                       - in_arr7(OPS_ACC8(0,0,0))*in_arr8(OPS_ACC9(0,0,0))

END SUBROUTINE math_MD_kernel_eqK

SUBROUTINE math_MD_kernel_eqL(out_arr1,out_arr2,in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1,out_arr2
    real(kind=dp), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4,in_arr5,in_arr6
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

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

    real(kind=dp), dimension(1) :: out_arr1,out_arr2
    real(kind=dp), dimension(1), intent(in) :: in_arr1,in_arr2,in_arr3,in_arr4
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = ( in_arr1(OPS_ACC3(0,0,0)) + in_arr2(OPS_ACC4(0,0,0)) + in_arr3(OPS_ACC5(0,0,0)) ) &
            * in_arr4(OPS_ACC6(0,0,0))

!   Y EQUATION
    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = out_arr1(OPS_ACC_MD1(ispec,0,0,0)) + fornow

!   DIFFUSION CORRECTION VELOCITY DIVERGENCE
    out_arr2(OPS_ACC2(0,0,0)) = out_arr2(OPS_ACC2(0,0,0)) + fornow

END SUBROUTINE math_MD_kernel_eqM

!A_multidim_lastIdx = A_multidim_lastIdx + A_multidim_loopIdx
SUBROUTINE math_MD_kernel_eqN(out_arr1,ispec,nspec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    integer, intent(in) :: ispec,nspec

    out_arr1(OPS_ACC_MD1(nspec,0,0,0)) = out_arr1(OPS_ACC_MD1(nspec,0,0,0)) + out_arr1(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqN

!A_multidim = 1.0 - A_multidim
SUBROUTINE math_MD_kernel_eqO(out_arr1,nspec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    integer, intent(in) :: nspec

    out_arr1(OPS_ACC_MD1(nspec,0,0,0)) = 1.0_dp - out_arr1(OPS_ACC_MD1(nspec,0,0,0))

END SUBROUTINE math_MD_kernel_eqO

!A = A + var*B_multidim
SUBROUTINE math_MD_kernel_eqP(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec
    
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + val1*in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqP

!A_multidim = var1 + B*(var2-var1)
SUBROUTINE math_MD_kernel_eqQ(out_arr1,in_arr1,val1,val2,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: val1,val2
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = val1 + in_arr1(OPS_ACC2(0,0,0))*(val2-val1)
    
END SUBROUTINE math_MD_kernel_eqQ

!A = A*max(B_multidim*var    zero)
SUBROUTINE math_MD_kernel_eqR(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow

    fornow = max(in_arr1(OPS_ACC_MD2(ispec,0,0,0))*val1, 0.0_dp)
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) * fornow

END SUBROUTINE math_MD_kernel_eqR
    
!A_multidim = A_multidim + B*var
SUBROUTINE math_MD_kernel_eqS(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = out_arr1(OPS_ACC_MD1(ispec,0,0,0)) + in_arr1(OPS_ACC2(0,0,0))*val1

END SUBROUTINE math_MD_kernel_eqS

!A_multidim = var
SUBROUTINE math_MD_kernel_eqT(out_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), intent(in) :: val1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = val1

END SUBROUTINE math_MD_kernel_eqT

!A_multidim = B*A_multidim
SUBROUTINE math_MD_kernel_eqU(out_arr1,in_arr1,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: out_arr1
    real(kind=dp), dimension(1), intent(in) :: in_arr1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = in_arr1(OPS_ACC2(0,0,0))*out_arr1(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE math_MD_kernel_eqU

SUBROUTINE math_MD_kernel_eqV(transp,store7,wmomix,trun,yrhs,drhs,condco,ovwmol,tdifgb,ncocon,ncocm1,nspec,nccfmx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: transp,store7,wmomix
    real(kind=dp), dimension(1), intent(in) :: trun,yrhs,drhs
    
    integer, intent(in) :: ncocon,ncocm1,nspec,nccfmx
    
    real(kind=dp), dimension(nspec), intent(in) :: ovwmol
    real(kind=dp), dimension(nccfmx,nspec), intent(in) :: condco
    
    real(kind=dp), intent(in) :: tdifgb

    real(kind=dp) :: fornow, combo1, combo2, combo3
    real(kind=dp), dimension(nspec) :: ctrans
    integer :: ispec,icp

!   CONDUCTIVITY FOR EACH SPECIES
    transp(OPS_ACC1(0,0,0)) = LOG(trun(OPS_ACC4(0,0,0))/tdifgb)
    DO ispec = 1, nspec
        fornow = condco(ncocon,ispec)
        DO icp = ncocm1,1,-1
            fornow = fornow*transp(OPS_ACC1(0,0,0)) + condco(icp,ispec)
        END DO
        ctrans(ispec) = EXP(fornow)
    END DO

!   COMBINATION RULE FOR CONDUCTIVITY
    combo1 = 0.0_dp
    combo2 = 0.0_dp
    combo3 = 0.0_dp
    DO ispec = 1, nspec
        fornow = yrhs(OPS_ACC_MD5(ispec,0,0,0))*ovwmol(ispec)
        combo1 = combo1 + fornow*ctrans(ispec)
        combo2 = combo2 + fornow/ctrans(ispec)
        combo3 = combo3 + fornow
    END DO

    combo3 = 1.0_dp/combo3
    combo1 = combo1*combo3
    combo2 = combo2*combo3
    store7(OPS_ACC2(0,0,0)) = 0.5_dp*(combo1 + 1.0_dp/combo2)
    wmomix(OPS_ACC3(0,0,0)) = drhs(OPS_ACC6(0,0,0))*combo3

END SUBROUTINE math_MD_kernel_eqV
    
SUBROUTINE math_MD_kernel_eqW(difmix,store7,transp,prun,yrhs,wmomix,drhs,diffco,ovwmol,pdifgb,dfctol,ispec,ncodif,ncodm1,nspec,ndcfmx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: difmix,store7
    real(kind=dp), dimension(1), intent(in) :: transp,prun,yrhs,wmomix,drhs

    integer, intent(in) :: ispec,ncodif,ncodm1,nspec,ndcfmx

    real(kind=dp), dimension(nspec), intent(in) :: ovwmol
    real(kind=dp), dimension(ndcfmx,nspec,nspec), intent(in) :: diffco
    real(kind=dp), intent(in) :: pdifgb,dfctol

    real(kind=dp) :: fornow, combo1, combo2
    real(kind=dp), dimension(nspec) :: ctrans
    integer :: jspec,icp

    DO jspec = 1, nspec
        fornow = diffco(ncodif,jspec,ispec)
        DO icp = ncodm1,1,-1
            fornow = fornow*transp(OPS_ACC3(0,0,0)) + diffco(icp,jspec,ispec)
        END DO
        ctrans(jspec) = EXP(fornow)*pdifgb/prun(OPS_ACC4(0,0,0))
    END DO

!   COMBINATION RULE FOR MASS DIFFUSIVITY
    combo1 = 0.0_dp
    combo2 = 0.0_dp
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

SUBROUTINE math_MD_kernel_eqX(tdrmix,trun,yrhs,wmomix,tdrcco,ovwmol,tdifgb,ncotdr,jspec,ispec,ncotm1,nspec,ndcfmx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: tdrmix
    real(kind=dp), dimension(1), intent(in) :: trun,yrhs,wmomix

    integer, intent(in) :: ncotdr,jspec,ispec,ncotm1,nspec,ndcfmx

    real(kind=dp), dimension(nspec), intent(in) :: ovwmol
    real(kind=dp), dimension(ndcfmx,nspec,nspec), intent(in) :: tdrcco
    real(kind=dp), intent(in) :: tdifgb

    real(kind=dp) :: fornow, combo2
    real(kind=dp), dimension(nspec) :: ctrans
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

SUBROUTINE math_MD_kernel_eqY(difmix,transp,yrhs,viscco,wilko1,wilko2,ovwmol,ncovis,ncovm1,nspec,nvcfmx)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: difmix
    real(kind=dp), dimension(1), intent(in) :: transp,yrhs

    integer, intent(in) :: ncovis,ncovm1,nspec,nvcfmx

    real(kind=dp), dimension(nspec), intent(in) :: ovwmol
    real(kind=dp), dimension(nvcfmx,nspec), intent(in) :: viscco
    real(kind=dp), dimension(nspec,nspec), intent(in) :: wilko1,wilko2

    real(kind=dp) :: fornow, combo1, combo2
    real(kind=dp), dimension(nspec) :: ctrans
    integer :: ispec, jspec, icp

!   VISCOSITY FOR EACH SPECIES
    DO ispec = 1, nspec
        fornow = viscco(ncovis,ispec)
        DO icp = ncovm1,1,-1
            fornow = fornow*transp(OPS_ACC2(0,0,0)) + viscco(icp,ispec)
        END DO
        ctrans(ispec) = EXP(fornow)
    END DO

!   COMBINATION RULE FOR VISCOSITY
    combo1 = 0.0_dp
    DO ispec = 1, nspec
        combo2 = 0.0_dp
        DO jspec = 1, nspec
            fornow = SQRT(ctrans(ispec)/ctrans(jspec))
            fornow = 1.0_dp + fornow*wilko2(jspec,ispec)
            fornow = wilko1(jspec,ispec)*fornow*fornow
            combo2 = combo2 + yrhs(OPS_ACC_MD3(jspec,0,0,0))*ovwmol(jspec)*fornow
        END DO
        fornow = ctrans(ispec)/combo2
        combo1 = combo1 + yrhs(OPS_ACC_MD3(ispec,0,0,0))*ovwmol(ispec)*fornow
    END DO

    difmix(OPS_ACC1(0,0,0)) = combo1

END SUBROUTINE math_MD_kernel_eqY

SUBROUTINE math_MD_kernel_eqZ(utmp,wtmp,trun,itndex,yrhs,amasch,ncpoly,ncpom1,ncenth,icoef1,icoef2,iindex,ispec,ncofmx,ntinmx,nspec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: utmp, wtmp
    real(kind=dp), dimension(1), intent(in) :: trun, yrhs
    integer, dimension(1), intent(in) :: itndex

    integer, intent(in) ::  iindex, icoef1, icoef2, ispec, ncofmx, ntinmx, nspec

    real(kind=dp), dimension(ncofmx,ntinmx,nspec), intent(in) :: amasch
    integer, dimension(ntinmx,nspec), intent(in) :: ncpoly, ncpom1, ncenth

    real(kind=dp) :: fornow
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

    real(kind=dp), dimension(1) :: store1
    real(kind=dp), dimension(1), intent(in) :: yrhs
    integer, intent(in) ::  ispec
    real(kind=dp), intent(in) :: ovwmol,scoef,ysmall,ydenom

    real(kind=dp) :: fornow

    fornow = yrhs(OPS_ACC_MD2(ispec,0,0,0))*ovwmol
    fornow = MAX(fornow,ysmall)
    fornow = EXP(scoef*LOG(fornow))
    fornow = fornow/(1.0_dp+ydenom*fornow)
    store1(OPS_ACC1(0,0,0)) = store1(OPS_ACC1(0,0,0))*fornow

END SUBROUTINE math_MD_kernel_eqAA

SUBROUTINE math_MD_kernel_eqAB(store2,trun,itndex,amolgb,ncpoly,ncpom1,ncenth,ncenpy,diffmu,nspimx,ntbase,ncofmx,ntinmx,nspec,nssmax,nstpmx,isspec,ispec,istep)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: store2
    real(kind=dp), dimension(1), intent(in) :: trun
    integer, dimension(1), intent(in) :: itndex
    integer, intent(in) ::  nspimx,ntbase,ncofmx,ntinmx,nspec,nssmax,nstpmx,isspec,ispec,istep

    real(kind=dp), dimension(ncofmx,ntinmx,nspec), intent(in) :: amolgb
    real(kind=dp), dimension(nssmax,nstpmx), intent(in) :: diffmu
    integer, dimension(ntinmx,nspec), intent(in) :: ncpoly,ncpom1,ncenth,ncenpy

    integer :: iindex,ipower,icoef2,icoef1,itint,icp
    real(kind=dp) :: gibbsp

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

SUBROUTINE math_MD_kernel_eqAC(itndex,trun,tinthi,ntint,nspimx,ntbase,nintmx,nspec)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: trun
    integer, dimension(1) :: itndex

    integer, intent(in) ::  nspimx,ntbase,nintmx,nspec

    real(kind=dp), dimension(nintmx,nspec), intent(in) :: tinthi
    integer, dimension(nspec), intent(in) :: ntint

    integer :: iindex,ispec,itint,ipower

    DO iindex = 1,nintmx
        itndex(OPS_ACC_MD1(iindex,0,0,0)) = 0
    END DO

    DO ispec = 1,nspec

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

SUBROUTINE math_MD_kernel_eqAD(erun,trun,itndex,yrun,amasch,ncpoly,ncpom1,ncenth,ispec,iindex,ipower,icoef1,icoef2,ncofmx,ntinmx,nspec)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: trun, yrun
    real(kind=dp), dimension(1) :: erun
    integer, dimension(1), intent(in) :: itndex
    
    integer, intent(in) ::  ispec,iindex,ipower,icoef1,icoef2,ncofmx,ntinmx,nspec

    real(kind=dp), dimension(ncofmx,ntinmx,nspec), intent(in) :: amasch
    integer, dimension(ntinmx,nspec), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=dp) :: fornow
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
