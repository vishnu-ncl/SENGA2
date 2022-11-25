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
