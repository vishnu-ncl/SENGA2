SUBROUTINE flamin_kernel_set_reaction_var(crin,deltag,clocat,cthick,idx)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: crin
    real(kind=8), intent(in) :: deltag,clocat,cthick
    integer, dimension(3), intent(in) :: idx

    real(kind=8) :: xcoord,argmnt

    real(kind=8) :: erfunc
    integer, parameter :: ncoeff = 5
    integer, parameter :: ncofm1 = ncoeff-1

    real(kind=8) :: ecoeff(ncoeff)
    real(kind=8) :: pcoeff
    real(kind=8) :: etotal,zvalue,tvalue
    integer :: ix,icoeff

!   SET THE COEFFICIENTS
    pcoeff = 0.3275911_8
    ecoeff(1) = 0.254829592_8
    ecoeff(2) =-0.284496736_8
    ecoeff(3) = 1.421413741_8
    ecoeff(4) =-1.453152027_8
    ecoeff(5) = 1.061405429_8

    ix = idx(1)
    xcoord = real(ix-1,kind=8)*deltag
    argmnt = (xcoord-clocat)/cthick

!   EVALUATE ERROR FUNCTION
    zvalue = abs(argmnt)
    tvalue = 1.0_8/(1.0_8+pcoeff*zvalue)

    etotal = ecoeff(ncoeff)
    DO icoeff = ncofm1,1,-1
        etotal = ecoeff(icoeff) + etotal*tvalue
    END DO
    etotal = etotal*tvalue

    erfunc = 1.0_8 - etotal*exp(-zvalue*zvalue)
    IF(argmnt < 0.0_8) erfunc = -erfunc

    crin(OPS_ACC1(0,0,0)) = 0.5_8*(1.0_8+erfunc)

END SUBROUTINE flamin_kernel_set_reaction_var

!A_multidim = var1 + B*(var2-var1)
SUBROUTINE flamin_kernel_eqA(out_arr1,in_arr1,val1,val2,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: val1,val2
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC_MD1(ispec,0,0,0)) = val1 + in_arr1(OPS_ACC2(0,0,0))*(val2-val1)

END SUBROUTINE flamin_kernel_eqA

!A_multidim_lastIdx = A_multidim_lastIdx + A_multidim_loopIdx
SUBROUTINE flamin_kernel_eqB(out_arr1,ispec,nspec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    integer, intent(in) :: ispec,nspec

    out_arr1(OPS_ACC_MD1(nspec,0,0,0)) = out_arr1(OPS_ACC_MD1(nspec,0,0,0)) + out_arr1(OPS_ACC_MD1(ispec,0,0,0))

END SUBROUTINE flamin_kernel_eqB

!A_multidim = 1.0 - A_multidim
SUBROUTINE flamin_kernel_eqC(out_arr1,nspec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    integer, intent(in) :: nspec

    out_arr1(OPS_ACC_MD1(nspec,0,0,0)) = 1.0_8 - out_arr1(OPS_ACC_MD1(nspec,0,0,0))

END SUBROUTINE flamin_kernel_eqC

!A = var1 + B_xdimonly*(var2-var1)
SUBROUTINE flamin_kernel_eqD(out_arr,in_arr1,var1,var2)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: var1,var2

    out_arr(OPS_ACC1(0,0,0)) = var1 + in_arr1(OPS_ACC2(0,0,0))*(var2-var1)

END SUBROUTINE flamin_kernel_eqD

!A = A + var*B_multidim
SUBROUTINE flamin_kernel_eqE(out_arr1,in_arr1,val1,ispec)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr1
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: val1
    integer, intent(in) :: ispec

    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0)) + val1*in_arr1(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE flamin_kernel_eqE

!A = var/(B*C)
SUBROUTINE flamin_kernel_eqF(out_arr, in_arr1, in_arr2, var1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=8), intent(in) :: var1

    out_arr(OPS_ACC1(0,0,0)) = var1/(in_arr1(OPS_ACC2(0,0,0))*in_arr2(OPS_ACC3(0,0,0)))

END SUBROUTINE flamin_kernel_eqF

!A = var/B
SUBROUTINE flamin_kernel_eqG(out_arr, in_arr1, var1)
    use data_types
    implicit none

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr1
    real(kind=8), intent(in) :: var1

    out_arr(OPS_ACC1(0,0,0)) = var1/in_arr1(OPS_ACC2(0,0,0))

END SUBROUTINE flamin_kernel_eqG
