SUBROUTINE flamin_kernel_set_reaction_var(crin,deltag,clocat,cthick,idx)

    real(kind=8), dimension(1) :: crin
    real(kind=8), intent(in) :: deltag,clocat,cthick
    integer(kind=4), dimension(3), intent(in) :: idx

    real(kind=8) :: xcoord,argmnt

    real(kind=8) :: erfunc
    integer(kind=4), parameter :: ncoeff = 5
    integer(kind=4), parameter :: ncofm1 = ncoeff-1

    real(kind=8) :: ecoeff(ncoeff)
    real(kind=8) :: pcoeff
    real(kind=8) :: etotal,zvalue,tvalue
    integer(kind=4) :: ix,icoeff

!   SET THE COEFFICIENTS
    pcoeff = 0.3275911_8
    ecoeff(1) = 0.254829592_8
    ecoeff(2) =-0.284496736_8
    ecoeff(3) = 1.421413741_8
    ecoeff(4) =-1.453152027_8
    ecoeff(5) = 1.061405429_8

    ix = idx(1)
    xcoord = REAL(ix-1,kind=8)*deltag
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
    IF (argmnt < 0.0_8) erfunc = -erfunc

    crin(OPS_ACC1(0,0,0)) = 0.5_8*(1.0_8+erfunc)

END SUBROUTINE flamin_kernel_set_reaction_var

!A = var1 + B*(var2-var1)
SUBROUTINE flamin_kernel_eqA(yrun,crin,yrinr,yrinp,ispec)

    real(kind=8), dimension(1) :: yrun
    real(kind=8), dimension(1), intent(in) :: crin
    real(kind=8), dimension(nspcmx), intent(in) ::  yrinr,yrinp
    integer(kind=4), intent(in) :: ispec

    yrun(OPS_ACC1(0,0,0)) = yrinr(ispec) + crin(OPS_ACC2(0,0,0))*(yrinp(ispec) - yrinr(ispec))

END SUBROUTINE flamin_kernel_eqA

!A = A + B
SUBROUTINE flamin_kernel_eqB(yrun_n,yrun_i)

    real(kind=8), dimension(1) :: yrun_n
    real(kind=8), dimension(1), intent(in) :: yrun_i

    yrun_n(OPS_ACC1(0,0,0)) = yrun_n(OPS_ACC1(0,0,0)) + yrun_i(OPS_ACC2(0,0,0))

END SUBROUTINE flamin_kernel_eqB

!A = 1.0 - A
SUBROUTINE flamin_kernel_eqC(yrun_n)

    real(kind=8), dimension(1) :: yrun_n

    yrun_n(OPS_ACC1(0,0,0)) = 1.0_8 - yrun_n(OPS_ACC1(0,0,0))

END SUBROUTINE flamin_kernel_eqC

!A = trinr + B_xdimonly*(trinp - trinr)
SUBROUTINE flamin_kernel_eqD(trun,crin,trinr,trinp)

    real(kind=8), dimension(1) :: trun
    real(kind=8), dimension(1), intent(in) :: crin
    real(kind=8), intent(in) :: trinr,trinp

    trun(OPS_ACC1(0,0,0)) = trinr + crin(OPS_ACC2(0,0,0))*(trinp - trinr)

END SUBROUTINE flamin_kernel_eqD

!A = A + var*B
SUBROUTINE flamin_kernel_eqE(store1,yrun,rgspec,ispec)

    real(kind=8), dimension(1) :: store1
    real(kind=8), dimension(1), intent(in) :: yrun
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), intent(in) :: ispec

    store1(OPS_ACC1(0,0,0)) = store1(OPS_ACC1(0,0,0)) + rgspec(ispec)*yrun(OPS_ACC2(0,0,0))

END SUBROUTINE flamin_kernel_eqE

!A = var/(B*C)
SUBROUTINE flamin_kernel_eqF(drun,store1,trun,prin)

    real(kind=8), dimension(1) :: drun
    real(kind=8), dimension(1), intent(in) :: store1, trun
    real(kind=8), intent(in) :: prin

    drun(OPS_ACC1(0,0,0)) = prin/(store1(OPS_ACC2(0,0,0))*trun(OPS_ACC3(0,0,0)))

END SUBROUTINE flamin_kernel_eqF

!A = var/B
SUBROUTINE flamin_kernel_eqG(urun,drun,flxmas)

    real(kind=8), dimension(1) :: urun
    real(kind=8), dimension(1), intent(in) :: drun
    real(kind=8), intent(in) :: flxmas

    urun(OPS_ACC1(0,0,0)) = flxmas/drun(OPS_ACC2(0,0,0))

END SUBROUTINE flamin_kernel_eqG
