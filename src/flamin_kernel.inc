SUBROUTINE flamin_kernel_defc(crin,deltag,clocat,cthick,idx)
  real(kind=8), dimension(1) :: crin
  real(kind=8), intent(in) :: deltag,clocat,cthick
  integer(kind=4),dimension(3),intent(in):: idx

  real(kind=8) :: xcoord,argmnt

  real(kind=8) :: erfunc
  integer, parameter :: ncoeff = 5
  integer, parameter :: ncofm1 = ncoeff-1

  real(kind=8) :: ecoeff(ncoeff)
  real(kind=8) :: pcoeff
  real(kind=8) :: etotal,zvalue,tvalue
  integer :: ix,icoeff

! SET THE COEFFICIENTS
  pcoeff = 0.3275911_8
  ecoeff(1) = 0.254829592_8
  ecoeff(2) =-0.284496736_8
  ecoeff(3) = 1.421413741_8
  ecoeff(4) =-1.453152027_8
  ecoeff(5) = 1.061405429_8

  ix = idx(1)
  xcoord = real(ix-1,kind=8)*deltag

  argmnt = (xcoord-clocat)/cthick

! EVALUATE ERROR FUNCTION
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
END SUBROUTINE flamin_kernel_defc

SUBROUTINE flamin_kernel_eqA(yrun,crin,yrinr,yrinp)

   real(kind=8), dimension(1) :: yrun
   real(kind=8), dimension(1), intent(in) :: crin
   real(kind=8), intent(in) :: yrinr,yrinp

   yrun(OPS_ACC1(0,0,0))=yrinr+crin(OPS_ACC2(0,0,0))*(yrinp-yrinr)
END SUBROUTINE flamin_kernel_eqA

SUBROUTINE flamin_kernel_eqB(yrun_n,yrun_i)
   real(kind=8), dimension(1) :: yrun_n
   real(kind=8), dimension(1), intent(in) :: yrun_i

   yrun_n(OPS_ACC1(0,0,0))=yrun_n(OPS_ACC1(0,0,0))+yrun_i(OPS_ACC2(0,0,0))
END SUBROUTINE flamin_kernel_eqB

SUBROUTINE flamin_kernel_eqC(yrun)
   real(kind=8), dimension(1) :: yrun

   yrun(OPS_ACC1(0,0,0)) = 1.0_8-yrun(OPS_ACC1(0,0,0))
END SUBROUTINE flamin_kernel_eqC

SUBROUTINE flamin_kernel_eqD(store1,yrun,rgspec)
   real(kind=8), dimension(1) :: store1
   real(kind=8), dimension(1), intent(in) :: yrun
   real(kind=8), intent(in) :: rgspec

   store1(OPS_ACC1(0,0,0)) = store1(OPS_ACC1(0,0,0))+rgspec*yrun(OPS_ACC2(0,0,0))
END SUBROUTINE flamin_kernel_eqD

SUBROUTINE flamin_kernel_eqE(drun,store1,trun,prin)
   real(kind=8), dimension(1) :: drun
   real(kind=8), dimension(1), intent(in) :: store1,trun
   real(kind=8), intent(in) :: prin

   drun(OPS_ACC1(0,0,0))=prin/(store1(OPS_ACC2(0,0,0))*trun(OPS_ACC3(0,0,0)))
END SUBROUTINE flamin_kernel_eqE

SUBROUTINE flamin_kernel_eqF(urun,drun,flxmas)
   real(kind=8), dimension(1) :: urun
   real(kind=8), dimension(1), intent(in) :: drun
   real(kind=8), intent(in) :: flxmas

   urun(OPS_ACC1(0,0,0)) = flxmas/(drun(OPS_ACC2(0,0,0)))
END SUBROUTINE flamin_kernel_eqF
