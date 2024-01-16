SUBROUTINE output_kernel_eqA(drun,utgv,vtgv,wtgv,tkes)
  use data_types
  implicit none

  real(kind=8),dimension(1),intent(in) :: drun,utgv,vtgv,wtgv
  real(kind=8) :: tkes

  tkes=tkes+(0.5d0*drun(OPS_ACC1(0,0,0))*(utgv(OPS_ACC1(0,0,0))*utgv(OPS_ACC1(0,0,0))+ &
     vtgv(OPS_ACC1(0,0,0))*vtgv(OPS_ACC1(0,0,0))+wtgv(OPS_ACC1(0,0,0))*wtgv(OPS_ACC1(0,0,0)))&
	  *(1.0d0/ovdelx)*(1.0d0/ovdelx)*(1.0d0/ovdelx))

END SUBROUTINE output_kernel_eqA
