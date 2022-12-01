SUBROUTINE radcal_kernel_addradiation(erhs,trun,store1)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: erhs
    real(kind=dp), dimension(1), intent(in) :: trun,store1
    real(kind=dp) :: fornow

    fornow = trun(OPS_ACC2(0,0,0))
    fornow = fornow*fornow*fornow*fornow

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
          - foursb*store1(OPS_ACC3(0,0,0))*(fornow - trfrth)
    
END SUBROUTINE radcal_kernel_addradiation

SUBROUTINE radcal_kernel_addspecies(store1,yrhs,trun,store2,rgspec,ispec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: store1
    real(kind=dp), dimension(1), intent(in) :: yrhs,trun,store2
    real(kind=dp), intent(in) :: rgspec
    integer, intent(in) :: ispec
    real(kind=dp) :: fornow
    
    fornow = yrhs(OPS_ACC_MD2(ispec,0,0,0))*rgspec*trun(OPS_ACC3(0,0,0))
    
    store1(OPS_ACC1(0,0,0)) = store1(OPS_ACC1(0,0,0)) + store2(OPS_ACC4(0,0,0))*fornow

END SUBROUTINE radcal_kernel_addspecies

SUBROUTINE radcal_kernel_meancoef(store2,trun,akprad,nkprad,nkprm1,jspec,ncfrmx,nspec)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: store2
    real(kind=dp), dimension(1), intent(in) :: trun

    integer, intent(in) :: jspec, ncfrmx, nspec

    real(kind=dp), dimension(ncfrmx,nspec), intent(in) :: akprad
    integer, dimension(nspec), intent(in) :: nkprad, nkprm1

    real(kind=dp) :: plspec,fornow
    integer :: icp

    fornow = trun(OPS_ACC2(0,0,0))
    plspec = akprad(nkprad(jspec),jspec)
    DO icp = nkprm1(jspec),1,-1
        plspec = plspec*fornow + akprad(icp,jspec)
    END DO
    store2(OPS_ACC1(0,0,0)) = plspec

END SUBROUTINE radcal_kernel_meancoef
