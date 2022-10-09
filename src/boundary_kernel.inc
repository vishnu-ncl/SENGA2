SUBROUTINE boundary_kernel_CPandGAS_xdir(transp,store7,drhs,strgx,strrx)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=dp), dimension(1) :: strgx,strrx

    strgx(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strrx(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_CPandGAS_xdir

SUBROUTINE boundary_kernel_CPandGAS_ydir(transp,store7,drhs,strgy,strry)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=dp), dimension(1) :: strgy,strry

    strgy(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strry(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_CPandGAS_ydir

SUBROUTINE boundary_kernel_CPandGAS_zdir(transp,store7,drhs,strgz,strrz)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=dp), dimension(1) :: strgz,strrz

    strgz(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strrz(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE boundary_kernel_CPandGAS_zdir

SUBROUTINE boundary_kernel_internalenergy_xdir(erhs,strex)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: erhs
    real(kind=dp), dimension(1) :: strex
    
    strex(OPS_ACC2(0,0,0)) = erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_internalenergy_xdir

SUBROUTINE boundary_kernel_internalenergy_ydir(erhs,strey)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: erhs
    real(kind=dp), dimension(1) :: strey

    strey(OPS_ACC2(0,0,0)) = erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_internalenergy_ydir

SUBROUTINE boundary_kernel_internalenergy_zdir(erhs,strez)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: erhs
    real(kind=dp), dimension(1) :: strez

    strez(OPS_ACC2(0,0,0)) = erhs(OPS_ACC1(0,0,0))

END SUBROUTINE boundary_kernel_internalenergy_zdir

SUBROUTINE boundary_kernel_temperature_xdir(trun,store1,strtx,bcltx)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: trun,store1
    real(kind=dp), dimension(1) :: strtx,bcltx

    strtx(OPS_ACC3(0,0,0)) = trun(OPS_ACC1(0,0,0))
    bcltx(OPS_ACC4(0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_temperature_xdir

SUBROUTINE boundary_kernel_temperature_ydir(trun,store1,strty,bclty)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: trun,store1
    real(kind=dp), dimension(1) :: strty,bclty

    strty(OPS_ACC3(0,0,0)) = trun(OPS_ACC1(0,0,0))
    bclty(OPS_ACC4(0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_temperature_ydir

SUBROUTINE boundary_kernel_temperature_zdir(trun,store1,strtz,bcltz)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: trun,store1
    real(kind=dp), dimension(1) :: strtz,bcltz

    strtz(OPS_ACC3(0,0,0)) = trun(OPS_ACC1(0,0,0))
    bcltz(OPS_ACC4(0,0,0)) = store1(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_temperature_zdir

SUBROUTINE boundary_kernel_density_xdir(drhs,store,strdx,bcl2x)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: drhs,store
    real(kind=dp), dimension(1) :: strdx,bcl2x

    strdx(OPS_ACC3(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    bcl2x(OPS_ACC4(0,0,0)) = store(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_density_xdir

SUBROUTINE boundary_kernel_density_ydir(drhs,store,strdy,bcl2y)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: drhs,store
    real(kind=dp), dimension(1) :: strdy,bcl2y

    strdy(OPS_ACC3(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    bcl2y(OPS_ACC4(0,0,0)) = store(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_density_ydir

SUBROUTINE boundary_kernel_density_zdir(drhs,store,strdz,bcl2z)
    use data_types
    implicit none

    real(kind=dp), dimension(1), intent(in) :: drhs,store
    real(kind=dp), dimension(1) :: strdz,bcl2z

    strdz(OPS_ACC3(0,0,0)) = drhs(OPS_ACC1(0,0,0))
    bcl2z(OPS_ACC4(0,0,0)) = store(OPS_ACC2(0,0,0))

END SUBROUTINE boundary_kernel_density_zdir
