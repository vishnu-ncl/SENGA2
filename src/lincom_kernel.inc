SUBROUTINE lincom_kernel_main(err_arr,run_arr,rhs_arr,rkerr,rklhs,rkrhs,irkstp)

    real(kind=8), dimension(1) :: err_arr, run_arr, rhs_arr
    integer(kind=4), intent(in) :: irkstp
    real(kind=8), dimension(nrkmax), intent(in) :: rkerr, rklhs, rkrhs
    real(kind=8) :: fornow

    err_arr(OPS_ACC1(0,0,0)) = err_arr(OPS_ACC1(0,0,0)) + rkerr(irkstp) * rhs_arr(OPS_ACC3(0,0,0))

    fornow = run_arr(OPS_ACC2(0,0,0))
    run_arr(OPS_ACC2(0,0,0)) = fornow + rklhs(irkstp) * rhs_arr(OPS_ACC3(0,0,0))
    rhs_arr(OPS_ACC3(0,0,0)) = fornow + rkrhs(irkstp) * rhs_arr(OPS_ACC3(0,0,0))

END SUBROUTINE lincom_kernel_main

SUBROUTINE lincom_kernel_eqA(yrun,yrhs,drhs)

    real(kind=8), dimension(1) :: yrun
    real(kind=8), dimension(1), intent(in) :: yrhs,drhs
    real(kind=8) :: temp1,temp2,temp3,temp4

    temp1 = yrhs(OPS_ACC2(1,0,0))/drhs(OPS_ACC3(1,0,0))
    temp2 = yrhs(OPS_ACC2(2,0,0))/drhs(OPS_ACC3(2,0,0))
    temp3 = yrhs(OPS_ACC2(3,0,0))/drhs(OPS_ACC3(3,0,0))
    temp4 = yrhs(OPS_ACC2(4,0,0))/drhs(OPS_ACC3(4,0,0))

    yrun(OPS_ACC1(0,0,0)) = (12.0_8/25.0_8)*(4.0_8*temp1 - 3.0_8*temp2 +  &
                            (4.0_8/3.0_8)*temp3 - (1.0_8/4.0_8*temp4))

    yrun(OPS_ACC1(0,0,0))=MIN(1.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0))=MAX(0.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0))=drhs(OPS_ACC3(0,0,0)) *yrun(OPS_ACC1(0,0,0))

END SUBROUTINE lincom_kernel_eqA

SUBROUTINE lincom_kernel_eqB(yrun,yrhs,drhs)

    real(kind=8), dimension(1) :: yrun
    real(kind=8), dimension(1), intent(in) :: yrhs,drhs
    real(kind=8) :: temp1,temp2,temp3,temp4

    temp1 = yrhs(OPS_ACC2(-1,0,0))/drhs(OPS_ACC3(-1,0,0))
    temp2 = yrhs(OPS_ACC2(-2,0,0))/drhs(OPS_ACC3(-2,0,0))
    temp3 = yrhs(OPS_ACC2(-3,0,0))/drhs(OPS_ACC3(-3,0,0))
    temp4 = yrhs(OPS_ACC2(-4,0,0))/drhs(OPS_ACC3(-4,0,0))

    yrun(OPS_ACC1(0,0,0)) = (12.0_8/25.0_8)*(4.0_8*temp1 - 3.0_8*temp2 +  &
                              (4.0_8/3.0_8)*temp3 - (1.0_8/4.0_8*temp4))

    yrun(OPS_ACC1(0,0,0)) = MIN(1.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0)) = MAX(0.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0)) = drhs(OPS_ACC3(0,0,0)) *yrun(OPS_ACC1(0,0,0))

END SUBROUTINE lincom_kernel_eqB

SUBROUTINE lincom_kernel_eqC(yrun,yrhs,drhs)

    real(kind=8), dimension(1) :: yrun
    real(kind=8), dimension(1), intent(in) :: yrhs,drhs
    real(kind=8) :: temp1,temp2,temp3,temp4

    temp1 = yrhs(OPS_ACC2(0,1,0))/drhs(OPS_ACC3(0,1,0))
    temp2 = yrhs(OPS_ACC2(0,2,0))/drhs(OPS_ACC3(0,2,0))
    temp3 = yrhs(OPS_ACC2(0,3,0))/drhs(OPS_ACC3(0,3,0))
    temp4 = yrhs(OPS_ACC2(0,4,0))/drhs(OPS_ACC3(0,4,0))

    yrun(OPS_ACC1(0,0,0)) = (12.0_8/25.0_8)*(4.0_8*temp1 - 3.0_8*temp2 +  &
                              (4.0_8/3.0_8)*temp3 - (1.0_8/4.0_8*temp4))

    yrun(OPS_ACC1(0,0,0)) = MIN(1.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0)) = MAX(0.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0)) = drhs(OPS_ACC3(0,0,0)) *yrun(OPS_ACC1(0,0,0))

END SUBROUTINE lincom_kernel_eqC

SUBROUTINE lincom_kernel_eqD(yrun,yrhs,drhs)

    real(kind=8), dimension(1) :: yrun
    real(kind=8), dimension(1), intent(in) :: yrhs,drhs
    real(kind=8) :: temp1,temp2,temp3,temp4

    temp1 = yrhs(OPS_ACC2(0,-1,0))/drhs(OPS_ACC3(0,-1,0))
    temp2 = yrhs(OPS_ACC2(0,-2,0))/drhs(OPS_ACC3(0,-2,0))
    temp3 = yrhs(OPS_ACC2(0,-3,0))/drhs(OPS_ACC3(0,-3,0))
    temp4 = yrhs(OPS_ACC2(0,-4,0))/drhs(OPS_ACC3(0,-4,0))

    yrun(OPS_ACC1(0,0,0)) = (12.0_8/25.0_8)*(4.0_8*temp1 - 3.0_8*temp2 +  &
                              (4.0_8/3.0_8)*temp3 - (1.0_8/4.0_8*temp4))

    yrun(OPS_ACC1(0,0,0)) = MIN(1.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0)) = MAX(0.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0)) = drhs(OPS_ACC3(0,0,0)) *yrun(OPS_ACC1(0,0,0))

END SUBROUTINE lincom_kernel_eqD

SUBROUTINE lincom_kernel_eqE(yrun,yrhs,drhs)

    real(kind=8), dimension(1) :: yrun
    real(kind=8), dimension(1), intent(in) :: yrhs,drhs
    real(kind=8) :: temp1,temp2,temp3,temp4

    temp1 = yrhs(OPS_ACC2(0,0,1))/drhs(OPS_ACC3(0,0,1))
    temp2 = yrhs(OPS_ACC2(0,0,2))/drhs(OPS_ACC3(0,0,2))
    temp3 = yrhs(OPS_ACC2(0,0,3))/drhs(OPS_ACC3(0,0,3))
    temp4 = yrhs(OPS_ACC2(0,0,4))/drhs(OPS_ACC3(0,0,4))

    yrun(OPS_ACC1(0,0,0)) = (12.0_8/25.0_8)*(4.0_8*temp1 - 3.0_8*temp2 +  &
            (4.0_8/3.0_8)*temp3 - (1.0_8/4.0_8*temp4))

    yrun(OPS_ACC1(0,0,0)) = MIN(1.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0)) = MAX(0.0_8,yrun(OPS_ACC1(0,0,0)))

    yrun(OPS_ACC1(0,0,0)) = drhs(OPS_ACC3(0,0,0)) *yrun(OPS_ACC1(0,0,0))

END SUBROUTINE lincom_kernel_eqE

SUBROUTINE lincom_kernel_eqF(yrun,yrhs,drhs)

    real(kind=8), dimension(1) :: yrun
    real(kind=8), dimension(1), intent(in) :: yrhs,drhs
    real(kind=8) :: temp1,temp2,temp3,temp4

    temp1 = yrhs(OPS_ACC2(0,0,-1))/drhs(OPS_ACC3(0,0,-1))
    temp2 = yrhs(OPS_ACC2(0,0,-2))/drhs(OPS_ACC3(0,0,-2))
    temp3 = yrhs(OPS_ACC2(0,0,-3))/drhs(OPS_ACC3(0,0,-3))
    temp4 = yrhs(OPS_ACC2(0,0,-4))/drhs(OPS_ACC3(0,0,-4))

    yrun(OPS_ACC1(0,0,0)) = (12.0_8/25.0_8)*(4.0_8*temp1 - 3.0_8*temp2 +  &
                        (4.0_8/3.0_8)*temp3 - (1.0_8/4.0_8*temp4))

    yrun(OPS_ACC1(0,0,0)) = MIN(1.0_8,yrun(OPS_ACC1(0,0,0)))
    yrun(OPS_ACC1(0,0,0)) = MAX(0.0_8,yrun(OPS_ACC1(0,0,0)))
    yrun(OPS_ACC1(0,0,0)) = drhs(OPS_ACC3(0,0,0)) *yrun(OPS_ACC1(0,0,0))

END SUBROUTINE lincom_kernel_eqF


