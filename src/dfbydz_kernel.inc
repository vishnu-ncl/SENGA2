SUBROUTINE dfbydz_kernel_null(fderiv)
    use data_types
    implicit none

    real(kind=dp), dimension(1) :: fderiv

    fderiv(OPS_ACC1(0,0,0)) = 0.0_dp

END SUBROUTINE dfbydz_kernel_null
