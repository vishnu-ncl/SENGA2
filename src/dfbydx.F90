SUBROUTINE dfbydx(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

    TYPE(ops_dat) :: functn, fderiv


!     ARGUMENTS
!     =========




!     LOCAL DATA
!     ==========
INTEGER :: istart,ifinis
INTEGER :: rangexyz(6)


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

istart = istal
ifinis = istol
IF(nendxl == nbound)    istart = istap5
IF(nendxr == nbound)    ifinis = istom5

!   =========================================================================

!   INTERIOR SCHEME
!   ===============

!   TENTH ORDER EXPLICIT DIFFERENCES
    rangexyz = (/istart,ifinis,jstal,jstol,kstal,kstol/)
    call ops_par_loop(dfbydx_kernel_interior, "dfbydx_interior_scheme", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(functn, 1, s3d_p500_to_m500_x, "real(dp)", OPS_READ),  &
                      ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!   =========================================================================

!   LH END
!   ======
    IF(nendxl == nbound)THEN
  
!   EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

        rangexyz = (/istal,istal,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_lhpoint_4th_onesided, "dfbydx_lh_4th_onesided", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_000_to_p400_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        rangexyz = (/istap1,istap1,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_lhpoint_4th_mixed, "dfbydx_lh_4th_mixed", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p300_to_m100_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        rangexyz = (/istap2,istap2,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_lhpoint_4th_centered, "dfbydx_lh_4th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p200_to_m200_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        rangexyz = (/istap3,istap3,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_lhpoint_6th_centered, "dfbydx_lh_6th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p300_to_m300_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        rangexyz = (/istap4,istap4,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_lhpoint_8th_centered, "dfbydx_lh_8th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p400_to_m400_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

    END IF

!   =========================================================================

!   RH END
!   ======
    IF(nendxr == nbound)THEN

        rangexyz = (/istom4,istom4,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_rhpoint_8th_centered, "dfbydx_rh_8th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p400_to_m400_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        rangexyz = (/istom3,istom3,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_rhpoint_6th_centered, "dfbydx_rh_6th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p300_to_m300_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        rangexyz = (/istom2,istom2,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_rhpoint_4th_centered, "dfbydx_rh_4th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p200_to_m200_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        rangexyz = (/istom1,istom1,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_rhpoint_4th_mixed, "dfbydx_rh_4th_mixed", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p100_to_m300_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        rangexyz = (/istol,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydx_kernel_rhpoint_4th_onesided, "dfbydx_rh_4th_onesided", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_000_to_m400_x, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

    END IF

!   =========================================================================

!   SCALING
!   =======
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(dfbydx_kernel_scaling, "dfbydx_scaling", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
!   =========================================================================


END SUBROUTINE dfbydx
