SUBROUTINE dfbydz(functn,fderiv)

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
INTEGER :: kstart,kfinis
INTEGER :: rangexyz(6)


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

kstart = kstal
kfinis = kstol
IF(nendzl == nbound)    kstart = kstap5
IF(nendzr == nbound)    kfinis = kstom5


!   =========================================================================
    if(nzglbl == 1) then
        rangexyz = (/istal,istol, jstal,jstol, kstal,kstol/)
        call ops_par_loop(dfbydz_kernel_null, "dfbydz_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
    else
!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz = (/istal,istol,jstal,jstol,kstart,kfinis/)
        call ops_par_loop(dfbydz_kernel_interior, "dfbydz_interior_scheme", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p005_to_m005_z, "real(dp)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

 
!       LH END
!       ======
        IF(nendzl == nbound)THEN

!           EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

            rangexyz = (/istal,istol,jstal,jstol,kstal,kstal/)
            call ops_par_loop(dfbydz_kernel_lhpoint_4th_onesided, "dfbydz_lh_4th_onesided", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_000_to_p004_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstal,jstol,kstap1,kstap1/)
            call ops_par_loop(dfbydz_kernel_lhpoint_4th_mixed, "dfbydz_lh_4th_mixed", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p003_to_m001_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstal,jstol,kstap2,kstap2/)
            call ops_par_loop(dfbydz_kernel_lhpoint_4th_centered, "dfbydz_lh_4th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p002_to_m002_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstal,jstol,kstap3,kstap3/)
            call ops_par_loop(dfbydz_kernel_lhpoint_6th_centered, "dfbydz_lh_6th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p003_to_m003_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstal,jstol,kstap4,kstap4/)
            call ops_par_loop(dfbydz_kernel_lhpoint_8th_centered, "dfbydz_lh_8th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p004_to_m004_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        END IF

!       RH END
!       ======
        IF(nendzr == nbound)THEN

            rangexyz = (/istal,istol,jstal,jstol,kstom4,kstom4/)
            call ops_par_loop(dfbydz_kernel_rhpoint_8th_centered, "dfbydz_rh_8th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p004_to_m004_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstal,jstol,kstom3,kstom3/)
            call ops_par_loop(dfbydz_kernel_rhpoint_6th_centered, "dfbydz_rh_6th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p003_to_m003_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstal,jstol,kstom2,kstom2/)
            call ops_par_loop(dfbydz_kernel_rhpoint_4th_centered, "dfbydz_rh_4th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p002_to_m002_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstal,jstol,kstom1,kstom1/)
            call ops_par_loop(dfbydz_kernel_rhpoint_4th_mixed, "dfbydz_rh_4th_mixed", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p001_to_m003_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstal,jstol,kstol,kstol/)
            call ops_par_loop(dfbydz_kernel_rhpoint_4th_onesided, "dfbydz_rh_4th_onesided", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_000_to_m004_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        END IF

!       =========================================================================

!       SCALING
!       =======
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(dfbydz_kernel_scaling, "dfbydz_scaling", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
!       =========================================================================


    end if


!   =========================================================================


END SUBROUTINE dfbydz
