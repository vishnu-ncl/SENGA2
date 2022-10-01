SUBROUTINE d2fdy2(functn,fderiv)
 
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
    INTEGER :: jstart,jfinis
    INTEGER :: rangexyz(6)


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

jstart = jstal
jfinis = jstol
IF(nendyl == nbound)    jstart = jstap5
IF(nendyr == nbound)    jfinis = jstom5


!   =========================================================================
    if(nyglbl == 1) then
        rangexyz = (/istal,istol, jstal,jstol, kstal,kstol/)
        call ops_par_loop(d2fdy2_kernel_null, "d2fdy2_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
    else
!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz = (/istal,istol,jstart,jfinis,kstal,kstol/)
        call ops_par_loop(d2fdy2_kernel_interior, "d2fdy2_interior_scheme", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p050_to_m050_y, "real(dp)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))


!       LH END
!       ======
        IF(nendyl == nbound)THEN

            rangexyz = (/istal,istol,jstal,jstal,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_lhpoint_4th_onesided, "d2fdy2_lh_4th_onesided", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_000_to_p050_y, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstap1,jstap1,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_lhpoint_4th_mixed, "d2fdy2_lh_4th_mixed", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p040_to_m010_y, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstap2,jstap2,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_lhpoint_4th_centered, "d2fdy2_lh_4th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p020_to_m020_y, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstap3,jstap3,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_lhpoint_6th_centered, "d2fdy2_lh_6th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p030_to_m030_y, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstap4,jstap4,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_lhpoint_8th_centered, "d2fdy2_lh_8th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p040_to_m040_y, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        END IF

!       RH END
!       ======
        IF(nendyr == nbound)THEN

            rangexyz = (/istal,istol,jstom4,jstom4,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_rhpoint_8th_centered, "d2fdy2_rh_8th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p040_to_m040_y, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstom3,jstom3,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_rhpoint_6th_centered, "d2fdy2_rh_6th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p030_to_m030_y, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstom2,jstom2,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_rhpoint_4th_centered, "d2fdy2_rh_4th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p020_to_m020_y, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstom1,jstom1,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_rhpoint_4th_mixed, "d2fdy2_rh_4th_mixed", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p010_to_m040_y, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz = (/istal,istol,jstol,jstol,kstal,kstol/)
            call ops_par_loop(d2fdy2_kernel_rhpoint_4th_onesided, "d2fdy2_rh_4th_onesided", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_000_to_m050_y, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        END IF

!       =========================================================================

!       SCALING
!       =======
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(d2fdy2_kernel_scaling, "d2fdy2_scaling", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
!       =========================================================================


    end if


!   =========================================================================


END SUBROUTINE d2fdy2
