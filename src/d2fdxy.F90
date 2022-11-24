SUBROUTINE d2fdxy(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

    TYPE(ops_dat) :: functn, fderiv

!   ARGUMENTS
!   =========

!   LOCAL DATA
!   ==========
    integer :: istart,ifinis,jstart,jfinis
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================


    if(nyglbl == 1) then
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(d2fdxy_kernel_null, "d2fdxy_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
    else
!       END CONDITIONS
!       ==============
        istart = istal
        ifinis = istol
        jstart = jstal
        jfinis = jstol

        IF(nendxl == nbound)    istart = istap5
        IF(nendxr == nbound)    ifinis = istom5
        IF(nendyl == nbound)    jstart = jstap5
        IF(nendyr == nbound)    jfinis = jstom5

!       =========================================================================

!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz = (/istart,ifinis,jstart,jfinis,kstal,kstol/)
        call ops_par_loop(d2fdxy_kernel_interior, "d2fdxy_interior", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p550_to_m550_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!       =========================================================================

!       LH END X-DIRECTION
!       ==================
        IF(nendxl == nbound) THEN
!           TAKE SECOND XY-DERIVATIVE IN X-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT

!           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
            rangexyz = (/istal,istal,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_4th_onesided, "d2fdxy_lh_xdir_4th_onesided", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p420_m020_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
            rangexyz = (/istap1,istap1,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_4th_mixed, "d2fdxy_lh_xdir_4th_mixed", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p320_m120_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!           LH POINT PLUS 2: 4TH ORDER CENTRED
            rangexyz = (/istap2,istap2,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_4th_centered, "d2fdxy_lh_xdir_4th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p220_m220_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!           LH POINT PLUS 3: 6TH ORDER CENTRED
            rangexyz = (/istap3,istap3,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_6th_centered, "d2fdxy_lh_xdir_6th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p330_m330_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!           LH POINT PLUS 4: 8TH ORDER CENTRED
            rangexyz = (/istap4,istap4,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_8th_centered, "d2fdxy_lh_xdir_8th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        END IF

!       =========================================================================

!       RH END X-DIRECTION
!       ==================
        IF(nendxr == nbound)THEN
!           TAKE SECOND XY-DERIVATIVE IN X-RIGHT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT

!           RH POINT MINUS 4: 8TH ORDER CENTRED
            rangexyz = (/istom4,istom4,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_8th_centered, "d2fdxy_rh_xdir_8th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!           RH POINT MINUS 3: 6TH ORDER CENTRED
            rangexyz = (/istom3,istom3,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_6th_centered, "d2fdxy_rh_xdir_6th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p330_m330_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!           RH POINT MINUS 2: 4TH ORDER CENTRED
            rangexyz = (/istom2,istom2,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_4th_centered, "d2fdxy_rh_xdir_4th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p220_m220_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
            rangexyz = (/istom1,istom1,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_4th_mixed, "d2fdxy_rh_xdir_4th_mixed", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p120_m320_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

!           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
            rangexyz = (/istol,istol,jstart,jfinis,kstal,kstol/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_4th_onesided, "d2fdxy_rh_xdir_4th_onesided", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p020_m420_mixed_xy, "real(dp)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        END IF




    end if !nyglbl == 1
 
!   =========================================================================

END SUBROUTINE d2fdxy
