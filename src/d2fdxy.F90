SUBROUTINE d2fdxy(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   D2FDXY
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   15-APR-2003:  RSC MODIFIED FOR SENGA2
!   10-OCT-2004:  RSC NULL VERSION

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES SECOND XY-DERIVATIVE OF SPECIFIED FUNCTION

!   *************************************************************************

!   ARGUMENTS
!   =========
    TYPE(ops_dat) :: functn, fderiv

!   LOCAL DATA
!   ==========
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

    IF (nysize == 1) THEN
        rangexyz = (/1,nxsize,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdxy_kernel_null, "d2fdxy_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

    ELSE
!       =========================================================================

!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES

        rangexyz = (/6,nxsize-5,6,nysize-5,1,nzsize/)
        call ops_par_loop(d2fdxy_kernel_interior, "d2fdxy_interior", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p550_to_m550_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       LH END X-DIRECTION
!       ==================
!           TAKE SECOND XY-DERIVATIVE IN X-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT

!           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
            rangexyz = (/1,1,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_4th_onesided, "d2fdxy_lh_xdir_4th_onesided", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p420_m020_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
            rangexyz = (/2,2,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_4th_mixed, "d2fdxy_lh_xdir_4th_mixed", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p320_m120_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           LH POINT PLUS 2: 4TH ORDER CENTRED
            rangexyz = (/3,3,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_4th_centered, "d2fdxy_lh_xdir_4th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p220_m220_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           LH POINT PLUS 3: 6TH ORDER CENTRED
            rangexyz = (/4,4,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_6th_centered, "d2fdxy_lh_xdir_6th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p330_m330_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           LH POINT PLUS 4: 8TH ORDER CENTRED
            rangexyz = (/5,5,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_lh_xdir_8th_centered, "d2fdxy_lh_xdir_8th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       RH END X-DIRECTION
!       ==================
!           TAKE SECOND XY-DERIVATIVE IN X-RIGHT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT

!           RH POINT MINUS 4: 8TH ORDER CENTRED
            rangexyz = (/nxsize-4,nxsize-4,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_8th_centered, "d2fdxy_rh_xdir_8th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           RH POINT MINUS 3: 6TH ORDER CENTRED
            rangexyz = (/nxsize-3,nxsize-3,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_6th_centered, "d2fdxy_rh_xdir_6th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p330_m330_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           RH POINT MINUS 2: 4TH ORDER CENTRED
            rangexyz = (/nxsize-2,nxsize-2,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_4th_centered, "d2fdxy_rh_xdir_4th_centred", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p220_m220_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
            rangexyz = (/nxsize-1,nxsize-1,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_4th_mixed, "d2fdxy_rh_xdir_4th_mixed", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p120_m320_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
            rangexyz = (/nxsize,nxsize,6,nysize-5,1,nzsize/)
            call ops_par_loop(d2fdxy_kernel_rh_xdir_4th_onesided, "d2fdxy_rh_xdir_4th_onesided", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p020_m420_mixed_xy, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))





    END IF
 
!   =========================================================================

END SUBROUTINE d2fdxy
