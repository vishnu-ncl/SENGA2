SUBROUTINE d2fdz2(functn,fderiv)
 
use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   D2FDZ2
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   11-APR-2003:  RSC MODIFIED FOR SENGA2
!   10-OCT-2004:  RSC NULL VERSION

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES SECOND Z-DERIVATIVE OF SPECIFIED FUNCTION
!   EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
!   EXPLICIT 8TH,6TH,4TH,2ND ORDER END CONDITIONS

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

    IF (nyglbl == 1) THEN
        rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
        call ops_par_loop(d2fdz2_kernel_null, "d2fdz2_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

    ELSE
!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES

        rangexyz = (/1,nxglbl,1,nyglbl,6,nzglbl-5/)
        call ops_par_loop(d2fdz2_kernel_interior, "d2fdz2_interior_scheme", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p005_to_m005_z, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       LH END
!       ======
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

            rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
            call ops_par_loop(d2fdz2_kernel_lhpoint_4th_onesided, "d2fdz2_lh_4th_onesided", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_000_to_p005_z, "real(8)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,nyglbl,2,2/)
            call ops_par_loop(d2fdz2_kernel_lhpoint_4th_mixed, "d2fdz2_lh_4th_mixed", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p004_to_m001_z, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,nyglbl,3,3/)
            call ops_par_loop(d2fdz2_kernel_lhpoint_4th_centered, "d2fdz2_lh_4th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p002_to_m002_z, "real(8)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,nyglbl,4,4/)
            call ops_par_loop(d2fdz2_kernel_lhpoint_6th_centered, "d2fdz2_lh_6th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p003_to_m003_z, "real(8)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,nyglbl,5,5/)
            call ops_par_loop(d2fdz2_kernel_lhpoint_8th_centered, "d2fdz2_lh_8th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p004_to_m004_z, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       RH END
!       ======
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl-4,nzglbl-4/)
            call ops_par_loop(d2fdz2_kernel_rhpoint_8th_centered, "d2fdz2_rh_8th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p004_to_m004_z, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl-3,nzglbl-3/)
            call ops_par_loop(d2fdz2_kernel_rhpoint_6th_centered, "d2fdz2_rh_6th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p003_to_m003_z, "real(8)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl-2,nzglbl-2/)
            call ops_par_loop(d2fdz2_kernel_rhpoint_4th_centered, "d2fdz2_rh_4th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p002_to_m002_z, "real(8)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl-1,nzglbl-1/)
            call ops_par_loop(d2fdz2_kernel_rhpoint_4th_mixed, "d2fdz2_rh_4th_mixed", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p001_to_m004_z, "real(8)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
            call ops_par_loop(d2fdz2_kernel_rhpoint_4th_onesided, "d2fdz2_rh_4th_onesided", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_000_to_m005_z, "real(8)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       SCALING
!       =======
        rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
        call ops_par_loop(d2fdz2_kernel_scaling, "d2fdz2_scaling", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_RW))

!   =========================================================================

    END IF

!   =========================================================================

END SUBROUTINE d2fdz2
