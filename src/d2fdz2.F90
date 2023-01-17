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
        rangexyz = (/1,nxsize, 1,nysize, 1,nzsize/)
        call ops_par_loop(d2fdz2_kernel_null, "d2fdz2_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
    ELSE
!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz(1) = 1
        rangexyz(2) = nxsize
        rangexyz(3) = 1
        rangexyz(4) = nysize

        rangexyz(5) = 1
        rangexyz(6) = nzsize
        IF(nendzl == nbound)    rangexyz(5) = 6
        IF(nendzr == nbound)    rangexyz(6) = nzsize-5
    
        call ops_par_loop(d2fdz2_kernel_interior, "d2fdz2_interior_scheme", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p005_to_m005_z, "real(dp)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))


!       LH END
!       ======
        IF(nendzl == nbound)THEN

            rangexyz(5) = 1
            rangexyz(6) = 1
            call ops_par_loop(d2fdz2_kernel_lhpoint_4th_onesided, "d2fdz2_lh_4th_onesided", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_000_to_p005_z, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz(5) = 2
            rangexyz(6) = 2
            call ops_par_loop(d2fdz2_kernel_lhpoint_4th_mixed, "d2fdz2_lh_4th_mixed", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p004_to_m001_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz(5) = 3
            rangexyz(6) = 3
            call ops_par_loop(d2fdz2_kernel_lhpoint_4th_centered, "d2fdz2_lh_4th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p002_to_m002_z, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz(5) = 4
            rangexyz(6) = 4
            call ops_par_loop(d2fdz2_kernel_lhpoint_6th_centered, "d2fdz2_lh_6th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p003_to_m003_z, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz(5) = 5
            rangexyz(6) = 5
            call ops_par_loop(d2fdz2_kernel_lhpoint_8th_centered, "d2fdz2_lh_8th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p004_to_m004_z, "real(dp)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        END IF

!       RH END
!       ======
        IF(nendzr == nbound)THEN

            rangexyz(5) = nzsize-4
            rangexyz(6) = nzsize-4
            call ops_par_loop(d2fdz2_kernel_rhpoint_8th_centered, "d2fdz2_rh_8th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p004_to_m004_z, "real(dp)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz(5) = nzsize-3
            rangexyz(6) = nzsize-3
            call ops_par_loop(d2fdz2_kernel_rhpoint_6th_centered, "d2fdz2_rh_6th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p003_to_m003_z, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz(5) = nzsize-2
            rangexyz(6) = nzsize-2
            call ops_par_loop(d2fdz2_kernel_rhpoint_4th_centered, "d2fdz2_rh_4th_centered", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p002_to_m002_z, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz(5) = nzsize-1
            rangexyz(6) = nzsize-1
            call ops_par_loop(d2fdz2_kernel_rhpoint_4th_mixed, "d2fdz2_rh_4th_mixed", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_p001_to_m004_z, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

            rangexyz(5) = nzsize
            rangexyz(6) = nzsize
            call ops_par_loop(d2fdz2_kernel_rhpoint_4th_onesided, "d2fdz2_rh_4th_onesided", senga_grid, 3, rangexyz,  &
                              ops_arg_dat(functn, 1, s3d_000_to_m005_z, "real(dp)", OPS_READ),  &
                              ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))

        END IF

!       =========================================================================

!       SCALING
!       =======
        rangexyz = (/1,nxsize,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdz2_kernel_scaling, "d2fdz2_scaling", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
!       =========================================================================

    end if

!   =========================================================================

END SUBROUTINE d2fdz2
