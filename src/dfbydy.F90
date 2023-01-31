SUBROUTINE dfbydy(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   DFBYDY
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
!   EVALUATES FIRST Y-DERIVATIVE OF SPECIFIED FUNCTION

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
        call ops_par_loop(dfbydy_kernel_null, "dfbydy_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))
    ELSE
!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz = (/1,nxsize,6,nysize-5,1,nzsize/)
        call ops_par_loop(dfbydy_kernel_interior, "dfbydy_interior_scheme", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p050_to_m050_y_no000, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================
    
!       LH END
!       ======
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

            rangexyz = (/1,nxsize,1,1,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_lhpoint_4th_onesided, "dfbydy_lh_4th_onesided", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_000_to_p040_y, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxsize,2,2,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_lhpoint_4th_mixed, "dfbydy_lh_4th_mixed", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p030_to_m010_y, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxsize,3,3,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_lhpoint_4th_centered, "dfbydy_lh_4th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p020_to_m020_y_no000, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxsize,4,4,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_lhpoint_6th_centered, "dfbydy_lh_6th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p030_to_m030_y_no000, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxsize,5,5,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_lhpoint_8th_centered, "dfbydy_lh_8th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p040_to_m040_y_no000, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       RH END
!       ======
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

            rangexyz = (/1,nxsize,nysize-4,nysize-4,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_rhpoint_8th_centered, "dfbydy_rh_8th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p040_to_m040_y_no000, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxsize,nysize-3,nysize-3,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_rhpoint_6th_centered, "dfbydy_rh_6th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p030_to_m030_y_no000, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxsize,nysize-2,nysize-2,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_rhpoint_4th_centered, "dfbydy_rh_4th_centered", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p020_to_m020_y_no000, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxsize,nysize-1,nysize-1,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_rhpoint_4th_mixed, "dfbydy_rh_4th_mixed", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_p010_to_m030_y, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxsize,nysize,nysize,1,nzsize/)
            call ops_par_loop(dfbydy_kernel_rhpoint_4th_onesided, "dfbydy_rh_4th_onesided", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(functn, 1, s3d_000_to_m040_y, "real(8)", OPS_READ),  &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       SCALING
!       =======
        rangexyz = (/1,nxsize,1,nysize,1,nzsize/)
        call ops_par_loop(dfbydy_kernel_scaling, "dfbydy_scaling", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_RW))

!   =========================================================================

    END IF

!   =========================================================================

END SUBROUTINE dfbydy
