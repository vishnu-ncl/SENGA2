SUBROUTINE d2fdx2(functn,fderiv)
 
use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   D2FDX2
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   06-APR-2003:  RSC MODIFIED FOR SENGA2

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES SECOND X-DERIVATIVE OF SPECIFIED FUNCTION
!   EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
!   EXPLICIT 8TH,6TH,4TH,4TH ORDER END CONDITIONS

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

!   INTERIOR SCHEME
!   ===============

!   TENTH ORDER EXPLICIT DIFFERENCES
    rangexyz(1) = 1
    rangexyz(2) = nxsize
    IF(nendxl == nbound) rangexyz(1) = 6
    IF(nendxr == nbound) rangexyz(2) = nxsize-5
 
    rangexyz(3) = 1
    rangexyz(4) = nysize
    rangexyz(5) = 1
    rangexyz(6) = nzsize
    
    call ops_par_loop(d2fdx2_kernel_interior, "d2fdx2_interior_scheme", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(functn, 1, s3d_p500_to_m500_x, "real(8)", OPS_READ),  &
                      ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!   =========================================================================

!   LH END
!   ======
    IF(nendxl == nbound)THEN

        rangexyz = (/1,1,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_lhpoint_4th_onesided, "d2fdx2_lh_4th_onesided", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_000_to_p500_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/2,2,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_lhpoint_4th_mixed, "d2fdx2_lh_4th_mixed", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p400_to_m100_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/3,3,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_lhpoint_4th_centered, "d2fdx2_lh_4th_centered", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p200_to_m200_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/4,4,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_lhpoint_6th_centered, "d2fdx2_lh_6th_centered", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p300_to_m300_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/5,5,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_lhpoint_8th_centered, "d2fdx2_lh_8th_centered", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p400_to_m400_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

    END IF

!   =========================================================================

!   RH END
!   ======
    IF(nendxr == nbound)THEN

        rangexyz = (/nxsize-4,nxsize-4,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_rhpoint_8th_centered, "d2fdx2_rh_8th_centered", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p400_to_m400_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/nxsize-3,nxsize-3,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_rhpoint_6th_centered, "d2fdx2_rh_6th_centered", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p300_to_m300_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/nxsize-2,nxsize-2,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_rhpoint_4th_centered, "d2fdx2_rh_4th_centered", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p200_to_m200_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/nxsize-1,nxsize-1,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_rhpoint_4th_mixed, "d2fdx2_rh_4th_mixed", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p100_to_m400_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/nxsize,nxsize,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdx2_kernel_rhpoint_4th_onesided, "d2fdx2_rh_4th_onesided", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_000_to_m500_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

    END IF

!   =========================================================================

!   SCALING
!   =======
    rangexyz = (/1,nxsize,1,nysize,1,nzsize/)
    call ops_par_loop(d2fdx2_kernel_scaling, "d2fdx2_scaling", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_RW))
!   =========================================================================

END SUBROUTINE d2fdx2
