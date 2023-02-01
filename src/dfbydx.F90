SUBROUTINE dfbydx(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   DFBYDX
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   28-MAR-2003:  RSC MODIFIED FOR SENGA2

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES FIRST X-DERIVATIVE OF SPECIFIED FUNCTION
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
    rangexyz = (/6,nxglbl-5,1,nyglbl,1,nzglbl/)
    call ops_par_loop(dfbydx_kernel_interior, "dfbydx_interior_scheme", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(functn, 1, s3d_p500_to_m500_x_no000, "real(8)", OPS_READ),  &
                      ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!   =========================================================================

!   LH END
!   ======
!   EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

        rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_lhpoint_4th_onesided, "dfbydx_lh_4th_onesided", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_000_to_p400_x, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/2,2,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_lhpoint_4th_mixed, "dfbydx_lh_4th_mixed", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p300_to_m100_x, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/3,3,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_lhpoint_4th_centered, "dfbydx_lh_4th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p200_to_m200_x_no000, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/4,4,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_lhpoint_6th_centered, "dfbydx_lh_6th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p300_to_m300_x_no000, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/5,5,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_lhpoint_8th_centered, "dfbydx_lh_8th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p400_to_m400_x_no000, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!   =========================================================================

!   RH END
!   ======
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

        rangexyz = (/nxglbl-4,nxglbl-4,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_rhpoint_8th_centered, "dfbydx_rh_8th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p400_to_m400_x_no000, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/nxglbl-3,nxglbl-3,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_rhpoint_6th_centered, "dfbydx_rh_6th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p300_to_m300_x_no000, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/nxglbl-2,nxglbl-2,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_rhpoint_4th_centered, "dfbydx_rh_4th_centered", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p200_to_m200_x_no000, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/nxglbl-1,nxglbl-1,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_rhpoint_4th_mixed, "dfbydx_rh_4th_mixed", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_p100_to_m300_x, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_rhpoint_4th_onesided, "dfbydx_rh_4th_onesided", senga_grid, 3, rangexyz,  &
                          ops_arg_dat(functn, 1, s3d_000_to_m400_x, "real(8)", OPS_READ),  &
                          ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!   =========================================================================

!   SCALING
!   =======
    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(dfbydx_kernel_scaling, "dfbydx_scaling", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_RW))
!   =========================================================================

END SUBROUTINE dfbydx
