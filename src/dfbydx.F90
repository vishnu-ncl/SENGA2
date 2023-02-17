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

        rangexyz = (/1,5,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_lhpoint, "dfbydx_lh", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p400_to_m400_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE), &
                        ops_arg_idx())
!   =========================================================================

!   RH END
!   ======
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT

        rangexyz = (/nxglbl-4,nxglbl,1,nyglbl,1,nzglbl/)
        call ops_par_loop(dfbydx_kernel_rhpoint, "dfbydx_rh", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p400_to_m400_x, "real(8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE), &
                        ops_arg_gbl(nxglbl, 1, "integer", OPS_READ), &
                        ops_arg_idx())
!   =========================================================================

!   SCALING
!   =======
    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(dfbydx_kernel_scaling, "dfbydx_scaling", senga_grid, 3, rangexyz,  &
                      ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_RW))
!   =========================================================================

END SUBROUTINE dfbydx
