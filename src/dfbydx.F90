SUBROUTINE dfbydx(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

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
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   INTERIOR SCHEME
!   ===============

!   TENTH ORDER EXPLICIT DIFFERENCES
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(dfbydx_kernel_main, "dfbydx_main_scheme", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(functn, 1, s3d_p500_to_m500_x, "real(kind=8)", OPS_READ),  &
                    ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(nxglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                    ops_arg_gbl(nendxl, 1, "integer(kind=4)", OPS_READ), &
                    ops_arg_gbl(nendxr, 1, "integer(kind=4)", OPS_READ), &
                    ops_arg_gbl(nbound_ops, 1, "integer(kind=4)", OPS_READ), &
                    ops_arg_idx())

!   =========================================================================

END SUBROUTINE dfbydx
