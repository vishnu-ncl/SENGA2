SUBROUTINE dfbydy(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

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
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

    IF (nyglbl == 1) THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(dfbydy_kernel_null, "dfbydy_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))
    ELSE
!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(dfbydy_kernel_main, "dfbydy_main_scheme", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(functn, 1, s3d_p050_to_m050_y, "real(kind=8)", OPS_READ),  &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_gbl(nyglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(nendyl, 1, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(nendyr, 1, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(nbound_ops, 1, "integer(kind=4)", OPS_READ), &
                        ops_arg_idx())

    END IF

!   =========================================================================

END SUBROUTINE dfbydy
