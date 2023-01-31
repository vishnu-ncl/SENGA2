SUBROUTINE d2fdyz(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   D2FDYZ
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   15-MAY-2003:  RSC MODIFIED FOR SENGA2
!   10-OCT-2004:  RSC NULL VERSION
!   31-DEC-2006:  RSC BUG FIX INDICES

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES SECOND YZ-DERIVATIVE OF SPECIFIED FUNCTION

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
        call ops_par_loop(d2fdyz_kernel_null, "d2fdyz_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

    ELSE
!       =========================================================================

!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES

        rangexyz = (/1,nxsize,6,nysize-5,6,nzsize-5/)

    END IF

!   =========================================================================

END SUBROUTINE d2fdyz
