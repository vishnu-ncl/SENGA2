SUBROUTINE d2fdxz(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   D2FDXZ
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   15-MAY-2003:  RSC MODIFIED FOR SENGA2
!   10-OCT-2004:  RSC NULL VERSION

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES SECOND XZ-DERIVATIVE OF SPECIFIED FUNCTION

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
        rangexyz = (/1,nxsize,1,nysize,1,nzsize/)
        call ops_par_loop(d2fdxz_kernel_null, "d2fdxz_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
    ELSE
!       =========================================================================

!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz(1) = 1
        rangexyz(2) = nxsize
        rangexyz(5) = 1
        rangexyz(6) = nzsize

        IF(nendxl == nbound)    rangexyz(1) = 6
        IF(nendxr == nbound)    rangexyz(2) = nxsize-5
        IF(nendzl == nbound)    rangexyz(5) = 6
        IF(nendzr == nbound)    rangexyz(6) = nzsize-5

        rangexyz(3) = 1
        rangexyz(4) = nysize



    END IF

!   =========================================================================

END SUBROUTINE d2fdxz
