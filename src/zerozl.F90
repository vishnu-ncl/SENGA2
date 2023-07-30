SUBROUTINE zerozl(farray)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   ARGUMENTS
!   =========
    TYPE(ops_dat) :: farray

!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

    rangexyz = [1,nxglbl,1,nyglbl,1,1]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                        ops_arg_dat(farray, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!   =========================================================================

END SUBROUTINE zerozl
