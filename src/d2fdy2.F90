SUBROUTINE d2fdy2(functn,fderiv)
 
use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

    TYPE(ops_dat) :: functn, fderiv


!     ARGUMENTS
!     =========




!     LOCAL DATA
!     ==========
    INTEGER :: jstart,jfinis
    INTEGER :: rangexyz(6)

!     BEGIN
!     =====

!     =========================================================================
    if(nyglbl == 1) then
        rangexyz = (/istal,istol, jstal,jstol, kstal,kstol/)
        call ops_par_loop(d2fdy2_kernel_null, "d2fdy2_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
    else

    end if

!     =========================================================================


END SUBROUTINE d2fdy2
