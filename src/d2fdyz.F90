SUBROUTINE d2fdyz(functn,fderiv)

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
INTEGER :: jstart,jfinis,kstart,kfinis
INTEGER :: rangexyz(6)


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============
jstart = jstal
jfinis = jstol
kstart = kstal
kfinis = kstol
IF(nendyl == nbound)    jstart = jstap5
IF(nendyr == nbound)    jfinis = jstom5
IF(nendzl == nbound)    kstart = kstap5
IF(nendzr == nbound)    kfinis = kstom5


    if(nyglbl == 1) then
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(d2fdyz_kernel_null, "d2fdyz_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(dp)", OPS_WRITE))
    else

    end if

!     =========================================================================



END SUBROUTINE d2fdyz
