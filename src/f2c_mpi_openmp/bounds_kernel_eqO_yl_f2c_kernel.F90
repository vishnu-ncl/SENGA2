! Auto-generated at 2026-04-28 18:43:32.311689 by ops-translator

MODULE BOUNDS_KERNEL_EQO_YL_MODULE

    USE OPS_FORTRAN_DECLARATIONS
    USE OPS_FORTRAN_RT_SUPPORT

    USE OPS_CONSTANTS
    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    INTERFACE

SUBROUTINE bounds_kernel_eqO_yl_host_c(name, args, nargs, dim, range, block) BIND(C,name='bounds_kernel_eqO_yl_host_c')

    USE, INTRINSIC :: ISO_C_BINDING
    import :: ops_block_core, ops_arg

    character(kind=c_char,len=1) :: name(*)
    type(c_ptr), value           :: args
    integer(kind=c_int), value   :: nargs
    integer(kind=c_int), value   :: dim
    type(c_ptr), value      :: range
    type(c_ptr), value      :: block

END SUBROUTINE bounds_kernel_eqO_yl_host_c

    END INTERFACE

    CONTAINS

SUBROUTINE bounds_kernel_eqO_yl_host( userSubroutine, block, dim, range, &
    opsArg1 &
    )

    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN), TARGET :: userSubroutine
    TYPE(ops_block), INTENT(IN) :: block
    INTEGER(KIND=4), INTENT(IN) :: dim
    INTEGER(KIND=4), DIMENSION(2*dim), INTENT(INOUT), TARGET :: range
    INTEGER(KIND=4), DIMENSION(2*dim), TARGET :: range_tmp

    TYPE(ops_arg), INTENT(IN) :: opsArg1

    TYPE(ops_arg), DIMENSION(1), TARGET :: opsArgArray
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: namelit

    namelit = "bounds_kernel_eqO_yl"

    opsArgArray(1) = opsArg1

    DO n_indx = 1, 3
        range_tmp(2*n_indx-1) = range(2*n_indx-1)-1
        range_tmp(2*n_indx)   = range(2*n_indx)
    END DO

    CALL bounds_kernel_eqO_yl_host_c(namelit//c_null_char, c_loc(opsArgArray), &
                                    1, dim, c_loc(range_tmp), &
                                    block%blockCptr)

END SUBROUTINE bounds_kernel_eqO_yl_host

END MODULE BOUNDS_KERNEL_EQO_YL_MODULE
