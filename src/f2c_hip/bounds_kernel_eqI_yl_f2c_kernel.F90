! Auto-generated at 2026-04-28 18:44:22.089017 by ops-translator

MODULE BOUNDS_KERNEL_EQI_YL_MODULE

    USE OPS_FORTRAN_DECLARATIONS
    USE OPS_FORTRAN_RT_SUPPORT

    USE OPS_CONSTANTS
    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    INTERFACE

SUBROUTINE bounds_kernel_eqI_yl_host_c(name, args, nargs, dim, range, block) BIND(C,name='bounds_kernel_eqI_yl_host_c')

    USE, INTRINSIC :: ISO_C_BINDING
    import :: ops_block_core, ops_arg

    character(kind=c_char,len=1) :: name(*)
    type(c_ptr), value           :: args
    integer(kind=c_int), value   :: nargs
    integer(kind=c_int), value   :: dim
    type(c_ptr), value      :: range
    type(c_ptr), value      :: block

END SUBROUTINE bounds_kernel_eqI_yl_host_c

    END INTERFACE

    CONTAINS

SUBROUTINE bounds_kernel_eqI_yl_host( userSubroutine, block, dim, range, &
    opsArg1, &
    opsArg2, &
    opsArg3, &
    opsArg4, &
    opsArg5, &
    opsArg6, &
    opsArg7, &
    opsArg8, &
    opsArg9, &
    opsArg10, &
    opsArg11, &
    opsArg12, &
    opsArg13 &
    )

    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN), TARGET :: userSubroutine
    TYPE(ops_block), INTENT(IN) :: block
    INTEGER(KIND=4), INTENT(IN) :: dim
    INTEGER(KIND=4), DIMENSION(2*dim), INTENT(INOUT), TARGET :: range
    INTEGER(KIND=4), DIMENSION(2*dim), TARGET :: range_tmp

    TYPE(ops_arg), INTENT(IN) :: opsArg1
    TYPE(ops_arg), INTENT(IN) :: opsArg2
    TYPE(ops_arg), INTENT(IN) :: opsArg3
    TYPE(ops_arg), INTENT(IN) :: opsArg4
    TYPE(ops_arg), INTENT(IN) :: opsArg5
    TYPE(ops_arg), INTENT(IN) :: opsArg6
    TYPE(ops_arg), INTENT(IN) :: opsArg7
    TYPE(ops_arg), INTENT(IN) :: opsArg8
    TYPE(ops_arg), INTENT(IN) :: opsArg9
    TYPE(ops_arg), INTENT(IN) :: opsArg10
    TYPE(ops_arg), INTENT(IN) :: opsArg11
    TYPE(ops_arg), INTENT(IN) :: opsArg12
    TYPE(ops_arg), INTENT(IN) :: opsArg13

    TYPE(ops_arg), DIMENSION(13), TARGET :: opsArgArray
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: namelit

    namelit = "bounds_kernel_eqI_yl"

    opsArgArray(1) = opsArg1
    opsArgArray(2) = opsArg2
    opsArgArray(3) = opsArg3
    opsArgArray(4) = opsArg4
    opsArgArray(5) = opsArg5
    opsArgArray(6) = opsArg6
    opsArgArray(7) = opsArg7
    opsArgArray(8) = opsArg8
    opsArgArray(9) = opsArg9
    opsArgArray(10) = opsArg10
    opsArgArray(11) = opsArg11
    opsArgArray(12) = opsArg12
    opsArgArray(13) = opsArg13

    DO n_indx = 1, 3
        range_tmp(2*n_indx-1) = range(2*n_indx-1)-1
        range_tmp(2*n_indx)   = range(2*n_indx)
    END DO

    CALL bounds_kernel_eqI_yl_host_c(namelit//c_null_char, c_loc(opsArgArray), &
                                    13, dim, c_loc(range_tmp), &
                                    block%blockCptr)

END SUBROUTINE bounds_kernel_eqI_yl_host

END MODULE BOUNDS_KERNEL_EQI_YL_MODULE
