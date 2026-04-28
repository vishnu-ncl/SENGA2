! Auto-generated at 2026-04-28 18:44:33.111347 by ops-translator

MODULE COPY_KERNEL_MDIM_TO_SDIM_MODULE

    USE OPS_FORTRAN_DECLARATIONS
    USE OPS_FORTRAN_RT_SUPPORT

    USE OPS_CONSTANTS
    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

!$OMP DECLARE TARGET(xdim1_copy_kernel_mdim_to_sdim)
    INTEGER(KIND=4) :: xdim1_copy_kernel_mdim_to_sdim
    INTEGER(KIND=4) :: xdim1_copy_kernel_mdim_to_sdim_h = -1
!$OMP DECLARE TARGET(ydim1_copy_kernel_mdim_to_sdim)
    INTEGER(KIND=4) :: ydim1_copy_kernel_mdim_to_sdim
    INTEGER(KIND=4) :: ydim1_copy_kernel_mdim_to_sdim_h = -1
#define OPS_ACC1(x,y,z) ((x) + (xdim1_copy_kernel_mdim_to_sdim*(y)) + (xdim1_copy_kernel_mdim_to_sdim*ydim1_copy_kernel_mdim_to_sdim*(z)) + 1)

!$OMP DECLARE TARGET(multi_d2)
    INTEGER(KIND=4) :: multi_d2
!$OMP DECLARE TARGET(xdim2_copy_kernel_mdim_to_sdim)
    INTEGER(KIND=4) :: xdim2_copy_kernel_mdim_to_sdim
    INTEGER(KIND=4) :: xdim2_copy_kernel_mdim_to_sdim_h = -1
!$OMP DECLARE TARGET(ydim2_copy_kernel_mdim_to_sdim)
    INTEGER(KIND=4) :: ydim2_copy_kernel_mdim_to_sdim
    INTEGER(KIND=4) :: ydim2_copy_kernel_mdim_to_sdim_h = -1
!$OMP DECLARE TARGET(zdim2_copy_kernel_mdim_to_sdim)
    INTEGER(KIND=4) :: zdim2_copy_kernel_mdim_to_sdim
    INTEGER(KIND=4) :: zdim2_copy_kernel_mdim_to_sdim_h = -1
#define OPS_ACC_MD2(d,x,y,z) ((x) + (xdim2_copy_kernel_mdim_to_sdim*(y)) + (xdim2_copy_kernel_mdim_to_sdim*ydim2_copy_kernel_mdim_to_sdim*(z)) + (xdim2_copy_kernel_mdim_to_sdim*ydim2_copy_kernel_mdim_to_sdim*zdim2_copy_kernel_mdim_to_sdim*(d-1)) + 1)

    CONTAINS

!   =============
!   User function
!   =============

SUBROUTINE copy_kernel_mdim_to_sdim(out_arr, in_arr, ispec)

    real(kind=8), dimension(1) :: out_arr
    real(kind=8), dimension(1), intent(in) :: in_arr
    integer(kind=4), intent(in) :: ispec

    out_arr(OPS_ACC1(0,0,0)) = in_arr(OPS_ACC_MD2(ispec,0,0,0))

END SUBROUTINE

#undef OPS_ACC1
#undef OPS_ACC_MD2

SUBROUTINE copy_kernel_mdim_to_sdim_wrap( &
    opsDat1Local, &
    opsDat2Local, &
    opsGblDat3Device, &
    dat1_base, &
    dat2_base, &
    dat3_base, &
    dat1_dim, &
    dat2_dim, &
    dat3_dim, &
    start_indx, &
    end_indx )

    INTEGER(KIND=4), VALUE :: dat1_dim
    REAL(KIND=8), DIMENSION(dat1_dim), INTENT(OUT) :: opsDat1Local
    INTEGER(KIND=4), VALUE :: dat1_base

    INTEGER(KIND=4), VALUE :: dat2_dim
    REAL(KIND=8), DIMENSION(dat2_dim), INTENT(IN) :: opsDat2Local
    INTEGER(KIND=4), VALUE :: dat2_base

    INTEGER(KIND=4), VALUE :: dat3_dim
    INTEGER(KIND=4), VALUE ::  opsGblDat3Device
    INTEGER(KIND=4), INTENT(IN) :: dat3_base

    INTEGER(KIND=4), DIMENSION(3), INTENT(IN) :: start_indx, end_indx

    INTEGER(KIND=4) :: n_x, n_y, n_z

    INTEGER(KIND=4) :: start_indx_1, end_indx_1
    INTEGER(KIND=4) :: start_indx_2, end_indx_2
    INTEGER(KIND=4) :: start_indx_3, end_indx_3

    start_indx_1 = start_indx(1)
    end_indx_1  = end_indx(1)
    start_indx_2 = start_indx(2)
    end_indx_2  = end_indx(2)
    start_indx_3 = start_indx(3)
    end_indx_3  = end_indx(3)

#ifdef _CRAYFTN
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(n_x,n_y,n_z)
#else
!$OMP TARGET TEAMS LOOP COLLAPSE(3) PRIVATE(n_x,n_y,n_z)
#endif
    DO n_z = 1, end_indx_3-start_indx_3+1
        DO n_y = 1, end_indx_2-start_indx_2+1
            DO n_x = 1, end_indx_1-start_indx_1+1

                CALL copy_kernel_mdim_to_sdim( &
                opsDat1Local(dat1_base + ((n_x-1)*1) + ((n_y-1)*xdim1_copy_kernel_mdim_to_sdim*1) + ((n_z-1)*ydim1_copy_kernel_mdim_to_sdim*xdim1_copy_kernel_mdim_to_sdim*1)), &
                opsDat2Local(dat2_base + ((n_x-1)*1) + ((n_y-1)*xdim2_copy_kernel_mdim_to_sdim*1) + ((n_z-1)*ydim2_copy_kernel_mdim_to_sdim*xdim2_copy_kernel_mdim_to_sdim*1)), &
                opsGblDat3Device &
               )

            END DO
        END DO
    END DO
#ifdef _CRAYFTN
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
!$OMP END TARGET TEAMS LOOP
#endif

END SUBROUTINE

!   ===============
!   Host subroutine
!   ===============
#ifndef OPS_LAZY
SUBROUTINE copy_kernel_mdim_to_sdim_host( userSubroutine, block, dim, range, &
    opsArg1, &
    opsArg2, &
    opsArg3 &
    )

    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: userSubroutine
    TYPE(ops_block), INTENT(IN) :: block
    INTEGER(KIND=4), INTENT(IN) :: dim
    INTEGER(KIND=4), DIMENSION(2*dim), INTENT(IN) :: range

    TYPE(ops_arg), INTENT(IN) :: opsArg1
    TYPE(ops_arg), INTENT(IN) :: opsArg2
    TYPE(ops_arg), INTENT(IN) :: opsArg3

    TYPE(ops_arg), DIMENSION(3) :: opsArgArray

#else
SUBROUTINE copy_kernel_mdim_to_sdim_host_execute( descPtr )

    TYPE(ops_kernel_descriptor), INTENT(IN) :: descPtr
    TYPE(ops_block) :: block
    INTEGER(KIND=C_INT) :: dim
    INTEGER(KIND=C_INT), POINTER, DIMENSION(:) :: range
    CHARACTER(KIND=C_CHAR), POINTER, DIMENSION(:) :: userSubroutine
    TYPE(ops_arg), POINTER, DIMENSION(:) :: opsArgArray

    TYPE(ops_arg) :: opsArg1
    TYPE(ops_arg) :: opsArg2
    TYPE(ops_arg) :: opsArg3

#endif

    REAL(KIND=8), DIMENSION(:), POINTER :: opsDat1Local
    INTEGER(KIND=4) :: opsDat1Cardinality
    INTEGER(KIND=4), DIMENSION(:), POINTER  :: dat1_size
    INTEGER(KIND=4) :: dat1_base
    INTEGER(KIND=4) :: xdim1, ydim1, zdim1

    REAL(KIND=8), DIMENSION(:), POINTER :: opsDat2Local
    INTEGER(KIND=4) :: opsDat2Cardinality
    INTEGER(KIND=4), DIMENSION(:), POINTER  :: dat2_size
    INTEGER(KIND=4) :: dat2_base
    INTEGER(KIND=4) :: xdim2, ydim2, zdim2

    INTEGER(KIND=4), DIMENSION(:), POINTER :: opsDat3Local
    INTEGER(KIND=4) :: opsDat3Cardinality
    INTEGER(KIND=4) :: dat3_base

    REAL(KIND=8) :: t1__, t2__, t3__
    REAL(KIND=4) :: transfer_total, transfer

    INTEGER(KIND=4) :: x_size, y_size, z_size

    INTEGER(KIND=4), DIMENSION(3) :: start_indx, end_indx
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: kernelName

    kernelName = "copy_kernel_mdim_to_sdim"

#ifdef OPS_LAZY
!   ==========================
!   Set from kernel descriptor
!   ==========================
    dim = descPtr%dim
    CALL c_f_pointer(descPtr%range, range, [2*dim])
    CALL c_f_pointer(descPtr%name, userSubroutine, [descPtr%name_len])
    block%blockCptr = descPtr%block
    CALL c_f_pointer(block%blockCptr, block%blockPtr)
    CALL c_f_pointer(descPtr%args, opsArgArray, [descPtr%nargs])

    opsArg1 = opsArgArray(1)
    opsArg2 = opsArgArray(2)
    opsArg3 = opsArgArray(3)
#else
    opsArgArray(1) = opsArg1
    opsArgArray(2) = opsArg2
    opsArgArray(3) = opsArg3
#endif

    CALL setKernelTime(566, kernelName//c_null_char, 0.0_8, 0.0_8, 0.0_4, 1)
    CALL ops_timers_core(t1__)

#if defined(OPS_MPI) && !defined(OPS_LAZY)
    IF ( getRange(block, start_indx, end_indx, range, 3) < 0 ) THEN
        RETURN
    END IF
#elif !defined(OPS_MPI)  && !defined(OPS_LAZY)
    DO n_indx = 1, 3
        start_indx(n_indx) = range(2*n_indx-1)
        end_indx  (n_indx) = range(2*n_indx)
    END DO
#else
    DO n_indx = 1, 3
        start_indx(n_indx) = range(2*n_indx-1) + 1
        end_indx  (n_indx) = range(2*n_indx)
    END DO
#endif

    x_size = MAX(0, end_indx(1)-start_indx(1)+1)
    y_size = MAX(0, end_indx(2)-start_indx(2)+1)
    z_size = MAX(0, end_indx(3)-start_indx(3)+1)

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg1), dat1_size, [dim])
    xdim1 = dat1_size(1)
    ydim1 = dat1_size(2)
    zdim1 = dat1_size(3)
    opsDat1Cardinality = opsArg1%dim * xdim1 * ydim1 * zdim1
    dat1_base = getDatBaseFromOpsArg3D(opsArg1, start_indx, 1)
    CALL c_f_pointer(opsArg1%data_d, opsDat1Local, [opsDat1Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg2), dat2_size, [dim])
    xdim2 = dat2_size(1)
    ydim2 = dat2_size(2)
    zdim2 = dat2_size(3)
    opsDat2Cardinality = opsArg2%dim * xdim2 * ydim2 * zdim2
    multi_d2 = getDatDimFromOpsArg(opsArg2) ! dimension of dat
!$OMP TARGET UPDATE TO(multi_d2)
    dat2_base = getDatBaseFromOpsArg3D(opsArg2, start_indx, multi_d2)
    CALL c_f_pointer(opsArg2%data_d, opsDat2Local, [opsDat2Cardinality])

    opsDat3Cardinality = opsArg3%dim
    CALL c_f_pointer(opsArg3%data, opsDat3Local, [opsDat3Cardinality])
    dat3_base = 1

    IF (&
         (xdim1 .NE. xdim1_copy_kernel_mdim_to_sdim_h) .OR. &
         (ydim1 .NE. ydim1_copy_kernel_mdim_to_sdim_h) .OR. &
         (xdim2 .NE. xdim2_copy_kernel_mdim_to_sdim_h) .OR. &
         (ydim2 .NE. ydim2_copy_kernel_mdim_to_sdim_h) .OR. &
         (zdim2 .NE. zdim2_copy_kernel_mdim_to_sdim_h) &
       ) THEN
            xdim1_copy_kernel_mdim_to_sdim = xdim1
            xdim1_copy_kernel_mdim_to_sdim_h = xdim1
!$OMP TARGET UPDATE TO(xdim1_copy_kernel_mdim_to_sdim)
            ydim1_copy_kernel_mdim_to_sdim = ydim1
            ydim1_copy_kernel_mdim_to_sdim_h = ydim1
!$OMP TARGET UPDATE TO(ydim1_copy_kernel_mdim_to_sdim)
            xdim2_copy_kernel_mdim_to_sdim = xdim2
            xdim2_copy_kernel_mdim_to_sdim_h = xdim2
!$OMP TARGET UPDATE TO(xdim2_copy_kernel_mdim_to_sdim)
            ydim2_copy_kernel_mdim_to_sdim = ydim2
            ydim2_copy_kernel_mdim_to_sdim_h = ydim2
            zdim2_copy_kernel_mdim_to_sdim = zdim2
            zdim2_copy_kernel_mdim_to_sdim_h = zdim2
!$OMP TARGET UPDATE TO(ydim2_copy_kernel_mdim_to_sdim,zdim2_copy_kernel_mdim_to_sdim)
    END IF

!   ==============
!   Halo exchanges
!   ==============
#ifndef OPS_LAZY
    CALL ops_H_D_exchanges_device(opsArgArray, 3)
    CALL ops_halo_exchanges(opsArgArray, 3, range, 3)
#endif

    CALL ops_timers_core(t2__)

!   ==============================
!   Call kernel wrapper subroutine
!   ==============================
    CALL copy_kernel_mdim_to_sdim_wrap( &
                        opsDat1Local, &
                        opsDat2Local, &
                        opsDat3Local(1), &
                        dat1_base, &
                        dat2_base, &
                        dat3_base, &
                        opsDat1Cardinality, &
                        opsDat2Cardinality, &
                        opsDat3Cardinality, &
                        start_indx, &
                        end_indx )

    CALL ops_timers_core(t3__)

#ifndef OPS_LAZY
    CALL ops_set_dirtybit_device(opsArgArray, 3)
    CALL ops_set_halo_dirtybit3(opsArg1, range, 3)
#endif

!   ========================
!   Timing and data movement
!   ========================
    transfer_total = 0.0_4
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg1, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg2, transfer)
    transfer_total = transfer_total + transfer

    CALL setKernelTime(566, kernelName//c_null_char, t3__-t2__, t2__-t1__, transfer_total, 0)

END SUBROUTINE

#ifdef OPS_LAZY
SUBROUTINE copy_kernel_mdim_to_sdim_host( userSubroutine, block, dim, range, &
    opsArg1, &
    opsArg2, &
    opsArg3 &
    )

    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN), TARGET :: userSubroutine
    TYPE(ops_block), INTENT(IN) :: block
    INTEGER(KIND=4), INTENT(IN) :: dim
    INTEGER(KIND=4), DIMENSION(2*dim), INTENT(INOUT), TARGET :: range
    INTEGER(KIND=4), DIMENSION(2*dim), TARGET :: range_tmp

    TYPE(ops_arg), INTENT(IN) :: opsArg1
    TYPE(ops_arg), INTENT(IN) :: opsArg2
    TYPE(ops_arg), INTENT(IN) :: opsArg3

    TYPE(ops_arg), DIMENSION(3), TARGET :: opsArgArray
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: namelit

    namelit = "copy_kernel_mdim_to_sdim"

    opsArgArray(1) = opsArg1
    opsArgArray(2) = opsArg2
    opsArgArray(3) = opsArg3

    DO n_indx = 1, 3
        range_tmp(2*n_indx-1) = range(2*n_indx-1)-1
        range_tmp(2*n_indx)   = range(2*n_indx)
    END DO

    CALL create_kerneldesc_and_enque(namelit//c_null_char, c_loc(opsArgArray), &
                                    3, 566, dim, 1, c_loc(range_tmp), &
                                    block%blockCptr, c_funloc(copy_kernel_mdim_to_sdim_host_execute))

END SUBROUTINE
#endif

END MODULE COPY_KERNEL_MDIM_TO_SDIM_MODULE
