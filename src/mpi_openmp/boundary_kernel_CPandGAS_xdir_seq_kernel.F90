! Auto-generated at 2024-04-21 01:07:42.385851 by ops-translator

MODULE BOUNDARY_KERNEL_CPANDGAS_XDIR_MODULE

    USE OPS_FORTRAN_DECLARATIONS
    USE OPS_FORTRAN_RT_SUPPORT

    USE OPS_CONSTANTS
    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    INTEGER(KIND=4) :: xdim1
    INTEGER(KIND=4) :: ydim1
    INTEGER(KIND=4) :: zdim1
#define OPS_ACC1(x,y,z) ((x) + (xdim1*(y)) + (xdim1*ydim1*(z)) + 1)

    INTEGER(KIND=4) :: xdim2
    INTEGER(KIND=4) :: ydim2
    INTEGER(KIND=4) :: zdim2
#define OPS_ACC2(x,y,z) ((x) + (xdim2*(y)) + (xdim2*ydim2*(z)) + 1)

    INTEGER(KIND=4) :: xdim3
    INTEGER(KIND=4) :: ydim3
    INTEGER(KIND=4) :: zdim3
#define OPS_ACC3(x,y,z) ((x) + (xdim3*(y)) + (xdim3*ydim3*(z)) + 1)

    INTEGER(KIND=4) :: xdim4
    INTEGER(KIND=4) :: ydim4
    INTEGER(KIND=4) :: zdim4
#define OPS_ACC4(x,y,z) ((x) + (xdim4*(y)) + (xdim4*ydim4*(z)) + 1)

    INTEGER(KIND=4) :: xdim5
    INTEGER(KIND=4) :: ydim5
    INTEGER(KIND=4) :: zdim5
#define OPS_ACC5(x,y,z) ((x) + (xdim5*(y)) + (xdim5*ydim5*(z)) + 1)

    CONTAINS

!   =============
!   User function
!   =============

!DEC$ ATTRIBUTES FORCEINLINE :: boundary_kernel_CPandGAS_xdir
SUBROUTINE boundary_kernel_CPandGAS_xdir(transp,store7,drhs,strgx,strrx)

    real(kind=8), dimension(1), intent(in) :: transp, store7, drhs
    real(kind=8), dimension(1) :: strgx,strrx

    strgx(OPS_ACC4(0,0,0)) = transp(OPS_ACC1(0,0,0))
    strrx(OPS_ACC5(0,0,0)) = store7(OPS_ACC2(0,0,0)) / drhs(OPS_ACC3(0,0,0))

END SUBROUTINE

#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5

SUBROUTINE boundary_kernel_CPandGAS_xdir_wrap( &
    opsDat1Local, &
    opsDat2Local, &
    opsDat3Local, &
    opsDat4Local, &
    opsDat5Local, &
    dat1_base, &
    dat2_base, &
    dat3_base, &
    dat4_base, &
    dat5_base, &
    start_indx, &
    end_indx )

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat1Local
    INTEGER(KIND=4), INTENT(IN) :: dat1_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat2Local
    INTEGER(KIND=4), INTENT(IN) :: dat2_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat3Local
    INTEGER(KIND=4), INTENT(IN) :: dat3_base

    REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: opsDat4Local
    INTEGER(KIND=4), INTENT(IN) :: dat4_base

    REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: opsDat5Local
    INTEGER(KIND=4), INTENT(IN) :: dat5_base

    INTEGER(KIND=4), DIMENSION(3), INTENT(IN) :: start_indx, end_indx

    INTEGER(KIND=4) :: n_x, n_y, n_z

    !$OMP PARALLEL DO PRIVATE(n_x,n_y,n_z)
    DO n_z = 1, end_indx(3)-start_indx(3)+1
        DO n_y = 1, end_indx(2)-start_indx(2)+1
            !$OMP SIMD
            DO n_x = 1, end_indx(1)-start_indx(1)+1

                CALL boundary_kernel_CPandGAS_xdir( &

                opsDat1Local(dat1_base + ((n_x-1)*1) + ((n_y-1)*xdim1*1) + ((n_z-1)*ydim1*xdim1*1)), &

                opsDat2Local(dat2_base + ((n_x-1)*1) + ((n_y-1)*xdim2*1) + ((n_z-1)*ydim2*xdim2*1)), &

                opsDat3Local(dat3_base + ((n_x-1)*1) + ((n_y-1)*xdim3*1) + ((n_z-1)*ydim3*xdim3*1)), &

                opsDat4Local(dat4_base + ((n_x-1)*0) + ((n_y-1)*xdim4*1) + ((n_z-1)*ydim4*xdim4*1)), &

                opsDat5Local(dat5_base + ((n_x-1)*0) + ((n_y-1)*xdim5*1) + ((n_z-1)*ydim5*xdim5*1)) &
               )

            END DO
            !$OMP END SIMD
        END DO
    END DO

END SUBROUTINE

!   ===============
!   Host subroutine
!   ===============
#ifndef OPS_LAZY
SUBROUTINE boundary_kernel_CPandGAS_xdir_host( userSubroutine, block, dim, range, &
    opsArg1, &
    opsArg2, &
    opsArg3, &
    opsArg4, &
    opsArg5 &
    )

    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: userSubroutine
    TYPE(ops_block), INTENT(IN) :: block
    INTEGER(KIND=4), INTENT(IN) :: dim
    INTEGER(KIND=4), DIMENSION(2*dim), INTENT(IN) :: range

    TYPE(ops_arg), INTENT(IN) :: opsArg1
    TYPE(ops_arg), INTENT(IN) :: opsArg2
    TYPE(ops_arg), INTENT(IN) :: opsArg3
    TYPE(ops_arg), INTENT(IN) :: opsArg4
    TYPE(ops_arg), INTENT(IN) :: opsArg5

    TYPE(ops_arg), DIMENSION(5) :: opsArgArray

#else
SUBROUTINE boundary_kernel_CPandGAS_xdir_host_execute( descPtr )

    TYPE(ops_kernel_descriptor), INTENT(IN) :: descPtr
    TYPE(ops_block) :: block
    INTEGER(KIND=C_INT) :: dim
    INTEGER(KIND=C_INT), POINTER, DIMENSION(:) :: range
    CHARACTER(KIND=C_CHAR), POINTER, DIMENSION(:) :: userSubroutine
    TYPE(ops_arg), POINTER, DIMENSION(:) :: opsArgArray

    TYPE(ops_arg) :: opsArg1
    TYPE(ops_arg) :: opsArg2
    TYPE(ops_arg) :: opsArg3
    TYPE(ops_arg) :: opsArg4
    TYPE(ops_arg) :: opsArg5

#endif

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat1Local
    INTEGER(KIND=4) :: opsDat1Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat1_size
    INTEGER(KIND=4) :: dat1_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat2Local
    INTEGER(KIND=4) :: opsDat2Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat2_size
    INTEGER(KIND=4) :: dat2_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat3Local
    INTEGER(KIND=4) :: opsDat3Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat3_size
    INTEGER(KIND=4) :: dat3_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat4Local
    INTEGER(KIND=4) :: opsDat4Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat4_size
    INTEGER(KIND=4) :: dat4_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat5Local
    INTEGER(KIND=4) :: opsDat5Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat5_size
    INTEGER(KIND=4) :: dat5_base

    REAL(KIND=8) :: t1__, t2__, t3__
    REAL(KIND=4) :: transfer_total, transfer

    INTEGER(KIND=4), DIMENSION(3) :: start_indx, end_indx
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: kernelName

    kernelName = "boundary_kernel_CPandGAS_xdir"

#ifdef OPS_LAZY
!   ==========================
!   Set from kernel descriptor
!   ==========================
    dim = descPtr%dim
    CALL c_f_pointer(descPtr%range, range, (/2*dim/))
    CALL c_f_pointer(descPtr%name, userSubroutine, (/descPtr%name_len/))
    block%blockCptr = descPtr%block
    CALL c_f_pointer(block%blockCptr, block%blockPtr)
    CALL c_f_pointer(descPtr%args, opsArgArray, (/descPtr%nargs/))

    opsArg1 = opsArgArray(1)
    opsArg2 = opsArgArray(2)
    opsArg3 = opsArgArray(3)
    opsArg4 = opsArgArray(4)
    opsArg5 = opsArgArray(5)
#else
    opsArgArray(1) = opsArg1
    opsArgArray(2) = opsArg2
    opsArgArray(3) = opsArg3
    opsArgArray(4) = opsArg4
    opsArgArray(5) = opsArg5
#endif

    CALL setKernelTime(193, kernelName//c_null_char, 0.0_8, 0.0_8, 0.0_4, 1)
    CALL ops_timers_core(t1__)

#if defined(OPS_MPI) && !defined(OPS_LAZY)
    IF ( getRange(block, start_indx, end_indx, range) < 0 ) THEN
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

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg1), dat1_size, (/dim/))
    xdim1 = dat1_size(1)
    ydim1 = dat1_size(2)
    zdim1 = dat1_size(3)
    opsDat1Cardinality = opsArg1%dim * xdim1 * ydim1 * zdim1
    dat1_base = getDatBaseFromOpsArg3D(opsArg1, start_indx, 1)
    CALL c_f_pointer(opsArg1%data, opsDat1Local, (/opsDat1Cardinality/))

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg2), dat2_size, (/dim/))
    xdim2 = dat2_size(1)
    ydim2 = dat2_size(2)
    zdim2 = dat2_size(3)
    opsDat2Cardinality = opsArg2%dim * xdim2 * ydim2 * zdim2
    dat2_base = getDatBaseFromOpsArg3D(opsArg2, start_indx, 1)
    CALL c_f_pointer(opsArg2%data, opsDat2Local, (/opsDat2Cardinality/))

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg3), dat3_size, (/dim/))
    xdim3 = dat3_size(1)
    ydim3 = dat3_size(2)
    zdim3 = dat3_size(3)
    opsDat3Cardinality = opsArg3%dim * xdim3 * ydim3 * zdim3
    dat3_base = getDatBaseFromOpsArg3D(opsArg3, start_indx, 1)
    CALL c_f_pointer(opsArg3%data, opsDat3Local, (/opsDat3Cardinality/))

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg4), dat4_size, (/dim/))
    xdim4 = dat4_size(1)
    ydim4 = dat4_size(2)
    zdim4 = dat4_size(3)
    opsDat4Cardinality = opsArg4%dim * xdim4 * ydim4 * zdim4
    dat4_base = getDatBaseFromOpsArg3D(opsArg4, start_indx, 1)
    CALL c_f_pointer(opsArg4%data, opsDat4Local, (/opsDat4Cardinality/))

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg5), dat5_size, (/dim/))
    xdim5 = dat5_size(1)
    ydim5 = dat5_size(2)
    zdim5 = dat5_size(3)
    opsDat5Cardinality = opsArg5%dim * xdim5 * ydim5 * zdim5
    dat5_base = getDatBaseFromOpsArg3D(opsArg5, start_indx, 1)
    CALL c_f_pointer(opsArg5%data, opsDat5Local, (/opsDat5Cardinality/))

!   ==============
!   Halo exchanges
!   ==============
#ifndef OPS_LAZY
    CALL ops_H_D_exchanges_host(opsArgArray, 5)
    CALL ops_halo_exchanges(opsArgArray, 5, range)
    CALL ops_H_D_exchanges_host(opsArgArray, 5)
#endif

    CALL ops_timers_core(t2__)

!   ==============================
!   Call kernel wrapper subroutine
!   ==============================
    CALL boundary_kernel_CPandGAS_xdir_wrap( &
                        opsDat1Local, &
                        opsDat2Local, &
                        opsDat3Local, &
                        opsDat4Local, &
                        opsDat5Local, &
                        dat1_base, &
                        dat2_base, &
                        dat3_base, &
                        dat4_base, &
                        dat5_base, &
                        start_indx, &
                        end_indx )

    CALL ops_timers_core(t3__)

#ifndef OPS_LAZY
    CALL ops_set_dirtybit_host(opsArgArray, 5)
    CALL ops_set_halo_dirtybit3(opsArg4, range)
    CALL ops_set_halo_dirtybit3(opsArg5, range)
#endif

!   ========================
!   Timing and data movement
!   ========================
    transfer_total = 0.0_4
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg1, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg2, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg3, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg4, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg5, transfer)
    transfer_total = transfer_total + transfer

    CALL setKernelTime(193, kernelName//c_null_char, t3__-t2__, t2__-t1__, transfer_total, 0)

END SUBROUTINE

#ifdef OPS_LAZY
SUBROUTINE boundary_kernel_CPandGAS_xdir_host( userSubroutine, block, dim, range, &
    opsArg1, &
    opsArg2, &
    opsArg3, &
    opsArg4, &
    opsArg5 &
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

    TYPE(ops_arg), DIMENSION(5), TARGET :: opsArgArray
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: namelit

    namelit = "boundary_kernel_CPandGAS_xdir"

    opsArgArray(1) = opsArg1
    opsArgArray(2) = opsArg2
    opsArgArray(3) = opsArg3
    opsArgArray(4) = opsArg4
    opsArgArray(5) = opsArg5

    DO n_indx = 1, 3
        range_tmp(2*n_indx-1) = range(2*n_indx-1)-1
        range_tmp(2*n_indx)   = range(2*n_indx)
    END DO

    CALL create_kerneldesc_and_enque(namelit//c_null_char, c_loc(opsArgArray), &
                                    5, 193, dim, 0, c_loc(range_tmp), &
                                    block%blockCptr, c_funloc(boundary_kernel_CPandGAS_xdir_host_execute))

END SUBROUTINE
#endif

END MODULE BOUNDARY_KERNEL_CPANDGAS_XDIR_MODULE
