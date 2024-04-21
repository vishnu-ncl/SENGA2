! Auto-generated at 2024-04-21 01:07:43.831834 by ops-translator

MODULE BOUNDS_KERNEL_EQQ_YR_MODULE

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

    INTEGER(KIND=4) :: xdim6
    INTEGER(KIND=4) :: ydim6
    INTEGER(KIND=4) :: zdim6
#define OPS_ACC6(x,y,z) ((x) + (xdim6*(y)) + (xdim6*ydim6*(z)) + 1)

    INTEGER(KIND=4) :: xdim7
    INTEGER(KIND=4) :: ydim7
    INTEGER(KIND=4) :: zdim7
#define OPS_ACC7(x,y,z) ((x) + (xdim7*(y)) + (xdim7*ydim7*(z)) + 1)

    INTEGER(KIND=4) :: xdim8
    INTEGER(KIND=4) :: ydim8
    INTEGER(KIND=4) :: zdim8
#define OPS_ACC8(x,y,z) ((x) + (xdim8*(y)) + (xdim8*ydim8*(z)) + 1)

    INTEGER(KIND=4) :: xdim9
    INTEGER(KIND=4) :: ydim9
    INTEGER(KIND=4) :: zdim9
#define OPS_ACC9(x,y,z) ((x) + (xdim9*(y)) + (xdim9*ydim9*(z)) + 1)

    INTEGER(KIND=4) :: xdim10
    INTEGER(KIND=4) :: ydim10
    INTEGER(KIND=4) :: zdim10
#define OPS_ACC10(x,y,z) ((x) + (xdim10*(y)) + (xdim10*ydim10*(z)) + 1)

    CONTAINS

!   =============
!   User function
!   =============

!DEC$ ATTRIBUTES FORCEINLINE :: bounds_kernel_eqQ_yr
SUBROUTINE bounds_kernel_eqQ_yr(erhs,yrhs,strdyr,strvyr,strhyr,stryyr,bclyyr,bcl1yr,bcl2yr,ova2yr,flag_pio_yr)

    real(kind=8), dimension(1) :: erhs,yrhs
    real(kind=8), dimension(1), intent(in) :: strdyr,strvyr,strhyr,stryyr,bclyyr,bcl1yr,bcl2yr,ova2yr
    integer(kind=4), intent(in) :: flag_pio_yr

    real(kind=8) :: fornow

    IF ( (strvyr(OPS_ACC4(0,0,0)) < 0.0_8) .AND. (flag_pio_yr==1) ) THEN

        fornow = bclyyr(OPS_ACC7(0,0,0))*strdyr(OPS_ACC3(0,0,0))

        erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0)) - fornow*strhyr(OPS_ACC5(0,0,0))

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - (bcl2yr(OPS_ACC9(0,0,0))+bcl1yr(OPS_ACC8(0,0,0))*ova2yr(OPS_ACC10(0,0,0)))*stryyr(OPS_ACC6(0,0,0)) &
                              - fornow

    ELSE

        yrhs(OPS_ACC2(0,0,0)) = yrhs(OPS_ACC2(0,0,0)) &
                              - bcl1yr(OPS_ACC8(0,0,0))*ova2yr(OPS_ACC10(0,0,0))*stryyr(OPS_ACC6(0,0,0))

    END IF

END SUBROUTINE

#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6
#undef OPS_ACC7
#undef OPS_ACC8
#undef OPS_ACC9
#undef OPS_ACC10

SUBROUTINE bounds_kernel_eqQ_yr_wrap( &
    opsDat1Local, &
    opsDat2Local, &
    opsDat3Local, &
    opsDat4Local, &
    opsDat5Local, &
    opsDat6Local, &
    opsDat7Local, &
    opsDat8Local, &
    opsDat9Local, &
    opsDat10Local, &
    opsDat11Local, &
    dat1_base, &
    dat2_base, &
    dat3_base, &
    dat4_base, &
    dat5_base, &
    dat6_base, &
    dat7_base, &
    dat8_base, &
    dat9_base, &
    dat10_base, &
    dat11_base, &
    start_indx, &
    end_indx )

    REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: opsDat1Local
    INTEGER(KIND=4), INTENT(IN) :: dat1_base

    REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: opsDat2Local
    INTEGER(KIND=4), INTENT(IN) :: dat2_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat3Local
    INTEGER(KIND=4), INTENT(IN) :: dat3_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat4Local
    INTEGER(KIND=4), INTENT(IN) :: dat4_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat5Local
    INTEGER(KIND=4), INTENT(IN) :: dat5_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat6Local
    INTEGER(KIND=4), INTENT(IN) :: dat6_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat7Local
    INTEGER(KIND=4), INTENT(IN) :: dat7_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat8Local
    INTEGER(KIND=4), INTENT(IN) :: dat8_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat9Local
    INTEGER(KIND=4), INTENT(IN) :: dat9_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat10Local
    INTEGER(KIND=4), INTENT(IN) :: dat10_base

    INTEGER(KIND=4), DIMENSION(*), INTENT(IN) :: opsDat11Local
    INTEGER(KIND=4), INTENT(IN) :: dat11_base

    INTEGER(KIND=4), DIMENSION(3), INTENT(IN) :: start_indx, end_indx

    INTEGER(KIND=4) :: n_x, n_y, n_z

    !$OMP PARALLEL DO PRIVATE(n_x,n_y,n_z)
    DO n_z = 1, end_indx(3)-start_indx(3)+1
        DO n_y = 1, end_indx(2)-start_indx(2)+1
            !$OMP SIMD
            DO n_x = 1, end_indx(1)-start_indx(1)+1

                CALL bounds_kernel_eqQ_yr( &

                opsDat1Local(dat1_base + ((n_x-1)*1) + ((n_y-1)*xdim1*1) + ((n_z-1)*ydim1*xdim1*1)), &

                opsDat2Local(dat2_base + ((n_x-1)*1) + ((n_y-1)*xdim2*1) + ((n_z-1)*ydim2*xdim2*1)), &

                opsDat3Local(dat3_base + ((n_x-1)*1) + ((n_y-1)*xdim3*0) + ((n_z-1)*ydim3*xdim3*1)), &

                opsDat4Local(dat4_base + ((n_x-1)*1) + ((n_y-1)*xdim4*0) + ((n_z-1)*ydim4*xdim4*1)), &

                opsDat5Local(dat5_base + ((n_x-1)*1) + ((n_y-1)*xdim5*0) + ((n_z-1)*ydim5*xdim5*1)), &

                opsDat6Local(dat6_base + ((n_x-1)*1) + ((n_y-1)*xdim6*0) + ((n_z-1)*ydim6*xdim6*1)), &

                opsDat7Local(dat7_base + ((n_x-1)*1) + ((n_y-1)*xdim7*0) + ((n_z-1)*ydim7*xdim7*1)), &

                opsDat8Local(dat8_base + ((n_x-1)*1) + ((n_y-1)*xdim8*0) + ((n_z-1)*ydim8*xdim8*1)), &

                opsDat9Local(dat9_base + ((n_x-1)*1) + ((n_y-1)*xdim9*0) + ((n_z-1)*ydim9*xdim9*1)), &

                opsDat10Local(dat10_base + ((n_x-1)*1) + ((n_y-1)*xdim10*0) + ((n_z-1)*ydim10*xdim10*1)), &
                opsDat11Local(dat11_base) &
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
SUBROUTINE bounds_kernel_eqQ_yr_host( userSubroutine, block, dim, range, &
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
    opsArg11 &
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
    TYPE(ops_arg), INTENT(IN) :: opsArg6
    TYPE(ops_arg), INTENT(IN) :: opsArg7
    TYPE(ops_arg), INTENT(IN) :: opsArg8
    TYPE(ops_arg), INTENT(IN) :: opsArg9
    TYPE(ops_arg), INTENT(IN) :: opsArg10
    TYPE(ops_arg), INTENT(IN) :: opsArg11

    TYPE(ops_arg), DIMENSION(11) :: opsArgArray

#else
SUBROUTINE bounds_kernel_eqQ_yr_host_execute( descPtr )

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
    TYPE(ops_arg) :: opsArg6
    TYPE(ops_arg) :: opsArg7
    TYPE(ops_arg) :: opsArg8
    TYPE(ops_arg) :: opsArg9
    TYPE(ops_arg) :: opsArg10
    TYPE(ops_arg) :: opsArg11

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

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat6Local
    INTEGER(KIND=4) :: opsDat6Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat6_size
    INTEGER(KIND=4) :: dat6_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat7Local
    INTEGER(KIND=4) :: opsDat7Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat7_size
    INTEGER(KIND=4) :: dat7_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat8Local
    INTEGER(KIND=4) :: opsDat8Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat8_size
    INTEGER(KIND=4) :: dat8_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat9Local
    INTEGER(KIND=4) :: opsDat9Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat9_size
    INTEGER(KIND=4) :: dat9_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat10Local
    INTEGER(KIND=4) :: opsDat10Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat10_size
    INTEGER(KIND=4) :: dat10_base

    INTEGER(KIND=4), POINTER, DIMENSION(:) :: opsDat11Local
    INTEGER(KIND=4) :: dat11_base

    REAL(KIND=8) :: t1__, t2__, t3__
    REAL(KIND=4) :: transfer_total, transfer

    INTEGER(KIND=4), DIMENSION(3) :: start_indx, end_indx
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: kernelName

    kernelName = "bounds_kernel_eqQ_yr"

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
    opsArg6 = opsArgArray(6)
    opsArg7 = opsArgArray(7)
    opsArg8 = opsArgArray(8)
    opsArg9 = opsArgArray(9)
    opsArg10 = opsArgArray(10)
    opsArg11 = opsArgArray(11)
#else
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
#endif

    CALL setKernelTime(411, kernelName//c_null_char, 0.0_8, 0.0_8, 0.0_4, 1)
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

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg6), dat6_size, (/dim/))
    xdim6 = dat6_size(1)
    ydim6 = dat6_size(2)
    zdim6 = dat6_size(3)
    opsDat6Cardinality = opsArg6%dim * xdim6 * ydim6 * zdim6
    dat6_base = getDatBaseFromOpsArg3D(opsArg6, start_indx, 1)
    CALL c_f_pointer(opsArg6%data, opsDat6Local, (/opsDat6Cardinality/))

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg7), dat7_size, (/dim/))
    xdim7 = dat7_size(1)
    ydim7 = dat7_size(2)
    zdim7 = dat7_size(3)
    opsDat7Cardinality = opsArg7%dim * xdim7 * ydim7 * zdim7
    dat7_base = getDatBaseFromOpsArg3D(opsArg7, start_indx, 1)
    CALL c_f_pointer(opsArg7%data, opsDat7Local, (/opsDat7Cardinality/))

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg8), dat8_size, (/dim/))
    xdim8 = dat8_size(1)
    ydim8 = dat8_size(2)
    zdim8 = dat8_size(3)
    opsDat8Cardinality = opsArg8%dim * xdim8 * ydim8 * zdim8
    dat8_base = getDatBaseFromOpsArg3D(opsArg8, start_indx, 1)
    CALL c_f_pointer(opsArg8%data, opsDat8Local, (/opsDat8Cardinality/))

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg9), dat9_size, (/dim/))
    xdim9 = dat9_size(1)
    ydim9 = dat9_size(2)
    zdim9 = dat9_size(3)
    opsDat9Cardinality = opsArg9%dim * xdim9 * ydim9 * zdim9
    dat9_base = getDatBaseFromOpsArg3D(opsArg9, start_indx, 1)
    CALL c_f_pointer(opsArg9%data, opsDat9Local, (/opsDat9Cardinality/))

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg10), dat10_size, (/dim/))
    xdim10 = dat10_size(1)
    ydim10 = dat10_size(2)
    zdim10 = dat10_size(3)
    opsDat10Cardinality = opsArg10%dim * xdim10 * ydim10 * zdim10
    dat10_base = getDatBaseFromOpsArg3D(opsArg10, start_indx, 1)
    CALL c_f_pointer(opsArg10%data, opsDat10Local, (/opsDat10Cardinality/))

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg11), opsDat11Local, (/opsArg11%dim/))
    dat11_base = 1

!   ==============
!   Halo exchanges
!   ==============
#ifndef OPS_LAZY
    CALL ops_H_D_exchanges_host(opsArgArray, 11)
    CALL ops_halo_exchanges(opsArgArray, 11, range)
    CALL ops_H_D_exchanges_host(opsArgArray, 11)
#endif

    CALL ops_timers_core(t2__)

!   ==============================
!   Call kernel wrapper subroutine
!   ==============================
    CALL bounds_kernel_eqQ_yr_wrap( &
                        opsDat1Local, &
                        opsDat2Local, &
                        opsDat3Local, &
                        opsDat4Local, &
                        opsDat5Local, &
                        opsDat6Local, &
                        opsDat7Local, &
                        opsDat8Local, &
                        opsDat9Local, &
                        opsDat10Local, &
                        opsDat11Local, &
                        dat1_base, &
                        dat2_base, &
                        dat3_base, &
                        dat4_base, &
                        dat5_base, &
                        dat6_base, &
                        dat7_base, &
                        dat8_base, &
                        dat9_base, &
                        dat10_base, &
                        dat11_base, &
                        start_indx, &
                        end_indx )

    CALL ops_timers_core(t3__)

#ifndef OPS_LAZY
    CALL ops_set_dirtybit_host(opsArgArray, 11)
    CALL ops_set_halo_dirtybit3(opsArg1, range)
    CALL ops_set_halo_dirtybit3(opsArg2, range)
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
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg6, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg7, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg8, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg9, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg10, transfer)
    transfer_total = transfer_total + transfer

    CALL setKernelTime(411, kernelName//c_null_char, t3__-t2__, t2__-t1__, transfer_total, 0)

END SUBROUTINE

#ifdef OPS_LAZY
SUBROUTINE bounds_kernel_eqQ_yr_host( userSubroutine, block, dim, range, &
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
    opsArg11 &
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

    TYPE(ops_arg), DIMENSION(11), TARGET :: opsArgArray
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: namelit

    namelit = "bounds_kernel_eqQ_yr"

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

    DO n_indx = 1, 3
        range_tmp(2*n_indx-1) = range(2*n_indx-1)-1
        range_tmp(2*n_indx)   = range(2*n_indx)
    END DO

    CALL create_kerneldesc_and_enque(namelit//c_null_char, c_loc(opsArgArray), &
                                    11, 411, dim, 0, c_loc(range_tmp), &
                                    block%blockCptr, c_funloc(bounds_kernel_eqQ_yr_host_execute))

END SUBROUTINE
#endif

END MODULE BOUNDS_KERNEL_EQQ_YR_MODULE
