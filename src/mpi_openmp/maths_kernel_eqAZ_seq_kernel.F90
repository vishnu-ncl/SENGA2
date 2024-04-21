! Auto-generated at 2024-04-21 01:07:44.681922 by ops-translator

MODULE MATHS_KERNEL_EQAZ_MODULE

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

!DEC$ ATTRIBUTES FORCEINLINE :: maths_kernel_eqAZ
SUBROUTINE maths_kernel_eqAZ(out_arr1,out_arr2,out_arr3,in_arr1,in_arr2,racnst,rncnst,reovrr,acfsri,bcfsrm,ovcsrm,dcfsri,ecfsri)

    real(kind=8), dimension(1) :: out_arr1, out_arr2, out_arr3
    real(kind=8), dimension(1), intent(in) :: in_arr1, in_arr2
    real(kind=8), intent(in) :: racnst,rncnst,reovrr,acfsri,bcfsrm,ovcsrm,dcfsri,ecfsri

    real(kind=8) :: preduc,fornow,fbroad

!   EVALUATE K0
    out_arr3(OPS_ACC3(0,0,0)) = racnst + rncnst*LOG(in_arr1(OPS_ACC4(0,0,0)))  &
                              - reovrr/in_arr1(OPS_ACC4(0,0,0))
    out_arr3(OPS_ACC3(0,0,0)) = EXP(out_arr3(OPS_ACC3(0,0,0)))

!   EVALUATE REDUCED PRESURE
    preduc = in_arr2(OPS_ACC5(0,0,0))*out_arr3(OPS_ACC3(0,0,0)) /out_arr1(OPS_ACC1(0,0,0))

!   EVALUATE FCENT
    fornow = bcfsrm/in_arr1(OPS_ACC4(0,0,0))
    fbroad = acfsri*EXP(fornow)
    fornow = in_arr1(OPS_ACC4(0,0,0))*ovcsrm
    fbroad = fbroad + EXP(fornow)
    fornow = LOG10(preduc)
    fornow = 1.0_8/(1.0_8+LOG10(fornow))
    fbroad = dcfsri*EXP(fornow*LOG(fbroad))
    fornow = fbroad*EXP(ecfsri*LOG(in_arr1(OPS_ACC4(0,0,0))))

!   EVALUATE UPDATED FORWARD RATE CONSTANT
    fornow = fbroad*preduc/(1.0_8+preduc)
    out_arr1(OPS_ACC1(0,0,0)) = out_arr1(OPS_ACC1(0,0,0))*fornow

!   RSC/ZN 08-AUG-2012 BUG FIX PRESSURE-DEPENDENT RATES
    out_arr2(OPS_ACC2(0,0,0)) = LOG(out_arr1(OPS_ACC1(0,0,0)))

END SUBROUTINE

#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5

SUBROUTINE maths_kernel_eqAZ_wrap( &
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
    opsDat12Local, &
    opsDat13Local, &
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
    dat12_base, &
    dat13_base, &
    start_indx, &
    end_indx )

    REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: opsDat1Local
    INTEGER(KIND=4), INTENT(IN) :: dat1_base

    REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: opsDat2Local
    INTEGER(KIND=4), INTENT(IN) :: dat2_base

    REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: opsDat3Local
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

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat11Local
    INTEGER(KIND=4), INTENT(IN) :: dat11_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat12Local
    INTEGER(KIND=4), INTENT(IN) :: dat12_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat13Local
    INTEGER(KIND=4), INTENT(IN) :: dat13_base

    INTEGER(KIND=4), DIMENSION(3), INTENT(IN) :: start_indx, end_indx

    INTEGER(KIND=4) :: n_x, n_y, n_z

    !$OMP PARALLEL DO PRIVATE(n_x,n_y,n_z)
    DO n_z = 1, end_indx(3)-start_indx(3)+1
        DO n_y = 1, end_indx(2)-start_indx(2)+1
            !$OMP SIMD
            DO n_x = 1, end_indx(1)-start_indx(1)+1

                CALL maths_kernel_eqAZ( &

                opsDat1Local(dat1_base + ((n_x-1)*1) + ((n_y-1)*xdim1*1) + ((n_z-1)*ydim1*xdim1*1)), &

                opsDat2Local(dat2_base + ((n_x-1)*1) + ((n_y-1)*xdim2*1) + ((n_z-1)*ydim2*xdim2*1)), &

                opsDat3Local(dat3_base + ((n_x-1)*1) + ((n_y-1)*xdim3*1) + ((n_z-1)*ydim3*xdim3*1)), &

                opsDat4Local(dat4_base + ((n_x-1)*1) + ((n_y-1)*xdim4*1) + ((n_z-1)*ydim4*xdim4*1)), &

                opsDat5Local(dat5_base + ((n_x-1)*1) + ((n_y-1)*xdim5*1) + ((n_z-1)*ydim5*xdim5*1)), &
                opsDat6Local(dat6_base), &
                opsDat7Local(dat7_base), &
                opsDat8Local(dat8_base), &
                opsDat9Local(dat9_base), &
                opsDat10Local(dat10_base), &
                opsDat11Local(dat11_base), &
                opsDat12Local(dat12_base), &
                opsDat13Local(dat13_base) &
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
SUBROUTINE maths_kernel_eqAZ_host( userSubroutine, block, dim, range, &
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
    TYPE(ops_arg), INTENT(IN) :: opsArg12
    TYPE(ops_arg), INTENT(IN) :: opsArg13

    TYPE(ops_arg), DIMENSION(13) :: opsArgArray

#else
SUBROUTINE maths_kernel_eqAZ_host_execute( descPtr )

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
    TYPE(ops_arg) :: opsArg12
    TYPE(ops_arg) :: opsArg13

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
    INTEGER(KIND=4) :: dat6_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat7Local
    INTEGER(KIND=4) :: dat7_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat8Local
    INTEGER(KIND=4) :: dat8_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat9Local
    INTEGER(KIND=4) :: dat9_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat10Local
    INTEGER(KIND=4) :: dat10_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat11Local
    INTEGER(KIND=4) :: dat11_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat12Local
    INTEGER(KIND=4) :: dat12_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat13Local
    INTEGER(KIND=4) :: dat13_base

    REAL(KIND=8) :: t1__, t2__, t3__
    REAL(KIND=4) :: transfer_total, transfer

    INTEGER(KIND=4), DIMENSION(3) :: start_indx, end_indx
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: kernelName

    kernelName = "maths_kernel_eqAZ"

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
    opsArg12 = opsArgArray(12)
    opsArg13 = opsArgArray(13)
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
    opsArgArray(12) = opsArg12
    opsArgArray(13) = opsArg13
#endif

    CALL setKernelTime(521, kernelName//c_null_char, 0.0_8, 0.0_8, 0.0_4, 1)
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

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg6), opsDat6Local, (/opsArg6%dim/))
    dat6_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg7), opsDat7Local, (/opsArg7%dim/))
    dat7_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg8), opsDat8Local, (/opsArg8%dim/))
    dat8_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg9), opsDat9Local, (/opsArg9%dim/))
    dat9_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg10), opsDat10Local, (/opsArg10%dim/))
    dat10_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg11), opsDat11Local, (/opsArg11%dim/))
    dat11_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg12), opsDat12Local, (/opsArg12%dim/))
    dat12_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg13), opsDat13Local, (/opsArg13%dim/))
    dat13_base = 1

!   ==============
!   Halo exchanges
!   ==============
#ifndef OPS_LAZY
    CALL ops_H_D_exchanges_host(opsArgArray, 13)
    CALL ops_halo_exchanges(opsArgArray, 13, range)
    CALL ops_H_D_exchanges_host(opsArgArray, 13)
#endif

    CALL ops_timers_core(t2__)

!   ==============================
!   Call kernel wrapper subroutine
!   ==============================
    CALL maths_kernel_eqAZ_wrap( &
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
                        opsDat12Local, &
                        opsDat13Local, &
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
                        dat12_base, &
                        dat13_base, &
                        start_indx, &
                        end_indx )

    CALL ops_timers_core(t3__)

#ifndef OPS_LAZY
    CALL ops_set_dirtybit_host(opsArgArray, 13)
    CALL ops_set_halo_dirtybit3(opsArg1, range)
    CALL ops_set_halo_dirtybit3(opsArg2, range)
    CALL ops_set_halo_dirtybit3(opsArg3, range)
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

    CALL setKernelTime(521, kernelName//c_null_char, t3__-t2__, t2__-t1__, transfer_total, 0)

END SUBROUTINE

#ifdef OPS_LAZY
SUBROUTINE maths_kernel_eqAZ_host( userSubroutine, block, dim, range, &
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

    namelit = "maths_kernel_eqAZ"

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

    CALL create_kerneldesc_and_enque(namelit//c_null_char, c_loc(opsArgArray), &
                                    13, 521, dim, 0, c_loc(range_tmp), &
                                    block%blockCptr, c_funloc(maths_kernel_eqAZ_host_execute))

END SUBROUTINE
#endif

END MODULE MATHS_KERNEL_EQAZ_MODULE
