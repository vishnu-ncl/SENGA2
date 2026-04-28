! Auto-generated at 2026-04-28 18:44:28.589681 by ops-translator

MODULE D2FDXY_KERNEL_EQI_MODULE

    USE OPS_FORTRAN_DECLARATIONS
    USE OPS_FORTRAN_RT_SUPPORT

    USE OPS_CONSTANTS
    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

!$OMP DECLARE TARGET(xdim1_d2fdxy_kernel_eqI)
    INTEGER(KIND=4) :: xdim1_d2fdxy_kernel_eqI
    INTEGER(KIND=4) :: xdim1_d2fdxy_kernel_eqI_h = -1
!$OMP DECLARE TARGET(ydim1_d2fdxy_kernel_eqI)
    INTEGER(KIND=4) :: ydim1_d2fdxy_kernel_eqI
    INTEGER(KIND=4) :: ydim1_d2fdxy_kernel_eqI_h = -1
#define OPS_ACC1(x,y,z) ((x) + (xdim1_d2fdxy_kernel_eqI*(y)) + (xdim1_d2fdxy_kernel_eqI*ydim1_d2fdxy_kernel_eqI*(z)) + 1)

!$OMP DECLARE TARGET(xdim2_d2fdxy_kernel_eqI)
    INTEGER(KIND=4) :: xdim2_d2fdxy_kernel_eqI
    INTEGER(KIND=4) :: xdim2_d2fdxy_kernel_eqI_h = -1
!$OMP DECLARE TARGET(ydim2_d2fdxy_kernel_eqI)
    INTEGER(KIND=4) :: ydim2_d2fdxy_kernel_eqI
    INTEGER(KIND=4) :: ydim2_d2fdxy_kernel_eqI_h = -1
#define OPS_ACC2(x,y,z) ((x) + (xdim2_d2fdxy_kernel_eqI*(y)) + (xdim2_d2fdxy_kernel_eqI*ydim2_d2fdxy_kernel_eqI*(z)) + 1)

    CONTAINS

!   =============
!   User function
!   =============

SUBROUTINE d2fdxy_kernel_eqI(functn, fderiv)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv

    real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd

    fdiffa = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(1,1,0)) + functn(OPS_ACC1(1,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(1,2,0)) + functn(OPS_ACC1(1,-2,0)))
    fdiffb = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-1,1,0)) + functn(OPS_ACC1(-1,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(-1,2,0)) + functn(OPS_ACC1(-1,-2,0)))
    fdiffc = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-2,1,0)) + functn(OPS_ACC1(-2,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(-2,2,0)) + functn(OPS_ACC1(-2,-2,0)))
    fdiffd = acofy1*(functn(OPS_ACC1(0,1,0)) - functn(OPS_ACC1(0,-1,0)) &
                   - functn(OPS_ACC1(-3,1,0)) + functn(OPS_ACC1(-3,-1,0))) &
            + bcofy1*(functn(OPS_ACC1(0,2,0)) - functn(OPS_ACC1(0,-2,0)) &
                    - functn(OPS_ACC1(-3,2,0)) + functn(OPS_ACC1(-3,-2,0)))

    fderiv(OPS_ACC2(0,0,0)) = acf2xy*fdiffa + bcf2xy*fdiffb  &
          + ccf2xy*fdiffc + dcf2xy*fdiffd

END SUBROUTINE

#undef OPS_ACC1
#undef OPS_ACC2

SUBROUTINE d2fdxy_kernel_eqI_wrap( &
    opsDat1Local, &
    opsDat2Local, &
    dat1_base, &
    dat2_base, &
    dat1_dim, &
    dat2_dim, &
    start_indx, &
    end_indx )

    INTEGER(KIND=4), VALUE :: dat1_dim
    REAL(KIND=8), DIMENSION(dat1_dim), INTENT(IN) :: opsDat1Local
    INTEGER(KIND=4), VALUE :: dat1_base

    INTEGER(KIND=4), VALUE :: dat2_dim
    REAL(KIND=8), DIMENSION(dat2_dim), INTENT(OUT) :: opsDat2Local
    INTEGER(KIND=4), VALUE :: dat2_base

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

                CALL d2fdxy_kernel_eqI( &
                opsDat1Local(dat1_base + ((n_x-1)*1) + ((n_y-1)*xdim1_d2fdxy_kernel_eqI*1) + ((n_z-1)*ydim1_d2fdxy_kernel_eqI*xdim1_d2fdxy_kernel_eqI*1)), &
                opsDat2Local(dat2_base + ((n_x-1)*1) + ((n_y-1)*xdim2_d2fdxy_kernel_eqI*1) + ((n_z-1)*ydim2_d2fdxy_kernel_eqI*xdim2_d2fdxy_kernel_eqI*1)) &
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
SUBROUTINE d2fdxy_kernel_eqI_host( userSubroutine, block, dim, range, &
    opsArg1, &
    opsArg2 &
    )

    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: userSubroutine
    TYPE(ops_block), INTENT(IN) :: block
    INTEGER(KIND=4), INTENT(IN) :: dim
    INTEGER(KIND=4), DIMENSION(2*dim), INTENT(IN) :: range

    TYPE(ops_arg), INTENT(IN) :: opsArg1
    TYPE(ops_arg), INTENT(IN) :: opsArg2

    TYPE(ops_arg), DIMENSION(2) :: opsArgArray

#else
SUBROUTINE d2fdxy_kernel_eqI_host_execute( descPtr )

    TYPE(ops_kernel_descriptor), INTENT(IN) :: descPtr
    TYPE(ops_block) :: block
    INTEGER(KIND=C_INT) :: dim
    INTEGER(KIND=C_INT), POINTER, DIMENSION(:) :: range
    CHARACTER(KIND=C_CHAR), POINTER, DIMENSION(:) :: userSubroutine
    TYPE(ops_arg), POINTER, DIMENSION(:) :: opsArgArray

    TYPE(ops_arg) :: opsArg1
    TYPE(ops_arg) :: opsArg2

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

    REAL(KIND=8) :: t1__, t2__, t3__
    REAL(KIND=4) :: transfer_total, transfer

    INTEGER(KIND=4) :: x_size, y_size, z_size

    INTEGER(KIND=4), DIMENSION(3) :: start_indx, end_indx
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: kernelName

    kernelName = "d2fdxy_kernel_eqI"

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
#else
    opsArgArray(1) = opsArg1
    opsArgArray(2) = opsArg2
#endif

    CALL setKernelTime(26, kernelName//c_null_char, 0.0_8, 0.0_8, 0.0_4, 1)
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
    dat2_base = getDatBaseFromOpsArg3D(opsArg2, start_indx, 1)
    CALL c_f_pointer(opsArg2%data_d, opsDat2Local, [opsDat2Cardinality])

    IF (&
         (xdim1 .NE. xdim1_d2fdxy_kernel_eqI_h) .OR. &
         (ydim1 .NE. ydim1_d2fdxy_kernel_eqI_h) .OR. &
         (xdim2 .NE. xdim2_d2fdxy_kernel_eqI_h) .OR. &
         (ydim2 .NE. ydim2_d2fdxy_kernel_eqI_h) &
       ) THEN
            xdim1_d2fdxy_kernel_eqI = xdim1
            xdim1_d2fdxy_kernel_eqI_h = xdim1
!$OMP TARGET UPDATE TO(xdim1_d2fdxy_kernel_eqI)
            ydim1_d2fdxy_kernel_eqI = ydim1
            ydim1_d2fdxy_kernel_eqI_h = ydim1
!$OMP TARGET UPDATE TO(ydim1_d2fdxy_kernel_eqI)
            xdim2_d2fdxy_kernel_eqI = xdim2
            xdim2_d2fdxy_kernel_eqI_h = xdim2
!$OMP TARGET UPDATE TO(xdim2_d2fdxy_kernel_eqI)
            ydim2_d2fdxy_kernel_eqI = ydim2
            ydim2_d2fdxy_kernel_eqI_h = ydim2
!$OMP TARGET UPDATE TO(ydim2_d2fdxy_kernel_eqI)
    END IF

!   ==============
!   Halo exchanges
!   ==============
#ifndef OPS_LAZY
    CALL ops_H_D_exchanges_device(opsArgArray, 2)
    CALL ops_halo_exchanges(opsArgArray, 2, range, 3)
#endif

    CALL ops_timers_core(t2__)

!   ==============================
!   Call kernel wrapper subroutine
!   ==============================
    CALL d2fdxy_kernel_eqI_wrap( &
                        opsDat1Local, &
                        opsDat2Local, &
                        dat1_base, &
                        dat2_base, &
                        opsDat1Cardinality, &
                        opsDat2Cardinality, &
                        start_indx, &
                        end_indx )

    CALL ops_timers_core(t3__)

#ifndef OPS_LAZY
    CALL ops_set_dirtybit_device(opsArgArray, 2)
    CALL ops_set_halo_dirtybit3(opsArg2, range, 3)
#endif

!   ========================
!   Timing and data movement
!   ========================
    transfer_total = 0.0_4
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg1, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg2, transfer)
    transfer_total = transfer_total + transfer

    CALL setKernelTime(26, kernelName//c_null_char, t3__-t2__, t2__-t1__, transfer_total, 0)

END SUBROUTINE

#ifdef OPS_LAZY
SUBROUTINE d2fdxy_kernel_eqI_host( userSubroutine, block, dim, range, &
    opsArg1, &
    opsArg2 &
    )

    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN), TARGET :: userSubroutine
    TYPE(ops_block), INTENT(IN) :: block
    INTEGER(KIND=4), INTENT(IN) :: dim
    INTEGER(KIND=4), DIMENSION(2*dim), INTENT(INOUT), TARGET :: range
    INTEGER(KIND=4), DIMENSION(2*dim), TARGET :: range_tmp

    TYPE(ops_arg), INTENT(IN) :: opsArg1
    TYPE(ops_arg), INTENT(IN) :: opsArg2

    TYPE(ops_arg), DIMENSION(2), TARGET :: opsArgArray
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: namelit

    namelit = "d2fdxy_kernel_eqI"

    opsArgArray(1) = opsArg1
    opsArgArray(2) = opsArg2

    DO n_indx = 1, 3
        range_tmp(2*n_indx-1) = range(2*n_indx-1)-1
        range_tmp(2*n_indx)   = range(2*n_indx)
    END DO

    CALL create_kerneldesc_and_enque(namelit//c_null_char, c_loc(opsArgArray), &
                                    2, 26, dim, 1, c_loc(range_tmp), &
                                    block%blockCptr, c_funloc(d2fdxy_kernel_eqI_host_execute))

END SUBROUTINE
#endif

END MODULE D2FDXY_KERNEL_EQI_MODULE
