! Auto-generated at 2024-04-21 01:07:41.114854 by ops-translator

MODULE D2FDZ2_KERNEL_MAIN_MODULE

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

    CONTAINS

!   =============
!   User function
!   =============

!DEC$ ATTRIBUTES FORCEINLINE :: d2fdz2_kernel_main
SUBROUTINE d2fdz2_kernel_main(functn, fderiv, nzglbl, nendzl, nendzr, nbound, idx)

    real(kind=8), dimension(1), intent(in) :: functn
    real(kind=8), dimension(1) :: fderiv
    real(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep
    real(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm,fdifem

    integer(kind=4), intent(in) :: nzglbl, nendzl, nendzr, nbound
    integer(kind=4), dimension(3), intent(in) :: idx
    integer(kind=4) :: kc, kstart, kfinis

    kstart = 1
    kfinis = nzglbl
    IF(nendzl == nbound) kstart = 6
    IF(nendzr == nbound) kfinis = nzglbl-5

    kc = idx(3)

    IF (kc >= kstart .and. kc <= kfinis) THEN
!       interior scheme
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
        fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))
        fdifep = functn(OPS_ACC1(0,0,5)) - functn(OPS_ACC1(0,0,0))
        fdifem = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-5))

        fderiv(OPS_ACC2(0,0,0)) = acofsz*(fdifap-fdifam) + bcofsz*(fdifbp-fdifbm)  &
                                + ccofsz*(fdifcp-fdifcm) + dcofsz*(fdifdp-fdifdm)  &
                                + ecofsz*(fdifep-fdifem)

    ELSE IF (kc == 1) THEN
!       lh point: 4th order one-sided
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,0,5)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs1z*fdifap + bcfs1z*fdifbp  &
                                + ccfs1z*fdifcp + dcfs1z*fdifdp  &
                                + ecfs1z*fdifep

    ELSE IF (kc == 2) THEN
!       lh point plus 1: 4th order mixed
        fdifap = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs2z*fdifap + bcfs2z*fdifbp  &
                                + ccfs2z*fdifcp + dcfs2z*fdifdp  &
                                + ecfs2z*fdifep

    ELSE IF (kc == 3) THEN
!       lh point plus 2: 4th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acfs3z*(fdifap-fdifam) + bcfs3z*(fdifbp-fdifbm)

    ELSE IF (kc == 4) THEN
!       lh point plus 3: 6th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acfs4z*(fdifap-fdifam)  &
                                + bcfs4z*(fdifbp-fdifbm) + ccfs4z*(fdifcp-fdifcm)

    ELSE IF (kc == 5) THEN
!       lh point plus 4: 8th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
        fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acfs5z*(fdifap-fdifam)  &
                                + bcfs5z*(fdifbp-fdifbm) + ccfs5z*(fdifcp-fdifcm)  &
                                + dcfs5z*(fdifdp-fdifdm)

    ELSE IF (kc == nzglbl-4) THEN
!       rh point minus 4: 8th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))
        fdifdp = functn(OPS_ACC1(0,0,4)) - functn(OPS_ACC1(0,0,0))
        fdifdm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-4))

        fderiv(OPS_ACC2(0,0,0)) = acfs5z*(fdifap-fdifam)  &
                                + bcfs5z*(fdifbp-fdifbm) + ccfs5z*(fdifcp-fdifcm)  &
                                + dcfs5z*(fdifdp-fdifdm)

    ELSE IF (kc == nzglbl-3) THEN
!       rh point minus 3: 6th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))
        fdifcp = functn(OPS_ACC1(0,0,3)) - functn(OPS_ACC1(0,0,0))
        fdifcm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-3))

        fderiv(OPS_ACC2(0,0,0)) = acfs4z*(fdifap-fdifam)  &
                                + bcfs4z*(fdifbp-fdifbm) + ccfs4z*(fdifcp-fdifcm)

    ELSE IF (kc == nzglbl-2) THEN
!       rh point minus 2: 4th order centered
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifam = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-1))
        fdifbp = functn(OPS_ACC1(0,0,2)) - functn(OPS_ACC1(0,0,0))
        fdifbm = functn(OPS_ACC1(0,0,0)) - functn(OPS_ACC1(0,0,-2))

        fderiv(OPS_ACC2(0,0,0)) = acfs3z*(fdifap-fdifam) + bcfs3z*(fdifbp-fdifbm)

    ELSE IF (kc == nzglbl-1) THEN
!       rh point minus 1: 4th order mixed
        fdifap = functn(OPS_ACC1(0,0,1)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,0,-2)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,0,-3)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,0,-4)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs2z*fdifap + bcfs2z*fdifbp  &
                                + ccfs2z*fdifcp + dcfs2z*fdifdp  &
                                + ecfs2z*fdifep

    ELSE IF (kc == nzglbl) THEN
!       rh point: 4th order one-sided
        fdifap = functn(OPS_ACC1(0,0,-1)) - functn(OPS_ACC1(0,0,0))
        fdifbp = functn(OPS_ACC1(0,0,-2)) - functn(OPS_ACC1(0,0,0))
        fdifcp = functn(OPS_ACC1(0,0,-3)) - functn(OPS_ACC1(0,0,0))
        fdifdp = functn(OPS_ACC1(0,0,-4)) - functn(OPS_ACC1(0,0,0))
        fdifep = functn(OPS_ACC1(0,0,-5)) - functn(OPS_ACC1(0,0,0))

        fderiv(OPS_ACC2(0,0,0)) = acfs1z*fdifap + bcfs1z*fdifbp  &
                                + ccfs1z*fdifcp + dcfs1z*fdifdp  &
                                + ecfs1z*fdifep

    END IF

!   scaling
    fderiv(OPS_ACC2(0,0,0)) = fderiv(OPS_ACC2(0,0,0))*ovdlz2

END SUBROUTINE

#undef OPS_ACC1
#undef OPS_ACC2

SUBROUTINE d2fdz2_kernel_main_wrap( &
    opsDat1Local, &
    opsDat2Local, &
    opsDat3Local, &
    opsDat4Local, &
    opsDat5Local, &
    opsDat6Local, &
    idx, &
    dat1_base, &
    dat2_base, &
    dat3_base, &
    dat4_base, &
    dat5_base, &
    dat6_base, &
    start_indx, &
    end_indx )

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat1Local
    INTEGER(KIND=4), INTENT(IN) :: dat1_base

    REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: opsDat2Local
    INTEGER(KIND=4), INTENT(IN) :: dat2_base

    INTEGER(KIND=4), DIMENSION(*), INTENT(IN) :: opsDat3Local
    INTEGER(KIND=4), INTENT(IN) :: dat3_base

    INTEGER(KIND=4), DIMENSION(*), INTENT(IN) :: opsDat4Local
    INTEGER(KIND=4), INTENT(IN) :: dat4_base

    INTEGER(KIND=4), DIMENSION(*), INTENT(IN) :: opsDat5Local
    INTEGER(KIND=4), INTENT(IN) :: dat5_base

    INTEGER(KIND=4), DIMENSION(*), INTENT(IN) :: opsDat6Local
    INTEGER(KIND=4), INTENT(IN) :: dat6_base

    INTEGER(KIND=4), DIMENSION(3), INTENT(IN) :: idx
    INTEGER(KIND=4), DIMENSION(3)             :: idx_local

    INTEGER(KIND=4), DIMENSION(3), INTENT(IN) :: start_indx, end_indx

    INTEGER(KIND=4) :: n_x, n_y, n_z

    DO n_z = 1, end_indx(3)-start_indx(3)+1
        idx_local(3) = idx(3) + n_z - 1
        DO n_y = 1, end_indx(2)-start_indx(2)+1
            idx_local(2) = idx(2) + n_y - 1
            !$OMP SIMD
            DO n_x = 1, end_indx(1)-start_indx(1)+1
                idx_local(1) = idx(1) + n_x - 1

                CALL d2fdz2_kernel_main( &

                opsDat1Local(dat1_base + ((n_x-1)*1) + ((n_y-1)*xdim1*1) + ((n_z-1)*ydim1*xdim1*1)), &

                opsDat2Local(dat2_base + ((n_x-1)*1) + ((n_y-1)*xdim2*1) + ((n_z-1)*ydim2*xdim2*1)), &
                opsDat3Local(dat3_base), &
                opsDat4Local(dat4_base), &
                opsDat5Local(dat5_base), &
                opsDat6Local(dat6_base), &
                idx_local &
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
SUBROUTINE d2fdz2_kernel_main_host( userSubroutine, block, dim, range, &
    opsArg1, &
    opsArg2, &
    opsArg3, &
    opsArg4, &
    opsArg5, &
    opsArg6, &
    opsArg7 &
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

    TYPE(ops_arg), DIMENSION(7) :: opsArgArray

#else
SUBROUTINE d2fdz2_kernel_main_host_execute( descPtr )

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

#endif

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat1Local
    INTEGER(KIND=4) :: opsDat1Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat1_size
    INTEGER(KIND=4) :: dat1_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat2Local
    INTEGER(KIND=4) :: opsDat2Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat2_size
    INTEGER(KIND=4) :: dat2_base

    INTEGER(KIND=4), POINTER, DIMENSION(:) :: opsDat3Local
    INTEGER(KIND=4) :: dat3_base

    INTEGER(KIND=4), POINTER, DIMENSION(:) :: opsDat4Local
    INTEGER(KIND=4) :: dat4_base

    INTEGER(KIND=4), POINTER, DIMENSION(:) :: opsDat5Local
    INTEGER(KIND=4) :: dat5_base

    INTEGER(KIND=4), POINTER, DIMENSION(:) :: opsDat6Local
    INTEGER(KIND=4) :: dat6_base

    REAL(KIND=8) :: t1__, t2__, t3__
    REAL(KIND=4) :: transfer_total, transfer

    INTEGER(KIND=4), DIMENSION(3) :: idx

    INTEGER(KIND=4), DIMENSION(3) :: start_indx, end_indx
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: kernelName

    kernelName = "d2fdz2_kernel_main"

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
#else
    opsArgArray(1) = opsArg1
    opsArgArray(2) = opsArg2
    opsArgArray(3) = opsArg3
    opsArgArray(4) = opsArg4
    opsArgArray(5) = opsArg5
    opsArgArray(6) = opsArg6
    opsArgArray(7) = opsArg7
#endif

    CALL setKernelTime(15, kernelName//c_null_char, 0.0_8, 0.0_8, 0.0_4, 1)
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

#ifdef OPS_MPI
    CALL getIdx(block, start_indx, idx)
#else
    idx(1) = start_indx(1)
    idx(2) = start_indx(2)
    idx(3) = start_indx(3)
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

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg3), opsDat3Local, (/opsArg3%dim/))
    dat3_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg4), opsDat4Local, (/opsArg4%dim/))
    dat4_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg5), opsDat5Local, (/opsArg5%dim/))
    dat5_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg6), opsDat6Local, (/opsArg6%dim/))
    dat6_base = 1

!   ==============
!   Halo exchanges
!   ==============
#ifndef OPS_LAZY
    CALL ops_H_D_exchanges_host(opsArgArray, 7)
    CALL ops_halo_exchanges(opsArgArray, 7, range)
    CALL ops_H_D_exchanges_host(opsArgArray, 7)
#endif

    CALL ops_timers_core(t2__)

!   ==============================
!   Call kernel wrapper subroutine
!   ==============================
    CALL d2fdz2_kernel_main_wrap( &
                        opsDat1Local, &
                        opsDat2Local, &
                        opsDat3Local, &
                        opsDat4Local, &
                        opsDat5Local, &
                        opsDat6Local, &
                        idx, &
                        dat1_base, &
                        dat2_base, &
                        dat3_base, &
                        dat4_base, &
                        dat5_base, &
                        dat6_base, &
                        start_indx, &
                        end_indx )

    CALL ops_timers_core(t3__)

#ifndef OPS_LAZY
    CALL ops_set_dirtybit_host(opsArgArray, 7)
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

    CALL setKernelTime(15, kernelName//c_null_char, t3__-t2__, t2__-t1__, transfer_total, 0)

END SUBROUTINE

#ifdef OPS_LAZY
SUBROUTINE d2fdz2_kernel_main_host( userSubroutine, block, dim, range, &
    opsArg1, &
    opsArg2, &
    opsArg3, &
    opsArg4, &
    opsArg5, &
    opsArg6, &
    opsArg7 &
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

    TYPE(ops_arg), DIMENSION(7), TARGET :: opsArgArray
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: namelit

    namelit = "d2fdz2_kernel_main"

    opsArgArray(1) = opsArg1
    opsArgArray(2) = opsArg2
    opsArgArray(3) = opsArg3
    opsArgArray(4) = opsArg4
    opsArgArray(5) = opsArg5
    opsArgArray(6) = opsArg6
    opsArgArray(7) = opsArg7

    DO n_indx = 1, 3
        range_tmp(2*n_indx-1) = range(2*n_indx-1)-1
        range_tmp(2*n_indx)   = range(2*n_indx)
    END DO

    CALL create_kerneldesc_and_enque(namelit//c_null_char, c_loc(opsArgArray), &
                                    7, 15, dim, 0, c_loc(range_tmp), &
                                    block%blockCptr, c_funloc(d2fdz2_kernel_main_host_execute))

END SUBROUTINE
#endif

END MODULE D2FDZ2_KERNEL_MAIN_MODULE
