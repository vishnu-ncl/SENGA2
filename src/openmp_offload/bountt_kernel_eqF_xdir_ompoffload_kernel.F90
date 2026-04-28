! Auto-generated at 2026-04-28 18:44:32.790339 by ops-translator

MODULE BOUNTT_KERNEL_EQF_XDIR_MODULE

    USE OPS_FORTRAN_DECLARATIONS
    USE OPS_FORTRAN_RT_SUPPORT

    USE OPS_CONSTANTS
    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

!$OMP DECLARE TARGET(xdim1_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: xdim1_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: xdim1_bountt_kernel_eqF_xdir_h = -1
!$OMP DECLARE TARGET(ydim1_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: ydim1_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: ydim1_bountt_kernel_eqF_xdir_h = -1
#define OPS_ACC1(x,y,z) ((x) + (xdim1_bountt_kernel_eqF_xdir*(y)) + (xdim1_bountt_kernel_eqF_xdir*ydim1_bountt_kernel_eqF_xdir*(z)) + 1)

!$OMP DECLARE TARGET(xdim2_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: xdim2_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: xdim2_bountt_kernel_eqF_xdir_h = -1
!$OMP DECLARE TARGET(ydim2_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: ydim2_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: ydim2_bountt_kernel_eqF_xdir_h = -1
#define OPS_ACC2(x,y,z) ((x) + (xdim2_bountt_kernel_eqF_xdir*(y)) + (xdim2_bountt_kernel_eqF_xdir*ydim2_bountt_kernel_eqF_xdir*(z)) + 1)

!$OMP DECLARE TARGET(xdim3_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: xdim3_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: xdim3_bountt_kernel_eqF_xdir_h = -1
!$OMP DECLARE TARGET(ydim3_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: ydim3_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: ydim3_bountt_kernel_eqF_xdir_h = -1
#define OPS_ACC3(x,y,z) ((x) + (xdim3_bountt_kernel_eqF_xdir*(y)) + (xdim3_bountt_kernel_eqF_xdir*ydim3_bountt_kernel_eqF_xdir*(z)) + 1)

!$OMP DECLARE TARGET(xdim4_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: xdim4_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: xdim4_bountt_kernel_eqF_xdir_h = -1
!$OMP DECLARE TARGET(ydim4_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: ydim4_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: ydim4_bountt_kernel_eqF_xdir_h = -1
#define OPS_ACC4(x,y,z) ((x) + (xdim4_bountt_kernel_eqF_xdir*(y)) + (xdim4_bountt_kernel_eqF_xdir*ydim4_bountt_kernel_eqF_xdir*(z)) + 1)

!$OMP DECLARE TARGET(xdim5_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: xdim5_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: xdim5_bountt_kernel_eqF_xdir_h = -1
!$OMP DECLARE TARGET(ydim5_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: ydim5_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: ydim5_bountt_kernel_eqF_xdir_h = -1
#define OPS_ACC5(x,y,z) ((x) + (xdim5_bountt_kernel_eqF_xdir*(y)) + (xdim5_bountt_kernel_eqF_xdir*ydim5_bountt_kernel_eqF_xdir*(z)) + 1)

!$OMP DECLARE TARGET(xdim6_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: xdim6_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: xdim6_bountt_kernel_eqF_xdir_h = -1
!$OMP DECLARE TARGET(ydim6_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: ydim6_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: ydim6_bountt_kernel_eqF_xdir_h = -1
#define OPS_ACC6(x,y,z) ((x) + (xdim6_bountt_kernel_eqF_xdir*(y)) + (xdim6_bountt_kernel_eqF_xdir*ydim6_bountt_kernel_eqF_xdir*(z)) + 1)

!$OMP DECLARE TARGET(xdim7_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: xdim7_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: xdim7_bountt_kernel_eqF_xdir_h = -1
!$OMP DECLARE TARGET(ydim7_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: ydim7_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: ydim7_bountt_kernel_eqF_xdir_h = -1
#define OPS_ACC7(x,y,z) ((x) + (xdim7_bountt_kernel_eqF_xdir*(y)) + (xdim7_bountt_kernel_eqF_xdir*ydim7_bountt_kernel_eqF_xdir*(z)) + 1)

!$OMP DECLARE TARGET(xdim8_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: xdim8_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: xdim8_bountt_kernel_eqF_xdir_h = -1
!$OMP DECLARE TARGET(ydim8_bountt_kernel_eqF_xdir)
    INTEGER(KIND=4) :: ydim8_bountt_kernel_eqF_xdir
    INTEGER(KIND=4) :: ydim8_bountt_kernel_eqF_xdir_h = -1
#define OPS_ACC8(x,y,z) ((x) + (xdim8_bountt_kernel_eqF_xdir*(y)) + (xdim8_bountt_kernel_eqF_xdir*ydim8_bountt_kernel_eqF_xdir*(z)) + 1)

    CONTAINS

!   =============
!   User function
!   =============

SUBROUTINE bountt_kernel_eqF_xdir(erhs,yrhs,yrun,yerr,itndex,drhs,strtx,stryx,amasch,rgspec,ncpoly,ncpom1,ncenth,ispec,icoef1,icoef2)

    real(kind=8), dimension(1) :: erhs,yrhs,yrun,yerr
    integer(kind=4), dimension(1), intent(in) :: itndex
    real(kind=8), dimension(1), intent(in) :: drhs,strtx,stryx

    integer(kind=4), intent(in) :: ispec,icoef1,icoef2

    real(kind=8), dimension(ncofmx,ntinmx,nspcmx), intent(in) :: amasch
    real(kind=8), dimension(nspcmx), intent(in) :: rgspec
    integer(kind=4), dimension(ntinmx,nspcmx), intent(in) :: ncpoly,ncpom1,ncenth

    real(kind=8) :: fornow
    integer(kind=4) :: itint,icp

    itint = 1 +MOD(itndex(OPS_ACC5(0,0,0)),icoef1)/icoef2
    fornow = amasch(ncpoly(itint,ispec),itint,ispec)

    DO icp = ncpom1(itint,ispec),1,-1
        fornow = fornow*strtx(OPS_ACC7(0,0,0)) + amasch(icp,itint,ispec)
    END DO

    fornow = amasch(ncenth(itint,ispec),itint,ispec)  &
              + fornow*strtx(OPS_ACC7(0,0,0))

    yrhs(OPS_ACC2(0,0,0)) = drhs(OPS_ACC6(0,0,0))*stryx(OPS_ACC8(0,0,0))

    yrun(OPS_ACC3(0,0,0)) = yrhs(OPS_ACC2(0,0,0))

    yerr(OPS_ACC4(0,0,0)) = 0.0_8

    erhs(OPS_ACC1(0,0,0)) = erhs(OPS_ACC1(0,0,0))  &
                      + (fornow-rgspec(ispec)*strtx(OPS_ACC7(0,0,0)))*yrhs(OPS_ACC2(0,0,0))

END SUBROUTINE

#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6
#undef OPS_ACC7
#undef OPS_ACC8

SUBROUTINE bountt_kernel_eqF_xdir_wrap( &
    opsDat1Local, &
    opsDat2Local, &
    opsDat3Local, &
    opsDat4Local, &
    opsDat5Local, &
    opsDat6Local, &
    opsDat7Local, &
    opsDat8Local, &
    opsGblDat9Device, &
    opsGblDat10Device, &
    opsGblDat11Device, &
    opsGblDat12Device, &
    opsGblDat13Device, &
    opsGblDat14Device, &
    opsGblDat15Device, &
    opsGblDat16Device, &
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
    dat14_base, &
    dat15_base, &
    dat16_base, &
    dat1_dim, &
    dat2_dim, &
    dat3_dim, &
    dat4_dim, &
    dat5_dim, &
    dat6_dim, &
    dat7_dim, &
    dat8_dim, &
    dat9_dim, &
    dat10_dim, &
    dat11_dim, &
    dat12_dim, &
    dat13_dim, &
    dat14_dim, &
    dat15_dim, &
    dat16_dim, &
    start_indx, &
    end_indx )

    INTEGER(KIND=4), VALUE :: dat1_dim
    REAL(KIND=8), DIMENSION(dat1_dim), INTENT(INOUT) :: opsDat1Local
    INTEGER(KIND=4), VALUE :: dat1_base

    INTEGER(KIND=4), VALUE :: dat2_dim
    REAL(KIND=8), DIMENSION(dat2_dim), INTENT(OUT) :: opsDat2Local
    INTEGER(KIND=4), VALUE :: dat2_base

    INTEGER(KIND=4), VALUE :: dat3_dim
    REAL(KIND=8), DIMENSION(dat3_dim), INTENT(OUT) :: opsDat3Local
    INTEGER(KIND=4), VALUE :: dat3_base

    INTEGER(KIND=4), VALUE :: dat4_dim
    REAL(KIND=8), DIMENSION(dat4_dim), INTENT(OUT) :: opsDat4Local
    INTEGER(KIND=4), VALUE :: dat4_base

    INTEGER(KIND=4), VALUE :: dat5_dim
    INTEGER(KIND=4), DIMENSION(dat5_dim), INTENT(IN) :: opsDat5Local
    INTEGER(KIND=4), VALUE :: dat5_base

    INTEGER(KIND=4), VALUE :: dat6_dim
    REAL(KIND=8), DIMENSION(dat6_dim), INTENT(IN) :: opsDat6Local
    INTEGER(KIND=4), VALUE :: dat6_base

    INTEGER(KIND=4), VALUE :: dat7_dim
    REAL(KIND=8), DIMENSION(dat7_dim), INTENT(IN) :: opsDat7Local
    INTEGER(KIND=4), VALUE :: dat7_base

    INTEGER(KIND=4), VALUE :: dat8_dim
    REAL(KIND=8), DIMENSION(dat8_dim), INTENT(IN) :: opsDat8Local
    INTEGER(KIND=4), VALUE :: dat8_base

    INTEGER(KIND=4), VALUE :: dat9_dim
    REAL(KIND=8), DIMENSION(dat9_dim), INTENT(IN) :: opsGblDat9Device
    INTEGER(KIND=4), INTENT(IN) :: dat9_base

    INTEGER(KIND=4), VALUE :: dat10_dim
    REAL(KIND=8), DIMENSION(dat10_dim), INTENT(IN) :: opsGblDat10Device
    INTEGER(KIND=4), INTENT(IN) :: dat10_base

    INTEGER(KIND=4), VALUE :: dat11_dim
    INTEGER(KIND=4), DIMENSION(dat11_dim), INTENT(IN) :: opsGblDat11Device
    INTEGER(KIND=4), INTENT(IN) :: dat11_base

    INTEGER(KIND=4), VALUE :: dat12_dim
    INTEGER(KIND=4), DIMENSION(dat12_dim), INTENT(IN) :: opsGblDat12Device
    INTEGER(KIND=4), INTENT(IN) :: dat12_base

    INTEGER(KIND=4), VALUE :: dat13_dim
    INTEGER(KIND=4), DIMENSION(dat13_dim), INTENT(IN) :: opsGblDat13Device
    INTEGER(KIND=4), INTENT(IN) :: dat13_base

    INTEGER(KIND=4), VALUE :: dat14_dim
    INTEGER(KIND=4), VALUE ::  opsGblDat14Device
    INTEGER(KIND=4), INTENT(IN) :: dat14_base

    INTEGER(KIND=4), VALUE :: dat15_dim
    INTEGER(KIND=4), VALUE ::  opsGblDat15Device
    INTEGER(KIND=4), INTENT(IN) :: dat15_base

    INTEGER(KIND=4), VALUE :: dat16_dim
    INTEGER(KIND=4), VALUE ::  opsGblDat16Device
    INTEGER(KIND=4), INTENT(IN) :: dat16_base

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

                CALL bountt_kernel_eqF_xdir( &
                opsDat1Local(dat1_base + ((n_x-1)*1) + ((n_y-1)*xdim1_bountt_kernel_eqF_xdir*1) + ((n_z-1)*ydim1_bountt_kernel_eqF_xdir*xdim1_bountt_kernel_eqF_xdir*1)), &
                opsDat2Local(dat2_base + ((n_x-1)*1) + ((n_y-1)*xdim2_bountt_kernel_eqF_xdir*1) + ((n_z-1)*ydim2_bountt_kernel_eqF_xdir*xdim2_bountt_kernel_eqF_xdir*1)), &
                opsDat3Local(dat3_base + ((n_x-1)*1) + ((n_y-1)*xdim3_bountt_kernel_eqF_xdir*1) + ((n_z-1)*ydim3_bountt_kernel_eqF_xdir*xdim3_bountt_kernel_eqF_xdir*1)), &
                opsDat4Local(dat4_base + ((n_x-1)*1) + ((n_y-1)*xdim4_bountt_kernel_eqF_xdir*1) + ((n_z-1)*ydim4_bountt_kernel_eqF_xdir*xdim4_bountt_kernel_eqF_xdir*1)), &
                opsDat5Local(dat5_base + ((n_x-1)*1) + ((n_y-1)*xdim5_bountt_kernel_eqF_xdir*1) + ((n_z-1)*ydim5_bountt_kernel_eqF_xdir*xdim5_bountt_kernel_eqF_xdir*1)), &
                opsDat6Local(dat6_base + ((n_x-1)*1) + ((n_y-1)*xdim6_bountt_kernel_eqF_xdir*1) + ((n_z-1)*ydim6_bountt_kernel_eqF_xdir*xdim6_bountt_kernel_eqF_xdir*1)), &
                opsDat7Local(dat7_base + ((n_x-1)*0) + ((n_y-1)*xdim7_bountt_kernel_eqF_xdir*1) + ((n_z-1)*ydim7_bountt_kernel_eqF_xdir*xdim7_bountt_kernel_eqF_xdir*1)), &
                opsDat8Local(dat8_base + ((n_x-1)*0) + ((n_y-1)*xdim8_bountt_kernel_eqF_xdir*1) + ((n_z-1)*ydim8_bountt_kernel_eqF_xdir*xdim8_bountt_kernel_eqF_xdir*1)), &
                opsGblDat9Device(1), &
                opsGblDat10Device(1), &
                opsGblDat11Device(1), &
                opsGblDat12Device(1), &
                opsGblDat13Device(1), &
                opsGblDat14Device, &
                opsGblDat15Device, &
                opsGblDat16Device &
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
SUBROUTINE bountt_kernel_eqF_xdir_host( userSubroutine, block, dim, range, &
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
    opsArg13, &
    opsArg14, &
    opsArg15, &
    opsArg16 &
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
    TYPE(ops_arg), INTENT(IN) :: opsArg14
    TYPE(ops_arg), INTENT(IN) :: opsArg15
    TYPE(ops_arg), INTENT(IN) :: opsArg16

    TYPE(ops_arg), DIMENSION(16) :: opsArgArray

#else
SUBROUTINE bountt_kernel_eqF_xdir_host_execute( descPtr )

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
    TYPE(ops_arg) :: opsArg14
    TYPE(ops_arg) :: opsArg15
    TYPE(ops_arg) :: opsArg16

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

    REAL(KIND=8), DIMENSION(:), POINTER :: opsDat3Local
    INTEGER(KIND=4) :: opsDat3Cardinality
    INTEGER(KIND=4), DIMENSION(:), POINTER  :: dat3_size
    INTEGER(KIND=4) :: dat3_base
    INTEGER(KIND=4) :: xdim3, ydim3, zdim3

    REAL(KIND=8), DIMENSION(:), POINTER :: opsDat4Local
    INTEGER(KIND=4) :: opsDat4Cardinality
    INTEGER(KIND=4), DIMENSION(:), POINTER  :: dat4_size
    INTEGER(KIND=4) :: dat4_base
    INTEGER(KIND=4) :: xdim4, ydim4, zdim4

    INTEGER(KIND=4), DIMENSION(:), POINTER :: opsDat5Local
    INTEGER(KIND=4) :: opsDat5Cardinality
    INTEGER(KIND=4), DIMENSION(:), POINTER  :: dat5_size
    INTEGER(KIND=4) :: dat5_base
    INTEGER(KIND=4) :: xdim5, ydim5, zdim5

    REAL(KIND=8), DIMENSION(:), POINTER :: opsDat6Local
    INTEGER(KIND=4) :: opsDat6Cardinality
    INTEGER(KIND=4), DIMENSION(:), POINTER  :: dat6_size
    INTEGER(KIND=4) :: dat6_base
    INTEGER(KIND=4) :: xdim6, ydim6, zdim6

    REAL(KIND=8), DIMENSION(:), POINTER :: opsDat7Local
    INTEGER(KIND=4) :: opsDat7Cardinality
    INTEGER(KIND=4), DIMENSION(:), POINTER  :: dat7_size
    INTEGER(KIND=4) :: dat7_base
    INTEGER(KIND=4) :: xdim7, ydim7, zdim7

    REAL(KIND=8), DIMENSION(:), POINTER :: opsDat8Local
    INTEGER(KIND=4) :: opsDat8Cardinality
    INTEGER(KIND=4), DIMENSION(:), POINTER  :: dat8_size
    INTEGER(KIND=4) :: dat8_base
    INTEGER(KIND=4) :: xdim8, ydim8, zdim8

    REAL(KIND=8), DIMENSION(:), POINTER :: opsDat9Local
    INTEGER(KIND=4) :: opsDat9Cardinality
    INTEGER(KIND=4) :: dat9_base

    REAL(KIND=8), DIMENSION(:), POINTER :: opsDat10Local
    INTEGER(KIND=4) :: opsDat10Cardinality
    INTEGER(KIND=4) :: dat10_base

    INTEGER(KIND=4), DIMENSION(:), POINTER :: opsDat11Local
    INTEGER(KIND=4) :: opsDat11Cardinality
    INTEGER(KIND=4) :: dat11_base

    INTEGER(KIND=4), DIMENSION(:), POINTER :: opsDat12Local
    INTEGER(KIND=4) :: opsDat12Cardinality
    INTEGER(KIND=4) :: dat12_base

    INTEGER(KIND=4), DIMENSION(:), POINTER :: opsDat13Local
    INTEGER(KIND=4) :: opsDat13Cardinality
    INTEGER(KIND=4) :: dat13_base

    INTEGER(KIND=4), DIMENSION(:), POINTER :: opsDat14Local
    INTEGER(KIND=4) :: opsDat14Cardinality
    INTEGER(KIND=4) :: dat14_base

    INTEGER(KIND=4), DIMENSION(:), POINTER :: opsDat15Local
    INTEGER(KIND=4) :: opsDat15Cardinality
    INTEGER(KIND=4) :: dat15_base

    INTEGER(KIND=4), DIMENSION(:), POINTER :: opsDat16Local
    INTEGER(KIND=4) :: opsDat16Cardinality
    INTEGER(KIND=4) :: dat16_base

    REAL(KIND=8) :: t1__, t2__, t3__
    REAL(KIND=4) :: transfer_total, transfer

    INTEGER(KIND=4) :: x_size, y_size, z_size

    INTEGER(KIND=4), DIMENSION(3) :: start_indx, end_indx
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: kernelName

    kernelName = "bountt_kernel_eqF_xdir"

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
    opsArg14 = opsArgArray(14)
    opsArg15 = opsArgArray(15)
    opsArg16 = opsArgArray(16)
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
    opsArgArray(14) = opsArg14
    opsArgArray(15) = opsArg15
    opsArgArray(16) = opsArg16
#endif

    CALL setKernelTime(495, kernelName//c_null_char, 0.0_8, 0.0_8, 0.0_4, 1)
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

    CALL ops_upload_gbls(opsArgArray, 16)

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

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg3), dat3_size, [dim])
    xdim3 = dat3_size(1)
    ydim3 = dat3_size(2)
    zdim3 = dat3_size(3)
    opsDat3Cardinality = opsArg3%dim * xdim3 * ydim3 * zdim3
    dat3_base = getDatBaseFromOpsArg3D(opsArg3, start_indx, 1)
    CALL c_f_pointer(opsArg3%data_d, opsDat3Local, [opsDat3Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg4), dat4_size, [dim])
    xdim4 = dat4_size(1)
    ydim4 = dat4_size(2)
    zdim4 = dat4_size(3)
    opsDat4Cardinality = opsArg4%dim * xdim4 * ydim4 * zdim4
    dat4_base = getDatBaseFromOpsArg3D(opsArg4, start_indx, 1)
    CALL c_f_pointer(opsArg4%data_d, opsDat4Local, [opsDat4Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg5), dat5_size, [dim])
    xdim5 = dat5_size(1)
    ydim5 = dat5_size(2)
    zdim5 = dat5_size(3)
    opsDat5Cardinality = opsArg5%dim * xdim5 * ydim5 * zdim5
    dat5_base = getDatBaseFromOpsArg3D(opsArg5, start_indx, 1)
    CALL c_f_pointer(opsArg5%data_d, opsDat5Local, [opsDat5Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg6), dat6_size, [dim])
    xdim6 = dat6_size(1)
    ydim6 = dat6_size(2)
    zdim6 = dat6_size(3)
    opsDat6Cardinality = opsArg6%dim * xdim6 * ydim6 * zdim6
    dat6_base = getDatBaseFromOpsArg3D(opsArg6, start_indx, 1)
    CALL c_f_pointer(opsArg6%data_d, opsDat6Local, [opsDat6Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg7), dat7_size, [dim])
    xdim7 = dat7_size(1)
    ydim7 = dat7_size(2)
    zdim7 = dat7_size(3)
    opsDat7Cardinality = opsArg7%dim * xdim7 * ydim7 * zdim7
    dat7_base = getDatBaseFromOpsArg3D(opsArg7, start_indx, 1)
    CALL c_f_pointer(opsArg7%data_d, opsDat7Local, [opsDat7Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg8), dat8_size, [dim])
    xdim8 = dat8_size(1)
    ydim8 = dat8_size(2)
    zdim8 = dat8_size(3)
    opsDat8Cardinality = opsArg8%dim * xdim8 * ydim8 * zdim8
    dat8_base = getDatBaseFromOpsArg3D(opsArg8, start_indx, 1)
    CALL c_f_pointer(opsArg8%data_d, opsDat8Local, [opsDat8Cardinality])

    opsDat9Cardinality = opsArg9%dim
    CALL c_f_pointer(opsArgArray(9)%data_d, opsDat9Local, [opsDat9Cardinality])
    dat9_base = 1

    opsDat10Cardinality = opsArg10%dim
    CALL c_f_pointer(opsArgArray(10)%data_d, opsDat10Local, [opsDat10Cardinality])
    dat10_base = 1

    opsDat11Cardinality = opsArg11%dim
    CALL c_f_pointer(opsArgArray(11)%data_d, opsDat11Local, [opsDat11Cardinality])
    dat11_base = 1

    opsDat12Cardinality = opsArg12%dim
    CALL c_f_pointer(opsArgArray(12)%data_d, opsDat12Local, [opsDat12Cardinality])
    dat12_base = 1

    opsDat13Cardinality = opsArg13%dim
    CALL c_f_pointer(opsArgArray(13)%data_d, opsDat13Local, [opsDat13Cardinality])
    dat13_base = 1

    opsDat14Cardinality = opsArg14%dim
    CALL c_f_pointer(opsArg14%data, opsDat14Local, [opsDat14Cardinality])
    dat14_base = 1

    opsDat15Cardinality = opsArg15%dim
    CALL c_f_pointer(opsArg15%data, opsDat15Local, [opsDat15Cardinality])
    dat15_base = 1

    opsDat16Cardinality = opsArg16%dim
    CALL c_f_pointer(opsArg16%data, opsDat16Local, [opsDat16Cardinality])
    dat16_base = 1

    IF (&
         (xdim1 .NE. xdim1_bountt_kernel_eqF_xdir_h) .OR. &
         (ydim1 .NE. ydim1_bountt_kernel_eqF_xdir_h) .OR. &
         (xdim2 .NE. xdim2_bountt_kernel_eqF_xdir_h) .OR. &
         (ydim2 .NE. ydim2_bountt_kernel_eqF_xdir_h) .OR. &
         (xdim3 .NE. xdim3_bountt_kernel_eqF_xdir_h) .OR. &
         (ydim3 .NE. ydim3_bountt_kernel_eqF_xdir_h) .OR. &
         (xdim4 .NE. xdim4_bountt_kernel_eqF_xdir_h) .OR. &
         (ydim4 .NE. ydim4_bountt_kernel_eqF_xdir_h) .OR. &
         (xdim5 .NE. xdim5_bountt_kernel_eqF_xdir_h) .OR. &
         (ydim5 .NE. ydim5_bountt_kernel_eqF_xdir_h) .OR. &
         (xdim6 .NE. xdim6_bountt_kernel_eqF_xdir_h) .OR. &
         (ydim6 .NE. ydim6_bountt_kernel_eqF_xdir_h) .OR. &
         (xdim7 .NE. xdim7_bountt_kernel_eqF_xdir_h) .OR. &
         (ydim7 .NE. ydim7_bountt_kernel_eqF_xdir_h) .OR. &
         (xdim8 .NE. xdim8_bountt_kernel_eqF_xdir_h) .OR. &
         (ydim8 .NE. ydim8_bountt_kernel_eqF_xdir_h) &
       ) THEN
            xdim1_bountt_kernel_eqF_xdir = xdim1
            xdim1_bountt_kernel_eqF_xdir_h = xdim1
!$OMP TARGET UPDATE TO(xdim1_bountt_kernel_eqF_xdir)
            ydim1_bountt_kernel_eqF_xdir = ydim1
            ydim1_bountt_kernel_eqF_xdir_h = ydim1
!$OMP TARGET UPDATE TO(ydim1_bountt_kernel_eqF_xdir)
            xdim2_bountt_kernel_eqF_xdir = xdim2
            xdim2_bountt_kernel_eqF_xdir_h = xdim2
!$OMP TARGET UPDATE TO(xdim2_bountt_kernel_eqF_xdir)
            ydim2_bountt_kernel_eqF_xdir = ydim2
            ydim2_bountt_kernel_eqF_xdir_h = ydim2
!$OMP TARGET UPDATE TO(ydim2_bountt_kernel_eqF_xdir)
            xdim3_bountt_kernel_eqF_xdir = xdim3
            xdim3_bountt_kernel_eqF_xdir_h = xdim3
!$OMP TARGET UPDATE TO(xdim3_bountt_kernel_eqF_xdir)
            ydim3_bountt_kernel_eqF_xdir = ydim3
            ydim3_bountt_kernel_eqF_xdir_h = ydim3
!$OMP TARGET UPDATE TO(ydim3_bountt_kernel_eqF_xdir)
            xdim4_bountt_kernel_eqF_xdir = xdim4
            xdim4_bountt_kernel_eqF_xdir_h = xdim4
!$OMP TARGET UPDATE TO(xdim4_bountt_kernel_eqF_xdir)
            ydim4_bountt_kernel_eqF_xdir = ydim4
            ydim4_bountt_kernel_eqF_xdir_h = ydim4
!$OMP TARGET UPDATE TO(ydim4_bountt_kernel_eqF_xdir)
            xdim5_bountt_kernel_eqF_xdir = xdim5
            xdim5_bountt_kernel_eqF_xdir_h = xdim5
!$OMP TARGET UPDATE TO(xdim5_bountt_kernel_eqF_xdir)
            ydim5_bountt_kernel_eqF_xdir = ydim5
            ydim5_bountt_kernel_eqF_xdir_h = ydim5
!$OMP TARGET UPDATE TO(ydim5_bountt_kernel_eqF_xdir)
            xdim6_bountt_kernel_eqF_xdir = xdim6
            xdim6_bountt_kernel_eqF_xdir_h = xdim6
!$OMP TARGET UPDATE TO(xdim6_bountt_kernel_eqF_xdir)
            ydim6_bountt_kernel_eqF_xdir = ydim6
            ydim6_bountt_kernel_eqF_xdir_h = ydim6
!$OMP TARGET UPDATE TO(ydim6_bountt_kernel_eqF_xdir)
            xdim7_bountt_kernel_eqF_xdir = xdim7
            xdim7_bountt_kernel_eqF_xdir_h = xdim7
!$OMP TARGET UPDATE TO(xdim7_bountt_kernel_eqF_xdir)
            ydim7_bountt_kernel_eqF_xdir = ydim7
            ydim7_bountt_kernel_eqF_xdir_h = ydim7
!$OMP TARGET UPDATE TO(ydim7_bountt_kernel_eqF_xdir)
            xdim8_bountt_kernel_eqF_xdir = xdim8
            xdim8_bountt_kernel_eqF_xdir_h = xdim8
!$OMP TARGET UPDATE TO(xdim8_bountt_kernel_eqF_xdir)
            ydim8_bountt_kernel_eqF_xdir = ydim8
            ydim8_bountt_kernel_eqF_xdir_h = ydim8
!$OMP TARGET UPDATE TO(ydim8_bountt_kernel_eqF_xdir)
    END IF

!   ==============
!   Halo exchanges
!   ==============
#ifndef OPS_LAZY
    CALL ops_H_D_exchanges_device(opsArgArray, 16)
    CALL ops_halo_exchanges(opsArgArray, 16, range, 3)
#endif

    CALL ops_timers_core(t2__)

!   ==============================
!   Call kernel wrapper subroutine
!   ==============================
    CALL bountt_kernel_eqF_xdir_wrap( &
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
                        opsDat14Local(1), &
                        opsDat15Local(1), &
                        opsDat16Local(1), &
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
                        dat14_base, &
                        dat15_base, &
                        dat16_base, &
                        opsDat1Cardinality, &
                        opsDat2Cardinality, &
                        opsDat3Cardinality, &
                        opsDat4Cardinality, &
                        opsDat5Cardinality, &
                        opsDat6Cardinality, &
                        opsDat7Cardinality, &
                        opsDat8Cardinality, &
                        opsDat9Cardinality, &
                        opsDat10Cardinality, &
                        opsDat11Cardinality, &
                        opsDat12Cardinality, &
                        opsDat13Cardinality, &
                        opsDat14Cardinality, &
                        opsDat15Cardinality, &
                        opsDat16Cardinality, &
                        start_indx, &
                        end_indx )

    CALL ops_timers_core(t3__)

#ifndef OPS_LAZY
    CALL ops_set_dirtybit_device(opsArgArray, 16)
    CALL ops_set_halo_dirtybit3(opsArg1, range, 3)
    CALL ops_set_halo_dirtybit3(opsArg2, range, 3)
    CALL ops_set_halo_dirtybit3(opsArg3, range, 3)
    CALL ops_set_halo_dirtybit3(opsArg4, range, 3)
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

    CALL setKernelTime(495, kernelName//c_null_char, t3__-t2__, t2__-t1__, transfer_total, 0)

END SUBROUTINE

#ifdef OPS_LAZY
SUBROUTINE bountt_kernel_eqF_xdir_host( userSubroutine, block, dim, range, &
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
    opsArg13, &
    opsArg14, &
    opsArg15, &
    opsArg16 &
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
    TYPE(ops_arg), INTENT(IN) :: opsArg14
    TYPE(ops_arg), INTENT(IN) :: opsArg15
    TYPE(ops_arg), INTENT(IN) :: opsArg16

    TYPE(ops_arg), DIMENSION(16), TARGET :: opsArgArray
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: namelit

    namelit = "bountt_kernel_eqF_xdir"

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
    opsArgArray(14) = opsArg14
    opsArgArray(15) = opsArg15
    opsArgArray(16) = opsArg16

    DO n_indx = 1, 3
        range_tmp(2*n_indx-1) = range(2*n_indx-1)-1
        range_tmp(2*n_indx)   = range(2*n_indx)
    END DO

    CALL create_kerneldesc_and_enque(namelit//c_null_char, c_loc(opsArgArray), &
                                    16, 495, dim, 1, c_loc(range_tmp), &
                                    block%blockCptr, c_funloc(bountt_kernel_eqF_xdir_host_execute))

END SUBROUTINE
#endif

END MODULE BOUNTT_KERNEL_EQF_XDIR_MODULE
