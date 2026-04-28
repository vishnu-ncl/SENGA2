! Auto-generated at 2026-04-28 18:43:12.234708 by ops-translator

MODULE BOUNDS_KERNEL_EQAB_XDIR_MODULE

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

    INTEGER(KIND=4) :: xdim11
    INTEGER(KIND=4) :: ydim11
    INTEGER(KIND=4) :: zdim11
#define OPS_ACC11(x,y,z) ((x) + (xdim11*(y)) + (xdim11*ydim11*(z)) + 1)

    INTEGER(KIND=4) :: xdim12
    INTEGER(KIND=4) :: ydim12
    INTEGER(KIND=4) :: zdim12
#define OPS_ACC12(x,y,z) ((x) + (xdim12*(y)) + (xdim12*ydim12*(z)) + 1)

    INTEGER(KIND=4) :: xdim13
    INTEGER(KIND=4) :: ydim13
    INTEGER(KIND=4) :: zdim13
#define OPS_ACC13(x,y,z) ((x) + (xdim13*(y)) + (xdim13*ydim13*(z)) + 1)

    INTEGER(KIND=4) :: xdim14
    INTEGER(KIND=4) :: ydim14
    INTEGER(KIND=4) :: zdim14
#define OPS_ACC14(x,y,z) ((x) + (xdim14*(y)) + (xdim14*ydim14*(z)) + 1)

    INTEGER(KIND=4) :: xdim15
    INTEGER(KIND=4) :: ydim15
    INTEGER(KIND=4) :: zdim15
#define OPS_ACC15(x,y,z) ((x) + (xdim15*(y)) + (xdim15*ydim15*(z)) + 1)

    INTEGER(KIND=4) :: xdim16
    INTEGER(KIND=4) :: ydim16
    INTEGER(KIND=4) :: zdim16
#define OPS_ACC16(x,y,z) ((x) + (xdim16*(y)) + (xdim16*ydim16*(z)) + 1)

    INTEGER(KIND=4) :: xdim17
    INTEGER(KIND=4) :: ydim17
    INTEGER(KIND=4) :: zdim17
#define OPS_ACC17(x,y,z) ((x) + (xdim17*(y)) + (xdim17*ydim17*(z)) + 1)

    INTEGER(KIND=4) :: xdim18
    INTEGER(KIND=4) :: ydim18
    INTEGER(KIND=4) :: zdim18
#define OPS_ACC18(x,y,z) ((x) + (xdim18*(y)) + (xdim18*ydim18*(z)) + 1)

    INTEGER(KIND=4) :: xdim19
    INTEGER(KIND=4) :: ydim19
    INTEGER(KIND=4) :: zdim19
#define OPS_ACC19(x,y,z) ((x) + (xdim19*(y)) + (xdim19*ydim19*(z)) + 1)

    INTEGER(KIND=4) :: xdim20
    INTEGER(KIND=4) :: ydim20
    INTEGER(KIND=4) :: zdim20
#define OPS_ACC20(x,y,z) ((x) + (xdim20*(y)) + (xdim20*ydim20*(z)) + 1)

    INTEGER(KIND=4) :: xdim21
    INTEGER(KIND=4) :: ydim21
    INTEGER(KIND=4) :: zdim21
#define OPS_ACC21(x,y,z) ((x) + (xdim21*(y)) + (xdim21*ydim21*(z)) + 1)

    CONTAINS

!   =============
!   User function
!   =============

#ifdef __INTEL_COMPILER
!DEC$ ATTRIBUTES FORCEINLINE :: bounds_kernel_eqAB_xdir
#endif
SUBROUTINE bounds_kernel_eqAB_xdir(bcl2x,bcl3x,bcl4x,bcl5x,bcl1x,strdx,acoux,ova2x,tt2x,tt3x,tt4x,tt5x,strrx,stlux,strux,stlvx,strvx,stlwx,strwx,stltx,strtx,xgdlen,nrieta2,nrieta3,nrieta4,nrieta5,m2max)

    real(kind=8), dimension(1) :: bcl2x,bcl3x,bcl4x,bcl5x
    real(kind=8), dimension(1), intent(in) :: bcl1x,strdx,acoux,ova2x,tt2x,tt3x,tt4x,tt5x
    real(kind=8), dimension(1), intent(in) :: strrx,stlux,strux,stlvx,strvx,stlwx,strwx,stltx,strtx
    real(kind=8), intent(in) :: xgdlen,nrieta2,nrieta3,nrieta4,nrieta5,m2max
    real(kind=8) :: fornow

!   OLD VALUE OF L's
    fornow = strdx(OPS_ACC6(0,0,0))*acoux(OPS_ACC7(0,0,0))*bcl1x(OPS_ACC5(0,0,0))

    bcl2x(OPS_ACC1(0,0,0)) = stlux(OPS_ACC14(0,0,0))  &
                  * (bcl2x(OPS_ACC1(0,0,0))-bcl5x(OPS_ACC4(0,0,0))*ova2x(OPS_ACC8(0,0,0)))

    bcl3x(OPS_ACC2(0,0,0)) = stlux(OPS_ACC14(0,0,0))*bcl3x(OPS_ACC2(0,0,0))

    bcl4x(OPS_ACC3(0,0,0)) = stlux(OPS_ACC14(0,0,0))*bcl4x(OPS_ACC3(0,0,0))

    bcl5x(OPS_ACC4(0,0,0)) = 0.5_8*(stlux(OPS_ACC14(0,0,0))+acoux(OPS_ACC7(0,0,0)))  &
                  * (bcl5x(OPS_ACC4(0,0,0))+fornow)

!   SUBTRACT FROM NEW VALUE OF L's
!   L1X UNCHANGED
    bcl2x(OPS_ACC1(0,0,0)) = nrieta2 * strdx(OPS_ACC6(0,0,0))  &
                  * strrx(OPS_ACC13(0,0,0))/(xgdlen*acoux(OPS_ACC7(0,0,0)))  &
                  * (stltx(OPS_ACC20(0,0,0))-strtx(OPS_ACC21(0,0,0)))  &
                  + tt2x(OPS_ACC9(0,0,0))*ova2x(OPS_ACC8(0,0,0))-bcl2x(OPS_ACC1(0,0,0))

    bcl3x(OPS_ACC2(0,0,0)) = nrieta3 * acoux(OPS_ACC7(0,0,0))  &
                  / xgdlen*(stlvx(OPS_ACC16(0,0,0))-strvx(OPS_ACC17(0,0,0))) +tt3x(OPS_ACC10(0,0,0))-bcl3x(OPS_ACC2(0,0,0))

    bcl4x(OPS_ACC3(0,0,0)) = nrieta4 * acoux(OPS_ACC7(0,0,0))  &
                  / xgdlen*(stlwx(OPS_ACC18(0,0,0))-strwx(OPS_ACC19(0,0,0))) +tt4x(OPS_ACC11(0,0,0))-bcl4x(OPS_ACC3(0,0,0))

    bcl5x(OPS_ACC4(0,0,0)) = nrieta5 *strdx(OPS_ACC6(0,0,0))*acoux(OPS_ACC7(0,0,0))  &
                  * acoux(OPS_ACC7(0,0,0))*(1.0_8-m2max)/(2.0_8*xgdlen)  &
                  * (stlux(OPS_ACC14(0,0,0))-strux(OPS_ACC15(0,0,0))) +0.5_8*tt5x(OPS_ACC12(0,0,0))  &
                  - bcl5x(OPS_ACC4(0,0,0))

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
#undef OPS_ACC11
#undef OPS_ACC12
#undef OPS_ACC13
#undef OPS_ACC14
#undef OPS_ACC15
#undef OPS_ACC16
#undef OPS_ACC17
#undef OPS_ACC18
#undef OPS_ACC19
#undef OPS_ACC20
#undef OPS_ACC21

SUBROUTINE bounds_kernel_eqAB_xdir_wrap( &
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
    opsDat14Local, &
    opsDat15Local, &
    opsDat16Local, &
    opsDat17Local, &
    opsDat18Local, &
    opsDat19Local, &
    opsDat20Local, &
    opsDat21Local, &
    opsDat22Local, &
    opsDat23Local, &
    opsDat24Local, &
    opsDat25Local, &
    opsDat26Local, &
    opsDat27Local, &
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
    dat17_base, &
    dat18_base, &
    dat19_base, &
    dat20_base, &
    dat21_base, &
    dat22_base, &
    dat23_base, &
    dat24_base, &
    dat25_base, &
    dat26_base, &
    dat27_base, &
    start_indx, &
    end_indx )

    REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: opsDat1Local
    INTEGER(KIND=4), INTENT(IN) :: dat1_base

    REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: opsDat2Local
    INTEGER(KIND=4), INTENT(IN) :: dat2_base

    REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: opsDat3Local
    INTEGER(KIND=4), INTENT(IN) :: dat3_base

    REAL(KIND=8), DIMENSION(*), INTENT(INOUT) :: opsDat4Local
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

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat14Local
    INTEGER(KIND=4), INTENT(IN) :: dat14_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat15Local
    INTEGER(KIND=4), INTENT(IN) :: dat15_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat16Local
    INTEGER(KIND=4), INTENT(IN) :: dat16_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat17Local
    INTEGER(KIND=4), INTENT(IN) :: dat17_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat18Local
    INTEGER(KIND=4), INTENT(IN) :: dat18_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat19Local
    INTEGER(KIND=4), INTENT(IN) :: dat19_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat20Local
    INTEGER(KIND=4), INTENT(IN) :: dat20_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat21Local
    INTEGER(KIND=4), INTENT(IN) :: dat21_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat22Local
    INTEGER(KIND=4), INTENT(IN) :: dat22_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat23Local
    INTEGER(KIND=4), INTENT(IN) :: dat23_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat24Local
    INTEGER(KIND=4), INTENT(IN) :: dat24_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat25Local
    INTEGER(KIND=4), INTENT(IN) :: dat25_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat26Local
    INTEGER(KIND=4), INTENT(IN) :: dat26_base

    REAL(KIND=8), DIMENSION(*), INTENT(IN) :: opsDat27Local
    INTEGER(KIND=4), INTENT(IN) :: dat27_base

    INTEGER(KIND=4), DIMENSION(3), INTENT(IN) :: start_indx, end_indx

    INTEGER(KIND=4) :: n_x, n_y, n_z

    !$OMP PARALLEL DO PRIVATE(n_x,n_y,n_z)
    DO n_z = 1, end_indx(3)-start_indx(3)+1
        DO n_y = 1, end_indx(2)-start_indx(2)+1
            !$OMP SIMD
            DO n_x = 1, end_indx(1)-start_indx(1)+1

#ifdef _CRAYFTN
                !DIR$ INLINE
#endif
                CALL bounds_kernel_eqAB_xdir( &

                opsDat1Local(dat1_base + ((n_x-1)*0) + ((n_y-1)*xdim1*1) + ((n_z-1)*ydim1*xdim1*1)), &

                opsDat2Local(dat2_base + ((n_x-1)*0) + ((n_y-1)*xdim2*1) + ((n_z-1)*ydim2*xdim2*1)), &

                opsDat3Local(dat3_base + ((n_x-1)*0) + ((n_y-1)*xdim3*1) + ((n_z-1)*ydim3*xdim3*1)), &

                opsDat4Local(dat4_base + ((n_x-1)*0) + ((n_y-1)*xdim4*1) + ((n_z-1)*ydim4*xdim4*1)), &

                opsDat5Local(dat5_base + ((n_x-1)*0) + ((n_y-1)*xdim5*1) + ((n_z-1)*ydim5*xdim5*1)), &

                opsDat6Local(dat6_base + ((n_x-1)*0) + ((n_y-1)*xdim6*1) + ((n_z-1)*ydim6*xdim6*1)), &

                opsDat7Local(dat7_base + ((n_x-1)*0) + ((n_y-1)*xdim7*1) + ((n_z-1)*ydim7*xdim7*1)), &

                opsDat8Local(dat8_base + ((n_x-1)*0) + ((n_y-1)*xdim8*1) + ((n_z-1)*ydim8*xdim8*1)), &

                opsDat9Local(dat9_base + ((n_x-1)*0) + ((n_y-1)*xdim9*1) + ((n_z-1)*ydim9*xdim9*1)), &

                opsDat10Local(dat10_base + ((n_x-1)*0) + ((n_y-1)*xdim10*1) + ((n_z-1)*ydim10*xdim10*1)), &

                opsDat11Local(dat11_base + ((n_x-1)*0) + ((n_y-1)*xdim11*1) + ((n_z-1)*ydim11*xdim11*1)), &

                opsDat12Local(dat12_base + ((n_x-1)*0) + ((n_y-1)*xdim12*1) + ((n_z-1)*ydim12*xdim12*1)), &

                opsDat13Local(dat13_base + ((n_x-1)*0) + ((n_y-1)*xdim13*1) + ((n_z-1)*ydim13*xdim13*1)), &

                opsDat14Local(dat14_base + ((n_x-1)*0) + ((n_y-1)*xdim14*1) + ((n_z-1)*ydim14*xdim14*1)), &

                opsDat15Local(dat15_base + ((n_x-1)*0) + ((n_y-1)*xdim15*1) + ((n_z-1)*ydim15*xdim15*1)), &

                opsDat16Local(dat16_base + ((n_x-1)*0) + ((n_y-1)*xdim16*1) + ((n_z-1)*ydim16*xdim16*1)), &

                opsDat17Local(dat17_base + ((n_x-1)*0) + ((n_y-1)*xdim17*1) + ((n_z-1)*ydim17*xdim17*1)), &

                opsDat18Local(dat18_base + ((n_x-1)*0) + ((n_y-1)*xdim18*1) + ((n_z-1)*ydim18*xdim18*1)), &

                opsDat19Local(dat19_base + ((n_x-1)*0) + ((n_y-1)*xdim19*1) + ((n_z-1)*ydim19*xdim19*1)), &

                opsDat20Local(dat20_base + ((n_x-1)*0) + ((n_y-1)*xdim20*1) + ((n_z-1)*ydim20*xdim20*1)), &

                opsDat21Local(dat21_base + ((n_x-1)*0) + ((n_y-1)*xdim21*1) + ((n_z-1)*ydim21*xdim21*1)), &
                opsDat22Local(dat22_base), &
                opsDat23Local(dat23_base), &
                opsDat24Local(dat24_base), &
                opsDat25Local(dat25_base), &
                opsDat26Local(dat26_base), &
                opsDat27Local(dat27_base) &
               )

            END DO
            !$OMP END SIMD
        END DO
    END DO

END SUBROUTINE bounds_kernel_eqAB_xdir_wrap

!   ===============
!   Host subroutine
!   ===============
#ifndef OPS_LAZY
SUBROUTINE bounds_kernel_eqAB_xdir_host( userSubroutine, block, dim, range, &
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
    opsArg16, &
    opsArg17, &
    opsArg18, &
    opsArg19, &
    opsArg20, &
    opsArg21, &
    opsArg22, &
    opsArg23, &
    opsArg24, &
    opsArg25, &
    opsArg26, &
    opsArg27 &
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
    TYPE(ops_arg), INTENT(IN) :: opsArg17
    TYPE(ops_arg), INTENT(IN) :: opsArg18
    TYPE(ops_arg), INTENT(IN) :: opsArg19
    TYPE(ops_arg), INTENT(IN) :: opsArg20
    TYPE(ops_arg), INTENT(IN) :: opsArg21
    TYPE(ops_arg), INTENT(IN) :: opsArg22
    TYPE(ops_arg), INTENT(IN) :: opsArg23
    TYPE(ops_arg), INTENT(IN) :: opsArg24
    TYPE(ops_arg), INTENT(IN) :: opsArg25
    TYPE(ops_arg), INTENT(IN) :: opsArg26
    TYPE(ops_arg), INTENT(IN) :: opsArg27

    TYPE(ops_arg), DIMENSION(27) :: opsArgArray

#else
SUBROUTINE bounds_kernel_eqAB_xdir_host_execute( descPtr )

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
    TYPE(ops_arg) :: opsArg17
    TYPE(ops_arg) :: opsArg18
    TYPE(ops_arg) :: opsArg19
    TYPE(ops_arg) :: opsArg20
    TYPE(ops_arg) :: opsArg21
    TYPE(ops_arg) :: opsArg22
    TYPE(ops_arg) :: opsArg23
    TYPE(ops_arg) :: opsArg24
    TYPE(ops_arg) :: opsArg25
    TYPE(ops_arg) :: opsArg26
    TYPE(ops_arg) :: opsArg27

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

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat11Local
    INTEGER(KIND=4) :: opsDat11Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat11_size
    INTEGER(KIND=4) :: dat11_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat12Local
    INTEGER(KIND=4) :: opsDat12Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat12_size
    INTEGER(KIND=4) :: dat12_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat13Local
    INTEGER(KIND=4) :: opsDat13Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat13_size
    INTEGER(KIND=4) :: dat13_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat14Local
    INTEGER(KIND=4) :: opsDat14Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat14_size
    INTEGER(KIND=4) :: dat14_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat15Local
    INTEGER(KIND=4) :: opsDat15Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat15_size
    INTEGER(KIND=4) :: dat15_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat16Local
    INTEGER(KIND=4) :: opsDat16Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat16_size
    INTEGER(KIND=4) :: dat16_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat17Local
    INTEGER(KIND=4) :: opsDat17Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat17_size
    INTEGER(KIND=4) :: dat17_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat18Local
    INTEGER(KIND=4) :: opsDat18Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat18_size
    INTEGER(KIND=4) :: dat18_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat19Local
    INTEGER(KIND=4) :: opsDat19Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat19_size
    INTEGER(KIND=4) :: dat19_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat20Local
    INTEGER(KIND=4) :: opsDat20Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat20_size
    INTEGER(KIND=4) :: dat20_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat21Local
    INTEGER(KIND=4) :: opsDat21Cardinality
    INTEGER(KIND=4), POINTER, DIMENSION(:)  :: dat21_size
    INTEGER(KIND=4) :: dat21_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat22Local
    INTEGER(KIND=4) :: dat22_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat23Local
    INTEGER(KIND=4) :: dat23_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat24Local
    INTEGER(KIND=4) :: dat24_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat25Local
    INTEGER(KIND=4) :: dat25_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat26Local
    INTEGER(KIND=4) :: dat26_base

    REAL(KIND=8), POINTER, DIMENSION(:) :: opsDat27Local
    INTEGER(KIND=4) :: dat27_base

    REAL(KIND=8) :: t1__, t2__, t3__
    REAL(KIND=4) :: transfer_total, transfer

    INTEGER(KIND=4), DIMENSION(3) :: start_indx, end_indx
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: kernelName

    kernelName = "bounds_kernel_eqAB_xdir"

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
    opsArg17 = opsArgArray(17)
    opsArg18 = opsArgArray(18)
    opsArg19 = opsArgArray(19)
    opsArg20 = opsArgArray(20)
    opsArg21 = opsArgArray(21)
    opsArg22 = opsArgArray(22)
    opsArg23 = opsArgArray(23)
    opsArg24 = opsArgArray(24)
    opsArg25 = opsArgArray(25)
    opsArg26 = opsArgArray(26)
    opsArg27 = opsArgArray(27)
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
    opsArgArray(17) = opsArg17
    opsArgArray(18) = opsArg18
    opsArgArray(19) = opsArg19
    opsArgArray(20) = opsArg20
    opsArgArray(21) = opsArg21
    opsArgArray(22) = opsArg22
    opsArgArray(23) = opsArg23
    opsArgArray(24) = opsArg24
    opsArgArray(25) = opsArg25
    opsArgArray(26) = opsArg26
    opsArgArray(27) = opsArg27
#endif

    CALL setKernelTime(356, kernelName//c_null_char, 0.0_8, 0.0_8, 0.0_4, 1)
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

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg1), dat1_size, [dim])
    xdim1 = dat1_size(1)
    ydim1 = dat1_size(2)
    zdim1 = dat1_size(3)
    opsDat1Cardinality = opsArg1%dim * xdim1 * ydim1 * zdim1
    dat1_base = getDatBaseFromOpsArg3D(opsArg1, start_indx, 1)
    CALL c_f_pointer(opsArg1%data, opsDat1Local, [opsDat1Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg2), dat2_size, [dim])
    xdim2 = dat2_size(1)
    ydim2 = dat2_size(2)
    zdim2 = dat2_size(3)
    opsDat2Cardinality = opsArg2%dim * xdim2 * ydim2 * zdim2
    dat2_base = getDatBaseFromOpsArg3D(opsArg2, start_indx, 1)
    CALL c_f_pointer(opsArg2%data, opsDat2Local, [opsDat2Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg3), dat3_size, [dim])
    xdim3 = dat3_size(1)
    ydim3 = dat3_size(2)
    zdim3 = dat3_size(3)
    opsDat3Cardinality = opsArg3%dim * xdim3 * ydim3 * zdim3
    dat3_base = getDatBaseFromOpsArg3D(opsArg3, start_indx, 1)
    CALL c_f_pointer(opsArg3%data, opsDat3Local, [opsDat3Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg4), dat4_size, [dim])
    xdim4 = dat4_size(1)
    ydim4 = dat4_size(2)
    zdim4 = dat4_size(3)
    opsDat4Cardinality = opsArg4%dim * xdim4 * ydim4 * zdim4
    dat4_base = getDatBaseFromOpsArg3D(opsArg4, start_indx, 1)
    CALL c_f_pointer(opsArg4%data, opsDat4Local, [opsDat4Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg5), dat5_size, [dim])
    xdim5 = dat5_size(1)
    ydim5 = dat5_size(2)
    zdim5 = dat5_size(3)
    opsDat5Cardinality = opsArg5%dim * xdim5 * ydim5 * zdim5
    dat5_base = getDatBaseFromOpsArg3D(opsArg5, start_indx, 1)
    CALL c_f_pointer(opsArg5%data, opsDat5Local, [opsDat5Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg6), dat6_size, [dim])
    xdim6 = dat6_size(1)
    ydim6 = dat6_size(2)
    zdim6 = dat6_size(3)
    opsDat6Cardinality = opsArg6%dim * xdim6 * ydim6 * zdim6
    dat6_base = getDatBaseFromOpsArg3D(opsArg6, start_indx, 1)
    CALL c_f_pointer(opsArg6%data, opsDat6Local, [opsDat6Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg7), dat7_size, [dim])
    xdim7 = dat7_size(1)
    ydim7 = dat7_size(2)
    zdim7 = dat7_size(3)
    opsDat7Cardinality = opsArg7%dim * xdim7 * ydim7 * zdim7
    dat7_base = getDatBaseFromOpsArg3D(opsArg7, start_indx, 1)
    CALL c_f_pointer(opsArg7%data, opsDat7Local, [opsDat7Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg8), dat8_size, [dim])
    xdim8 = dat8_size(1)
    ydim8 = dat8_size(2)
    zdim8 = dat8_size(3)
    opsDat8Cardinality = opsArg8%dim * xdim8 * ydim8 * zdim8
    dat8_base = getDatBaseFromOpsArg3D(opsArg8, start_indx, 1)
    CALL c_f_pointer(opsArg8%data, opsDat8Local, [opsDat8Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg9), dat9_size, [dim])
    xdim9 = dat9_size(1)
    ydim9 = dat9_size(2)
    zdim9 = dat9_size(3)
    opsDat9Cardinality = opsArg9%dim * xdim9 * ydim9 * zdim9
    dat9_base = getDatBaseFromOpsArg3D(opsArg9, start_indx, 1)
    CALL c_f_pointer(opsArg9%data, opsDat9Local, [opsDat9Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg10), dat10_size, [dim])
    xdim10 = dat10_size(1)
    ydim10 = dat10_size(2)
    zdim10 = dat10_size(3)
    opsDat10Cardinality = opsArg10%dim * xdim10 * ydim10 * zdim10
    dat10_base = getDatBaseFromOpsArg3D(opsArg10, start_indx, 1)
    CALL c_f_pointer(opsArg10%data, opsDat10Local, [opsDat10Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg11), dat11_size, [dim])
    xdim11 = dat11_size(1)
    ydim11 = dat11_size(2)
    zdim11 = dat11_size(3)
    opsDat11Cardinality = opsArg11%dim * xdim11 * ydim11 * zdim11
    dat11_base = getDatBaseFromOpsArg3D(opsArg11, start_indx, 1)
    CALL c_f_pointer(opsArg11%data, opsDat11Local, [opsDat11Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg12), dat12_size, [dim])
    xdim12 = dat12_size(1)
    ydim12 = dat12_size(2)
    zdim12 = dat12_size(3)
    opsDat12Cardinality = opsArg12%dim * xdim12 * ydim12 * zdim12
    dat12_base = getDatBaseFromOpsArg3D(opsArg12, start_indx, 1)
    CALL c_f_pointer(opsArg12%data, opsDat12Local, [opsDat12Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg13), dat13_size, [dim])
    xdim13 = dat13_size(1)
    ydim13 = dat13_size(2)
    zdim13 = dat13_size(3)
    opsDat13Cardinality = opsArg13%dim * xdim13 * ydim13 * zdim13
    dat13_base = getDatBaseFromOpsArg3D(opsArg13, start_indx, 1)
    CALL c_f_pointer(opsArg13%data, opsDat13Local, [opsDat13Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg14), dat14_size, [dim])
    xdim14 = dat14_size(1)
    ydim14 = dat14_size(2)
    zdim14 = dat14_size(3)
    opsDat14Cardinality = opsArg14%dim * xdim14 * ydim14 * zdim14
    dat14_base = getDatBaseFromOpsArg3D(opsArg14, start_indx, 1)
    CALL c_f_pointer(opsArg14%data, opsDat14Local, [opsDat14Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg15), dat15_size, [dim])
    xdim15 = dat15_size(1)
    ydim15 = dat15_size(2)
    zdim15 = dat15_size(3)
    opsDat15Cardinality = opsArg15%dim * xdim15 * ydim15 * zdim15
    dat15_base = getDatBaseFromOpsArg3D(opsArg15, start_indx, 1)
    CALL c_f_pointer(opsArg15%data, opsDat15Local, [opsDat15Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg16), dat16_size, [dim])
    xdim16 = dat16_size(1)
    ydim16 = dat16_size(2)
    zdim16 = dat16_size(3)
    opsDat16Cardinality = opsArg16%dim * xdim16 * ydim16 * zdim16
    dat16_base = getDatBaseFromOpsArg3D(opsArg16, start_indx, 1)
    CALL c_f_pointer(opsArg16%data, opsDat16Local, [opsDat16Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg17), dat17_size, [dim])
    xdim17 = dat17_size(1)
    ydim17 = dat17_size(2)
    zdim17 = dat17_size(3)
    opsDat17Cardinality = opsArg17%dim * xdim17 * ydim17 * zdim17
    dat17_base = getDatBaseFromOpsArg3D(opsArg17, start_indx, 1)
    CALL c_f_pointer(opsArg17%data, opsDat17Local, [opsDat17Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg18), dat18_size, [dim])
    xdim18 = dat18_size(1)
    ydim18 = dat18_size(2)
    zdim18 = dat18_size(3)
    opsDat18Cardinality = opsArg18%dim * xdim18 * ydim18 * zdim18
    dat18_base = getDatBaseFromOpsArg3D(opsArg18, start_indx, 1)
    CALL c_f_pointer(opsArg18%data, opsDat18Local, [opsDat18Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg19), dat19_size, [dim])
    xdim19 = dat19_size(1)
    ydim19 = dat19_size(2)
    zdim19 = dat19_size(3)
    opsDat19Cardinality = opsArg19%dim * xdim19 * ydim19 * zdim19
    dat19_base = getDatBaseFromOpsArg3D(opsArg19, start_indx, 1)
    CALL c_f_pointer(opsArg19%data, opsDat19Local, [opsDat19Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg20), dat20_size, [dim])
    xdim20 = dat20_size(1)
    ydim20 = dat20_size(2)
    zdim20 = dat20_size(3)
    opsDat20Cardinality = opsArg20%dim * xdim20 * ydim20 * zdim20
    dat20_base = getDatBaseFromOpsArg3D(opsArg20, start_indx, 1)
    CALL c_f_pointer(opsArg20%data, opsDat20Local, [opsDat20Cardinality])

    CALL c_f_pointer(getDatSizeFromOpsArg(opsArg21), dat21_size, [dim])
    xdim21 = dat21_size(1)
    ydim21 = dat21_size(2)
    zdim21 = dat21_size(3)
    opsDat21Cardinality = opsArg21%dim * xdim21 * ydim21 * zdim21
    dat21_base = getDatBaseFromOpsArg3D(opsArg21, start_indx, 1)
    CALL c_f_pointer(opsArg21%data, opsDat21Local, [opsDat21Cardinality])

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg22), opsDat22Local, [opsArg22%dim])
    dat22_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg23), opsDat23Local, [opsArg23%dim])
    dat23_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg24), opsDat24Local, [opsArg24%dim])
    dat24_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg25), opsDat25Local, [opsArg25%dim])
    dat25_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg26), opsDat26Local, [opsArg26%dim])
    dat26_base = 1

    CALL c_f_pointer(getGblPtrFromOpsArg(opsArg27), opsDat27Local, [opsArg27%dim])
    dat27_base = 1

!   ==============
!   Halo exchanges
!   ==============
#ifndef OPS_LAZY
    CALL ops_H_D_exchanges_host(opsArgArray, 27)
    CALL ops_halo_exchanges(opsArgArray, 27, range, 3)
    CALL ops_H_D_exchanges_host(opsArgArray, 27)
#endif

    CALL ops_timers_core(t2__)

!   ==============================
!   Call kernel wrapper subroutine
!   ==============================
    CALL bounds_kernel_eqAB_xdir_wrap( &
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
                        opsDat14Local, &
                        opsDat15Local, &
                        opsDat16Local, &
                        opsDat17Local, &
                        opsDat18Local, &
                        opsDat19Local, &
                        opsDat20Local, &
                        opsDat21Local, &
                        opsDat22Local, &
                        opsDat23Local, &
                        opsDat24Local, &
                        opsDat25Local, &
                        opsDat26Local, &
                        opsDat27Local, &
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
                        dat17_base, &
                        dat18_base, &
                        dat19_base, &
                        dat20_base, &
                        dat21_base, &
                        dat22_base, &
                        dat23_base, &
                        dat24_base, &
                        dat25_base, &
                        dat26_base, &
                        dat27_base, &
                        start_indx, &
                        end_indx )

    CALL ops_timers_core(t3__)

#ifndef OPS_LAZY
    CALL ops_set_dirtybit_host(opsArgArray, 27)
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
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg9, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg10, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg11, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg12, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg13, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg14, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg15, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg16, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg17, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg18, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg19, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg20, transfer)
    transfer_total = transfer_total + transfer
    CALL ops_compute_transfer(3, start_indx, end_indx, opsArg21, transfer)
    transfer_total = transfer_total + transfer

    CALL setKernelTime(356, kernelName//c_null_char, t3__-t2__, t2__-t1__, transfer_total, 0)

END SUBROUTINE

#ifdef OPS_LAZY
SUBROUTINE bounds_kernel_eqAB_xdir_host( userSubroutine, block, dim, range, &
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
    opsArg16, &
    opsArg17, &
    opsArg18, &
    opsArg19, &
    opsArg20, &
    opsArg21, &
    opsArg22, &
    opsArg23, &
    opsArg24, &
    opsArg25, &
    opsArg26, &
    opsArg27 &
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
    TYPE(ops_arg), INTENT(IN) :: opsArg17
    TYPE(ops_arg), INTENT(IN) :: opsArg18
    TYPE(ops_arg), INTENT(IN) :: opsArg19
    TYPE(ops_arg), INTENT(IN) :: opsArg20
    TYPE(ops_arg), INTENT(IN) :: opsArg21
    TYPE(ops_arg), INTENT(IN) :: opsArg22
    TYPE(ops_arg), INTENT(IN) :: opsArg23
    TYPE(ops_arg), INTENT(IN) :: opsArg24
    TYPE(ops_arg), INTENT(IN) :: opsArg25
    TYPE(ops_arg), INTENT(IN) :: opsArg26
    TYPE(ops_arg), INTENT(IN) :: opsArg27

    TYPE(ops_arg), DIMENSION(27), TARGET :: opsArgArray
    INTEGER(KIND=4) :: n_indx
    CHARACTER(LEN=40) :: namelit

    namelit = "bounds_kernel_eqAB_xdir"

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
    opsArgArray(17) = opsArg17
    opsArgArray(18) = opsArg18
    opsArgArray(19) = opsArg19
    opsArgArray(20) = opsArg20
    opsArgArray(21) = opsArg21
    opsArgArray(22) = opsArg22
    opsArgArray(23) = opsArg23
    opsArgArray(24) = opsArg24
    opsArgArray(25) = opsArg25
    opsArgArray(26) = opsArg26
    opsArgArray(27) = opsArg27

    DO n_indx = 1, 3
        range_tmp(2*n_indx-1) = range(2*n_indx-1)-1
        range_tmp(2*n_indx)   = range(2*n_indx)
    END DO

    CALL create_kerneldesc_and_enque(namelit//c_null_char, c_loc(opsArgArray), &
                                    27, 356, dim, 0, c_loc(range_tmp), &
                                    block%blockCptr, c_funloc(bounds_kernel_eqAB_xdir_host_execute))

END SUBROUTINE
#endif

END MODULE BOUNDS_KERNEL_EQAB_XDIR_MODULE
