SUBROUTINE d2fdxy(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   D2FDXY
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   15-APR-2003:  RSC MODIFIED FOR SENGA2
!   10-OCT-2004:  RSC NULL VERSION

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES SECOND XY-DERIVATIVE OF SPECIFIED FUNCTION

!   *************************************************************************

!   ARGUMENTS
!   =========
    TYPE(ops_dat) :: functn, fderiv

!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)
    integer(kind=4) :: istart,ifinis,jstart,jfinis

!   BEGIN
!   =====

!   =========================================================================

    IF (nxglbl==1 .or. nyglbl==1) THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(d2fdxy_kernel_null, "d2fdxy_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    ELSE
!       =========================================================================

!       END CONDITIONS
!       ==============

        istart = 1
        ifinis = nxglbl
        jstart = 1
        jfinis = nyglbl
        IF(nendxl == nbound) istart = 6
        IF(nendxr == nbound) ifinis = nxglbl-5
        IF(nendyl == nbound) jstart = 6
        IF(nendyr == nbound) jfinis = nyglbl-5

!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES

        rangexyz = [istart,ifinis,jstart,jfinis,1,nzglbl]
        call ops_par_loop(d2fdxy_kernel_interior, "d2fdxy_interior", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p550_m550_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!       =========================================================================

!       LH END X-DIRECTION
!       ==================
        IF(nendxl == nbound) THEN
!           TAKE SECOND XY-DERIVATIVE IN X-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT

!           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
            rangexyz = [1,1,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqA, "d2fdxy_kernel_eqA", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p420_m020_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
            rangexyz = [2,2,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqB, "d2fdxy_kernel_eqB", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p320_m120_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           LH POINT PLUS 2: 4TH ORDER CENTRED
            rangexyz = [3,3,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqC, "d2fdxy_kernel_eqC", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p220_m220_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           LH POINT PLUS 3: 6TH ORDER CENTRED
            rangexyz = [4,4,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqD, "d2fdxy_kernel_eqD", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p330_m330_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           LH POINT PLUS 4: 8TH ORDER CENTRED
            rangexyz = [5,5,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqE, "d2fdxy_kernel_eqE", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        END IF
!       =========================================================================

!       RH END X-DIRECTION
!       ==================
        IF(nendxr == nbound) THEN
!           TAKE SECOND XY-DERIVATIVE IN X-RIGHT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT

!           RH POINT MINUS 4: 8TH ORDER CENTRED
            rangexyz = [nxglbl-4,nxglbl-4,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqF, "d2fdxy_kernel_eqF", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           RH POINT MINUS 3: 6TH ORDER CENTRED
            rangexyz = [nxglbl-3,nxglbl-3,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqG, "d2fdxy_kernel_eqG", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p330_m330_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           RH POINT MINUS 2: 4TH ORDER CENTRED
            rangexyz = [nxglbl-2,nxglbl-2,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqH, "d2fdxy_kernel_eqH", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p220_m220_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
            rangexyz = [nxglbl-1,nxglbl-1,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqI, "d2fdxy_kernel_eqI", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p120_m320_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
            rangexyz = [nxglbl,nxglbl,jstart,jfinis,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqJ, "d2fdxy_kernel_eqJ", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p020_m420_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        END IF

!       =========================================================================

!       LH END Y-DIRECTION
!       ==================
        IF(nendyl == nbound) THEN
!           TAKE SECOND XY-DERIVATIVE IN Y-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [istart,ifinis,1,1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqK, "d2fdxy_kernel_eqK", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p240_m200_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,2,2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqL, "d2fdxy_kernel_eqL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p230_m210_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,3,3,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqM, "d2fdxy_kernel_eqM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p220_m220_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,4,4,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqN, "d2fdxy_kernel_eqN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p330_m330_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,5,5,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqO, "d2fdxy_kernel_eqO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))


!           LH IN X LH IN Y CORNER
!           ======================
            IF(nendxl == nbound) THEN
            rangexyz = [1,1,1,1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqP, "d2fdxy_kernel_eqP", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p440_p000_mixed_xy, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,2,2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqQ, "d2fdxy_kernel_eqQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p330_m110_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,1,2,2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqR, "d2fdxy_kernel_eqR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p430_m010_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,1,1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqS, "d2fdxy_kernel_eqS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p340_m100_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,1,1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqT, "d2fdxy_kernel_eqT", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p240_m200_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,2,2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqU, "d2fdxy_kernel_eqU", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p230_m210_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,1,3,5,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqV, "d2fdxy_kernel_eqV", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p420_m020_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,3,5,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqW, "d2fdxy_kernel_eqW", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p320_m120_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,3,5,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqX, "d2fdxy_kernel_eqX", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_idx())

            END IF

!           RH IN X LH IN Y CORNER
!           ======================
            IF(nendxr == nbound) THEN
            rangexyz = [nxglbl,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqY, "d2fdxy_kernel_eqY", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p040_m400_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,2,2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqZ, "d2fdxy_kernel_eqZ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p100_m300_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,2,2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAA, "d2fdxy_kernel_eqAA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p030_m410_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAB, "d2fdxy_kernel_eqAB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p140_m300_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAC, "d2fdxy_kernel_eqAC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p240_m200_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,2,2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAD, "d2fdxy_kernel_eqAD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p230_m210_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,3,5,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAE, "d2fdxy_kernel_eqAE", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p020_m420_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,3,5,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAF, "d2fdxy_kernel_eqAF", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p120_m320_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,3,5,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAG, "d2fdxy_kernel_eqAG", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(nxglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_idx())

            END IF

        END IF

!       =========================================================================

!       RH END Y-DIRECTION
!       ==================
        IF(nendyr == nbound) THEN
!           TAKE SECOND XY-DERIVATIVE IN Y-RIGHT INNER HALO
!           EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [istart,ifinis,nyglbl-4,nyglbl-4,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAH, "d2fdxy_kernel_eqAH", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,nyglbl-3,nyglbl-3,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAI, "d2fdxy_kernel_eqAI", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p330_m330_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,nyglbl-2,nyglbl-2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAJ, "d2fdxy_kernel_eqAJ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p220_m220_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,nyglbl-1,nyglbl-1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAK, "d2fdxy_kernel_eqAK", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p210_m230_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAL, "d2fdxy_kernel_eqAL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p200_m240_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           LH IN X RH IN Y CORNER
!           ======================
            IF(nendxl == nbound) THEN
            rangexyz = [1,1,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAM, "d2fdxy_kernel_eqAM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p400_p040_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,nyglbl-1,nyglbl-1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAN, "d2fdxy_kernel_eqAN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p300_m100_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,1,nyglbl-1,nyglbl-1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAO, "d2fdxy_kernel_eqAO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p410_p030_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAP, "d2fdxy_kernel_eqAP", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p300_m140_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,nyglbl-1,nyglbl-1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAQ, "d2fdxy_kernel_eqAQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p210_m230_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAR, "d2fdxy_kernel_eqAR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p200_m240_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,1,nyglbl-4,nyglbl-2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAS, "d2fdxy_kernel_eqAS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p420_m020_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,nyglbl-4,nyglbl-2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAT, "d2fdxy_kernel_eqAT", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p320_m120_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,nyglbl-4,nyglbl-2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAU, "d2fdxy_kernel_eqAU", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(nyglbl_ops, 1 , "integer(kind=4)", OPS_READ), &
                            ops_arg_idx())

            END IF

!           RH IN X RH IN Y CORNER
!           ======================
            IF(nendxr == nbound) THEN
            rangexyz = [nxglbl,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAV, "d2fdxy_kernel_eqAV", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p000_m440_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,nyglbl-1,nyglbl-1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAW, "d2fdxy_kernel_eqAW", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p110_m330_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,nyglbl-1,nyglbl-1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAX, "d2fdxy_kernel_eqAX", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p010_m430_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAY, "d2fdxy_kernel_eqAY", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p100_m340_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,nyglbl-1,nyglbl-1,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqAZ, "d2fdxy_kernel_eqAZ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p210_m230_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqBA, "d2fdxy_kernel_eqBA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p200_m240_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,nyglbl-4,nyglbl-2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqBB, "d2fdxy_kernel_eqBB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p020_m420_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,nyglbl-4,nyglbl-2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqBC, "d2fdxy_kernel_eqBC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p120_m320_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,nyglbl-4,nyglbl-2,1,nzglbl]
            call ops_par_loop(d2fdxy_kernel_eqBD, "d2fdxy_kernel_eqBD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p440_m440_mixed_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(nxglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_gbl(nyglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_idx())

            END IF

        END IF

!       =========================================================================

!       SCALING
!       =======
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(d2fdxy_kernel_scaling, "d2fdxy_scaling", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_RW))

    END IF

!   =========================================================================

END SUBROUTINE d2fdxy













