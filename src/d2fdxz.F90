SUBROUTINE d2fdxz(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   D2FDXZ
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   15-MAY-2003:  RSC MODIFIED FOR SENGA2
!   10-OCT-2004:  RSC NULL VERSION

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES SECOND XZ-DERIVATIVE OF SPECIFIED FUNCTION

!   *************************************************************************

!   ARGUMENTS
!   =========
    TYPE(ops_dat) :: functn, fderiv

!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)
    integer(kind=4) :: istart,ifinis,kstart,kfinis

!   BEGIN
!   =====

!   =========================================================================

    IF (nxglbl==1 .or. nzglbl==1) THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(d2fdxz_kernel_null, "d2fdxz_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    ELSE
!       =========================================================================

!       END CONDITIONS
!       ==============

        istart = 1
        ifinis = nxglbl
        kstart = 1
        kfinis = nzglbl
        IF(nendxl == nbound) istart = 6
        IF(nendxr == nbound) ifinis = nxglbl-5
        IF(nendzl == nbound) kstart = 6
        IF(nendzr == nbound) kfinis = nzglbl-5

!       INTERIOR SCHEME
!       ===============
!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz = [istart,ifinis,1,nyglbl,kstart,kfinis]
        call ops_par_loop(d2fdxz_kernel_interior, "d2fdxz_interior", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p505_m505_mixed_xz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!       =========================================================================

!       LH END X-DIRECTION
!       ==================
        IF(nendxl == nbound) THEN
!           TAKE SECOND XZ-DERIVATIVE IN X-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [1,1,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqA, "d2fdxz_kernel_eqA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p402_m002_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqB, "d2fdxz_kernel_eqB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p302_m102_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,3,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqC, "d2fdxz_kernel_eqC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p202_m202_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [4,4,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqD, "d2fdxz_kernel_eqD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m303_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [5,5,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqE, "d2fdxz_kernel_eqE", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        END IF

!       =========================================================================

!       RH END X-DIRECTION
!       ==================
        IF(nendxr == nbound) THEN
!           TAKE SECOND XZ-DERIVATIVE IN X-RIGHT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [nxglbl-4,nxglbl-4,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqF, "d2fdxz_kernel_eqF", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-3,nxglbl-3,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqG, "d2fdxz_kernel_eqG", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m303_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-2,nxglbl-2,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqH, "d2fdxz_kernel_eqH", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p202_m202_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqI, "d2fdxz_kernel_eqI", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p102_m302_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_eqJ, "d2fdxz_kernel_eqJ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m402_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        END IF

!       =========================================================================

!       LH END Z-DIRECTION
!       ==================
        IF(nendxl == nbound) THEN
!           TAKE SECOND XZ-DERIVATIVE IN Z-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [istart,ifinis,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_eqK, "d2fdxz_kernel_eqK", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p204_m200_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_eqL, "d2fdxz_kernel_eqL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p203_m201_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,3,3]
            call ops_par_loop(d2fdxz_kernel_eqM, "d2fdxz_kernel_eqM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p202_m202_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,4,4]
            call ops_par_loop(d2fdxz_kernel_eqN, "d2fdxz_kernel_eqN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m303_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,5,5]
            call ops_par_loop(d2fdxz_kernel_eqO, "d2fdxz_kernel_eqO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           LH IN X LH IN Z CORNER
!           ======================
            IF(nendxl == nbound) THEN
            rangexyz = [1,1,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_eqP, "d2fdxz_kernel_eqP", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_p000_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_eqQ, "d2fdxz_kernel_eqQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m101_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,1,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_eqR, "d2fdxz_kernel_eqR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p403_m001_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_eqS, "d2fdxz_kernel_eqS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p304_m100_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_eqT, "d2fdxz_kernel_eqT", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p204_m200_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_eqU, "d2fdxz_kernel_eqU", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p203_m201_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,1,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_eqV, "d2fdxz_kernel_eqV", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p402_m002_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_eqW, "d2fdxz_kernel_eqW", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p302_m102_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_eqX, "d2fdxz_kernel_eqX", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_idx())

            END IF

!           RH IN X LH IN Z CORNER
!           ======================
            IF(nendxr == nbound)THEN
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_eqY, "d2fdxz_kernel_eqY", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p004_m400_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_eqZ, "d2fdxz_kernel_eqZ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p100_m300_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_eqAA, "d2fdxz_kernel_eqAA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p003_m401_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_eqAB, "d2fdxz_kernel_eqAB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p104_m300_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_eqAC, "d2fdxz_kernel_eqAC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p204_m200_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_eqAD, "d2fdxz_kernel_eqAD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p203_m201_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_eqAE, "d2fdxz_kernel_eqAE", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m402_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_eqAF, "d2fdxz_kernel_eqAF", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p102_m302_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_eqAG, "d2fdxz_kernel_eqAG", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(nxglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_idx())

            END IF

        END IF

!       =========================================================================

!       RH END Z-DIRECTION
!       ==================
        IF(nendzr == nbound) THEN
!           TAKE SECOND XZ-DERIVATIVE IN Z-RIGHT INNER HALO
!           EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [istart,ifinis,1,nyglbl,nzglbl-4,nzglbl-4]
            call ops_par_loop(d2fdxz_kernel_eqAH, "d2fdxz_kernel_eqAH", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,nzglbl-3,nzglbl-3]
            call ops_par_loop(d2fdxz_kernel_eqAI, "d2fdxz_kernel_eqAI", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m303_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,nzglbl-2,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_eqAJ, "d2fdxz_kernel_eqAJ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p202_m202_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_eqAK, "d2fdxz_kernel_eqAK", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p201_m203_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_eqAL, "d2fdxz_kernel_eqAL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p200_m204_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           LH IN X RH IN Z CORNER
!           ======================
            IF(nendxl == nbound) THEN
            rangexyz = [1,1,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_eqAM, "d2fdxz_kernel_eqAM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p400_p004_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_eqAN, "d2fdxz_kernel_eqAN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p300_m100_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,1,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_eqAO, "d2fdxz_kernel_eqAO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p401_p003_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_eqAP, "d2fdxz_kernel_eqAP", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p300_m104_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_eqAQ, "d2fdxz_kernel_eqAQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p201_m203_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_eqAR, "d2fdxz_kernel_eqAR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p200_m204_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,1,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_eqAS, "d2fdxz_kernel_eqAS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p402_m002_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_eqAT, "d2fdxz_kernel_eqAT", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p302_m102_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_eqAU, "d2fdxz_kernel_eqAU", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(nzglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_idx())

            END IF

!           RH IN X RH IN Z CORNER
!           ======================
            IF(nendxr == nbound) THEN
            rangexyz = [nxglbl,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_eqAV, "d2fdxz_kernel_eqAV", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p000_m404_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_eqAW, "d2fdxz_kernel_eqAW", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p101_m303_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_eqAX, "d2fdxz_kernel_eqAX", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p001_m403_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_eqAY, "d2fdxz_kernel_eqAY", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p100_m304_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_eqAZ, "d2fdxz_kernel_eqAZ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p201_m203_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_eqBA, "d2fdxz_kernel_eqBA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p200_m204_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_eqBB, "d2fdxz_kernel_eqBB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m402_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_eqBC, "d2fdxz_kernel_eqBC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p102_m302_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_eqBD, "d2fdxz_kernel_eqBD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(nxglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_gbl(nzglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_idx())

            END IF

        END IF

!       =========================================================================

!       SCALING
!       =======
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(d2fdxz_kernel_scaling, "d2fdxz_scaling", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_RW))

    END IF

!   =========================================================================

END SUBROUTINE d2fdxz
