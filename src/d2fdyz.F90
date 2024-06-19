SUBROUTINE d2fdyz(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   D2FDYZ
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   01-AUG-1996:  CREATED
!   15-MAY-2003:  RSC MODIFIED FOR SENGA2
!   10-OCT-2004:  RSC NULL VERSION
!   31-DEC-2006:  RSC BUG FIX INDICES

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES SECOND YZ-DERIVATIVE OF SPECIFIED FUNCTION

!   *************************************************************************

!   ARGUMENTS
!   =========
    TYPE(ops_dat) :: functn, fderiv

!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)
    integer(kind=4) :: jstart,jfinis,kstart,kfinis

!   BEGIN
!   =====

!   =========================================================================

    IF (nyglbl==1 .or. nzglbl==1) THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(d2fdyz_kernel_null, "d2fdyz_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    ELSE
!       =========================================================================

!       END CONDITIONS
!       ==============

        jstart = 1
        jfinis = nyglbl
        kstart = 1
        kfinis = nzglbl
        IF(nendyl == nbound) jstart = 6
        IF(nendyr == nbound) jfinis = nyglbl-5
        IF(nendzl == nbound) kstart = 6
        IF(nendzr == nbound) kfinis = nzglbl-5

!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz = [1,nxglbl,jstart,jfinis,kstart,kfinis]
        call ops_par_loop(d2fdyz_kernel_interior, "d2fdyz_interior", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p055_m055_mixed_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!       =========================================================================

!       LH END Y-DIRECTION
!       ==================
        IF(nendyl == nbound) THEN
!           TAKE SECOND YZ-DERIVATIVE IN Y-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [1,nxglbl,1,1,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqA, "d2fdyz_kernel_eqA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p042_m002_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,2,2,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqB, "d2fdyz_kernel_eqB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p032_m012_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,3,3,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqC, "d2fdyz_kernel_eqC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p022_m022_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,4,4,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqD, "d2fdyz_kernel_eqD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m033_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,5,5,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqE, "d2fdyz_kernel_eqE", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        END IF

!       =========================================================================

!       RH END Y-DIRECTION
!       ==================
        IF(nendyr == nbound) THEN
!           TAKE SECOND YZ-DERIVATIVE IN Y-RIGHT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [1,nxglbl,nyglbl-4,nyglbl-4,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqF, "d2fdyz_kernel_eqF", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-3,nyglbl-3,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqG, "d2fdyz_kernel_eqG", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m033_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-2,nyglbl-2,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqH, "d2fdyz_kernel_eqH", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p022_m022_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-1,nyglbl-1,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqI, "d2fdyz_kernel_eqI", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p012_m032_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdyz_kernel_eqJ, "d2fdyz_kernel_eqJ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m042_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

        END IF

!       =========================================================================

!       LH END Z-DIRECTION
!       ==================
        IF(nendzl == nbound) THEN
!           TAKE SECOND XZ-DERIVATIVE IN Z-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [1,nxglbl,jstart,jfinis,1,1]
            call ops_par_loop(d2fdyz_kernel_eqK, "d2fdyz_kernel_eqK", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p024_m020_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,jstart,jfinis,2,2]
            call ops_par_loop(d2fdyz_kernel_eqL, "d2fdyz_kernel_eqL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p023_m021_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,jstart,jfinis,3,3]
            call ops_par_loop(d2fdyz_kernel_eqM, "d2fdyz_kernel_eqM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p022_m022_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,jstart,jfinis,4,4]
            call ops_par_loop(d2fdyz_kernel_eqN, "d2fdyz_kernel_eqN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m033_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,jstart,jfinis,5,5]
            call ops_par_loop(d2fdyz_kernel_eqO, "d2fdyz_kernel_eqO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           LH IN Y LH IN Z CORNER
!           ======================
            IF(nendyl == nbound) THEN
            rangexyz = [1,nxglbl,1,1,1,1]
            call ops_par_loop(d2fdyz_kernel_eqP, "d2fdyz_kernel_eqP", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_p000_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,2,2,2,2]
            call ops_par_loop(d2fdyz_kernel_eqQ, "d2fdyz_kernel_eqQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m011_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,1,1,2,2]
            call ops_par_loop(d2fdyz_kernel_eqR, "d2fdyz_kernel_eqR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p043_m001_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,2,2,1,1]
            call ops_par_loop(d2fdyz_kernel_eqS, "d2fdyz_kernel_eqS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p034_m010_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,3,5,1,1]
            call ops_par_loop(d2fdyz_kernel_eqT, "d2fdyz_kernel_eqT", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p024_m020_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,3,5,2,2]
            call ops_par_loop(d2fdyz_kernel_eqU, "d2fdyz_kernel_eqU", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p023_m021_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,1,1,3,5]
            call ops_par_loop(d2fdyz_kernel_eqV, "d2fdyz_kernel_eqV", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p042_m002_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,2,2,3,5]
            call ops_par_loop(d2fdyz_kernel_eqW, "d2fdyz_kernel_eqW", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p032_m012_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,3,5,3,5]
            call ops_par_loop(d2fdyz_kernel_eqX, "d2fdyz_kernel_eqX", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_idx())

            END IF

!           RH IN Y LH IN Z CORNER
!           ======================
            IF(nendyr == nbound) THEN
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,1]
            call ops_par_loop(d2fdyz_kernel_eqY, "d2fdyz_kernel_eqY", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p004_m040_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-1,nyglbl-1,2,2]
            call ops_par_loop(d2fdyz_kernel_eqZ, "d2fdyz_kernel_eqZ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p010_m030_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl,nyglbl,2,2]
            call ops_par_loop(d2fdyz_kernel_eqAA, "d2fdyz_kernel_eqAA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p003_m041_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-1,nyglbl-1,1,1]
            call ops_par_loop(d2fdyz_kernel_eqAB, "d2fdyz_kernel_eqAB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p014_m030_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-4,nyglbl-2,1,1]
            call ops_par_loop(d2fdyz_kernel_eqAC, "d2fdyz_kernel_eqAC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p024_m020_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-4,nyglbl-2,2,2]
            call ops_par_loop(d2fdyz_kernel_eqAD, "d2fdyz_kernel_eqAD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p023_m021_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl,nyglbl,3,5]
            call ops_par_loop(d2fdyz_kernel_eqAE, "d2fdyz_kernel_eqAE", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m042_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-1,nyglbl-1,3,5]
            call ops_par_loop(d2fdyz_kernel_eqAF, "d2fdyz_kernel_eqAF", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p012_m032_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-4,nyglbl-2,3,5]
            call ops_par_loop(d2fdyz_kernel_eqAG, "d2fdyz_kernel_eqAG", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(nyglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_idx())

            END IF

        END IF

!       =========================================================================

!       RH END Z-DIRECTION
!       ==================
        IF(nendzr == nbound) THEN
!           TAKE SECOND YZ-DERIVATIVE IN Z-RIGHT INNER HALO
!           EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [1,nxglbl,jstart,jfinis,nzglbl-4,nzglbl-4]
            call ops_par_loop(d2fdyz_kernel_eqAH, "d2fdyz_kernel_eqAH", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,jstart,jfinis,nzglbl-3,nzglbl-3]
            call ops_par_loop(d2fdyz_kernel_eqAI, "d2fdyz_kernel_eqAI", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m033_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,jstart,jfinis,nzglbl-2,nzglbl-2]
            call ops_par_loop(d2fdyz_kernel_eqAJ, "d2fdyz_kernel_eqAJ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p022_m022_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,jstart,jfinis,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdyz_kernel_eqAK, "d2fdyz_kernel_eqAK", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p021_m023_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,jstart,jfinis,nzglbl,nzglbl]
            call ops_par_loop(d2fdyz_kernel_eqAL, "d2fdyz_kernel_eqAL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p020_m024_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!           LH IN Y RH IN Z CORNER
!           ======================
            IF(nendyl == nbound)THEN
            rangexyz = [1,nxglbl,1,1,nzglbl,nzglbl]
            call ops_par_loop(d2fdyz_kernel_eqAM, "d2fdyz_kernel_eqAM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p040_p004_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,2,2,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdyz_kernel_eqAN, "d2fdyz_kernel_eqAN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p030_m010_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,1,1,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdyz_kernel_eqAO, "d2fdyz_kernel_eqAO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p041_p003_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,2,2,nzglbl,nzglbl]
            call ops_par_loop(d2fdyz_kernel_eqAP, "d2fdyz_kernel_eqAP", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p030_m014_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,3,5,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdyz_kernel_eqAQ, "d2fdyz_kernel_eqAQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p021_m023_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,3,5,nzglbl,nzglbl]
            call ops_par_loop(d2fdyz_kernel_eqAR, "d2fdyz_kernel_eqAR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p020_m024_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,1,1,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdyz_kernel_eqAS, "d2fdyz_kernel_eqAS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p042_m002_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,2,2,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdyz_kernel_eqAT, "d2fdyz_kernel_eqAT", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p032_m012_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,3,5,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdyz_kernel_eqAU, "d2fdyz_kernel_eqAU", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(nzglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_idx())

            END IF

!           RH IN Y RH IN Z CORNER
!           ======================
            IF(nendyr == nbound) THEN
            rangexyz = [1,nxglbl,nyglbl,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdyz_kernel_eqAV, "d2fdyz_kernel_eqAV", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p000_m044_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-1,nyglbl-1,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdyz_kernel_eqAW, "d2fdyz_kernel_eqAW", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p011_m033_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdyz_kernel_eqAX, "d2fdyz_kernel_eqAX", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p001_m043_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-1,nyglbl-1,nzglbl,nzglbl]
            call ops_par_loop(d2fdyz_kernel_eqAY, "d2fdyz_kernel_eqAY", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p010_m034_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-4,nyglbl-2,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdyz_kernel_eqAZ, "d2fdyz_kernel_eqAZ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p021_m023_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-4,nyglbl-2,nzglbl,nzglbl]
            call ops_par_loop(d2fdyz_kernel_eqBA, "d2fdyz_kernel_eqBA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p020_m024_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdyz_kernel_eqBB, "d2fdyz_kernel_eqBB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m042_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-1,nyglbl-1,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdyz_kernel_eqBC, "d2fdyz_kernel_eqBC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p012_m032_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE))

            rangexyz = [1,nxglbl,nyglbl-4,nyglbl-2,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdyz_kernel_eqBD, "d2fdyz_kernel_eqBD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(nyglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_gbl(nzglbl_ops, 1, "integer(kind=4)", OPS_READ), &
                            ops_arg_idx())

            END IF

        END IF

!       =========================================================================

!       SCALING
!       =======
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(d2fdyz_kernel_scaling, "d2fdyz_scaling", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(kind=8)", OPS_RW))
    END IF

!   =========================================================================

END SUBROUTINE d2fdyz
