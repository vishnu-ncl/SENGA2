SUBROUTINE d2fdyz(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
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
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

    IF (nyglbl == 1) THEN
        rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
        call ops_par_loop(d2fdyz_kernel_null, "d2fdyz_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

    ELSE
!       =========================================================================

!       INTERIOR SCHEME
!       ===============

!       TENTH ORDER EXPLICIT DIFFERENCES
        rangexyz = (/1,nxglbl,6,nyglbl-5,6,nzglbl-5/)
        call ops_par_loop(d2fdyz_kernel_interior, "d2fdyz_interior", senga_grid, 3, rangexyz, &
                        ops_arg_dat(functn, 1, s3d_p055_m055_mixed_yz, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       LH END Y-DIRECTION
!       ==================
!           TAKE SECOND YZ-DERIVATIVE IN Y-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = (/1,nxglbl,1,1,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_lh_ydir_4th_onesided, "d2fdyz_lh_ydir_4th_onesided", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p042_m002_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,2,2,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_lh_ydir_4th_mixed, "d2fdyz_lh_ydir_4th_mixed", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p032_m012_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,3,3,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_lh_ydir_4th_centred, "d2fdyz_lh_ydir_4th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p022_m022_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,4,4,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_lh_ydir_6th_centred, "d2fdyz_lh_ydir_6th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m033_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,5,5,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_lh_ydir_8th_centred, "d2fdyz_lh_ydir_8th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       RH END Y-DIRECTION
!       ==================
!           TAKE SECOND YZ-DERIVATIVE IN Y-RIGHT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-4,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_rh_ydir_8th_centred, "d2fdyz_rh_ydir_8th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-3,nyglbl-3,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_rh_ydir_6th_centred, "d2fdyz_rh_ydir_6th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m033_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-2,nyglbl-2,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_rh_ydir_4th_centred, "d2fdyz_rh_ydir_4th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p022_m022_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-1,nyglbl-1,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_rh_ydir_4th_mixed, "d2fdyz_rh_ydir_4th_mixed", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p012_m032_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl,nyglbl,6,nzglbl-5/)
            call ops_par_loop(d2fdyz_kernel_rh_ydir_4th_onesided, "d2fdyz_rh_ydir_4th_onesided", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m042_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       LH END Z-DIRECTION
!       ==================
!           TAKE SECOND XZ-DERIVATIVE IN Z-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = (/1,nxglbl,6,nyglbl-5,1,1/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_4th_onesided, "d2fdyz_lh_zdir_4th_onesided", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p024_m020_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,6,nyglbl-5,2,2/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_4th_mixed, "d2fdyz_lh_zdir_4th_mixed", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p023_m021_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,6,nyglbl-5,3,3/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_4th_centred, "d2fdyz_lh_zdir_4th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p022_m022_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,6,nyglbl-5,4,4/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_6th_centred, "d2fdyz_lh_zdir_6th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m033_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,6,nyglbl-5,5,5/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_8th_centred, "d2fdyz_lh_zdir_8th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           LH IN Y LH IN Z CORNER
!           ======================
            rangexyz = (/1,nxglbl,1,1,1,1/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqA, "d2fdyz_lh_zdir_corner_eqA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_p000_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,2,2,2,2/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqB, "d2fdyz_lh_zdir_corner_eqB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m011_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,1,2,2/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqC, "d2fdyz_lh_zdir_corner_eqC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p043_m001_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,2,2,1,1/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqD, "d2fdyz_lh_zdir_corner_eqD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p034_m010_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,3,5,1,1/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqE, "d2fdyz_lh_zdir_corner_eqE", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p024_m020_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,3,5,2,2/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqF, "d2fdyz_lh_zdir_corner_eqF", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p023_m021_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,1,3,5/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqG, "d2fdyz_lh_zdir_corner_eqG", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p042_m002_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,2,2,3,5/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqH, "d2fdyz_lh_zdir_corner_eqH", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p032_m012_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,3,5,3,5/)


            rangexyz = (/1,nxglbl,4,5,4,5/)


            rangexyz = (/1,nxglbl,5,5,5,5/)


!           RH IN Y LH IN Z CORNER
!           ======================
            rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,1/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqL, "d2fdyz_lh_zdir_corner_eqL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p004_m040_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-1,nyglbl-1,2,2/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqM, "d2fdyz_lh_zdir_corner_eqM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p010_m030_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl,nyglbl,2,2/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqN, "d2fdyz_lh_zdir_corner_eqN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p003_m041_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-1,nyglbl-1,1,1/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqO, "d2fdyz_lh_zdir_corner_eqO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p014_m030_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-2,1,1/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqP, "d2fdyz_lh_zdir_corner_eqP", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p024_m020_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-2,2,2/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqQ, "d2fdyz_lh_zdir_corner_eqQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p023_m021_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl,nyglbl,3,5/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqR, "d2fdyz_lh_zdir_corner_eqR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m042_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-1,nyglbl-1,3,5/)
            call ops_par_loop(d2fdyz_kernel_lh_zdir_corner_eqS, "d2fdyz_lh_zdir_corner_eqS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p012_m032_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-2,3,5/)


            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-3,4,5/)


            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-4,5,5/)


!       =========================================================================

!       RH END Z-DIRECTION
!       ==================
!           TAKE SECOND YZ-DERIVATIVE IN Z-RIGHT INNER HALO
!           EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = (/1,nxglbl,6,nyglbl-5,nzglbl-4,nzglbl-4/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_8th_centred, "d2fdyz_rh_zdir_8th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p044_m044_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,6,nyglbl-5,nzglbl-3,nzglbl-3/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_6th_centred, "d2fdyz_rh_zdir_6th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p033_m033_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,6,nyglbl-5,nzglbl-2,nzglbl-2/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_4th_centred, "d2fdyz_rh_zdir_4th_centred", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p022_m022_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,6,nyglbl-5,nzglbl-1,nzglbl-1/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_4th_mixed, "d2fdyz_rh_zdir_4th_mixed", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p021_m023_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,6,nyglbl-5,nzglbl,nzglbl/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_4th_onesided, "d2fdyz_rh_zdir_4th_onesided", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p020_m024_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           LH IN Y RH IN Z CORNER
!           ======================
            rangexyz = (/1,nxglbl,1,1,nzglbl,nzglbl/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqA, "d2fdyz_rh_zdir_corner_eqA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p040_p004_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,2,2,nzglbl-1,nzglbl-1/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqB, "d2fdyz_rh_zdir_corner_eqB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p030_m010_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,1,nzglbl-1,nzglbl-1/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqC, "d2fdyz_rh_zdir_corner_eqC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p041_p003_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,2,2,nzglbl,nzglbl/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqD, "d2fdyz_rh_zdir_corner_eqD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p030_m014_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,3,5,nzglbl-1,nzglbl-1/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqE, "d2fdyz_rh_zdir_corner_eqE", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p021_m023_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,3,5,nzglbl,nzglbl/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqF, "d2fdyz_rh_zdir_corner_eqF", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p020_m024_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,1,1,nzglbl-4,nzglbl-2/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqG, "d2fdyz_rh_zdir_corner_eqG", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p042_m002_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,2,2,nzglbl-4,nzglbl-2/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqH, "d2fdyz_rh_zdir_corner_eqH", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p032_m012_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,3,5,nzglbl-4,nzglbl-2/)


            rangexyz = (/1,nxglbl,4,5,nzglbl-4,nzglbl-3/)


            rangexyz = (/1,nxglbl,5,5,nzglbl-4,nzglbl-4/)


!           RH IN Y RH IN Z CORNER
!           ======================
            rangexyz = (/1,nxglbl,nyglbl,nyglbl,nzglbl,nzglbl/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqL, "d2fdyz_rh_zdir_corner_eqL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p000_m044_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-1,nyglbl-1,nzglbl-1,nzglbl-1/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqM, "d2fdyz_rh_zdir_corner_eqM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p011_m033_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl,nyglbl,nzglbl-1,nzglbl-1/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqN, "d2fdyz_rh_zdir_corner_eqN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p001_m043_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-1,nyglbl-1,nzglbl,nzglbl/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqO, "d2fdyz_rh_zdir_corner_eqO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p010_m034_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-2,nzglbl-1,nzglbl-1/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqP, "d2fdyz_rh_zdir_corner_eqP", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p021_m023_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-2,nzglbl,nzglbl/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqQ, "d2fdyz_rh_zdir_corner_eqQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p020_m024_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl,nyglbl,nzglbl-4,nzglbl-2/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqR, "d2fdyz_rh_zdir_corner_eqR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m042_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-1,nyglbl-1,nzglbl-4,nzglbl-2/)
            call ops_par_loop(d2fdyz_kernel_rh_zdir_corner_eqS, "d2fdyz_rh_zdir_corner_eqS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p012_m032_mixed_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-2,nzglbl-4,nzglbl-2/)


            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-3,nzglbl-4,nzglbl-3/)


            rangexyz = (/1,nxglbl,nyglbl-4,nyglbl-4,nzglbl-4,nzglbl-4/)


!       =========================================================================

!       SCALING
!       =======
        rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
        call ops_par_loop(d2fdyz_kernel_scaling, "d2fdyz_scaling", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_RW))
    END IF

!   =========================================================================

END SUBROUTINE d2fdyz
