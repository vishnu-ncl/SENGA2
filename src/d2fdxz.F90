SUBROUTINE d2fdxz(functn,fderiv)

use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
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
    integer(4) :: rangexyz(6)
    integer(4) :: istart,ifinis,kstart,kfinis

!   BEGIN
!   =====

!   =========================================================================

    IF (nxglbl==1 .or. nzglbl==1) THEN
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(d2fdxz_kernel_null, "d2fdxz_null", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

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
                        ops_arg_dat(functn, 1, s3d_p505_m505_mixed_xz, "real(8)", OPS_READ), &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!       =========================================================================

!       LH END X-DIRECTION
!       ==================
        IF(nendxl == nbound) THEN
!           TAKE SECOND XZ-DERIVATIVE IN X-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [1,1,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_lh_xdir_4th_onesided, "d2fdxz_lh_xdir_4th_onesided", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p402_m002_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_lh_xdir_4th_mixed, "d2fdxz_lh_xdir_4th_mixed", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p302_m102_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [3,3,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_lh_xdir_4th_centered, "d2fdxz_lh_xdir_4th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p202_m202_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [4,4,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_lh_xdir_6th_centered, "d2fdxz_lh_xdir_6th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m303_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [5,5,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_lh_xdir_8th_centered, "d2fdxz_lh_xdir_8th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

        END IF

!       =========================================================================

!       RH END X-DIRECTION
!       ==================
        IF(nendxr == nbound) THEN
!           TAKE SECOND XZ-DERIVATIVE IN X-RIGHT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [nxglbl-4,nxglbl-4,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_rh_xdir_8th_centered, "d2fdxz_rh_xdir_8th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-3,nxglbl-3,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_rh_xdir_6th_centered, "d2fdxz_rh_xdir_6th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m303_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-2,nxglbl-2,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_rh_xdir_4th_centered, "d2fdxz_rh_xdir_4th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p202_m202_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_rh_xdir_4th_mixed, "d2fdxz_rh_xdir_4th_mixed", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p102_m302_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,kstart,kfinis]
            call ops_par_loop(d2fdxz_kernel_rh_xdir_4th_onesided, "d2fdxz_rh_xdir_4th_onesided", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m402_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))            

        END IF

!       =========================================================================

!       LH END Z-DIRECTION
!       ==================
        IF(nendxl == nbound) THEN
!           TAKE SECOND XZ-DERIVATIVE IN Z-LEFT INNER HALO
!           EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
            rangexyz = [istart,ifinis,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_4th_onesided, "d2fdxz_lh_zdir_4th_onesided", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p204_m200_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_4th_mixed, "d2fdxz_lh_zdir_4th_mixed", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p203_m201_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,3,3]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_4th_centered, "d2fdxz_lh_zdir_4th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p202_m202_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,4,4]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_6th_centered, "d2fdxz_lh_zdir_6th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m303_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,5,5]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_8th_centered, "d2fdxz_lh_zdir_8th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           LH IN X LH IN Z CORNER
!           ======================
            IF(nendxl == nbound) THEN
            rangexyz = [1,1,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqA, "d2fdxz_lh_zdir_corner_eqA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_p000_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))            
            
            rangexyz = [2,2,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqB, "d2fdxz_lh_zdir_corner_eqB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m101_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [1,1,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqC, "d2fdxz_lh_zdir_corner_eqC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p403_m001_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))    

            rangexyz = [2,2,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqD, "d2fdxz_lh_zdir_corner_eqD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p304_m100_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqE, "d2fdxz_lh_zdir_corner_eqE", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p204_m200_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqF, "d2fdxz_lh_zdir_corner_eqF", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p203_m201_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [1,1,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqG, "d2fdxz_lh_zdir_corner_eqG", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p402_m002_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqH, "d2fdxz_lh_zdir_corner_eqH", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p302_m102_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqI, "d2fdxz_lh_zdir_corner_eqI", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_idx())

            END IF

!           RH IN X LH IN Z CORNER
!           ======================
            IF(nendxr == nbound)THEN
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqL, "d2fdxz_lh_zdir_corner_eqL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p004_m400_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqM, "d2fdxz_lh_zdir_corner_eqM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p100_m300_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqN, "d2fdxz_lh_zdir_corner_eqN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p003_m401_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqO, "d2fdxz_lh_zdir_corner_eqO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p104_m300_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,1,1]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqP, "d2fdxz_lh_zdir_corner_eqP", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p204_m200_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,2,2]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqQ, "d2fdxz_lh_zdir_corner_eqQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p203_m201_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqR, "d2fdxz_lh_zdir_corner_eqR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m402_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqS, "d2fdxz_lh_zdir_corner_eqS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p102_m302_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,3,5]
            call ops_par_loop(d2fdxz_kernel_lh_zdir_corner_eqT, "d2fdxz_lh_zdir_corner_eqT", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_gbl(nxglbl, 1, "integer(4)", OPS_READ), &
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
            call ops_par_loop(d2fdxz_kernel_rh_zdir_8th_centered, "d2fdxz_rh_zdir_8th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,nzglbl-3,nzglbl-3]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_6th_centered, "d2fdxz_rh_zdir_6th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p303_m303_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,nzglbl-2,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_4th_centered, "d2fdxz_rh_zdir_4th_centered", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p202_m202_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_4th_mixed, "d2fdxz_rh_zdir_4th_mixed", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p201_m203_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [istart,ifinis,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_4th_onesided, "d2fdxz_rh_zdir_4th_onesided", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p200_m204_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

!           LH IN X RH IN Z CORNER
!           ======================
            IF(nendxl == nbound) THEN
            rangexyz = [1,1,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqA, "d2fdxz_rh_zdir_corner_eqA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p400_p004_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqB, "d2fdxz_rh_zdir_corner_eqB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p300_m100_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [1,1,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqC, "d2fdxz_rh_zdir_corner_eqC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p401_p003_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqD, "d2fdxz_rh_zdir_corner_eqD", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p300_m104_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqE, "d2fdxz_rh_zdir_corner_eqE", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p201_m203_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqF, "d2fdxz_rh_zdir_corner_eqF", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p200_m204_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [1,1,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqG, "d2fdxz_rh_zdir_corner_eqG", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p402_m002_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [2,2,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqH, "d2fdxz_rh_zdir_corner_eqH", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p302_m102_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [3,5,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqI, "d2fdxz_rh_zdir_corner_eqI", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_gbl(nzglbl, 1, "integer(4)", OPS_READ), &
                            ops_arg_idx())

            END IF

!           RH IN X RH IN Z CORNER
!           ======================
            IF(nendxr == nbound) THEN
            rangexyz = [nxglbl,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqL, "d2fdxz_rh_zdir_corner_eqL", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p000_m404_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqM, "d2fdxz_rh_zdir_corner_eqM", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p101_m303_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqN, "d2fdxz_rh_zdir_corner_eqN", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p001_m403_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqO, "d2fdxz_rh_zdir_corner_eqO", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p100_m304_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,nzglbl-1,nzglbl-1]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqP, "d2fdxz_rh_zdir_corner_eqP", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p201_m203_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqQ, "d2fdxz_rh_zdir_corner_eqQ", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p200_m204_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl,nxglbl,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqR, "d2fdxz_rh_zdir_corner_eqR", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p002_m402_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-1,nxglbl-1,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqS, "d2fdxz_rh_zdir_corner_eqS", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p102_m302_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE))

            rangexyz = [nxglbl-4,nxglbl-2,1,nyglbl,nzglbl-4,nzglbl-2]
            call ops_par_loop(d2fdxz_kernel_rh_zdir_corner_eqT, "d2fdxz_rh_zdir_corner_eqT", senga_grid, 3, rangexyz, &
                            ops_arg_dat(functn, 1, s3d_p404_m404_mixed_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_gbl(nxglbl, 1, "integer(4)", OPS_READ), &
                            ops_arg_gbl(nzglbl, 1, "integer(4)", OPS_READ), &
                            ops_arg_idx())

            END IF

        END IF

!       =========================================================================

!       SCALING
!       =======
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(d2fdxz_kernel_scaling, "d2fdxz_scaling", senga_grid, 3, rangexyz, &
                        ops_arg_dat(fderiv, 1, s3d_000, "real(8)", OPS_RW))

    END IF

!   =========================================================================

END SUBROUTINE d2fdxz
