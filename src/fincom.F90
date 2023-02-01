SUBROUTINE fincom

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga
 
!   *************************************************************************

!   FINCOM
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   15-JAN-2003:  CREATED
!   08-AUG-2012:  RSC EVALUATE ALL SPECIES

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   COMPUTES FINAL SOLUTION VALUES IN ERK SCHEME
!   BY DOING A LINEAR COMBINATION OF LEFT- AND RIGHT-HAND SIDES

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------


!   LOCAL DATA
!   ==========
    integer :: ispec
    integer :: rangexyz(6)

!   -------------------------------------------------------------------------

!   BEGIN
!   =====

!   =========================================================================

!   FINAL ERK SUBSTEP
!    =================

!   -------------------------------------------------------------------------
!   NOTE: ALL ERK ERROR ARRAYS ARE INITIALISED TO ZERO IN SUBROUTINE ADAPTT
!   -------------------------------------------------------------------------

!   DENSITY
!   ----------
    rangexyz(1) = 1
    IF (nsbcxl == nsbci3) rangexyz(1) = 2
    rangexyz(2) = nxglbl
    IF (nsbcxr == nsbci3) rangexyz(2) = nxglbl-1
    rangexyz(3) = 1
    IF (nsbcyl == nsbci3) rangexyz(3) = 2
    rangexyz(4) = nyglbl
    IF (nsbcyr == nsbci3) rangexyz(4) = nyglbl-1
    rangexyz(5) = 1
    IF (nsbczl == nsbci3) rangexyz(5) = 2
    rangexyz(6) = nzglbl
    IF (nsbczr == nsbci3) rangexyz(6) = nzglbl-1

    call ops_par_loop(fincom_kernel_main, "fincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_derr, 1, s3d_000, "real(8)", OPS_INC), &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_gbl(rkerr(nrkstp), 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(rklhs(nrkstp), 1, "real(8)", OPS_READ))

!   -------------------------------------------------------------------------
!   U VELOCITY
!   ----------
    rangexyz(1) = 1
    IF (nsbcxl == nsbci2 .or. nsbcxl == nsbci3 .or. nsbcxl == nsbcw1 .or. nsbcxl == nsbcw2) rangexyz(1) = 2
    rangexyz(2) = nxglbl
    IF (nsbcxr == nsbci2 .or. nsbcxr == nsbci3 .or. nsbcxr == nsbcw1 .or. nsbcxr == nsbcw2) rangexyz(2) = nxglbl-1
    rangexyz(3) = 1
    IF (nsbcyl == nsbci2 .or. nsbcyl == nsbci3 .or. nsbcyl == nsbcw1 .or. nsbcyl == nsbcw2) rangexyz(3) = 2
    rangexyz(4) = nyglbl
    IF (nsbcyr == nsbci2 .or. nsbcyr == nsbci3 .or. nsbcyr == nsbcw1 .or. nsbcyr == nsbcw2) rangexyz(4) = nyglbl-1
    rangexyz(5) = 1
    IF (nsbczl == nsbci2 .or. nsbczl == nsbci3 .or. nsbczl == nsbcw1 .or. nsbczl == nsbcw2) rangexyz(5) = 2
    rangexyz(6) = nzglbl
    IF (nsbczr == nsbci2 .or. nsbczr == nsbci3 .or. nsbczr == nsbcw1 .or. nsbczr == nsbcw2) rangexyz(6) = nzglbl-1

    call ops_par_loop(fincom_kernel_main, "fincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_uerr, 1, s3d_000, "real(8)", OPS_INC), &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_gbl(rkerr(nrkstp), 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(rklhs(nrkstp), 1, "real(8)", OPS_READ))

!   -------------------------------------------------------------------------
!   V-VELOCITY
!   ----------
    call ops_par_loop(fincom_kernel_main, "fincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_verr, 1, s3d_000, "real(8)", OPS_INC), &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_gbl(rkerr(nrkstp), 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(rklhs(nrkstp), 1, "real(8)", OPS_READ))

!   -------------------------------------------------------------------------
!   W-VELOCITY
!   ----------
    call ops_par_loop(fincom_kernel_main, "fincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_werr, 1, s3d_000, "real(8)", OPS_INC), &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_gbl(rkerr(nrkstp), 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(rklhs(nrkstp), 1, "real(8)", OPS_READ))

!   -------------------------------------------------------------------------
!   STAGNATION INTERNAL ENERGY
!   --------------------------
    rangexyz(1) = 1
    IF (nsbcxl == nsbci2 .or. nsbcxl == nsbcw2) rangexyz(1) = 2
    rangexyz(2) = nxglbl
    IF (nsbcxr == nsbci2 .or. nsbcxr == nsbcw2) rangexyz(2) = nxglbl-1
    rangexyz(3) = 1
    IF (nsbcyl == nsbci2 .or. nsbcyl == nsbcw2) rangexyz(3) = 2
    rangexyz(4) = nyglbl
    IF (nsbcyr == nsbci2 .or. nsbcyr == nsbcw2) rangexyz(4) = nyglbl-1
    rangexyz(5) = 1
    IF (nsbczl == nsbci2 .or. nsbczl == nsbcw2) rangexyz(5) = 2
    rangexyz(6) = nzglbl
    IF (nsbczr == nsbci2 .or. nsbczr == nsbcw2) rangexyz(6) = nzglbl-1

    call ops_par_loop(fincom_kernel_main, "fincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_eerr, 1, s3d_000, "real(8)", OPS_INC), &
                    ops_arg_dat(d_erun, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_gbl(rkerr(nrkstp), 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(rklhs(nrkstp), 1, "real(8)", OPS_READ))

!   -------------------------------------------------------------------------
!   SPECIES MASS FRACTIONS
!   ----------------------
!   RSC 08-AUG-2012 EVALUATE ALL SPECIES
!   DO ISPEC = 1,NSPM1
    rangexyz(1) = 1
    IF (nsbcxl == nsbci2 .or. nsbcxl == nsbci3) rangexyz(1) = 2
    rangexyz(2) = nxglbl
    IF (nsbcxr == nsbci2 .or. nsbcxr == nsbci3) rangexyz(2) = nxglbl-1
    rangexyz(3) = 1
    IF (nsbcyl == nsbci2 .or. nsbcyl == nsbci3) rangexyz(3) = 2
    rangexyz(4) = nyglbl
    IF (nsbcyr == nsbci2 .or. nsbcyr == nsbci3) rangexyz(4) = nyglbl-1
    rangexyz(5) = 1
    IF (nsbczl == nsbci2 .or. nsbczl == nsbci3) rangexyz(5) = 2
    rangexyz(6) = nzglbl
    IF (nsbczr == nsbci2 .or. nsbczr == nsbci3) rangexyz(6) = nzglbl-1

    DO ispec = 1,nspec

        call ops_par_loop(fincom_kernel_MD, "fincom mulit-dim", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_yerr, 9, s3d_000, "real(8)", OPS_INC), &
                        ops_arg_dat(d_yrun, 9, s3d_000, "real(8)", OPS_RW), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_RW), &
                        ops_arg_gbl(rkerr(irkstp), 1, "real(8)", OPS_READ), &
                        ops_arg_gbl(rklhs(irkstp), 1, "real(8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

    END DO

!   -------------------------------------------------------------------------

END SUBROUTINE fincom
