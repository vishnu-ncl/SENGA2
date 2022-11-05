SUBROUTINE lincom

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   LINCOM
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
!   COMPUTES INTERMEDIATE SOLUTION VALUES IN ERK SCHEME
!   BY DOING LINEAR COMBINATIONS OF LEFT- AND RIGHT-HAND SIDES

!   *************************************************************************


!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------


!   LOCAL DATA
!   ==========
    real(kind=dp) :: fornow
    integer :: ispec
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   ERK SUBSTEP
!   ===========

!   -------------------------------------------------------------------------
!   NOTE: ALL ERK ERROR ARRAYS ARE INITIALISED TO ZERO IN SUBROUTINE ADAPTT
!   -------------------------------------------------------------------------

!   DENSITY
!   -------
    rangexyz = (/istald,istold,jstald,jstold,kstald,kstold/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_derr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!   -------------------------------------------------------------------------
!   U-VELOCITY
!   ----------
    rangexyz = (/istalu,istolu,jstalu,jstolu,kstalu,kstolu/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!   -------------------------------------------------------------------------
!   V-VELOCITY
!   ----------
    rangexyz = (/istalv,istolv,jstalv,jstolv,kstalv,kstolv/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!   -------------------------------------------------------------------------
!   W-VELOCITY
!   ----------
    rangexyz = (/istalw,istolw,jstalw,jstolw,kstalw,kstolw/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!   -------------------------------------------------------------------------
!   STAGNATION INTERNAL ENERGY
!   --------------------------
    rangexyz = (/istale,istole,jstale,jstole,kstale,kstole/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!   -------------------------------------------------------------------------
!   SPECIES MASS FRACTIONS
!   ----------------------
!   RSC 08-AUG-2012 EVALUATE ALL SPECIES
!   DO ISPEC = 1,NSPM1
    DO ispec = 1,nspec

        rangexyz = (/istaly,istoly,jstaly,jstoly,kstaly,kstoly/)
        call ops_par_loop(lincom_kernel_MD, "lincom mulit-dim", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_yerr, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrun, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

    END DO

!   -------------------------------------------------------------------------

END SUBROUTINE lincom
