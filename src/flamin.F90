SUBROUTINE flamin

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   FLAMIN
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   28-DEC-2003:  CREATED
!   08-JAN-2005:  RSC INITIAL 1D LAMINAR FLAME PROFILE

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   SETS INITIAL THERMOCHEMICAL FIELD
!   1D LAMINAR FLAME PROFILE (LEFT OR RIGHT FACING)
!   SPECIAL FOR 21 STEP HYDROGEN MECHAMISM

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   POPULATE DATA FROM THE FLAMIN INDATA GENERATED FROM STANDALONE ROUTINE
!   ======================================================================
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_drun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_urun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_vrun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_wrun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_trun_dump, 1, s3d_000, "real(kind=8)", OPS_READ))

    DO ispec = 1,nspcmx
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_yrun_dump(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
    END DO

!   =========================================================================

END SUBROUTINE flamin
