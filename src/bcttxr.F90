SUBROUTINE bcttxr

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCTTXR
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   30-DEC-2003:  CREATED
!   09-MAY-2015:  RSC MODIFIED FOR ISOTHERMAL WALL

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR TEMPERATURE
!   AND ITS TIME DERIVATIVE

!   X-DIRECTION RIGHT-HAND END

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

!   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!   =========================================================================

!   EVALUATE AND RETURN STRTXR,DTDTXR
    rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(bcdt_kernel_xdir, "bcdt_kernel_xdir", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                    ops_arg_dat(d_dtdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(trin, 1, "real(kind=8)", OPS_READ))

!   =========================================================================

!   ISOTHERMAL WALL
    IF(nsbcxr == nsbcw2) THEN
        rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(bcdt_kernel_xdir_eqA, "bcdt_kernel_xdir", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE),  &
                        ops_arg_dat(d_dtdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_gbl(rxrprm, nbcprr, "real(kind=8)", OPS_READ))

    END IF

!   =========================================================================

END SUBROUTINE bcttxr
