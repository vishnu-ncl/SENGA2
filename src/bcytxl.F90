SUBROUTINE bcytxl

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BCYTXL
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   30-DEC-2003:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
!   AND THEIR TIME DERIVATIVES

!   X-DIRECTION LEFT-HAND END

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    integer(kind=4) :: ispec
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

!   =========================================================================

!   EVALUATE AND RETURN STRYXL,DYDTXL
    DO ispec = 1,nspec
        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(bcyt_kernel_xdir_eqA, "bcyt_kernel_xdir_eqA", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dydtxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(yrin, nspcmx, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END DO

!   VM: SYNTHETIC SCALAR INFLOW
!   VM: NXLPRM(2)=1 IMPLIES THAT THE SCALAR SYTHETIC DIGITAL FILTERING IS ON
    IF ((nxlprm(2)==1).AND.(nxlprm(1)==4).AND.(ngbcxl==12)) THEN
        DO ispec = 1,nspec
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bcyt_kernel_xdir_eqB, "bcyt_kernel_xdir_eqB", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dydtxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_yinf2(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_RW), &
                            ops_arg_dat(d_yinf1(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(yrin, nspcmx, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(tstep, 1, "real(kind=8)", OPS_READ), &
                            ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
        END DO

        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(set_zero_kernel_xdir, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_totyxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE))
        DO ispec = 1,nspec-1
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bcyt_kernel_xdir_eqC, "bcyt_kernel_xdir_eqC", senga_grid, 3, rangexyz, &
                            ops_arg_dat(d_totyxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_INC), &
                            ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))
        END DO

        rangexyz = [1,1,1,nyglbl,1,nzglbl]
        call ops_par_loop(bcyt_kernel_xdir_eqD, "bcyt_kernel_xdir_eqD", senga_grid, 3, rangexyz, &
                        ops_arg_dat(d_stryxl(nspec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_dydtxl(nspec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_yinf2(nspec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_yinf1(nspec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_totyxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(tstep, 1, "real(kind=8)", OPS_READ))

    END IF

!   =========================================================================

END SUBROUTINE bcytxl
