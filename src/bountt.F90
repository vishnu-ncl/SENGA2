SUBROUTINE bountt

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BOUNTT
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   29-SEP-2003:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   SYNCHRONISES THE TIME-DEPENDENT BOUNDARY CONDITIONS

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    real(kind=8) :: fornow
    integer(kind=4) :: jc,kc
    integer(kind=4) :: ispec
    integer(kind=4) :: iindex,ipower,icoef1,icoef2
    integer(kind=4) :: rangexyz(6)
    real(kind=8) :: rxlprm_1, rxrprm_1, rylprm_1, ryrprm_1, rzlprm_1, rzrprm_1

!   BEGIN
!   =====

!   =========================================================================

!   SYNCHRONISE AT CURRENT TIME STEP
!   --------------------------------
    irkstp = 1
    !KA   FIX INFLOW BC
    btime  = tstep
    fupelc = .true.
!   =========================================================================

!   X-DIRECTION LEFT-HAND END
!   -------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fxltrb) call bcutxl

!   LOCAL BC SUPPORT
    IF(fxlcnv) THEN

!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE

!       =======================================================================

        IF(nsbcxl == nsbci2) THEN

!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxl

!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttxl

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_xdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ),  &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxl

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bountt_kernel_eqF_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcxl == nsbci3) THEN

!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           SET DENSITY AND TIME DERIVATIVE
            call bcdtxl

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxl

!           CONSERVATIVE VARIABLES
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqC_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxl

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(bountt_kernel_eqE_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_stryxl(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcxl == nsbcw1)THEN

!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_xl, "bounds_kernel_eqAF_xl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_p400_x, "real(kind=8)", OPS_RW))

            DO iindex = 1,nintmx
                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(set_zero_kernel_int, "set_zero itndex", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqE_xdir, "boundt_kernel_eqE_xdir", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           NC SUGGESTION FOR WALL:
!           IF ISSUES STILL PERSISTS, COPY LINCOM LAST STATEMENT ADDED FOR
!           ENFORCING DY/DN=0,HERE
!           CONSERVATIVE VARIABLES
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_xdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcxl == nsbcw2)THEN

!           WALL BC No 2
!           NO-SLIP WALL - ISOTHERMAL

            rxlprm_1 = rxlprm(1)
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqH_xdir, "bountt_kernel_eqH_xdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dtdtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(rxlprm_1, 1, "real(kind=8)", OPS_READ))

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_xdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,1,1,nyglbl,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_xdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,1,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

    END IF
!   X-DIRECTION LEFT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   X-DIRECTION RIGHT-HAND END
!   --------------------------
    IF(fxrcnv) THEN

!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE

!       =======================================================================

        IF(nsbcxr == nsbci2) THEN

!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxr

!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttxr

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_xdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ),  &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxr

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bountt_kernel_eqF_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcxr == nsbci3) THEN

!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           SET DENSITY AND TIME DERIVATIVE
            call bcdtxr

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxr

!           CONSERVATIVE VARIABLES
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqC_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxr

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(bountt_kernel_eqE_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_stryxr(ispec), 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcxr == nsbcw1) THEN

!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_xr, "bounds_kernel_eqAF_xr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_m400_x, "real(kind=8)", OPS_RW))

            DO iindex = 1,nintmx
                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(set_zero_kernel_int, "set_zero itndex", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqE_xdir, "boundt_kernel_eqE_xdir", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_xdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcxr == nsbcw2) THEN

!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***

            rxrprm_1 = rxrprm(1)
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqH_xdir, "bountt_kernel_eqH_xdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dtdtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(rxrprm_1, 1, "real(kind=8)", OPS_READ))

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_xdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_xdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [nxglbl,nxglbl,1,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

    END IF
!   X-DIRECTION RIGHT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Y-DIRECTION LEFT-HAND END
!   -------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fyltrb) call bcutyl

!   LOCAL BC SUPPORT
    IF(fylcnv) THEN

!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE

!       =======================================================================

        IF(nsbcyl == nsbci2) THEN

!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyl

!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttyl

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_ydir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ),  &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyl

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bountt_kernel_eqF_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcyl == nsbci3) THEN

!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           SET DENSITY AND TIME DERIVATIVE
            call bcdtyl

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyl

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqC_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyl

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(bountt_kernel_eqE_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_stryyl(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcyl == nsbcw1) THEN

!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_yl, "bounds_kernel_eqAF_yl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_p040_y, "real(kind=8)", OPS_RW))

            DO iindex = 1,nintmx
                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(set_zero_kernel_int, "set_zero itndex", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqE_ydir, "boundt_kernel_eqE_ydir", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_ydir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcyl == nsbcw2)THEN

!           WALL BC No 2
!           NO-SLIP WALL - ISOTHERMAL

            rylprm_1 = rylprm(1)
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqH_ydir, "bountt_kernel_eqH_ydir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dtdtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(rylprm_1, 1, "real(kind=8)", OPS_READ))

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_ydir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,1,1,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_ydir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,1,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

    END IF

!   Y-DIRECTION LEFT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Y-DIRECTION RIGHT-HAND END
!   --------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fyrtrb) call bcutyr

!   LOCAL BC SUPPORT
    IF(fyrcnv) THEN

!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE

!       =======================================================================

        IF(nsbcyr == nsbci2) THEN

!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyr

!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttyr

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_ydir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

             call ops_par_loop(bountt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ),  &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyr

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bountt_kernel_eqF_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcyr == nsbci3) THEN

!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           SET DENSITY AND TIME DERIVATIVE
            call bcdtyr

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyr

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqC_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyr

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(bountt_kernel_eqE_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_stryyr(ispec), 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbcyr == nsbcw1) THEN

!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_yr, "bounds_kernel_eqAF_yr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_m040_y, "real(kind=8)", OPS_RW))

            DO iindex = 1,nintmx
                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(set_zero_kernel_int, "set_zero itndex", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqE_ydir, "boundt_kernel_eqE_ydir", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_ydir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcyr == nsbcw2) THEN

!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***

            ryrprm_1 = ryrprm(1)
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqH_ydir, "bountt_kernel_eqH_ydir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dtdtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(ryrprm_1, 1, "real(kind=8)", OPS_READ))

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_ydir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_ydir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,nyglbl,nyglbl,1,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

    END IF
!   Y-DIRECTION RIGHT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Z-DIRECTION LEFT-HAND END
!   -------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fzltrb) call bcutzl

!   LOCAL BC SUPPORT
    IF(fzlcnv) THEN

!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE

!       =======================================================================

        IF(nsbczl == nsbci2) THEN

!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzl

!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttzl

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_zdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bountt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

             call ops_par_loop(bountt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ),  &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzl

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bountt_kernel_eqF_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbczl == nsbci3) THEN

!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           SET DENSITY AND TIME DERIVATIVE
            call bcdtzl

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzl

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bountt_kernel_eqC_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzl

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(bountt_kernel_eqE_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_stryzl(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbczl == nsbcw1) THEN

!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bounds_kernel_eqAF_zl, "bounds_kernel_eqAF_zl", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_p004_z, "real(kind=8)", OPS_RW))

            DO iindex = 1,nintmx
                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(set_zero_kernel_int, "set_zero itndex", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(boundt_kernel_eqE_zdir, "boundt_kernel_eqE_zdir", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(boundt_kernel_eqG_zdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbczl == nsbcw2)THEN

!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***

            rzlprm_1 = rzlprm(1)
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bountt_kernel_eqH_zdir, "bountt_kernel_eqH_zdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dtdtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(rzlprm_1, 1, "real(kind=8)", OPS_READ))

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_zdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,1,nyglbl,1,1]
                call ops_par_loop(boundt_kernel_eqG_zdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,1,1]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

    END IF
!   Z-DIRECTION LEFT-HAND END

!   =========================================================================
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   =========================================================================

!   Z-DIRECTION RIGHT-HAND END
!   --------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fzrtrb) call bcutzr

!   LOCAL BC SUPPORT
    IF(fzrcnv) THEN

!       =======================================================================

!       OUTFLOW BOUNDARY CONDITIONS
!       ---------------------------

!       OUTFLOW BC No 1
!       SUBSONIC NON-REFLECTING OUTFLOW
!       WITH OPTION TO SET PRESSURE AT INFINITY
!       REQUIRES NO ACTION HERE

!       =======================================================================

!       INFLOW BOUNDARY CONDITIONS
!       --------------------------

!       INFLOW BC No 1
!       SUBSONIC NON-REFLECTING LAMINAR INFLOW
!       REQUIRES NO ACTION HERE

!       =======================================================================

        IF(nsbczr == nsbci2) THEN

!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzr

!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttzr

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_zdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bountt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            call ops_par_loop(bountt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ),  &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzr

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bountt_kernel_eqF_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ), &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbczr == nsbci3) THEN

!           INFLOW BC No 3
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED DENSITY

!           SET DENSITY AND TIME DERIVATIVE
            call bcdtzr

!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzr

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bountt_kernel_eqC_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_derr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzr

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(bountt_kernel_eqE_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_yerr(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                                ops_arg_dat(d_stryzr(ispec), 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ))

            END DO

        END IF

!       =======================================================================

        IF(nsbczr == nsbcw1) THEN

!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bounds_kernel_eqAF_zr, "bounds_kernel_eqAF_zr", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_trun, 1, s3d_000_to_m004_z, "real(kind=8)", OPS_RW))

            DO iindex = 1,nintmx
                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(set_zero_kernel_int, "set_zero itndex", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(boundt_kernel_eqE_zdir, "boundt_kernel_eqE_zdir", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_zdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbczr == nsbcw2) THEN

!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
!           *** RSC 10-APRIL-2005 CODING CHECKED BUT BC UNTESTED ***

            rzrprm_1 = rzrprm(1)
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bountt_kernel_eqH_zdir, "bountt_kernel_eqH_zdir", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_dat(d_dtdtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_WRITE), &
                            ops_arg_gbl(rzrprm_1, 1, "real(kind=8)", OPS_READ))

!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            DO iindex = 1,nintmx
                call ops_par_loop(set_zero_kernel_int, "set_zero", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_WRITE))
            END DO

            DO ispec = 1,nspec
!               SET THE TEMPERATURE INTERVAL INDEX
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1

                call ops_par_loop(boundt_kernel_eqE_zdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ), &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_INC),  &
                                ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ntint, nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ipower, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
            END DO

!           CONSERVATIVE VARIABLES
            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bountt_kernel_eqG_xyz, "bountt_kernel_eqG_xyz", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_uerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_verr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_werr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(kind=8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
                call ops_par_loop(boundt_kernel_eqG_zdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ),  &
                                ops_arg_dat(d_itndex(iindex), 1, s3d_000, "integer(kind=4)", OPS_READ),  &
                                ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(kind=8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer(kind=4)", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer(kind=4)", OPS_READ))

            END DO

            rangexyz = [1,nxglbl,1,nyglbl,nzglbl,nzglbl]
            call ops_par_loop(bountt_kernel_eqD, "init values", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erun, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_eerr, 1, s3d_000, "real(kind=8)", OPS_WRITE),  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_READ))

        END IF

!       =======================================================================

    END IF
!   Z-DIRECTION RIGHT-HAND END

!   =========================================================================

END SUBROUTINE bountt
