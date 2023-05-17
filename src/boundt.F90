SUBROUTINE boundt

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   BOUNDT
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   29-SEP-2003:  CREATED
!   10-MAY-2015:  RSC WALL BCS UPDATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   APPLIES BOUNDARY CONDITIONS TO PRIMITIVE VARIABLES

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    real(kind=8) :: fornow
    integer :: jc,kc
    integer :: ispec
    integer :: iindex,ipower,icoef1,icoef2
    integer :: itint,icp
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================
!KA   FIX INFLOW BC
!KA   RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
    btime  = rktim(irkstp)
    fupelc = .false.

!   X-DIRECTION LEFT-HAND END
!   -------------------------

!   GLOBAL BC SUPPORT
!   TURBULENT INFLOW VELOCITY FIELD
    IF(fxltrb) call bcutxl

!   LOCAL BC SUPPORT
    IF(fxlcnv)THEN
  
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
            rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_xdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxl
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
      
!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase
      
                rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqF_xdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_RW),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_stryxl, 9, s3d_000_strid3d_yz, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
            rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqC_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxl
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqD_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_WRITE), &
                                ops_arg_dat(d_stryxl, 9, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                                ops_arg_dat(d_strdxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbcxl == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxl
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

        END IF

!       =======================================================================

        IF(nsbcxl == nsbcw2) THEN

!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxl
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttxl
    
!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_xdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/1,1,1,nyglbl,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqG_xdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_strtxl, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
  
        IF(nsbcxr == nsbci2)THEN
    
!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttxr
    
!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_xdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqF_xdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_RW),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_stryxr, 9, s3d_000_strid3d_yz, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
            rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqC_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytxr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqD_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_WRITE), &
                                ops_arg_dat(d_stryxr, 9, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                                ops_arg_dat(d_strdxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbcxr == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxr
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbcxr == nsbcw2) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutxr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttxr
    
!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_xdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_xdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase
      
                rangexyz = (/nxglbl,nxglbl,1,nyglbl,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqG_xdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_strtxr, 1, s3d_000_strid3d_yz, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
            rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_ydir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyl

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqF_ydir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_RW),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_stryyl, 9, s3d_000_strid3d_xz, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
            rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqC_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyl
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqD_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_WRITE), &
                                ops_arg_dat(d_stryyl, 9, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                                ops_arg_dat(d_strdyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbcyl == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyl
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbcyl == nsbcw2)THEN
    
!           WALL BC No 2
!           NO-SLIP WALL - ISOTHERMAL
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyl
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttyl
    
!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_ydir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/1,nxglbl,1,1,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqG_ydir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_strtyl, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
            rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_ydir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))
    
!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqF_ydir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_RW),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_stryyr, 9, s3d_000_strid3d_xz, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
            rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqC_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytyr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
                 call ops_par_loop(boundt_kernel_eqD_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_WRITE), &
                                ops_arg_dat(d_stryyr, 9, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                                ops_arg_dat(d_strdyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbcyr == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyr
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

        END IF
  
!       =======================================================================
  
        IF(nsbcyr == nsbcw2) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutyr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttyr
    
!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_ydir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_ydir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/1,nxglbl,nyglbl,nyglbl,1,nzglbl/)
                call ops_par_loop(boundt_kernel_eqG_ydir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_strtyr, 1, s3d_000_strid3d_xz, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
            rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
            call ops_par_loop(boundt_kernel_eqE_zdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
            call ops_par_loop(boundt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzl

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
                call ops_par_loop(boundt_kernel_eqF_zdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_RW),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_stryzl, 9, s3d_000_strid3d_xy, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
            rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
            call ops_par_loop(boundt_kernel_eqC_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzl
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
                call ops_par_loop(boundt_kernel_eqD_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_WRITE), &
                                ops_arg_dat(d_stryzl, 9, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                                ops_arg_dat(d_strdzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbczl == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzl
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
            call ops_par_loop(boundt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))
    
        END IF
  
!       =======================================================================
  
        IF(nsbczl == nsbcw2) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzl
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttzl
    
!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
            call ops_par_loop(boundt_kernel_eqE_zdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
            call ops_par_loop(boundt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/1,nxglbl,1,nyglbl,1,1/)
                call ops_par_loop(boundt_kernel_eqG_zdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_strtzl, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
  
        IF(nsbczr == nsbci2)THEN
    
!           INFLOW BC No 2
!           SUBSONIC REFLECTING INFLOW WITH SPECIFIED TEMPERATURE
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttzr
    
!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_zdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzr

!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
                call ops_par_loop(boundt_kernel_eqF_zdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_RW),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_stryzr, 9, s3d_000_strid3d_xy, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

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
            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
            call ops_par_loop(boundt_kernel_eqC_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_WRITE), &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

!           SET MASS FRACTIONS AND TIME DERIVATIVES
            call bcytzr
    
!           CONSERVATIVE VARIABLES
            DO ispec = 1,nspec
                rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
                call ops_par_loop(boundt_kernel_eqD_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_WRITE), &
                                ops_arg_dat(d_stryzr, 9, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                                ops_arg_dat(d_strdzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ))

            END DO

        END IF

!       =======================================================================
  
        IF(nsbczr == nsbcw1) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ADIABATIC
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzr
    
!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))
    
        END IF
  
!       =======================================================================
  
        IF(nsbczr == nsbcw2) THEN
    
!           WALL BC No 1
!           NO-SLIP WALL - ISOTHERMAL
    
!           SET VELOCITY COMPONENTS AND TIME DERIVATIVES
            call bcutzr
    
!           SET TEMPERATURE AND TIME DERIVATIVE
            call bcttzr
    
!           SET TEMPERATURE INTERVAL INDEX
            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
            call ops_par_loop(boundt_kernel_eqE_zdir, "SET TEMPERATURE INTERVAL INDEX", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_RW),  &
                            ops_arg_gbl(tinthi, ntinmx*nspcmx, "real(8)", OPS_READ), &
                            ops_arg_gbl(ntint, nspcmx, "integer", OPS_READ))

!           CONSERVATIVE VARIABLES
            rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
            call ops_par_loop(boundt_kernel_eqA_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_urhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_vrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_wrhs, 1, s3d_000, "real(8)", OPS_WRITE),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ), &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

            call ops_par_loop(boundt_kernel_eqB_zdir, "CONSERVATIVE VARIABLES", senga_grid, 3, rangexyz,  &
                            ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_RW),  &
                            ops_arg_dat(d_drhs, 1, s3d_000, "real(8)", OPS_READ),  &
                            ops_arg_dat(d_struzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strvzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ), &
                            ops_arg_dat(d_strwzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ))

            DO ispec = 1,nspec

!               TEMPERATURE INTERVAL INDEXING
                iindex = 1 + (ispec-1)/nspimx
                ipower = ispec - (iindex-1)*nspimx - 1
                icoef2 = ntbase**ipower
                icoef1 = icoef2*ntbase

                rangexyz = (/1,nxglbl,1,nyglbl,nzglbl,nzglbl/)
                call ops_par_loop(boundt_kernel_eqG_zdir, "TEMPERATURE INTERVAL INDEXING", senga_grid, 3, rangexyz,  &
                                ops_arg_dat(d_erhs, 1, s3d_000, "real(8)", OPS_INC),  &
                                ops_arg_dat(d_yrhs, 9, s3d_000, "real(8)", OPS_READ),  &
                                ops_arg_dat(d_itndex, 2, s3d_000, "integer", OPS_READ),  &
                                ops_arg_dat(d_strtzr, 1, s3d_000_strid3d_xy, "real(8)", OPS_READ),  &
                                ops_arg_gbl(amasch, ncofmx*ntinmx*nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(rgspec, nspcmx, "real(8)", OPS_READ), &
                                ops_arg_gbl(ncpoly, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncpom1, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ncenth, ntinmx*nspcmx, "integer", OPS_READ), &
                                ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                                ops_arg_gbl(iindex, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef1, 1, "integer", OPS_READ), &
                                ops_arg_gbl(icoef2, 1, "integer", OPS_READ))

            END DO

        END IF

    END IF
!   Z-DIRECTION RIGHT-HAND END

!   =========================================================================

END SUBROUTINE boundt
