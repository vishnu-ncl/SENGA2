SUBROUTINE flamin

    use OPS_Fortran_Reference
    use OPS_Fortran_hdf5_Declarations
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

!   PARAMETERS
!   ==========
!   ESTIMATED FLAME LOCATION AND THICKNESS
    real(kind=8), parameter :: clocat = 0.0025_8, cthick = 0.0005_8

!   FUNCTION
!   ========
    real(kind=8) :: erfunc
    EXTERNAL erfunc

!   LOCAL DATA
!   ==========
    real(kind=8) :: u0
    real(kind=8) :: deltagx,deltagy,deltagz
    
    real(kind=8) :: angfrx,angfry,angfrz
    real(kind=8), parameter :: cpsi = 3.0_8, cfo = 0.0556_8
    real(kind=8) :: rpsi,rgspec_ispec

    integer(kind=4) :: ispec,il
    integer(kind=4) :: rangexyz(6)

    real(kind=4) :: xr(nl),tr(nl)
    real(kind=4) :: yrunr(nl,nspcmx)
    real(kind=8) :: xl(nl),theta(nl)
    real(kind=8) :: yrunl(nl,nspcmx)

    character(len=60) :: fname
    character(len=3) :: pnxhdf
    parameter(pnxhdf = '.h5')
    fname = 'output/timestep_flamin'//pnxhdf

!   BEGIN
!   =====

!   =========================================================================

!   SPECIFY INITIAL THERMOCHEMICAL FIELD HERE
!   =========================================

!   SET PRODUCT TEMPERATURE
!   -----------------------
!   REACTANT TEMPERATURE SET IN CONTROL FILE
    u0 = 4.0_8
    rpsi = 1.0_8/8.0_8
    rpsi = rpsi*xgdlen

!   Global indexing
!   ---------------
    deltagx = xgdlen/(REAL(nxglbl-1,kind=8))
    deltagy = ygdlen/(REAL(nyglbl-1,kind=8))
    deltagz = zgdlen/(REAL(nzglbl-1,kind=8))

!   SET THE VELOCITY PROFILE FOR TGV
!   --------------------------------
    angfrx = 8.0_8*atan(1.0_8)/xgdlen
    angfry = 8.0_8*atan(1.0_8)/ygdlen
    angfrz = 8.0_8*atan(1.0_8)/zgdlen

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(flamin_kernel_set_velocity_tgv, "SET THE VELOCITY PROFILE FOR TGV",  senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_psi, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_gbl(prin, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(drin, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(u0, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(deltagx, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(deltagy, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(deltagz, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(angfrx, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(angfry, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(angfrz, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(xgdlen, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(cpsi, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(rpsi, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_idx())

!   SET SPECIES MASS FRACTION PROFILES
!   ----------------------------------
    DO ispec = 1,nspcmx
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                          ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE))
    END DO

!    OPEN(UNIT=14,FILE='TGV_species.dat',STATUS='unknown')
!    DO il = 1, nl
!!       Read xr(il) and the nspec values for yrunr(il,:)
!        READ(14,*) xr(il), (yrunr(il,ispec), ispec=1,nspec)
!
!!       Convert to REAL if needed (often not strictly necessary in modern compilers)
!        xl(il) = REAL(xr(il),kind=8)
!        DO ispec=1,nspec
!            yrunl(il,ispec) = REAL(yrunr(il,ispec),kind=8)
!        END DO
!    END DO
!    CLOSE(14)

!    DO ispec = 1,nspcmx
!        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
!        call ops_par_loop(copy_kernel_sdim_to_mdim, "A_multidim(ispec) = B", senga_grid, 3, rangexyz,  &
!                        ops_arg_dat(d_yrhs_mdim, 9, s3d_000, "real(kind=8)", OPS_WRITE), &
!                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
!                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
!    END DO
!
!    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
!    call ops_par_loop(flamin_kernel_eqA, "flamin_kernel_eqA", senga_grid, 3, rangexyz,  &
!                    ops_arg_dat(d_yrhs_mdim, 9, s3d_000, "real(kind=8)", OPS_RW), &
!                    ops_arg_gbl(yrunl, nl*nspcmx, "real(kind=8)", OPS_READ), &
!                    ops_arg_gbl(xl, nl, "real(kind=8)", OPS_READ), &
!                    ops_arg_gbl(deltagx, 1, "real(kind=8)", OPS_READ), &
!                    ops_arg_gbl(nspcmx, 1, "integer(kind=4)", OPS_READ), &
!                    ops_arg_idx())
!
!    DO ispec = 1,nspcmx
!        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
!        call ops_par_loop(copy_kernel_mdim_to_sdim, "A = B_multidim(ispec)", senga_grid, 3, rangexyz,  &
!                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
!                        ops_arg_dat(d_yrhs_mdim, 9, s3d_000, "real(kind=8)", OPS_READ), &
!                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))
!    END DO

    DO ispec = 1,nspcmx
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_yrun_dump(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
    END DO

    OPEN(UNIT=16,FILE='t_initial.dat',STATUS='unknown')
    DO il = 1,nl
        READ(16,*)xr(il),tr(il)
    END DO
    CLOSE(16)

    DO il=1,nl
        xl(il) = REAL(xr(il),kind=8)
        theta(il) = REAL(tr(il),kind=8)
    END DO

!   SET TEMPERATURE PROFILE
!   -----------------------
    call ops_par_loop(flamin_kernel_eqD, "flamin_kernel_eqD",  senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_RW), &
                    ops_arg_gbl(xl, nl, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(theta, nl, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(deltagx, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_idx())

!   SET DENSITY PROFILE ASSUMING CONSTANT PRESSURE
!   -------------------
!   PRESSURE SET IN CONTROL FILE
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    DO ispec = 1,nspcmx
        rgspec_ispec = rgspec(ispec)
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(flamin_kernel_eqB, "flamin_kernel_eqB", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rgspec_ispec, 1, "real(kind=8)", OPS_READ))
    END DO

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(flamin_kernel_eqE, "flamin_kernel_eqE", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(prin, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(trin, 1, "real(kind=8)", OPS_READ))

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(flamin_kernel_eqF, "flamin_kernel_eqF",  senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(prin, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(u0, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(deltagx, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(deltagy, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(deltagz, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(angfrx, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(angfry, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_gbl(angfrz, 1, "real(kind=8)", OPS_READ), &
                    ops_arg_idx())

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(flamin_kernel_eqC, "flamin_kernel_eqC", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!   =========================================================================

    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d2prun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_prun, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_par_loop(copy_kernel, "copy", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d2trun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ))

    call ops_fetch_block_hdf5_file(senga_grid, trim(fname))

    DO ispec = 1,nspcmx
        call ops_fetch_dat_hdf5_file(d_yrun(ispec), trim(fname))
    END DO

    call ops_fetch_dat_hdf5_file(d2trun, trim(fname))
    call ops_fetch_dat_hdf5_file(d2prun, trim(fname))

9000  FORMAT(a)

END SUBROUTINE flamin
