SUBROUTINE flamin
 
    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
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
    real(kind=8) :: clocat,cthick
    PARAMETER(clocat = 0.0025_8, cthick = 0.0005_8)

!C     PINCH OF HYDROGEN ATOM
!      real(kind=8) HPINCH,HLOCAT,HTHICK
!      PARAMETER(HPINCH = 1.0D-10, HLOCAT = 2.5D-3, HTHICK = 1.0D-4)
!C     PINCH OF HYDROGEN MOLECULE
!      real(kind=8) H2PNCH,H2LOCT,H2THCK
!      PARAMETER(H2PNCH = 1.0D-6, H2LOCT = 2.5D-3, H2THCK = 2.5D-4)

!   FUNCTION
!   ========
    real(kind=8) :: erfunc
    EXTERNAL erfunc

!   LOCAL DATA
!   ==========
    real(kind=8) :: yrinr(nspcmx),yrinp(nspcmx)
    real(kind=8) :: trinr,trinp
    real(kind=8) :: deltag,xcoord,argmnt
    real(kind=8) :: flxmas
    integer :: icproc
    integer :: ix
    integer :: ic,jc,kc
    integer :: ispec
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   SPECIFY INITIAL THERMOCHEMICAL FIELD HERE
!   =========================================

!   SET PRODUCT TEMPERATURE
!   -----------------------
!   REACTANT TEMPERATURE SET IN CONTROL FILE
    trinr = trin
!   TRINP = 2330.96554
    trinp = 2200.0_8

!   SET SPECIES MASS FRACTIONS
!   --------------------------
!   OVERRIDE MASS FRACTION VALUES SET IN CONTROL FILE

!   REACTANTS
    yrinr(1) = 0.0199886_8     !2.8312571D-2
    yrinr(2) = 0.2286239_8     !2.26500566D-1
    DO ispec = 3,nspm1
        yrinr(ispec) = zero
    END DO

    yrinr(nspec) = zero
    DO ispec = 1,nspm1
        yrinr(nspec) = yrinr(nspec) + yrinr(ispec)
    END DO
    yrinr(nspec) = one - yrinr(nspec)

!   PRODUCTS
    yrinp(1) = zero
    yrinp(2) = 0.0685323_8     !ZERO
    yrinp(3) = 0.1798974_8     !2.54716981D-1
    DO ispec = 4,nspm1
        yrinp(ispec) = zero
    END DO
    yrinp(nspec) = zero
    DO ispec = 1,nspm1
        yrinp(nspec) = yrinp(nspec) + yrinp(ispec)
    END DO
    yrinp(nspec) = one - yrinp(nspec)

!   WRITE TO REPORT FILE
    IF(iproc == 0) THEN
  
        WRITE(ncrept,*)
        WRITE(ncrept,*)'FLAMIN: reactant mass fractions:'
        DO ispec = 1,nspec
            WRITE(ncrept,'(I5,1PE15.7)')ispec,yrinr(ispec)
        END DO
        WRITE(ncrept,*)
  
        WRITE(ncrept,*)'FLAMIN: product mass fractions:'
        DO ispec = 1,nspec
            WRITE(ncrept,'(I5,1PE15.7)')ispec,yrinp(ispec)
        END DO
        WRITE(ncrept,*)
  
        WRITE(ncrept,*)'FLAMIN: reactant and product temperatures:'
        WRITE(ncrept,'(2(1PE15.7))')trinr,trinp
        WRITE(ncrept,*)
  
    END IF

!   GLOBAL INDEXING
!   ---------------
    deltag = xgdlen/(REAL(nxglbl-1,kind=8))

!   SET REACTION PROGRESS VARIABLE PROFILE
!   --------------------------------------
!   SIMPLE 1D LEFT-FACING ERROR FUNCTION PROFILE
    rangexyz = (/1,nxglbl,1,1,1,1/)
    call ops_par_loop(flamin_kernel_set_reaction_var, "SIMPLE 1D LEFT-FACING ERROR FUNCTION PROFILE", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_crin, 1, s3d_000_strid3d_x, "real(8)", OPS_WRITE), &
                    ops_arg_gbl(deltag, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(clocat, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(cthick, 1, "real(8)", OPS_READ), &
                    ops_arg_idx())

!   SET SPECIES MASS FRACTION PROFILES
!   ----------------------------------
    DO ispec = 1, nspm1
        rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
        call ops_par_loop(flamin_kernel_eqA, "SET SPECIES MASS FRACTION PROFILES", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrun, 2, s3d_000, "real(8)", OPS_WRITE), &
                        ops_arg_dat(d_crin, 1, s3d_000_strid3d_x, "real(8)", OPS_READ), &
                        ops_arg_gbl(yrinr(ispec), 1, "real(8)", OPS_READ), &
                        ops_arg_gbl(yrinp(ispec), 1, "real(8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

    END DO

    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(set_zero_kernel_MD, "set zero multi-dim", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_yrun, 2, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_gbl(nspec, 1, "integer", OPS_READ))

    DO ispec = 1, nspm1
        rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
        call ops_par_loop(flamin_kernel_eqB, "A_multidim = A_multidim + A_multidim", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_yrun, 2, s3d_000, "real(8)", OPS_RW), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ), &
                        ops_arg_gbl(nspec, 1, "integer", OPS_READ))

    END DO

    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(flamin_kernel_eqC, "A_multidim = 1.0 - A_multidim", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_yrun, 2, s3d_000, "real(8)", OPS_RW), &
                    ops_arg_gbl(nspec, 1, "integer", OPS_READ))

!   SET TEMPERATURE PROFILE
!   -----------------------
    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(flamin_kernel_eqD, "A = var1 + B_xdimonly*(var2-var1)", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_dat(d_crin, 1, s3d_000_strid3d_x, "real(8)", OPS_READ), &
                    ops_arg_gbl(trinr, 1, "real(8)", OPS_READ), &
                    ops_arg_gbl(trinp, 1, "real(8)", OPS_READ))

!   SET DENSITY PROFILE ASSUMING CONSTANT PRESSURE
!   -------------------
!   PRESSURE SET IN CONTROL FILE

    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(set_zero_kernel, "set zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(8)", OPS_WRITE))

    DO ispec = 1,nspec
        rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
        call ops_par_loop(flamin_kernel_eqE, "A = A + var*B_multidim", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(8)", OPS_INC), &
                        ops_arg_dat(d_yrun, 2, s3d_000, "real(8)", OPS_READ), &
                        ops_arg_gbl(rgspec(ispec), 1, "real(8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

    END DO

    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(flamin_kernel_eqF, "A = var/(B*C)", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_gbl(prin, 1, "real(8)", OPS_READ))

!   SET VELOCITY PROFILE ASSUMING CONSTANT MASS FLUX
!   --------------------
!   INITIAL (INLET) VEOCITY SET IN CONTROL FILE
    flxmas = drin*urin
    rangexyz = (/1,nxglbl,1,nyglbl,1,nzglbl/)
    call ops_par_loop(flamin_kernel_eqG, "A = var/B", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(8)", OPS_WRITE), &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(8)", OPS_READ), &
                    ops_arg_gbl(flxmas, 1, "real(8)", OPS_READ))

!   =========================================================================

9000  FORMAT(a)

END SUBROUTINE flamin
