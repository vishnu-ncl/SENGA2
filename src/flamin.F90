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
    real(kind=8) :: trinr,trinp
    real(kind=8) :: yrinr(nspcmx),yrinp(nspcmx)
    real(kind=8) :: deltag
    real(kind=8) :: flxmas
    integer(kind=4) :: ic,jx,kc,ispec
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   SPECIFY INITIAL THERMOCHEMICAL FIELD HERE
!   =========================================

!   SET PRODUCT TEMPERATURE
!   -----------------------
!   REACTANT TEMPERATURE SET IN CONTROL FILE
    trinr = trin
    trinp = 2200.0d0

!   SET SPECIES MASS FRACTION
!   -------------------------
!   OVERRIDE MASS FRACTION VALUES SET IN CONTROL FILE

!   REACTANT
    yrinr(1) = 1.99886d-2
    yrinr(2) = 2.286239d-1
    do ispec = 3,nspm1
      yrinr(ispec) = zero
    enddo

    yrinr(nspec) = zero
    do ispec = 1,nspm1
      yrinr(nspec) = yrinr(nspec) + yrinr(ispec)
    enddo
    yrinr(nspec) = one - yrinr(nspec)

!   PRODUCTS
    yrinp(1) = zero
    yrinp(2) = 6.85323d-2
    yrinp(3) = 1.798974d-1
    do ispec = 4,nspm1
      yrinp(ispec) = zero
    enddo
    yrinp(nspec) = zero
    do ispec = 1,nspm1
      yrinp(nspec) = yrinp(nspec) + yrinp(ispec)
    enddo
    yrinp(nspec) = one - yrinp(nspec)

!   WRITE TO REPORT FILE
    if(iproc.eq.0)then

!      OPEN(UNIT=NCREPT,FILE=FNREPT,STATUS='OLD',FORM='FORMATTED')

!!     GO TO EOF
!1000    CONTINUE
!        READ(NCREPT,9000,END=1010)
!        GOTO 1000
!1010    BACKSPACE(NCREPT)

      write(ncrept,*)
      write(ncrept,*)'flamin: reactant mass fractions:'
      do ispec = 1,nspec
        write(ncrept,'(I5,1PE15.7)')ispec,yrinr(ispec)
      enddo
      write(ncrept,*)

      write(ncrept,*)'flamin: product mass fractions:'
      do ispec = 1,nspec
        write(ncrept,'(I5,1PE15.7)')ispec,yrinp(ispec)
      enddo
      write(ncrept,*)

      write(ncrept,*)'flamin: reactant and product temperatures:'
      write(ncrept,'(2(1PE15.7))')trinr,trinp
      write(ncrept,*)

!      CLOSE(NCREPT)

    endif

!   GLOBAL INDEXING
!   ---------------
    deltag = xgdlen/(real(nxglbl-1))

!   SET REACTION PROGRESS VARIABLE PROFILE
!   --------------------------------------
!   SIMPLE 1D LEFT-FACING ERROR FUNCTION PROFILE
    rangexyz = [1,nxglbl,1,1,1,1]
    call ops_par_loop(flamin_kernel_defc, "Calculating reaction progress variable", senga_grid, 3, rangexyz, &
            ops_arg_dat(d_crin,   1, s3d_000_strid3d_x, "real(kind=8)", OPS_WRITE), &
            ops_arg_gbl(deltag, 1, "real(kind=8)", OPS_READ), &
            ops_arg_gbl(clocat, 1, "real(kind=8)", OPS_READ), &
            ops_arg_gbl(cthick, 1, "real(kind=8)", OPS_READ), &
            ops_arg_idx())

!   SET SPECIES MASS FRACTION PROFILES
!   ----------------------------------
    do ispec = 1, nspm1
      rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]

      call ops_par_loop(flamin_kernel_eqA, "calculate streamwise profile", senga_grid, 3, rangexyz, &
              ops_arg_dat(d_yrun(ispec),       1, s3d_000, "real(kind=8)", OPS_WRITE), &
              ops_arg_dat(d_crin,       1, s3d_000_strid3d_x, "real(kind=8)", OPS_READ), &
              ops_arg_gbl(yrinr(ispec), 1, "real(kind=8)", OPS_READ), &
              ops_arg_gbl(yrinp(ispec), 1, "real(kind=8)", OPS_READ))

    enddo

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
            ops_arg_dat(d_yrun(nspec), 1, s3d_000, "real(kind=8)", OPS_WRITE)) 

    do ispec=1,nspm1
      call ops_par_loop(flamin_kernel_eqB, "Summing mass fractions", senga_grid, 3, rangexyz, &
              ops_arg_dat(d_yrun(nspec), 1, s3d_000, "real(kind=8)", OPS_WRITE), &
              ops_arg_dat(d_yrun(ispec), 1, s3d_000, "real(kind=8)", OPS_READ))
    end do

    call ops_par_loop(flamin_kernel_eqC, "A = 1-A", senga_grid, 3, rangexyz, &
            ops_arg_dat(d_yrun(nspec), 1, s3d_000, "real(kind=8)", OPS_RW))

    call ops_par_loop(flamin_kernel_eqA, "calculate streamwise profile", senga_grid, 3, rangexyz, &
            ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
            ops_arg_dat(d_crin, 1, s3d_000_strid3d_x, "real(kind=8)", OPS_READ), &
            ops_arg_gbl(trinr , 1, "real(kind=8)", OPS_READ), &
            ops_arg_gbl(trinp , 1, "real(kind=8)", OPS_READ))

!   SET DENSITY PROFILE ASSUMING CONSTANT PRESSURE
!   -------------------
!   PRESSURE SET IN CONTROL FILE

    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz, &
            ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_WRITE))

    do ispec=1,nspec
      call ops_par_loop(flamin_kernel_eqD, "calculate R", senga_grid, 3, rangexyz, &
              ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_RW), &
              ops_arg_dat(d_yrun(ispec),1, s3d_000, "real(kind=8)", OPS_READ), &
              ops_arg_gbl(rgspec(ispec), 1, "real(kind=8)",OPS_READ))
    end do

    call ops_par_loop(flamin_kernel_eqE, "Calculate density", senga_grid, 3, rangexyz, &
            ops_arg_dat(d_drun,   1, s3d_000, "real(kind=8)", OPS_WRITE), &
            ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ), &
            ops_arg_dat(d_trun,   1, s3d_000, "real(kind=8)", OPS_READ), &
            ops_arg_gbl(prin,     1, "real(kind=8)", OPS_READ))

!   SET VELOCITY PROFILE ASSUMING CONSTANT MASS FLUX
!   --------------------
!   INITIAL (INLET) VEOCITY SET IN CONTROL FILE
    flxmas = drin*urin
    call ops_par_loop(flamin_kernel_eqF, "Calculate velocity", senga_grid, 3, rangexyz, &
            ops_arg_dat(d_urun, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
            ops_arg_dat(d_drun, 1, s3d_000, "real(kind=8)", OPS_READ), &
            ops_arg_gbl(flxmas, 1, "real(kind=8)", OPS_READ))

!   =========================================================================

9000  FORMAT(a)

END SUBROUTINE flamin
