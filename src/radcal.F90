SUBROUTINE radcal

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   RADCAL
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   14-JUL-2013:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   RADIATION TREATMENT
!   USING OPTICALLY THIN ASSUMPTION: Ju et al: JFM 342, 315-334, 1997.
!   AFTER TOM DUNSTAN 2012

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    real(kind=8) :: plspec,fornow
    integer(kind=4) :: ispec,jspec,icp
    integer(kind=4) :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   BUILD THE PLANCK MEAN ABSORPTION COEFFICIENT OF THE MIXTURE
!   -----------------------------------------------------------

!   INITIALISE THE ACCUMULATOR
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(set_zero_kernel, "set_zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_WRITE))

!   -------------------------------------------------------------------------

!   RUN THROUGH ALL RADIATING SPECIES
    DO jspec = 1, nsprad

!       PLANCK MEAN ABSORPTION COEFFICIENT OF EACH SPECIES
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(radcal_kernel_meancoef, "PLANCK MEAN ABSORPTION COEF", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_WRITE), &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(akprad, ncfrmx*nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(nkprad, nspcmx, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(nkprm1, nspcmx, "integer(kind=4)", OPS_READ), &
                        ops_arg_gbl(jspec, 1, "integer(kind=4)", OPS_READ))

!       SPECIES ID
        ispec = nsprid(jspec)

!       ADD THE SPECIES CONTRIBUTION
        rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
        call ops_par_loop(radcal_kernel_addspecies, "ADD THE SPECIES CONTRIBUTION", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_INC), &
                        ops_arg_dat(d_yrhs(ispec), 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(rgspec, nspcmx, "real(kind=8)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer(kind=4)", OPS_READ))

    END DO

!   =========================================================================

!   INCLUDE THE RADIATION TERM IN THE ENERGY EQUATION
    rangexyz = [1,nxglbl,1,nyglbl,1,nzglbl]
    call ops_par_loop(radcal_kernel_addradiation, "INCLUDE THE RADIATION TERM IN THE ENERGY EQUATION", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(kind=8)", OPS_INC), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(kind=8)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(kind=8)", OPS_READ))

!   =========================================================================

END SUBROUTINE radcal
