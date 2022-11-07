SUBROUTINE radcal

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
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
    real(kind=dp) :: plspec,fornow
    integer :: ic,jc,kc,ispec,jspec,icp
    integer :: rangexyz(6)

!   BEGIN
!   =====

!   =========================================================================

!   BUILD THE PLANCK MEAN ABSORPTION COEFFICIENT OF THE MIXTURE
!   -----------------------------------------------------------

!   INITIALISE THE ACCUMULATOR
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(set_zero_kernel, "set zero", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_WRITE))

!   -------------------------------------------------------------------------

!   RUN THROUGH ALL RADIATING SPECIES
    DO jspec = 1, nsprad
  
!       PLANCK MEAN ABSORPTION COEFFICIENT OF EACH SPECIES
        DO kc = kstal,kstol
        DO jc = jstal,jstol
            DO ic = istal,istol
        
                fornow = trun(ic,jc,kc)
                plspec = akprad(nkprad(jspec),jspec)
                DO icp = nkprm1(jspec),1,-1
                    plspec = plspec*fornow + akprad(icp,jspec)
                END DO
                store2(ic,jc,kc) = plspec
        
            END DO
        END DO
        END DO
  
!       SPECIES ID
        ispec = nsprid(jspec)
  
!       ADD THE SPECIES CONTRIBUTION
        rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
        call ops_par_loop(radcal_kernel_addspecies, "ADD THE SPECIES CONTRIBUTION", senga_grid, 3, rangexyz,  &
                        ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_WRITE), &
                        ops_arg_dat(d_yrhs, 9, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_dat(d_store2, 1, s3d_000, "real(dp)", OPS_READ), &
                        ops_arg_gbl(rgspec(ispec), 1, "real(dp)", OPS_READ), &
                        ops_arg_gbl(ispec, 1, "integer", OPS_READ))

    END DO

!   =========================================================================

!   INCLUDE THE RADIATION TERM IN THE ENERGY EQUATION
    rangexyz = (/istal,istol,jstal,jstol,kstal,kstol/)
    call ops_par_loop(radcal_kernel_addradiation, "INCLUDE THE RADIATION TERM IN THE ENERGY EQUATION", senga_grid, 3, rangexyz,  &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_trun, 1, s3d_000, "real(dp)", OPS_READ), &
                    ops_arg_dat(d_store1, 1, s3d_000, "real(dp)", OPS_READ))

!   =========================================================================

END SUBROUTINE radcal
