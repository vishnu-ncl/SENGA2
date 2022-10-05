SUBROUTINE lincom

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:26:06

!     *************************************************************************

!     LINCOM
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     15-JAN-2003:  CREATED
!     08-AUG-2012:  RSC EVALUATE ALL SPECIES

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     COMPUTES INTERMEDIATE SOLUTION VALUES IN ERK SCHEME
!     BY DOING LINEAR COMBINATIONS OF LEFT- AND RIGHT-HAND SIDES

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
real(kind=dp) :: fornow
INTEGER :: ic,jc,kc,ispec
integer :: rangexyz(6)



!     BEGIN
!     =====

!     =========================================================================

!     ERK SUBSTEP
!     ===========

!     -------------------------------------------------------------------------
!     NOTE: ALL ERK ERROR ARRAYS ARE INITIALISED TO ZERO IN SUBROUTINE ADAPTT
!     -------------------------------------------------------------------------

!     DENSITY
!     -------
    rangexyz = (/istald,istold,jstald,jstold,kstald,kstold/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_derr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_drun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_drhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!     -------------------------------------------------------------------------

!     U-VELOCITY
!     ----------
    rangexyz = (/istalu,istolu,jstalu,jstolu,kstalu,kstolu/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_uerr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_urun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_urhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!     -------------------------------------------------------------------------

!     V-VELOCITY
!     ----------
    rangexyz = (/istalv,istolv,jstalv,jstolv,kstalv,kstolv/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_verr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_vrun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_vrhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!     -------------------------------------------------------------------------

!     W-VELOCITY
!     ----------
    rangexyz = (/istalw,istolw,jstalw,jstolw,kstalw,kstolw/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_werr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_wrun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_wrhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!     -------------------------------------------------------------------------

!     STAGNATION INTERNAL ENERGY
!     --------------------------
    rangexyz = (/istale,istole,jstale,jstole,kstale,kstole/)
    call ops_par_loop(lincom_kernel_main, "lincom_main", senga_grid, 3, rangexyz, &
                    ops_arg_dat(d_eerr, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_erun, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_dat(d_erhs, 1, s3d_000, "real(dp)", OPS_WRITE), &
                    ops_arg_gbl(rkerr(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rklhs(irkstp), 1, "real(dp)", OPS_READ), &
                    ops_arg_gbl(rkrhs(irkstp), 1, "real(dp)", OPS_READ))

!     -------------------------------------------------------------------------

!     SPECIES MASS FRACTIONS
!     ----------------------
!     RSC 08-AUG-2012 EVALUATE ALL SPECIES
!      DO ISPEC = 1,NSPM1
DO ispec = 1,nspec
  
  DO kc = kstaly,kstoly
    DO jc = jstaly,jstoly
      DO ic = istaly,istoly
        
        yerr(ic,jc,kc,ispec) = yerr(ic,jc,kc,ispec)  &
            + rkerr(irkstp)*yrhs(ic,jc,kc,ispec)
        
        fornow = yrun(ic,jc,kc,ispec)
        yrun(ic,jc,kc,ispec) = fornow + rklhs(irkstp)*yrhs(ic,jc,kc,ispec)
        yrhs(ic,jc,kc,ispec) = fornow + rkrhs(irkstp)*yrhs(ic,jc,kc,ispec)
        
      END DO
    END DO
  END DO
  
END DO

!     -------------------------------------------------------------------------

!C     NTH SPECIES
!      DO KC = KSTALY,KSTOLY
!        DO JC = JSTALY,JSTOLY
!          DO IC = ISTALY,ISTOLY

!            YRUN(IC,JC,KC,NSPEC) = ZERO
!            YRHS(IC,JC,KC,NSPEC) = ZERO

!          ENDDO
!        ENDDO
!      ENDDO

!      DO ISPEC = 1,NSPM1
!        DO KC = KSTALY,KSTOLY
!          DO JC = JSTALY,JSTOLY
!            DO IC = ISTALY,ISTOLY

!              YRUN(IC,JC,KC,NSPEC) = YRUN(IC,JC,KC,NSPEC)
!     +                             + YRUN(IC,JC,KC,ISPEC)
!              YRHS(IC,JC,KC,NSPEC) = YRHS(IC,JC,KC,NSPEC)
!     +                             + YRHS(IC,JC,KC,ISPEC)

!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO

!      DO KC = KSTALY,KSTOLY
!        DO JC = JSTALY,JSTOLY
!          DO IC = ISTALY,ISTOLY

!            YRUN(IC,JC,KC,NSPEC)
!     +        = DRUN(IC,JC,KC)*(ONE-YRUN(IC,JC,KC,NSPEC)/DRUN(IC,JC,KC))

!            YRHS(IC,JC,KC,NSPEC)
!     +        = DRHS(IC,JC,KC)*(ONE-YRHS(IC,JC,KC,NSPEC)/DRHS(IC,JC,KC))

!          ENDDO
!        ENDDO
!      ENDDO

!     =========================================================================


END SUBROUTINE lincom
