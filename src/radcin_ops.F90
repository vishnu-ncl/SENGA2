
! Auto-generated at 2026-04-28 18:43:09.519739 by ops-translator
SUBROUTINE radcin

  USE ops_fortran_declarations
  USE ops_fortran_rt_support

  USE OPS_CONSTANTS
  USE, INTRINSIC :: ISO_C_BINDING

  USE com_senga
  USE com_ops_senga

  !   *************************************************************************

  !   RADCIN
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
  !   USING OPTICALLY-THIN ASSUMPTION: Ju et al: JFM 342, 315-334, 1997.
  !   INITIALISES DATA

  !   *************************************************************************

  !   GLOBAL DATA
  !   ===========
  !   -------------------------------------------------------------------------
  !   -------------------------------------------------------------------------

  !   LOCAL DATA
  !   ==========
  REAL(KIND = 8) :: pcount
  INTEGER(KIND = 4) :: ispec, jspec, icoeff
  INTEGER(KIND = 4) :: nrdspc, ncount
  INTEGER(KIND = 4) :: iroot
  CHARACTER(LEN = 10) :: spcrad(nspcmx)


  !   BEGIN
  !   =====

  !   =========================================================================

  !   CONTROL FLAGS
  !   DEFAULTS: RADIATION TREATMENT NOT ACTIVATED
  flradn = .FALSE.

    !   =========================================================================

    !   CHECK SWITCH FOR RADIATION TREATMENT
    IF (nfradn == 1) THEN

    !       =======================================================================

    !       DATA IS READ BY THE LOWEST-RANKED PROCESSOR AND BROADCAST TO THE OTHERS
    !       -----------------------------------------------------------------------

    !       SET ID OF ROOT (BROADCASTING) PROCESSOR
    iroot = 0

    !       INITIALISE THE BROADCAST COUNTER
    pcount = zero

      !       INITIALISE THE BROADCAST ARRAY
      DO ncount = 1, nparay
      parray(ncount) = zero
    END DO

      !       =======================================================================

      !       LOWEST-RANKED PROCESSOR
      IF (iproc == 0) THEN

      !       =======================================================================

      !           READ THE RADIATION DATA FILE
      !           ----------------------------
      !           ---------------------------------------------------------------------

      OPEN(UNIT = ncradn, FILE = fnradn, STATUS = 'OLD', FORM = 'FORMATTED')

      !           ---------------------------------------------------------------------

      !           READ AND IGNORE THE HEADER
      READ(ncradn, *)
      READ(ncradn, *)
      READ(ncradn, *)
      READ(ncradn, *)
      READ(ncradn, *)

      !           ---------------------------------------------------------------------

      !           READ SPECIES LIST AND CHECK AGAINST CHEMICAL DATA
      READ(ncradn, *)
      READ(ncradn, '(I5)') nrdspc
      IF (nrdspc /= nspec) THEN
        WRITE(ncradn, *)
        WRITE(ncradn, *) 'Warning: RADCIN: mismatch in no. of species'
        WRITE(ncradn, *) 'CHEMIN: ', nspec, '; RADCIN: ', nrdspc
        WRITE(ncradn, *)
      END IF

        DO jspec = 1, nspec
        READ(ncradn, '(I5,3X,A)') ispec, spcrad(ispec)
        IF (spcsym(ispec) /= spcrad(ispec)) THEN
          WRITE(ncrept, *)
          WRITE(ncrept, *) 'Warning: RADCIN: mismatch in species symbol'
          WRITE(ncrept, *) 'CHEMIN: ', spcsym(ispec), '; RADCIN: ', spcrad(ispec)
          WRITE(ncrept, *)
        END IF
      END DO

      !           RADIATION DATA
      READ(ncradn, *)
      READ(ncradn, *) nsprad
      DO ispec = 1, nsprad
        READ(ncradn, *) jspec, nkprad(ispec)
        READ(ncradn, *) (akprad(icoeff, ispec), icoeff = 1, nkprad(ispec))
        nsprid(ispec) = jspec
      END DO

      !           ---------------------------------------------------------------------

      !           READ END-OF-FILE LINE
      READ(ncradn, *)
      CLOSE(UNIT = ncradn)

      !           =====================================================================

      !           COMPRESS THE DATA
      !           -----------------

      !           RADIATION DATA
      ncount = 1
      DO ispec = 1, nspec
        ncount = ncount + 1
        parray(ncount) = REAL(nsprad, kind = 8)
        DO icoeff = 1, nsprad
          ncount = ncount + 1
          parray(ncount) = akprad(icoeff, ispec)
        END DO
        ncount = ncount + 1
        parray(ncount) = REAL(nsprid(ispec), kind = 8)
      END DO

      !           ---------------------------------------------------------------------

      !           BROADCAST COUNTER
      pcount = REAL(ncount, kind = 8)

        !           CHECK BROADCAST COUNTER AGAINST BROADCAST ARRAY SIZE
        IF (ncount > nparay) THEN
        WRITE(ncrept, *)
        WRITE(ncrept, *) 'Warning: RADCIN: broadcast array size too small'
        WRITE(ncrept, *) 'Actual: ', ncount, '; Available: ', nparay
        WRITE(ncrept, *)
      END IF

      !           =====================================================================

    END IF

    !       =======================================================================

    !       BROADCAST (OR RECEIVE) THE COUNTER
    !       ----------------------------------
    CALL p_bcst(pcount, 1, 1, iroot)
    ncount = NINT(pcount)

    !       BROADCAST (OR RECEIVE) THE DATA
    !       -------------------------------
    CALL p_bcst(parray, nparay, ncount, iroot)

      !       =======================================================================

      !       NOT THE LOWEST-RANKED PROCESSOR
      IF (iproc /= 0) THEN

      !           =====================================================================

      !           UNCOMPRESS THE DATA
      !           -------------------

      !           RADIATION DATA
      ncount = 1
      DO ispec = 1, nspec
        ncount = ncount + 1
        nsprad = NINT(parray(ncount))
        DO icoeff = 1, nsprad
          ncount = ncount + 1
          akprad(icoeff, ispec) = parray(ncount)
        END DO
        ncount = ncount + 1
        nsprid(ispec) = NINT(parray(ncount))
      END DO

      !           =====================================================================

    END IF

      !       =======================================================================

      !       EVALUATE DERIVED QUANTITIES
      !       ===========================

      !       COEFFICIENT COUNTERS
      DO ispec = 1, nspec
      nkprm1(ispec) = nkprad(ispec) - 1
    END DO

    !       =======================================================================

    !       CONSTANTS
    foursb = four * stefbo
    trfrth = trefrn * trefrn * trefrn * trefrn
#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("foursb", 1, "real(kind=8)", foursb)
#endif

#ifdef OPS_F2C_INTEROP
CALL ops_decl_const("trfrth", 1, "real(kind=8)", trfrth)
#endif
#ifdef OPS_WITH_OMPOFFLOADFOR
!$OMP TARGET UPDATE TO(foursb)
!$OMP TARGET UPDATE TO(trfrth)
#endif

!       =======================================================================

    !       CONTROL FLAGS
    flradn = .TRUE.

    !       =======================================================================

  END IF

  !   =========================================================================

END SUBROUTINE radcin