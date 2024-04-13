SUBROUTINE diffin

!   *************************************************************************

!   DIFFIN
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   06-JAN-2013:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   INITIALISES DATA FOR MULTI-SPECIES MIXTURE AVERAGED MOLECULAR TRANSPORT

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
    use com_senga
!   -------------------------------------------------------------------------

!   PARAMETER
!   =========
    real(kind=8), parameter :: wmltol = 0.001_8

!   LOCAL DATA
!   ==========
    real(kind=8) :: pcount,fornow
    integer(kind=4) :: ispec,jspec,kspec,icoeff
    integer(kind=4) :: ncount,ndspec
    integer(kind=4) :: iroot
    CHARACTER (LEN=10) :: spcdif(nspcmx)

!   BEGIN
!   =====

!   =========================================================================

!   CONTROL FLAGS
!   DEFAULTS: MIXTURE AVERAGED TRANSPORT NOT ACTIVATED
    flmavt = .false.
    flmixw = .false.
    flmixp = .false.
    flmixt = .false.
    DO ispec = 1, nspec
        flmsor(ispec) = .false.
        flmduf(ispec) = .false.
        flmtdr(ispec) = .false.
    END DO

!   =========================================================================

!   CHECK SWITCH FOR MIXTURE AVERAGED TRANSPORT
    IF(nfmavt == 1) THEN

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
        IF(iproc == 0) THEN

!           =======================================================================

!           READ THE MOLECULAR TRANSPORT DATA FILE
!           --------------------------------------
!           ---------------------------------------------------------------------

            OPEN(UNIT=ncdiff,FILE=fndiff,STATUS='OLD',FORM='FORMATTED')

!           ---------------------------------------------------------------------

!           READ AND IGNORE THE HEADER
            READ(ncdiff,*)
            READ(ncdiff,*)
            READ(ncdiff,*)
            READ(ncdiff,*)
            READ(ncdiff,*)

!           ---------------------------------------------------------------------

!           READ SPECIES LIST AND CHECK AGAINST CHEMICAL DATA
            READ(ncdiff,*)
            READ(ncdiff,'(I5)')ndspec
            IF(ndspec /= nspec) THEN
                WRITE(ncrept,*)
                WRITE(ncrept,*)'Warning: DIFFIN: mismatch in no. of species'
                WRITE(ncrept,*)'CHEMIN: ',nspec,'; DIFFIN: ',ndspec
                WRITE(ncrept,*)
            END IF

            DO jspec = 1, nspec
                READ(ncdiff,'(I5,3X,A)')ispec,spcdif(ispec)
                IF(spcsym(ispec) /= spcdif(ispec)) THEN
                    WRITE(ncrept,*)
                    WRITE(ncrept,*)'Warning: DIFFIN: mismatch in species symbol'
                    WRITE(ncrept,*)'CHEMIN: ',spcsym(ispec), '; DIFFIN: ',spcdif(ispec)
                    WRITE(ncrept,*)
                END IF
            END DO

!           MOLECULAR TRANSPORT DATA
            READ(ncdiff,*)
            READ(ncdiff,*)pdifgb,tdifgb

            READ(ncdiff,*)
            DO ispec = 1,nspec
                READ(ncdiff,*)jspec,ncovis
                IF(jspec /= ispec) THEN
                    WRITE(ncrept,*)
                    WRITE(ncrept,*)'Warning: DIFFIN: species ID mismatch'
                    WRITE(ncrept,*)'Viscosity data'
                    WRITE(ncrept,*)'Expected: ',ispec,'; Input: ',jspec
                END IF
                READ(ncdiff,*) (viscco(icoeff,ispec),icoeff=1,ncovis)
            END DO

            READ(ncdiff,*)
            DO ispec = 1,nspec
                READ(ncdiff,*)jspec,ncocon
                IF(jspec /= ispec) THEN
                    WRITE(ncrept,*)
                    WRITE(ncrept,*)'Warning: DIFFIN: species ID mismatch'
                    WRITE(ncrept,*)'Thermal conductivity data'
                    WRITE(ncrept,*)'Expected: ',ispec,'; Input: ',jspec
                END IF
                READ(ncdiff,*) (condco(icoeff,ispec),icoeff=1,ncocon)
            END DO

            READ(ncdiff,*)
            DO ispec = 1,nspec
                READ(ncdiff,*)jspec,ncodif
                IF(jspec /= ispec) THEN
                    WRITE(ncrept,*)
                    WRITE(ncrept,*)'Warning: DIFFIN: species ID mismatch'
                    WRITE(ncrept,*)'Binary diffusion coefficient data ISPEC'
                    WRITE(ncrept,*)'Expected: ',ispec,'; Input: ',jspec
                END IF
                DO jspec = 1,ispec
                    READ(ncdiff,*)kspec, (diffco(icoeff,jspec,ispec),icoeff=1,ncodif)
                    IF(kspec /= jspec) THEN
                        WRITE(ncrept,*)
                        WRITE(ncrept,*)'Warning: DIFFIN: species ID mismatch'
                        WRITE(ncrept,*)'Binary diffusion coefficient data JSPEC'
                        WRITE(ncrept,*)'Expected: ',jspec,'; Input: ',kspec
                    END IF
                END DO
            END DO

            READ(ncdiff,*)
            DO ispec = 1,nspec
                READ(ncdiff,*)jspec,ncotdr
                IF(jspec /= ispec) THEN
                    WRITE(ncrept,*)
                    WRITE(ncrept,*)'Warning: DIFFIN: species ID mismatch'
                    WRITE(ncrept,*)'Thermal diffusion coefficient data ISPEC'
                    WRITE(ncrept,*)'Expected: ',ispec,'; Input: ',jspec
                END IF
                DO jspec = 1,ispec-1
                    READ(ncdiff,*)kspec, (tdrcco(icoeff,jspec,ispec),icoeff=1,ncotdr)
                    IF(kspec /= jspec) THEN
                        WRITE(ncrept,*)
                        WRITE(ncrept,*)'Warning: DIFFIN: species ID mismatch'
                        WRITE(ncrept,*)'Thermal diffusion coefficient data JSPEC'
                        WRITE(ncrept,*)'Expected: ',ispec,'; Input: ',jspec
                    END IF
                END DO
            END DO

!           ---------------------------------------------------------------------

!           READ END-OF-FILE LINE
            READ(ncdiff,*)
            CLOSE(ncdiff)

!           =====================================================================

!           COMPRESS THE DATA
!           -----------------

!           MOLECULAR TRANSPORT DATA
            ncount = 1
            parray(ncount) = pdifgb
            ncount = ncount + 1
            parray(ncount) = tdifgb
            DO ispec = 1,nspec
                ncount = ncount + 1
                parray(ncount) = REAL(ncovis,kind=8)
                DO icoeff = 1, ncovis
                    ncount = ncount + 1
                    parray(ncount) = viscco(icoeff,ispec)
                END DO
                ncount = ncount + 1
                parray(ncount) = REAL(ncocon,kind=8)
                DO icoeff = 1, ncocon
                    ncount = ncount + 1
                    parray(ncount) = condco(icoeff,ispec)
                END DO
                ncount = ncount + 1
                parray(ncount) = REAL(ncodif,kind=8)
                DO jspec = 1, ispec
                    DO icoeff = 1, ncodif
                        ncount = ncount + 1
                        parray(ncount) = diffco(icoeff,jspec,ispec)
                    END DO
                END DO
                ncount = ncount + 1
                parray(ncount) = REAL(ncotdr,kind=8)
                DO jspec = 1, ispec-1
                    DO icoeff = 1, ncotdr
                        ncount = ncount + 1
                        parray(ncount) = tdrcco(icoeff,jspec,ispec)
                    END DO
                END DO
            END DO

!           ---------------------------------------------------------------------

!           BROADCAST COUNTER
            pcount = REAL(ncount,kind=8)

!           CHECK BROADCAST COUNTER AGAINST BROADCAST ARRAY SIZE
            IF(ncount > nparay) THEN
                WRITE(ncrept,*)
                WRITE(ncrept,*)'Warning: DIFFIN: broadcast array size too small'
                WRITE(ncrept,*)'Actual: ',ncount,'; Available: ',nparay
                WRITE(ncrept,*)
            END IF

!           =====================================================================

        END IF

!       =======================================================================

!       BROADCAST (OR RECEIVE) THE COUNTER
!       ----------------------------------
        call p_bcst(pcount,1,1,iroot)
        ncount = nint(pcount)

!       BROADCAST (OR RECEIVE) THE DATA
!       -------------------------------
        call p_bcst(parray,nparay,ncount,iroot)

!       =======================================================================

!       NOT THE LOWEST-RANKED PROCESSOR
        IF(iproc /= 0) THEN

!           =====================================================================

!           UNCOMPRESS THE DATA
!           -------------------

!           MOLECULAR TRANSPORT DATA
            ncount = 1
            pdifgb = parray(ncount)
            ncount = ncount + 1
            tdifgb = parray(ncount)
            DO ispec = 1,nspec
                ncount = ncount + 1
                ncovis = nint(parray(ncount))
                DO icoeff = 1, ncovis
                    ncount = ncount + 1
                    viscco(icoeff,ispec) = parray(ncount)
                END DO
                ncount = ncount + 1
                ncocon = nint(parray(ncount))
                DO icoeff = 1, ncocon
                    ncount = ncount + 1
                    condco(icoeff,ispec) = parray(ncount)
                END DO
                ncount = ncount + 1
                ncodif = nint(parray(ncount))
                DO jspec = 1, ispec
                    DO icoeff = 1, ncodif
                        ncount = ncount + 1
                        diffco(icoeff,jspec,ispec) = parray(ncount)
                    END DO
                END DO
                ncount = ncount + 1
                ncotdr = nint(parray(ncount))
                DO jspec = 1, ispec-1
                    DO icoeff = 1, ncotdr
                        ncount = ncount + 1
                        tdrcco(icoeff,jspec,ispec) = parray(ncount)
                    END DO
                END DO
            END DO

!           =====================================================================

        END IF

!       =======================================================================

!       EVALUATE DERIVED QUANTITIES
!       ===========================

!       COEFFICIENT COUNTERS
        ncovm1 = ncovis - 1
        ncocm1 = ncocon - 1
        ncodm1 = ncodif - 1
        ncotm1 = ncotdr - 1

!       FILL OUT THE DIFFUSION COEFFICIENT MATRIX
!       MATRIX IS SYMMETRIC
        DO ispec = 1, nspec
            DO jspec = ispec+1, nspec
                DO icoeff = 1, ncodif
                    fornow = diffco(icoeff,ispec,jspec)
                    diffco(icoeff,jspec,ispec) = fornow
                END DO
            END DO
        END DO

!       FILL OUT THE THERMAL DIFFUSION COEFFICIENT MATRIX
!       RSC 26-MAY-2023
        DO ispec = 1, nspec
!           LEADING DIAGONAL IS ZERO
            DO icoeff = 1, ncotdr
                tdrcco(icoeff,ispec,ispec) = zero
            END DO

!           MATRIX IS ANTISYMMETRIC
!           ELEMENTS OF THE MATRIX TRANSPOSE ARE STORED
            DO jspec = 1, ispec-1
                DO icoeff = 1, ncotdr
                    fornow = tdrcco(icoeff,jspec,ispec)
                    tdrcco(icoeff,jspec,ispec) = -fornow
                END DO
            END DO

            DO jspec = ispec+1, nspec
                DO icoeff = 1, ncotdr
                    fornow = tdrcco(icoeff,ispec,jspec)
                    tdrcco(icoeff,jspec,ispec) = fornow
                END DO
            END DO
        END DO

!       PRECOMPUTE ELEMENTS OF WILKES COMBINATION RULE FOR VISCOSITY
!       NOTE THAT THE WILKES MATRIX IS NOT SYMMETRIC
!       ELEMENTS OF THE MATRIX TRANSPOSE ARE STORED
        DO ispec = 1, nspec
            DO jspec = 1, nspec
                fornow = wmolar(ispec)/wmolar(jspec)
                wilko2(jspec,ispec) = one/SQRT(SQRT(fornow))
                fornow = eight*(one + fornow)
                wilko1(jspec,ispec) = one/SQRT(fornow)
            END DO
        END DO

!       =======================================================================

!       CONTROL FLAGS
        flmavt = .true.
        IF(nfmixw == 1) flmixw = .true.
        IF(nfmixp == 1) flmixp = .true.

        IF(nfmsor == 1) THEN
            flmixt = .true.
            DO ispec = 1, nspec
                IF(wmolar(ispec) <= (wmltdr+wmltol)) THEN
                    flmsor(ispec) = .true.
                    flmtdr(ispec) = .true.
                END IF
            END DO
        END IF

        IF(nfmduf == 1) THEN
            flmixt = .true.
            DO ispec = 1, nspec
                IF(wmolar(ispec) <= (wmltdr+wmltol)) THEN
                    flmduf(ispec) = .true.
                    flmtdr(ispec) = .true.
                END IF
            END DO
        END IF

!       =======================================================================

    END IF

!   =========================================================================

END SUBROUTINE diffin
