SUBROUTINE chemin

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use com_senga
    use com_ops_senga

!   *************************************************************************

!   CHEMIN
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   08-JAN-2003:  CREATED
!   29-DEC-2006:  RSC BUG FIX INDEXING (THANKS TO MDLLP)

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   INITIALISES DATA FOR MULTI-SPECIES MULTI-STEP CHEMISTRY

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ==========
    real(kind=8) :: pcount
    real(kind=8) :: giblet
    integer(kind=4) :: ispec,istep,ibody,ilind,itroe,isrif
    integer(kind=4) :: jspec,jstep,jbody,jlind,jtroe,jsrif
    integer(kind=4) :: itint,icp
    integer(kind=4) :: ncount
    integer(kind=4) :: iroot


!   BEGIN
!   =====

!   =========================================================================

!   GLOBAL DATA
!   -----------
!   NSPEC total no of species
!   NSTEP total no of steps
!   NBODY total no of third bodies
!   NGIBB total no of Gibbs steps
!   NLIND total no of Lindemann steps
!   NTROE total no of Troe steps
!   NSRIF total no of SRI steps
!   NSPCMX max no of species
!   NSTPMX max no of steps
!   NSSMAX max length of step species-list
!   NRSMAX max length of step reactant- and product-lists
!   NBDYMX max no of third bodies
!   NLLMAX max no of Lindemann steps
!   RPARAM(1,NSTPMX) pre-exponential factor ln(A)
!   RPARAM(2,NSTPMX) temperature exponent n
!   RPARAM(3,NSTPMX) activation energy E/R0
!   NSSPEC(NSSMAX,NSTPMX) is the step species-list
!   NRSPEC(NRSMAX,NSTPMX) is the step reactant-list
!   NPSPEC(NRSMAX,NSTPMX) is the step product-list
!   NRCPEC(NRSMAX,NSTPMX) is the step coefficient reactant-list
!   NPCPEC(NRSMAX,NSTPMX) is the step coefficient product-list
!   NSSLEN(NSTPMX) contains the length of the species-list for each step
!   NRSLEN(NSTPMX) contains the length of the reactant-list for each step
!   NPSLEN(NSTPMX) contains the length of the product-list for each step
!   NRCLEN(NSTPMX) contains the length of the reactant coefficient-list for
!                  each step
!   NPCLEN(NSTPMX) contains the length of the reactant coefficient-list for
!                  each step
!   WMOLAR(NSPCMX) molar mass of species
!   OVWMOL(NSPCMX) reciprocal of molar mass of species
!   DIFFMU(NSSMAX,NSTPMX) is the species delta-mu for each step
!   DIFFMW(NSSMAX,NSTPMX) is the species delta-mu (times molar mass)
!                         for each step
!   CRSPEC(NSSMAX,NSTPMX) is the reactant stoichiometric coefficient-list
!                         for each step
!   CPSPEC(NSSMAX,NSTPMX) is the product stoichiometric coefficient-list
!                         for each step
!   MBLIST(NSTPMX) contains 0 if no third body in a step
!                  else contains the index number of the third body
!   MGLIST(NSTPMX) contains 0 if no Gibbs function evaluation in a step
!                  else contains the index number of the Gibbs step
!   MLLIST(NSTPMX) contains 0 if no Lindemann rate in a step
!                  else contains the index number of the Lindemann step
!   MTLIST(NSTPMX) contains 0 if no Troe form rate in a step
!                  else contains the index number of the Troe form step
!   MSLIST(NSTPMX) contains 0 if no SRI form rate in a step
!                  else contains the index number of the SRI form step
!   EFFY3B(NSPCMX,NBDYMX) contains the third-body efficiencies
!   RCLIND(1,MLLIST(NSTPMX)) Lindemann rate pre-exponential factor ln(A)
!   RCLIND(2,MLLIST(NSTPMX)) Lindemann rate temperature exponent n
!   RCLIND(3,MLLIST(NSTPMX)) Lindemann rate activation energy E/R0
!   RCLIND(4,MLLIST(NSTPMX)) Lindemann rate broadening factor
!   RCTROE(1,MTLIST(NSTPMX)) Troe form rate pre-exponential factor ln(A)
!   RCTROE(2,MTLIST(NSTPMX)) Troe form rate temperature exponent n
!   RCTROE(3,MTLIST(NSTPMX)) Troe form rate activation energy E/R0
!   RCTROE(4,MTLIST(NSTPMX)) Troe form rate temperature parameter alpha
!   RCTROE(5,MTLIST(NSTPMX)) Troe form rate temperature parameter no 1*
!   RCTROE(6,MTLIST(NSTPMX)) Troe form rate temperature parameter no 2*
!   RCTROE(7,MTLIST(NSTPMX)) Troe form rate temperature parameter no 3*
!   RCTROE(8,MTLIST(NSTPMX)) Troe form rate c-const no 1
!   RCTROE(9,MTLIST(NSTPMX)) Troe form rate c-const no 2
!   RCTROE(10,MTLIST(NSTPMX)) Troe form rate n-const no 1
!   RCTROE(11,MTLIST(NSTPMX)) Troe form rate n-const no 2
!   RCTROE(12,MTLIST(NSTPMX)) Troe form rate d-const
!   RCSRIF(1,MSLIST(NSTPMX)) SRI form rate pre-exponential factor ln(A)
!   RCSRIF(2,MSLIST(NSTPMX)) SRI form rate temperature exponent n
!   RCSRIF(3,MSLIST(NSTPMX)) SRI form rate activation energy E/R0
!   RCSRIF(4,MSLIST(NSTPMX)) SRI form rate parameter a
!   RCSRIF(5,MSLIST(NSTPMX)) SRI form rate parameter b
!   RCSRIF(6,MSLIST(NSTPMX)) SRI form rate parameter c
!   RCSRIF(7,MSLIST(NSTPMX)) SRI form rate parameter d
!   RCSRIF(8,MSLIST(NSTPMX)) SRI form rate parameter e
!   =========================================================================

!   DATA IS READ BY THE LOWEST-RANKED PROCESSOR AND BROADCAST TO THE OTHERS
!   -----------------------------------------------------------------------

!   SET ID OF ROOT (BROADCASTING) PROCESSOR
    iroot = 0

!   INITIALISE THE BROADCAST COUNTER
    pcount = zero

!   INITIALISE THE BROADCAST ARRAY
    DO ncount = 1, nparay
        parray(ncount) = zero
    END DO

!   =========================================================================

!   LOWEST-RANKED PROCESSOR
    IF(iproc == 0)THEN

!       =======================================================================

!       READ THE CHEMICAL DATA FILE
!       ---------------------------
!       -----------------------------------------------------------------------

        OPEN(UNIT=ncchem,FILE=fnchem,STATUS='OLD',FORM='FORMATTED')

!           -----------------------------------------------------------------------

!           READ AND IGNORE THE HEADER
            READ(ncchem,*)
            READ(ncchem,*)
            READ(ncchem,*)
            READ(ncchem,*)
            READ(ncchem,*)

!           -----------------------------------------------------------------------

!           SPECIES LIST
            READ(ncchem,*)
            READ(ncchem,'(I5)')nspec
            DO jspec = 1, nspec
                READ(ncchem,'(I5,3X,A)')ispec,spcsym(ispec)
            END DO

!           SPECIES DATA
            READ(ncchem,*)
            READ(ncchem,*)prefgb
            DO jspec = 1,nspec
                READ(ncchem,'(I5,1PE12.4)')ispec,wmolar(ispec)
                READ(ncchem,'(1PE12.4)')clewis(ispec)
                READ(ncchem,'(I5)')ntint(ispec)
                DO itint = 1,ntint(ispec)
                    READ(ncchem,'(2(1PE12.4),I5)')tintlo(itint,ispec),  &
                                        tinthi(itint,ispec),ncofcp(itint,ispec)
                    DO icp = 1,ncofcp(itint,ispec)
                        READ(ncchem,'(1PE15.7)')amolcp(icp,itint,ispec)
                    END DO
                END DO
            END DO

!           -----------------------------------------------------------------------

!           STEP RATE DATA
            READ(ncchem,*)
            READ(ncchem,'(I5)')nstep
            DO jstep = 1, nstep
                READ(ncchem,'(I5,3(1PE12.4))')istep,(rparam(icp,istep),icp=1,3)
            END DO

!           -----------------------------------------------------------------------

!           STEP SPECIES-LIST
            READ(ncchem,*)
            DO jstep = 1,nstep
                READ(ncchem,'(2I5)')istep,nsslen(istep)
                DO jspec = 1, nsslen(istep)
                    READ(ncchem,'(2I5,I8)')istep,ispec,nsspec(ispec,istep)
                END DO
            END DO

!           STEP REACTANT-LIST
            READ(ncchem,*)
            DO jstep = 1,nstep
                READ(ncchem,'(2I5)')istep,nrslen(istep)
                DO jspec = 1, nrslen(istep)
                    READ(ncchem,'(2I5,I8)')istep,ispec,nrspec(ispec,istep)
                END DO
            END DO

!           STEP PRODUCT-LIST
            READ(ncchem,*)
            DO jstep = 1,nstep
                READ(ncchem,'(2I5)')istep,npslen(istep)
                DO jspec = 1, npslen(istep)
                    READ(ncchem,'(2I5,I8)')istep,ispec,npspec(ispec,istep)
                END DO
            END DO

!           STEP REACTANT COEFFICIENT-LIST
            READ(ncchem,*)
            DO jstep = 1,nstep
                READ(ncchem,'(2I5)')istep,nrclen(istep)
                DO jspec = 1, nrclen(istep)
                    READ(ncchem,'(2I5,I8,1PE12.4)')istep,ispec,  &
                                nrcpec(ispec,istep),crspec(ispec,istep)
                END DO
            END DO

!           STEP PRODUCT COEFFICIENT-LIST
            READ(ncchem,*)
            DO jstep = 1,nstep
                READ(ncchem,'(2I5)')istep,npclen(istep)
                DO jspec = 1, npclen(istep)
                    READ(ncchem,'(2I5,I8,1PE12.4)')istep,ispec,  &
                                npcpec(ispec,istep),cpspec(ispec,istep)
                END DO
            END DO

!           SPECIES DELTA-LIST
            READ(ncchem,*)
            DO jstep = 1,nstep
                DO jspec = 1,nsslen(jstep)
                    READ(ncchem,'(2I5,F5.1)')istep,ispec,diffmu(ispec,istep)
                END DO
            END DO

!           -----------------------------------------------------------------------

!           THIRD-BODY LIST
            READ(ncchem,*)
            READ(ncchem,'(I5)')nbody
            DO jbody = 1, nbody
                READ(ncchem,'(I5,3X,A)')ibody,bdysym(ibody)
            END DO

!           THIRD-BODY STEP-LIST
            READ(ncchem,*)
            IF(nbody > 0) THEN
                DO jstep = 1, nstep
                    READ(ncchem,'(I5,I8)')istep,mblist(istep)
                END DO
            END IF

!           THIRD-BODY EFFICIENCIES
            READ(ncchem,*)
            DO jbody = 1,nbody
                DO jspec = 1,nspec
                    READ(ncchem,'(2I5,1PE12.4)')ibody,ispec,effy3b(ispec,ibody)
                END DO
            END DO

!           -----------------------------------------------------------------------

!           GIBBS STEP-LIST
            READ(ncchem,*)
            READ(ncchem,'(I5)')ngibb
            IF(ngibb > 0) THEN
                DO jstep = 1, nstep
                    READ(ncchem,'(I5,I8)')istep,mglist(istep)
                END DO
            END IF

!           -----------------------------------------------------------------------

!           LINDEMANN STEP-LIST
            READ(ncchem,*)
            READ(ncchem,'(I5)')nlind
            IF(nlind > 0) THEN
                DO jstep = 1, nstep
                    READ(ncchem,'(I5,I8)')istep,mllist(istep)
                END DO
            END IF

!           LINDEMANN STEP RATE DATA
            READ(ncchem,*)
            DO jlind = 1, nlind
                READ(ncchem,'(I5,4(1PE12.4))')ilind,(rclind(icp,ilind),icp=1,4)
            END DO

!           -----------------------------------------------------------------------

!           TROE FORM STEP-LIST
            READ(ncchem,*)
            READ(ncchem,'(I5)')ntroe
            IF(ntroe > 0) THEN
                DO jstep = 1, nstep
                    READ(ncchem,'(I5,I8)')istep,mtlist(istep)
                END DO
            END IF

!           TROE FORM STEP RATE DATA
            READ(ncchem,*)
            DO jtroe = 1, ntroe
                READ(ncchem,'(I5,6(1PE12.4))')itroe,(rctroe(icp,itroe),icp=1,6)
                READ(ncchem,'(I5,6(1PE12.4))')itroe,(rctroe(icp,itroe),icp=7,12)
            END DO

!           -----------------------------------------------------------------------

!           SRI FORM STEP-LIST
            READ(ncchem,*)
            READ(ncchem,'(I5)')nsrif
            IF(nsrif > 0) THEN
                DO jstep = 1, nstep
                    READ(ncchem,'(I5,I8)')istep,mslist(istep)
                END DO
            END IF

!           SRI FORM STEP RATE DATA
            READ(ncchem,*)
            DO jsrif = 1, nsrif
                READ(ncchem,'(I5,3(1PE12.4))')isrif,(rcsrif(icp,isrif),icp=1,3)
                READ(ncchem,'(I5,5(1PE12.4))')isrif,(rcsrif(icp,isrif),icp=4,8)
            END DO

!           -----------------------------------------------------------------------

!           READ END-OF-FILE LINE
            READ(ncchem,*)
        CLOSE(ncchem)

!       =======================================================================

!       COMPRESS THE DATA
!       -----------------

!       NO OF SPECIES
        ncount = 1
        parray(ncount) = REAL(nspec,kind=8)

!       SPECIES STRINGS
!       RSC 29-DEC-2006 BUG FIX INDEXING
        DO ispec = 1,nspec
            DO icp = 1, nspstr
                ncount = ncount + 1
                parray(ncount) = REAL(ICHAR(spcsym(ispec)(icp:icp)),kind=8)
            END DO
        END DO

!       SPECIES DATA
        ncount = ncount + 1
        parray(ncount) = prefgb
        DO ispec = 1,nspec
            ncount = ncount + 1
            parray(ncount) = wmolar(ispec)
            ncount = ncount + 1
            parray(ncount) = clewis(ispec)
            ncount = ncount + 1
            parray(ncount) = REAL(ntint(ispec),kind=8)
            DO itint = 1,ntint(ispec)
                ncount = ncount + 1
                parray(ncount) = tintlo(itint,ispec)
                ncount = ncount + 1
                parray(ncount) = tinthi(itint,ispec)
                ncount = ncount + 1
                parray(ncount) = REAL(ncofcp(itint,ispec),kind=8)
                DO icp = 1,ncofcp(itint,ispec)
                    ncount = ncount + 1
                    parray(ncount) = amolcp(icp,itint,ispec)
                END DO
            END DO
        END DO

!       -----------------------------------------------------------------------

!       NO OF STEPS
        ncount = ncount + 1
        parray(ncount) = REAL(nstep,kind=8)

!       STEP RATE DATA
        DO istep = 1, nstep
            DO icp = 1,3
                ncount = ncount + 1
                parray(ncount) = rparam(icp,istep)
            END DO
        END DO

!       STEP SPECIES-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            parray(ncount) = REAL(nsslen(istep),kind=8)
            DO ispec = 1, nsslen(istep)
                ncount = ncount + 1
                parray(ncount) = REAL(nsspec(ispec,istep),kind=8)
            END DO
        END DO

!       STEP REACTANT-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            parray(ncount) = REAL(nrslen(istep),kind=8)
            DO ispec = 1, nrslen(istep)
                ncount = ncount + 1
                parray(ncount) = REAL(nrspec(ispec,istep),kind=8)
            END DO
        END DO

!       STEP PRODUCT-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            parray(ncount) = REAL(npslen(istep),kind=8)
            DO ispec = 1, npslen(istep)
                ncount = ncount + 1
                parray(ncount) = REAL(npspec(ispec,istep),kind=8)
            END DO
        END DO

!       STEP REACTANT COEFFICIENT-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            parray(ncount) = REAL(nrclen(istep),kind=8)
            DO ispec = 1, nrclen(istep)
                ncount = ncount + 1
                parray(ncount) = REAL(nrcpec(ispec,istep),kind=8)
                ncount = ncount + 1
                parray(ncount) = crspec(ispec,istep)
            END DO
        END DO

!       STEP PRODUCT COEFFICIENT-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            parray(ncount) = REAL(npclen(istep),kind=8)
            DO ispec = 1, npclen(istep)
                ncount = ncount + 1
                parray(ncount) = REAL(npcpec(ispec,istep),kind=8)
                ncount = ncount + 1
                parray(ncount) = cpspec(ispec,istep)
            END DO
        END DO

!       SPECIES DELTA-LIST
        DO istep = 1,nstep
            DO ispec = 1, nsslen(istep)
                ncount = ncount + 1
                parray(ncount) = diffmu(ispec,istep)
            END DO
        END DO

!       -----------------------------------------------------------------------

!       THIRD-BODY LIST
        ncount = ncount + 1
        parray(ncount) = REAL(nbody,kind=8)
        DO ibody = 1, nbody
            DO icp = 1, nspstr
                ncount = ncount + 1
                parray(ncount) = REAL(ICHAR(bdysym(ibody)(icp:icp)),kind=8)
            END DO
        END DO

!       THIRD-BODY STEP-LIST
        IF(nbody > 0)THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                parray(ncount) = REAL(mblist(istep),kind=8)
            END DO
        END IF

!       THIRD-BODY EFFICIENCIES
        DO ibody = 1, nbody
            DO ispec = 1,nspec
                ncount = ncount + 1
                parray(ncount) = effy3b(ispec,ibody)
            END DO
        END DO

!       -----------------------------------------------------------------------

!       GIBBS STEP-LIST
        ncount = ncount + 1
        parray(ncount) = REAL(ngibb,kind=8)
        IF(ngibb > 0) THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                parray(ncount) = REAL(mglist(istep),kind=8)
            END DO
        END IF

!       -----------------------------------------------------------------------

!       LINDEMANN STEP-LIST
        ncount = ncount + 1
        parray(ncount) = REAL(nlind,kind=8)
        IF(nlind > 0) THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                parray(ncount) = REAL(mllist(istep),kind=8)
            END DO
        END IF

!       LINDEMANN STEP RATE DATA
        DO ilind = 1, nlind
            DO icp = 1,4
                ncount = ncount + 1
                parray(ncount) = rclind(icp,ilind)
            END DO
        END DO

!       -----------------------------------------------------------------------

!       TROE FORM STEP-LIST
        ncount = ncount + 1
        parray(ncount) = REAL(ntroe,kind=8)
        IF(ntroe > 0) THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                parray(ncount) = REAL(mtlist(istep),kind=8)
            END DO
        END IF

!       TROE FORM STEP RATE DATA
        DO itroe = 1, ntroe
            DO icp = 1,12
                ncount = ncount + 1
                parray(ncount) = rctroe(icp,itroe)
            END DO
        END DO

!       -----------------------------------------------------------------------

!       SRI FORM STEP-LIST
        ncount = ncount + 1
        parray(ncount) = REAL(nsrif,kind=8)
        IF(nsrif > 0) THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                parray(ncount) = REAL(mslist(istep),kind=8)
            END DO
        END IF

!       SRI FORM STEP RATE DATA
        DO isrif = 1, nsrif
            DO icp = 1,8
                ncount = ncount + 1
                parray(ncount) = rcsrif(icp,isrif)
            END DO
        END DO

!       -----------------------------------------------------------------------

!       BROADCAST COUNTER
        pcount = REAL(ncount,kind=8)

!       CHECK BROADCAST COUNTER AGAINST BROADCAST ARRAY SIZE
        IF(ncount > nparay) THEN
            WRITE(ncrept,*)
            WRITE(ncrept,*)'Warning: CHEMIN: broadcast array size too small'
            WRITE(ncrept,*)'Actual: ',ncount,'; Available: ',nparay
            WRITE(ncrept,*)
        END IF

!       =======================================================================

    END IF

!   =========================================================================

!   BROADCAST (OR RECEIVE) THE COUNTER
!   ----------------------------------
    call p_bcst(pcount,1,1,iroot)
    ncount = nint(pcount)

!   BROADCAST (OR RECEIVE) THE DATA
!   -------------------------------
    call p_bcst(parray,nparay,ncount,iroot)

!   =========================================================================

!   NOT THE LOWEST-RANKED PROCESSOR
    IF(iproc /= 0) THEN

!       =======================================================================

!       UNCOMPRESS THE DATA
!       -------------------

!       NO OF SPECIES
        ncount = 1
        nspec = nint(parray(ncount))

!       SPECIES STRINGS
!       RSC 29-DEC-2006 BUG FIX INDEXING
        DO ispec = 1,nspec
            DO icp = 1, nspstr
                ncount = ncount + 1
                spcsym(ispec)(icp:icp) = CHAR(nint(parray(ncount)))
            END DO
        END DO

!       SPECIES DATA
        ncount = ncount + 1
        prefgb = parray(ncount)
        DO ispec = 1,nspec
            ncount = ncount + 1
            wmolar(ispec) = parray(ncount)
            ncount = ncount + 1
            clewis(ispec) = parray(ncount)
            ncount = ncount + 1
            ntint(ispec) = nint(parray(ncount))
            DO itint = 1,ntint(ispec)
                ncount = ncount + 1
                tintlo(itint,ispec) = parray(ncount)
                ncount = ncount + 1
                tinthi(itint,ispec) = parray(ncount)
                ncount = ncount + 1
                ncofcp(itint,ispec) = nint(parray(ncount))
                DO icp = 1,ncofcp(itint,ispec)
                    ncount = ncount + 1
                    amolcp(icp,itint,ispec) = parray(ncount)
                END DO
            END DO
        END DO

!       -----------------------------------------------------------------------

!       NO OF STEPS
        ncount = ncount + 1
        nstep = nint(parray(ncount))

!       STEP RATE DATA
        DO istep = 1, nstep
            DO icp = 1,3
                ncount = ncount + 1
                rparam(icp,istep) = parray(ncount)
            END DO
        END DO

!       STEP SPECIES-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            nsslen(istep) = nint(parray(ncount))
            DO ispec = 1, nsslen(istep)
                ncount = ncount + 1
                nsspec(ispec,istep) = nint(parray(ncount))
            END DO
        END DO

!       STEP REACTANT-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            nrslen(istep) = nint(parray(ncount))
            DO ispec = 1, nrslen(istep)
                ncount = ncount + 1
                nrspec(ispec,istep) = nint(parray(ncount))
            END DO
        END DO

!       STEP PRODUCT-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            npslen(istep) = nint(parray(ncount))
            DO ispec = 1, npslen(istep)
                ncount = ncount + 1
                npspec(ispec,istep) = nint(parray(ncount))
            END DO
        END DO

!       STEP REACTANT COEFFICIENT-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            nrclen(istep) = nint(parray(ncount))
            DO ispec = 1, nrclen(istep)
                ncount = ncount + 1
                nrcpec(ispec,istep) = nint(parray(ncount))
                ncount = ncount + 1
                crspec(ispec,istep) = parray(ncount)
            END DO
        END DO

!       STEP PRODUCT COEFFICIENT-LIST
        DO istep = 1,nstep
            ncount = ncount + 1
            npclen(istep) = nint(parray(ncount))
            DO ispec = 1, npclen(istep)
                ncount = ncount + 1
                npcpec(ispec,istep) = nint(parray(ncount))
                ncount = ncount + 1
                cpspec(ispec,istep) = parray(ncount)
            END DO
        END DO

!       SPECIES DELTA-LIST
        DO istep = 1,nstep
            DO ispec = 1, nsslen(istep)
                ncount = ncount + 1
                diffmu(ispec,istep) = parray(ncount)
            END DO
        END DO

!       -----------------------------------------------------------------------

!       THIRD-BODY LIST
        ncount = ncount + 1
        nbody = nint(parray(ncount))
        DO ibody = 1, nbody
            DO icp = 1, nspstr
                ncount = ncount + 1
                bdysym(ibody)(icp:icp) = CHAR(nint(parray(ncount)))
            END DO
        END DO

!       THIRD-BODY STEP-LIST
        IF(nbody > 0) THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                mblist(istep) = nint(parray(ncount))
            END DO
        END IF

!       THIRD-BODY EFFICIENCIES
        DO ibody = 1, nbody
            DO ispec = 1,nspec
                ncount = ncount + 1
                effy3b(ispec,ibody) = parray(ncount)
            END DO
        END DO

!       -----------------------------------------------------------------------

!       GIBBS STEP-LIST
        ncount = ncount + 1
        ngibb = nint(parray(ncount))
        IF(ngibb > 0) THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                mglist(istep) = nint(parray(ncount))
            END DO
        END IF

!       -----------------------------------------------------------------------

!       LINDEMANN STEP-LIST
        ncount = ncount + 1
        nlind = nint(parray(ncount))
        IF(nlind > 0) THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                mllist(istep) = nint(parray(ncount))
            END DO
        END IF

!       LINDEMANN STEP RATE DATA
        DO ilind = 1, nlind
            DO icp = 1,4
                ncount = ncount + 1
                rclind(icp,ilind) = parray(ncount)
            END DO
        END DO

!       -----------------------------------------------------------------------

!       TROE FORM STEP-LIST
        ncount = ncount + 1
        ntroe = nint(parray(ncount))
        IF(ntroe > 0) THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                mtlist(istep) = nint(parray(ncount))
            END DO
        END IF

!       TROE FORM STEP RATE DATA
        DO itroe = 1, ntroe
            DO icp = 1,12
                ncount = ncount + 1
                rctroe(icp,itroe) = parray(ncount)
            END DO
        END DO

!       -----------------------------------------------------------------------

!       SRI FORM STEP-LIST
        ncount = ncount + 1
        nsrif = nint(parray(ncount))
        IF(nsrif > 0) THEN
            DO istep = 1, nstep
                ncount = ncount + 1
                mslist(istep) = nint(parray(ncount))
            END DO
        END IF

!       SRI FORM STEP RATE DATA
        DO isrif = 1, nsrif
            DO icp = 1,8
                ncount = ncount + 1
                rcsrif(icp,isrif) = parray(ncount)
            END DO
        END DO

!       =======================================================================

    END IF

!   =========================================================================

!   EVALUATE DERIVED QUANTITIES
!   ===========================

!   NO OF ACTIVE SPECIES
    nspm1 = nspec - 1

!   -------------------------------------------------------------------------

!   CONVERT RATE PARAMETERS
    DO istep = 1,nstep
        rparam(1,istep) = LOG(rparam(1,istep))
        rparam(3,istep) = rparam(3,istep)/rguniv
    END DO

!   LINDEMANN STEP RATE DATA
    DO ilind = 1, nlind
        rclind(1,ilind) = LOG(rclind(1,ilind))
        rclind(3,ilind) = rclind(3,ilind)/rguniv
    END DO

!   TROE FORM STEP RATE DATA
    DO itroe = 1, ntroe
        rctroe(1,itroe) = LOG(rctroe(1,itroe))
        rctroe(3,itroe) = rctroe(3,itroe)/rguniv
        rctroe(5,itroe) = -one/rctroe(5,itroe)      !       VM: sign change based on Khalil
        rctroe(7,itroe) = -one/rctroe(7,itroe)
    END DO

!   SRI FORM STEP RATE DATA
    DO isrif = 1, nsrif
        rcsrif(1,isrif) = LOG(rcsrif(1,isrif))
        rcsrif(3,isrif) = rcsrif(3,isrif)/rguniv
        rcsrif(5,isrif) = -rcsrif(5,isrif)
        rcsrif(6,isrif) = -one/rcsrif(6,isrif)
    END DO

!   -------------------------------------------------------------------------

!   STOICHIOMETRIC COEFFICIENTS TIMES MOLAR MASS
    DO istep = 1,nstep
        DO ispec = 1,nsslen(istep)
            jspec = nsspec(ispec,istep)
            diffmw(ispec,istep) = diffmu(ispec,istep)*wmolar(jspec)
        END DO
    END DO

!   -------------------------------------------------------------------------

!   RECIPROCAL OF MOLAR MASS
!   SPECIFIC GAS CONSTANT
    DO ispec = 1,nspec
        ovwmol(ispec) = one/wmolar(ispec)
        rgspec(ispec) = rguniv*ovwmol(ispec)
    END DO

!   -------------------------------------------------------------------------

!   SPECIFIC HEAT CAPACITY
!   ----------------------

!   NUMBER OF CP POLYNOMIAL COEFFICIENTS
!   NUMBER OF CP POLYNOMIAL COEFFICIENTS MINUS ONE
!   INDEX NUMBERS OF ENTHALPY AND ENTROPY COEFFICIENTS
    DO ispec = 1,nspec
        DO itint = 1,ntint(ispec)
            ncpoly(itint,ispec) = ncofcp(itint,ispec)-2
            ncpom1(itint,ispec) = ncpoly(itint,ispec)-1
            ncenth(itint,ispec) = ncofcp(itint,ispec)-1
            ncenpy(itint,ispec) = ncofcp(itint,ispec)
        END DO
    END DO

!   SPECIFIC HEAT CAPACITY CP PER UNIT MASS
    DO ispec = 1,nspec
        DO itint = 1,ntint(ispec)
            DO icp = 1,ncofcp(itint,ispec)
                amascp(icp,itint,ispec) = amolcp(icp,itint,ispec)*rguniv*ovwmol(ispec)
            END DO
        END DO
    END DO

!   COEFFICIENTS FOR TEMPERATURE
    DO ispec = 1,nspec
        DO itint = 1,ntint(ispec)
            amasct(1,itint,ispec) = amascp(1,itint,ispec) - rgspec(ispec)
            DO icp = 2,ncpoly(itint,ispec)
                amasct(icp,itint,ispec) = amascp(icp,itint,ispec)/REAL(icp,kind=8)
            END DO
            amasct(ncenth(itint,ispec),itint,ispec) = zero
            amasct(ncenpy(itint,ispec),itint,ispec) = zero
        END DO
    END DO

!   COEFFICIENTS FOR ENTHALPY PER UNIT MASS
    DO ispec = 1,nspec
        DO itint = 1,ntint(ispec)
            amasch(1,itint,ispec) = amascp(1,itint,ispec)
            DO icp = 2,ncpoly(itint,ispec)
                amasch(icp,itint,ispec) = amascp(icp,itint,ispec)/REAL(icp,kind=8)
            END DO
            amasch(ncenth(itint,ispec),itint,ispec)  &
                = amascp(ncenth(itint,ispec),itint,ispec)
            amasch(ncenpy(itint,ispec),itint,ispec) = zero
        END DO
    END DO

!   COEFFICIENTS FOR ENTROPY PER UNIT MASS
    DO ispec = 1,nspec
        DO itint = 1,ntint(ispec)
            amascs(1,itint,ispec) = amascp(1,itint,ispec)
            DO icp = 2,ncpoly(itint,ispec)
                amascs(icp,itint,ispec) = amascp(icp,itint,ispec) /REAL(icp-1,kind=8)
            END DO
            amascs(ncenth(itint,ispec),itint,ispec) = zero
            amascs(ncenpy(itint,ispec),itint,ispec)  &
                = amascp(ncenpy(itint,ispec),itint,ispec)
        END DO
    END DO

!   COEFFICIENTS FOR GIBBS FUNCTION PER MOLE
!   ACTUALLY GIBBS/(R^0 T) WITH PRESSURE TERM
    giblet = LOG(prefgb/rguniv)
    DO ispec = 1,nspec
        DO itint = 1,ntint(ispec)
            amolgb(1,itint,ispec)= amolcp(ncenpy(itint,ispec),itint,ispec)  &
                - amolcp(1,itint,ispec) + giblet
            DO icp = 2,ncpoly(itint,ispec)
                amolgb(icp,itint,ispec) = amolcp(icp,itint,ispec) /REAL(icp*(icp-1),kind=8)
            END DO
            amolgb(ncenth(itint,ispec),itint,ispec)  &
                = amolcp(ncenth(itint,ispec),itint,ispec)
            amolgb(ncenpy(itint,ispec),itint,ispec) = amolcp(1,itint,ispec) - one
        END DO
    END DO

!   SPECIFIC HEAT CAPACITY CP PER MOLE
    DO ispec = 1,nspec
        DO itint = 1,ntint(ispec)
            DO icp = 1,ncofcp(itint,ispec)
                amolcp(icp,itint,ispec) = amolcp(icp,itint,ispec)*rguniv
            END DO
        END DO
    END DO

!   -------------------------------------------------------------------------

!   RECIPROCAL OF LEWIS NUMBER
    DO ispec = 1,nspec
        olewis(ispec) = one/clewis(ispec)
    END DO

!   -------------------------------------------------------------------------

!   TRANSPORT COEFFICIENTS
    rlamda = 0.70_8

!   CONDUCTIVITY COEFFICIENT
    alamda = alamdc*EXP(-rlamda*LOG(tlamda))

!   =========================================================================

    call ops_decl_const("rlamda", 1, "real(kind=8)", rlamda)
    call ops_decl_const("alamda", 1, "real(kind=8)", alamda)

!   =========================================================================

END SUBROUTINE chemin
