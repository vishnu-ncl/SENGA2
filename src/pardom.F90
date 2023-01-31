SUBROUTINE pardom
 
    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

    use data_types
    use com_senga
    use com_ops_senga

!   *************************************************************************

!   PARDOM
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   01-OCT-1996:  CREATED
!   09-MAR-2003:  RSC UPDATED FOR SENGA2
!   27-SEP-2009:  RSC UPDATE INDEXING FOR FILTERS

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   CARRIES OUT PARALLEL DOMAIN DECOMPOSITION

!   *************************************************************************

!   GLOBAL DATA
!   ===========
!   -------------------------------------------------------------------------
!   -------------------------------------------------------------------------

!   LOCAL DATA
!   ===========
    INTEGER :: nextra
    INTEGER :: icproc

!   BEGIN
!   =====

!   =========================================================================

!   INITIALISE PARALLEL MESSAGE PASSING
!   ===================================
!   INITIALISE
   CALL p_init

!   GET THE NUMBER OF PROCESSORS NPROC
!   AND THE LOCAL PROCESSOR ID IPROC: IPROC=0,NPROC-1
#ifdef OPS_MPI
    call p_size(nproc,iproc)
#else
    nproc = 1
    iproc = 0
#endif

    write(*, '(a,i5,a,i5)') "Senga2 nproc: ",nproc,"  iproc: ",iproc

!   =========================================================================

!   GET XYZ COORDINATES OF LOCAL PROCESSOR
!   ======================================
!   IPROC = IXPROC+IYPROC*NXPROC+IZPROC*NXPROC*NYPROC
    ixproc = MOD(iproc,nxproc)
    iyproc = iproc/nxproc
    izproc = iyproc/nyproc
    iyproc = MOD(iyproc,nyproc)

!   =========================================================================

!   ASSIGN GRID POINTS TO PROCESSORS
!   ================================
!   STATIC LOAD-BALANCING:
!   LOAD IS AUTOMATICALLY BALANCED IF NXGLBL IS EXACTLY DIVISIBLE BY NXPROC
!   OTHERWISE ANY REMAINDER IS SPLIT EQUALLY BETWEEN THE HIGHER-RANKED
!   PROCESSORS IN EACH DIRECTION
    nxnode = nxglbl/nxproc
    nextra = MOD(nxglbl,nxproc)
    nextra = nxproc-nextra
    nxprm1 = nxproc-1
    DO icproc = 0,nxprm1
        npmapx(icproc) = nxnode
        IF(icproc >= nextra) npmapx(icproc) = nxnode+1
    END DO
   nxnode = npmapx(ixproc)
!    nxnode = nxsize

    nynode = nyglbl/nyproc
    nextra = MOD(nyglbl,nyproc)
    nextra = nyproc-nextra
    nyprm1 = nyproc-1
    DO icproc = 0,nyprm1
        npmapy(icproc) = nynode
        IF(icproc >= nextra) npmapy(icproc) = nynode+1
   END DO
   nynode = npmapy(iyproc)
!    nynode = nysize

    nznode = nzglbl/nzproc
    nextra = MOD(nzglbl,nzproc)
    nextra = nzproc-nextra
    nzprm1 = nzproc-1
    DO icproc = 0,nzprm1
        npmapz(icproc) = nznode
        IF(icproc >= nextra) npmapz(icproc) = nznode+1
    END DO
   nznode = npmapz(izproc)
!    nznode = nzsize

!   =========================================================================

!   MAP PROCESSORS ALONG CURRENT ROW
!   ================================
    DO icproc = 0,nxprm1
        nprocx(icproc) = icproc+nxproc*(iyproc+nyproc*izproc)
    END DO
    DO icproc = 0,nyprm1
        nprocy(icproc) = ixproc+nxproc*(icproc+nyproc*izproc)
    END DO
    DO icproc = 0,nzprm1
        nprocz(icproc) = ixproc+nxproc*(iyproc+nyproc*icproc)
    END DO

!   =========================================================================

!   IDENTIFY NEAREST NEIGHBOURS FOR MESSAGE-PASSING
!   ===============================================
    ixprom = ixproc-1
    IF(ixprom < 0) ixprom = nxproc-1
    ixprop = ixproc+1
    IF(ixprop >= nxproc) ixprop = 0
    ixprom = ixprom+nxproc*(iyproc+izproc*nyproc)
    ixprop = ixprop+nxproc*(iyproc+izproc*nyproc)

    iyprom = iyproc-1
    IF(iyprom < 0) iyprom = nyproc-1
    iyprop = iyproc+1
    IF(iyprop >= nyproc) iyprop = 0
    iyprom = ixproc+nxproc*(iyprom+izproc*nyproc)
    iyprop = ixproc+nxproc*(iyprop+izproc*nyproc)

    izprom = izproc-1
    IF(izprom < 0) izprom = nzproc-1
    izprop = izproc+1
    IF(izprop >= nzproc) izprop = 0
    izprom = ixproc+nxproc*(iyproc+izprom*nyproc)
    izprop = ixproc+nxproc*(iyproc+izprop*nyproc)

!   =========================================================================

!   SET MESSAGE TAGS
!   EACH TAG UNIQUELY IDENTIFIES THE SENDING AND RECEIVING PROCESS
!   FOR EACH PROCESSOR FOR EACH DIRECTION THERE ARE FOUR TAGS
!   CORRESPONDING TO SEND AND RECEIVE LEFT AND RIGHT
    itgxsl = iproc*nproc + ixprom
    itgxrl = ixprom*nproc + iproc
    itgxsr = iproc*nproc + ixprop
    itgxrr = ixprop*nproc + iproc
    itgysl = iproc*nproc + iyprom
    itgyrl = iyprom*nproc + iproc
    itgysr = iproc*nproc + iyprop
    itgyrr = iyprop*nproc + iproc
    itgzsl = iproc*nproc + izprom
    itgzrl = izprom*nproc + iproc
    itgzsr = iproc*nproc + izprop
    itgzrr = izprop*nproc + iproc

!   =========================================================================

END SUBROUTINE pardom
