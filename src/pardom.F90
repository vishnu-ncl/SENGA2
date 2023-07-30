SUBROUTINE pardom

    use OPS_Fortran_Reference

    use OPS_CONSTANTS
    use, intrinsic :: ISO_C_BINDING

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
    integer(kind=4) :: nextra
    integer(kind=4) :: icproc

!   BEGIN
!   =====

!   =========================================================================

!   INITIALISE PARALLEL MESSAGE PASSING
!   ===================================
!   INITIALISE
#ifndef OPS_MPI
   CALL p_init
#endif

!   GET THE NUMBER OF PROCESSORS NPROC
!   AND THE LOCAL PROCESSOR ID IPROC: IPROC=0,NPROC-1
#ifdef OPS_MPI
    call p_size(nproc,iproc)
#else
    nproc = 1
    iproc = 0
#endif

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

    nynode = nyglbl/nyproc
    nextra = MOD(nyglbl,nyproc)
    nextra = nyproc-nextra
    nyprm1 = nyproc-1
    DO icproc = 0,nyprm1
        npmapy(icproc) = nynode
        IF(icproc >= nextra) npmapy(icproc) = nynode+1
   END DO
   nynode = npmapy(iyproc)

    nznode = nzglbl/nzproc
    nextra = MOD(nzglbl,nzproc)
    nextra = nzproc-nextra
    nzprm1 = nzproc-1
    DO icproc = 0,nzprm1
        npmapz(icproc) = nznode
        IF(icproc >= nextra) npmapz(icproc) = nznode+1
    END DO
   nznode = npmapz(izproc)

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

END SUBROUTINE pardom
