SUBROUTINE contin

! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:25:00

!     *************************************************************************

!     CONTIN
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     30-DEC-2003:  CREATED
!     11-JUL-2009:  RSC ADD A DUMP FORMAT SWITCH
!     01-MAY-2013:  RSC MIXTURE AVERAGED TRANSPORT
!     21-JUL-2013:  RSC RADIATION TREATMENT

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     INITIALISES RUN CONTROL DATA

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_espect
use com_senga
!     -------------------------------------------------------------------------


!     PARAMETER
!     =========
!     RSC 01-MAY-2013
CHARACTER (LEN=32) :: strmol
PARAMETER(strmol = 'Molecular transport control data')


!     LOCAL DATA
!     ==========
real(kind=8) :: pcount
integer(kind=4) :: ncount
integer(kind=4) :: ispec,jspec
integer(kind=4) :: ic
integer(kind=4) :: iroot
!     RSC 01-MAY-2013
CHARACTER (LEN=50) :: eoflin


!     BEGIN
!     =====

!     =========================================================================

!     DATA IS READ BY THE LOWEST-RANKED PROCESSOR AND BROADCAST TO THE OTHERS
!     -----------------------------------------------------------------------

!     SET ID OF ROOT (BROADCASTING) PROCESSOR
iroot = 0

!     INITIALISE THE BROADCAST COUNTER
pcount = zero

!     INITIALISE THE BROADCAST ARRAY
DO ncount = 1,nparay
  parray(ncount) = zero
END DO

!     =========================================================================

!     LOWEST-RANKED PROCESSOR
IF(iproc == 0)THEN

!       =======================================================================

!       READ THE RUN CONTROL DATA FILE
!       ------------------------------
!       -----------------------------------------------------------------------

  OPEN(UNIT=nccont,FILE=fncont,STATUS='OLD',FORM='FORMATTED')

!       -----------------------------------------------------------------------

!       READ AND IGNORE THE HEADER
  READ(nccont,*)
  READ(nccont,*)
  READ(nccont,*)
  READ(nccont,*)
  READ(nccont,*)
!       AND THE SUBSEQUENT BLANK LINE
  READ(nccont,*)

!       -----------------------------------------------------------------------

!       GRID SIZE DATA
!       --------------
!       REQUIRED GLOBAL GRID SIZE X,Y,Z
  READ(nccont,*)
  READ(nccont,*)xgdlen,ygdlen,zgdlen
  READ(nccont,*)

!       REQUIRED GLOBAL GRID SIZE NX,NY,NZ
  READ(nccont,*)
  READ(nccont,*)nxgreq,nygreq,nzgreq
  READ(nccont,*)

!       REQUIRED NUMBER OF PROCESSORS X,Y,Z
  READ(nccont,*)
  READ(nccont,*)nxpreq,nypreq,nzpreq
  READ(nccont,*)

!       TIME STEPPING DATA
!       ------------------
!       TIME STEP, START STEP, NO OF STEPS, STEP SWITCH
  READ(nccont,*)
  READ(nccont,*)tstep,ntime1,ntime,nstpsw
  READ(nccont,*)

!       INTERVALS BETWEEN DUMPS,REPORTS,STATISTICS
!       ------------------------------------------
!       DUMP OUTPUT FORMAT
!       RSC 11-JUL-2009
  READ(nccont,*)
  READ(nccont,*)ntdump,ntrept,ntstat,ndofmt
  READ(nccont,*)

!       COLD START SWITCH (0=COLD START, 1=RESTART)
!       -----------------
!       DUMP INPUT FORMAT
!       RSC 11-JUL-2009
  READ(nccont,*)
  READ(nccont,*)ncdmpi,ndifmt
  READ(nccont,*)

!       INITIAL TURBULENCE (0=INITIAL TURBULENCE OFF, 1=ON)
!       ------------------
!       RANDOM SEED, SPECTRUM PARAMETERS (FOUR)
  READ(nccont,*)
  READ(nccont,*)
  READ(nccont,*)inturb,inseed,(sparam(ic),ic=1,nsparm)
  READ(nccont,*)

!       FLAME START OPTION (0=NO FLAME, 1=FLAME)
!       ------------------
  READ(nccont,*)
  READ(nccont,*)inflam
  READ(nccont,*)

!       DEFAULT INITIAL CONDITIONS
!       --------------------------
  READ(nccont,*)
  READ(nccont,*)
  READ(nccont,*)prin,trin,urin,vrin,wrin
  READ(nccont,*)
  READ(nccont,*)nspreq
  DO ispec = 1, nspreq
    READ(nccont,*)jspec,yrin(ispec)
  END DO
  READ(nccont,*)

!       GLOBAL BOUNDARY CONDITION TYPES
!       -------------------------------
  READ(nccont,*)
  READ(nccont,*)
  READ(nccont,*)
  READ(nccont,*)
  READ(nccont,*)ngbcxl,(nxlprm(ic),ic=1,nbcpri), (rxlprm(ic),ic=1,nbcprr)
! FY - ADDED FOR READING OF EXTRA VARIABLES FOR NON-REFLECTING INFLOW
  IF(ngbcxl == nsbci4)THEN
    READ(nccont,*)(rxlprc(ic),ic=1,nbcprc)
  END IF
  READ(nccont,*)ngbcxr,(nxrprm(ic),ic=1,nbcpri), (rxrprm(ic),ic=1,nbcprr)
  READ(nccont,*)ngbcyl,(nylprm(ic),ic=1,nbcpri), (rylprm(ic),ic=1,nbcprr)
  READ(nccont,*)ngbcyr,(nyrprm(ic),ic=1,nbcpri), (ryrprm(ic),ic=1,nbcprr)
  READ(nccont,*)ngbczl,(nzlprm(ic),ic=1,nbcpri), (rzlprm(ic),ic=1,nbcprr)
  READ(nccont,*)ngbczr,(nzrprm(ic),ic=1,nbcpri), (rzrprm(ic),ic=1,nbcprr)
  READ(nccont,*)

!       -----------------------------------------------------------------------

!       RSC 01-MAY-2013
!       READ MOLECULAR TRANSPORT CONTROL DATA
!       RSC 21-JUL-2013
!       READ RADIATION CONTROL DATA
  nfmavt = 0
  nfmixw = 0
  nfmixp = 0
  nfmsor = 0
  nfmduf = 0
  wmltdr = 0.0_8
  nfradn = 0
  trefrn = 3.0_8
  READ(nccont,'(A)')eoflin
  IF(eoflin(1:32) == strmol)THEN
    READ(nccont,*)
    READ(nccont,*)nfmavt,nfmixw,nfmixp,nfmsor,nfmduf,wmltdr
    READ(nccont,*)
    READ(nccont,*)
    READ(nccont,*)
    READ(nccont,*)nfradn,trefrn
    READ(nccont,*)
  END IF

!       READ END-OF-FILE LINE
  READ(nccont,'(A)')eoflin

  CLOSE(nccont)

!       =======================================================================

!       COMPRESS THE DATA
!       -----------------

!       REQUIRED GLOBAL DOMAIN SIZE (x,y,z)
  ncount = 1
  parray(ncount) = xgdlen
  ncount = ncount + 1
  parray(ncount) = ygdlen
  ncount = ncount + 1
  parray(ncount) = zgdlen

!       TIME STEP, START STEP, NO OF STEPS, STEP SWITCH
  ncount = ncount + 1
  parray(ncount) = tstep
  ncount = ncount + 1
  parray(ncount) = REAL(ntime1,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(ntime,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(nstpsw,kind=8)

!       INTERVALS BETWEEN DUMPS,REPORTS,STATISTICS; DUMP OUTPUT SWITCH
!       RSC 11-JUL-2009
  ncount = ncount + 1
  parray(ncount) = REAL(ntdump,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(ntrept,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(ntstat,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(ndofmt,kind=8)

!       COLD START SWITCH; DUMP INPUT SWITCH
!       RSC 11-JUL-2009
  ncount = ncount + 1
  parray(ncount) = REAL(ncdmpi,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(ndifmt,kind=8)

!       INITIAL TURBULENCE
  ncount = ncount + 1
  parray(ncount) = REAL(inturb,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(inseed,kind=8)
  DO ispec = 1, nsparm
    ncount = ncount + 1
    parray(ncount) = sparam(ispec)
  END DO

!       FLAME START OPTION
  ncount = ncount + 1
  parray(ncount) = REAL(inflam,kind=8)

!       DEFAULT INITIAL CONDITIONS
  ncount = ncount + 1
  parray(ncount) = prin
  ncount = ncount + 1
  parray(ncount) = trin
  ncount = ncount + 1
  parray(ncount) = urin
  ncount = ncount + 1
  parray(ncount) = vrin
  ncount = ncount + 1
  parray(ncount) = wrin
  DO ispec = 1, nspcmx
    ncount = ncount + 1
    parray(ncount) = yrin(ispec)
  END DO

!       GLOBAL BOUNDARY CONDITION TYPES
  ncount = ncount + 1
  parray(ncount) = REAL(ngbcxl,kind=8)
  DO ic = 1, nbcpri
    ncount = ncount + 1
    parray(ncount) = REAL(nxlprm(ic),kind=8)
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    parray(ncount) = rxlprm(ic)
  END DO
  IF(ngbcxl == nsbci4)THEN ! FY - ADDED FOR NSBCI4
    DO ic = 1, nbcprc
      ncount = ncount + 1
      parray(ncount) = rxlprc(ic)
    END DO
  END IF
  ncount = ncount + 1
  parray(ncount) = REAL(ngbcxr,kind=8)
  DO ic = 1, nbcpri
    ncount = ncount + 1
    parray(ncount) = REAL(nxrprm(ic),kind=8)
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    parray(ncount) = rxrprm(ic)
  END DO
  ncount = ncount + 1
  parray(ncount) = REAL(ngbcyl,kind=8)
  DO ic = 1, nbcpri
    ncount = ncount + 1
    parray(ncount) = REAL(nylprm(ic),kind=8)
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    parray(ncount) = rylprm(ic)
  END DO
  ncount = ncount + 1
  parray(ncount) = REAL(ngbcyr,kind=8)
  DO ic = 1, nbcpri
    ncount = ncount + 1
    parray(ncount) = REAL(nyrprm(ic),kind=8)
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    parray(ncount) = ryrprm(ic)
  END DO
  ncount = ncount + 1
  parray(ncount) = REAL(ngbczl,kind=8)
  DO ic = 1, nbcpri
    ncount = ncount + 1
    parray(ncount) = REAL(nzlprm(ic),kind=8)
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    parray(ncount) = rzlprm(ic)
  END DO
  ncount = ncount + 1
  parray(ncount) = REAL(ngbczr,kind=8)
  DO ic = 1, nbcpri
    ncount = ncount + 1
    parray(ncount) = REAL(nzrprm(ic),kind=8)
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    parray(ncount) = rzrprm(ic)
  END DO

!       MIXTURE AVERAGED TRANSPORT
!       RSC 05-MAY-2013
  ncount = ncount + 1
  parray(ncount) = REAL(nfmavt,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(nfmixw,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(nfmixp,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(nfmsor,kind=8)
  ncount = ncount + 1
  parray(ncount) = REAL(nfmduf,kind=8)
  ncount = ncount + 1
  parray(ncount) = wmltdr

!       RADIATION TREATMENT
!       RSC 21-JUL-2013
  ncount = ncount + 1
  parray(ncount) = REAL(nfradn,kind=8)
  ncount = ncount + 1
  parray(ncount) = trefrn

!       BROADCAST COUNTER
  pcount = REAL(ncount,kind=8)

!       CHECK BROADCAST COUNTER AGAINST BROADCAST ARRAY SIZE
  IF(ncount > nparay)THEN
    WRITE(ncrept,*)'Warning: INCONT: broadcast array size too small'
  END IF

!       =======================================================================

END IF

!     =========================================================================

!     BROADCAST (OR RECEIVE) THE COUNTER
!     ----------------------------------
CALL p_bcst(pcount,1,1,iroot)
ncount = nint(pcount)

!     BROADCAST (OR RECEIVE) THE DATA
!     -------------------------------
CALL p_bcst(parray,nparay,ncount,iroot)

!     =========================================================================

!     NOT THE LOWEST-RANKED PROCESSOR
IF(iproc /= 0)THEN

!       =======================================================================

!       UNCOMPRESS THE DATA
!       -------------------

!       REQUIRED GLOBAL GRID SIZE X,Y,Z
  ncount = 1
  xgdlen = parray(ncount)
  ncount = ncount + 1
  ygdlen = parray(ncount)
  ncount = ncount + 1
  zgdlen = parray(ncount)

!       TIME STEP, START STEP, NO OF STEPS, STEP SWITCH
  ncount = ncount + 1
  tstep = parray(ncount)
  ncount = ncount + 1
  ntime1 = nint(parray(ncount))
  ncount = ncount + 1
  ntime = nint(parray(ncount))
  ncount = ncount + 1
  nstpsw = nint(parray(ncount))

!       INTERVALS BETWEEN DUMPS,REPORTS,STATISTICS; DUMP OUTPUT SWITCH
!       RSC 11-JUL-2009
  ncount = ncount + 1
  ntdump = nint(parray(ncount))
  ncount = ncount + 1
  ntrept = nint(parray(ncount))
  ncount = ncount + 1
  ntstat = nint(parray(ncount))
  ncount = ncount + 1
  ndofmt = nint(parray(ncount))

!       COLD START SWITCH; DUMP INPUT SWITCH
!       RSC 11-JUL-2009
  ncount = ncount + 1
  ncdmpi = nint(parray(ncount))
  ncount = ncount + 1
  ndifmt = nint(parray(ncount))

!       INITIAL TURBULENCE
  ncount = ncount + 1
  inturb = nint(parray(ncount))
  ncount = ncount + 1
  inseed = nint(parray(ncount))
  DO ispec = 1, nsparm
    ncount = ncount + 1
    sparam(ispec) = parray(ncount)
  END DO

!       FLAME START OPTION
  ncount = ncount + 1
  inflam = nint(parray(ncount))

!       DEFAULT INITIAL CONDITIONS
  ncount = ncount + 1
  prin = parray(ncount)
  ncount = ncount + 1
  trin = parray(ncount)
  ncount = ncount + 1
  urin = parray(ncount)
  ncount = ncount + 1
  vrin = parray(ncount)
  ncount = ncount + 1
  wrin = parray(ncount)
  DO ispec = 1, nspcmx
    ncount = ncount + 1
    yrin(ispec) = parray(ncount)
  END DO

!       GLOBAL BOUNDARY CONDITION TYPES
  ncount = ncount + 1
  ngbcxl = nint(parray(ncount))
  DO ic = 1, nbcpri
    ncount = ncount + 1
    nxlprm(ic) = nint(parray(ncount))
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    rxlprm(ic) = parray(ncount)
  END DO
  IF(ngbcxl == nsbci4)THEN ! FY - ADDED FOR NSBCI4
    DO ic = 1, nbcprc
      ncount = ncount + 1
      rxlprc(ic) = parray(ncount)
    END DO
  END IF
  ncount = ncount + 1
  ngbcxr = nint(parray(ncount))
  DO ic = 1, nbcpri
    ncount = ncount + 1
    nxrprm(ic) = nint(parray(ncount))
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    rxrprm(ic) = parray(ncount)
  END DO
  ncount = ncount + 1
  ngbcyl = nint(parray(ncount))
  DO ic = 1, nbcpri
    ncount = ncount + 1
    nylprm(ic) = nint(parray(ncount))
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    rylprm(ic) = parray(ncount)
  END DO
  ncount = ncount + 1
  ngbcyr = nint(parray(ncount))
  DO ic = 1, nbcpri
    ncount = ncount + 1
    nyrprm(ic) = nint(parray(ncount))
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    ryrprm(ic) = parray(ncount)
  END DO
  ncount = ncount + 1
  ngbczl = nint(parray(ncount))
  DO ic = 1, nbcpri
    ncount = ncount + 1
    nzlprm(ic) = nint(parray(ncount))
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    rzlprm(ic) = parray(ncount)
  END DO
  ncount = ncount + 1
  ngbczr = nint(parray(ncount))
  DO ic = 1, nbcpri
    ncount = ncount + 1
    nzrprm(ic) = nint(parray(ncount))
  END DO
  DO ic = 1, nbcprr
    ncount = ncount + 1
    rzrprm(ic) = parray(ncount)
  END DO

!       MIXTURE AVERAGED TRANSPORT
!       RSC 05-MAY-2013
  ncount = ncount + 1
  nfmavt = nint(parray(ncount))
  ncount = ncount + 1
  nfmixw = nint(parray(ncount))
  ncount = ncount + 1
  nfmixp = nint(parray(ncount))
  ncount = ncount + 1
  nfmsor = nint(parray(ncount))
  ncount = ncount + 1
  nfmduf = nint(parray(ncount))
  ncount = ncount + 1
  wmltdr = parray(ncount)

!       RADIATION TREATMENT
!       RSC 21-JUL-2013
  ncount = ncount + 1
  nfradn = nint(parray(ncount))
  ncount = ncount + 1
  trefrn = parray(ncount)

!       =======================================================================

END IF

!     =========================================================================


RETURN
END SUBROUTINE contin
