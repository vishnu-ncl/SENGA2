MODULE com_senga

!     PHYDIM-------------------------------------------------------------------
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-02  Time: 14:22:03

!     PHYSICAL DIMENSIONS OF ARRAYS
!     -----------------------------
!     NOTE: ALL ARRAY SIZES MUST BE CONSISTENT
!     NXSIZE MUST BE >= NXGLBL/NXPROC
!     NYSIZE MUST BE >= NYGLBL/NYPROC
!     NZSIZE MUST BE >= NZGLBL/NZPROC
!     WITH AN EXTRA ALLOWANCE FOR ANY REMAINDER

    

!     GLOBAL GRID SIZE
real(kind=8) :: start_time, finish_time, total_time
real(kind=8) :: start_comm_time, finish_comm_time, total_comm_time

INTEGER :: nxglbl,nyglbl,nzglbl
PARAMETER(nxglbl=384, nyglbl=384, nzglbl=384)
INTEGER :: ngzmax
!     SET NGZMAX=MAX(NXGLBL,NYGLBL,NZGLBL)
PARAMETER(ngzmax=nxglbl)

!     NUMBER OF PROCESSORS
INTEGER :: nxproc,nyproc,nzproc
PARAMETER(nxproc=8, nyproc=4, nzproc=4)
INTEGER :: nprmax
!     SET NPRMAX=MAX(NXPROC,NYPROC,NZPROC)
PARAMETER(nprmax=nxproc)

!     LOCAL GRID SIZE
INTEGER :: nxsize,nysize,nzsize
PARAMETER(nxsize=48, nysize=96, nzsize=96)
INTEGER :: nszmax
!     SET NSZMAX=MAX(NXSIZE,NYSIZE,NZSIZE)
PARAMETER(nszmax=nxsize)

!     SIZE OF HALO
INTEGER :: nhalox,nhaloy,nhaloz
PARAMETER(nhalox=5,nhaloy=5,nhaloz=5)

!     SIZE OF PARALLEL TRANSFER ARRAY
!     NPARAY MUST BE >= MAX(NHALOX,NHALOY,NHALOZ)
!                            *MAX((NXSIZE+2*NHALOX)*(NYSIZE+2*NHALOY),
!                                 (NXSIZE+2*NHALOX)*(NZSIZE+2*NHALOZ),
!                                 (NYSIZE+2*NHALOY)*(NZSIZE+2*NHALOZ))
!     AND ALSO LARGE ENOUGH FOR PARALLEL BROADCAST OF CHEMICAL DATA
!     AND ALSO LARGE ENOUGH FOR PARALLEL TRANSFER OF INITIAL TURBULENCE DATA
INTEGER :: nparay
PARAMETER(nparay=8000000)

!     LOCAL BIG ARRAY SIZE
INTEGER :: nxbigl,nxbigr,nybigl,nybigr,nzbigl,nzbigr
PARAMETER(nxbigl=1-nhalox, nxbigr=nxsize+nhalox,  &
    nybigl=1-nhaloy, nybigr=nysize+nhaloy, nzbigl=1-nhaloz, nzbigr=nzsize+nhaloz)

!     PHYDIM-------------------------------------------------------------------
!     IFTURB-------------------------------------------------------------------

!     DATA FOR INITIAL TURBULENCE FIELD

!     NUMBER OF PENCILS
INTEGER :: npenmx
PARAMETER(npenmx=1)

REAL(kind=8) :: fftrow(2*ngzmax,3,npenmx), ftpart(2*nszmax,3,npenmx),  &
    fftinx(2*ngzmax)

COMMON/ifturb/fftrow,ftpart,fftinx

!     IFTURB-------------------------------------------------------------------
!     CHEMIC-------------------------------------------------------------------

!     PARAMETERS
!     ==========
!     MAX NO OF SPECIES, NO OF STEPS
INTEGER :: nspcmx,nstpmx
PARAMETER(nspcmx=2, nstpmx=1)

!     THERMODYNAMIC DATA
!     MAX NO OF TEMPERATURE INTERVALS, THERMO POLYNOMIAL COEFFICIENTS
INTEGER :: ntinmx,ncofmx
PARAMETER(ntinmx=2, ncofmx=7)
!     MAX NO OF TEMPERATURE COEFFICIENTS, DITTO MINUS ONE
INTEGER :: nctmax,nctmm1
PARAMETER(nctmax=5, nctmm1=nctmax-1)

!     TEMPERATURE INTERVAL INDEXING
!     NTBASE = NUMBER BASE FOR INDEXING:
!     MUST BE A POWER OF TWO >= MAX NO OF TEMPERATURE INTERVALS PER SPECIES
!     NSPIMX = MAX NO OF SPECIES STORED PER SINGLE (32-BIT) SIGNED INTEGER:
!     MUST BE SET EQUAL TO 31 DIV LOG2(NTBASE)
!     NINTMX = NO OF INTEGERS REQUIRED PER SPATIAL POINT:
!     MUST BE SET EQUAL TO (1 + NSPCMX DIV NSPIMX)
INTEGER :: nspimx,ntbase,nintmx
PARAMETER(nspimx=15, ntbase=4, nintmx=2)

!     CHEMICAL RATE DATA
!     MAX NO OF THIRD BODIES
INTEGER :: nbdymx
PARAMETER(nbdymx=10)
!     MAX SIZE OF STEP SPECIES-LIST, STEP REACTANT-LIST
INTEGER :: nssmax,nrsmax
PARAMETER(nssmax=10, nrsmax=10)
!     MAX NO OF LINDEMANN STEPS
INTEGER :: nllmax
PARAMETER(nllmax=10)

!     TRANSPORT COEFFICIENTS
REAL(kind=8) :: alamdc,rlamda,tlamda
PARAMETER(alamdc=0.0000258_8, rlamda=0.70_8, tlamda=298.0_8)
REAL(kind=8) :: prantl
PARAMETER(prantl=0.70_8)

!     UNIVERSAL GAS CONSTANT
REAL(kind=8) :: rguniv
PARAMETER(rguniv=8314.20_8)


!     CHEMICAL DATA
!     =============
REAL(kind=8) :: amolcp(ncofmx,ntinmx,nspcmx)
REAL(kind=8) :: amolgb(ncofmx,ntinmx,nspcmx)
REAL(kind=8) :: amascp(ncofmx,ntinmx,nspcmx)
REAL(kind=8) :: amasct(ncofmx,ntinmx,nspcmx)
REAL(kind=8) :: amasch(ncofmx,ntinmx,nspcmx)
REAL(kind=8) :: amascs(ncofmx,ntinmx,nspcmx)
REAL(kind=8) :: tintlo(ntinmx,nspcmx),tinthi(ntinmx,nspcmx)
REAL(kind=8) :: diffmu(nssmax,nstpmx),diffmw(nssmax,nstpmx)
REAL(kind=8) :: crspec(nssmax,nstpmx),cpspec(nssmax,nstpmx)
REAL(kind=8) :: effy3b(nspcmx,nbdymx)
REAL(kind=8) :: rparam(3,nstpmx)
REAL(kind=8) :: rclind(4,nllmax)
REAL(kind=8) :: rctroe(12,nllmax)
REAL(kind=8) :: rcsrif(8,nllmax)
REAL(kind=8) :: wmolar(nspcmx),ovwmol(nspcmx),rgspec(nspcmx)
REAL(kind=8) :: clewis(nspcmx),olewis(nspcmx)
REAL(kind=8) :: prefgb,alamda

INTEGER :: nsspec(nssmax,nstpmx)
INTEGER :: nrspec(nssmax,nstpmx),npspec(nssmax,nstpmx)
INTEGER :: nrcpec(nssmax,nstpmx),npcpec(nssmax,nstpmx)
INTEGER :: ncofcp(ntinmx,nspcmx)
INTEGER :: ncpoly(ntinmx,nspcmx),ncpom1(ntinmx,nspcmx)
INTEGER :: ncenth(ntinmx,nspcmx),ncenpy(ntinmx,nspcmx)
INTEGER :: nsslen(nstpmx)
INTEGER :: nrslen(nstpmx),npslen(nstpmx)
INTEGER :: nrclen(nstpmx),npclen(nstpmx)
INTEGER :: mblist(nstpmx)
INTEGER :: mglist(nstpmx)
INTEGER :: mllist(nstpmx)
INTEGER :: mtlist(nstpmx)
INTEGER :: mslist(nstpmx)
INTEGER :: ntint(nspcmx)
INTEGER :: nspec,nspm1,nstep,nbody,ngibb,nlind,ntroe,nsrif

!     MAX LENGTH OF SPECIES SYMBOL STRINGS
INTEGER :: nspstr
PARAMETER(nspstr=10)
!     SPECIES AND THIRD BODY STRINGS MUST BE DECLARED AS CHARACTER*NSPSTR
CHARACTER (LEN=10) :: spcsym(nspcmx)
CHARACTER (LEN=10) :: bdysym(nspcmx)

COMMON/chemic/amolcp,amolgb,amascp,amasct,amasch,amascs, tintlo,tinthi,  &
    diffmu,diffmw,crspec,cpspec,effy3b, rparam,rclind,rctroe,rcsrif,  &
    wmolar,ovwmol,rgspec,clewis,olewis, prefgb,alamda,  &
    nsspec,nrspec,npspec,nrcpec,npcpec, ncofcp,ncpoly,ncpom1,ncenth,ncenpy,  &
    nsslen,nrslen,npslen,nrclen,npclen, mblist,mglist,mllist,mtlist,mslist,  &
    ntint, nspec,nspm1,nstep,nbody,ngibb,nlind,ntroe,nsrif,  &
    spcsym,bdysym

!     CHEMIC-------------------------------------------------------------------
!     DIFFUS-------------------------------------------------------------------

!     MIXTURE AVERAGED TRANSPORT

!     MASS FRACTION TOLERANCE
REAL(kind=8) :: dfctol
PARAMETER(dfctol = 0.0000000000010_8)

!     MAX NUMBERS OF POLYNOMIAL COEFFICIENTS
INTEGER :: ndcfmx,nvcfmx,nccfmx
PARAMETER(ndcfmx = 4, nvcfmx = 4, nccfmx = 4)

!     MOLECULAR TRANSPORT DATA
REAL(kind=8) :: diffco(ndcfmx,nspcmx,nspcmx)
REAL(kind=8) :: tdrcco(ndcfmx,nspcmx,nspcmx)
REAL(kind=8) :: wilko1(nspcmx,nspcmx),wilko2(nspcmx,nspcmx)
REAL(kind=8) :: viscco(nvcfmx,nspcmx)
REAL(kind=8) :: condco(nccfmx,nspcmx)
REAL(kind=8) :: pdifgb,tdifgb
REAL(kind=8) :: wmltdr

INTEGER :: ncodif,ncotdr,ncovis,ncocon
INTEGER :: ncodm1,ncotm1,ncovm1,ncocm1
INTEGER :: nfmavt,nfmixw,nfmixp,nfmsor,nfmduf

!     CONTROL FLAGS
LOGICAL :: flmsor(nspcmx),flmduf(nspcmx),flmtdr(nspcmx)
LOGICAL :: flmavt,flmixw,flmixp,flmixt

COMMON/diffus/diffco,tdrcco,wilko1,wilko2,viscco,condco,  &
    pdifgb,tdifgb,wmltdr, ncodif,ncotdr,ncovis,ncocon,  &
    ncodm1,ncotm1,ncovm1,ncocm1, nfmavt,nfmixw,nfmixp,nfmsor,nfmduf,  &
    flmsor,flmduf,flmtdr,flmavt,flmixw,flmixp,flmixt

!     DIFFUS-------------------------------------------------------------------
!     RADIAT-------------------------------------------------------------------

!     RADIATION TREATMENT

REAL(kind=8) :: stefbo
PARAMETER(stefbo = 0.00000005670373210_8)

INTEGER :: ncfrmx
PARAMETER(ncfrmx = 6)

REAL(kind=8) :: akprad(ncfrmx,nspcmx)
REAL(kind=8) :: foursb,trefrn,trfrth

INTEGER :: nsprid(nspcmx),nkprad(nspcmx),nkprm1(nspcmx)
INTEGER :: nsprad,nfradn

LOGICAL :: flradn

COMMON/radiat/akprad,foursb,trefrn,trfrth, nsprid,nkprad,nkprm1,nsprad,  &
    nfradn, flradn

!     RADIAT-------------------------------------------------------------------
!     PARAMS-------------------------------------------------------------------

!     PHYSICAL AND NUMERICAL PARAMETERS

!     MASS FRACTION TOLERANCE
REAL(kind=8) :: ytoler
PARAMETER(ytoler = 1.0E-10)

!     NUMBERS
REAL(kind=8) :: zero,one,two,three,four,eight,  &
    half,thrf,thrd,tthd,fthd,qrtr
PARAMETER(zero=0.0_8, one=1.0_8, two=2.0_8, three=3.0_8,  &
    four=4.0_8, eight=8.0_8)
PARAMETER(half=one/two,  thrf=three*half, thrd=one/three,  &
    tthd=two*thrd, fthd=two*tthd,   qrtr=one/four)

!     LOCATION FOR REFERENCE PRESSURE VALUE
INTEGER :: ipref,jpref,kpref
PARAMETER(ipref=1,jpref=1,kpref=1)

!     PHYSICAL AND NUMERICAL DATA
REAL(kind=8) :: yrin(nspcmx)
REAL(kind=8) :: prin,trin,drin,urin,vrin,wrin,erin
REAL(kind=8) :: xgdlen,ygdlen,zgdlen
REAL(kind=8) :: etime,tstep,deltax,deltay,deltaz,pi,clnten
INTEGER :: itime,ntime,ntime1,ntime2,nstpsw,nsaved, inturb,inflam,inseed,  &
    nxgreq,nygreq,nzgreq,nxpreq,nypreq,nzpreq,nspreq,  &
    ntdump,ntrept,niters,itstat,ntstat,idflag,  &
    ncont1,ncont2,ncont3,ncont4,ncont5,ncont6,  &
    nccont,ncchem,ncdiff,ncradn,ncrept,ncstat,ncdmpi,ncdmpo, ndifmt,ndofmt
CHARACTER (LEN=60) :: fncont,fnchem,fndiff,fnradn,fnrept,fnstat,fndmpo(2)
LOGICAL :: fladpt

COMMON/params/ yrin,prin,trin,drin,urin,vrin,wrin,erin,  &
    xgdlen,ygdlen,zgdlen, etime,tstep,deltax,deltay,deltaz,pi,clnten,  &
    itime,ntime,ntime1,ntime2,nstpsw,nsaved, inturb,inflam,inseed,  &
    nxgreq,nygreq,nzgreq,nxpreq,nypreq,nzpreq,nspreq,  &
    ntdump,ntrept,niters,itstat,ntstat,idflag,  &
    ncont1,ncont2,ncont3,ncont4,ncont5,ncont6,  &
    nccont,ncchem,ncdiff,ncradn,ncrept,ncstat,ncdmpi,ncdmpo, ndifmt,ndofmt,  &
    fncont,fnchem,fndiff,fnradn,fnrept,fnstat,fndmpo,fladpt

!     PARAMS-------------------------------------------------------------------
!     DFDIFF-------------------------------------------------------------------

!     COEFFICIENTS AND END CONDITIONS FOR SPATIAL DIFFERENCING SCHEMES

!     DERIVATIVE BC STATUS
INTEGER :: nobc,nbound,nperi
PARAMETER(nobc=0,nbound=1,nperi=2)

!     COEFFICIENTS OF SPATIAL DIFFERENCING SCHEMES
REAL(kind=8) :: acoeff,bcoeff,ccoeff,dcoeff,ecoeff,  &
    acoffx,bcoffx,ccoffx,dcoffx,ecoffx, acoffy,bcoffy,ccoffy,dcoffy,ecoffy,  &
    acoffz,bcoffz,ccoffz,dcoffz,ecoffz, acoef1,bcoef1,ccoef1,dcoef1,  &
    acof1x,bcof1x,ccof1x,dcof1x, acof1y,bcof1y,ccof1y,dcof1y,  &
    acof1z,bcof1z,ccof1z,dcof1z, acoef2,bcoef2,ccoef2,dcoef2,  &
    acof2x,bcof2x,ccof2x,dcof2x, acof2y,bcof2y,ccof2y,dcof2y,  &
    acof2z,bcof2z,ccof2z,dcof2z, acoef3,bcoef3,  &
    acof3x,bcof3x,acof3y,bcof3y,acof3z,bcof3z, acoef4,bcoef4,ccoef4,  &
    acof4x,bcof4x,ccof4x, acof4y,bcof4y,ccof4y,  &
    acof4z,bcof4z,ccof4z, acoef5,bcoef5,ccoef5,dcoef5,  &
    acof5x,bcof5x,ccof5x,dcof5x, acof5y,bcof5y,ccof5y,dcof5y,  &
    acof5z,bcof5z,ccof5z,dcof5z, acoefs,bcoefs,ccoefs,dcoefs,ecoefs,  &
    acofsx,bcofsx,ccofsx,dcofsx,ecofsx, acofsy,bcofsy,ccofsy,dcofsy,ecofsy,  &
    acofsz,bcofsz,ccofsz,dcofsz,ecofsz, acofs1,bcofs1,ccofs1,dcofs1,ecofs1,  &
    acfs1x,bcfs1x,ccfs1x,dcfs1x,ecfs1x, acfs1y,bcfs1y,ccfs1y,dcfs1y,ecfs1y,  &
    acfs1z,bcfs1z,ccfs1z,dcfs1z,ecfs1z, acofs2,bcofs2,ccofs2,dcofs2,ecofs2,  &
    acfs2x,bcfs2x,ccfs2x,dcfs2x,ecfs2x, acfs2y,bcfs2y,ccfs2y,dcfs2y,ecfs2y,  &
    acfs2z,bcfs2z,ccfs2z,dcfs2z,ecfs2z, acofs3,bcofs3,  &
    acfs3x,bcfs3x,acfs3y,bcfs3y,acfs3z,bcfs3z, acofs4,bcofs4,ccofs4,  &
    acfs4x,bcfs4x,ccfs4x, acfs4y,bcfs4y,ccfs4y,  &
    acfs4z,bcfs4z,ccfs4z, acofs5,bcofs5,ccofs5,dcofs5,  &
    acfs5x,bcfs5x,ccfs5x,dcfs5x, acfs5y,bcfs5y,ccfs5y,dcfs5y,  &
    acfs5z,bcfs5z,ccfs5z,dcfs5z, acoefx,bcoefx,ccoefx,dcoefx,ecoefx,  &
    acofxy,bcofxy,ccofxy,dcofxy,ecofxy, acofxz,bcofxz,ccofxz,dcofxz,ecofxz,  &
    acofyz,bcofyz,ccofyz,dcofyz,ecofyz,  &
    acofx1,bcofx1,acofy1,bcofy1,acofz1,bcofz1, acf1xy,bcf1xy,ccf1xy,dcf1xy,  &
    acf1xz,bcf1xz,ccf1xz,dcf1xz, acf1yz,bcf1yz,ccf1yz,dcf1yz,  &
    acf2xy,bcf2xy,ccf2xy,dcf2xy, acf2xz,bcf2xz,ccf2xz,dcf2xz,  &
    acf2yz,bcf2yz,ccf2yz,dcf2yz, acf3xy,bcf3xy,  &
    acf3xz,bcf3xz, acf3yz,bcf3yz,  &
    acf4xy,bcf4xy,ccf4xy, acf4xz,bcf4xz,ccf4xz,  &
    acf4yz,bcf4yz,ccf4yz, acf5xy,bcf5xy,ccf5xy,dcf5xy,  &
    acf5xz,bcf5xz,ccf5xz,dcf5xz, acf5yz,bcf5yz,ccf5yz,dcf5yz,  &
    acofc1,bcofc1,ccofc1,dcofc1, acofc2,bcofc2,ccofc2,dcofc2,  &
    acc1xy,bcc1xy,ccc1xy,dcc1xy, acc1xz,bcc1xz,ccc1xz,dcc1xz,  &
    acc1yz,bcc1yz,ccc1yz,dcc1yz, acc2xy,bcc2xy,ccc2xy,dcc2xy,  &
    acc2xz,bcc2xz,ccc2xz,dcc2xz, acc2yz,bcc2yz,ccc2yz,dcc2yz,  &
    ovdelx,ovdely,ovdelz,ovdlx2,ovdly2,ovdlz2

!     SPATIAL DERIVATIVE END CONDITIONS
INTEGER :: nendxl,nendxr,nendyl,nendyr,nendzl,nendzr

COMMON/dfdiff/acoeff,bcoeff,ccoeff,dcoeff,ecoeff,  &
    acoffx,bcoffx,ccoffx,dcoffx,ecoffx, acoffy,bcoffy,ccoffy,dcoffy,ecoffy,  &
    acoffz,bcoffz,ccoffz,dcoffz,ecoffz, acoef1,bcoef1,ccoef1,dcoef1,  &
    acof1x,bcof1x,ccof1x,dcof1x, acof1y,bcof1y,ccof1y,dcof1y,  &
    acof1z,bcof1z,ccof1z,dcof1z, acoef2,bcoef2,ccoef2,dcoef2,  &
    acof2x,bcof2x,ccof2x,dcof2x, acof2y,bcof2y,ccof2y,dcof2y,  &
    acof2z,bcof2z,ccof2z,dcof2z, acoef3,bcoef3,  &
    acof3x,bcof3x,acof3y,bcof3y,acof3z,bcof3z, acoef4,bcoef4,ccoef4,  &
    acof4x,bcof4x,ccof4x, acof4y,bcof4y,ccof4y,  &
    acof4z,bcof4z,ccof4z, acoef5,bcoef5,ccoef5,dcoef5,  &
    acof5x,bcof5x,ccof5x,dcof5x, acof5y,bcof5y,ccof5y,dcof5y,  &
    acof5z,bcof5z,ccof5z,dcof5z, acoefs,bcoefs,ccoefs,dcoefs,ecoefs,  &
    acofsx,bcofsx,ccofsx,dcofsx,ecofsx, acofsy,bcofsy,ccofsy,dcofsy,ecofsy,  &
    acofsz,bcofsz,ccofsz,dcofsz,ecofsz, acofs1,bcofs1,ccofs1,dcofs1,ecofs1,  &
    acfs1x,bcfs1x,ccfs1x,dcfs1x,ecfs1x, acfs1y,bcfs1y,ccfs1y,dcfs1y,ecfs1y,  &
    acfs1z,bcfs1z,ccfs1z,dcfs1z,ecfs1z, acofs2,bcofs2,ccofs2,dcofs2,ecofs2,  &
    acfs2x,bcfs2x,ccfs2x,dcfs2x,ecfs2x, acfs2y,bcfs2y,ccfs2y,dcfs2y,ecfs2y,  &
    acfs2z,bcfs2z,ccfs2z,dcfs2z,ecfs2z, acofs3,bcofs3,  &
    acfs3x,bcfs3x,acfs3y,bcfs3y,acfs3z,bcfs3z, acofs4,bcofs4,ccofs4,  &
    acfs4x,bcfs4x,ccfs4x, acfs4y,bcfs4y,ccfs4y,  &
    acfs4z,bcfs4z,ccfs4z, acofs5,bcofs5,ccofs5,dcofs5,  &
    acfs5x,bcfs5x,ccfs5x,dcfs5x, acfs5y,bcfs5y,ccfs5y,dcfs5y,  &
    acfs5z,bcfs5z,ccfs5z,dcfs5z, acoefx,bcoefx,ccoefx,dcoefx,ecoefx,  &
    acofxy,bcofxy,ccofxy,dcofxy,ecofxy, acofxz,bcofxz,ccofxz,dcofxz,ecofxz,  &
    acofyz,bcofyz,ccofyz,dcofyz,ecofyz,  &
    acofx1,bcofx1,acofy1,bcofy1,acofz1,bcofz1, acf1xy,bcf1xy,ccf1xy,dcf1xy,  &
    acf1xz,bcf1xz,ccf1xz,dcf1xz, acf1yz,bcf1yz,ccf1yz,dcf1yz,  &
    acf2xy,bcf2xy,ccf2xy,dcf2xy, acf2xz,bcf2xz,ccf2xz,dcf2xz,  &
    acf2yz,bcf2yz,ccf2yz,dcf2yz, acf3xy,bcf3xy,  &
    acf3xz,bcf3xz, acf3yz,bcf3yz,  &
    acf4xy,bcf4xy,ccf4xy, acf4xz,bcf4xz,ccf4xz,  &
    acf4yz,bcf4yz,ccf4yz, acf5xy,bcf5xy,ccf5xy,dcf5xy,  &
    acf5xz,bcf5xz,ccf5xz,dcf5xz, acf5yz,bcf5yz,ccf5yz,dcf5yz,  &
    acofc1,bcofc1,ccofc1,dcofc1, acofc2,bcofc2,ccofc2,dcofc2,  &
    acc1xy,bcc1xy,ccc1xy,dcc1xy, acc1xz,bcc1xz,ccc1xz,dcc1xz,  &
    acc1yz,bcc1yz,ccc1yz,dcc1yz, acc2xy,bcc2xy,ccc2xy,dcc2xy,  &
    acc2xz,bcc2xz,ccc2xz,dcc2xz, acc2yz,bcc2yz,ccc2yz,dcc2yz,  &
    ovdelx,ovdely,ovdelz,ovdlx2,ovdly2,ovdlz2,  &
    nendxl,nendxr,nendyl,nendyr,nendzl,nendzr

!     DFDIFF-------------------------------------------------------------------
!     EMSTRT-------------------------------------------------------------------

!     DATA FOR MESH STRETCHING
REAL(kind=8) :: gcmreg(nszmax),gcmstr(nszmax),  &
    dgdhat(nszmax),dgdhsq(nszmax),d2gdh2(nszmax)

COMMON/emstrt/gcmreg,gcmstr,dgdhat,dgdhsq,d2gdh2

!     EMSTRT-------------------------------------------------------------------
!     FILCOM-------------------------------------------------------------------

!     SPATIAL FILTER COEFFICIENTS
REAL(kind=8) :: facofx,fbcofx,fccofx,fdcofx,fecofx,ffcofx,fgcofx,  &
    facofy,fbcofy,fccofy,fdcofy,fecofy,ffcofy,fgcofy,  &
    facofz,fbcofz,fccofz,fdcofz,fecofz,ffcofz,fgcofz,  &
    facf1x,fbcf1x,fccf1x,fdcf1x,fecf1x,ffcf1x,fgcf1x,  &
    facf1y,fbcf1y,fccf1y,fdcf1y,fecf1y,ffcf1y,fgcf1y,  &
    facf1z,fbcf1z,fccf1z,fdcf1z,fecf1z,ffcf1z,fgcf1z,  &
    facf2x,fbcf2x,fccf2x,fdcf2x,fecf2x,ffcf2x,fgcf2x, fhcf2x,  &
    facf2y,fbcf2y,fccf2y,fdcf2y,fecf2y,ffcf2y,fgcf2y, fhcf2y,  &
    facf2z,fbcf2z,fccf2z,fdcf2z,fecf2z,ffcf2z,fgcf2z, fhcf2z,  &
    facf3x,fbcf3x,fccf3x,fdcf3x,fecf3x,ffcf3x,fgcf3x, fhcf3x,ficf3x,  &
    facf3y,fbcf3y,fccf3y,fdcf3y,fecf3y,ffcf3y,fgcf3y, fhcf3y,ficf3y,  &
    facf3z,fbcf3z,fccf3z,fdcf3z,fecf3z,ffcf3z,fgcf3z, fhcf3z,ficf3z,  &
    facf4x,fbcf4x,fccf4x,fdcf4x,fecf4x,ffcf4x,fgcf4x, fhcf4x,ficf4x,fjcf4x,  &
    facf4y,fbcf4y,fccf4y,fdcf4y,fecf4y,ffcf4y,fgcf4y, fhcf4y,ficf4y,fjcf4y,  &
    facf4z,fbcf4z,fccf4z,fdcf4z,fecf4z,ffcf4z,fgcf4z, fhcf4z,ficf4z,fjcf4z,  &
    facf5x,fbcf5x,fccf5x,fdcf5x,fecf5x,ffcf5x,fgcf5x,  &
    fhcf5x,ficf5x,fjcf5x,fkcf5x,  &
    facf5y,fbcf5y,fccf5y,fdcf5y,fecf5y,ffcf5y,fgcf5y,  &
    fhcf5y,ficf5y,fjcf5y,fkcf5y,  &
    facf5z,fbcf5z,fccf5z,fdcf5z,fecf5z,ffcf5z,fgcf5z,  &
    fhcf5z,ficf5z,fjcf5z,fkcf5z,  &
    facf6x,fbcf6x,fccf6x,fdcf6x,fecf6x,ffcf6x,fgcf6x,  &
    fhcf6x,ficf6x,fjcf6x,fkcf6x,flcf6x,  &
    facf6y,fbcf6y,fccf6y,fdcf6y,fecf6y,ffcf6y,fgcf6y,  &
    fhcf6y,ficf6y,fjcf6y,fkcf6y,flcf6y,  &
    facf6z,fbcf6z,fccf6z,fdcf6z,fecf6z,ffcf6z,fgcf6z,  &
    fhcf6z,ficf6z,fjcf6z,fkcf6z,flcf6z

COMMON/filcom/facofx,fbcofx,fccofx,fdcofx,fecofx,ffcofx,fgcofx,  &
    facofy,fbcofy,fccofy,fdcofy,fecofy,ffcofy,fgcofy,  &
    facofz,fbcofz,fccofz,fdcofz,fecofz,ffcofz,fgcofz,  &
    facf1x,fbcf1x,fccf1x,fdcf1x,fecf1x,ffcf1x,fgcf1x,  &
    facf1y,fbcf1y,fccf1y,fdcf1y,fecf1y,ffcf1y,fgcf1y,  &
    facf1z,fbcf1z,fccf1z,fdcf1z,fecf1z,ffcf1z,fgcf1z,  &
    facf2x,fbcf2x,fccf2x,fdcf2x,fecf2x,ffcf2x,fgcf2x, fhcf2x,  &
    facf2y,fbcf2y,fccf2y,fdcf2y,fecf2y,ffcf2y,fgcf2y, fhcf2y,  &
    facf2z,fbcf2z,fccf2z,fdcf2z,fecf2z,ffcf2z,fgcf2z, fhcf2z,  &
    facf3x,fbcf3x,fccf3x,fdcf3x,fecf3x,ffcf3x,fgcf3x, fhcf3x,ficf3x,  &
    facf3y,fbcf3y,fccf3y,fdcf3y,fecf3y,ffcf3y,fgcf3y, fhcf3y,ficf3y,  &
    facf3z,fbcf3z,fccf3z,fdcf3z,fecf3z,ffcf3z,fgcf3z, fhcf3z,ficf3z,  &
    facf4x,fbcf4x,fccf4x,fdcf4x,fecf4x,ffcf4x,fgcf4x, fhcf4x,ficf4x,fjcf4x,  &
    facf4y,fbcf4y,fccf4y,fdcf4y,fecf4y,ffcf4y,fgcf4y, fhcf4y,ficf4y,fjcf4y,  &
    facf4z,fbcf4z,fccf4z,fdcf4z,fecf4z,ffcf4z,fgcf4z, fhcf4z,ficf4z,fjcf4z,  &
    facf5x,fbcf5x,fccf5x,fdcf5x,fecf5x,ffcf5x,fgcf5x,  &
    fhcf5x,ficf5x,fjcf5x,fkcf5x,  &
    facf5y,fbcf5y,fccf5y,fdcf5y,fecf5y,ffcf5y,fgcf5y,  &
    fhcf5y,ficf5y,fjcf5y,fkcf5y,  &
    facf5z,fbcf5z,fccf5z,fdcf5z,fecf5z,ffcf5z,fgcf5z,  &
    fhcf5z,ficf5z,fjcf5z,fkcf5z,  &
    facf6x,fbcf6x,fccf6x,fdcf6x,fecf6x,ffcf6x,fgcf6x,  &
    fhcf6x,ficf6x,fjcf6x,fkcf6x,flcf6x,  &
    facf6y,fbcf6y,fccf6y,fdcf6y,fecf6y,ffcf6y,fgcf6y,  &
    fhcf6y,ficf6y,fjcf6y,fkcf6y,flcf6y,  &
    facf6z,fbcf6z,fccf6z,fdcf6z,fecf6z,ffcf6z,fgcf6z,  &
    fhcf6z,ficf6z,fjcf6z,fkcf6z,flcf6z

!     FILCOM-------------------------------------------------------------------
!     RUNKUT-------------------------------------------------------------------

!     RUNGE-KUTTA PARAMETERS
INTEGER :: nrkmax
PARAMETER(nrkmax=5)

REAL(kind=8) :: rklhs(nrkmax),rkrhs(nrkmax)
REAL(kind=8) :: rkerr(nrkmax),rktim(nrkmax)
REAL(kind=8) :: ctmult,ctalph,ctbeta,ctgama
REAL(kind=8) :: errtol,errlow,trmin,trmax,tsmin,tsmax
REAL(kind=8) :: errold,errldr,btime
REAL(kind=8) :: erdnrm,erunrm,ervnrm,erwnrm,erenrm
REAL(kind=8) :: erynrm(nspcmx)
INTEGER :: irkstp,nrkstp,nrksm1
INTEGER :: inderr
LOGICAL :: fupelc

COMMON/runkut/rklhs,rkrhs,rkerr,rktim, ctmult,ctalph,ctbeta,ctgama,  &
    errtol,errlow,trmin,trmax,tsmin,tsmax, errold,errldr,btime,  &
    erdnrm,erunrm,ervnrm,erwnrm,erenrm,erynrm, irkstp,nrkstp,nrksm1,inderr,fupelc

!     RUNKUT-------------------------------------------------------------------
!     RKNORM-------------------------------------------------------------------

!     RUNGE-KUTTA SUBSTEP ERROR NORMS
REAL(kind=8) :: erdrhs(nrkmax),erurhs(nrkmax),ervrhs(nrkmax),  &
    erwrhs(nrkmax),ererhs(nrkmax)
REAL(kind=8) :: eryrhs(nspcmx,nrkmax)

COMMON/rknorm/erdrhs,erurhs,ervrhs,erwrhs,ererhs,eryrhs

!     RKNORM-------------------------------------------------------------------
!     PRVARS-------------------------------------------------------------------

!     PRINCIPAL VARIABLES (ALL STANDARD SIZE ARRAYS)
REAL(kind=8) :: drun(nxsize,nysize,nzsize),urun(nxsize,nysize,nzsize),  &
    vrun(nxsize,nysize,nzsize),wrun(nxsize,nysize,nzsize),  &
    erun(nxsize,nysize,nzsize), yrun(nxsize,nysize,nzsize,nspcmx)

COMMON/prvars/drun,urun,vrun,wrun,erun,yrun

!     PRVARS-------------------------------------------------------------------
!     TIMRHS-------------------------------------------------------------------

!     TIME-STEPPING RHS TERMS (ALL BIGGER SIZE ARRAYS)
REAL(kind=8) :: drhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    urhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    vrhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    wrhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    erhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    yrhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr,nspcmx)

COMMON/timrhs/drhs,urhs,vrhs,wrhs,erhs,yrhs

!     TIMRHS-------------------------------------------------------------------
!     TIMERR-------------------------------------------------------------------

!     TIME-STEPPING ERROR ARRAYS (ALL STANDARD SIZE ARRAYS)
REAL(kind=8) :: derr(nxsize,nysize,nzsize),uerr(nxsize,nysize,nzsize),  &
    verr(nxsize,nysize,nzsize),werr(nxsize,nysize,nzsize),  &
    eerr(nxsize,nysize,nzsize), yerr(nxsize,nysize,nzsize,nspcmx)

COMMON/timerr/derr,uerr,verr,werr,eerr,yerr

!     TIMERR-------------------------------------------------------------------
!     WORKSP-------------------------------------------------------------------

!     WORKSPACE (STANDARD SIZE ARRAYS)
REAL(kind=8) ::  &
    store1(nxsize,nysize,nzsize),store2(nxsize,nysize,nzsize),  &
    store3(nxsize,nysize,nzsize),store4(nxsize,nysize,nzsize),  &
    store5(nxsize,nysize,nzsize),store6(nxsize,nysize,nzsize),  &
    divm(nxsize,nysize,nzsize), rate(nxsize,nysize,nzsize,nspcmx),  &
    rrte(nxsize,nysize,nzsize,nspcmx)

COMMON/worksp/store1,store2,store3,store4,store5,store6,divm,rate, rrte

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     WORKSPACE (DIFFUSIVE CORRECTION VELOCITY, STANDARD SIZE ARRAYS)
REAL(kind=8) :: ucor(nxsize,nysize,nzsize),vcor(nxsize,nysize,nzsize),  &
    wcor(nxsize,nysize,nzsize)

COMMON/worksc/ucor,vcor,wcor
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     WORKSPACE (MIXTURE AVERAGED TRANSPORT, STANDARD SIZE ARRAYS)
REAL(kind=8) :: wd1x(nxsize,nysize,nzsize),wd1y(nxsize,nysize,nzsize),  &
    wd1z(nxsize,nysize,nzsize),  &
    wd2x(nxsize,nysize,nzsize),wd2y(nxsize,nysize,nzsize),  &
    wd2z(nxsize,nysize,nzsize),  &
    pd1x(nxsize,nysize,nzsize),pd1y(nxsize,nysize,nzsize),  &
    pd1z(nxsize,nysize,nzsize),  &
    pd2x(nxsize,nysize,nzsize),pd2y(nxsize,nysize,nzsize),  &
    pd2z(nxsize,nysize,nzsize),  &
    td1x(nxsize,nysize,nzsize),td1y(nxsize,nysize,nzsize),  &
    td1z(nxsize,nysize,nzsize),  &
    td2x(nxsize,nysize,nzsize),td2y(nxsize,nysize,nzsize),  &
    td2z(nxsize,nysize,nzsize)

COMMON/workst/wd1x,wd1y,wd1z,wd2x,wd2y,wd2z, pd1x,pd1y,pd1z,pd2x,pd2y,pd2z,  &
    td1x,td1y,td1z,td2x,td2y,td2z

!     WORKSPACE (MIXTURE AVERAGED TRANSPORT, BIGGER SIZE ARRAYS)
REAL(kind=8) :: wmomix(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    difmix(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    tdrmix(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)

COMMON/workbt/wmomix,difmix,tdrmix
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     WORKSPACE (BIGGER SIZE ARRAYS)
REAL(kind=8) :: utmp(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    vtmp(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    wtmp(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    prun(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    trun(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    transp(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    store7(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)

COMMON/worksb/utmp,vtmp,wtmp,prun,trun,transp,store7

!     WORKSPACE (TEMPERATURE INDEXING)
INTEGER :: itndex(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr,nintmx)

COMMON/worksi/itndex

!     WORKSPACE (PARALLEL TRANSFER ARRAY)
REAL(kind=8) :: parray(nparay)

COMMON/workpt/parray

!     WORKSP-------------------------------------------------------------------
!     NSBCCL-------------------------------------------------------------------

!     BC TYPE PARAMETERS
INTEGER :: nsnull,nsperi, nsbci1,nsbci2,nsbci3,nsbci4,  &
    nsbco1,nsbco2,nsbco3,nsbco4, nsbcw1,nsbcw2,nsbcw3,nsbcw4
PARAMETER(nsnull=0, nsperi=1, nsbci1=11,nsbci2=12,nsbci3=13,nsbci4=14,  &
    nsbco1=21,nsbco2=22,nsbco3=23,nsbco4=24,  &
    nsbcw1=31,nsbcw2=32,nsbcw3=33,nsbcw4=34)

!     BC NUMBER OF PARAMETERS
INTEGER :: nbcpri,nbcprr
PARAMETER(nbcpri=4, nbcprr=4)

!     GLOBAL BC TYPE VARIABLES
INTEGER :: ngbcxl,ngbcxr,ngbcyl,ngbcyr,ngbczl,ngbczr
COMMON/nglbct/ngbcxl,ngbcxr,ngbcyl,ngbcyr,ngbczl,ngbczr

!     LOCAL BC TYPE VARIABLES
INTEGER :: nsbcxl,nsbcxr,nsbcyl,nsbcyr,nsbczl,nsbczr
COMMON/nsbcct/nsbcxl,nsbcxr,nsbcyl,nsbcyr,nsbczl,nsbczr

!     BC FLAGS
LOGICAL :: fxlcnv,fxlvsn,fxlvst,fxlcon,fxladb,fxldif, fxlcnw,fxladw,fxldfw,  &
    fxrcnv,fxrvsn,fxrvst,fxrcon,fxradb,fxrdif, fxrcnw,fxradw,fxrdfw,  &
    fylcnv,fylvsn,fylvst,fylcon,fyladb,fyldif, fylcnw,fyladw,fyldfw,  &
    fyrcnv,fyrvsn,fyrvst,fyrcon,fyradb,fyrdif, fyrcnw,fyradw,fyrdfw,  &
    fzlcnv,fzlvsn,fzlvst,fzlcon,fzladb,fzldif, fzlcnw,fzladw,fzldfw,  &
    fzrcnv,fzrvsn,fzrvst,fzrcon,fzradb,fzrdif, fzrcnw,fzradw,fzrdfw,  &
    fxltrb,fxrtrb,fyltrb,fyrtrb,fzltrb,fzrtrb
COMMON/nsbccf/fxlcnv,fxlvsn,fxlvst,fxlcon,fxladb,fxldif,  &
    fxlcnw,fxladw,fxldfw, fxrcnv,fxrvsn,fxrvst,fxrcon,fxradb,fxrdif,  &
    fxrcnw,fxradw,fxrdfw, fylcnv,fylvsn,fylvst,fylcon,fyladb,fyldif,  &
    fylcnw,fyladw,fyldfw, fyrcnv,fyrvsn,fyrvst,fyrcon,fyradb,fyrdif,  &
    fyrcnw,fyradw,fyrdfw, fzlcnv,fzlvsn,fzlvst,fzlcon,fzladb,fzldif,  &
    fzlcnw,fzladw,fzldfw, fzrcnv,fzrvsn,fzrvst,fzrcon,fzradb,fzrdif,  &
    fzrcnw,fzradw,fzrdfw, fxltrb,fxrtrb,fyltrb,fyrtrb,fzltrb,fzrtrb

!     BC DATA
REAL(kind=8) :: rxlprm(nbcprr),rxrprm(nbcprr),  &
    rylprm(nbcprr),ryrprm(nbcprr), rzlprm(nbcprr),rzrprm(nbcprr)
REAL(kind=8) :: cobcxl,cobcxr,cobcyl,cobcyr,cobczl,cobczr,  &
    pinfxl,pinfxr,pinfyl,pinfyr,pinfzl,pinfzr
REAL(kind=8) :: slocxl,slocxr,slocyl,slocyr,sloczl,sloczr,  &
    elocxl,elocxr,elocyl,elocyr,eloczl,eloczr,  &
    bvelxl,bvelxr,bvelyl,bvelyr,bvelzl,bvelzr,  &
    svelxl,svelxr,svelyl,svelyr,svelzl,svelzr,  &
    scauxl,scauxr,scauyl,scauyr,scauzl,scauzr,  &
    scduxl,scduxr,scduyl,scduyr,scduzl,scduzr, tpovxg,tpovyg,tpovzg
INTEGER :: nxlprm(nbcpri),nxrprm(nbcpri), nylprm(nbcpri),nyrprm(nbcpri),  &
    nzlprm(nbcpri),nzrprm(nbcpri)
INTEGER :: istaxl,istoxl,istayl,istoyl,istazl,istozl,  &
    istaxr,istoxr,istayr,istoyr,istazr,istozr,  &
    kminxl,kminxr,kminyl,kminyr,kminzl,kminzr
INTEGER :: nctixl
CHARACTER (LEN=60) :: fntixl,fntcxl
LOGICAL :: fllixl,fllixr,flliyl,flliyr,fllizl,fllizr,  &
    fltrxl,fltrxr,fltryl,fltryr,fltrzl,fltrzr
COMMON/nsbccp/rxlprm,rxrprm,rylprm,ryrprm,rzlprm,rzrprm,  &
    cobcxl,cobcxr,cobcyl,cobcyr,cobczl,cobczr,  &
    pinfxl,pinfxr,pinfyl,pinfyr,pinfzl,pinfzr,  &
    slocxl,slocxr,slocyl,slocyr,sloczl,sloczr,  &
    elocxl,elocxr,elocyl,elocyr,eloczl,eloczr,  &
    bvelxl,bvelxr,bvelyl,bvelyr,bvelzl,bvelzr,  &
    svelxl,svelxr,svelyl,svelyr,svelzl,svelzr,  &
    scauxl,scauxr,scauyl,scauyr,scauzl,scauzr,  &
    scduxl,scduxr,scduyl,scduyr,scduzl,scduzr, tpovxg,tpovyg,tpovzg,  &
    nxlprm,nxrprm,nylprm,nyrprm,nzlprm,nzrprm,  &
    istaxl,istoxl,istayl,istoyl,istazl,istozl,  &
    istaxr,istoxr,istayr,istoyr,istazr,istozr,  &
    kminxl,kminxr,kminyl,kminyr,kminzl,kminzr, nctixl,  &
    fntixl,fntcxl, fllixl,fllixr,flliyl,flliyr,fllizl,fllizr,  &
    fltrxl,fltrxr,fltryl,fltryr,fltrzl,fltrzr

!     INFLOW VELOCITY FIELD DATA
!     REQUIRED ONLY FOR TURBULENT INFLOW
REAL(kind=8) :: ufxl(nxsize,nysize,nzsize), vfxl(nxsize,nysize,nzsize),  &
    wfxl(nxsize,nysize,nzsize)
COMMON/nbcinv/ufxl,vfxl,wfxl

!     WALL BC DIFFERENCING DATA
INTEGER :: ncbcsz
PARAMETER(ncbcsz=5)
REAL(kind=8) :: acbcxl(ncbcsz),acbcxr(ncbcsz),  &
    acbcyl(ncbcsz),acbcyr(ncbcsz), acbczl(ncbcsz),acbczr(ncbcsz)
COMMON/bcdifw/acbcxl,acbcxr,acbcyl,acbcyr,acbczl,acbczr

!     X-DIRECTION LEFT-HAND END
REAL(kind=8) ::  &
    bclyxl(nysize,nzsize,nspcmx),stryxl(nysize,nzsize,nspcmx),  &
    dydtxl(nysize,nzsize,nspcmx),ratexl(nysize,nzsize,nspcmx),  &
    strhxl(nysize,nzsize,nspcmx), bcl1xl(nysize,nzsize),bcl2xl(nysize,nzsize),  &
    bcl3xl(nysize,nzsize),bcl4xl(nysize,nzsize),  &
    bcl5xl(nysize,nzsize),bcltxl(nysize,nzsize),  &
    struxl(nysize,nzsize),strvxl(nysize,nzsize),  &
    strwxl(nysize,nzsize),strpxl(nysize,nzsize),  &
    strdxl(nysize,nzsize),strtxl(nysize,nzsize),  &
    strexl(nysize,nzsize),strgxl(nysize,nzsize), strrxl(nysize,nzsize),  &
    dudtxl(nysize,nzsize),dvdtxl(nysize,nzsize),  &
    dwdtxl(nysize,nzsize),dtdtxl(nysize,nzsize), dddtxl(nysize,nzsize),  &
    acouxl(nysize,nzsize),ova2xl(nysize,nzsize),  &
    gam1xl(nysize,nzsize),ovgmxl(nysize,nzsize),  &
    sydtxl(nysize,nzsize),sorpxl(nysize,nzsize)
COMMON/nbccxl/bclyxl,stryxl,dydtxl,ratexl,strhxl,  &
    bcl1xl,bcl2xl,bcl3xl,bcl4xl,bcl5xl,bcltxl,  &
    struxl,strvxl,strwxl,strpxl,strdxl,strtxl, strexl,strgxl,strrxl,  &
    dudtxl,dvdtxl,dwdtxl,dtdtxl,dddtxl, acouxl,ova2xl,gam1xl,ovgmxl,sydtxl,sorpxl

!     X-DIRECTION RIGHT-HAND END
REAL(kind=8) ::  &
    bclyxr(nysize,nzsize,nspcmx),stryxr(nysize,nzsize,nspcmx),  &
    dydtxr(nysize,nzsize,nspcmx),ratexr(nysize,nzsize,nspcmx),  &
    strhxr(nysize,nzsize,nspcmx), bcl1xr(nysize,nzsize),bcl2xr(nysize,nzsize),  &
    bcl3xr(nysize,nzsize),bcl4xr(nysize,nzsize),  &
    bcl5xr(nysize,nzsize),bcltxr(nysize,nzsize),  &
    struxr(nysize,nzsize),strvxr(nysize,nzsize),  &
    strwxr(nysize,nzsize),strpxr(nysize,nzsize),  &
    strdxr(nysize,nzsize),strtxr(nysize,nzsize),  &
    strexr(nysize,nzsize),strgxr(nysize,nzsize), strrxr(nysize,nzsize),  &
    dudtxr(nysize,nzsize),dvdtxr(nysize,nzsize),  &
    dwdtxr(nysize,nzsize),dtdtxr(nysize,nzsize), dddtxr(nysize,nzsize),  &
    acouxr(nysize,nzsize),ova2xr(nysize,nzsize),  &
    gam1xr(nysize,nzsize),ovgmxr(nysize,nzsize),  &
    sydtxr(nysize,nzsize),sorpxr(nysize,nzsize)
COMMON/nbccxr/bclyxr,stryxr,dydtxr,ratexr,strhxr,  &
    bcl1xr,bcl2xr,bcl3xr,bcl4xr,bcl5xr,bcltxr,  &
    struxr,strvxr,strwxr,strpxr,strdxr,strtxr, strexr,strgxr,strrxr,  &
    dudtxr,dvdtxr,dwdtxr,dtdtxr,dddtxr, acouxr,ova2xr,gam1xr,ovgmxr,sydtxr,sorpxr


!     Y-DIRECTION LEFT-HAND END
REAL(kind=8) ::  &
    bclyyl(nxsize,nzsize,nspcmx),stryyl(nxsize,nzsize,nspcmx),  &
    dydtyl(nxsize,nzsize,nspcmx),rateyl(nxsize,nzsize,nspcmx),  &
    strhyl(nxsize,nzsize,nspcmx), bcl1yl(nxsize,nzsize),bcl2yl(nxsize,nzsize),  &
    bcl3yl(nxsize,nzsize),bcl4yl(nxsize,nzsize),  &
    bcl5yl(nxsize,nzsize),bcltyl(nxsize,nzsize),  &
    struyl(nxsize,nzsize),strvyl(nxsize,nzsize),  &
    strwyl(nxsize,nzsize),strpyl(nxsize,nzsize),  &
    strdyl(nxsize,nzsize),strtyl(nxsize,nzsize),  &
    streyl(nxsize,nzsize),strgyl(nxsize,nzsize), strryl(nxsize,nzsize),  &
    dudtyl(nxsize,nzsize),dvdtyl(nxsize,nzsize),  &
    dwdtyl(nxsize,nzsize),dtdtyl(nxsize,nzsize), dddtyl(nxsize,nzsize),  &
    acouyl(nxsize,nzsize),ova2yl(nxsize,nzsize),  &
    gam1yl(nxsize,nzsize),ovgmyl(nxsize,nzsize),  &
    sydtyl(nxsize,nzsize),sorpyl(nxsize,nzsize)
COMMON/nbccyl/bclyyl,stryyl,dydtyl,rateyl,strhyl,  &
    bcl1yl,bcl2yl,bcl3yl,bcl4yl,bcl5yl,bcltyl,  &
    struyl,strvyl,strwyl,strpyl,strdyl,strtyl, streyl,strgyl,strryl,  &
    dudtyl,dvdtyl,dwdtyl,dtdtyl,dddtyl, acouyl,ova2yl,gam1yl,ovgmyl,sydtyl,sorpyl

!     Y-DIRECTION RIGHT-HAND END
REAL(kind=8) ::  &
    bclyyr(nxsize,nzsize,nspcmx),stryyr(nxsize,nzsize,nspcmx),  &
    dydtyr(nxsize,nzsize,nspcmx),rateyr(nxsize,nzsize,nspcmx),  &
    strhyr(nxsize,nzsize,nspcmx), bcl1yr(nxsize,nzsize),bcl2yr(nxsize,nzsize),  &
    bcl3yr(nxsize,nzsize),bcl4yr(nxsize,nzsize),  &
    bcl5yr(nxsize,nzsize),bcltyr(nxsize,nzsize),  &
    struyr(nxsize,nzsize),strvyr(nxsize,nzsize),  &
    strwyr(nxsize,nzsize),strpyr(nxsize,nzsize),  &
    strdyr(nxsize,nzsize),strtyr(nxsize,nzsize),  &
    streyr(nxsize,nzsize),strgyr(nxsize,nzsize), strryr(nxsize,nzsize),  &
    dudtyr(nxsize,nzsize),dvdtyr(nxsize,nzsize),  &
    dwdtyr(nxsize,nzsize),dtdtyr(nxsize,nzsize), dddtyr(nxsize,nzsize),  &
    acouyr(nxsize,nzsize),ova2yr(nxsize,nzsize),  &
    gam1yr(nxsize,nzsize),ovgmyr(nxsize,nzsize),  &
    sydtyr(nxsize,nzsize),sorpyr(nxsize,nzsize)
COMMON/nbccyr/bclyyr,stryyr,dydtyr,rateyr,strhyr,  &
    bcl1yr,bcl2yr,bcl3yr,bcl4yr,bcl5yr,bcltyr,  &
    struyr,strvyr,strwyr,strpyr,strdyr,strtyr, streyr,strgyr,strryr,  &
    dudtyr,dvdtyr,dwdtyr,dtdtyr,dddtyr, acouyr,ova2yr,gam1yr,ovgmyr,sydtyr,sorpyr


!     Z-DIRECTION LEFT-HAND END
REAL(kind=8) ::  &
    bclyzl(nxsize,nysize,nspcmx),stryzl(nxsize,nysize,nspcmx),  &
    dydtzl(nxsize,nysize,nspcmx),ratezl(nxsize,nysize,nspcmx),  &
    strhzl(nxsize,nysize,nspcmx), bcl1zl(nxsize,nysize),bcl2zl(nxsize,nysize),  &
    bcl3zl(nxsize,nysize),bcl4zl(nxsize,nysize),  &
    bcl5zl(nxsize,nysize),bcltzl(nxsize,nysize),  &
    struzl(nxsize,nysize),strvzl(nxsize,nysize),  &
    strwzl(nxsize,nysize),strpzl(nxsize,nysize),  &
    strdzl(nxsize,nysize),strtzl(nxsize,nysize),  &
    strezl(nxsize,nysize),strgzl(nxsize,nysize), strrzl(nxsize,nysize),  &
    dudtzl(nxsize,nysize),dvdtzl(nxsize,nysize),  &
    dwdtzl(nxsize,nysize),dtdtzl(nxsize,nysize), dddtzl(nxsize,nysize),  &
    acouzl(nxsize,nysize),ova2zl(nxsize,nysize),  &
    gam1zl(nxsize,nysize),ovgmzl(nxsize,nysize),  &
    sydtzl(nxsize,nysize),sorpzl(nxsize,nysize)
COMMON/nbcczl/bclyzl,stryzl,dydtzl,ratezl,strhzl,  &
    bcl1zl,bcl2zl,bcl3zl,bcl4zl,bcl5zl,bcltzl,  &
    struzl,strvzl,strwzl,strpzl,strdzl,strtzl, strezl,strgzl,strrzl,  &
    dudtzl,dvdtzl,dwdtzl,dtdtzl,dddtzl, acouzl,ova2zl,gam1zl,ovgmzl,sydtzl,sorpzl

!     Z-DIRECTION RIGHT-HAND END
REAL(kind=8) ::  &
    bclyzr(nxsize,nysize,nspcmx),stryzr(nxsize,nysize,nspcmx),  &
    dydtzr(nxsize,nysize,nspcmx),ratezr(nxsize,nysize,nspcmx),  &
    strhzr(nxsize,nysize,nspcmx), bcl1zr(nxsize,nysize),bcl2zr(nxsize,nysize),  &
    bcl3zr(nxsize,nysize),bcl4zr(nxsize,nysize),  &
    bcl5zr(nxsize,nysize),bcltzr(nxsize,nysize),  &
    struzr(nxsize,nysize),strvzr(nxsize,nysize),  &
    strwzr(nxsize,nysize),strpzr(nxsize,nysize),  &
    strdzr(nxsize,nysize),strtzr(nxsize,nysize),  &
    strezr(nxsize,nysize),strgzr(nxsize,nysize), strrzr(nxsize,nysize),  &
    dudtzr(nxsize,nysize),dvdtzr(nxsize,nysize),  &
    dwdtzr(nxsize,nysize),dtdtzr(nxsize,nysize), dddtzr(nxsize,nysize),  &
    acouzr(nxsize,nysize),ova2zr(nxsize,nysize),  &
    gam1zr(nxsize,nysize),ovgmzr(nxsize,nysize),  &
    sydtzr(nxsize,nysize),sorpzr(nxsize,nysize)
COMMON/nbcczr/bclyzr,stryzr,dydtzr,ratezr,strhzr,  &
    bcl1zr,bcl2zr,bcl3zr,bcl4zr,bcl5zr,bcltzr,  &
    struzr,strvzr,strwzr,strpzr,strdzr,strtzr, strezr,strgzr,strrzr,  &
    dudtzr,dvdtzr,dwdtzr,dtdtzr,dddtzr, acouzr,ova2zr,gam1zr,ovgmzr,sydtzr,sorpzr

!     NSBCCL-------------------------------------------------------------------
!     DOMDEC-------------------------------------------------------------------

!     PARALLEL DOMAIN DECOMPOSITION DATA
INTEGER :: npmapx(0:nxproc),npmapy(0:nyproc),npmapz(0:nzproc),  &
    nprocx(0:nxproc),nprocy(0:nyproc),nprocz(0:nzproc),  &
    nxnode,nynode,nznode,nxnbig,nynbig,nznbig, nxprm1,nyprm1,nzprm1,  &
    iproc,nproc,ixproc,iyproc,izproc,  &
    ixprom,ixprop,iyprom,iyprop,izprom,izprop,  &
    istal, istol, jstal, jstol, kstal, kstol,  &
    istap1,istom1,jstap1,jstom1,kstap1,kstom1,  &
    istap2,istom2,jstap2,jstom2,kstap2,kstom2,  &
    istap3,istom3,jstap3,jstom3,kstap3,kstom3,  &
    istap4,istom4,jstap4,jstom4,kstap4,kstom4,  &
    istap5,istom5,jstap5,jstom5,kstap5,kstom5,  &
    istap6,istom6,jstap6,jstom6,kstap6,kstom6,  &
    istap7,istom7,jstap7,jstom7,kstap7,kstom7,  &
    istap8,istom8,jstap8,jstom8,kstap8,kstom8,  &
    istap9,istom9,jstap9,jstom9,kstap9,kstom9,  &
    istapa,istoma,jstapa,jstoma,kstapa,kstoma,  &
    istapb,istomb,jstapb,jstomb,kstapb,kstomb,  &
    istald,istold,jstald,jstold,kstald,kstold,  &
    istalu,istolu,jstalu,jstolu,kstalu,kstolu,  &
    istalv,istolv,jstalv,jstolv,kstalv,kstolv,  &
    istalw,istolw,jstalw,jstolw,kstalw,kstolw,  &
    istale,istole,jstale,jstole,kstale,kstole,  &
    istaly,istoly,jstaly,jstoly,kstaly,kstoly,  &
    istalt,istolt,jstalt,jstolt,kstalt,kstolt,  &
    istab, istob, jstab, jstob, kstab, kstob,  &
    istali,istoli,jstali,jstoli,kstali,kstoli,  &
    istari,istori,jstari,jstori,kstari,kstori,  &
    istalo,istolo,jstalo,jstolo,kstalo,kstolo,  &
    istaro,istoro,jstaro,jstoro,kstaro,kstoro,  &
    istaw, istow, jstaw, jstow, kstaw, kstow,  &
    itgxsl,itgxrl,itgysl,itgyrl,itgzsl,itgzrl,  &
    itgxsr,itgxrr,itgysr,itgyrr,itgzsr,itgzrr

LOGICAL :: prgoxl,prgoxr,prgoyl,prgoyr,prgozl,prgozr, proddx,proddy,proddz

COMMON/domdec/npmapx,npmapy,npmapz,nprocx,nprocy,nprocz,  &
    nxnode,nynode,nznode,nxnbig,nynbig,nznbig, nxprm1,nyprm1,nzprm1,  &
    iproc,nproc,ixproc,iyproc,izproc,  &
    ixprom,ixprop,iyprom,iyprop,izprom,izprop,  &
    istal, istol, jstal, jstol, kstal, kstol,  &
    istap1,istom1,jstap1,jstom1,kstap1,kstom1,  &
    istap2,istom2,jstap2,jstom2,kstap2,kstom2,  &
    istap3,istom3,jstap3,jstom3,kstap3,kstom3,  &
    istap4,istom4,jstap4,jstom4,kstap4,kstom4,  &
    istap5,istom5,jstap5,jstom5,kstap5,kstom5,  &
    istap6,istom6,jstap6,jstom6,kstap6,kstom6,  &
    istap7,istom7,jstap7,jstom7,kstap7,kstom7,  &
    istap8,istom8,jstap8,jstom8,kstap8,kstom8,  &
    istap9,istom9,jstap9,jstom9,kstap9,kstom9,  &
    istapa,istoma,jstapa,jstoma,kstapa,kstoma,  &
    istapb,istomb,jstapb,jstomb,kstapb,kstomb,  &
    istald,istold,jstald,jstold,kstald,kstold,  &
    istalu,istolu,jstalu,jstolu,kstalu,kstolu,  &
    istalv,istolv,jstalv,jstolv,kstalv,kstolv,  &
    istalw,istolw,jstalw,jstolw,kstalw,kstolw,  &
    istale,istole,jstale,jstole,kstale,kstole,  &
    istaly,istoly,jstaly,jstoly,kstaly,kstoly,  &
    istalt,istolt,jstalt,jstolt,kstalt,kstolt,  &
    istab, istob, jstab, jstob, kstab, kstob,  &
    istali,istoli,jstali,jstoli,kstali,kstoli,  &
    istari,istori,jstari,jstori,kstari,kstori,  &
    istalo,istolo,jstalo,jstolo,kstalo,kstolo,  &
    istaro,istoro,jstaro,jstoro,kstaro,kstoro,  &
    istaw, istow, jstaw, jstow, kstaw, kstow,  &
    itgxsl,itgxrl,itgysl,itgyrl,itgzsl,itgzrl,  &
    itgxsr,itgxrr,itgysr,itgyrr,itgzsr,itgzrr,  &
    prgoxl,prgoxr,prgoyl,prgoyr,prgozl,prgozr, proddx,proddy,proddz

!     DOMDEC-------------------------------------------------------------------
!     STATIS-------------------------------------------------------------------

!     DATA FOR ON-LINE STATISTICS

INTEGER :: nstore
PARAMETER(nstore=32)

REAL(kind=8) :: umax(nstore),umin(nstore),vmax(nstore),vmin(nstore),  &
    wmax(nstore),wmin(nstore), ubar(nstore),vbar(nstore),wbar(nstore),  &
    uvar(nstore),vvar(nstore),wvar(nstore),  &
    tke(nstore),dissip(nstore),dissom(nstore),  &
    disot1(nstore),disot2(nstore),disot3(nstore),  &
    residx(nxsize,nstore),residy(nysize,nstore), residz(nzsize,nstore),  &
    xlenth(nstore),ylenth(nstore),zlenth(nstore),tkecor(nstore),  &
    uveloc(nstore),vveloc(nstore),wveloc(nstore),  &
    skewdx(nstore),skewdy(nstore),skewdz(nstore),  &
    tayscx(nstore),tayscy(nstore),tayscz(nstore),  &
    etalen(nstore),taueta(nstore),veleta(nstore),  &
    glorat,fsturb(nstore),fmassb(nstore),  &
    tint(nstore),uprm(nstore),ubarg(nstore),vbarg(nstore),  &
    wbarg(nstore),uvarg(nstore),vvarg(nstore),wvarg(nstore), dissipg(nstore)
INTEGER :: itstim(nstore)

COMMON/statis/umax,umin,vmax,vmin,wmax,wmin, ubar,vbar,wbar,uvar,vvar,wvar,  &
    tke,dissip,dissom,disot1,disot2,disot3, residx,residy,residz,  &
    xlenth,ylenth,zlenth,tkecor, uveloc,vveloc,wveloc,  &
    skewdx,skewdy,skewdz, tayscx,tayscy,tayscz,  &
    etalen,taueta,veleta, glorat,fsturb,fmassb,  &
    tint,uprm,ubarg,vbarg,wbarg, uvarg,vvarg,wvarg,dissipg,  &
    itstim

!     STATIS-------------------------------------------------------------------
!     umod record data for tgv------------------------------
REAL(kind=8) :: dutgvdx(nxsize,nysize,nzsize)
REAL(kind=8) :: dutgvdy(nxsize,nysize,nzsize)
REAL(kind=8) :: dutgvdz(nxsize,nysize,nzsize)
REAL(kind=8) :: dvtgvdx(nxsize,nysize,nzsize)
REAL(kind=8) :: dvtgvdy(nxsize,nysize,nzsize)
REAL(kind=8) :: dvtgvdz(nxsize,nysize,nzsize)
REAL(kind=8) :: dwtgvdx(nxsize,nysize,nzsize)
REAL(kind=8) :: dwtgvdy(nxsize,nysize,nzsize)
REAL(kind=8) :: dwtgvdz(nxsize,nysize,nzsize)

common/tgv/dutgvdx,dutgvdy,dutgvdz,dvtgvdx,dvtgvdy,dvtgvdz, &
           dwtgvdx,dwtgvdy,dwtgvdz

END MODULE com_senga

