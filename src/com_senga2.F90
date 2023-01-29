MODULE com_senga

!     PHYDIM-------------------------------------------------------------------
 
!     PHYSICAL DIMENSIONS OF ARRAYS
!     -----------------------------
!     NOTE: ALL ARRAY SIZES MUST BE CONSISTENT
!     NXSIZE MUST BE >= NXGLBL/NXPROC
!     NYSIZE MUST BE >= NYGLBL/NYPROC
!     NZSIZE MUST BE >= NZGLBL/NZPROC
!     WITH AN EXTRA ALLOWANCE FOR ANY REMAINDER

    use data_types

!     GLOBAL GRID SIZE
INTEGER :: nxglbl,nyglbl,nzglbl
PARAMETER(nxglbl=2000, nyglbl=1, nzglbl=1)
INTEGER :: ngzmax
!     SET NGZMAX=MAX(NXGLBL,NYGLBL,NZGLBL)
PARAMETER(ngzmax=nxglbl)

!     NUMBER OF PROCESSORS
INTEGER :: nxproc,nyproc,nzproc
PARAMETER(nxproc=1, nyproc=1, nzproc=1)
INTEGER :: nprmax
!     SET NPRMAX=MAX(NXPROC,NYPROC,NZPROC)
PARAMETER(nprmax=nxproc)

!     LOCAL GRID SIZE
INTEGER :: nxsize,nysize,nzsize
PARAMETER(nxsize=2000, nysize=1, nzsize=1)
INTEGER :: nszmax
!     SET NSZMAX=MAX(NXSIZE,NYSIZE,NZSIZE)
PARAMETER(nszmax=nxsize)

!     SIZE OF HALO
INTEGER :: nhalox,nhaloy,nhaloz
PARAMETER(nhalox=5,nhaloy=0,nhaloz=0)

!     SIZE OF PARALLEL TRANSFER ARRAY
!     NPARAY MUST BE >= MAX(NHALOX,NHALOY,NHALOZ)
!                            *MAX((NXSIZE+2*NHALOX)*(NYSIZE+2*NHALOY),
!                                 (NXSIZE+2*NHALOX)*(NZSIZE+2*NHALOZ),
!                                 (NYSIZE+2*NHALOY)*(NZSIZE+2*NHALOZ))
!     AND ALSO LARGE ENOUGH FOR PARALLEL BROADCAST OF CHEMICAL DATA
!     AND ALSO LARGE ENOUGH FOR PARALLEL TRANSFER OF INITIAL TURBULENCE DATA
INTEGER :: nparay
PARAMETER(nparay=80000)

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

REAL(KIND=8) :: fftrow(2*ngzmax,3,npenmx), ftpart(2*nszmax,3,npenmx),  &
    fftinx(2*ngzmax)

COMMON/ifturb/fftrow,ftpart,fftinx

!     IFTURB-------------------------------------------------------------------
!     CHEMIC-------------------------------------------------------------------

!     PARAMETERS
!     ==========
!     MAX NO OF SPECIES, NO OF STEPS
INTEGER :: nspcmx,nstpmx
PARAMETER(nspcmx=9, nstpmx=21)

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

!     UNIVERSAL GAS CONSTANT
REAL(KIND=8) :: rguniv
PARAMETER(rguniv=8314.20_8)

!     CHEMICAL DATA
!     =============
REAL(KIND=8) :: amolcp(ncofmx,ntinmx,nspcmx)
REAL(KIND=8) :: amolgb(ncofmx,ntinmx,nspcmx)
REAL(KIND=8) :: amascp(ncofmx,ntinmx,nspcmx)
REAL(KIND=8) :: amasct(ncofmx,ntinmx,nspcmx)
REAL(KIND=8) :: amasch(ncofmx,ntinmx,nspcmx)
REAL(KIND=8) :: amascs(ncofmx,ntinmx,nspcmx)
REAL(KIND=8) :: tintlo(ntinmx,nspcmx),tinthi(ntinmx,nspcmx)
REAL(KIND=8) :: diffmu(nssmax,nstpmx),diffmw(nssmax,nstpmx)
REAL(KIND=8) :: crspec(nssmax,nstpmx),cpspec(nssmax,nstpmx)
REAL(KIND=8) :: effy3b(nspcmx,nbdymx)
REAL(KIND=8) :: rparam(3,nstpmx)
REAL(KIND=8) :: rclind(4,nllmax)
REAL(KIND=8) :: rctroe(12,nllmax)
REAL(KIND=8) :: rcsrif(8,nllmax)
REAL(KIND=8) :: wmolar(nspcmx),ovwmol(nspcmx),rgspec(nspcmx)
REAL(KIND=8) :: clewis(nspcmx),olewis(nspcmx)
REAL(KIND=8) :: prefgb

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
    wmolar,ovwmol,rgspec,clewis,olewis, prefgb,  &
    nsspec,nrspec,npspec,nrcpec,npcpec, ncofcp,ncpoly,ncpom1,ncenth,ncenpy,  &
    nsslen,nrslen,npslen,nrclen,npclen, mblist,mglist,mllist,mtlist,mslist,  &
    ntint, nspec,nspm1,nstep,nbody,ngibb,nlind,ntroe,nsrif,  &
    spcsym,bdysym

!     CHEMIC-------------------------------------------------------------------
!     DIFFUS-------------------------------------------------------------------

!     MIXTURE AVERAGED TRANSPORT

!     MASS FRACTION TOLERANCE
REAL(KIND=8) :: dfctol
PARAMETER(dfctol = 0.0000000000010_8)

!     MAX NUMBERS OF POLYNOMIAL COEFFICIENTS
INTEGER :: ndcfmx,nvcfmx,nccfmx
PARAMETER(ndcfmx = 4, nvcfmx = 4, nccfmx = 4)

!     MOLECULAR TRANSPORT DATA
REAL(KIND=8) :: diffco(ndcfmx,nspcmx,nspcmx)
REAL(KIND=8) :: tdrcco(ndcfmx,nspcmx,nspcmx)
REAL(KIND=8) :: wilko1(nspcmx,nspcmx),wilko2(nspcmx,nspcmx)
REAL(KIND=8) :: viscco(nvcfmx,nspcmx)
REAL(KIND=8) :: condco(nccfmx,nspcmx)
REAL(KIND=8) :: pdifgb,tdifgb
REAL(KIND=8) :: wmltdr

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

REAL(KIND=8) :: stefbo
PARAMETER(stefbo = 0.00000005670373210_8)

INTEGER :: ncfrmx
PARAMETER(ncfrmx = 6)

REAL(KIND=8) :: akprad(ncfrmx,nspcmx)
REAL(KIND=8) :: trefrn

INTEGER :: nsprid(nspcmx),nkprad(nspcmx),nkprm1(nspcmx)
INTEGER :: nsprad,nfradn

LOGICAL :: flradn

COMMON/radiat/akprad,trefrn,nsprid,nkprad,nkprm1,nsprad,  &
    nfradn, flradn

!     RADIAT-------------------------------------------------------------------
!     PARAMS-------------------------------------------------------------------

!     PHYSICAL AND NUMERICAL PARAMETERS

!     MASS FRACTION TOLERANCE
REAL(KIND=8) :: ytoler
PARAMETER(ytoler = 1.0E-10)

!     NUMBERS
REAL(KIND=8) :: zero,one,two,three,four,eight,  &
    half,thrf,thrd,tthd,fthd,qrtr
PARAMETER(zero=0.0_8, one=1.0_8, two=2.0_8, three=3.0_8,  &
    four=4.0_8, eight=8.0_8)
PARAMETER(half=one/two,  thrf=three*half, thrd=one/three,  &
    tthd=two*thrd, fthd=two*tthd,   qrtr=one/four)

!     LOCATION FOR REFERENCE PRESSURE VALUE
INTEGER :: ipref,jpref,kpref
PARAMETER(ipref=1,jpref=1,kpref=1)

!     PHYSICAL AND NUMERICAL DATA
REAL(KIND=8) :: yrin(nspcmx)
REAL(KIND=8) :: prin,trin,drin,urin,vrin,wrin,erin
REAL(KIND=8) :: xgdlen,ygdlen,zgdlen
REAL(KIND=8) :: etime,tstep,deltax,deltay,deltaz,pi,clnten
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
REAL(KIND=8) :: acoeff,bcoeff,ccoeff,dcoeff,ecoeff,  &
    acoef1,bcoef1,ccoef1,dcoef1,  &
    acoef2,bcoef2,ccoef2,dcoef2,  &
    acoef3,bcoef3,  &
    acoef4,bcoef4,ccoef4,  &
    acoef5,bcoef5,ccoef5,dcoef5,  &
    acoefs,bcoefs,ccoefs,dcoefs,ecoefs,  &
    acofs1,bcofs1,ccofs1,dcofs1,ecofs1,  &
    acofs2,bcofs2,ccofs2,dcofs2,ecofs2,  &
    acofs3,bcofs3,  &
    acofs4,bcofs4,ccofs4,  &
    acofs5,bcofs5,ccofs5,dcofs5,  &
    acoefx,bcoefx,ccoefx,dcoefx,ecoefx,  &
    acofxz,bcofxz,ccofxz,dcofxz,ecofxz,  &
    acofyz,bcofyz,ccofyz,dcofyz,ecofyz,  &
    acf1xz,bcf1xz,ccf1xz,dcf1xz, acf1yz,bcf1yz,ccf1yz,dcf1yz,  &
    acf2xz,bcf2xz,ccf2xz,dcf2xz, &
    acf2yz,bcf2yz,ccf2yz,dcf2yz, &
    acf3xz,bcf3xz, acf3yz,bcf3yz,  &
    acf4xz,bcf4xz,ccf4xz, &
    acf4yz,bcf4yz,ccf4yz, &
    acf5xz,bcf5xz,ccf5xz,dcf5xz, acf5yz,bcf5yz,ccf5yz,dcf5yz,  &
    acofc1,bcofc1,ccofc1,dcofc1, acofc2,bcofc2,ccofc2,dcofc2,  &
    acc1xz,bcc1xz,ccc1xz,dcc1xz, &
    acc1yz,bcc1yz,ccc1yz,dcc1yz, &
    acc2xz,bcc2xz,ccc2xz,dcc2xz, acc2yz,bcc2yz,ccc2yz,dcc2yz

!     SPATIAL DERIVATIVE END CONDITIONS
INTEGER :: nendxl,nendxr,nendyl,nendyr,nendzl,nendzr

COMMON/dfdiff/acoeff,bcoeff,ccoeff,dcoeff,ecoeff,  &
    acoef1,bcoef1,ccoef1,dcoef1,  &
    acoef2,bcoef2,ccoef2,dcoef2,  &
    acoef3,bcoef3,  &
    acoef4,bcoef4,ccoef4,  &
    acoef5,bcoef5,ccoef5,dcoef5,  &
    acoefs,bcoefs,ccoefs,dcoefs,ecoefs,  &
    acofs1,bcofs1,ccofs1,dcofs1,ecofs1,  &
    acofs2,bcofs2,ccofs2,dcofs2,ecofs2,  &
    acofs3,bcofs3,  &
    acofs4,bcofs4,ccofs4,  &
    acofs5,bcofs5,ccofs5,dcofs5,  &
    acoefx,bcoefx,ccoefx,dcoefx,ecoefx,  &
    acofxz,bcofxz,ccofxz,dcofxz,ecofxz,  &
    acofyz,bcofyz,ccofyz,dcofyz,ecofyz,  &
    acf1xz,bcf1xz,ccf1xz,dcf1xz, acf1yz,bcf1yz,ccf1yz,dcf1yz,  &
    acf2xz,bcf2xz,ccf2xz,dcf2xz, &
    acf2yz,bcf2yz,ccf2yz,dcf2yz, &
    acf3xz,bcf3xz, acf3yz,bcf3yz,  &
    acf4xz,bcf4xz,ccf4xz, &
    acf4yz,bcf4yz,ccf4yz, &
    acf5xz,bcf5xz,ccf5xz,dcf5xz, acf5yz,bcf5yz,ccf5yz,dcf5yz,  &
    acofc1,bcofc1,ccofc1,dcofc1, acofc2,bcofc2,ccofc2,dcofc2,  &
    acc1xz,bcc1xz,ccc1xz,dcc1xz, &
    acc1yz,bcc1yz,ccc1yz,dcc1yz, &
    acc2xz,bcc2xz,ccc2xz,dcc2xz, acc2yz,bcc2yz,ccc2yz,dcc2yz,  &
    nendxl,nendxr,nendyl,nendyr,nendzl,nendzr

!     DFDIFF-------------------------------------------------------------------
!     EMSTRT-------------------------------------------------------------------

!     DATA FOR MESH STRETCHING
REAL(KIND=8) :: gcmreg(nszmax),gcmstr(nszmax),  &
    dgdhat(nszmax),dgdhsq(nszmax),d2gdh2(nszmax)

COMMON/emstrt/gcmreg,gcmstr,dgdhat,dgdhsq,d2gdh2

!     EMSTRT-------------------------------------------------------------------
!     FILCOM-------------------------------------------------------------------

!     SPATIAL FILTER COEFFICIENTS
REAL(KIND=8) :: facofx,fbcofx,fccofx,fdcofx,fecofx,ffcofx,fgcofx,  &
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

REAL(KIND=8) :: rklhs(nrkmax),rkrhs(nrkmax)
REAL(KIND=8) :: rkerr(nrkmax),rktim(nrkmax)
REAL(KIND=8) :: ctmult,ctalph,ctbeta,ctgama
REAL(KIND=8) :: errtol,errlow,trmin,trmax,tsmin,tsmax
REAL(KIND=8) :: errold,errldr,btime
REAL(KIND=8) :: erdnrm,erunrm,ervnrm,erwnrm,erenrm
REAL(KIND=8) :: erynrm(nspcmx)
INTEGER :: irkstp,nrkstp,nrksm1
INTEGER :: inderr
LOGICAL :: fupelc

COMMON/runkut/rklhs,rkrhs,rkerr,rktim, ctmult,ctalph,ctbeta,ctgama,  &
    errtol,errlow,trmin,trmax,tsmin,tsmax, errold,errldr,btime,  &
    erdnrm,erunrm,ervnrm,erwnrm,erenrm,erynrm, irkstp,nrkstp,nrksm1,inderr,fupelc

!     RUNKUT-------------------------------------------------------------------
!     RKNORM-------------------------------------------------------------------

!     RUNGE-KUTTA SUBSTEP ERROR NORMS
REAL(KIND=8) :: erdrhs(nrkmax),erurhs(nrkmax),ervrhs(nrkmax),  &
    erwrhs(nrkmax),ererhs(nrkmax)
REAL(KIND=8) :: eryrhs(nspcmx,nrkmax)

COMMON/rknorm/erdrhs,erurhs,ervrhs,erwrhs,ererhs,eryrhs

!     PRINCIPAL VARIABLES (ALL STANDARD SIZE ARRAYS)
real(kind=8), dimension(:,:,:), allocatable :: drun,urun,vrun,wrun,erun
real(kind=8), dimension(:,:,:,:), allocatable :: yrun

!     WORKSPACE (STANDARD SIZE ARRAYS)
real(kind=8), dimension(:,:,:,:), allocatable :: rrte

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     WORKSPACE (BIGGER SIZE ARRAYS)
real(kind=8), dimension(:,:,:), allocatable :: utmp,vtmp,wtmp,prun,trun

!     WORKSPACE (PARALLEL TRANSFER ARRAY)
REAL(KIND=8) :: parray(nparay)

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
REAL(KIND=8) :: rxlprm(nbcprr),rxrprm(nbcprr),  &
    rylprm(nbcprr),ryrprm(nbcprr), rzlprm(nbcprr),rzrprm(nbcprr)
REAL(KIND=8) :: cobcxl,cobcxr,cobcyl,cobcyr,cobczl,cobczr,  &
    pinfxl,pinfxr,pinfyl,pinfyr,pinfzl,pinfzr
REAL(KIND=8) :: slocxl,slocxr,slocyl,slocyr,sloczl,sloczr,  &
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
real(kind=8), dimension(:,:,:), allocatable :: ufxl,vfxl,wfxl

!     WALL BC DIFFERENCING DATA
INTEGER :: ncbcsz
PARAMETER(ncbcsz=5)
REAL(KIND=8) :: acbcxl(ncbcsz),acbcxr(ncbcsz),  &
    acbcyl(ncbcsz),acbcyr(ncbcsz), acbczl(ncbcsz),acbczr(ncbcsz)
COMMON/bcdifw/acbcxl,acbcxr,acbcyl,acbcyr,acbczl,acbczr

!     X-DIRECTION LEFT-HAND END
real(kind=8), dimension(:,:,:), allocatable :: struxl,strvxl,strwxl,dudtxl,dvdtxl,dwdtxl

!     NSBCCL-------------------------------------------------------------------
!     DOMDEC-------------------------------------------------------------------

!     PARALLEL DOMAIN DECOMPOSITION DATA
INTEGER :: npmapx(0:nxproc),npmapy(0:nyproc),npmapz(0:nzproc),  &
    nprocx(0:nxproc),nprocy(0:nyproc),nprocz(0:nzproc),  &
    nxnode,nynode,nznode,nxnbig,nynbig,nznbig, nxprm1,nyprm1,nzprm1,  &
    iproc,nproc,ixproc,iyproc,izproc,  &
    ixprom,ixprop,iyprom,iyprop,izprom,izprop,  &
    itgxsl,itgxrl,itgysl,itgyrl,itgzsl,itgzrl,  &
    itgxsr,itgxrr,itgysr,itgyrr,itgzsr,itgzrr

LOGICAL :: prgoxl,prgoxr,prgoyl,prgoyr,prgozl,prgozr, proddx,proddy,proddz

COMMON/domdec/npmapx,npmapy,npmapz,nprocx,nprocy,nprocz,  &
    nxnode,nynode,nznode,nxnbig,nynbig,nznbig, nxprm1,nyprm1,nzprm1,  &
    iproc,nproc,ixproc,iyproc,izproc,  &
    ixprom,ixprop,iyprom,iyprop,izprom,izprop,  &
    itgxsl,itgxrl,itgysl,itgyrl,itgzsl,itgzrl,  &
    itgxsr,itgxrr,itgysr,itgyrr,itgzsr,itgzrr,  &
    prgoxl,prgoxr,prgoyl,prgoyr,prgozl,prgozr, proddx,proddy,proddz

!     DOMDEC-------------------------------------------------------------------
!     STATIS-------------------------------------------------------------------

!     DATA FOR ON-LINE STATISTICS

INTEGER :: nstore
PARAMETER(nstore=32)

REAL(KIND=8) :: umax(nstore),umin(nstore),vmax(nstore),vmin(nstore),  &
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
END MODULE com_senga

