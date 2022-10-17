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

    use data_types

!     GLOBAL GRID SIZE
INTEGER :: nxglbl,nyglbl,nzglbl
PARAMETER(nxglbl=2000, nyglbl=1, nzglbl=1)
INTEGER :: ngzmax
!     SET NGZMAX=MAX(NXGLBL,NYGLBL,NZGLBL)
PARAMETER(ngzmax=nxglbl)

!     NUMBER OF PROCESSORS
INTEGER :: nxproc,nyproc,nzproc
PARAMETER(nxproc=40, nyproc=1, nzproc=1)
INTEGER :: nprmax
!     SET NPRMAX=MAX(NXPROC,NYPROC,NZPROC)
PARAMETER(nprmax=nxproc)

!     LOCAL GRID SIZE
INTEGER :: nxsize,nysize,nzsize
PARAMETER(nxsize=50, nysize=1, nzsize=1)
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

REAL(KIND=dp) :: fftrow(2*ngzmax,3,npenmx), ftpart(2*nszmax,3,npenmx),  &
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

!     TRANSPORT COEFFICIENTS
REAL(KIND=dp) :: alamdc,rlamda,tlamda
PARAMETER(alamdc=0.0000258_dp, rlamda=0.70_dp, tlamda=298.0_dp)
REAL(KIND=dp) :: prantl
PARAMETER(prantl=0.70_dp)

!     UNIVERSAL GAS CONSTANT
REAL(KIND=dp) :: rguniv
PARAMETER(rguniv=8314.20_dp)


!     CHEMICAL DATA
!     =============
REAL(KIND=dp) :: amolcp(ncofmx,ntinmx,nspcmx)
REAL(KIND=dp) :: amolgb(ncofmx,ntinmx,nspcmx)
REAL(KIND=dp) :: amascp(ncofmx,ntinmx,nspcmx)
REAL(KIND=dp) :: amasct(ncofmx,ntinmx,nspcmx)
REAL(KIND=dp) :: amasch(ncofmx,ntinmx,nspcmx)
REAL(KIND=dp) :: amascs(ncofmx,ntinmx,nspcmx)
REAL(KIND=dp) :: tintlo(ntinmx,nspcmx),tinthi(ntinmx,nspcmx)
REAL(KIND=dp) :: diffmu(nssmax,nstpmx),diffmw(nssmax,nstpmx)
REAL(KIND=dp) :: crspec(nssmax,nstpmx),cpspec(nssmax,nstpmx)
REAL(KIND=dp) :: effy3b(nspcmx,nbdymx)
REAL(KIND=dp) :: rparam(3,nstpmx)
REAL(KIND=dp) :: rclind(4,nllmax)
REAL(KIND=dp) :: rctroe(12,nllmax)
REAL(KIND=dp) :: rcsrif(8,nllmax)
REAL(KIND=dp) :: wmolar(nspcmx),ovwmol(nspcmx),rgspec(nspcmx)
REAL(KIND=dp) :: clewis(nspcmx),olewis(nspcmx)
REAL(KIND=dp) :: prefgb,alamda

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
REAL(KIND=dp) :: dfctol
PARAMETER(dfctol = 0.0000000000010_dp)

!     MAX NUMBERS OF POLYNOMIAL COEFFICIENTS
INTEGER :: ndcfmx,nvcfmx,nccfmx
PARAMETER(ndcfmx = 4, nvcfmx = 4, nccfmx = 4)

!     MOLECULAR TRANSPORT DATA
REAL(KIND=dp) :: diffco(ndcfmx,nspcmx,nspcmx)
REAL(KIND=dp) :: tdrcco(ndcfmx,nspcmx,nspcmx)
REAL(KIND=dp) :: wilko1(nspcmx,nspcmx),wilko2(nspcmx,nspcmx)
REAL(KIND=dp) :: viscco(nvcfmx,nspcmx)
REAL(KIND=dp) :: condco(nccfmx,nspcmx)
REAL(KIND=dp) :: pdifgb,tdifgb
REAL(KIND=dp) :: wmltdr

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

REAL(KIND=dp) :: stefbo
PARAMETER(stefbo = 0.00000005670373210_dp)

INTEGER :: ncfrmx
PARAMETER(ncfrmx = 6)

REAL(KIND=dp) :: akprad(ncfrmx,nspcmx)
REAL(KIND=dp) :: foursb,trefrn,trfrth

INTEGER :: nsprid(nspcmx),nkprad(nspcmx),nkprm1(nspcmx)
INTEGER :: nsprad,nfradn

LOGICAL :: flradn

COMMON/radiat/akprad,foursb,trefrn,trfrth, nsprid,nkprad,nkprm1,nsprad,  &
    nfradn, flradn

!     RADIAT-------------------------------------------------------------------
!     PARAMS-------------------------------------------------------------------

!     PHYSICAL AND NUMERICAL PARAMETERS

!     MASS FRACTION TOLERANCE
REAL(KIND=dp) :: ytoler
PARAMETER(ytoler = 1.0E-10)

!     NUMBERS
REAL(KIND=dp) :: zero,one,two,three,four,eight,  &
    half,thrf,thrd,tthd,fthd,qrtr
PARAMETER(zero=0.0_dp, one=1.0_dp, two=2.0_dp, three=3.0_dp,  &
    four=4.0_dp, eight=8.0_dp)
PARAMETER(half=one/two,  thrf=three*half, thrd=one/three,  &
    tthd=two*thrd, fthd=two*tthd,   qrtr=one/four)

!     LOCATION FOR REFERENCE PRESSURE VALUE
INTEGER :: ipref,jpref,kpref
PARAMETER(ipref=1,jpref=1,kpref=1)

!     PHYSICAL AND NUMERICAL DATA
REAL(KIND=dp) :: yrin(nspcmx)
REAL(KIND=dp) :: prin,trin,drin,urin,vrin,wrin,erin
REAL(KIND=dp) :: xgdlen,ygdlen,zgdlen
REAL(KIND=dp) :: etime,tstep,deltax,deltay,deltaz,pi,clnten
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
REAL(KIND=dp) :: acoeff,bcoeff,ccoeff,dcoeff,ecoeff,  &
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
    nendxl,nendxr,nendyl,nendyr,nendzl,nendzr

!     DFDIFF-------------------------------------------------------------------
!     EMSTRT-------------------------------------------------------------------

!     DATA FOR MESH STRETCHING
REAL(KIND=dp) :: gcmreg(nszmax),gcmstr(nszmax),  &
    dgdhat(nszmax),dgdhsq(nszmax),d2gdh2(nszmax)

COMMON/emstrt/gcmreg,gcmstr,dgdhat,dgdhsq,d2gdh2

!     EMSTRT-------------------------------------------------------------------
!     FILCOM-------------------------------------------------------------------

!     SPATIAL FILTER COEFFICIENTS
REAL(KIND=dp) :: facofx,fbcofx,fccofx,fdcofx,fecofx,ffcofx,fgcofx,  &
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

REAL(KIND=dp) :: rklhs(nrkmax),rkrhs(nrkmax)
REAL(KIND=dp) :: rkerr(nrkmax),rktim(nrkmax)
REAL(KIND=dp) :: ctmult,ctalph,ctbeta,ctgama
REAL(KIND=dp) :: errtol,errlow,trmin,trmax,tsmin,tsmax
REAL(KIND=dp) :: errold,errldr,btime
REAL(KIND=dp) :: erdnrm,erunrm,ervnrm,erwnrm,erenrm
REAL(KIND=dp) :: erynrm(nspcmx)
INTEGER :: irkstp,nrkstp,nrksm1
INTEGER :: inderr
LOGICAL :: fupelc

COMMON/runkut/rklhs,rkrhs,rkerr,rktim, ctmult,ctalph,ctbeta,ctgama,  &
    errtol,errlow,trmin,trmax,tsmin,tsmax, errold,errldr,btime,  &
    erdnrm,erunrm,ervnrm,erwnrm,erenrm,erynrm, irkstp,nrkstp,nrksm1,inderr,fupelc

!     RUNKUT-------------------------------------------------------------------
!     RKNORM-------------------------------------------------------------------

!     RUNGE-KUTTA SUBSTEP ERROR NORMS
REAL(KIND=dp) :: erdrhs(nrkmax),erurhs(nrkmax),ervrhs(nrkmax),  &
    erwrhs(nrkmax),ererhs(nrkmax)
REAL(KIND=dp) :: eryrhs(nspcmx,nrkmax)

COMMON/rknorm/erdrhs,erurhs,ervrhs,erwrhs,ererhs,eryrhs

!     RKNORM-------------------------------------------------------------------
!     PRVARS-------------------------------------------------------------------

!     PRINCIPAL VARIABLES (ALL STANDARD SIZE ARRAYS)
REAL(KIND=dp) :: drun(nxsize,nysize,nzsize),urun(nxsize,nysize,nzsize),  &
    vrun(nxsize,nysize,nzsize),wrun(nxsize,nysize,nzsize),  &
    erun(nxsize,nysize,nzsize), yrun(nspcmx,nxsize,nysize,nzsize)

COMMON/prvars/drun,urun,vrun,wrun,erun,yrun

!     PRVARS-------------------------------------------------------------------
!     TIMRHS-------------------------------------------------------------------

!     TIME-STEPPING RHS TERMS (ALL BIGGER SIZE ARRAYS)
REAL(KIND=dp) :: drhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    urhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    vrhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    wrhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    erhs(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    yrhs(nspcmx,nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)

COMMON/timrhs/drhs,urhs,vrhs,wrhs,erhs,yrhs

!     TIMRHS-------------------------------------------------------------------
!     TIMERR-------------------------------------------------------------------

!     TIME-STEPPING ERROR ARRAYS (ALL STANDARD SIZE ARRAYS)
REAL(KIND=dp) :: derr(nxsize,nysize,nzsize),uerr(nxsize,nysize,nzsize),  &
    verr(nxsize,nysize,nzsize),werr(nxsize,nysize,nzsize),  &
    eerr(nxsize,nysize,nzsize), yerr(nspcmx,nxsize,nysize,nzsize)

COMMON/timerr/derr,uerr,verr,werr,eerr,yerr

!     TIMERR-------------------------------------------------------------------
!     WORKSP-------------------------------------------------------------------

!     WORKSPACE (STANDARD SIZE ARRAYS)
REAL(KIND=dp) ::  &
    store1(nxsize,nysize,nzsize),store2(nxsize,nysize,nzsize),  &
    store3(nxsize,nysize,nzsize),store4(nxsize,nysize,nzsize),  &
    store5(nxsize,nysize,nzsize),store6(nxsize,nysize,nzsize),  &
    divm(nxsize,nysize,nzsize), rate(nspcmx,nxsize,nysize,nzsize),  &
    rrte(nspcmx,nxsize,nysize,nzsize)

COMMON/worksp/store1,store2,store3,store4,store5,store6,divm,rate, rrte

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     WORKSPACE (DIFFUSIVE CORRECTION VELOCITY, STANDARD SIZE ARRAYS)
REAL(KIND=dp) :: ucor(nxsize,nysize,nzsize),vcor(nxsize,nysize,nzsize),  &
    wcor(nxsize,nysize,nzsize)

COMMON/worksc/ucor,vcor,wcor
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     WORKSPACE (MIXTURE AVERAGED TRANSPORT, STANDARD SIZE ARRAYS)
REAL(KIND=dp) :: wd1x(nxsize,nysize,nzsize),wd1y(nxsize,nysize,nzsize),  &
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
REAL(KIND=dp) :: wmomix(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    difmix(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
    tdrmix(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)

COMMON/workbt/wmomix,difmix,tdrmix
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     WORKSPACE (BIGGER SIZE ARRAYS)
REAL(KIND=dp) :: utmp(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr),  &
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
REAL(KIND=dp) :: parray(nparay)

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
REAL(KIND=dp) :: rxlprm(nbcprr),rxrprm(nbcprr),  &
    rylprm(nbcprr),ryrprm(nbcprr), rzlprm(nbcprr),rzrprm(nbcprr)
REAL(KIND=dp) :: cobcxl,cobcxr,cobcyl,cobcyr,cobczl,cobczr,  &
    pinfxl,pinfxr,pinfyl,pinfyr,pinfzl,pinfzr
REAL(KIND=dp) :: slocxl,slocxr,slocyl,slocyr,sloczl,sloczr,  &
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
REAL(KIND=dp) :: ufxl(nxsize,nysize,nzsize), vfxl(nxsize,nysize,nzsize),  &
    wfxl(nxsize,nysize,nzsize)
COMMON/nbcinv/ufxl,vfxl,wfxl

!     WALL BC DIFFERENCING DATA
INTEGER :: ncbcsz
PARAMETER(ncbcsz=5)
REAL(KIND=dp) :: acbcxl(ncbcsz),acbcxr(ncbcsz),  &
    acbcyl(ncbcsz),acbcyr(ncbcsz), acbczl(ncbcsz),acbczr(ncbcsz)
COMMON/bcdifw/acbcxl,acbcxr,acbcyl,acbcyr,acbczl,acbczr

!     X-DIRECTION LEFT-HAND END
REAL(KIND=dp) ::  &
    bclyxl(nspcmx,1,nysize,nzsize),stryxl(nspcmx,1,nysize,nzsize),  &
    dydtxl(nspcmx,1,nysize,nzsize),ratexl(nspcmx,1,nysize,nzsize),  &
    strhxl(nspcmx,1,nysize,nzsize), bcl1xl(1,nysize,nzsize),bcl2xl(1,nysize,nzsize),  &
    bcl3xl(1,nysize,nzsize),bcl4xl(1,nysize,nzsize),  &
    bcl5xl(1,nysize,nzsize),bcltxl(1,nysize,nzsize),  &
    struxl(1,nysize,nzsize),strvxl(1,nysize,nzsize),  &
    strwxl(1,nysize,nzsize),strpxl(1,nysize,nzsize),  &
    strdxl(1,nysize,nzsize),strtxl(1,nysize,nzsize),  &
    strexl(1,nysize,nzsize),strgxl(1,nysize,nzsize), strrxl(1,nysize,nzsize),  &
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
REAL(KIND=dp) ::  &
    bclyxr(nspcmx,1,nysize,nzsize),stryxr(nspcmx,1,nysize,nzsize),  &
    dydtxr(nspcmx,1,nysize,nzsize),ratexr(nspcmx,1,nysize,nzsize),  &
    strhxr(nspcmx,1,nysize,nzsize), bcl1xr(1,nysize,nzsize),bcl2xr(1,nysize,nzsize),  &
    bcl3xr(1,nysize,nzsize),bcl4xr(1,nysize,nzsize),  &
    bcl5xr(1,nysize,nzsize),bcltxr(1,nysize,nzsize),  &
    struxr(1,nysize,nzsize),strvxr(1,nysize,nzsize),  &
    strwxr(1,nysize,nzsize),strpxr(1,nysize,nzsize),  &
    strdxr(1,nysize,nzsize),strtxr(1,nysize,nzsize),  &
    strexr(1,nysize,nzsize),strgxr(1,nysize,nzsize), strrxr(1,nysize,nzsize),  &
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
REAL(KIND=dp) ::  &
    bclyyl(nspcmx,nxsize,1,nzsize),stryyl(nspcmx,nxsize,1,nzsize),  &
    dydtyl(nspcmx,nxsize,1,nzsize),rateyl(nspcmx,nxsize,1,nzsize),  &
    strhyl(nspcmx,nxsize,1,nzsize), bcl1yl(nxsize,1,nzsize),bcl2yl(nxsize,1,nzsize),  &
    bcl3yl(nxsize,1,nzsize),bcl4yl(nxsize,1,nzsize),  &
    bcl5yl(nxsize,1,nzsize),bcltyl(nxsize,1,nzsize),  &
    struyl(nxsize,1,nzsize),strvyl(nxsize,1,nzsize),  &
    strwyl(nxsize,1,nzsize),strpyl(nxsize,1,nzsize),  &
    strdyl(nxsize,1,nzsize),strtyl(nxsize,1,nzsize),  &
    streyl(nxsize,1,nzsize),strgyl(nxsize,1,nzsize), strryl(nxsize,1,nzsize),  &
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
REAL(KIND=dp) ::  &
    bclyyr(nspcmx,nxsize,1,nzsize),stryyr(nspcmx,nxsize,1,nzsize),  &
    dydtyr(nspcmx,nxsize,1,nzsize),rateyr(nspcmx,nxsize,1,nzsize),  &
    strhyr(nspcmx,nxsize,1,nzsize), bcl1yr(nxsize,1,nzsize),bcl2yr(nxsize,1,nzsize),  &
    bcl3yr(nxsize,1,nzsize),bcl4yr(nxsize,1,nzsize),  &
    bcl5yr(nxsize,1,nzsize),bcltyr(nxsize,1,nzsize),  &
    struyr(nxsize,1,nzsize),strvyr(nxsize,1,nzsize),  &
    strwyr(nxsize,1,nzsize),strpyr(nxsize,1,nzsize),  &
    strdyr(nxsize,1,nzsize),strtyr(nxsize,1,nzsize),  &
    streyr(nxsize,1,nzsize),strgyr(nxsize,1,nzsize), strryr(nxsize,1,nzsize),  &
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
REAL(KIND=dp) ::  &
    bclyzl(nspcmx,nxsize,nysize,1),stryzl(nspcmx,nxsize,nysize,1),  &
    dydtzl(nspcmx,nxsize,nysize,1),ratezl(nspcmx,nxsize,nysize,1),  &
    strhzl(nspcmx,nxsize,nysize,1), bcl1zl(nxsize,nysize,1),bcl2zl(nxsize,nysize,1),  &
    bcl3zl(nxsize,nysize,1),bcl4zl(nxsize,nysize,1),  &
    bcl5zl(nxsize,nysize,1),bcltzl(nxsize,nysize,1),  &
    struzl(nxsize,nysize,1),strvzl(nxsize,nysize,1),  &
    strwzl(nxsize,nysize,1),strpzl(nxsize,nysize,1),  &
    strdzl(nxsize,nysize,1),strtzl(nxsize,nysize,1),  &
    strezl(nxsize,nysize,1),strgzl(nxsize,nysize,1), strrzl(nxsize,nysize,1),  &
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
REAL(KIND=dp) ::  &
    bclyzr(nspcmx,nxsize,nysize,1),stryzr(nspcmx,nxsize,nysize,1),  &
    dydtzr(nspcmx,nxsize,nysize,1),ratezr(nspcmx,nxsize,nysize,1),  &
    strhzr(nspcmx,nxsize,nysize,1), bcl1zr(nxsize,nysize,1),bcl2zr(nxsize,nysize,1),  &
    bcl3zr(nxsize,nysize,1),bcl4zr(nxsize,nysize,1),  &
    bcl5zr(nxsize,nysize,1),bcltzr(nxsize,nysize,1),  &
    struzr(nxsize,nysize,1),strvzr(nxsize,nysize,1),  &
    strwzr(nxsize,nysize,1),strpzr(nxsize,nysize,1),  &
    strdzr(nxsize,nysize,1),strtzr(nxsize,nysize,1),  &
    strezr(nxsize,nysize,1),strgzr(nxsize,nysize,1), strrzr(nxsize,nysize,1),  &
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

REAL(KIND=dp) :: umax(nstore),umin(nstore),vmax(nstore),vmin(nstore),  &
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

