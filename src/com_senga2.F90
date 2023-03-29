MODULE com_senga

!     PHYDIM-------------------------------------------------------------------
 
!     PHYSICAL DIMENSIONS OF ARRAYS
!     -----------------------------
!     NOTE: ALL ARRAY SIZES MUST BE CONSISTENT
!     NXSIZE MUST BE >= NXGLBL/NXPROC
!     NYSIZE MUST BE >= NYGLBL/NYPROC
!     NZSIZE MUST BE >= NZGLBL/NZPROC
!     WITH AN EXTRA ALLOWANCE FOR ANY REMAINDER
    
    use OPS_CONSTANTS

!     GLOBAL GRID SIZE
integer :: nxglbl,nyglbl,nzglbl
PARAMETER(nxglbl=16, nyglbl=16, nzglbl=16)
integer :: ngzmax
!     SET NGZMAX=MAX(NXGLBL,NYGLBL,NZGLBL)
PARAMETER(ngzmax=nxglbl)

!     NUMBER OF PROCESSORS
integer :: nxproc,nyproc,nzproc
PARAMETER(nxproc=1, nyproc=1, nzproc=1)
integer :: nprmax
!     SET NPRMAX=MAX(NXPROC,NYPROC,NZPROC)
PARAMETER(nprmax=nxproc)

!     LOCAL GRID SIZE
integer :: nxsize,nysize,nzsize
PARAMETER(nxsize=16, nysize=16, nzsize=16)
integer :: nszmax
!     SET NSZMAX=MAX(NXSIZE,NYSIZE,NZSIZE)
PARAMETER(nszmax=nxsize)

!     SIZE OF HALO
integer :: nhalox,nhaloy,nhaloz
PARAMETER(nhalox=5,nhaloy=5,nhaloz=5)

!     SIZE OF PARALLEL TRANSFER ARRAY
!     NPARAY MUST BE >= MAX(NHALOX,NHALOY,NHALOZ)
!                            *MAX((NXSIZE+2*NHALOX)*(NYSIZE+2*NHALOY),
!                                 (NXSIZE+2*NHALOX)*(NZSIZE+2*NHALOZ),
!                                 (NYSIZE+2*NHALOY)*(NZSIZE+2*NHALOZ))
!     AND ALSO LARGE ENOUGH FOR PARALLEL BROADCAST OF CHEMICAL DATA
!     AND ALSO LARGE ENOUGH FOR PARALLEL TRANSFER OF INITIAL TURBULENCE DATA
integer :: nparay
PARAMETER(nparay=8000000)

!     LOCAL BIG ARRAY SIZE
integer :: nxbigl,nxbigr,nybigl,nybigr,nzbigl,nzbigr
PARAMETER(nxbigl=1-nhalox, nxbigr=nxsize+nhalox,  &
    nybigl=1-nhaloy, nybigr=nysize+nhaloy, nzbigl=1-nhaloz, nzbigr=nzsize+nhaloz)

!     PHYDIM-------------------------------------------------------------------
!     IFTURB-------------------------------------------------------------------

!     DATA FOR INITIAL TURBULENCE FIELD

!     NUMBER OF PENCILS
integer :: npenmx
PARAMETER(npenmx=1)

real(kind=8) :: fftrow(2*ngzmax,3,npenmx), ftpart(2*nszmax,3,npenmx),  &
    fftinx(2*ngzmax)

COMMON/ifturb/fftrow,ftpart,fftinx

!     IFTURB-------------------------------------------------------------------
!     CHEMIC-------------------------------------------------------------------

!     PARAMETERS
!     ==========
!     MAX NO OF SPECIES, NO OF STEPS
integer :: nspcmx,nstpmx
PARAMETER(nspcmx=2, nstpmx=1)

!     THERMODYNAMIC DATA
!     MAX NO OF TEMPERATURE INTERVALS, THERMO POLYNOMIAL COEFFICIENTS
integer :: ntinmx,ncofmx
PARAMETER(ntinmx=2, ncofmx=7)
!     MAX NO OF TEMPERATURE COEFFICIENTS, DITTO MINUS ONE
integer :: nctmax,nctmm1
PARAMETER(nctmax=5, nctmm1=nctmax-1)

!     TEMPERATURE INTERVAL INDEXING
!     NTBASE = NUMBER BASE FOR INDEXING:
!     MUST BE A POWER OF TWO >= MAX NO OF TEMPERATURE INTERVALS PER SPECIES
!     NSPIMX = MAX NO OF SPECIES STORED PER SINGLE (32-BIT) SIGNED integer:
!     MUST BE SET EQUAL TO 31 DIV LOG2(NTBASE)
!     NINTMX = NO OF integerS REQUIRED PER SPATIAL POINT:
!     MUST BE SET EQUAL TO (1 + NSPCMX DIV NSPIMX)
integer :: nspimx,ntbase,nintmx
PARAMETER(nspimx=15, ntbase=4, nintmx=2)

!     CHEMICAL RATE DATA
!     MAX NO OF THIRD BODIES
integer :: nbdymx
PARAMETER(nbdymx=10)
!     MAX SIZE OF STEP SPECIES-LIST, STEP REACTANT-LIST
integer :: nssmax,nrsmax
PARAMETER(nssmax=10, nrsmax=10)
!     MAX NO OF LINDEMANN STEPS
integer :: nllmax
PARAMETER(nllmax=10)

!     UNIVERSAL GAS CONSTANT
real(kind=8) :: rguniv
PARAMETER(rguniv=8314.20_8)

!     CHEMICAL DATA
!     =============
real(kind=8) :: amolcp(ncofmx,ntinmx,nspcmx)
real(kind=8) :: amolgb(ncofmx,ntinmx,nspcmx)
real(kind=8) :: amascp(ncofmx,ntinmx,nspcmx)
real(kind=8) :: amasct(ncofmx,ntinmx,nspcmx)
real(kind=8) :: amasch(ncofmx,ntinmx,nspcmx)
real(kind=8) :: amascs(ncofmx,ntinmx,nspcmx)
real(kind=8) :: tintlo(ntinmx,nspcmx),tinthi(ntinmx,nspcmx)
real(kind=8) :: diffmu(nssmax,nstpmx),diffmw(nssmax,nstpmx)
real(kind=8) :: crspec(nssmax,nstpmx),cpspec(nssmax,nstpmx)
real(kind=8) :: effy3b(nspcmx,nbdymx)
real(kind=8) :: rparam(3,nstpmx)
real(kind=8) :: rclind(4,nllmax)
real(kind=8) :: rctroe(12,nllmax)
real(kind=8) :: rcsrif(8,nllmax)
real(kind=8) :: wmolar(nspcmx),ovwmol(nspcmx),rgspec(nspcmx)
real(kind=8) :: clewis(nspcmx),olewis(nspcmx)
real(kind=8) :: prefgb

integer :: nsspec(nssmax,nstpmx)
integer :: nrspec(nssmax,nstpmx),npspec(nssmax,nstpmx)
integer :: nrcpec(nssmax,nstpmx),npcpec(nssmax,nstpmx)
integer :: ncofcp(ntinmx,nspcmx)
integer :: ncpoly(ntinmx,nspcmx),ncpom1(ntinmx,nspcmx)
integer :: ncenth(ntinmx,nspcmx),ncenpy(ntinmx,nspcmx)
integer :: nsslen(nstpmx)
integer :: nrslen(nstpmx),npslen(nstpmx)
integer :: nrclen(nstpmx),npclen(nstpmx)
integer :: mblist(nstpmx)
integer :: mglist(nstpmx)
integer :: mllist(nstpmx)
integer :: mtlist(nstpmx)
integer :: mslist(nstpmx)
integer :: ntint(nspcmx)
integer :: nspec,nspm1,nstep,nbody,ngibb,nlind,ntroe,nsrif

!     MAX LENGTH OF SPECIES SYMBOL STRINGS
integer :: nspstr
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
real(kind=8) :: dfctol
PARAMETER(dfctol = 0.0000000000010_8)

!     MAX NUMBERS OF POLYNOMIAL COEFFICIENTS
integer :: ndcfmx,nvcfmx,nccfmx
PARAMETER(ndcfmx = 4, nvcfmx = 4, nccfmx = 4)

!     MOLECULAR TRANSPORT DATA
real(kind=8) :: diffco(ndcfmx,nspcmx,nspcmx)
real(kind=8) :: tdrcco(ndcfmx,nspcmx,nspcmx)
real(kind=8) :: wilko1(nspcmx,nspcmx),wilko2(nspcmx,nspcmx)
real(kind=8) :: viscco(nvcfmx,nspcmx)
real(kind=8) :: condco(nccfmx,nspcmx)
real(kind=8) :: pdifgb,tdifgb
real(kind=8) :: wmltdr

integer :: ncodif,ncotdr,ncovis,ncocon
integer :: ncodm1,ncotm1,ncovm1,ncocm1
integer :: nfmavt,nfmixw,nfmixp,nfmsor,nfmduf

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

real(kind=8) :: stefbo
PARAMETER(stefbo = 0.00000005670373210_8)

integer :: ncfrmx
PARAMETER(ncfrmx = 6)

real(kind=8) :: akprad(ncfrmx,nspcmx)
real(kind=8) :: trefrn

integer :: nsprid(nspcmx),nkprad(nspcmx),nkprm1(nspcmx)
integer :: nsprad,nfradn

LOGICAL :: flradn

COMMON/radiat/akprad,trefrn,nsprid,nkprad,nkprm1,nsprad,  &
    nfradn, flradn

!     RADIAT-------------------------------------------------------------------
!     PARAMS-------------------------------------------------------------------

!     PHYSICAL AND NUMERICAL PARAMETERS

!     MASS FRACTION TOLERANCE
real(kind=8) :: ytoler
PARAMETER(ytoler = 1.0E-10)

!     NUMBERS
real(kind=8) :: zero,one,two,three,four,eight,  &
    half,thrf,thrd,tthd,fthd,qrtr
PARAMETER(zero=0.0_8, one=1.0_8, two=2.0_8, three=3.0_8,  &
    four=4.0_8, eight=8.0_8)
PARAMETER(half=one/two,  thrf=three*half, thrd=one/three,  &
    tthd=two*thrd, fthd=two*tthd,   qrtr=one/four)

!     LOCATION FOR REFERENCE PRESSURE VALUE
integer :: ipref,jpref,kpref
PARAMETER(ipref=1,jpref=1,kpref=1)

!     PHYSICAL AND NUMERICAL DATA
real(kind=8) :: yrin(nspcmx)
real(kind=8) :: prin,trin,drin,urin,vrin,wrin,erin
real(kind=8) :: xgdlen,ygdlen,zgdlen
real(kind=8) :: etime,tstep,deltax,deltay,deltaz,pi,clnten
integer :: itime,ntime,ntime1,ntime2,nstpsw,nsaved, inturb,inflam,inseed,  &
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
integer :: nobc,nbound,nperi
PARAMETER(nobc=0,nbound=1,nperi=2)

!     COEFFICIENTS OF SPATIAL DIFFERENCING SCHEMES
real(kind=8) :: acoeff,bcoeff,ccoeff,dcoeff,ecoeff,  &
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
    acofc1,bcofc1,ccofc1,dcofc1, acofc2,bcofc2,ccofc2,dcofc2

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
    acofc1,bcofc1,ccofc1,dcofc1, acofc2,bcofc2,ccofc2,dcofc2, &
    nendxl,nendxr,nendyl,nendyr,nendzl,nendzr

!     DFDIFF-------------------------------------------------------------------
!     EMSTRT-------------------------------------------------------------------

!     DATA FOR MESH STRETCHING
real(kind=8) :: gcmreg(nszmax),gcmstr(nszmax),  &
    dgdhat(nszmax),dgdhsq(nszmax),d2gdh2(nszmax)

COMMON/emstrt/gcmreg,gcmstr,dgdhat,dgdhsq,d2gdh2

!     FILCOM-------------------------------------------------------------------
!     RUNKUT-------------------------------------------------------------------

!     RUNGE-KUTTA PARAMETERS
integer :: nrkmax
PARAMETER(nrkmax=5)

real(kind=8) :: rklhs(nrkmax),rkrhs(nrkmax)
real(kind=8) :: rkerr(nrkmax),rktim(nrkmax)
real(kind=8) :: ctmult,ctalph,ctbeta,ctgama
real(kind=8) :: errtol,errlow,trmin,trmax,tsmin,tsmax
real(kind=8) :: errold,errldr,btime
real(kind=8) :: erdnrm,erunrm,ervnrm,erwnrm,erenrm
real(kind=8) :: erynrm(nspcmx)
integer :: irkstp,nrkstp,nrksm1
integer :: inderr
LOGICAL :: fupelc

COMMON/runkut/rklhs,rkrhs,rkerr,rktim, ctmult,ctalph,ctbeta,ctgama,  &
    errtol,errlow,trmin,trmax,tsmin,tsmax, errold,errldr,btime,  &
    erdnrm,erunrm,ervnrm,erwnrm,erenrm,erynrm, irkstp,nrkstp,nrksm1,inderr,fupelc

!     RUNKUT-------------------------------------------------------------------
!     RKNORM-------------------------------------------------------------------

!     RUNGE-KUTTA SUBSTEP ERROR NORMS
real(kind=8) :: erdrhs(nrkmax),erurhs(nrkmax),ervrhs(nrkmax),  &
    erwrhs(nrkmax),ererhs(nrkmax)
real(kind=8) :: eryrhs(nspcmx,nrkmax)

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
real(kind=8) :: parray(nparay)

COMMON/workpt/parray

!     WORKSP-------------------------------------------------------------------
!     NSBCCL-------------------------------------------------------------------

!     BC TYPE PARAMETERS
integer :: nsnull,nsperi, nsbci1,nsbci2,nsbci3,nsbci4,  &
    nsbco1,nsbco2,nsbco3,nsbco4, nsbcw1,nsbcw2,nsbcw3,nsbcw4
PARAMETER(nsnull=0, nsperi=1, nsbci1=11,nsbci2=12,nsbci3=13,nsbci4=14,  &
    nsbco1=21,nsbco2=22,nsbco3=23,nsbco4=24,  &
    nsbcw1=31,nsbcw2=32,nsbcw3=33,nsbcw4=34)

!     BC NUMBER OF PARAMETERS
integer :: nbcpri,nbcprr
PARAMETER(nbcpri=4, nbcprr=4)

!     GLOBAL BC TYPE VARIABLES
integer :: ngbcxl,ngbcxr,ngbcyl,ngbcyr,ngbczl,ngbczr
COMMON/nglbct/ngbcxl,ngbcxr,ngbcyl,ngbcyr,ngbczl,ngbczr

!     LOCAL BC TYPE VARIABLES
integer :: nsbcxl,nsbcxr,nsbcyl,nsbcyr,nsbczl,nsbczr
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
real(kind=8) :: rxlprm(nbcprr),rxrprm(nbcprr),  &
    rylprm(nbcprr),ryrprm(nbcprr), rzlprm(nbcprr),rzrprm(nbcprr)
real(kind=8) :: cobcxl,cobcxr,cobcyl,cobcyr,cobczl,cobczr,  &
    pinfxl,pinfxr,pinfyl,pinfyr,pinfzl,pinfzr
real(kind=8) :: slocxl,slocxr,slocyl,slocyr,sloczl,sloczr,  &
    elocxl,elocxr,elocyl,elocyr,eloczl,eloczr,  &
    bvelxl,bvelxr,bvelyl,bvelyr,bvelzl,bvelzr,  &
    svelxl,svelxr,svelyl,svelyr,svelzl,svelzr,  &
    scauxl,scauxr,scauyl,scauyr,scauzl,scauzr,  &
    scduxl,scduxr,scduyl,scduyr,scduzl,scduzr, tpovxg,tpovyg,tpovzg
integer :: nxlprm(nbcpri),nxrprm(nbcpri), nylprm(nbcpri),nyrprm(nbcpri),  &
    nzlprm(nbcpri),nzrprm(nbcpri)
integer :: istaxl,istoxl,istayl,istoyl,istazl,istozl,  &
    istaxr,istoxr,istayr,istoyr,istazr,istozr,  &
    kminxl,kminxr,kminyl,kminyr,kminzl,kminzr
integer :: nctixl
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
integer :: ncbcsz
PARAMETER(ncbcsz=5)
real(kind=8) :: acbcxl(ncbcsz),acbcxr(ncbcsz),  &
    acbcyl(ncbcsz),acbcyr(ncbcsz), acbczl(ncbcsz),acbczr(ncbcsz)
COMMON/bcdifw/acbcxl,acbcxr,acbcyl,acbcyr,acbczl,acbczr

!     X-DIRECTION LEFT-HAND END
real(kind=8), dimension(:,:,:), allocatable :: struxl,strvxl,strwxl,dudtxl,dvdtxl,dwdtxl

!     NSBCCL-------------------------------------------------------------------
!     DOMDEC-------------------------------------------------------------------

!     PARALLEL DOMAIN DECOMPOSITION DATA
integer :: npmapx(0:nxproc),npmapy(0:nyproc),npmapz(0:nzproc),  &
    nprocx(0:nxproc),nprocy(0:nyproc),nprocz(0:nzproc),  &
    nxnode,nynode,nznode, nxprm1,nyprm1,nzprm1,  &
    iproc,nproc,ixproc,iyproc,izproc

COMMON/domdec/npmapx,npmapy,npmapz,nprocx,nprocy,nprocz,  &
    nxnode,nynode,nznode, nxprm1,nyprm1,nzprm1,  &
    iproc,nproc,ixproc,iyproc,izproc

!     DOMDEC-------------------------------------------------------------------
!     STATIS-------------------------------------------------------------------

!     DATA FOR ON-LINE STATISTICS

integer :: nstore
PARAMETER(nstore=32)

real(kind=8) :: umax(nstore),umin(nstore),vmax(nstore),vmin(nstore),  &
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
integer :: itstim(nstore)

COMMON/statis/umax,umin,vmax,vmin,wmax,wmin, ubar,vbar,wbar,uvar,vvar,wvar,  &
    tke,dissip,dissom,disot1,disot2,disot3, residx,residy,residz,  &
    xlenth,ylenth,zlenth,tkecor, uveloc,vveloc,wveloc,  &
    skewdx,skewdy,skewdz, tayscx,tayscy,tayscz,  &
    etalen,taueta,veleta, glorat,fsturb,fmassb,  &
    tint,uprm,ubarg,vbarg,wbarg, uvarg,vvarg,wvarg,dissipg,  &
    itstim

!     STATIS-------------------------------------------------------------------
END MODULE com_senga

