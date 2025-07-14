MODULE com_senga

!   PHYDIM-------------------------------------------------------------------

!   PHYSICAL DIMENSIONS OF ARRAYS
!   -----------------------------
!   NOTE: ALL ARRAY SIZES MUST BE CONSISTENT
!   NXSIZE MUST BE >= NXGLBL/NXPROC
!   NYSIZE MUST BE >= NYGLBL/NYPROC
!   NZSIZE MUST BE >= NZGLBL/NZPROC
!   WITH AN EXTRA ALLOWANCE FOR ANY REMAINDER
    use OPS_CONSTANTS
    implicit none

!   GLOBAL GRID SIZE
    integer(kind=4), parameter :: nxglbl=512, nyglbl=256, nzglbl=256
!   SET NGZMAX=MAX(NXGLBL,NYGLBL,NZGLBL)
    integer(kind=4), parameter :: ngzmax=nxglbl
    integer(kind=4), target :: nxglbl_ops=nxglbl, nyglbl_ops=nyglbl, nzglbl_ops=nzglbl

!   NUMBER OF PROCESSORS - DO NOT CHANGE - KEEP AS 1 ONLY
    integer(kind=4), parameter :: nxproc=1, nyproc=1, nzproc=1
!   SET NPRMAX=MAX(NXPROC,NYPROC,NZPROC)
    integer(kind=4), parameter :: nprmax=nxproc

!   LOCAL GRID SIZE - MENTION SAME VALUE AS GLOBAL GRID SIZE
    integer(kind=4), parameter :: nxsize=512, nysize=256, nzsize=256
!   SET NSZMAX=MAX(NXSIZE,NYSIZE,NZSIZE)
    integer(kind=4), parameter :: nszmax=nxsize

!   SIZE OF HALO
    integer(kind=4), parameter :: nhalox=5,nhaloy=5,nhaloz=5

!   SIZE OF PARALLEL TRANSFER ARRAY
!   NPARAY MUST BE >= MAX(NHALOX,NHALOY,NHALOZ)
!                          *MAX((NXSIZE+2*NHALOX)*(NYSIZE+2*NHALOY),
!                               (NXSIZE+2*NHALOX)*(NZSIZE+2*NHALOZ),
!                               (NYSIZE+2*NHALOY)*(NZSIZE+2*NHALOZ))
!   AND ALSO LARGE ENOUGH FOR PARALLEL BROADCAST OF CHEMICAL DATA
!   AND ALSO LARGE ENOUGH FOR PARALLEL TRANSFER OF INITIAL TURBULENCE DATA
    integer(kind=4), parameter :: nparay=8000000

!   LOCAL BIG ARRAY SIZE
    integer(kind=4), parameter :: nxbigl=1-nhalox, nxbigr=nxsize+nhalox,  &
                             nybigl=1-nhaloy, nybigr=nysize+nhaloy,  &
                             nzbigl=1-nhaloz, nzbigr=nzsize+nhaloz

!   INFLOW TURBULENCE GENERATOR BY SYNTHETIC DIGITAL FILTERING
    integer(kind=4), parameter :: nfy=60,nfz=60,lnx=30,lny=30,lnz=30
!   NFY AND NFZ ARE TWICE OF LNY AND LNZ

!   PHYDIM-------------------------------------------------------------------
!   IFTURB-------------------------------------------------------------------

!   DATA FOR INITIAL TURBULENCE FIELD

!   NUMBER OF PENCILS
    integer(kind=4), parameter :: npenmx=1

    real(kind=8) :: fftrow(2*ngzmax,3,npenmx), ftpart(2*nszmax,3,npenmx),  &
               fftinx(2*ngzmax)

!   IFTURB-------------------------------------------------------------------
!   CHEMIC-------------------------------------------------------------------

!   PARAMETERS
!   ==========
!   CHEMICAL RATE DATA
!   MAX NO OF THIRD BODIES
    integer(kind=4), parameter :: nbdymx=10
!   MAX NO OF LINDEMANN STEPS
    integer(kind=4), parameter :: nllmax=10

!   TRANSPORT COEFFICIENTS
    real(kind=8), parameter :: alamdc=0.0000258_8, tlamda=298.0_8, prantl=0.70_8
    real(kind=8), target :: prantl_ops=prantl

!   UNIVERSAL GAS CONSTANT
    real(kind=8), parameter :: rguniv=8314.20_8

!   CHEMICAL DATA
!   =============
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

    integer(kind=4) :: nsspec(nssmax,nstpmx)
    integer(kind=4) :: nrspec(nssmax,nstpmx),npspec(nssmax,nstpmx)
    integer(kind=4) :: nrcpec(nssmax,nstpmx),npcpec(nssmax,nstpmx)
    integer(kind=4) :: ncofcp(ntinmx,nspcmx)
    integer(kind=4) :: ncpoly(ntinmx,nspcmx),ncpom1(ntinmx,nspcmx)
    integer(kind=4) :: ncenth(ntinmx,nspcmx),ncenpy(ntinmx,nspcmx)
    integer(kind=4) :: nsslen(nstpmx)
    integer(kind=4) :: nrslen(nstpmx),npslen(nstpmx)
    integer(kind=4) :: nrclen(nstpmx),npclen(nstpmx)
    integer(kind=4) :: mblist(nstpmx)
    integer(kind=4) :: mglist(nstpmx)
    integer(kind=4) :: mllist(nstpmx)
    integer(kind=4) :: mtlist(nstpmx)
    integer(kind=4) :: mslist(nstpmx)
    integer(kind=4) :: ntint(nspcmx)
    integer(kind=4) :: nspec,nspm1,nstep,nbody,ngibb,nlind,ntroe,nsrif

!   MAX LENGTH OF SPECIES SYMBOL STRINGS
    integer(kind=4), parameter :: nspstr=10
!   SPECIES AND THIRD BODY STRINGS MUST BE DECLARED AS CHARACTER*NSPSTR
    character (len=10) :: spcsym(nspcmx)
    character (len=10) :: bdysym(nspcmx)

!   CHEMIC-------------------------------------------------------------------
!   DIFFUS-------------------------------------------------------------------

!   MIXTURE AVERAGED TRANSPORT

!   MASS FRACTION TOLERANCE
    real(kind=8), parameter :: dfctol = 0.0000000000010_8
    real(kind=8), target :: dfctol_ops=dfctol

!   MOLECULAR TRANSPORT DATA
    real(kind=8) :: diffco(ndcfmx,nspcmx,nspcmx)
    real(kind=8) :: tdrcco(ndcfmx,nspcmx,nspcmx)
    real(kind=8) :: wilko1(nspcmx,nspcmx),wilko2(nspcmx,nspcmx)
    real(kind=8) :: viscco(nvcfmx,nspcmx)
    real(kind=8) :: condco(nccfmx,nspcmx)
    real(kind=8) :: pdifgb,tdifgb
    real(kind=8) :: wmltdr

    integer(kind=4) :: ncodif,ncotdr,ncovis,ncocon
    integer(kind=4) :: ncodm1,ncotm1,ncovm1,ncocm1
    integer(kind=4) :: nfmavt,nfmixw,nfmixp,nfmsor,nfmduf

!   CONTROL FLAGS
    logical :: flmsor(nspcmx),flmduf(nspcmx),flmtdr(nspcmx)
    logical :: flmavt,flmixw,flmixp,flmixt

!   DIFFUS-------------------------------------------------------------------
!   RADIAT-------------------------------------------------------------------

!   RADIATION TREATMENT
    real(kind=8),parameter :: stefbo = 0.00000005670373210_8

    real(kind=8) :: akprad(ncfrmx,nspcmx)
    real(kind=8) :: trefrn

    integer(kind=4) :: nsprid(nspcmx),nkprad(nspcmx),nkprm1(nspcmx)
    integer(kind=4) :: nsprad,nfradn

    logical :: flradn

!   RADIAT-------------------------------------------------------------------
!   PARAMS-------------------------------------------------------------------

!   PHYSICAL AND NUMERICAL PARAMETERS

!   MASS FRACTION TOLERANCE
    real(kind=8), parameter :: ytoler = 1.0E-10

!   NUMBERS
    real(kind=8), parameter :: zero=0.0_8, one=1.0_8, two=2.0_8, three=3.0_8,  &
                          four=4.0_8, eight=8.0_8, &
                          half=one/two,  thrf=three*half, thrd=one/three,  &
                          tthd=two*thrd, fthd=two*tthd,   qrtr=one/four

    real(kind=8), target :: half_ops=half, tthd_ops=tthd, fthd_ops=fthd

!   PHYSICAL AND NUMERICAL DATA
    real(kind=8) :: yrin(nspcmx)
    real(kind=8) :: prin,trin,drin,urin,vrin,wrin,erin
    real(kind=8) :: xgdlen,ygdlen,zgdlen
    real(kind=8) :: etime,tstep,deltax,deltay,deltaz,pi,clnten
    integer(kind=4) :: itime,ntime,ntime1,ntime2,nstpsw,nsaved, inturb,inflam,inseed,  &
                  nxgreq,nygreq,nzgreq,nxpreq,nypreq,nzpreq,nspreq,  &
                  ntdump,ntrept,niters,itstat,ntstat,idflag,  &
                  ncont1,ncont2,ncont3,ncont4,ncont5,ncont6,  &
                  nccont,ncchem,ncdiff,ncradn,ncrept,ncstat,ncdmpi,ncdmpo, ndifmt,ndofmt
    character (len=60) :: fncont,fnchem,fndiff,fnradn,fnrept,fnstat,fndmpo(2)
    logical :: fladpt

!   PARAMS-------------------------------------------------------------------
!   DFDIFF-------------------------------------------------------------------

!   COEFFICIENTS AND END CONDITIONS FOR SPATIAL DIFFERENCING SCHEMES

!   DERIVATIVE BC STATUS
    integer(kind=4), parameter :: nobc=0, nbound=1, nperi=2
    integer(kind=4), target :: nbound_ops=nbound

!   COEFFICIENTS OF SPATIAL DIFFERENCING SCHEMES
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

!   SPATIAL DERIVATIVE END CONDITIONS
    integer(kind=4) :: nendxl,nendxr,nendyl,nendyr,nendzl,nendzr

!   DFDIFF-------------------------------------------------------------------
!   EMSTRT-------------------------------------------------------------------

!   DATA FOR MESH STRETCHING
    real(kind=8) :: gcmreg(nszmax),gcmstr(nszmax),  &
               dgdhat(nszmax),dgdhsq(nszmax),d2gdh2(nszmax)

!   FILCOM-------------------------------------------------------------------
!   RUNKUT-------------------------------------------------------------------

!   RUNGE-KUTTA PARAMETERS
    real(kind=8) :: rklhs(nrkmax),rkrhs(nrkmax)
    real(kind=8) :: rkerr(nrkmax),rktim(nrkmax)
    real(kind=8) :: ctmult,ctalph,ctbeta,ctgama
    real(kind=8) :: errtol,errlow,trmin,trmax,tsmin,tsmax
    real(kind=8) :: errold,errldr,btime
    real(kind=8) :: erdnrm,erunrm,ervnrm,erwnrm,erenrm
    real(kind=8) :: erynrm(nspcmx)
    integer(kind=4) :: irkstp,nrkstp,nrksm1
    integer(kind=4) :: inderr
    logical :: fupelc

!   RUNKUT-------------------------------------------------------------------
!   RKNORM-------------------------------------------------------------------

!   RUNGE-KUTTA SUBSTEP ERROR NORMS
    real(kind=8) :: erdrhs(nrkmax),erurhs(nrkmax),ervrhs(nrkmax),  &
               erwrhs(nrkmax),ererhs(nrkmax)
    real(kind=8) :: eryrhs(nspcmx,nrkmax)

!   PRINCIPAL VARIABLES (ALL STANDARD SIZE ARRAYS)
    real(kind=8), dimension(:,:,:), allocatable :: drun,urun,vrun,wrun,erun
    real(kind=8), dimension(:,:,:,:), allocatable :: yrun

!   WORKSPACE (STANDARD SIZE ARRAYS)
    real(kind=8), dimension(:,:,:,:), allocatable :: rrte

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!   WORKSPACE (BIGGER SIZE ARRAYS)
    real(kind=8), dimension(:,:,:), allocatable :: utmp,vtmp,wtmp,prun,trun

!   WORKSPACE (PARALLEL TRANSFER ARRAY)
    real(kind=8) :: parray(nparay)

!   WORKSP-------------------------------------------------------------------
!   NSBCCL-------------------------------------------------------------------

!   BC TYPE PARAMETERS
    integer(kind=4), parameter :: nsnull=0, nsperi=1, nsbci1=11,nsbci2=12,nsbci3=13,nsbci4=14,  &
                             nsbco1=21,nsbco2=22,nsbco3=23,nsbco4=24,  &
                             nsbcw1=31,nsbcw2=32,nsbcw3=33,nsbcw4=34

!   GLOBAL BC TYPE VARIABLES
    integer(kind=4) :: ngbcxl,ngbcxr,ngbcyl,ngbcyr,ngbczl,ngbczr

!   LOCAL BC TYPE VARIABLES
    integer(kind=4) :: nsbcxl,nsbcxr,nsbcyl,nsbcyr,nsbczl,nsbczr

!   BC FLAGS
    logical :: fxlcnv,fxlvsn,fxlvst,fxlcon,fxladb,fxldif, fxlcnw,fxladw,fxldfw,  &
               fxrcnv,fxrvsn,fxrvst,fxrcon,fxradb,fxrdif, fxrcnw,fxradw,fxrdfw,  &
               fylcnv,fylvsn,fylvst,fylcon,fyladb,fyldif, fylcnw,fyladw,fyldfw,  &
               fyrcnv,fyrvsn,fyrvst,fyrcon,fyradb,fyrdif, fyrcnw,fyradw,fyrdfw,  &
               fzlcnv,fzlvsn,fzlvst,fzlcon,fzladb,fzldif, fzlcnw,fzladw,fzldfw,  &
               fzrcnv,fzrvsn,fzrvst,fzrcon,fzradb,fzrdif, fzrcnw,fzradw,fzrdfw,  &
               fxltrb,fxrtrb,fyltrb,fyrtrb,fzltrb,fzrtrb

!   BC DATA
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
    integer(kind=4) :: nxlprm(nbcpri),nxrprm(nbcpri), nylprm(nbcpri),nyrprm(nbcpri),  &
                  nzlprm(nbcpri),nzrprm(nbcpri)
    integer(kind=4) :: istaxl,istoxl,istayl,istoyl,istazl,istozl,  &
                  istaxr,istoxr,istayr,istoyr,istazr,istozr,  &
                  kminxl,kminxr,kminyl,kminyr,kminzl,kminzr
    integer(kind=4) :: nctixl
    character (len=60) :: fntixl,fntcxl
    logical :: fllixl,fllixr,flliyl,flliyr,fllizl,fllizr,  &
               fltrxl,fltrxr,fltryl,fltryr,fltrzl,fltrzr

!   INFLOW VELOCITY FIELD DATA
!   REQUIRED ONLY FOR TURBULENT INFLOW
    real(kind=8), dimension(:,:,:), allocatable :: ufxl,vfxl,wfxl

!   WALL BC DIFFERENCING DATA
    real(kind=8) :: acbcxl(ncbcsz),acbcxr(ncbcsz),  &
               acbcyl(ncbcsz),acbcyr(ncbcsz), acbczl(ncbcsz),acbczr(ncbcsz)

!   X-DIRECTION LEFT-HAND END
    integer(kind=4) :: intran
    real(kind=8) :: passer
    real(kind=8), dimension(:,:,:), allocatable :: struxl,strvxl,strwxl,dudtxl,dvdtxl,dwdtxl

!   NSBCCL-------------------------------------------------------------------
!   DOMDEC-------------------------------------------------------------------

!   PARALLEL DOMAIN DECOMPOSITION DATA
    integer(kind=4) :: npmapx(0:nxproc),npmapy(0:nyproc),npmapz(0:nzproc),  &
                  nprocx(0:nxproc),nprocy(0:nyproc),nprocz(0:nzproc),  &
                  nxnode,nynode,nznode, nxprm1,nyprm1,nzprm1,  &
                  iproc,nproc,ixproc,iyproc,izproc

!   DOMDEC-------------------------------------------------------------------
!   STATIS-------------------------------------------------------------------

!   DATA FOR ON-LINE STATISTICS

    integer(kind=4), parameter :: nstore=32

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
    integer(kind=4) :: itstim(nstore)

!   STATIS-------------------------------------------------------------------

!   FY - FOR NON REFLECTING INFLOW
    integer(kind=4), parameter :: nbcprc=9
    real(kind=8) :: rxlprc(nbcprc)

!   -------------------------------------------------------------------------

END MODULE com_senga

