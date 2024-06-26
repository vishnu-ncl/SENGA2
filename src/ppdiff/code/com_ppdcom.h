C     PPCCOM--------------------------------------------------------------------

C     PARAMETERS
C     ==========
      DOUBLE PRECISION AVOGNO,BOLTZC
      PARAMETER(AVOGNO = 6.02214129D26, BOLTZC = 1.3806488D-23)

      DOUBLE PRECISION DIMTOL
      PARAMETER(DIMTOL = 1.0D-35)

      DOUBLE PRECISION TFITLO,TFITHI
      PARAMETER(TFITLO = 2.0D2, TFITHI = 3.0D3)

      DOUBLE PRECISION ZERO, ONE, TWO, THREE
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0)

      DOUBLE PRECISION HALF
      PARAMETER(HALF = 5.0D-1)

      DOUBLE PRECISION VSFACT,DIFACT
      PARAMETER(VSFACT = 5.0D0/1.6D1, DIFACT = 3.0D0/1.6D1)

      DOUBLE PRECISION TREFZD
      PARAMETER(TREFZD = 2.98D2)

      INTEGER NCINPT,NCREPT,NCSPEC,NCTHMO,NCDIFF,NCCOLL,NCOUTP
      PARAMETER(NCINPT=5, NCREPT=6,
     +          NCSPEC=1, NCTHMO=2, NCDIFF=3, NCCOLL=4, NCOUTP=7)

      INTEGER IENULL,IEINFO,IEWARN,IERROR,IESEVR,IEFATL
      PARAMETER(IENULL=0,
     +          IEINFO=1, IEWARN=2, IERROR=3, IESEVR=4, IEFATL=5)

      INTEGER ICHNUL,ICHSTR,ICHEND,ICHINT,ICHSNT,ICHRNO,ICHFNO
      PARAMETER(ICHNUL=0,
     +      ICHSTR=1, ICHEND=2, ICHINT=3, ICHSNT=4, ICHRNO=5, ICHFNO=6)

      INTEGER NTINMX,NCOFMX
      PARAMETER(NTINMX=2, NCOFMX=7)

      INTEGER NTRDMX,NDELMX
      PARAMETER(NTRDMX=40, NDELMX=10)

      INTEGER NFITMX,NFCOMX
      PARAMETER(NFITMX=40, NFCOMX=5)

      CHARACTER*50 BLANKS
      PARAMETER(
     +   BLANKS='                                                  ')


C     GLOBAL DATA
C     ===========
      DOUBLE PRECISION ACOFCP(NCOFMX,NTINMX,NSPCMX)
      DOUBLE PRECISION TINTLO(NTINMX,NSPCMX),TINTHI(NTINMX,NSPCMX)
      DOUBLE PRECISION OMEG11(NTRDMX,NDELMX),OMEG22(NTRDMX,NDELMX)
      DOUBLE PRECISION ASTAR(NTRDMX,NDELMX),BSTAR(NTRDMX,NDELMX)
      DOUBLE PRECISION CSTAR(NTRDMX,NDELMX)
      DOUBLE PRECISION TRED11(NTRDMX),TRED22(NTRDMX)
      DOUBLE PRECISION TREDAA(NTRDMX),TREDBB(NTRDMX),TREDCC(NTRDMX)
      DOUBLE PRECISION DELR11(NDELMX),DELR22(NDELMX)
      DOUBLE PRECISION DELRAA(NDELMX),DELRBB(NDELMX),DELRCC(NDELMX)
      DOUBLE PRECISION WMOLAR(NSPCMX),EMOMAS(NSPCMX),CLEWIS(NSPCMX)
      DOUBLE PRECISION EPSOKB(NSPCMX),SIGMAD(NSPCMX)
      DOUBLE PRECISION DIMOMU(NSPCMX),POLALF(NSPCMX),ZEDROT(NSPCMX)
      DOUBLE PRECISION RGUNIV
      DOUBLE PRECISION CLJP2K,CLJD2M,CDM2SI,CPL2M3
      DOUBLE PRECISION PI,TWOPI,ROOTPI,ROOTWO
      INTEGER NCOFCP(NTINMX,NSPCMX)
      INTEGER NTINT(NSPCMX)
      INTEGER IMGEOM(NSPCMX)
      INTEGER LENNAM(NSPCMX)
      INTEGER NUMSPC(NSPCMX)
      INTEGER NTRD11,NTRD22,NTRDAA,NTRDBB,NTRDCC
      INTEGER NDEL11,NDEL22,NDELAA,NDELBB,NDELCC
 
      CHARACTER*50 SPCNAM(NSPCMX)
      CHARACTER*80 FNSPEC,FNTHMO,FNDIFF,FNCOLL,FNOUTP

      LOGICAL FPOLAR(NSPCMX)

      COMMON/PPCCOM/ACOFCP,TINTLO,TINTHI,
     +              OMEG11,OMEG22,ASTAR,BSTAR,CSTAR,
     +              TRED11,TRED22,TREDAA,TREDBB,TREDCC,
     +              DELR11,DELR22,DELRAA,DELRBB,DELRCC,
     +              WMOLAR,EMOMAS,CLEWIS,
     +              EPSOKB,SIGMAD,DIMOMU,POLALF,ZEDROT,
     +              RGUNIV,
     +              CLJP2K,CLJD2M,CDM2SI,CPL2M3,
     +              PI,TWOPI,ROOTPI,ROOTWO,
     +              NCOFCP,NTINT,
     +              IMGEOM,
     +              NTRD11,NTRD22,NTRDAA,NTRDBB,NTRDCC,
     +              NDEL11,NDEL22,NDELAA,NDELBB,NDELCC,
     +              NUMSPC,LENNAM,SPCNAM,
     +              FNSPEC,FNTHMO,FNDIFF,FNCOLL,FNOUTP,
     +              FPOLAR

C     PPCCOM--------------------------------------------------------------------
