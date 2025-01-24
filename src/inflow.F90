SUBROUTINE inflow
 
! Code converted using TO_F90 by Alan Miller
! Date: 2024-11-08  Time: 05:07:37

!     *************************************************************************

!     INFLOW
!     ======

!     AUTHOR
!     ------
!     M.KLEIN  --   UNIVERSITÄT DER BUNDESWEHR MÜNCHEN

!     CHANGE RECORD
!     -------------
!     27-MAR-2015:  CREATED
!     10-APR-2019:  CHANGED FOR SENGA2 (M.A.)

!     DESCRIPTION
!     -----------

!     *************************************************************************
!      IMPLICIT NONE


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!INCLUDE 'mpif.h'
!INCLUDE 'com_senga2.h'
!     -------------------------------------------------------------------------


!     LOCAL DATA
!     ==========
REAL :: ran1        ! Random number Generator
DOUBLE PRECISION :: lenx,leny,lenz,theta,umean,delx,phi,tstp,urms

DOUBLE PRECISION :: norm,ay(-nfy:nfy),az(-nfz:nfz)  ! filter coefficients

LOGICAL :: periy,periz!,TRADIT
PARAMETER (periy=.true.,periz=.false.)!,TRADIT=.TRUE.)

INTEGER :: yinfl,yinfr,zinfl,zinfr,nylcl,nzlcl      ! rectangular inflow area
PARAMETER (yinfl=1,yinfr=nyglbl, zinfl=1,zinfr=nzglbl)

DOUBLE PRECISION :: ufilt(yinfl:nysize,zinfl:nzsize),  &
    vfilt(yinfl:nysize,zinfl:nzsize),  &
    wfilt(yinfl:nysize,zinfl:nzsize),  &
    yfilt(yinfl:nysize,zinfl:nzsize,nspec)

DOUBLE PRECISION :: ufold(yinfl:nysize,zinfl-nfz:nzsize+nfz),  &
    vfold(yinfl:nysize,zinfl-nfz:nzsize+nfz),  &
    wfold(yinfl:nysize,zinfl-nfz:nzsize+nfz),  &
    yfold(yinfl:nysize,zinfl-nfz:nzsize+nfz,nspec)

DOUBLE PRECISION :: urand(yinfl-nfy:yinfr+nfy,zinfl-nfz:zinfr+nfz),  &
    vrand(yinfl-nfy:yinfr+nfy,zinfl-nfz:zinfr+nfz),  &
    wrand(yinfl-nfy:yinfr+nfy,zinfl-nfz:zinfr+nfz),  &
    yrand(yinfl-nfy:yinfr+nfy,zinfl-nfz:zinfr+nfz,nspec)

DOUBLE PRECISION :: umeanglbl,sumumeanglbl

INTEGER :: ic,jc,kc,jc2,kc2,yp,zp
INTEGER :: jfiltstart,kfiltstart
INTEGER :: tspc(2),nspc,ispec

!     VM: CALCULATING USTEAD
DOUBLE PRECISION :: rad,denom,sumdenom
DOUBLE PRECISION :: zdist,ydist
DOUBLE PRECISION :: yref(nspec)
INTEGER :: kg,jg
INTEGER :: yoffset,zoffset,ymidpnt,zmidpnt
!     VM: DEBUG
INTEGER :: snap,icproc
CHARACTER (LEN=40) :: fname
CHARACTER (LEN=4) :: psnap

!     LNX, LNY, LNZ: LENGTH SCALES IN TERMS OF DELX
!     THETA: INDUCED TIME SCALE
!     PHI: REQUIRED WEIGHTING PARAMETER
!     NFY,NFZ: FILTER SIZE: NF=2*LN
delx = xgdlen/REAL(nxglbl-1)
umean = rxlprm(1)
tstp = tstep
urms = rxlprm(2)
yrms = 0.050D0
lenx = DBLE(lnx)
leny = DBLE(lny)
lenz = DBLE(lnz)
theta=lenx*delx/umean
phi=1.0-tstp/theta

zp=INT(iproc/(nxproc*nyproc))
yp=INT((iproc-zp*nyproc*nxproc)/nxproc)

jfiltstart=yp*nysize
kfiltstart=zp*nzsize

tspc = (/1,2/)
nspc = 2
DO ispec=1,nspec
  yref(ispec)=0.0D0
END DO
yref(1)=1.0D0
yref(2)=0.233D0

!     ------------------------------------------------------------------

!C    URMS=MAX(5.0D0,3.0D0+FLOAT(ITIME)/10000.D0*1.0D0)
!     CHECK RAMP UP OF VELOCITY
!      IF (NXLPRM(2).EQ.1)THEN
!         TO BE DONE
!      END IF
pi = four*ATAN(1.0D0)
!     CALCULATE FILTER COEFFICIENTS
norm=0.0
DO jc=-nfy,nfy
  ay(jc)=EXP(-pi*jc*jc/(2*lny*lny))
  norm=norm+ay(jc)**2
END DO
norm=SQRT(norm)
DO jc=-nfy,nfy
  ay(jc)=ay(jc)/norm
END DO

norm=0.0
DO kc=-nfz,nfz
  az(kc)=EXP(-pi*kc*kc/(2*lnz*lnz))
  norm=norm+az(kc)**2
END DO
norm=SQRT(norm)
DO kc=-nfz,nfz
  az(kc)=az(kc)/norm
END DO

intran = itime

!     INITIALIZE RANDOM ARRAYS
DO jc=yinfl-nfy,yinfr+nfy
  DO kc=zinfl-nfz,zinfr+nfz
    urand(jc,kc)=ran1(intran)
    vrand(jc,kc)=ran1(intran)
    wrand(jc,kc)=ran1(intran)
  END DO
END DO
IF(nxlprm(2)==1)THEN
  DO jc=yinfl-nfy,yinfr+nfy
    DO kc=zinfl-nfz,zinfr+nfz
      yrand(jc,kc,:)=0.0D0
      DO ispc = 1,nspc
        ispec = tspc(ispc)
        yrand(jc,kc,ispec)=ran1(intran)
      END DO
    END DO
  END DO
END IF

IF (periy) THEN   ! overwrite in a periodic manner
  DO jc=-nfy+1,0
    DO kc=zinfl-nfz,zinfr+nfz
      urand(jc,kc)=urand(jc+yinfr,kc)
      vrand(jc,kc)=vrand(jc+yinfr,kc)
      wrand(jc,kc)=wrand(jc+yinfr,kc)
    END DO
  END DO
  DO jc=yinfr+1,yinfr+nfy
    DO kc=zinfl-nfz,zinfr+nfz
      urand(jc,kc)=urand(jc-yinfr,kc)
      vrand(jc,kc)=vrand(jc-yinfr,kc)
      wrand(jc,kc)=wrand(jc-yinfr,kc)
    END DO
  END DO
  IF(nxlprm(2)==1)THEN
    DO jc=-nfy+1,0
      DO kc=zinfl-nfz,zinfr+nfz
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yrand(jc,kc,ispec)=yrand(jc+yinfr,kc,ispec)
        END DO
      END DO
    END DO
    DO jc=yinfr+1,yinfr+nfy
      DO kc=zinfl-nfz,zinfr+nfz
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yrand(jc,kc,ispec)=yrand(jc-yinfr,kc,ispec)
        END DO
      END DO
    END DO
  END IF
END IF

IF (periz) THEN   ! overwrite in a periodic manner
  DO kc=-nfz+1,0
    DO jc=yinfl-nfy,yinfr+nfy
      urand(jc,kc)=urand(jc,kc+zinfr)
      vrand(jc,kc)=vrand(jc,kc+zinfr)
      wrand(jc,kc)=wrand(jc,kc+zinfr)
    END DO
  END DO
  DO kc=zinfr+1,zinfr+nfz
    DO jc=yinfl-nfy,yinfr+nfy
      urand(jc,kc)=urand(jc,kc-zinfr)
      vrand(jc,kc)=vrand(jc,kc-zinfr)
      wrand(jc,kc)=wrand(jc,kc-zinfr)
    END DO
  END DO
  IF(nxlprm(2)==1)THEN
    DO kc=-nfz+1,0
      DO jc=yinfl-nfy,yinfr+nfy
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yrand(jc,kc,ispec)=yrand(jc,kc+zinfr,ispec)
        END DO
      END DO
    END DO
    DO kc=zinfr+1,zinfr+nfz
      DO jc=yinfl-nfy,yinfr+nfy
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yrand(jc,kc,ispec)=yrand(jc,kc-zinfr,ispec)
        END DO
      END DO
    END DO
  END IF
END IF

ufilt(:,:)=0.0
vfilt(:,:)=0.0
wfilt(:,:)=0.0
yfilt(:,:,:)=0.0
ufold(:,:)=0.0
vfold(:,:)=0.0
wfold(:,:)=0.0
yfold(:,:,:)=0.0


!     What is TRADIT?
DO jc=jstal,jstol
  DO kc=kstal-nfz,kstol+nfz
    DO jc2=-nfy,nfy
      
      ufold(jc,kc)=ufold(jc,kc)+  &
          urand(jc+jfiltstart+jc2,kc+kfiltstart)*ay(jc2)
      vfold(jc,kc)=vfold(jc,kc)+  &
          vrand(jc+jfiltstart+jc2,kc+kfiltstart)*ay(jc2)
      wfold(jc,kc)=wfold(jc,kc)+  &
          wrand(jc+jfiltstart+jc2,kc+kfiltstart)*ay(jc2)
    END DO
  END DO
END DO
IF(nxlprm(2)==1)THEN
  DO jc=jstal,jstol
    DO kc=kstal-nfz,kstol+nfz
      DO jc2=-nfy,nfy
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yfold(jc,kc,ispec)=yfold(jc,kc,ispec)+  &
              yrand(jc+jfiltstart+jc2,kc+kfiltstart,ispec)*ay(jc2)
        END DO
      END DO
    END DO
  END DO
END IF

DO jc=jstal,jstol
  DO kc=kstal,kstol
    DO kc2=-nfz,nfz
      ufilt(jc,kc)=ufilt(jc,kc)+ ufold(jc,kc+kc2)*az(kc2)
      vfilt(jc,kc)=vfilt(jc,kc)+ vfold(jc,kc+kc2)*az(kc2)
      wfilt(jc,kc)=wfilt(jc,kc)+ wfold(jc,kc+kc2)*az(kc2)
    END DO
  END DO
END DO
IF(nxlprm(2)==1)THEN
  DO jc=jstal,jstol
    DO kc=kstal,kstol
      DO kc2=-nfz,nfz
        DO ispc = 1,nspc
          ispec = tspc(ispc)
          yfilt(jc,kc,ispec)=yfilt(jc,kc,ispec)+  &
              yfold(jc,kc+kc2,ispec)*az(kc2)
        END DO
      END DO
    END DO
  END DO
END IF

!     Turbulence intensity
ufilt=ufilt*SQRT(urms**two*(one-phi**two))
vfilt=vfilt*SQRT(urms**two*(one-phi**two))
wfilt=wfilt*SQRT(urms**two*(one-phi**two))
IF(nxlprm(2)==1)THEN
  DO jc=jstal,jstol
    DO kc=kstal,kstol
      DO ispc = 1,nspc
        ispec = tspc(ispc)
        yfilt(jc,kc,ispec)=yfilt(jc,kc,ispec)*  &
            SQRT((yrms*yref(ispec))**two*(one-phi**two))
      END DO
    END DO
  END DO
END IF

uinf1=uinf2
vinf1=vinf2
winf1=winf2
IF(nxlprm(2)==1)THEN
  DO jc=jstal,jstol
    DO kc=kstal,kstol
      DO ispc = 1,nspc
        ispec = tspc(ispc)
        yinf1(jc,kc,ispec)=yinf2(jc,kc,ispec)
      END DO
    END DO
  END DO
END IF


DO jc=jstal,jstol
  DO kc=kstal,kstol
    uinf2(jc,kc)=ufilt(jc,kc)+phi*uinf1(jc,kc)
    vinf2(jc,kc)=vfilt(jc,kc)+phi*vinf1(jc,kc)
    winf2(jc,kc)=wfilt(jc,kc)+phi*winf1(jc,kc)
  END DO
END DO
IF(nxlprm(2)==1)THEN
  DO jc=jstal,jstol
    DO kc=kstal,kstol
      DO ispc = 1,nspc
        ispec = tspc(ispc)
        yinf2(jc,kc,ispec)=yfilt(jc,kc,ispec)+ phi*yinf1(jc,kc,ispec)
      END DO
    END DO
  END DO
END IF

!     VM: CALCULATING USTEAD

yoffset = 0
DO icproc = 0, iyproc-1
  yoffset = yoffset + npmapy(icproc)
END DO
zoffset = 0
ymidpnt = nyglbl/2
zmidpnt = 1
IF (nyglbl==1) THEN
  deltay = 0.0D0
ELSE
  deltay = ygdlen/DBLE(nyglbl-1)
END IF
IF (nzglbl==1) THEN
  deltaz = 0.0D0
ELSE
  deltaz = zgdlen/DBLE(nzglbl-1)
END IF

DO kc = kstal,kstol
  DO jc = jstal,jstol
    kg = zoffset+kc
    jg = yoffset+jc
    ydist = deltay*DBLE((jg-ymidpnt))
    zdist = deltaz*DBLE((kg-zmidpnt))
    rad=SQRT(ydist**2+zdist**2)
    ustead(jc,kc)=1.0D0!0.5D0*(1.0D0-TANH((RAD-RXLPRM(3)*YGDLEN)
!     +                  /(0.001*YGDLEN)))
  END DO
END DO

DO kc=kstal,kstol
  DO jc=jstal,jstol
    uinf2(jc,kc)=uinf2(jc,kc)*ustead(jc,kc)
    vinf2(jc,kc)=vinf2(jc,kc)*ustead(jc,kc)
    winf2(jc,kc)=0.0D0!WINF2(JC,KC)*USTEAD(JC,KC)
  END DO
END DO
IF(nxlprm(2)==1)THEN
  DO kc=kstal,kstol
    DO jc=jstal,jstol
      DO ispc = 1,nspc
        ispec = tspc(ispc)
        yinf2(jc,kc,ispec)=yinf2(jc,kc,ispec)*ustead(jc,kc)
      END DO
    END DO
  END DO
END IF

! Calculate the average speed in the entire domain, which is then
! subtracted from u_new
sumumeanglbl=zero
sumdenom=zero
denom = zero
umeanglbl=zero
DO jc=jstal,jstol
  DO kc=kstal,kstol
    umeanglbl=umeanglbl+uinf2(jc,kc)
    denom = denom+ustead(jc,kc)
  END DO
END DO
!      UMEANGLBL=UMEANGLBL/DENOM!FLOAT(NYSIZE*NZSIZE)
CALL mpi_allreduce(umeanglbl,sumumeanglbl,1,  &
    mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
CALL mpi_allreduce(denom,sumdenom,1,  &
    mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
uinf2=uinf2-sumumeanglbl/(sumdenom)

sumumeanglbl=zero
sumdenom=zero
denom = zero
umeanglbl=zero
DO jc=jstal,jstol
  DO kc=kstal,kstol
    umeanglbl=umeanglbl+vinf2(jc,kc)
    denom = denom+ustead(jc,kc)
  END DO
END DO
!      UMEANGLBL=UMEANGLBL/DENOM!FLOAT(NYSIZE*NZSIZE)
CALL mpi_allreduce(umeanglbl,sumumeanglbl,1,  &
    mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
CALL mpi_allreduce(denom,sumdenom,1,  &
    mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
vinf2=vinf2-sumumeanglbl/(sumdenom)

sumumeanglbl=zero
sumdenom=zero
denom = zero
umeanglbl=zero
DO jc=jstal,jstol
  DO kc=kstal,kstol
    umeanglbl=umeanglbl+winf2(jc,kc)
    denom = denom+ustead(jc,kc)
  END DO
END DO
!      UMEANGLBL=UMEANGLBL/DENOM!FLOAT(NYSIZE*NZSIZE)
CALL mpi_allreduce(umeanglbl,sumumeanglbl,1,  &
    mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
CALL mpi_allreduce(denom,sumdenom,1,  &
    mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
winf2=winf2-sumumeanglbl/(sumdenom)

!     -------------------------------------------------------------------------

END SUBROUTINE inflow

!     ==================================================================

FUNCTION ran1(idum)

INTEGER, INTENT(OUT)                     :: idum

REAL :: ran1
INTEGER, PARAMETER :: ia=16807
INTEGER, PARAMETER :: im=2147483647
REAL, PARAMETER :: am=1./im
INTEGER, PARAMETER :: iq=127773
INTEGER, PARAMETER :: ir=2836
INTEGER, PARAMETER :: ntab=32
INTEGER, PARAMETER :: ndiv=1+(im-1)/ntab
REAL, PARAMETER :: eps=1.2E-7
REAL, PARAMETER :: rnmx=1.-eps
INTEGER :: j,k,iv(ntab),iy
SAVE iv,iy
DATA iv /ntab*0/, iy /0/

IF (idum <= 0.OR.iy == 0) THEN
  idum=MAX(-idum,1)
  DO  j=ntab+8,1,-1
    k=idum/iq
    idum=ia*(idum-k*iq)-ir*k
    IF (idum < 0) idum=idum+im
    IF (j <= ntab) iv(j)=idum
  END DO
  iy=iv(1)
END IF
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
IF (idum < 0) idum=idum+im
j=1+iy/ndiv
iy=iv(j)
iv(j)=idum
ran1=MIN(am*iy,rnmx)
ran1=(ran1*2.0-1.0)/0.577
RETURN
END FUNCTION ran1
!     ==================================================================

