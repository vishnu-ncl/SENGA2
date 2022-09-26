SUBROUTINE turbin
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 15:26:33

!     *************************************************************************

!     TURBIN
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     28-MAR-1997:  CREATED
!     20-OCT-1999:  KWJ PARALLEL FFT IMPLEMENTATION
!     24-MAY-2003:  RSC UPDATED FOR SENGA2
!     08-SEP-2012:  RSC/RACG BUG FIX SPECIAL CASE OF ZERO 12-PLANE VECTOR

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     INITIAL TURBULENCE GENERATOR
!     COMPUTES A HOMOGENEOUS ISOTROPIC TURBULENT VELOCITY FIELD
!     SATISFYING CONTINUITY, ISOTROPY AND TOTAL ENERGY CONDITIONS

!     FIELD IS BUILT IN FOURIER SPACE ACCORDING TO ROGALLO (NASA TM 81315)

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     FUNCTIONS
!     =========
real(kind=dp) :: espect,ranuni
EXTERNAL espect,ranuni


!     PARAMETERS
!     ==========
real(kind=dp) :: vectol,tolimg
PARAMETER(vectol=0.00001_dp, tolimg=1.0E-6)


!     LOCAL DATA
!     ==========
real(kind=dp) :: vectk1,vectk2,vectk3,vksize,ovksiz
real(kind=dp) :: velmag,vfactr,plnmag
real(kind=dp) :: ovplmg,aziang,cosazi,sinazi
real(kind=dp) :: phang1,phang2,cosph1,sinph1,cosph2,sinph2
real(kind=dp) :: vaterm,vbterm
real(kind=dp) :: tklodd
real(kind=dp) :: twopi,ovtopi
real(kind=dp) :: ubart,vbart,wbart,uvart,vvart,wvart,tket
real(kind=dp) :: ubartt,vbartt,wbartt,uvartt,vvartt,wvartt,tketot
real(kind=dp) :: ubartl,vbartl,wbartl,uvartl,vvartl,wvartl,tketl
real(kind=dp) :: ubartg,vbartg,wbartg,uvartg,vvartg,wvartg,tketg
real(kind=dp) :: udev,vdev,wdev,faclav,facgav
INTEGER :: ic,jc,kc,ix,jx,kx,icproc
INTEGER :: igofst,jgofst,kgofst,igofm1,jgofm1,kgofm1
INTEGER :: igstal,jgstal,kgstal,igstol,jgstol,kgstol
INTEGER :: nodblx,nodbly,nodblz
INTEGER :: ngoddx,ngoddy,ngoddz
INTEGER :: iseed
LOGICAL :: flagim
LOGICAL :: flrani,flranj,flrank


!     BEGIN
!     =====

!     =========================================================================

!     INDEXING
!     --------

!     SET ODD-NUMBER INDICATORS
ngoddx = MOD(nxglbl,2)
ngoddy = MOD(nyglbl,2)
ngoddz = MOD(nzglbl,2)

!     SET ODDBALL WAVENUMBER INDICES
nodblx = nxglbl/2 - 1 + ngoddx
nodbly = nyglbl/2 - 1 + ngoddy
nodblz = nzglbl/2 - 1 + ngoddz

!     PHYSICAL-SPACE GLOBAL INDEX OFFSETS
igofst = 0
DO icproc = 0, ixproc-1
  igofst = igofst + npmapx(icproc)
END DO
igstal = igofst+1
igstol = igofst + npmapx(ixproc)
igofm1 = igofst-1

jgofst = 0
DO icproc = 0, iyproc-1
  jgofst = jgofst + npmapy(icproc)
END DO
jgstal = jgofst+1
jgstol = jgofst + npmapy(iyproc)
jgofm1 = jgofst-1

kgofst = 0
DO icproc = 0, izproc-1
  kgofst = kgofst + npmapz(icproc)
END DO
kgstal = kgofst+1
kgstol = kgofst + npmapz(izproc)
kgofm1 = kgofst-1

!     =========================================================================

!     SET OTHER FACTORS
twopi = two*pi
ovtopi = one/twopi
tklodd = zero

!     =========================================================================

!     INITIALISE THE RANDOM NUMBER GENERATOR
!     --------------------------------------
IF(inseed >= 0)THEN
  
!       LOCAL INITIALISATION
!       --------------------
!       RANDOM SEED MUST BE DIFFERENT FOR EVERY PROCESSOR
  iseed = inseed + iproc
  CALL ranini(iseed)
  
!       SWEEP THROUGH LOCAL PHYSICAL SPACE
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      DO ic = istal,istol
        
!             GET AND SAVE THREE RANDOM NUMBERS
        utmp(ic,jc,kc) = ranuni(iseed)
        vtmp(ic,jc,kc) = ranuni(iseed)
        wtmp(ic,jc,kc) = ranuni(iseed)
        
      END DO
    END DO
  END DO
  
ELSE
  
!       GLOBAL INITIALISATION
!       ---------------------
!       RANDOM SEED MUST BE IDENTICAL FOR EVERY PROCESSOR
  iseed = inseed
  CALL ranini(iseed)
  
!       SWEEP THROUGH GLOBAL PHYSICAL SPACE
  DO kc = 1,nzglbl
    DO jc = 1,nyglbl
      DO ic = 1,nxglbl
        
!             GET THREE RANDOM NUMBERS
        aziang = ranuni(iseed)
        phang1 = ranuni(iseed)
        phang2 = ranuni(iseed)
        
!             CHECK PHYSICAL-SPACE INDEX RANGE FOR THIS PROCESSOR
        flrani = (ic >= igstal).AND.(ic <= igstol)
        flranj = (jc >= jgstal).AND.(jc <= jgstol)
        flrank = (kc >= kgstal).AND.(kc <= kgstol)
        
        IF(flrani.AND.flranj.AND.flrank)THEN
          
!               SET LOCAL INDEXING
          ix = ic - igofst
          jx = jc - jgofst
          kx = kc - kgofst
          
!               SAVE THE THREE RANDOM NUMBERS
          utmp(ix,jx,kx) = aziang
          vtmp(ix,jx,kx) = phang1
          wtmp(ix,jx,kx) = phang2
          
        END IF
        
      END DO
    END DO
  END DO
  
END IF

!     =========================================================================

!     INITIALISE THE ENERGY SPECTRUM
!     ------------------------------
CALL espini

!     =========================================================================

!     EVALUATE VELOCITY COMPONENTS IN FOURIER SPACE
!     ---------------------------------------------

!     SWEEP THROUGH LOCAL PHYSICAL SPACE
!     ----------------------------------
DO kc = kstal,kstol
  
!       FOURIER-SPACE GLOBAL INDEXING
  kx = kgofm1 + kc
  IF(kx > nodblz)kx = kx-nzglbl
  
  DO jc = jstal,jstol
    
!         FOURIER-SPACE GLOBAL INDEXING
    jx = jgofm1 + jc
    IF(jx > nodbly)jx = jx-nyglbl
    
    DO ic = istal,istol
      
!           FOURIER-SPACE GLOBAL INDEXING
      ix = igofm1 + ic
      IF(ix > nodblx)ix = ix-nxglbl
      
!           CHECK LOCATION IN FOURIER SPACE
      flrank = (kx > 0)
      flranj = (kx == 0).AND.(jx > 0)
      flrani = (kx == 0).AND.(jx == 0).AND.(ix >= 0)
      
      IF(flrani.OR.flranj.OR.flrank)THEN
        
!             IN UPPER-CENTRAL HALF OF FOURIER SPACE
        
!             EVALUATE THE WAVENUMBER VECTOR COMPONENTS
        vectk1 = REAL(ix,kind=dp)
        vectk2 = REAL(jx,kind=dp)
        vectk3 = REAL(kx,kind=dp)
        vksize = SQRT(vectk1*vectk1+vectk2*vectk2+vectk3*vectk3)
!             SPECIAL CASE OF ZERO VECTOR
        ovksiz = one
        IF(vksize > vectol)ovksiz = one/vksize
        
!             WAVENUMBER VECTOR MAGNITUDE IN THE 12-PLANE
        plnmag = SQRT(vectk1*vectk1 + vectk2*vectk2)
!             SPECIAL CASE OF ZERO 12-PLANE VECTOR
!             RSC/RACG 08-SEP-2012 BUG FIX
!              OVPLMG = ONE
!              IF(PLNMAG.GT.VECTOL)OVPLMG = ONE/PLNMAG
        IF(plnmag < vectol)THEN
!               AZIMUTHAL ANGLE IS ARBITRARY/STOCHASTIC: CHOOSE K1=0 K2=1
          ovplmg = one
          vectk1 = zero
          vectk2 = one
        ELSE
          ovplmg = one/plnmag
        END IF
        
!             EVALUATE THE ENERGY SPECTRUM FUNCTION
!             ENERGY SPECTRUM SPECIFIED EXTERNALLY IN FUNCTION ESPECT
        velmag = espect(vksize)
        
!             EVALUATE THE FOURIER-VELOCITY MAGNITUDE
        velmag = SQRT(velmag*ovtopi)*ovksiz
        
!             GENERATE THE (RANDOM) AZIMUTH ANGLE
!             RANGE IS ZERO TO 2*PI
        aziang = utmp(ic,jc,kc)
        aziang = twopi*aziang
        cosazi = COS(aziang)
        sinazi = SIN(aziang)
        
!             GENERATE THE (RANDOM) PHASE ANGLES
!             RANGE IS -PI TO PI
        phang1 = vtmp(ic,jc,kc)
        phang1 = pi*(two*phang1-one)
        cosph1 = COS(phang1)
        sinph1 = SIN(phang1)
        phang2 = wtmp(ic,jc,kc)
        phang2 = pi*(two*phang2-one)
        cosph2 = COS(phang2)
        sinph2 = SIN(phang2)
        
!             EVALUATE COMMON FACTOR FOR ALL COMPONENTS
        vfactr = velmag*ovksiz*ovplmg
        
!             U-COMPONENT
        vaterm = vfactr*vksize*vectk2*cosazi
        vbterm = vfactr*vectk1*vectk3*sinazi
!             REAL PART
        urun(ic,jc,kc) = vaterm*cosph1 + vbterm*cosph2
!             IMAGINARY PART
        utmp(ic,jc,kc) = vaterm*sinph1 + vbterm*sinph2
        
!             V-COMPONENT
        vaterm = vfactr*vksize*vectk1*cosazi
        vbterm = vfactr*vectk2*vectk3*sinazi
!             REAL PART
        vrun(ic,jc,kc) = vbterm*cosph2 - vaterm*cosph1
!             IMAGINARY PART
        vtmp(ic,jc,kc) = vbterm*sinph2 - vaterm*sinph1
        
!             W-COMPONENT
!             VATERM IS ZERO
        vbterm = -vfactr*plnmag*plnmag*sinazi
!             REAL PART
        wrun(ic,jc,kc) = vbterm*cosph2
!             IMAGINARY PART
        wtmp(ic,jc,kc) = vbterm*sinph2
        
      ELSE
        
!             IN LOWER-CENTRAL HALF OF FOURIER SPACE
        
!             SET VELOCITY COMPONENT VALUES TO ZERO FOR NOW
        urun(ic,jc,kc) = zero
        utmp(ic,jc,kc) = zero
        vrun(ic,jc,kc) = zero
        vtmp(ic,jc,kc) = zero
        wrun(ic,jc,kc) = zero
        wtmp(ic,jc,kc) = zero
        
      END IF
      
    END DO
  END DO
END DO

!     =========================================================================

!     CHECK ENERGY CONTENT
!     --------------------
tket = zero
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      tket = tket + urun(ic,jc,kc)*urun(ic,jc,kc)  &
          + utmp(ic,jc,kc)*utmp(ic,jc,kc) + vrun(ic,jc,kc)*vrun(ic,jc,kc)  &
          + vtmp(ic,jc,kc)*vtmp(ic,jc,kc) + wrun(ic,jc,kc)*wrun(ic,jc,kc)  &
          + wtmp(ic,jc,kc)*wtmp(ic,jc,kc)
      
    END DO
  END DO
END DO

!     SUM OVER ALL PROCESSORS
CALL p_summ(tket,tketot)

!     REPORT
IF(iproc == 0)THEN
  WRITE(ncrept,*)'Fourier-space turbulence kinetic energy'
  WRITE(ncrept,'(1PE12.4)')tketot
END IF

!     =========================================================================

!     CARRY OUT AN INVERSE FOURIER TRANSFORM
!     ======================================
!     WITH CONJUGATE ANTISYMMETRY
CALL turbft

!     =========================================================================

!     CHECK THAT TRANSFORMED DATA IS REAL
!     -----------------------------------

DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      flagim = .false.
      IF(utmp(ic,jc,kc) > tolimg)flagim = .true.
      IF(vtmp(ic,jc,kc) > tolimg)flagim = .true.
      IF(wtmp(ic,jc,kc) > tolimg)flagim = .true.
      IF(flagim)THEN
        WRITE(6,*)'Warning: Imaginary part too large',iproc
        WRITE(6,'(3I5,3(1PE12.4))')ic,jc,kc,  &
            utmp(ic,jc,kc),vtmp(ic,jc,kc),wtmp(ic,jc,kc)
      END IF
      
    END DO
  END DO
END DO

!     =========================================================================

!     CHECK TURBULENCE STATISTICS
!     ---------------------------

!     AVERAGING FACTORS
faclav = one/REAL(nxnode,kind=dp)/REAL(nynode,kind=dp)/REAL(nznode,kind=dp)
facgav = one/REAL(nxglbl,kind=dp)/REAL(nyglbl,kind=dp)/REAL(nzglbl,kind=dp)


!     VELOCITY MEANS
!     --------------
!     SUM OVER LOCAL VALUES
ubart = zero
vbart = zero
wbart = zero
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      ubart = ubart + urun(ic,jc,kc)
      vbart = vbart + vrun(ic,jc,kc)
      wbart = wbart + wrun(ic,jc,kc)
      
    END DO
  END DO
END DO

!     LOCAL AVERAGES
ubartl = ubart*faclav
vbartl = vbart*faclav
wbartl = wbart*faclav

!     SUM OVER ALL PROCESSORS
CALL p_summ(ubart,ubartt)
CALL p_summ(vbart,vbartt)
CALL p_summ(wbart,wbartt)

!     GLOBAL AVERAGES
ubartg = ubartt*facgav
vbartg = vbartt*facgav
wbartg = wbartt*facgav


!     VELOCITY VARIANCES
!     ------------------
!     SUM OVER LOCAL VALUES
uvart = zero
vvart = zero
wvart = zero
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      udev = urun(ic,jc,kc)-ubartg
      uvart = uvart + udev*udev
      vdev = vrun(ic,jc,kc)-vbartg
      vvart = vvart + vdev*vdev
      wdev = wrun(ic,jc,kc)-wbartg
      wvart = wvart + wdev*wdev
      
    END DO
  END DO
END DO

!     TURBULENCE TOTAL ENERGY
tket   = half*(uvart+vvart+wvart)

!     LOCAL AVERAGES
uvartl = uvart*faclav
vvartl = vvart*faclav
wvartl = wvart*faclav
tketl  = half*(uvartl+vvartl+wvartl)

!     SUM OVER ALL PROCESSORS
CALL p_summ(uvart,uvartt)
CALL p_summ(vvart,vvartt)
CALL p_summ(wvart,wvartt)
CALL p_summ(tket,tketot)

!     GLOBAL AVERAGES
uvartg = uvartt*facgav
vvartg = vvartt*facgav
wvartg = wvartt*facgav
tketg  = tketot*facgav

!     REPORT
IF(iproc == 0)THEN
  WRITE(ncrept,*)'Physical-space turbulence kinetic energy'
  WRITE(ncrept,'(1PE12.4)')tketg
  WRITE(ncrept,*)
  WRITE(ncrept,*)'Velocity means:'
  WRITE(ncrept,'(3(1PE12.4))')ubartg,vbartg,wbartg
  WRITE(ncrept,*)'Velocity variances:'
  WRITE(ncrept,'(3(1PE12.4))')uvartg,vvartg,wvartg
  WRITE(ncrept,*)
END IF

!     =========================================================================


RETURN
END SUBROUTINE turbin
