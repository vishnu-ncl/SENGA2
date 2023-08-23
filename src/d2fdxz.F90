SUBROUTINE d2fdxz(functn,fderiv)

! Code converted using TO_F90 by Alan Miller
! Date: 2022-11-09  Time: 14:51:24

!     *************************************************************************

!     D2FDXZ
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     15-APR-2003:  RSC MODIFIED FOR SENGA2

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES SECOND XZ-DERIVATIVE OF SPECIFIED FUNCTION
!     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
!     EXPLICIT 8TH,6TH,4TH,4TH,4TH COMPATIBLE ORDER END CONDITIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------

!     ARGUMENTS
!     =========

real(kind=8),INTENT(IN)             :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
real(kind=8),INTENT(OUT)            :: fderiv(nxsize,nysize,nzsize)




!     LOCAL DATA
!     ==========
real(kind=8):: fstora(3,3),fstorb(3,3),fstorc(2,2)
real(kind=8):: fdiffa,fdiffb,fdiffc,fdiffd,fdiffe
INTEGER :: ic,jc,kc
INTEGER :: is,ks,ism1,ksm1
INTEGER :: istart,ifinis,kstart,kfinis
INTEGER :: icm5,icm4,icm3,icm2,icm1,iccc,icp1,icp2,icp3,icp4,icp5
INTEGER :: kcm5,kcm4,kcm3,kcm2,kcm1,kccc,kcp1,kcp2,kcp3,kcp4,kcp5


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

istart = istal
ifinis = istol
kstart = kstal
kfinis = kstol
IF(nendxl == nbound)istart = istap5
IF(nendxr == nbound)ifinis = istom5
IF(nendzl == nbound)kstart = kstap5
IF(nendzr == nbound)kfinis = kstom5

!     =========================================================================

!     INTERIOR SCHEME
!     ===============

!     TENTH ORDER EXPLICIT DIFFERENCES

kcm4 = kstart-5
kcm3 = kstart-4
kcm2 = kstart-3
kcm1 = kstart-2
kccc = kstart-1
kcp1 = kstart
kcp2 = kstart+1
kcp3 = kstart+2
kcp4 = kstart+3
kcp5 = kstart+4

DO kc = kstart,kfinis
  
  kcm5 = kcm4
  kcm4 = kcm3
  kcm3 = kcm2
  kcm2 = kcm1
  kcm1 = kccc
  kccc = kcp1
  kcp1 = kcp2
  kcp2 = kcp3
  kcp3 = kcp4
  kcp4 = kcp5
  kcp5 = kc+5
  
  DO jc = jstal,jstol
    
    icm4 = istart-5
    icm3 = istart-4
    icm2 = istart-3
    icm1 = istart-2
    iccc = istart-1
    icp1 = istart
    icp2 = istart+1
    icp3 = istart+2
    icp4 = istart+3
    icp5 = istart+4
    
    DO ic = istart,ifinis
      
      icm5 = icm4
      icm4 = icm3
      icm3 = icm2
      icm2 = icm1
      icm1 = iccc
      iccc = icp1
      icp1 = icp2
      icp2 = icp3
      icp3 = icp4
      icp4 = icp5
      icp5 = ic+5
      
      fdiffa = functn(icp1,jc,kcp1) - functn(icp1,jc,kcm1)  &
          - functn(icm1,jc,kcp1) + functn(icm1,jc,kcm1)
      fdiffb = functn(icp2,jc,kcp2) - functn(icp2,jc,kcm2)  &
          - functn(icm2,jc,kcp2) + functn(icm2,jc,kcm2)
      fdiffc = functn(icp3,jc,kcp3) - functn(icp3,jc,kcm3)  &
          - functn(icm3,jc,kcp3) + functn(icm3,jc,kcm3)
      fdiffd = functn(icp4,jc,kcp4) - functn(icp4,jc,kcm4)  &
          - functn(icm4,jc,kcp4) + functn(icm4,jc,kcm4)
      fdiffe = functn(icp5,jc,kcp5) - functn(icp5,jc,kcm5)  &
          - functn(icm5,jc,kcp5) + functn(icm5,jc,kcm5)
      
      fderiv(ic,jc,kc) = acofxz*fdiffa + bcofxz*fdiffb  &
          + ccofxz*fdiffc + dcofxz*fdiffd  &
          + ecofxz*fdiffe
      
    END DO
    
  END DO
  
END DO

!     =========================================================================

!     LH END X-DIRECTION
!     ==================
IF(nendxl == nbound)THEN
  
!       TAKE SECOND XZ-DERIVATIVE IN X-LEFT INNER HALO
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  kcm3 = kstart-4
  kcm2 = kstart-3
  kcm1 = kstart-2
  kccc = kstart-1
  kcp1 = kstart
  kcp2 = kstart+1
  kcp3 = kstart+2
  kcp4 = kstart+3
  
  DO kc = kstart,kfinis
    
    kcm4 = kcm3
    kcm3 = kcm2
    kcm2 = kcm1
    kcm1 = kccc
    kccc = kcp1
    kcp1 = kcp2
    kcp2 = kcp3
    kcp3 = kcp4
    kcp4 = kc+4
    
    DO jc = jstal,jstol
      
!           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofz1*(functn(istap1,jc,kcp1) - functn(istap1,jc,kcm1)  &
          - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
          + bcofz1*(functn(istap1,jc,kcp2) - functn(istap1,jc,kcm2)  &
          - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
      fdiffb = acofz1*(functn(istap2,jc,kcp1) - functn(istap2,jc,kcm1)  &
          - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
          + bcofz1*(functn(istap2,jc,kcp2) - functn(istap2,jc,kcm2)  &
          - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
      fdiffc = acofz1*(functn(istap3,jc,kcp1) - functn(istap3,jc,kcm1)  &
          - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
          + bcofz1*(functn(istap3,jc,kcp2) - functn(istap3,jc,kcm2)  &
          - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
      fdiffd = acofz1*(functn(istap4,jc,kcp1) - functn(istap4,jc,kcm1)  &
          - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
          + bcofz1*(functn(istap4,jc,kcp2) - functn(istap4,jc,kcm2)  &
          - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
      fderiv(istal,jc,kc) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofz1*(functn(istal,jc,kcp1)  - functn(istal,jc,kcm1)  &
          - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
          + bcofz1*(functn(istal,jc,kcp2)  - functn(istal,jc,kcm2)  &
          - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
      fdiffb = acofz1*(functn(istap2,jc,kcp1) - functn(istap2,jc,kcm1)  &
          - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
          + bcofz1*(functn(istap2,jc,kcp2) - functn(istap2,jc,kcm2)  &
          - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
      fdiffc = acofz1*(functn(istap3,jc,kcp1) - functn(istap3,jc,kcm1)  &
          - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
          + bcofz1*(functn(istap3,jc,kcp2) - functn(istap3,jc,kcm2)  &
          - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
      fdiffd = acofz1*(functn(istap4,jc,kcp1) - functn(istap4,jc,kcm1)  &
          - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
          + bcofz1*(functn(istap4,jc,kcp2) - functn(istap4,jc,kcm2)  &
          - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
      fderiv(istap1,jc,kc) = acf2xz*fdiffa + bcf2xz*fdiffb  &
          + ccf2xz*fdiffc + dcf2xz*fdiffd
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdiffa = functn(istap3,jc,kcp1) - functn(istap3,jc,kcm1)  &
          - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1)
      fdiffb = functn(istap4,jc,kcp2) - functn(istap4,jc,kcm2)  &
          - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2)
      fderiv(istap2,jc,kc) = acf3xz*fdiffa + bcf3xz*fdiffb
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdiffa = functn(istap4,jc,kcp1) - functn(istap4,jc,kcm1)  &
          - functn(istap2,jc,kcp1) + functn(istap2,jc,kcm1)
      fdiffb = functn(istap5,jc,kcp2) - functn(istap5,jc,kcm2)  &
          - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2)
      fdiffc = functn(istap6,jc,kcp3) - functn(istap6,jc,kcm3)  &
          - functn(istal,jc,kcp3)  + functn(istal,jc,kcm3)
      fderiv(istap3,jc,kc) = acf4xz*fdiffa + bcf4xz*fdiffb  &
          + ccf4xz*fdiffc
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdiffa = functn(istap5,jc,kcp1) - functn(istap5,jc,kcm1)  &
          - functn(istap3,jc,kcp1) + functn(istap3,jc,kcm1)
      fdiffb = functn(istap6,jc,kcp2) - functn(istap6,jc,kcm2)  &
          - functn(istap2,jc,kcp2) + functn(istap2,jc,kcm2)
      fdiffc = functn(istap7,jc,kcp3) - functn(istap7,jc,kcm3)  &
          - functn(istap1,jc,kcp3) + functn(istap1,jc,kcm3)
      fdiffd = functn(istap8,jc,kcp4) - functn(istap8,jc,kcm4)  &
          - functn(istal,jc,kcp4)  + functn(istal,jc,kcm4)
      fderiv(istap4,jc,kc) = acf5xz*fdiffa + bcf5xz*fdiffb  &
          + ccf5xz*fdiffc + dcf5xz*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     RH END X-DIRECTION
!     ==================
IF(nendxr == nbound)THEN
  
!       TAKE SECOND XZ-DERIVATIVE IN X-RIGHT INNER HALO
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  kcm3 = kstart-4
  kcm2 = kstart-3
  kcm1 = kstart-2
  kccc = kstart-1
  kcp1 = kstart
  kcp2 = kstart+1
  kcp3 = kstart+2
  kcp4 = kstart+3
  
  DO kc = kstart,kfinis
    
    kcm4 = kcm3
    kcm3 = kcm2
    kcm2 = kcm1
    kcm1 = kccc
    kccc = kcp1
    kcp1 = kcp2
    kcp2 = kcp3
    kcp3 = kcp4
    kcp4 = kc+4
    
    DO jc = jstal,jstol
      
!           RH POINT MINUS 4: 8TH ORDER CENTRED
      fdiffa = functn(istom3,jc,kcp1) - functn(istom3,jc,kcm1)  &
          - functn(istom5,jc,kcp1) + functn(istom5,jc,kcm1)
      fdiffb = functn(istom2,jc,kcp2) - functn(istom2,jc,kcm2)  &
          - functn(istom6,jc,kcp2) + functn(istom6,jc,kcm2)
      fdiffc = functn(istom1,jc,kcp3) - functn(istom1,jc,kcm3)  &
          - functn(istom7,jc,kcp3) + functn(istom7,jc,kcm3)
      fdiffd = functn(istol,jc,kcp4)  - functn(istol,jc,kcm4)  &
          - functn(istom8,jc,kcp4) + functn(istom8,jc,kcm4)
      fderiv(istom4,jc,kc) = acf5xz*fdiffa + bcf5xz*fdiffb  &
          + ccf5xz*fdiffc + dcf5xz*fdiffd
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdiffa = functn(istom2,jc,kcp1) - functn(istom2,jc,kcm1)  &
          - functn(istom4,jc,kcp1) + functn(istom4,jc,kcm1)
      fdiffb = functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
          - functn(istom5,jc,kcp2) + functn(istom5,jc,kcm2)
      fdiffc = functn(istol,jc,kcp3)  - functn(istol,jc,kcm3)  &
          - functn(istom6,jc,kcp3) + functn(istom6,jc,kcm3)
      fderiv(istom3,jc,kc) = acf4xz*fdiffa + bcf4xz*fdiffb  &
          + ccf4xz*fdiffc
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdiffa = functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
          - functn(istom3,jc,kcp1) + functn(istom3,jc,kcm1)
      fdiffb = functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
          - functn(istom4,jc,kcp2) + functn(istom4,jc,kcm2)
      fderiv(istom2,jc,kc) = acf3xz*fdiffa + bcf3xz*fdiffb
      
!           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
          - functn(istol,jc,kcp1)  + functn(istol,jc,kcm1))  &
          + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
          - functn(istol,jc,kcp2)  + functn(istol,jc,kcm2))
      fdiffb = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
          - functn(istom2,jc,kcp1) + functn(istom2,jc,kcm1))  &
          + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
          - functn(istom2,jc,kcp2) + functn(istom2,jc,kcm2))
      fdiffc = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
          - functn(istom3,jc,kcp1) + functn(istom3,jc,kcm1))  &
          + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
          - functn(istom3,jc,kcp2) + functn(istom3,jc,kcm2))
      fdiffd = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
          - functn(istom4,jc,kcp1) + functn(istom4,jc,kcm1))  &
          + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
          - functn(istom4,jc,kcp2) + functn(istom4,jc,kcm2))
      fderiv(istom1,jc,kc) = acf2xz*fdiffa + bcf2xz*fdiffb  &
          + ccf2xz*fdiffc + dcf2xz*fdiffd
      
!           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
          - functn(istom1,jc,kcp1) + functn(istom1,jc,kcm1))  &
          + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
          - functn(istom1,jc,kcp2) + functn(istom1,jc,kcm2))
      fdiffb = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
          - functn(istom2,jc,kcp1) + functn(istom2,jc,kcm1))  &
          + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
          - functn(istom2,jc,kcp2) + functn(istom2,jc,kcm2))
      fdiffc = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
          - functn(istom3,jc,kcp1) + functn(istom3,jc,kcm1))  &
          + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
          - functn(istom3,jc,kcp2) + functn(istom3,jc,kcm2))
      fdiffd = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
          - functn(istom4,jc,kcp1) + functn(istom4,jc,kcm1))  &
          + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
          - functn(istom4,jc,kcp2) + functn(istom4,jc,kcm2))
      fderiv(istol,jc,kc) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     LH END Z-DIRECTION
!     ==================
IF(nendzl == nbound)THEN
  
!       TAKE SECOND XZ-DERIVATIVE IN Z-LEFT INNER HALO
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  DO jc = jstal,jstol
    
    icm3 = istart-4
    icm2 = istart-3
    icm1 = istart-2
    iccc = istart-1
    icp1 = istart
    icp2 = istart+1
    icp3 = istart+2
    icp4 = istart+3
    
    DO ic = istart,ifinis
      
      icm4 = icm3
      icm3 = icm2
      icm2 = icm1
      icm1 = iccc
      iccc = icp1
      icp1 = icp2
      icp2 = icp3
      icp3 = icp4
      icp4 = ic+4
      
!           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofx1*(functn(icp1,jc,kstap1) - functn(icm1,jc,kstap1)  &
          - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
          + bcofx1*(functn(icp2,jc,kstap1) - functn(icm2,jc,kstap1)  &
          - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
      fdiffb = acofx1*(functn(icp1,jc,kstap2) - functn(icm1,jc,kstap2)  &
          - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
          + bcofx1*(functn(icp2,jc,kstap2) - functn(icm2,jc,kstap2)  &
          - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
      fdiffc = acofx1*(functn(icp1,jc,kstap3) - functn(icm1,jc,kstap3)  &
          - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
          + bcofx1*(functn(icp2,jc,kstap3) - functn(icm2,jc,kstap3)  &
          - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
      fdiffd = acofx1*(functn(icp1,jc,kstap4) - functn(icm1,jc,kstap4)  &
          - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
          + bcofx1*(functn(icp2,jc,kstap4) - functn(icm2,jc,kstap4)  &
          - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
      fderiv(ic,jc,kstal) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofx1*(functn(icp1,jc,kstal)  - functn(icm1,jc,kstal)  &
          - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
          + bcofx1*(functn(icp2,jc,kstal)  - functn(icm2,jc,kstal)  &
          - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
      fdiffb = acofx1*(functn(icp1,jc,kstap2) - functn(icm1,jc,kstap2)  &
          - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
          + bcofx1*(functn(icp2,jc,kstap2) - functn(icm2,jc,kstap2)  &
          - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
      fdiffc = acofx1*(functn(icp1,jc,kstap3) - functn(icm1,jc,kstap3)  &
          - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
          + bcofx1*(functn(icp2,jc,kstap3) - functn(icm2,jc,kstap3)  &
          - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
      fdiffd = acofx1*(functn(icp1,jc,kstap4) - functn(icm1,jc,kstap4)  &
          - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
          + bcofx1*(functn(icp2,jc,kstap4) - functn(icm2,jc,kstap4)  &
          - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
      fderiv(ic,jc,kstap1) = acf2xz*fdiffa + bcf2xz*fdiffb  &
          + ccf2xz*fdiffc + dcf2xz*fdiffd
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdiffa = functn(icp1,jc,kstap3) - functn(icm1,jc,kstap3)  &
          - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1)
      fdiffb = functn(icp2,jc,kstap4) - functn(icm2,jc,kstap4)  &
          - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal)
      fderiv(ic,jc,kstap2) = acf3xz*fdiffa + bcf3xz*fdiffb
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdiffa = functn(icp1,jc,kstap4) - functn(icm1,jc,kstap4)  &
          - functn(icp1,jc,kstap2) + functn(icm1,jc,kstap2)
      fdiffb = functn(icp2,jc,kstap5) - functn(icm2,jc,kstap5)  &
          - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1)
      fdiffc = functn(icp3,jc,kstap6) - functn(icm3,jc,kstap6)  &
          - functn(icp3,jc,kstal)  + functn(icm3,jc,kstal)
      fderiv(ic,jc,kstap3) = acf4xz*fdiffa + bcf4xz*fdiffb  &
          + ccf4xz*fdiffc
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdiffa = functn(icp1,jc,kstap5) - functn(icm1,jc,kstap5)  &
          - functn(icp1,jc,kstap3) + functn(icm1,jc,kstap3)
      fdiffb = functn(icp2,jc,kstap6) - functn(icm2,jc,kstap6)  &
          - functn(icp2,jc,kstap2) + functn(icm2,jc,kstap2)
      fdiffc = functn(icp3,jc,kstap7) - functn(icm3,jc,kstap7)  &
          - functn(icp3,jc,kstap1) + functn(icm3,jc,kstap1)
      fdiffd = functn(icp4,jc,kstap8) - functn(icm4,jc,kstap8)  &
          - functn(icp4,jc,kstal)  + functn(icm4,jc,kstal)
      fderiv(ic,jc,kstap4) = acf5xz*fdiffa + bcf5xz*fdiffb  &
          + ccf5xz*fdiffc + dcf5xz*fdiffd
      
    END DO
  END DO
  
!       LH IN X LH IN Z CORNER
!       ======================
  IF(nendxl == nbound)THEN
    
    DO jc = jstal,jstol
      
!           LH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(istap1,jc,kstap1) - functn(istap1,jc,kstal)  &
          - functn(istal,jc,kstap1)  + functn(istal,jc,kstal)
      fdiffb = functn(istap2,jc,kstap2) - functn(istap2,jc,kstal)  &
          - functn(istal,jc,kstap2)  + functn(istal,jc,kstal)
      fdiffc = functn(istap3,jc,kstap3) - functn(istap3,jc,kstal)  &
          - functn(istal,jc,kstap3)  + functn(istal,jc,kstal)
      fdiffd = functn(istap4,jc,kstap4) - functn(istap4,jc,kstal)  &
          - functn(istal,jc,kstap4)  + functn(istal,jc,kstal)
      fderiv(istal,jc,kstal) = acc1xz*fdiffa + bcc1xz*fdiffb  &
          + ccc1xz*fdiffc + dcc1xz*fdiffd
      
!           LH+1 LH+1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(istal,jc,kstal)   - functn(istal,jc,kstap1)  &
          - functn(istap1,jc,kstal)  + functn(istap1,jc,kstap1)
      fdiffb = functn(istap2,jc,kstap2) - functn(istap2,jc,kstap1)  &
          - functn(istap1,jc,kstap2) + functn(istap1,jc,kstap1)
      fdiffc = functn(istap3,jc,kstap3) - functn(istap3,jc,kstap1)  &
          - functn(istap1,jc,kstap3) + functn(istap1,jc,kstap1)
      fdiffd = functn(istap4,jc,kstap4) - functn(istap4,jc,kstap1)  &
          - functn(istap1,jc,kstap4) + functn(istap1,jc,kstap1)
      fderiv(istap1,jc,kstap1) = acc2xz*fdiffa + bcc2xz*fdiffb  &
          + ccc2xz*fdiffc + dcc2xz*fdiffd
      
!           LH LH+1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xz*(functn(istap1,jc,kstal)  - functn(istap1,jc,kstap1)  &
          - functn(istal,jc,kstal)   + functn(istal,jc,kstap1))  &
          + bcf2xz*(functn(istap1,jc,kstap2) - functn(istap1,jc,kstap1)  &
          - functn(istal,jc,kstap2)  + functn(istal,jc,kstap1))  &
          + ccf2xz*(functn(istap1,jc,kstap3) - functn(istap1,jc,kstap1)  &
          - functn(istal,jc,kstap3)  + functn(istal,jc,kstap1))  &
          + dcf2xz*(functn(istap1,jc,kstap4) - functn(istap1,jc,kstap1)  &
          - functn(istal,jc,kstap4)  + functn(istal,jc,kstap1))
      fdiffb = acf2xz*(functn(istap2,jc,kstal)  - functn(istap2,jc,kstap1)  &
          - functn(istal,jc,kstal)   + functn(istal,jc,kstap1))  &
          + bcf2xz*(functn(istap2,jc,kstap2) - functn(istap2,jc,kstap1)  &
          - functn(istal,jc,kstap2)  + functn(istal,jc,kstap1))  &
          + ccf2xz*(functn(istap2,jc,kstap3) - functn(istap2,jc,kstap1)  &
          - functn(istal,jc,kstap3)  + functn(istal,jc,kstap1))  &
          + dcf2xz*(functn(istap2,jc,kstap4) - functn(istap2,jc,kstap1)  &
          - functn(istal,jc,kstap4)  + functn(istal,jc,kstap1))
      fdiffc = acf2xz*(functn(istap3,jc,kstal)  - functn(istap3,jc,kstap1)  &
          - functn(istal,jc,kstal)   + functn(istal,jc,kstap1))  &
          + bcf2xz*(functn(istap3,jc,kstap2) - functn(istap3,jc,kstap1)  &
          - functn(istal,jc,kstap2)  + functn(istal,jc,kstap1))  &
          + ccf2xz*(functn(istap3,jc,kstap3) - functn(istap3,jc,kstap1)  &
          - functn(istal,jc,kstap3)  + functn(istal,jc,kstap1))  &
          + dcf2xz*(functn(istap3,jc,kstap4) - functn(istap3,jc,kstap1)  &
          - functn(istal,jc,kstap4)  + functn(istal,jc,kstap1))
      fdiffd = acf2xz*(functn(istap4,jc,kstal)  - functn(istap4,jc,kstap1)  &
          - functn(istal,jc,kstal)   + functn(istal,jc,kstap1))  &
          + bcf2xz*(functn(istap4,jc,kstap2) - functn(istap4,jc,kstap1)  &
          - functn(istal,jc,kstap2)  + functn(istal,jc,kstap1))  &
          + ccf2xz*(functn(istap4,jc,kstap3) - functn(istap4,jc,kstap1)  &
          - functn(istal,jc,kstap3)  + functn(istal,jc,kstap1))  &
          + dcf2xz*(functn(istap4,jc,kstap4) - functn(istap4,jc,kstap1)  &
          - functn(istal,jc,kstap4)  + functn(istal,jc,kstap1))
      fderiv(istal,jc,kstap1) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           LH+1 LH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xz*(functn(istal,jc,kstap1)  - functn(istap1,jc,kstap1)  &
          - functn(istal,jc,kstal)   + functn(istap1,jc,kstal))  &
          + bcf2xz*(functn(istap2,jc,kstap1) - functn(istap1,jc,kstap1)  &
          - functn(istap2,jc,kstal)  + functn(istap1,jc,kstal))  &
          + ccf2xz*(functn(istap3,jc,kstap1) - functn(istap1,jc,kstap1)  &
          - functn(istap3,jc,kstal)  + functn(istap1,jc,kstal))  &
          + dcf2xz*(functn(istap4,jc,kstap1) - functn(istap1,jc,kstap1)  &
          - functn(istap4,jc,kstal)  + functn(istap1,jc,kstal))
      fdiffb = acf2xz*(functn(istal,jc,kstap2)  - functn(istap1,jc,kstap2)  &
          - functn(istal,jc,kstal)   + functn(istap1,jc,kstal))  &
          + bcf2xz*(functn(istap2,jc,kstap2) - functn(istap1,jc,kstap2)  &
          - functn(istap2,jc,kstal)  + functn(istap1,jc,kstal))  &
          + ccf2xz*(functn(istap3,jc,kstap2) - functn(istap1,jc,kstap2)  &
          - functn(istap3,jc,kstal)  + functn(istap1,jc,kstal))  &
          + dcf2xz*(functn(istap4,jc,kstap2) - functn(istap1,jc,kstap2)  &
          - functn(istap4,jc,kstal)  + functn(istap1,jc,kstal))
      fdiffc = acf2xz*(functn(istal,jc,kstap3)  - functn(istap1,jc,kstap3)  &
          - functn(istal,jc,kstal)   + functn(istap1,jc,kstal))  &
          + bcf2xz*(functn(istap2,jc,kstap3) - functn(istap1,jc,kstap3)  &
          - functn(istap2,jc,kstal)  + functn(istap1,jc,kstal))  &
          + ccf2xz*(functn(istap3,jc,kstap3) - functn(istap1,jc,kstap3)  &
          - functn(istap3,jc,kstal)  + functn(istap1,jc,kstal))  &
          + dcf2xz*(functn(istap4,jc,kstap3) - functn(istap1,jc,kstap3)  &
          - functn(istap4,jc,kstal)  + functn(istap1,jc,kstal))
      fdiffd = acf2xz*(functn(istal,jc,kstap4)  - functn(istap1,jc,kstap4)  &
          - functn(istal,jc,kstal)   + functn(istap1,jc,kstal))  &
          + bcf2xz*(functn(istap2,jc,kstap4) - functn(istap1,jc,kstap4)  &
          - functn(istap2,jc,kstal)  + functn(istap1,jc,kstal))  &
          + ccf2xz*(functn(istap3,jc,kstap4) - functn(istap1,jc,kstap4)  &
          - functn(istap3,jc,kstal)  + functn(istap1,jc,kstal))  &
          + dcf2xz*(functn(istap4,jc,kstap4) - functn(istap1,jc,kstap4)  &
          - functn(istap4,jc,kstal)  + functn(istap1,jc,kstal))
      fderiv(istap1,jc,kstal) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           LH EDGE IN Z
      DO ic = istap2,istap4
        
        icm2 = ic-2
        icm1 = ic-1
        icp1 = ic+1
        icp2 = ic+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofx1*(functn(icp1,jc,kstap1) - functn(icm1,jc,kstap1)  &
            - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
            + bcofx1*(functn(icp2,jc,kstap1) - functn(icm2,jc,kstap1)  &
            - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
        fdiffb = acofx1*(functn(icp1,jc,kstap2) - functn(icm1,jc,kstap2)  &
            - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
            + bcofx1*(functn(icp2,jc,kstap2) - functn(icm2,jc,kstap2)  &
            - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
        fdiffc = acofx1*(functn(icp1,jc,kstap3) - functn(icm1,jc,kstap3)  &
            - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
            + bcofx1*(functn(icp2,jc,kstap3) - functn(icm2,jc,kstap3)  &
            - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
        fdiffd = acofx1*(functn(icp1,jc,kstap4) - functn(icm1,jc,kstap4)  &
            - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
            + bcofx1*(functn(icp2,jc,kstap4) - functn(icm2,jc,kstap4)  &
            - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
        fderiv(ic,jc,kstal) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofx1*(functn(icp1,jc,kstal)  - functn(icm1,jc,kstal)  &
            - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
            + bcofx1*(functn(icp2,jc,kstal)  - functn(icm2,jc,kstal)  &
            - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
        fdiffb = acofx1*(functn(icp1,jc,kstap2) - functn(icm1,jc,kstap2)  &
            - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
            + bcofx1*(functn(icp2,jc,kstap2) - functn(icm2,jc,kstap2)  &
            - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
        fdiffc = acofx1*(functn(icp1,jc,kstap3) - functn(icm1,jc,kstap3)  &
            - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
            + bcofx1*(functn(icp2,jc,kstap3) - functn(icm2,jc,kstap3)  &
            - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
        fdiffd = acofx1*(functn(icp1,jc,kstap4) - functn(icm1,jc,kstap4)  &
            - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
            + bcofx1*(functn(icp2,jc,kstap4) - functn(icm2,jc,kstap4)  &
            - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
        fderiv(ic,jc,kstap1) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd
        
      END DO
      
!           LH EDGE IN X
      DO kc = kstap2,kstap4
        
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofz1*(functn(istap1,jc,kcp1) - functn(istap1,jc,kcm1)  &
            - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
            + bcofz1*(functn(istap1,jc,kcp2) - functn(istap1,jc,kcm2)  &
            - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
        fdiffb = acofz1*(functn(istap2,jc,kcp1) - functn(istap2,jc,kcm1)  &
            - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
            + bcofz1*(functn(istap2,jc,kcp2) - functn(istap2,jc,kcm2)  &
            - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
        fdiffc = acofz1*(functn(istap3,jc,kcp1) - functn(istap3,jc,kcm1)  &
            - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
            + bcofz1*(functn(istap3,jc,kcp2) - functn(istap3,jc,kcm2)  &
            - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
        fdiffd = acofz1*(functn(istap4,jc,kcp1) - functn(istap4,jc,kcm1)  &
            - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
            + bcofz1*(functn(istap4,jc,kcp2) - functn(istap4,jc,kcm2)  &
            - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
        fderiv(istal,jc,kc) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofz1*(functn(istal,jc,kcp1)  - functn(istal,jc,kcm1)  &
            - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
            + bcofz1*(functn(istal,jc,kcp2)  - functn(istal,jc,kcm2)  &
            - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
        fdiffb = acofz1*(functn(istap2,jc,kcp1) - functn(istap2,jc,kcm1)  &
            - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
            + bcofz1*(functn(istap2,jc,kcp2) - functn(istap2,jc,kcm2)  &
            - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
        fdiffc = acofz1*(functn(istap3,jc,kcp1) - functn(istap3,jc,kcm1)  &
            - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
            + bcofz1*(functn(istap3,jc,kcp2) - functn(istap3,jc,kcm2)  &
            - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
        fdiffd = acofz1*(functn(istap4,jc,kcp1) - functn(istap4,jc,kcm1)  &
            - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
            + bcofz1*(functn(istap4,jc,kcp2) - functn(istap4,jc,kcm2)  &
            - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
        fderiv(istap1,jc,kc) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      ks = 0
      DO kc = kstap2,kstap4
        
        ks = ks+1
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
        is = 0
        DO ic = istap2,istap4
          
          is = is+1
          icm2 = ic-2
          icm1 = ic-1
          icp1 = ic+1
          icp2 = ic+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(icp1,jc,kcp1) - functn(icp1,jc,kcm1)  &
              - functn(icm1,jc,kcp1) + functn(icm1,jc,kcm1)
          fdiffb = functn(icp2,jc,kcp2) - functn(icp2,jc,kcm2)  &
              - functn(icm2,jc,kcp2) + functn(icm2,jc,kcm2)
          fderiv(ic,jc,kc) = acf3xz*fdiffa + bcf3xz*fdiffb
          fstora(is,ks) = fdiffa
          fstorb(is,ks) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      ks = 1
      DO kc = kstap3,kstap4
        
        ksm1 = ks
        ks = ks+1
        kcm3 = kc-3
        kcp3 = kc+3
        
        is = 1
        DO ic = istap3,istap4
          
          ism1 = is
          is = is+1
          icm3 = ic-3
          icp3 = ic+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(icp3,jc,kcp3) - functn(icp3,jc,kcm3)  &
              - functn(icm3,jc,kcp3) + functn(icm3,jc,kcm3)
          fderiv(ic,jc,kc) = acf4xz*fstora(is,ks) + bcf4xz*fstorb(is,ks)  &
              + ccf4xz*fdiffc
          fstorc(ism1,ksm1) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      ks = 3
      is = 3
      ksm1 = 2
      ism1 = 2
      kc = kstap4
      ic = istap4
      kcm4 = kc-4
      kcp4 = kc+4
      icm4 = ic-4
      icp4 = ic+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(icp4,jc,kcp4) - functn(icp4,jc,kcm4)  &
          - functn(icm4,jc,kcp4) + functn(icm4,jc,kcm4)
      fderiv(ic,jc,kc) = acf5xz*fstora(is,ks) + bcf5xz*fstorb(is,ks)  &
          + ccf5xz*fstorc(ism1,ksm1) + dcf5xz*fdiffd
      
    END DO
    
  END IF
  
!       RH IN X LH IN Z CORNER
!       ======================
  IF(nendxr == nbound)THEN
    
    DO jc = jstal,jstol
      
!           RH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(istol,jc,kstap1)  - functn(istol,jc,kstal)  &
          - functn(istom1,jc,kstap1) + functn(istom1,jc,kstal)
      fdiffb = functn(istol,jc,kstap2)  - functn(istol,jc,kstal)  &
          - functn(istom2,jc,kstap2) + functn(istom2,jc,kstal)
      fdiffc = functn(istol,jc,kstap3)  - functn(istol,jc,kstal)  &
          - functn(istom3,jc,kstap3) + functn(istom3,jc,kstal)
      fdiffd = functn(istol,jc,kstap4)  - functn(istol,jc,kstal)  &
          - functn(istom4,jc,kstap4) + functn(istom4,jc,kstal)
      fderiv(istol,jc,kstal) = acc1xz*fdiffa + bcc1xz*fdiffb  &
          + ccc1xz*fdiffc + dcc1xz*fdiffd
      
!           RH-1 LH+1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(istom1,jc,kstal)  - functn(istom1,jc,kstap1)  &
          - functn(istol,jc,kstal)   + functn(istol,jc,kstap1)
      fdiffb = functn(istom1,jc,kstap2) - functn(istom1,jc,kstap1)  &
          - functn(istom2,jc,kstap2) + functn(istom2,jc,kstap1)
      fdiffc = functn(istom1,jc,kstap3) - functn(istom1,jc,kstap1)  &
          - functn(istom3,jc,kstap3) + functn(istom3,jc,kstap1)
      fdiffd = functn(istom1,jc,kstap4) - functn(istom1,jc,kstap1)  &
          - functn(istom4,jc,kstap4) + functn(istom4,jc,kstap1)
      fderiv(istom1,jc,kstap1) = acc2xz*fdiffa + bcc2xz*fdiffb  &
          + ccc2xz*fdiffc + dcc2xz*fdiffd
      
!           RH LH+1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xz*(functn(istol,jc,kstal)   - functn(istol,jc,kstap1)  &
          - functn(istom1,jc,kstal)  + functn(istom1,jc,kstap1))  &
          + bcf2xz*(functn(istol,jc,kstap2)  - functn(istol,jc,kstap1)  &
          - functn(istom1,jc,kstap2) + functn(istom1,jc,kstap1))  &
          + ccf2xz*(functn(istol,jc,kstap3)  - functn(istol,jc,kstap1)  &
          - functn(istom1,jc,kstap3) + functn(istom1,jc,kstap1))  &
          + dcf2xz*(functn(istol,jc,kstap4)  - functn(istol,jc,kstap1)  &
          - functn(istom1,jc,kstap4) + functn(istom1,jc,kstap1))
      fdiffb = acf2xz*(functn(istol,jc,kstal)   - functn(istol,jc,kstap1)  &
          - functn(istom2,jc,kstal)  + functn(istom2,jc,kstap1))  &
          + bcf2xz*(functn(istol,jc,kstap2)  - functn(istol,jc,kstap1)  &
          - functn(istom2,jc,kstap2) + functn(istom2,jc,kstap1))  &
          + ccf2xz*(functn(istol,jc,kstap3)  - functn(istol,jc,kstap1)  &
          - functn(istom2,jc,kstap3) + functn(istom2,jc,kstap1))  &
          + dcf2xz*(functn(istol,jc,kstap4)  - functn(istol,jc,kstap1)  &
          - functn(istom2,jc,kstap4) + functn(istom2,jc,kstap1))
      fdiffc = acf2xz*(functn(istol,jc,kstal)   - functn(istol,jc,kstap1)  &
          - functn(istom3,jc,kstal)  + functn(istom3,jc,kstap1))  &
          + bcf2xz*(functn(istol,jc,kstap2)  - functn(istol,jc,kstap1)  &
          - functn(istom3,jc,kstap2) + functn(istom3,jc,kstap1))  &
          + ccf2xz*(functn(istol,jc,kstap3)  - functn(istol,jc,kstap1)  &
          - functn(istom3,jc,kstap3) + functn(istom3,jc,kstap1))  &
          + dcf2xz*(functn(istol,jc,kstap4)  - functn(istol,jc,kstap1)  &
          - functn(istom3,jc,kstap4) + functn(istom3,jc,kstap1))
      fdiffd = acf2xz*(functn(istol,jc,kstal)   - functn(istol,jc,kstap1)  &
          - functn(istom4,jc,kstal)  + functn(istom4,jc,kstap1))  &
          + bcf2xz*(functn(istol,jc,kstap2)  - functn(istol,jc,kstap1)  &
          - functn(istom4,jc,kstap2) + functn(istom4,jc,kstap1))  &
          + ccf2xz*(functn(istol,jc,kstap3)  - functn(istol,jc,kstap1)  &
          - functn(istom4,jc,kstap3) + functn(istom4,jc,kstap1))  &
          + dcf2xz*(functn(istol,jc,kstap4)  - functn(istol,jc,kstap1)  &
          - functn(istom4,jc,kstap4) + functn(istom4,jc,kstap1))
      fderiv(istol,jc,kstap1) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           RH-1 LH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xz*(functn(istom1,jc,kstap1) - functn(istol,jc,kstap1)  &
          - functn(istom1,jc,kstal)  + functn(istol,jc,kstal))  &
          + bcf2xz*(functn(istom1,jc,kstap1) - functn(istom2,jc,kstap1)  &
          - functn(istom1,jc,kstal)  + functn(istom2,jc,kstal))  &
          + ccf2xz*(functn(istom1,jc,kstap1) - functn(istom3,jc,kstap1)  &
          - functn(istom1,jc,kstal)  + functn(istom3,jc,kstal))  &
          + dcf2xz*(functn(istom1,jc,kstap1) - functn(istom4,jc,kstap1)  &
          - functn(istom1,jc,kstal)  + functn(istom4,jc,kstal))
      fdiffb = acf2xz*(functn(istom1,jc,kstap2) - functn(istol,jc,kstap2)  &
          - functn(istom1,jc,kstal)  + functn(istol,jc,kstal))  &
          + bcf2xz*(functn(istom1,jc,kstap2) - functn(istom2,jc,kstap2)  &
          - functn(istom1,jc,kstal)  + functn(istom2,jc,kstal))  &
          + ccf2xz*(functn(istom1,jc,kstap2) - functn(istom3,jc,kstap2)  &
          - functn(istom1,jc,kstal)  + functn(istom3,jc,kstal))  &
          + dcf2xz*(functn(istom1,jc,kstap2) - functn(istom4,jc,kstap2)  &
          - functn(istom1,jc,kstal)  + functn(istom4,jc,kstal))
      fdiffc = acf2xz*(functn(istom1,jc,kstap3) - functn(istol,jc,kstap3)  &
          - functn(istom1,jc,kstal)  + functn(istol,jc,kstal))  &
          + bcf2xz*(functn(istom1,jc,kstap3) - functn(istom2,jc,kstap3)  &
          - functn(istom1,jc,kstal)  + functn(istom2,jc,kstal))  &
          + ccf2xz*(functn(istom1,jc,kstap3) - functn(istom3,jc,kstap3)  &
          - functn(istom1,jc,kstal)  + functn(istom3,jc,kstal))  &
          + dcf2xz*(functn(istom1,jc,kstap3) - functn(istom4,jc,kstap3)  &
          - functn(istom1,jc,kstal)  + functn(istom4,jc,kstal))
      fdiffd = acf2xz*(functn(istom1,jc,kstap4) - functn(istol,jc,kstap4)  &
          - functn(istom1,jc,kstal)  + functn(istol,jc,kstal))  &
          + bcf2xz*(functn(istom1,jc,kstap4) - functn(istom2,jc,kstap4)  &
          - functn(istom1,jc,kstal)  + functn(istom2,jc,kstal))  &
          + ccf2xz*(functn(istom1,jc,kstap4) - functn(istom3,jc,kstap4)  &
          - functn(istom1,jc,kstal)  + functn(istom3,jc,kstal))  &
          + dcf2xz*(functn(istom1,jc,kstap4) - functn(istom4,jc,kstap4)  &
          - functn(istom1,jc,kstal)  + functn(istom4,jc,kstal))
      fderiv(istom1,jc,kstal) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           LH EDGE IN Z
      DO ic = istom4,istom2
        
        icm2 = ic-2
        icm1 = ic-1
        icp1 = ic+1
        icp2 = ic+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofx1*(functn(icp1,jc,kstap1) - functn(icm1,jc,kstap1)  &
            - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
            + bcofx1*(functn(icp2,jc,kstap1) - functn(icm2,jc,kstap1)  &
            - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
        fdiffb = acofx1*(functn(icp1,jc,kstap2) - functn(icm1,jc,kstap2)  &
            - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
            + bcofx1*(functn(icp2,jc,kstap2) - functn(icm2,jc,kstap2)  &
            - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
        fdiffc = acofx1*(functn(icp1,jc,kstap3) - functn(icm1,jc,kstap3)  &
            - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
            + bcofx1*(functn(icp2,jc,kstap3) - functn(icm2,jc,kstap3)  &
            - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
        fdiffd = acofx1*(functn(icp1,jc,kstap4) - functn(icm1,jc,kstap4)  &
            - functn(icp1,jc,kstal)  + functn(icm1,jc,kstal))  &
            + bcofx1*(functn(icp2,jc,kstap4) - functn(icm2,jc,kstap4)  &
            - functn(icp2,jc,kstal)  + functn(icm2,jc,kstal))
        fderiv(ic,jc,kstal) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofx1*(functn(icp1,jc,kstal)  - functn(icm1,jc,kstal)  &
            - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
            + bcofx1*(functn(icp2,jc,kstal)  - functn(icm2,jc,kstal)  &
            - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
        fdiffb = acofx1*(functn(icp1,jc,kstap2) - functn(icm1,jc,kstap2)  &
            - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
            + bcofx1*(functn(icp2,jc,kstap2) - functn(icm2,jc,kstap2)  &
            - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
        fdiffc = acofx1*(functn(icp1,jc,kstap3) - functn(icm1,jc,kstap3)  &
            - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
            + bcofx1*(functn(icp2,jc,kstap3) - functn(icm2,jc,kstap3)  &
            - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
        fdiffd = acofx1*(functn(icp1,jc,kstap4) - functn(icm1,jc,kstap4)  &
            - functn(icp1,jc,kstap1) + functn(icm1,jc,kstap1))  &
            + bcofx1*(functn(icp2,jc,kstap4) - functn(icm2,jc,kstap4)  &
            - functn(icp2,jc,kstap1) + functn(icm2,jc,kstap1))
        fderiv(ic,jc,kstap1) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd
        
      END DO
      
!           RH EDGE IN X
      DO kc = kstap2,kstap4
        
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
            - functn(istom1,jc,kcp1) + functn(istom1,jc,kcm1))  &
            + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
            - functn(istom1,jc,kcp2) + functn(istom1,jc,kcm2))
        fdiffb = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
            - functn(istom2,jc,kcp1) + functn(istom2,jc,kcm1))  &
            + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
            - functn(istom2,jc,kcp2) + functn(istom2,jc,kcm2))
        fdiffc = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
            - functn(istom3,jc,kcp1) + functn(istom3,jc,kcm1))  &
            + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
            - functn(istom3,jc,kcp2) + functn(istom3,jc,kcm2))
        fdiffd = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
            - functn(istom4,jc,kcp1) + functn(istom4,jc,kcm1))  &
            + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
            - functn(istom4,jc,kcp2) + functn(istom4,jc,kcm2))
        fderiv(istol,jc,kc) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
            - functn(istol,jc,kcp1)  + functn(istol,jc,kcm1))  &
            + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
            - functn(istol,jc,kcp2)  + functn(istol,jc,kcm2))
        fdiffb = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
            - functn(istom2,jc,kcp1) + functn(istom2,jc,kcm1))  &
            + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
            - functn(istom2,jc,kcp2) + functn(istom2,jc,kcm2))
        fdiffc = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
            - functn(istom3,jc,kcp1) + functn(istom3,jc,kcm1))  &
            + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
            - functn(istom3,jc,kcp2) + functn(istom3,jc,kcm2))
        fdiffd = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
            - functn(istom4,jc,kcp1) + functn(istom4,jc,kcm1))  &
            + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
            - functn(istom4,jc,kcp2) + functn(istom4,jc,kcm2))
        fderiv(istom1,jc,kc) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      ks = 0
      DO kc = kstap2,kstap4
        
        ks = ks+1
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
        is = 0
        DO ic = istom4,istom2
          
          is = is+1
          icm2 = ic-2
          icm1 = ic-1
          icp1 = ic+1
          icp2 = ic+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(icp1,jc,kcp1) - functn(icp1,jc,kcm1)  &
              - functn(icm1,jc,kcp1) + functn(icm1,jc,kcm1)
          fdiffb = functn(icp2,jc,kcp2) - functn(icp2,jc,kcm2)  &
              - functn(icm2,jc,kcp2) + functn(icm2,jc,kcm2)
          fderiv(ic,jc,kc) = acf3xz*fdiffa + bcf3xz*fdiffb
          fstora(is,ks) = fdiffa
          fstorb(is,ks) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      ks = 1
      DO kc = kstap3,kstap4
        
        ksm1 = ks
        ks = ks+1
        kcm3 = kc-3
        kcp3 = kc+3
        
        is = 0
        DO ic = istom4,istom3
          
          is = is+1
          icm3 = ic-3
          icp3 = ic+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(icp3,jc,kcp3) - functn(icp3,jc,kcm3)  &
              - functn(icm3,jc,kcp3) + functn(icm3,jc,kcm3)
          fderiv(ic,jc,kc) = acf4xz*fstora(is,ks) + bcf4xz*fstorb(is,ks)  &
              + ccf4xz*fdiffc
          fstorc(is,ksm1) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      ks = 3
      is = 1
      ksm1 = 2
      kc = kstap4
      ic = istom4
      kcm4 = kc-4
      kcp4 = kc+4
      icm4 = ic-4
      icp4 = ic+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(icp4,jc,kcp4) - functn(icp4,jc,kcm4)  &
          - functn(icm4,jc,kcp4) + functn(icm4,jc,kcm4)
      fderiv(ic,jc,kc) = acf5xz*fstora(is,ks) + bcf5xz*fstorb(is,ks)  &
          + ccf5xz*fstorc(is,ksm1) + dcf5xz*fdiffd
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     RH END Z-DIRECTION
!     ==================
IF(nendzr == nbound)THEN
  
!       TAKE SECOND XZ-DERIVATIVE IN Z-RIGHT INNER HALO
!       EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  DO jc = jstal,jstol
    
    icm3 = istart-4
    icm2 = istart-3
    icm1 = istart-2
    iccc = istart-1
    icp1 = istart
    icp2 = istart+1
    icp3 = istart+2
    icp4 = istart+3
    
    DO ic = istart,ifinis
      
      icm4 = icm3
      icm3 = icm2
      icm2 = icm1
      icm1 = iccc
      iccc = icp1
      icp1 = icp2
      icp2 = icp3
      icp3 = icp4
      icp4 = ic+4
      
!           RH POINT MINUS 4: 8TH ORDER CENTRED
      fdiffa = functn(icp1,jc,kstom3) - functn(icm1,jc,kstom3)  &
          - functn(icp1,jc,kstom5) + functn(icm1,jc,kstom5)
      fdiffb = functn(icp2,jc,kstom2) - functn(icm2,jc,kstom2)  &
          - functn(icp2,jc,kstom6) + functn(icm2,jc,kstom6)
      fdiffc = functn(icp3,jc,kstom1) - functn(icm3,jc,kstom1)  &
          - functn(icp3,jc,kstom7) + functn(icm3,jc,kstom7)
      fdiffd = functn(icp4,jc,kstol)  - functn(icm4,jc,kstol)  &
          - functn(icp4,jc,kstom8) + functn(icm4,jc,kstom8)
      fderiv(ic,jc,kstom4) = acf5xz*fdiffa + bcf5xz*fdiffb  &
          + ccf5xz*fdiffc + dcf5xz*fdiffd
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdiffa = functn(icp1,jc,kstom2) - functn(icm1,jc,kstom2)  &
          - functn(icp1,jc,kstom4) + functn(icm1,jc,kstom4)
      fdiffb = functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
          - functn(icp2,jc,kstom5) + functn(icm2,jc,kstom5)
      fdiffc = functn(icp3,jc,kstol)  - functn(icm3,jc,kstol)  &
          - functn(icp3,jc,kstom6) + functn(icm3,jc,kstom6)
      fderiv(ic,jc,kstom3) = acf4xz*fdiffa + bcf4xz*fdiffb  &
          + ccf4xz*fdiffc
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdiffa = functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
          - functn(icp1,jc,kstom3) + functn(icm1,jc,kstom3)
      fdiffb = functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
          - functn(icp2,jc,kstom4) + functn(icm2,jc,kstom4)
      fderiv(ic,jc,kstom2) = acf3xz*fdiffa + bcf3xz*fdiffb
      
!           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
          - functn(icp1,jc,kstol)  + functn(icm1,jc,kstol))  &
          + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
          - functn(icp2,jc,kstol)  + functn(icm2,jc,kstol))
      fdiffb = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
          - functn(icp1,jc,kstom2) + functn(icm1,jc,kstom2))  &
          + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
          - functn(icp2,jc,kstom2) + functn(icm2,jc,kstom2))
      fdiffc = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
          - functn(icp1,jc,kstom3) + functn(icm1,jc,kstom3))  &
          + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
          - functn(icp2,jc,kstom3) + functn(icm2,jc,kstom3))
      fdiffd = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
          - functn(icp1,jc,kstom4) + functn(icm1,jc,kstom4))  &
          + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
          - functn(icp2,jc,kstom4) + functn(icm2,jc,kstom4))
      fderiv(ic,jc,kstom1) = acf2xz*fdiffa + bcf2xz*fdiffb  &
          + ccf2xz*fdiffc + dcf2xz*fdiffd
      
!           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
          - functn(icp1,jc,kstom1) + functn(icm1,jc,kstom1))  &
          + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
          - functn(icp2,jc,kstom1) + functn(icm2,jc,kstom1))
      fdiffb = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
          - functn(icp1,jc,kstom2) + functn(icm1,jc,kstom2))  &
          + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
          - functn(icp2,jc,kstom2) + functn(icm2,jc,kstom2))
      fdiffc = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
          - functn(icp1,jc,kstom3) + functn(icm1,jc,kstom3))  &
          + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
          - functn(icp2,jc,kstom3) + functn(icm2,jc,kstom3))
      fdiffd = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
          - functn(icp1,jc,kstom4) + functn(icm1,jc,kstom4))  &
          + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
          - functn(icp2,jc,kstom4) + functn(icm2,jc,kstom4))
      fderiv(ic,jc,kstol) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
    END DO
  END DO
  
!       LH IN X RH IN Z CORNER
!       ======================
  IF(nendxl == nbound)THEN
    
    DO jc = jstal,jstol
      
!           LH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(istap1,jc,kstol) - functn(istap1,jc,kstom1)  &
          - functn(istal,jc,kstol)  + functn(istal,jc,kstom1)
      fdiffb = functn(istap2,jc,kstol) - functn(istap2,jc,kstom2)  &
          - functn(istal,jc,kstol)  + functn(istal,jc,kstom2)
      fdiffc = functn(istap3,jc,kstol) - functn(istap3,jc,kstom3)  &
          - functn(istal,jc,kstol)  + functn(istal,jc,kstom3)
      fdiffd = functn(istap4,jc,kstol) - functn(istap4,jc,kstom4)  &
          - functn(istal,jc,kstol)  + functn(istal,jc,kstom4)
      fderiv(istal,jc,kstol) = acc1xz*fdiffa + bcc1xz*fdiffb  &
          + ccc1xz*fdiffc + dcc1xz*fdiffd
      
!           LH+1 RH-1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(istal,jc,kstom1)  - functn(istal,jc,kstol)  &
          - functn(istap1,jc,kstom1) + functn(istap1,jc,kstol)
      fdiffb = functn(istap2,jc,kstom1) - functn(istap2,jc,kstom2)  &
          - functn(istap1,jc,kstom1) + functn(istap1,jc,kstom2)
      fdiffc = functn(istap3,jc,kstom1) - functn(istap3,jc,kstom3)  &
          - functn(istap1,jc,kstom1) + functn(istap1,jc,kstom3)
      fdiffd = functn(istap4,jc,kstom1) - functn(istap4,jc,kstom4)  &
          - functn(istap1,jc,kstom1) + functn(istap1,jc,kstom4)
      fderiv(istap1,jc,kstom1) = acc2xz*fdiffa + bcc2xz*fdiffb  &
          + ccc2xz*fdiffc + dcc2xz*fdiffd
      
!           LH RH-1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xz*(functn(istap1,jc,kstom1) - functn(istap1,jc,kstol)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstol))  &
          + bcf2xz*(functn(istap1,jc,kstom1) - functn(istap1,jc,kstom2)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom2))  &
          + ccf2xz*(functn(istap1,jc,kstom1) - functn(istap1,jc,kstom3)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom3))  &
          + dcf2xz*(functn(istap1,jc,kstom1) - functn(istap1,jc,kstom4)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom4))
      fdiffb = acf2xz*(functn(istap2,jc,kstom1) - functn(istap2,jc,kstol)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstol))  &
          + bcf2xz*(functn(istap2,jc,kstom1) - functn(istap2,jc,kstom2)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom2))  &
          + ccf2xz*(functn(istap2,jc,kstom1) - functn(istap2,jc,kstom3)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom3))  &
          + dcf2xz*(functn(istap2,jc,kstom1) - functn(istap2,jc,kstom4)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom4))
      fdiffc = acf2xz*(functn(istap3,jc,kstom1) - functn(istap3,jc,kstol)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstol))  &
          + bcf2xz*(functn(istap3,jc,kstom1) - functn(istap3,jc,kstom2)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom2))  &
          + ccf2xz*(functn(istap3,jc,kstom1) - functn(istap3,jc,kstom3)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom3))  &
          + dcf2xz*(functn(istap3,jc,kstom1) - functn(istap3,jc,kstom4)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom4))
      fdiffd = acf2xz*(functn(istap4,jc,kstom1) - functn(istap4,jc,kstol)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstol))  &
          + bcf2xz*(functn(istap4,jc,kstom1) - functn(istap4,jc,kstom2)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom2))  &
          + ccf2xz*(functn(istap4,jc,kstom1) - functn(istap4,jc,kstom3)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom3))  &
          + dcf2xz*(functn(istap4,jc,kstom1) - functn(istap4,jc,kstom4)  &
          - functn(istal,jc,kstom1)  + functn(istal,jc,kstom4))
      fderiv(istal,jc,kstom1) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           LH+1 RH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xz*(functn(istal,jc,kstol)   - functn(istap1,jc,kstol)  &
          - functn(istal,jc,kstom1)  + functn(istap1,jc,kstom1))  &
          + bcf2xz*(functn(istap2,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap2,jc,kstom1) + functn(istap1,jc,kstom1))  &
          + ccf2xz*(functn(istap3,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap3,jc,kstom1) + functn(istap1,jc,kstom1))  &
          + dcf2xz*(functn(istap4,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap4,jc,kstom1) + functn(istap1,jc,kstom1))
      fdiffb = acf2xz*(functn(istal,jc,kstol)   - functn(istap1,jc,kstol)  &
          - functn(istal,jc,kstom2)  + functn(istap1,jc,kstom2))  &
          + bcf2xz*(functn(istap2,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap2,jc,kstom2) + functn(istap1,jc,kstom2))  &
          + ccf2xz*(functn(istap3,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap3,jc,kstom2) + functn(istap1,jc,kstom2))  &
          + dcf2xz*(functn(istap4,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap4,jc,kstom2) + functn(istap1,jc,kstom2))
      fdiffc = acf2xz*(functn(istal,jc,kstol)   - functn(istap1,jc,kstol)  &
          - functn(istal,jc,kstom3)  + functn(istap1,jc,kstom3))  &
          + bcf2xz*(functn(istap2,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap2,jc,kstom3) + functn(istap1,jc,kstom3))  &
          + ccf2xz*(functn(istap3,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap3,jc,kstom3) + functn(istap1,jc,kstom3))  &
          + dcf2xz*(functn(istap4,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap4,jc,kstom3) + functn(istap1,jc,kstom3))
      fdiffd = acf2xz*(functn(istal,jc,kstol)   - functn(istap1,jc,kstol)  &
          - functn(istal,jc,kstom4)  + functn(istap1,jc,kstom4))  &
          + bcf2xz*(functn(istap2,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap2,jc,kstom4) + functn(istap1,jc,kstom4))  &
          + ccf2xz*(functn(istap3,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap3,jc,kstom4) + functn(istap1,jc,kstom4))  &
          + dcf2xz*(functn(istap4,jc,kstol)  - functn(istap1,jc,kstol)  &
          - functn(istap4,jc,kstom4) + functn(istap1,jc,kstom4))
      fderiv(istap1,jc,kstol) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           RH EDGE IN Z
      DO ic = istap2,istap4
        
        icm2 = ic-2
        icm1 = ic-1
        icp1 = ic+1
        icp2 = ic+2
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
            - functn(icp1,jc,kstol)  + functn(icm1,jc,kstol))  &
            + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
            - functn(icp2,jc,kstol)  + functn(icm2,jc,kstol))
        fdiffb = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
            - functn(icp1,jc,kstom2) + functn(icm1,jc,kstom2))  &
            + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
            - functn(icp2,jc,kstom2) + functn(icm2,jc,kstom2))
        fdiffc = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
            - functn(icp1,jc,kstom3) + functn(icm1,jc,kstom3))  &
            + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
            - functn(icp2,jc,kstom3) + functn(icm2,jc,kstom3))
        fdiffd = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
            - functn(icp1,jc,kstom4) + functn(icm1,jc,kstom4))  &
            + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
            - functn(icp2,jc,kstom4) + functn(icm2,jc,kstom4))
        fderiv(ic,jc,kstom1) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
            - functn(icp1,jc,kstom1) + functn(icm1,jc,kstom1))  &
            + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
            - functn(icp2,jc,kstom1) + functn(icm2,jc,kstom1))
        fdiffb = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
            - functn(icp1,jc,kstom2) + functn(icm1,jc,kstom2))  &
            + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
            - functn(icp2,jc,kstom2) + functn(icm2,jc,kstom2))
        fdiffc = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
            - functn(icp1,jc,kstom3) + functn(icm1,jc,kstom3))  &
            + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
            - functn(icp2,jc,kstom3) + functn(icm2,jc,kstom3))
        fdiffd = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
            - functn(icp1,jc,kstom4) + functn(icm1,jc,kstom4))  &
            + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
            - functn(icp2,jc,kstom4) + functn(icm2,jc,kstom4))
        fderiv(ic,jc,kstol) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd
        
      END DO
      
!           LH EDGE IN X
      DO kc = kstom4,kstom2
        
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofz1*(functn(istap1,jc,kcp1) - functn(istap1,jc,kcm1)  &
            - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
            + bcofz1*(functn(istap1,jc,kcp2) - functn(istap1,jc,kcm2)  &
            - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
        fdiffb = acofz1*(functn(istap2,jc,kcp1) - functn(istap2,jc,kcm1)  &
            - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
            + bcofz1*(functn(istap2,jc,kcp2) - functn(istap2,jc,kcm2)  &
            - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
        fdiffc = acofz1*(functn(istap3,jc,kcp1) - functn(istap3,jc,kcm1)  &
            - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
            + bcofz1*(functn(istap3,jc,kcp2) - functn(istap3,jc,kcm2)  &
            - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
        fdiffd = acofz1*(functn(istap4,jc,kcp1) - functn(istap4,jc,kcm1)  &
            - functn(istal,jc,kcp1)  + functn(istal,jc,kcm1))  &
            + bcofz1*(functn(istap4,jc,kcp2) - functn(istap4,jc,kcm2)  &
            - functn(istal,jc,kcp2)  + functn(istal,jc,kcm2))
        fderiv(istal,jc,kc) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofz1*(functn(istal,jc,kcp1)  - functn(istal,jc,kcm1)  &
            - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
            + bcofz1*(functn(istal,jc,kcp2)  - functn(istal,jc,kcm2)  &
            - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
        fdiffb = acofz1*(functn(istap2,jc,kcp1) - functn(istap2,jc,kcm1)  &
            - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
            + bcofz1*(functn(istap2,jc,kcp2) - functn(istap2,jc,kcm2)  &
            - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
        fdiffc = acofz1*(functn(istap3,jc,kcp1) - functn(istap3,jc,kcm1)  &
            - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
            + bcofz1*(functn(istap3,jc,kcp2) - functn(istap3,jc,kcm2)  &
            - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
        fdiffd = acofz1*(functn(istap4,jc,kcp1) - functn(istap4,jc,kcm1)  &
            - functn(istap1,jc,kcp1) + functn(istap1,jc,kcm1))  &
            + bcofz1*(functn(istap4,jc,kcp2) - functn(istap4,jc,kcm2)  &
            - functn(istap1,jc,kcp2) + functn(istap1,jc,kcm2))
        fderiv(istap1,jc,kc) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      ks = 0
      DO kc = kstom4,kstom2
        
        ks = ks+1
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
        is = 0
        DO ic = istap2,istap4
          
          is = is+1
          icm2 = ic-2
          icm1 = ic-1
          icp1 = ic+1
          icp2 = ic+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(icp1,jc,kcp1) - functn(icp1,jc,kcm1)  &
              - functn(icm1,jc,kcp1) + functn(icm1,jc,kcm1)
          fdiffb = functn(icp2,jc,kcp2) - functn(icp2,jc,kcm2)  &
              - functn(icm2,jc,kcp2) + functn(icm2,jc,kcm2)
          fderiv(ic,jc,kc) = acf3xz*fdiffa + bcf3xz*fdiffb
          fstora(is,ks) = fdiffa
          fstorb(is,ks) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      ks = 0
      DO kc = kstom4,kstom3
        
        ks = ks+1
        kcm3 = kc-3
        kcp3 = kc+3
        
        is = 1
        DO ic = istap3,istap4
          
          ism1 = is
          is = is+1
          icm3 = ic-3
          icp3 = ic+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(icp3,jc,kcp3) - functn(icp3,jc,kcm3)  &
              - functn(icm3,jc,kcp3) + functn(icm3,jc,kcm3)
          fderiv(ic,jc,kc) = acf4xz*fstora(is,ks) + bcf4xz*fstorb(is,ks)  &
              + ccf4xz*fdiffc
          fstorc(ism1,ks) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      ks = 1
      is = 3
      ism1 = 2
      kc = kstom4
      ic = istap4
      kcm4 = kc-4
      kcp4 = kc+4
      icm4 = ic-4
      icp4 = ic+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(icp4,jc,kcp4) - functn(icp4,jc,kcm4)  &
          - functn(icm4,jc,kcp4) + functn(icm4,jc,kcm4)
      fderiv(ic,jc,kc) = acf5xz*fstora(is,ks) + bcf5xz*fstorb(is,ks)  &
          + ccf5xz*fstorc(ism1,ks) + dcf5xz*fdiffd
      
    END DO
    
  END IF
  
!       RH IN X RH IN Z CORNER
!       ======================
  IF(nendxr == nbound)THEN
    
    DO jc = jstal,jstol
      
!           RH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(istom1,jc,kstom1) - functn(istom1,jc,kstol)  &
          - functn(istol,jc,kstom1)  + functn(istol,jc,kstol)
      fdiffb = functn(istom2,jc,kstom2) - functn(istom2,jc,kstol)  &
          - functn(istol,jc,kstom2)  + functn(istol,jc,kstol)
      fdiffc = functn(istom3,jc,kstom3) - functn(istom3,jc,kstol)  &
          - functn(istol,jc,kstom3)  + functn(istol,jc,kstol)
      fdiffd = functn(istom4,jc,kstom4) - functn(istom4,jc,kstol)  &
          - functn(istol,jc,kstom4)  + functn(istol,jc,kstol)
      fderiv(istol,jc,kstol) = acc1xz*fdiffa + bcc1xz*fdiffb  &
          + ccc1xz*fdiffc + dcc1xz*fdiffd
      
!           RH-1 RH-1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(istol,jc,kstol)   - functn(istol,jc,kstom1)  &
          - functn(istom1,jc,kstol)  + functn(istom1,jc,kstom1)
      fdiffb = functn(istom2,jc,kstom2) - functn(istom2,jc,kstom1)  &
          - functn(istom1,jc,kstom2) + functn(istom1,jc,kstom1)
      fdiffc = functn(istom3,jc,kstom3) - functn(istom3,jc,kstom1)  &
          - functn(istom1,jc,kstom3) + functn(istom1,jc,kstom1)
      fdiffd = functn(istom4,jc,kstom4) - functn(istom4,jc,kstom1)  &
          - functn(istom1,jc,kstom4) + functn(istom1,jc,kstom1)
      fderiv(istom1,jc,kstom1) = acc2xz*fdiffa + bcc2xz*fdiffb  &
          + ccc2xz*fdiffc + dcc2xz*fdiffd
      
!           RH RH-1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xz*(functn(istom1,jc,kstol)  - functn(istom1,jc,kstom1)  &
          - functn(istol,jc,kstol)   + functn(istol,jc,kstom1))  &
          + bcf2xz*(functn(istom1,jc,kstom2) - functn(istom1,jc,kstom1)  &
          - functn(istol,jc,kstom2)  + functn(istol,jc,kstom1))  &
          + ccf2xz*(functn(istom1,jc,kstom3) - functn(istom1,jc,kstom1)  &
          - functn(istol,jc,kstom3)  + functn(istol,jc,kstom1))  &
          + dcf2xz*(functn(istom1,jc,kstom4) - functn(istom1,jc,kstom1)  &
          - functn(istol,jc,kstom4)  + functn(istol,jc,kstom1))
      fdiffb = acf2xz*(functn(istom2,jc,kstol)  - functn(istom2,jc,kstom1)  &
          - functn(istol,jc,kstol)   + functn(istol,jc,kstom1))  &
          + bcf2xz*(functn(istom2,jc,kstom2) - functn(istom2,jc,kstom1)  &
          - functn(istol,jc,kstom2)  + functn(istol,jc,kstom1))  &
          + ccf2xz*(functn(istom2,jc,kstom3) - functn(istom2,jc,kstom1)  &
          - functn(istol,jc,kstom3)  + functn(istol,jc,kstom1))  &
          + dcf2xz*(functn(istom2,jc,kstom4) - functn(istom2,jc,kstom1)  &
          - functn(istol,jc,kstom4)  + functn(istol,jc,kstom1))
      fdiffc = acf2xz*(functn(istom3,jc,kstol)  - functn(istom3,jc,kstom1)  &
          - functn(istol,jc,kstol)   + functn(istol,jc,kstom1))  &
          + bcf2xz*(functn(istom3,jc,kstom2) - functn(istom3,jc,kstom1)  &
          - functn(istol,jc,kstom2)  + functn(istol,jc,kstom1))  &
          + ccf2xz*(functn(istom3,jc,kstom3) - functn(istom3,jc,kstom1)  &
          - functn(istol,jc,kstom3)  + functn(istol,jc,kstom1))  &
          + dcf2xz*(functn(istom3,jc,kstom4) - functn(istom3,jc,kstom1)  &
          - functn(istol,jc,kstom4)  + functn(istol,jc,kstom1))
      fdiffd = acf2xz*(functn(istom4,jc,kstol)  - functn(istom4,jc,kstom1)  &
          - functn(istol,jc,kstol)   + functn(istol,jc,kstom1))  &
          + bcf2xz*(functn(istom4,jc,kstom2) - functn(istom4,jc,kstom1)  &
          - functn(istol,jc,kstom2)  + functn(istol,jc,kstom1))  &
          + ccf2xz*(functn(istom4,jc,kstom3) - functn(istom4,jc,kstom1)  &
          - functn(istol,jc,kstom3)  + functn(istol,jc,kstom1))  &
          + dcf2xz*(functn(istom4,jc,kstom4) - functn(istom4,jc,kstom1)  &
          - functn(istol,jc,kstom4)  + functn(istol,jc,kstom1))
      fderiv(istol,jc,kstom1) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           RH+1 RH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xz*(functn(istol,jc,kstom1)  - functn(istom1,jc,kstom1)  &
          - functn(istol,jc,kstol)   + functn(istom1,jc,kstol))  &
          + bcf2xz*(functn(istom2,jc,kstom1) - functn(istom1,jc,kstom1)  &
          - functn(istom2,jc,kstol)  + functn(istom1,jc,kstol))  &
          + ccf2xz*(functn(istom3,jc,kstom1) - functn(istom1,jc,kstom1)  &
          - functn(istom3,jc,kstol)  + functn(istom1,jc,kstol))  &
          + dcf2xz*(functn(istom4,jc,kstom1) - functn(istom1,jc,kstom1)  &
          - functn(istom4,jc,kstol)  + functn(istom1,jc,kstol))
      fdiffb = acf2xz*(functn(istol,jc,kstom2)  - functn(istom1,jc,kstom2)  &
          - functn(istol,jc,kstol)   + functn(istom1,jc,kstol))  &
          + bcf2xz*(functn(istom2,jc,kstom2) - functn(istom1,jc,kstom2)  &
          - functn(istom2,jc,kstol)  + functn(istom1,jc,kstol))  &
          + ccf2xz*(functn(istom3,jc,kstom2) - functn(istom1,jc,kstom2)  &
          - functn(istom3,jc,kstol)  + functn(istom1,jc,kstol))  &
          + dcf2xz*(functn(istom4,jc,kstom2) - functn(istom1,jc,kstom2)  &
          - functn(istom4,jc,kstol)  + functn(istom1,jc,kstol))
      fdiffc = acf2xz*(functn(istol,jc,kstom3)  - functn(istom1,jc,kstom3)  &
          - functn(istol,jc,kstol)   + functn(istom1,jc,kstol))  &
          + bcf2xz*(functn(istom2,jc,kstom3) - functn(istom1,jc,kstom3)  &
          - functn(istom2,jc,kstol)  + functn(istom1,jc,kstol))  &
          + ccf2xz*(functn(istom3,jc,kstom3) - functn(istom1,jc,kstom3)  &
          - functn(istom3,jc,kstol)  + functn(istom1,jc,kstol))  &
          + dcf2xz*(functn(istom4,jc,kstom3) - functn(istom1,jc,kstom3)  &
          - functn(istom4,jc,kstol)  + functn(istom1,jc,kstol))
      fdiffd = acf2xz*(functn(istol,jc,kstom4)  - functn(istom1,jc,kstom4)  &
          - functn(istol,jc,kstol)   + functn(istom1,jc,kstol))  &
          + bcf2xz*(functn(istom2,jc,kstom4) - functn(istom1,jc,kstom4)  &
          - functn(istom2,jc,kstol)  + functn(istom1,jc,kstol))  &
          + ccf2xz*(functn(istom3,jc,kstom4) - functn(istom1,jc,kstom4)  &
          - functn(istom3,jc,kstol)  + functn(istom1,jc,kstol))  &
          + dcf2xz*(functn(istom4,jc,kstom4) - functn(istom1,jc,kstom4)  &
          - functn(istom4,jc,kstol)  + functn(istom1,jc,kstol))
      fderiv(istom1,jc,kstol) = acf1xz*fdiffa + bcf1xz*fdiffb  &
          + ccf1xz*fdiffc + dcf1xz*fdiffd
      
!           RH EDGE IN Z
      DO ic = istom4,istom2
        
        icm2 = ic-2
        icm1 = ic-1
        icp1 = ic+1
        icp2 = ic+2
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
            - functn(icp1,jc,kstol)  + functn(icm1,jc,kstol))  &
            + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
            - functn(icp2,jc,kstol)  + functn(icm2,jc,kstol))
        fdiffb = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
            - functn(icp1,jc,kstom2) + functn(icm1,jc,kstom2))  &
            + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
            - functn(icp2,jc,kstom2) + functn(icm2,jc,kstom2))
        fdiffc = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
            - functn(icp1,jc,kstom3) + functn(icm1,jc,kstom3))  &
            + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
            - functn(icp2,jc,kstom3) + functn(icm2,jc,kstom3))
        fdiffd = acofx1*(functn(icp1,jc,kstom1) - functn(icm1,jc,kstom1)  &
            - functn(icp1,jc,kstom4) + functn(icm1,jc,kstom4))  &
            + bcofx1*(functn(icp2,jc,kstom1) - functn(icm2,jc,kstom1)  &
            - functn(icp2,jc,kstom4) + functn(icm2,jc,kstom4))
        fderiv(ic,jc,kstom1) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
            - functn(icp1,jc,kstom1) + functn(icm1,jc,kstom1))  &
            + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
            - functn(icp2,jc,kstom1) + functn(icm2,jc,kstom1))
        fdiffb = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
            - functn(icp1,jc,kstom2) + functn(icm1,jc,kstom2))  &
            + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
            - functn(icp2,jc,kstom2) + functn(icm2,jc,kstom2))
        fdiffc = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
            - functn(icp1,jc,kstom3) + functn(icm1,jc,kstom3))  &
            + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
            - functn(icp2,jc,kstom3) + functn(icm2,jc,kstom3))
        fdiffd = acofx1*(functn(icp1,jc,kstol)  - functn(icm1,jc,kstol)  &
            - functn(icp1,jc,kstom4) + functn(icm1,jc,kstom4))  &
            + bcofx1*(functn(icp2,jc,kstol)  - functn(icm2,jc,kstol)  &
            - functn(icp2,jc,kstom4) + functn(icm2,jc,kstom4))
        fderiv(ic,jc,kstol) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd
        
      END DO
      
!           RH EDGE IN X
      DO kc = kstom4,kstom2
        
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
            - functn(istom1,jc,kcp1) + functn(istom1,jc,kcm1))  &
            + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
            - functn(istom1,jc,kcp2) + functn(istom1,jc,kcm2))
        fdiffb = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
            - functn(istom2,jc,kcp1) + functn(istom2,jc,kcm1))  &
            + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
            - functn(istom2,jc,kcp2) + functn(istom2,jc,kcm2))
        fdiffc = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
            - functn(istom3,jc,kcp1) + functn(istom3,jc,kcm1))  &
            + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
            - functn(istom3,jc,kcp2) + functn(istom3,jc,kcm2))
        fdiffd = acofz1*(functn(istol,jc,kcp1)  - functn(istol,jc,kcm1)  &
            - functn(istom4,jc,kcp1) + functn(istom4,jc,kcm1))  &
            + bcofz1*(functn(istol,jc,kcp2)  - functn(istol,jc,kcm2)  &
            - functn(istom4,jc,kcp2) + functn(istom4,jc,kcm2))
        fderiv(istol,jc,kc) = acf1xz*fdiffa + bcf1xz*fdiffb  &
            + ccf1xz*fdiffc + dcf1xz*fdiffd
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
            - functn(istol,jc,kcp1)  + functn(istol,jc,kcm1))  &
            + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
            - functn(istol,jc,kcp2)  + functn(istol,jc,kcm2))
        fdiffb = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
            - functn(istom2,jc,kcp1) + functn(istom2,jc,kcm1))  &
            + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
            - functn(istom2,jc,kcp2) + functn(istom2,jc,kcm2))
        fdiffc = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
            - functn(istom3,jc,kcp1) + functn(istom3,jc,kcm1))  &
            + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
            - functn(istom3,jc,kcp2) + functn(istom3,jc,kcm2))
        fdiffd = acofz1*(functn(istom1,jc,kcp1) - functn(istom1,jc,kcm1)  &
            - functn(istom4,jc,kcp1) + functn(istom4,jc,kcm1))  &
            + bcofz1*(functn(istom1,jc,kcp2) - functn(istom1,jc,kcm2)  &
            - functn(istom4,jc,kcp2) + functn(istom4,jc,kcm2))
        fderiv(istom1,jc,kc) = acf2xz*fdiffa + bcf2xz*fdiffb  &
            + ccf2xz*fdiffc + dcf2xz*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      ks = 0
      DO kc = kstom4,kstom2
        
        ks = ks+1
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
        is = 0
        DO ic = istom4,istom2
          
          is = is+1
          icm2 = ic-2
          icm1 = ic-1
          icp1 = ic+1
          icp2 = ic+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(icp1,jc,kcp1) - functn(icp1,jc,kcm1)  &
              - functn(icm1,jc,kcp1) + functn(icm1,jc,kcm1)
          fdiffb = functn(icp2,jc,kcp2) - functn(icp2,jc,kcm2)  &
              - functn(icm2,jc,kcp2) + functn(icm2,jc,kcm2)
          fderiv(ic,jc,kc) = acf3xz*fdiffa + bcf3xz*fdiffb
          fstora(is,ks) = fdiffa
          fstorb(is,ks) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      ks = 0
      DO kc = kstom4,kstom3
        
        ks = ks+1
        kcm3 = kc-3
        kcp3 = kc+3
        
        is = 0
        DO ic = istom4,istom3
          
          is = is+1
          icm3 = ic-3
          icp3 = ic+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(icp3,jc,kcp3) - functn(icp3,jc,kcm3)  &
              - functn(icm3,jc,kcp3) + functn(icm3,jc,kcm3)
          fderiv(ic,jc,kc) = acf4xz*fstora(is,ks) + bcf4xz*fstorb(is,ks)  &
              + ccf4xz*fdiffc
          fstorc(is,ks) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      ks = 1
      is = 1
      kc = kstom4
      ic = istom4
      kcm4 = kc-4
      kcp4 = kc+4
      icm4 = ic-4
      icp4 = ic+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(icp4,jc,kcp4) - functn(icp4,jc,kcm4)  &
          - functn(icm4,jc,kcp4) + functn(icm4,jc,kcm4)
      fderiv(ic,jc,kc) = acf5xz*fstora(is,ks) + bcf5xz*fstorb(is,ks)  &
          + ccf5xz*fstorc(is,ks) + dcf5xz*fdiffd
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     SCALING
!     =======
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fderiv(ic,jc,kc) = fderiv(ic,jc,kc)*ovdelx*ovdelz
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE d2fdxz
