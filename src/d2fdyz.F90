SUBROUTINE d2fdyz(functn,fderiv)

! Code converted using TO_F90 by Alan Miller
! Date: 2022-11-09  Time: 14:51:28

!     *************************************************************************

!     D2FDYZ
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
!     EVALUATES SECOND YZ-DERIVATIVE OF SPECIFIED FUNCTION
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
INTEGER :: js,ks,jsm1,ksm1
INTEGER :: jstart,jfinis,kstart,kfinis
INTEGER :: jcm5,jcm4,jcm3,jcm2,jcm1,jccc,jcp1,jcp2,jcp3,jcp4,jcp5
INTEGER :: kcm5,kcm4,kcm3,kcm2,kcm1,kccc,kcp1,kcp2,kcp3,kcp4,kcp5


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

jstart = jstal
jfinis = jstol
kstart = kstal
kfinis = kstol
IF(nendyl == nbound)jstart = jstap5
IF(nendyr == nbound)jfinis = jstom5
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
  
  jcm4 = jstart-5
  jcm3 = jstart-4
  jcm2 = jstart-3
  jcm1 = jstart-2
  jccc = jstart-1
  jcp1 = jstart
  jcp2 = jstart+1
  jcp3 = jstart+2
  jcp4 = jstart+3
  jcp5 = jstart+4
  
  DO jc = jstart,jfinis
    
    jcm5 = jcm4
    jcm4 = jcm3
    jcm3 = jcm2
    jcm2 = jcm1
    jcm1 = jccc
    jccc = jcp1
    jcp1 = jcp2
    jcp2 = jcp3
    jcp3 = jcp4
    jcp4 = jcp5
    jcp5 = jc+5
    
    DO ic = istal,istol
      
      fdiffa = functn(ic,jcp1,kcp1) - functn(ic,jcp1,kcm1)  &
          - functn(ic,jcm1,kcp1) + functn(ic,jcm1,kcm1)
      fdiffb = functn(ic,jcp2,kcp2) - functn(ic,jcp2,kcm2)  &
          - functn(ic,jcm2,kcp2) + functn(ic,jcm2,kcm2)
      fdiffc = functn(ic,jcp3,kcp3) - functn(ic,jcp3,kcm3)  &
          - functn(ic,jcm3,kcp3) + functn(ic,jcm3,kcm3)
      fdiffd = functn(ic,jcp4,kcp4) - functn(ic,jcp4,kcm4)  &
          - functn(ic,jcm4,kcp4) + functn(ic,jcm4,kcm4)
      fdiffe = functn(ic,jcp5,kcp5) - functn(ic,jcp5,kcm5)  &
          - functn(ic,jcm5,kcp5) + functn(ic,jcm5,kcm5)
      
      fderiv(ic,jc,kc) = acofyz*fdiffa + bcofyz*fdiffb  &
          + ccofyz*fdiffc + dcofyz*fdiffd  &
          + ecofyz*fdiffe
      
    END DO
    
  END DO
  
END DO

!     =========================================================================

!     LH END Y-DIRECTION
!     ==================
IF(nendyl == nbound)THEN
  
!       TAKE SECOND YZ-DERIVATIVE IN Y-LEFT INNER HALO
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
    
    DO ic = istal,istol
      
!           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofz1*(functn(ic,jstap1,kcp1) - functn(ic,jstap1,kcm1)  &
          - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
          + bcofz1*(functn(ic,jstap1,kcp2) - functn(ic,jstap1,kcm2)  &
          - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
      fdiffb = acofz1*(functn(ic,jstap2,kcp1) - functn(ic,jstap2,kcm1)  &
          - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
          + bcofz1*(functn(ic,jstap2,kcp2) - functn(ic,jstap2,kcm2)  &
          - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
      fdiffc = acofz1*(functn(ic,jstap3,kcp1) - functn(ic,jstap3,kcm1)  &
          - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
          + bcofz1*(functn(ic,jstap3,kcp2) - functn(ic,jstap3,kcm2)  &
          - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
      fdiffd = acofz1*(functn(ic,jstap4,kcp1) - functn(ic,jstap4,kcm1)  &
          - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
          + bcofz1*(functn(ic,jstap4,kcp2) - functn(ic,jstap4,kcm2)  &
          - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
      fderiv(ic,jstal,kc) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofz1*(functn(ic,jstal,kcp1)  - functn(ic,jstal,kcm1)  &
          - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
          + bcofz1*(functn(ic,jstal,kcp2)  - functn(ic,jstal,kcm2)  &
          - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
      fdiffb = acofz1*(functn(ic,jstap2,kcp1) - functn(ic,jstap2,kcm1)  &
          - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
          + bcofz1*(functn(ic,jstap2,kcp2) - functn(ic,jstap2,kcm2)  &
          - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
      fdiffc = acofz1*(functn(ic,jstap3,kcp1) - functn(ic,jstap3,kcm1)  &
          - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
          + bcofz1*(functn(ic,jstap3,kcp2) - functn(ic,jstap3,kcm2)  &
          - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
      fdiffd = acofz1*(functn(ic,jstap4,kcp1) - functn(ic,jstap4,kcm1)  &
          - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
          + bcofz1*(functn(ic,jstap4,kcp2) - functn(ic,jstap4,kcm2)  &
          - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
      fderiv(ic,jstap1,kc) = acf2yz*fdiffa + bcf2yz*fdiffb  &
          + ccf2yz*fdiffc + dcf2yz*fdiffd
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdiffa = functn(ic,jstap3,kcp1) - functn(ic,jstap3,kcm1)  &
          - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1)
      fdiffb = functn(ic,jstap4,kcp2) - functn(ic,jstap4,kcm2)  &
          - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2)
      fderiv(ic,jstap2,kc) = acf3yz*fdiffa + bcf3yz*fdiffb
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdiffa = functn(ic,jstap4,kcp1) - functn(ic,jstap4,kcm1)  &
          - functn(ic,jstap2,kcp1) + functn(ic,jstap2,kcm1)
      fdiffb = functn(ic,jstap5,kcp2) - functn(ic,jstap5,kcm2)  &
          - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2)
      fdiffc = functn(ic,jstap6,kcp3) - functn(ic,jstap6,kcm3)  &
          - functn(ic,jstal,kcp3)  + functn(ic,jstal,kcm3)
      fderiv(ic,jstap3,kc) = acf4yz*fdiffa + bcf4yz*fdiffb  &
          + ccf4yz*fdiffc
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdiffa = functn(ic,jstap5,kcp1) - functn(ic,jstap5,kcm1)  &
          - functn(ic,jstap3,kcp1) + functn(ic,jstap3,kcm1)
      fdiffb = functn(ic,jstap6,kcp2) - functn(ic,jstap6,kcm2)  &
          - functn(ic,jstap2,kcp2) + functn(ic,jstap2,kcm2)
      fdiffc = functn(ic,jstap7,kcp3) - functn(ic,jstap7,kcm3)  &
          - functn(ic,jstap1,kcp3) + functn(ic,jstap1,kcm3)
      fdiffd = functn(ic,jstap8,kcp4) - functn(ic,jstap8,kcm4)  &
          - functn(ic,jstal,kcp4)  + functn(ic,jstal,kcm4)
      fderiv(ic,jstap4,kc) = acf5yz*fdiffa + bcf5yz*fdiffb  &
          + ccf5yz*fdiffc + dcf5yz*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     RH END Y-DIRECTION
!     ==================
IF(nendyr == nbound)THEN
  
!       TAKE SECOND YZ-DERIVATIVE IN Y-RIGHT INNER HALO
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
    
    DO ic = istal,istol
      
!           RH POINT MINUS 4: 8TH ORDER CENTRED
      fdiffa = functn(ic,jstom3,kcp1) - functn(ic,jstom3,kcm1)  &
          - functn(ic,jstom5,kcp1) + functn(ic,jstom5,kcm1)
      fdiffb = functn(ic,jstom2,kcp2) - functn(ic,jstom2,kcm2)  &
          - functn(ic,jstom6,kcp2) + functn(ic,jstom6,kcm2)
      fdiffc = functn(ic,jstom1,kcp3) - functn(ic,jstom1,kcm3)  &
          - functn(ic,jstom7,kcp3) + functn(ic,jstom7,kcm3)
      fdiffd = functn(ic,jstol,kcp4)  - functn(ic,jstol,kcm4)  &
          - functn(ic,jstom8,kcp4) + functn(ic,jstom8,kcm4)
      fderiv(ic,jstom4,kc) = acf5yz*fdiffa + bcf5yz*fdiffb  &
          + ccf5yz*fdiffc + dcf5yz*fdiffd
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdiffa = functn(ic,jstom2,kcp1) - functn(ic,jstom2,kcm1)  &
          - functn(ic,jstom4,kcp1) + functn(ic,jstom4,kcm1)
      fdiffb = functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
          - functn(ic,jstom5,kcp2) + functn(ic,jstom5,kcm2)
      fdiffc = functn(ic,jstol,kcp3)  - functn(ic,jstol,kcm3)  &
          - functn(ic,jstom6,kcp3) + functn(ic,jstom6,kcm3)
      fderiv(ic,jstom3,kc) = acf4yz*fdiffa + bcf4yz*fdiffb  &
          + ccf4yz*fdiffc
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdiffa = functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
          - functn(ic,jstom3,kcp1) + functn(ic,jstom3,kcm1)
      fdiffb = functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
          - functn(ic,jstom4,kcp2) + functn(ic,jstom4,kcm2)
      fderiv(ic,jstom2,kc) = acf3yz*fdiffa + bcf3yz*fdiffb
      
!           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
          - functn(ic,jstol,kcp1)  + functn(ic,jstol,kcm1))  &
          + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
          - functn(ic,jstol,kcp2)  + functn(ic,jstol,kcm2))
      fdiffb = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
          - functn(ic,jstom2,kcp1) + functn(ic,jstom2,kcm1))  &
          + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
          - functn(ic,jstom2,kcp2) + functn(ic,jstom2,kcm2))
      fdiffc = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
          - functn(ic,jstom3,kcp1) + functn(ic,jstom3,kcm1))  &
          + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
          - functn(ic,jstom3,kcp2) + functn(ic,jstom3,kcm2))
      fdiffd = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
          - functn(ic,jstom4,kcp1) + functn(ic,jstom4,kcm1))  &
          + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
          - functn(ic,jstom4,kcp2) + functn(ic,jstom4,kcm2))
      fderiv(ic,jstom1,kc) = acf2yz*fdiffa + bcf2yz*fdiffb  &
          + ccf2yz*fdiffc + dcf2yz*fdiffd
      
!           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
          - functn(ic,jstom1,kcp1) + functn(ic,jstom1,kcm1))  &
          + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
          - functn(ic,jstom1,kcp2) + functn(ic,jstom1,kcm2))
      fdiffb = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
          - functn(ic,jstom2,kcp1) + functn(ic,jstom2,kcm1))  &
          + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
          - functn(ic,jstom2,kcp2) + functn(ic,jstom2,kcm2))
      fdiffc = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
          - functn(ic,jstom3,kcp1) + functn(ic,jstom3,kcm1))  &
          + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
          - functn(ic,jstom3,kcp2) + functn(ic,jstom3,kcm2))
      fdiffd = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
          - functn(ic,jstom4,kcp1) + functn(ic,jstom4,kcm1))  &
          + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
          - functn(ic,jstom4,kcp2) + functn(ic,jstom4,kcm2))
      fderiv(ic,jstol,kc) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     LH END Z-DIRECTION
!     ==================
IF(nendzl == nbound)THEN
  
!       TAKE SECOND XZ-DERIVATIVE IN Z-LEFT INNER HALO
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  jcm3 = jstart-4
  jcm2 = jstart-3
  jcm1 = jstart-2
  jccc = jstart-1
  jcp1 = jstart
  jcp2 = jstart+1
  jcp3 = jstart+2
  jcp4 = jstart+3
  
  DO jc = jstart,jfinis
    
    jcm4 = jcm3
    jcm3 = jcm2
    jcm2 = jcm1
    jcm1 = jccc
    jccc = jcp1
    jcp1 = jcp2
    jcp2 = jcp3
    jcp3 = jcp4
    jcp4 = jc+4
    
    DO ic = istal,istol
      
!           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofy1*(functn(ic,jcp1,kstap1) - functn(ic,jcm1,kstap1)  &
          - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
          + bcofy1*(functn(ic,jcp2,kstap1) - functn(ic,jcm2,kstap1)  &
          - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
      fdiffb = acofy1*(functn(ic,jcp1,kstap2) - functn(ic,jcm1,kstap2)  &
          - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
          + bcofy1*(functn(ic,jcp2,kstap2) - functn(ic,jcm2,kstap2)  &
          - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
      fdiffc = acofy1*(functn(ic,jcp1,kstap3) - functn(ic,jcm1,kstap3)  &
          - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
          + bcofy1*(functn(ic,jcp2,kstap3) - functn(ic,jcm2,kstap3)  &
          - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
      fdiffd = acofy1*(functn(ic,jcp1,kstap4) - functn(ic,jcm1,kstap4)  &
          - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
          + bcofy1*(functn(ic,jcp2,kstap4) - functn(ic,jcm2,kstap4)  &
          - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
      fderiv(ic,jc,kstal) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofy1*(functn(ic,jcp1,kstal)  - functn(ic,jcm1,kstal)  &
          - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
          + bcofy1*(functn(ic,jcp2,kstal)  - functn(ic,jcm2,kstal)  &
          - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
      fdiffb = acofy1*(functn(ic,jcp1,kstap2) - functn(ic,jcm1,kstap2)  &
          - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
          + bcofy1*(functn(ic,jcp2,kstap2) - functn(ic,jcm2,kstap2)  &
          - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
      fdiffc = acofy1*(functn(ic,jcp1,kstap3) - functn(ic,jcm1,kstap3)  &
          - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
          + bcofy1*(functn(ic,jcp2,kstap3) - functn(ic,jcm2,kstap3)  &
          - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
      fdiffd = acofy1*(functn(ic,jcp1,kstap4) - functn(ic,jcm1,kstap4)  &
          - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
          + bcofy1*(functn(ic,jcp2,kstap4) - functn(ic,jcm2,kstap4)  &
          - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
      fderiv(ic,jc,kstap1) = acf2yz*fdiffa + bcf2yz*fdiffb  &
          + ccf2yz*fdiffc + dcf2yz*fdiffd
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdiffa = functn(ic,jcp1,kstap3) - functn(ic,jcm1,kstap3)  &
          - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1)
      fdiffb = functn(ic,jcp2,kstap4) - functn(ic,jcm2,kstap4)  &
          - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal)
      fderiv(ic,jc,kstap2) = acf3yz*fdiffa + bcf3yz*fdiffb
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdiffa = functn(ic,jcp1,kstap4) - functn(ic,jcm1,kstap4)  &
          - functn(ic,jcp1,kstap2) + functn(ic,jcm1,kstap2)
      fdiffb = functn(ic,jcp2,kstap5) - functn(ic,jcm2,kstap5)  &
          - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1)
      fdiffc = functn(ic,jcp3,kstap6) - functn(ic,jcm3,kstap6)  &
          - functn(ic,jcp3,kstal)  + functn(ic,jcm3,kstal)
      fderiv(ic,jc,kstap3) = acf4yz*fdiffa + bcf4yz*fdiffb  &
          + ccf4yz*fdiffc
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdiffa = functn(ic,jcp1,kstap5) - functn(ic,jcm1,kstap5)  &
          - functn(ic,jcp1,kstap3) + functn(ic,jcm1,kstap3)
      fdiffb = functn(ic,jcp2,kstap6) - functn(ic,jcm2,kstap6)  &
          - functn(ic,jcp2,kstap2) + functn(ic,jcm2,kstap2)
      fdiffc = functn(ic,jcp3,kstap7) - functn(ic,jcm3,kstap7)  &
          - functn(ic,jcp3,kstap1) + functn(ic,jcm3,kstap1)
      fdiffd = functn(ic,jcp4,kstap8) - functn(ic,jcm4,kstap8)  &
          - functn(ic,jcp4,kstal)  + functn(ic,jcm4,kstal)
      fderiv(ic,jc,kstap4) = acf5yz*fdiffa + bcf5yz*fdiffb  &
          + ccf5yz*fdiffc + dcf5yz*fdiffd
      
    END DO
  END DO
  
!       LH IN Y LH IN Z CORNER
!       ======================
  IF(nendyl == nbound)THEN
    
    DO ic = istal,istol
      
!           LH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(ic,jstap1,kstap1) - functn(ic,jstap1,kstal)  &
          - functn(ic,jstal,kstap1)  + functn(ic,jstal,kstal)
      fdiffb = functn(ic,jstap2,kstap2) - functn(ic,jstap2,kstal)  &
          - functn(ic,jstal,kstap2)  + functn(ic,jstal,kstal)
      fdiffc = functn(ic,jstap3,kstap3) - functn(ic,jstap3,kstal)  &
          - functn(ic,jstal,kstap3)  + functn(ic,jstal,kstal)
      fdiffd = functn(ic,jstap4,kstap4) - functn(ic,jstap4,kstal)  &
          - functn(ic,jstal,kstap4)  + functn(ic,jstal,kstal)
      fderiv(ic,jstal,kstal) = acc1yz*fdiffa + bcc1yz*fdiffb  &
          + ccc1yz*fdiffc + dcc1yz*fdiffd
      
!           LH+1 LH+1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(ic,jstal,kstal)   - functn(ic,jstal,kstap1)  &
          - functn(ic,jstap1,kstal)  + functn(ic,jstap1,kstap1)
      fdiffb = functn(ic,jstap2,kstap2) - functn(ic,jstap2,kstap1)  &
          - functn(ic,jstap1,kstap2) + functn(ic,jstap1,kstap1)
      fdiffc = functn(ic,jstap3,kstap3) - functn(ic,jstap3,kstap1)  &
          - functn(ic,jstap1,kstap3) + functn(ic,jstap1,kstap1)
      fdiffd = functn(ic,jstap4,kstap4) - functn(ic,jstap4,kstap1)  &
          - functn(ic,jstap1,kstap4) + functn(ic,jstap1,kstap1)
      fderiv(ic,jstap1,kstap1) = acc2yz*fdiffa + bcc2yz*fdiffb  &
          + ccc2yz*fdiffc + dcc2yz*fdiffd
      
!           LH LH+1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2yz*(functn(ic,jstap1,kstal)  - functn(ic,jstap1,kstap1)  &
          - functn(ic,jstal,kstal)   + functn(ic,jstal,kstap1))  &
          + bcf2yz*(functn(ic,jstap1,kstap2) - functn(ic,jstap1,kstap1)  &
          - functn(ic,jstal,kstap2)  + functn(ic,jstal,kstap1))  &
          + ccf2yz*(functn(ic,jstap1,kstap3) - functn(ic,jstap1,kstap1)  &
          - functn(ic,jstal,kstap3)  + functn(ic,jstal,kstap1))  &
          + dcf2yz*(functn(ic,jstap1,kstap4) - functn(ic,jstap1,kstap1)  &
          - functn(ic,jstal,kstap4)  + functn(ic,jstal,kstap1))
      fdiffb = acf2yz*(functn(ic,jstap2,kstal)  - functn(ic,jstap2,kstap1)  &
          - functn(ic,jstal,kstal)   + functn(ic,jstal,kstap1))  &
          + bcf2yz*(functn(ic,jstap2,kstap2) - functn(ic,jstap2,kstap1)  &
          - functn(ic,jstal,kstap2)  + functn(ic,jstal,kstap1))  &
          + ccf2yz*(functn(ic,jstap2,kstap3) - functn(ic,jstap2,kstap1)  &
          - functn(ic,jstal,kstap3)  + functn(ic,jstal,kstap1))  &
          + dcf2yz*(functn(ic,jstap2,kstap4) - functn(ic,jstap2,kstap1)  &
          - functn(ic,jstal,kstap4)  + functn(ic,jstal,kstap1))
      fdiffc = acf2yz*(functn(ic,jstap3,kstal)  - functn(ic,jstap3,kstap1)  &
          - functn(ic,jstal,kstal)   + functn(ic,jstal,kstap1))  &
          + bcf2yz*(functn(ic,jstap3,kstap2) - functn(ic,jstap3,kstap1)  &
          - functn(ic,jstal,kstap2)  + functn(ic,jstal,kstap1))  &
          + ccf2yz*(functn(ic,jstap3,kstap3) - functn(ic,jstap3,kstap1)  &
          - functn(ic,jstal,kstap3)  + functn(ic,jstal,kstap1))  &
          + dcf2yz*(functn(ic,jstap3,kstap4) - functn(ic,jstap3,kstap1)  &
          - functn(ic,jstal,kstap4)  + functn(ic,jstal,kstap1))
      fdiffd = acf2yz*(functn(ic,jstap4,kstal)  - functn(ic,jstap4,kstap1)  &
          - functn(ic,jstal,kstal)   + functn(ic,jstal,kstap1))  &
          + bcf2yz*(functn(ic,jstap4,kstap2) - functn(ic,jstap4,kstap1)  &
          - functn(ic,jstal,kstap2)  + functn(ic,jstal,kstap1))  &
          + ccf2yz*(functn(ic,jstap4,kstap3) - functn(ic,jstap4,kstap1)  &
          - functn(ic,jstal,kstap3)  + functn(ic,jstal,kstap1))  &
          + dcf2yz*(functn(ic,jstap4,kstap4) - functn(ic,jstap4,kstap1)  &
          - functn(ic,jstal,kstap4)  + functn(ic,jstal,kstap1))
      fderiv(ic,jstal,kstap1) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           LH+1 LH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2yz*(functn(ic,jstal,kstap1)  - functn(ic,jstap1,kstap1)  &
          - functn(ic,jstal,kstal)   + functn(ic,jstap1,kstal))  &
          + bcf2yz*(functn(ic,jstap2,kstap1) - functn(ic,jstap1,kstap1)  &
          - functn(ic,jstap2,kstal)  + functn(ic,jstap1,kstal))  &
          + ccf2yz*(functn(ic,jstap3,kstap1) - functn(ic,jstap1,kstap1)  &
          - functn(ic,jstap3,kstal)  + functn(ic,jstap1,kstal))  &
          + dcf2yz*(functn(ic,jstap4,kstap1) - functn(ic,jstap1,kstap1)  &
          - functn(ic,jstap4,kstal)  + functn(ic,jstap1,kstal))
      fdiffb = acf2yz*(functn(ic,jstal,kstap2)  - functn(ic,jstap1,kstap2)  &
          - functn(ic,jstal,kstal)   + functn(ic,jstap1,kstal))  &
          + bcf2yz*(functn(ic,jstap2,kstap2) - functn(ic,jstap1,kstap2)  &
          - functn(ic,jstap2,kstal)  + functn(ic,jstap1,kstal))  &
          + ccf2yz*(functn(ic,jstap3,kstap2) - functn(ic,jstap1,kstap2)  &
          - functn(ic,jstap3,kstal)  + functn(ic,jstap1,kstal))  &
          + dcf2yz*(functn(ic,jstap4,kstap2) - functn(ic,jstap1,kstap2)  &
          - functn(ic,jstap4,kstal)  + functn(ic,jstap1,kstal))
      fdiffc = acf2yz*(functn(ic,jstal,kstap3)  - functn(ic,jstap1,kstap3)  &
          - functn(ic,jstal,kstal)   + functn(ic,jstap1,kstal))  &
          + bcf2yz*(functn(ic,jstap2,kstap3) - functn(ic,jstap1,kstap3)  &
          - functn(ic,jstap2,kstal)  + functn(ic,jstap1,kstal))  &
          + ccf2yz*(functn(ic,jstap3,kstap3) - functn(ic,jstap1,kstap3)  &
          - functn(ic,jstap3,kstal)  + functn(ic,jstap1,kstal))  &
          + dcf2yz*(functn(ic,jstap4,kstap3) - functn(ic,jstap1,kstap3)  &
          - functn(ic,jstap4,kstal)  + functn(ic,jstap1,kstal))
      fdiffd = acf2yz*(functn(ic,jstal,kstap4)  - functn(ic,jstap1,kstap4)  &
          - functn(ic,jstal,kstal)   + functn(ic,jstap1,kstal))  &
          + bcf2yz*(functn(ic,jstap2,kstap4) - functn(ic,jstap1,kstap4)  &
          - functn(ic,jstap2,kstal)  + functn(ic,jstap1,kstal))  &
          + ccf2yz*(functn(ic,jstap3,kstap4) - functn(ic,jstap1,kstap4)  &
          - functn(ic,jstap3,kstal)  + functn(ic,jstap1,kstal))  &
          + dcf2yz*(functn(ic,jstap4,kstap4) - functn(ic,jstap1,kstap4)  &
          - functn(ic,jstap4,kstal)  + functn(ic,jstap1,kstal))
      fderiv(ic,jstap1,kstal) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           LH EDGE IN Z
      DO jc = jstap2,jstap4
        
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofy1*(functn(ic,jcp1,kstap1) - functn(ic,jcm1,kstap1)  &
            - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
            + bcofy1*(functn(ic,jcp2,kstap1) - functn(ic,jcm2,kstap1)  &
            - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
        fdiffb = acofy1*(functn(ic,jcp1,kstap2) - functn(ic,jcm1,kstap2)  &
            - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
            + bcofy1*(functn(ic,jcp2,kstap2) - functn(ic,jcm2,kstap2)  &
            - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
        fdiffc = acofy1*(functn(ic,jcp1,kstap3) - functn(ic,jcm1,kstap3)  &
            - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
            + bcofy1*(functn(ic,jcp2,kstap3) - functn(ic,jcm2,kstap3)  &
            - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
        fdiffd = acofy1*(functn(ic,jcp1,kstap4) - functn(ic,jcm1,kstap4)  &
            - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
            + bcofy1*(functn(ic,jcp2,kstap4) - functn(ic,jcm2,kstap4)  &
            - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
        fderiv(ic,jc,kstal) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofy1*(functn(ic,jcp1,kstal)  - functn(ic,jcm1,kstal)  &
            - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
            + bcofy1*(functn(ic,jcp2,kstal)  - functn(ic,jcm2,kstal)  &
            - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
        fdiffb = acofy1*(functn(ic,jcp1,kstap2) - functn(ic,jcm1,kstap2)  &
            - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
            + bcofy1*(functn(ic,jcp2,kstap2) - functn(ic,jcm2,kstap2)  &
            - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
        fdiffc = acofy1*(functn(ic,jcp1,kstap3) - functn(ic,jcm1,kstap3)  &
            - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
            + bcofy1*(functn(ic,jcp2,kstap3) - functn(ic,jcm2,kstap3)  &
            - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
        fdiffd = acofy1*(functn(ic,jcp1,kstap4) - functn(ic,jcm1,kstap4)  &
            - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
            + bcofy1*(functn(ic,jcp2,kstap4) - functn(ic,jcm2,kstap4)  &
            - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
        fderiv(ic,jc,kstap1) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd
        
      END DO
      
!           LH EDGE IN Y
      DO kc = kstap2,kstap4
        
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofz1*(functn(ic,jstap1,kcp1) - functn(ic,jstap1,kcm1)  &
            - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
            + bcofz1*(functn(ic,jstap1,kcp2) - functn(ic,jstap1,kcm2)  &
            - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
        fdiffb = acofz1*(functn(ic,jstap2,kcp1) - functn(ic,jstap2,kcm1)  &
            - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
            + bcofz1*(functn(ic,jstap2,kcp2) - functn(ic,jstap2,kcm2)  &
            - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
        fdiffc = acofz1*(functn(ic,jstap3,kcp1) - functn(ic,jstap3,kcm1)  &
            - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
            + bcofz1*(functn(ic,jstap3,kcp2) - functn(ic,jstap3,kcm2)  &
            - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
        fdiffd = acofz1*(functn(ic,jstap4,kcp1) - functn(ic,jstap4,kcm1)  &
            - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
            + bcofz1*(functn(ic,jstap4,kcp2) - functn(ic,jstap4,kcm2)  &
            - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
        fderiv(ic,jstal,kc) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofz1*(functn(ic,jstal,kcp1)  - functn(ic,jstal,kcm1)  &
            - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
            + bcofz1*(functn(ic,jstal,kcp2)  - functn(ic,jstal,kcm2)  &
            - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
        fdiffb = acofz1*(functn(ic,jstap2,kcp1) - functn(ic,jstap2,kcm1)  &
            - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
            + bcofz1*(functn(ic,jstap2,kcp2) - functn(ic,jstap2,kcm2)  &
            - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
        fdiffc = acofz1*(functn(ic,jstap3,kcp1) - functn(ic,jstap3,kcm1)  &
            - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
            + bcofz1*(functn(ic,jstap3,kcp2) - functn(ic,jstap3,kcm2)  &
            - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
        fdiffd = acofz1*(functn(ic,jstap4,kcp1) - functn(ic,jstap4,kcm1)  &
            - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
            + bcofz1*(functn(ic,jstap4,kcp2) - functn(ic,jstap4,kcm2)  &
            - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
        fderiv(ic,jstap1,kc) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      ks = 0
      DO kc = kstap2,kstap4
        
        ks = ks+1
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
        js = 0
        DO jc = jstap2,jstap4
          
          js = js+1
          jcm2 = jc-2
          jcm1 = jc-1
          jcp1 = jc+1
          jcp2 = jc+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(ic,jcp1,kcp1) - functn(ic,jcp1,kcm1)  &
              - functn(ic,jcm1,kcp1) + functn(ic,jcm1,kcm1)
          fdiffb = functn(ic,jcp2,kcp2) - functn(ic,jcp2,kcm2)  &
              - functn(ic,jcm2,kcp2) + functn(ic,jcm2,kcm2)
          fderiv(ic,jc,kc) = acf3yz*fdiffa + bcf3yz*fdiffb
          fstora(js,ks) = fdiffa
          fstorb(js,ks) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      ks = 1
      DO kc = kstap3,kstap4
        
        ksm1 = ks
        ks = ks+1
        kcm3 = kc-3
        kcp3 = kc+3
        
        js = 1
        DO jc = jstap3,jstap4
          
          jsm1 = js
          js = js+1
          jcm3 = jc-3
          jcp3 = jc+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(ic,jcp3,kcp3) - functn(ic,jcp3,kcm3)  &
              - functn(ic,jcm3,kcp3) + functn(ic,jcm3,kcm3)
          fderiv(ic,jc,kc) = acf4yz*fstora(js,ks) + bcf4yz*fstorb(js,ks)  &
              + ccf4yz*fdiffc
          fstorc(jsm1,ksm1) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      ks = 3
      js = 3
      ksm1 = 2
      jsm1 = 2
      kc = kstap4
      jc = jstap4
      kcm4 = kc-4
      kcp4 = kc+4
      jcm4 = jc-4
      jcp4 = jc+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(ic,jcp4,kcp4) - functn(ic,jcp4,kcm4)  &
          - functn(ic,jcm4,kcp4) + functn(ic,jcm4,kcm4)
      fderiv(ic,jc,kc) = acf5yz*fstora(js,ks) + bcf5yz*fstorb(js,ks)  &
          + ccf5yz*fstorc(jsm1,ksm1) + dcf5yz*fdiffd
      
    END DO
    
  END IF
  
!       RH IN Y LH IN Z CORNER
!       ======================
  IF(nendyr == nbound)THEN
    
    DO ic = istal,istol
      
!           RH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(ic,jstol,kstap1)  - functn(ic,jstol,kstal)  &
          - functn(ic,jstom1,kstap1) + functn(ic,jstom1,kstal)
      fdiffb = functn(ic,jstol,kstap2)  - functn(ic,jstol,kstal)  &
          - functn(ic,jstom2,kstap2) + functn(ic,jstom2,kstal)
      fdiffc = functn(ic,jstol,kstap3)  - functn(ic,jstol,kstal)  &
          - functn(ic,jstom3,kstap3) + functn(ic,jstom3,kstal)
      fdiffd = functn(ic,jstol,kstap4)  - functn(ic,jstol,kstal)  &
          - functn(ic,jstom4,kstap4) + functn(ic,jstom4,kstal)
      fderiv(ic,jstol,kstal) = acc1yz*fdiffa + bcc1yz*fdiffb  &
          + ccc1yz*fdiffc + dcc1yz*fdiffd
      
!           RH-1 LH+1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(ic,jstom1,kstal)  - functn(ic,jstom1,kstap1)  &
          - functn(ic,jstol,kstal)   + functn(ic,jstol,kstap1)
      fdiffb = functn(ic,jstom1,kstap2) - functn(ic,jstom1,kstap1)  &
          - functn(ic,jstom2,kstap2) + functn(ic,jstom2,kstap1)
      fdiffc = functn(ic,jstom1,kstap3) - functn(ic,jstom1,kstap1)  &
          - functn(ic,jstom3,kstap3) + functn(ic,jstom3,kstap1)
      fdiffd = functn(ic,jstom1,kstap4) - functn(ic,jstom1,kstap1)  &
          - functn(ic,jstom4,kstap4) + functn(ic,jstom4,kstap1)
      fderiv(ic,jstom1,kstap1) = acc2yz*fdiffa + bcc2yz*fdiffb  &
          + ccc2yz*fdiffc + dcc2yz*fdiffd
      
!           RH LH+1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2yz*(functn(ic,jstol,kstal)   - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom1,kstap1))  &
          + bcf2yz*(functn(ic,jstol,kstap2)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom1,kstap2) + functn(ic,jstom1,kstap1))  &
          + ccf2yz*(functn(ic,jstol,kstap3)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom1,kstap3) + functn(ic,jstom1,kstap1))  &
          + dcf2yz*(functn(ic,jstol,kstap4)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom1,kstap4) + functn(ic,jstom1,kstap1))
      fdiffb = acf2yz*(functn(ic,jstol,kstal)   - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom2,kstal)  + functn(ic,jstom2,kstap1))  &
          + bcf2yz*(functn(ic,jstol,kstap2)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom2,kstap2) + functn(ic,jstom2,kstap1))  &
          + ccf2yz*(functn(ic,jstol,kstap3)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom2,kstap3) + functn(ic,jstom2,kstap1))  &
          + dcf2yz*(functn(ic,jstol,kstap4)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom2,kstap4) + functn(ic,jstom2,kstap1))
      fdiffc = acf2yz*(functn(ic,jstol,kstal)   - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom3,kstal)  + functn(ic,jstom3,kstap1))  &
          + bcf2yz*(functn(ic,jstol,kstap2)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom3,kstap2) + functn(ic,jstom3,kstap1))  &
          + ccf2yz*(functn(ic,jstol,kstap3)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom3,kstap3) + functn(ic,jstom3,kstap1))  &
          + dcf2yz*(functn(ic,jstol,kstap4)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom3,kstap4) + functn(ic,jstom3,kstap1))
      fdiffd = acf2yz*(functn(ic,jstol,kstal)   - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom4,kstal)  + functn(ic,jstom4,kstap1))  &
          + bcf2yz*(functn(ic,jstol,kstap2)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom4,kstap2) + functn(ic,jstom4,kstap1))  &
          + ccf2yz*(functn(ic,jstol,kstap3)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom4,kstap3) + functn(ic,jstom4,kstap1))  &
          + dcf2yz*(functn(ic,jstol,kstap4)  - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom4,kstap4) + functn(ic,jstom4,kstap1))
      fderiv(ic,jstol,kstap1) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           RH-1 LH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2yz*(functn(ic,jstom1,kstap1) - functn(ic,jstol,kstap1)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstol,kstal))  &
          + bcf2yz*(functn(ic,jstom1,kstap1) - functn(ic,jstom2,kstap1)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom2,kstal))  &
          + ccf2yz*(functn(ic,jstom1,kstap1) - functn(ic,jstom3,kstap1)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom3,kstal))  &
          + dcf2yz*(functn(ic,jstom1,kstap1) - functn(ic,jstom4,kstap1)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom4,kstal))
      fdiffb = acf2yz*(functn(ic,jstom1,kstap2) - functn(ic,jstol,kstap2)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstol,kstal))  &
          + bcf2yz*(functn(ic,jstom1,kstap2) - functn(ic,jstom2,kstap2)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom2,kstal))  &
          + ccf2yz*(functn(ic,jstom1,kstap2) - functn(ic,jstom3,kstap2)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom3,kstal))  &
          + dcf2yz*(functn(ic,jstom1,kstap2) - functn(ic,jstom4,kstap2)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom4,kstal))
      fdiffc = acf2yz*(functn(ic,jstom1,kstap3) - functn(ic,jstol,kstap3)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstol,kstal))  &
          + bcf2yz*(functn(ic,jstom1,kstap3) - functn(ic,jstom2,kstap3)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom2,kstal))  &
          + ccf2yz*(functn(ic,jstom1,kstap3) - functn(ic,jstom3,kstap3)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom3,kstal))  &
          + dcf2yz*(functn(ic,jstom1,kstap3) - functn(ic,jstom4,kstap3)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom4,kstal))
      fdiffd = acf2yz*(functn(ic,jstom1,kstap4) - functn(ic,jstol,kstap4)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstol,kstal))  &
          + bcf2yz*(functn(ic,jstom1,kstap4) - functn(ic,jstom2,kstap4)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom2,kstal))  &
          + ccf2yz*(functn(ic,jstom1,kstap4) - functn(ic,jstom3,kstap4)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom3,kstal))  &
          + dcf2yz*(functn(ic,jstom1,kstap4) - functn(ic,jstom4,kstap4)  &
          - functn(ic,jstom1,kstal)  + functn(ic,jstom4,kstal))
      fderiv(ic,jstom1,kstal) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           LH EDGE IN Z
      DO jc = jstom4,jstom2
        
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofy1*(functn(ic,jcp1,kstap1) - functn(ic,jcm1,kstap1)  &
            - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
            + bcofy1*(functn(ic,jcp2,kstap1) - functn(ic,jcm2,kstap1)  &
            - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
        fdiffb = acofy1*(functn(ic,jcp1,kstap2) - functn(ic,jcm1,kstap2)  &
            - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
            + bcofy1*(functn(ic,jcp2,kstap2) - functn(ic,jcm2,kstap2)  &
            - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
        fdiffc = acofy1*(functn(ic,jcp1,kstap3) - functn(ic,jcm1,kstap3)  &
            - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
            + bcofy1*(functn(ic,jcp2,kstap3) - functn(ic,jcm2,kstap3)  &
            - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
        fdiffd = acofy1*(functn(ic,jcp1,kstap4) - functn(ic,jcm1,kstap4)  &
            - functn(ic,jcp1,kstal)  + functn(ic,jcm1,kstal))  &
            + bcofy1*(functn(ic,jcp2,kstap4) - functn(ic,jcm2,kstap4)  &
            - functn(ic,jcp2,kstal)  + functn(ic,jcm2,kstal))
        fderiv(ic,jc,kstal) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofy1*(functn(ic,jcp1,kstal)  - functn(ic,jcm1,kstal)  &
            - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
            + bcofy1*(functn(ic,jcp2,kstal)  - functn(ic,jcm2,kstal)  &
            - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
        fdiffb = acofy1*(functn(ic,jcp1,kstap2) - functn(ic,jcm1,kstap2)  &
            - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
            + bcofy1*(functn(ic,jcp2,kstap2) - functn(ic,jcm2,kstap2)  &
            - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
        fdiffc = acofy1*(functn(ic,jcp1,kstap3) - functn(ic,jcm1,kstap3)  &
            - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
            + bcofy1*(functn(ic,jcp2,kstap3) - functn(ic,jcm2,kstap3)  &
            - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
        fdiffd = acofy1*(functn(ic,jcp1,kstap4) - functn(ic,jcm1,kstap4)  &
            - functn(ic,jcp1,kstap1) + functn(ic,jcm1,kstap1))  &
            + bcofy1*(functn(ic,jcp2,kstap4) - functn(ic,jcm2,kstap4)  &
            - functn(ic,jcp2,kstap1) + functn(ic,jcm2,kstap1))
        fderiv(ic,jc,kstap1) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd
        
      END DO
      
!           RH EDGE IN Y
      DO kc = kstap2,kstap4
        
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
            - functn(ic,jstom1,kcp1) + functn(ic,jstom1,kcm1))  &
            + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
            - functn(ic,jstom1,kcp2) + functn(ic,jstom1,kcm2))
        fdiffb = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
            - functn(ic,jstom2,kcp1) + functn(ic,jstom2,kcm1))  &
            + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
            - functn(ic,jstom2,kcp2) + functn(ic,jstom2,kcm2))
        fdiffc = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
            - functn(ic,jstom3,kcp1) + functn(ic,jstom3,kcm1))  &
            + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
            - functn(ic,jstom3,kcp2) + functn(ic,jstom3,kcm2))
        fdiffd = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
            - functn(ic,jstom4,kcp1) + functn(ic,jstom4,kcm1))  &
            + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
            - functn(ic,jstom4,kcp2) + functn(ic,jstom4,kcm2))
        fderiv(ic,jstol,kc) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
            - functn(ic,jstol,kcp1)  + functn(ic,jstol,kcm1))  &
            + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
            - functn(ic,jstol,kcp2)  + functn(ic,jstol,kcm2))
        fdiffb = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
            - functn(ic,jstom2,kcp1) + functn(ic,jstom2,kcm1))  &
            + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
            - functn(ic,jstom2,kcp2) + functn(ic,jstom2,kcm2))
        fdiffc = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
            - functn(ic,jstom3,kcp1) + functn(ic,jstom3,kcm1))  &
            + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
            - functn(ic,jstom3,kcp2) + functn(ic,jstom3,kcm2))
        fdiffd = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
            - functn(ic,jstom4,kcp1) + functn(ic,jstom4,kcm1))  &
            + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
            - functn(ic,jstom4,kcp2) + functn(ic,jstom4,kcm2))
        fderiv(ic,jstom1,kc) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      ks = 0
      DO kc = kstap2,kstap4
        
        ks = ks+1
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
        js = 0
        DO jc = jstom4,jstom2
          
          js = js+1
          jcm2 = jc-2
          jcm1 = jc-1
          jcp1 = jc+1
          jcp2 = jc+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(ic,jcp1,kcp1) - functn(ic,jcp1,kcm1)  &
              - functn(ic,jcm1,kcp1) + functn(ic,jcm1,kcm1)
          fdiffb = functn(ic,jcp2,kcp2) - functn(ic,jcp2,kcm2)  &
              - functn(ic,jcm2,kcp2) + functn(ic,jcm2,kcm2)
          fderiv(ic,jc,kc) = acf3yz*fdiffa + bcf3yz*fdiffb
          fstora(js,ks) = fdiffa
          fstorb(js,ks) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      ks = 1
      DO kc = kstap3,kstap4
        
        ksm1 = ks
        ks = ks+1
        kcm3 = kc-3
        kcp3 = kc+3
        
        js = 0
        DO jc = jstom4,jstom3
          
          js = js+1
          jcm3 = jc-3
          jcp3 = jc+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(ic,jcp3,kcp3) - functn(ic,jcp3,kcm3)  &
              - functn(ic,jcm3,kcp3) + functn(ic,jcm3,kcm3)
          fderiv(ic,jc,kc) = acf4yz*fstora(js,ks) + bcf4yz*fstorb(js,ks)  &
              + ccf4yz*fdiffc
          fstorc(js,ksm1) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      ks = 3
      js = 1
      ksm1 = 2
      kc = kstap4
      jc = jstom4
      kcm4 = kc-4
      kcp4 = kc+4
      jcm4 = jc-4
      jcp4 = jc+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(ic,jcp4,kcp4) - functn(ic,jcp4,kcm4)  &
          - functn(ic,jcm4,kcp4) + functn(ic,jcm4,kcm4)
      fderiv(ic,jc,kc) = acf5yz*fstora(js,ks) + bcf5yz*fstorb(js,ks)  &
          + ccf5yz*fstorc(js,ksm1) + dcf5yz*fdiffd
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     RH END Z-DIRECTION
!     ==================
IF(nendzr == nbound)THEN
  
!       TAKE SECOND YZ-DERIVATIVE IN Z-RIGHT INNER HALO
!       EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  
  jcm3 = jstart-4
  jcm2 = jstart-3
  jcm1 = jstart-2
  jccc = jstart-1
  jcp1 = jstart
  jcp2 = jstart+1
  jcp3 = jstart+2
  jcp4 = jstart+3
  
  DO jc = jstart,jfinis
    
    jcm4 = jcm3
    jcm3 = jcm2
    jcm2 = jcm1
    jcm1 = jccc
    jccc = jcp1
    jcp1 = jcp2
    jcp2 = jcp3
    jcp3 = jcp4
    jcp4 = jc+4
    
    DO ic = istal,istol
      
!           RH POINT MINUS 4: 8TH ORDER CENTRED
      fdiffa = functn(ic,jcp1,kstom3) - functn(ic,jcm1,kstom3)  &
          - functn(ic,jcp1,kstom5) + functn(ic,jcm1,kstom5)
      fdiffb = functn(ic,jcp2,kstom2) - functn(ic,jcm2,kstom2)  &
          - functn(ic,jcp2,kstom6) + functn(ic,jcm2,kstom6)
      fdiffc = functn(ic,jcp3,kstom1) - functn(ic,jcm3,kstom1)  &
          - functn(ic,jcp3,kstom7) + functn(ic,jcm3,kstom7)
      fdiffd = functn(ic,jcp4,kstol)  - functn(ic,jcm4,kstol)  &
          - functn(ic,jcp4,kstom8) + functn(ic,jcm4,kstom8)
      fderiv(ic,jc,kstom4) = acf5yz*fdiffa + bcf5yz*fdiffb  &
          + ccf5yz*fdiffc + dcf5yz*fdiffd
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdiffa = functn(ic,jcp1,kstom2) - functn(ic,jcm1,kstom2)  &
          - functn(ic,jcp1,kstom4) + functn(ic,jcm1,kstom4)
      fdiffb = functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
          - functn(ic,jcp2,kstom5) + functn(ic,jcm2,kstom5)
      fdiffc = functn(ic,jcp3,kstol)  - functn(ic,jcm3,kstol)  &
          - functn(ic,jcp3,kstom6) + functn(ic,jcm3,kstom6)
      fderiv(ic,jc,kstom3) = acf4yz*fdiffa + bcf4yz*fdiffb  &
          + ccf4yz*fdiffc
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdiffa = functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
          - functn(ic,jcp1,kstom3) + functn(ic,jcm1,kstom3)
      fdiffb = functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
          - functn(ic,jcp2,kstom4) + functn(ic,jcm2,kstom4)
      fderiv(ic,jc,kstom2) = acf3yz*fdiffa + bcf3yz*fdiffb
      
!           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
          - functn(ic,jcp1,kstol)  + functn(ic,jcm1,kstol))  &
          + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
          - functn(ic,jcp2,kstol)  + functn(ic,jcm2,kstol))
      fdiffb = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
          - functn(ic,jcp1,kstom2) + functn(ic,jcm1,kstom2))  &
          + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
          - functn(ic,jcp2,kstom2) + functn(ic,jcm2,kstom2))
      fdiffc = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
          - functn(ic,jcp1,kstom3) + functn(ic,jcm1,kstom3))  &
          + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
          - functn(ic,jcp2,kstom3) + functn(ic,jcm2,kstom3))
      fdiffd = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
          - functn(ic,jcp1,kstom4) + functn(ic,jcm1,kstom4))  &
          + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
          - functn(ic,jcp2,kstom4) + functn(ic,jcm2,kstom4))
      fderiv(ic,jc,kstom1) = acf2yz*fdiffa + bcf2yz*fdiffb  &
          + ccf2yz*fdiffc + dcf2yz*fdiffd
      
!           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
          - functn(ic,jcp1,kstom1) + functn(ic,jcm1,kstom1))  &
          + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
          - functn(ic,jcp2,kstom1) + functn(ic,jcm2,kstom1))
      fdiffb = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
          - functn(ic,jcp1,kstom2) + functn(ic,jcm1,kstom2))  &
          + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
          - functn(ic,jcp2,kstom2) + functn(ic,jcm2,kstom2))
      fdiffc = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
          - functn(ic,jcp1,kstom3) + functn(ic,jcm1,kstom3))  &
          + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
          - functn(ic,jcp2,kstom3) + functn(ic,jcm2,kstom3))
      fdiffd = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
          - functn(ic,jcp1,kstom4) + functn(ic,jcm1,kstom4))  &
          + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
          - functn(ic,jcp2,kstom4) + functn(ic,jcm2,kstom4))
      fderiv(ic,jc,kstol) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
    END DO
  END DO
  
!       LH IN Y RH IN Z CORNER
!       ======================
  IF(nendyl == nbound)THEN
    
    DO ic = istal,istol
      
!           LH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(ic,jstap1,kstol) - functn(ic,jstap1,kstom1)  &
          - functn(ic,jstal,kstol)  + functn(ic,jstal,kstom1)
      fdiffb = functn(ic,jstap2,kstol) - functn(ic,jstap2,kstom2)  &
          - functn(ic,jstal,kstol)  + functn(ic,jstal,kstom2)
      fdiffc = functn(ic,jstap3,kstol) - functn(ic,jstap3,kstom3)  &
          - functn(ic,jstal,kstol)  + functn(ic,jstal,kstom3)
      fdiffd = functn(ic,jstap4,kstol) - functn(ic,jstap4,kstom4)  &
          - functn(ic,jstal,kstol)  + functn(ic,jstal,kstom4)
      fderiv(ic,jstal,kstol) = acc1yz*fdiffa + bcc1yz*fdiffb  &
          + ccc1yz*fdiffc + dcc1yz*fdiffd
      
!           LH+1 RH-1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(ic,jstal,kstom1)  - functn(ic,jstal,kstol)  &
          - functn(ic,jstap1,kstom1) + functn(ic,jstap1,kstol)
      fdiffb = functn(ic,jstap2,kstom1) - functn(ic,jstap2,kstom2)  &
          - functn(ic,jstap1,kstom1) + functn(ic,jstap1,kstom2)
      fdiffc = functn(ic,jstap3,kstom1) - functn(ic,jstap3,kstom3)  &
          - functn(ic,jstap1,kstom1) + functn(ic,jstap1,kstom3)
      fdiffd = functn(ic,jstap4,kstom1) - functn(ic,jstap4,kstom4)  &
          - functn(ic,jstap1,kstom1) + functn(ic,jstap1,kstom4)
      fderiv(ic,jstap1,kstom1) = acc2yz*fdiffa + bcc2yz*fdiffb  &
          + ccc2yz*fdiffc + dcc2yz*fdiffd
      
!           LH RH-1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2yz*(functn(ic,jstap1,kstom1) - functn(ic,jstap1,kstol)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstol))  &
          + bcf2yz*(functn(ic,jstap1,kstom1) - functn(ic,jstap1,kstom2)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom2))  &
          + ccf2yz*(functn(ic,jstap1,kstom1) - functn(ic,jstap1,kstom3)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom3))  &
          + dcf2yz*(functn(ic,jstap1,kstom1) - functn(ic,jstap1,kstom4)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom4))
      fdiffb = acf2yz*(functn(ic,jstap2,kstom1) - functn(ic,jstap2,kstol)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstol))  &
          + bcf2yz*(functn(ic,jstap2,kstom1) - functn(ic,jstap2,kstom2)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom2))  &
          + ccf2yz*(functn(ic,jstap2,kstom1) - functn(ic,jstap2,kstom3)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom3))  &
          + dcf2yz*(functn(ic,jstap2,kstom1) - functn(ic,jstap2,kstom4)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom4))
      fdiffc = acf2yz*(functn(ic,jstap3,kstom1) - functn(ic,jstap3,kstol)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstol))  &
          + bcf2yz*(functn(ic,jstap3,kstom1) - functn(ic,jstap3,kstom2)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom2))  &
          + ccf2yz*(functn(ic,jstap3,kstom1) - functn(ic,jstap3,kstom3)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom3))  &
          + dcf2yz*(functn(ic,jstap3,kstom1) - functn(ic,jstap3,kstom4)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom4))
      fdiffd = acf2yz*(functn(ic,jstap4,kstom1) - functn(ic,jstap4,kstol)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstol))  &
          + bcf2yz*(functn(ic,jstap4,kstom1) - functn(ic,jstap4,kstom2)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom2))  &
          + ccf2yz*(functn(ic,jstap4,kstom1) - functn(ic,jstap4,kstom3)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom3))  &
          + dcf2yz*(functn(ic,jstap4,kstom1) - functn(ic,jstap4,kstom4)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstal,kstom4))
      fderiv(ic,jstal,kstom1) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           LH+1 RH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2yz*(functn(ic,jstal,kstol)   - functn(ic,jstap1,kstol)  &
          - functn(ic,jstal,kstom1)  + functn(ic,jstap1,kstom1))  &
          + bcf2yz*(functn(ic,jstap2,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap2,kstom1) + functn(ic,jstap1,kstom1))  &
          + ccf2yz*(functn(ic,jstap3,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap3,kstom1) + functn(ic,jstap1,kstom1))  &
          + dcf2yz*(functn(ic,jstap4,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap4,kstom1) + functn(ic,jstap1,kstom1))
      fdiffb = acf2yz*(functn(ic,jstal,kstol)   - functn(ic,jstap1,kstol)  &
          - functn(ic,jstal,kstom2)  + functn(ic,jstap1,kstom2))  &
          + bcf2yz*(functn(ic,jstap2,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap2,kstom2) + functn(ic,jstap1,kstom2))  &
          + ccf2yz*(functn(ic,jstap3,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap3,kstom2) + functn(ic,jstap1,kstom2))  &
          + dcf2yz*(functn(ic,jstap4,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap4,kstom2) + functn(ic,jstap1,kstom2))
      fdiffc = acf2yz*(functn(ic,jstal,kstol)   - functn(ic,jstap1,kstol)  &
          - functn(ic,jstal,kstom3)  + functn(ic,jstap1,kstom3))  &
          + bcf2yz*(functn(ic,jstap2,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap2,kstom3) + functn(ic,jstap1,kstom3))  &
          + ccf2yz*(functn(ic,jstap3,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap3,kstom3) + functn(ic,jstap1,kstom3))  &
          + dcf2yz*(functn(ic,jstap4,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap4,kstom3) + functn(ic,jstap1,kstom3))
      fdiffd = acf2yz*(functn(ic,jstal,kstol)   - functn(ic,jstap1,kstol)  &
          - functn(ic,jstal,kstom4)  + functn(ic,jstap1,kstom4))  &
          + bcf2yz*(functn(ic,jstap2,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap2,kstom4) + functn(ic,jstap1,kstom4))  &
          + ccf2yz*(functn(ic,jstap3,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap3,kstom4) + functn(ic,jstap1,kstom4))  &
          + dcf2yz*(functn(ic,jstap4,kstol)  - functn(ic,jstap1,kstol)  &
          - functn(ic,jstap4,kstom4) + functn(ic,jstap1,kstom4))
      fderiv(ic,jstap1,kstol) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           RH EDGE IN Z
      DO jc = jstap2,jstap4
        
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
            - functn(ic,jcp1,kstol)  + functn(ic,jcm1,kstol))  &
            + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
            - functn(ic,jcp2,kstol)  + functn(ic,jcm2,kstol))
        fdiffb = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
            - functn(ic,jcp1,kstom2) + functn(ic,jcm1,kstom2))  &
            + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
            - functn(ic,jcp2,kstom2) + functn(ic,jcm2,kstom2))
        fdiffc = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
            - functn(ic,jcp1,kstom3) + functn(ic,jcm1,kstom3))  &
            + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
            - functn(ic,jcp2,kstom3) + functn(ic,jcm2,kstom3))
        fdiffd = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
            - functn(ic,jcp1,kstom4) + functn(ic,jcm1,kstom4))  &
            + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
            - functn(ic,jcp2,kstom4) + functn(ic,jcm2,kstom4))
        fderiv(ic,jc,kstom1) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
            - functn(ic,jcp1,kstom1) + functn(ic,jcm1,kstom1))  &
            + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
            - functn(ic,jcp2,kstom1) + functn(ic,jcm2,kstom1))
        fdiffb = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
            - functn(ic,jcp1,kstom2) + functn(ic,jcm1,kstom2))  &
            + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
            - functn(ic,jcp2,kstom2) + functn(ic,jcm2,kstom2))
        fdiffc = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
            - functn(ic,jcp1,kstom3) + functn(ic,jcm1,kstom3))  &
            + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
            - functn(ic,jcp2,kstom3) + functn(ic,jcm2,kstom3))
        fdiffd = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
            - functn(ic,jcp1,kstom4) + functn(ic,jcm1,kstom4))  &
            + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
            - functn(ic,jcp2,kstom4) + functn(ic,jcm2,kstom4))
        fderiv(ic,jc,kstol) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd
        
      END DO
      
!           LH EDGE IN Y
      DO kc = kstom4,kstom2
        
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofz1*(functn(ic,jstap1,kcp1) - functn(ic,jstap1,kcm1)  &
            - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
            + bcofz1*(functn(ic,jstap1,kcp2) - functn(ic,jstap1,kcm2)  &
            - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
        fdiffb = acofz1*(functn(ic,jstap2,kcp1) - functn(ic,jstap2,kcm1)  &
            - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
            + bcofz1*(functn(ic,jstap2,kcp2) - functn(ic,jstap2,kcm2)  &
            - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
        fdiffc = acofz1*(functn(ic,jstap3,kcp1) - functn(ic,jstap3,kcm1)  &
            - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
            + bcofz1*(functn(ic,jstap3,kcp2) - functn(ic,jstap3,kcm2)  &
            - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
        fdiffd = acofz1*(functn(ic,jstap4,kcp1) - functn(ic,jstap4,kcm1)  &
            - functn(ic,jstal,kcp1)  + functn(ic,jstal,kcm1))  &
            + bcofz1*(functn(ic,jstap4,kcp2) - functn(ic,jstap4,kcm2)  &
            - functn(ic,jstal,kcp2)  + functn(ic,jstal,kcm2))
        fderiv(ic,jstal,kc) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofz1*(functn(ic,jstal,kcp1)  - functn(ic,jstal,kcm1)  &
            - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
            + bcofz1*(functn(ic,jstal,kcp2)  - functn(ic,jstal,kcm2)  &
            - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
        fdiffb = acofz1*(functn(ic,jstap2,kcp1) - functn(ic,jstap2,kcm1)  &
            - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
            + bcofz1*(functn(ic,jstap2,kcp2) - functn(ic,jstap2,kcm2)  &
            - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
        fdiffc = acofz1*(functn(ic,jstap3,kcp1) - functn(ic,jstap3,kcm1)  &
            - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
            + bcofz1*(functn(ic,jstap3,kcp2) - functn(ic,jstap3,kcm2)  &
            - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
        fdiffd = acofz1*(functn(ic,jstap4,kcp1) - functn(ic,jstap4,kcm1)  &
            - functn(ic,jstap1,kcp1) + functn(ic,jstap1,kcm1))  &
            + bcofz1*(functn(ic,jstap4,kcp2) - functn(ic,jstap4,kcm2)  &
            - functn(ic,jstap1,kcp2) + functn(ic,jstap1,kcm2))
        fderiv(ic,jstap1,kc) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      ks = 0
      DO kc = kstom4,kstom2
        
        ks = ks+1
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
        js = 0
        DO jc = jstap2,jstap4
          
          js = js+1
          jcm2 = jc-2
          jcm1 = jc-1
          jcp1 = jc+1
          jcp2 = jc+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(ic,jcp1,kcp1) - functn(ic,jcp1,kcm1)  &
              - functn(ic,jcm1,kcp1) + functn(ic,jcm1,kcm1)
          fdiffb = functn(ic,jcp2,kcp2) - functn(ic,jcp2,kcm2)  &
              - functn(ic,jcm2,kcp2) + functn(ic,jcm2,kcm2)
          fderiv(ic,jc,kc) = acf3yz*fdiffa + bcf3yz*fdiffb
          fstora(js,ks) = fdiffa
          fstorb(js,ks) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      ks = 0
      DO kc = kstom4,kstom3
        
        ks = ks+1
        kcm3 = kc-3
        kcp3 = kc+3
        
        js = 1
        DO jc = jstap3,jstap4
          
          jsm1 = js
          js = js+1
          jcm3 = jc-3
          jcp3 = jc+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(ic,jcp3,kcp3) - functn(ic,jcp3,kcm3)  &
              - functn(ic,jcm3,kcp3) + functn(ic,jcm3,kcm3)
          fderiv(ic,jc,kc) = acf4yz*fstora(js,ks) + bcf4yz*fstorb(js,ks)  &
              + ccf4yz*fdiffc
          fstorc(jsm1,ks) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      ks = 1
      js = 3
      jsm1 = 2
      kc = kstom4
      jc = jstap4
      kcm4 = kc-4
      kcp4 = kc+4
      jcm4 = jc-4
      jcp4 = jc+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(ic,jcp4,kcp4) - functn(ic,jcp4,kcm4)  &
          - functn(ic,jcm4,kcp4) + functn(ic,jcm4,kcm4)
      fderiv(ic,jc,kc) = acf5yz*fstora(js,ks) + bcf5yz*fstorb(js,ks)  &
          + ccf5yz*fstorc(jsm1,ks) + dcf5yz*fdiffd
      
    END DO
    
  END IF
  
!       RH IN Y RH IN Z CORNER
!       ======================
  IF(nendyr == nbound)THEN
    
    DO ic = istal,istol
      
!           RH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(ic,jstom1,kstom1) - functn(ic,jstom1,kstol)  &
          - functn(ic,jstol,kstom1)  + functn(ic,jstol,kstol)
      fdiffb = functn(ic,jstom2,kstom2) - functn(ic,jstom2,kstol)  &
          - functn(ic,jstol,kstom2)  + functn(ic,jstol,kstol)
      fdiffc = functn(ic,jstom3,kstom3) - functn(ic,jstom3,kstol)  &
          - functn(ic,jstol,kstom3)  + functn(ic,jstol,kstol)
      fdiffd = functn(ic,jstom4,kstom4) - functn(ic,jstom4,kstol)  &
          - functn(ic,jstol,kstom4)  + functn(ic,jstol,kstol)
      fderiv(ic,jstol,kstol) = acc1yz*fdiffa + bcc1yz*fdiffb  &
          + ccc1yz*fdiffc + dcc1yz*fdiffd
      
!           RH-1 RH-1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(ic,jstol,kstol)   - functn(ic,jstol,kstom1)  &
          - functn(ic,jstom1,kstol)  + functn(ic,jstom1,kstom1)
      fdiffb = functn(ic,jstom2,kstom2) - functn(ic,jstom2,kstom1)  &
          - functn(ic,jstom1,kstom2) + functn(ic,jstom1,kstom1)
      fdiffc = functn(ic,jstom3,kstom3) - functn(ic,jstom3,kstom1)  &
          - functn(ic,jstom1,kstom3) + functn(ic,jstom1,kstom1)
      fdiffd = functn(ic,jstom4,kstom4) - functn(ic,jstom4,kstom1)  &
          - functn(ic,jstom1,kstom4) + functn(ic,jstom1,kstom1)
      fderiv(ic,jstom1,kstom1) = acc2yz*fdiffa + bcc2yz*fdiffb  &
          + ccc2yz*fdiffc + dcc2yz*fdiffd
      
!           RH RH-1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2yz*(functn(ic,jstom1,kstol)  - functn(ic,jstom1,kstom1)  &
          - functn(ic,jstol,kstol)   + functn(ic,jstol,kstom1))  &
          + bcf2yz*(functn(ic,jstom1,kstom2) - functn(ic,jstom1,kstom1)  &
          - functn(ic,jstol,kstom2)  + functn(ic,jstol,kstom1))  &
          + ccf2yz*(functn(ic,jstom1,kstom3) - functn(ic,jstom1,kstom1)  &
          - functn(ic,jstol,kstom3)  + functn(ic,jstol,kstom1))  &
          + dcf2yz*(functn(ic,jstom1,kstom4) - functn(ic,jstom1,kstom1)  &
          - functn(ic,jstol,kstom4)  + functn(ic,jstol,kstom1))
      fdiffb = acf2yz*(functn(ic,jstom2,kstol)  - functn(ic,jstom2,kstom1)  &
          - functn(ic,jstol,kstol)   + functn(ic,jstol,kstom1))  &
          + bcf2yz*(functn(ic,jstom2,kstom2) - functn(ic,jstom2,kstom1)  &
          - functn(ic,jstol,kstom2)  + functn(ic,jstol,kstom1))  &
          + ccf2yz*(functn(ic,jstom2,kstom3) - functn(ic,jstom2,kstom1)  &
          - functn(ic,jstol,kstom3)  + functn(ic,jstol,kstom1))  &
          + dcf2yz*(functn(ic,jstom2,kstom4) - functn(ic,jstom2,kstom1)  &
          - functn(ic,jstol,kstom4)  + functn(ic,jstol,kstom1))
      fdiffc = acf2yz*(functn(ic,jstom3,kstol)  - functn(ic,jstom3,kstom1)  &
          - functn(ic,jstol,kstol)   + functn(ic,jstol,kstom1))  &
          + bcf2yz*(functn(ic,jstom3,kstom2) - functn(ic,jstom3,kstom1)  &
          - functn(ic,jstol,kstom2)  + functn(ic,jstol,kstom1))  &
          + ccf2yz*(functn(ic,jstom3,kstom3) - functn(ic,jstom3,kstom1)  &
          - functn(ic,jstol,kstom3)  + functn(ic,jstol,kstom1))  &
          + dcf2yz*(functn(ic,jstom3,kstom4) - functn(ic,jstom3,kstom1)  &
          - functn(ic,jstol,kstom4)  + functn(ic,jstol,kstom1))
      fdiffd = acf2yz*(functn(ic,jstom4,kstol)  - functn(ic,jstom4,kstom1)  &
          - functn(ic,jstol,kstol)   + functn(ic,jstol,kstom1))  &
          + bcf2yz*(functn(ic,jstom4,kstom2) - functn(ic,jstom4,kstom1)  &
          - functn(ic,jstol,kstom2)  + functn(ic,jstol,kstom1))  &
          + ccf2yz*(functn(ic,jstom4,kstom3) - functn(ic,jstom4,kstom1)  &
          - functn(ic,jstol,kstom3)  + functn(ic,jstol,kstom1))  &
          + dcf2yz*(functn(ic,jstom4,kstom4) - functn(ic,jstom4,kstom1)  &
          - functn(ic,jstol,kstom4)  + functn(ic,jstol,kstom1))
      fderiv(ic,jstol,kstom1) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           RH+1 RH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2yz*(functn(ic,jstol,kstom1)  - functn(ic,jstom1,kstom1)  &
          - functn(ic,jstol,kstol)   + functn(ic,jstom1,kstol))  &
          + bcf2yz*(functn(ic,jstom2,kstom1) - functn(ic,jstom1,kstom1)  &
          - functn(ic,jstom2,kstol)  + functn(ic,jstom1,kstol))  &
          + ccf2yz*(functn(ic,jstom3,kstom1) - functn(ic,jstom1,kstom1)  &
          - functn(ic,jstom3,kstol)  + functn(ic,jstom1,kstol))  &
          + dcf2yz*(functn(ic,jstom4,kstom1) - functn(ic,jstom1,kstom1)  &
          - functn(ic,jstom4,kstol)  + functn(ic,jstom1,kstol))
      fdiffb = acf2yz*(functn(ic,jstol,kstom2)  - functn(ic,jstom1,kstom2)  &
          - functn(ic,jstol,kstol)   + functn(ic,jstom1,kstol))  &
          + bcf2yz*(functn(ic,jstom2,kstom2) - functn(ic,jstom1,kstom2)  &
          - functn(ic,jstom2,kstol)  + functn(ic,jstom1,kstol))  &
          + ccf2yz*(functn(ic,jstom3,kstom2) - functn(ic,jstom1,kstom2)  &
          - functn(ic,jstom3,kstol)  + functn(ic,jstom1,kstol))  &
          + dcf2yz*(functn(ic,jstom4,kstom2) - functn(ic,jstom1,kstom2)  &
          - functn(ic,jstom4,kstol)  + functn(ic,jstom1,kstol))
      fdiffc = acf2yz*(functn(ic,jstol,kstom3)  - functn(ic,jstom1,kstom3)  &
          - functn(ic,jstol,kstol)   + functn(ic,jstom1,kstol))  &
          + bcf2yz*(functn(ic,jstom2,kstom3) - functn(ic,jstom1,kstom3)  &
          - functn(ic,jstom2,kstol)  + functn(ic,jstom1,kstol))  &
          + ccf2yz*(functn(ic,jstom3,kstom3) - functn(ic,jstom1,kstom3)  &
          - functn(ic,jstom3,kstol)  + functn(ic,jstom1,kstol))  &
          + dcf2yz*(functn(ic,jstom4,kstom3) - functn(ic,jstom1,kstom3)  &
          - functn(ic,jstom4,kstol)  + functn(ic,jstom1,kstol))
      fdiffd = acf2yz*(functn(ic,jstol,kstom4)  - functn(ic,jstom1,kstom4)  &
          - functn(ic,jstol,kstol)   + functn(ic,jstom1,kstol))  &
          + bcf2yz*(functn(ic,jstom2,kstom4) - functn(ic,jstom1,kstom4)  &
          - functn(ic,jstom2,kstol)  + functn(ic,jstom1,kstol))  &
          + ccf2yz*(functn(ic,jstom3,kstom4) - functn(ic,jstom1,kstom4)  &
          - functn(ic,jstom3,kstol)  + functn(ic,jstom1,kstol))  &
          + dcf2yz*(functn(ic,jstom4,kstom4) - functn(ic,jstom1,kstom4)  &
          - functn(ic,jstom4,kstol)  + functn(ic,jstom1,kstol))
      fderiv(ic,jstom1,kstol) = acf1yz*fdiffa + bcf1yz*fdiffb  &
          + ccf1yz*fdiffc + dcf1yz*fdiffd
      
!           RH EDGE IN Z
      DO jc = jstom4,jstom2
        
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
            - functn(ic,jcp1,kstol)  + functn(ic,jcm1,kstol))  &
            + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
            - functn(ic,jcp2,kstol)  + functn(ic,jcm2,kstol))
        fdiffb = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
            - functn(ic,jcp1,kstom2) + functn(ic,jcm1,kstom2))  &
            + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
            - functn(ic,jcp2,kstom2) + functn(ic,jcm2,kstom2))
        fdiffc = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
            - functn(ic,jcp1,kstom3) + functn(ic,jcm1,kstom3))  &
            + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
            - functn(ic,jcp2,kstom3) + functn(ic,jcm2,kstom3))
        fdiffd = acofy1*(functn(ic,jcp1,kstom1) - functn(ic,jcm1,kstom1)  &
            - functn(ic,jcp1,kstom4) + functn(ic,jcm1,kstom4))  &
            + bcofy1*(functn(ic,jcp2,kstom1) - functn(ic,jcm2,kstom1)  &
            - functn(ic,jcp2,kstom4) + functn(ic,jcm2,kstom4))
        fderiv(ic,jc,kstom1) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
            - functn(ic,jcp1,kstom1) + functn(ic,jcm1,kstom1))  &
            + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
            - functn(ic,jcp2,kstom1) + functn(ic,jcm2,kstom1))
        fdiffb = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
            - functn(ic,jcp1,kstom2) + functn(ic,jcm1,kstom2))  &
            + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
            - functn(ic,jcp2,kstom2) + functn(ic,jcm2,kstom2))
        fdiffc = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
            - functn(ic,jcp1,kstom3) + functn(ic,jcm1,kstom3))  &
            + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
            - functn(ic,jcp2,kstom3) + functn(ic,jcm2,kstom3))
        fdiffd = acofy1*(functn(ic,jcp1,kstol)  - functn(ic,jcm1,kstol)  &
            - functn(ic,jcp1,kstom4) + functn(ic,jcm1,kstom4))  &
            + bcofy1*(functn(ic,jcp2,kstol)  - functn(ic,jcm2,kstol)  &
            - functn(ic,jcp2,kstom4) + functn(ic,jcm2,kstom4))
        fderiv(ic,jc,kstol) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd
        
      END DO
      
!           RH EDGE IN Y
      DO kc = kstom4,kstom2
        
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
            - functn(ic,jstom1,kcp1) + functn(ic,jstom1,kcm1))  &
            + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
            - functn(ic,jstom1,kcp2) + functn(ic,jstom1,kcm2))
        fdiffb = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
            - functn(ic,jstom2,kcp1) + functn(ic,jstom2,kcm1))  &
            + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
            - functn(ic,jstom2,kcp2) + functn(ic,jstom2,kcm2))
        fdiffc = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
            - functn(ic,jstom3,kcp1) + functn(ic,jstom3,kcm1))  &
            + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
            - functn(ic,jstom3,kcp2) + functn(ic,jstom3,kcm2))
        fdiffd = acofz1*(functn(ic,jstol,kcp1)  - functn(ic,jstol,kcm1)  &
            - functn(ic,jstom4,kcp1) + functn(ic,jstom4,kcm1))  &
            + bcofz1*(functn(ic,jstol,kcp2)  - functn(ic,jstol,kcm2)  &
            - functn(ic,jstom4,kcp2) + functn(ic,jstom4,kcm2))
        fderiv(ic,jstol,kc) = acf1yz*fdiffa + bcf1yz*fdiffb  &
            + ccf1yz*fdiffc + dcf1yz*fdiffd
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
            - functn(ic,jstol,kcp1)  + functn(ic,jstol,kcm1))  &
            + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
            - functn(ic,jstol,kcp2)  + functn(ic,jstol,kcm2))
        fdiffb = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
            - functn(ic,jstom2,kcp1) + functn(ic,jstom2,kcm1))  &
            + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
            - functn(ic,jstom2,kcp2) + functn(ic,jstom2,kcm2))
        fdiffc = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
            - functn(ic,jstom3,kcp1) + functn(ic,jstom3,kcm1))  &
            + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
            - functn(ic,jstom3,kcp2) + functn(ic,jstom3,kcm2))
        fdiffd = acofz1*(functn(ic,jstom1,kcp1) - functn(ic,jstom1,kcm1)  &
            - functn(ic,jstom4,kcp1) + functn(ic,jstom4,kcm1))  &
            + bcofz1*(functn(ic,jstom1,kcp2) - functn(ic,jstom1,kcm2)  &
            - functn(ic,jstom4,kcp2) + functn(ic,jstom4,kcm2))
        fderiv(ic,jstom1,kc) = acf2yz*fdiffa + bcf2yz*fdiffb  &
            + ccf2yz*fdiffc + dcf2yz*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      ks = 0
      DO kc = kstom4,kstom2
        
        ks = ks+1
        kcm2 = kc-2
        kcm1 = kc-1
        kcp1 = kc+1
        kcp2 = kc+2
        
        js = 0
        DO jc = jstom4,jstom2
          
          js = js+1
          jcm2 = jc-2
          jcm1 = jc-1
          jcp1 = jc+1
          jcp2 = jc+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(ic,jcp1,kcp1) - functn(ic,jcp1,kcm1)  &
              - functn(ic,jcm1,kcp1) + functn(ic,jcm1,kcm1)
          fdiffb = functn(ic,jcp2,kcp2) - functn(ic,jcp2,kcm2)  &
              - functn(ic,jcm2,kcp2) + functn(ic,jcm2,kcm2)
          fderiv(ic,jc,kc) = acf3yz*fdiffa + bcf3yz*fdiffb
          fstora(js,ks) = fdiffa
          fstorb(js,ks) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      ks = 0
      DO kc = kstom4,kstom3
        
        ks = ks+1
        kcm3 = kc-3
        kcp3 = kc+3
        
        js = 0
        DO jc = jstom4,jstom3
          
          js = js+1
          jcm3 = jc-3
          jcp3 = jc+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(ic,jcp3,kcp3) - functn(ic,jcp3,kcm3)  &
              - functn(ic,jcm3,kcp3) + functn(ic,jcm3,kcm3)
          fderiv(ic,jc,kc) = acf4yz*fstora(js,ks) + bcf4yz*fstorb(js,ks)  &
              + ccf4yz*fdiffc
          fstorc(js,ks) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      ks = 1
      js = 1
      kc = kstom4
      jc = jstom4
      kcm4 = kc-4
      kcp4 = kc+4
      jcm4 = jc-4
      jcp4 = jc+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(ic,jcp4,kcp4) - functn(ic,jcp4,kcm4)  &
          - functn(ic,jcm4,kcp4) + functn(ic,jcm4,kcm4)
      fderiv(ic,jc,kc) = acf5yz*fstora(js,ks) + bcf5yz*fstorb(js,ks)  &
          + ccf5yz*fstorc(js,ks) + dcf5yz*fdiffd
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     SCALING
!     =======
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fderiv(ic,jc,kc) = fderiv(ic,jc,kc)*ovdely*ovdelz
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE d2fdyz
