SUBROUTINE d2fdxy(functn,fderiv)

! Code converted using TO_F90 by Alan Miller
! Date: 2022-11-09  Time: 14:51:19

!     *************************************************************************

!     D2FDXY
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
!     EVALUATES SECOND XY-DERIVATIVE OF SPECIFIED FUNCTION
!     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
!     EXPLICIT 8TH,6TH,4TH,4TH,4TH COMPATIBLE ORDER END CONDITIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------

!     ARGUMENTS
!     =========

real(kind=dp), INTENT(IN)             :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
real(kind=dp), INTENT(OUT)            :: fderiv(nxsize,nysize,nzsize)




!     LOCAL DATA
!     ==========
real(kind=dp) :: fstora(3,3),fstorb(3,3),fstorc(2,2)
real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd,fdiffe
INTEGER :: ic,jc,kc
INTEGER :: is,js,ism1,jsm1
INTEGER :: istart,ifinis,jstart,jfinis
INTEGER :: icm5,icm4,icm3,icm2,icm1,iccc,icp1,icp2,icp3,icp4,icp5
INTEGER :: jcm5,jcm4,jcm3,jcm2,jcm1,jccc,jcp1,jcp2,jcp3,jcp4,jcp5


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

istart = istal
ifinis = istol
jstart = jstal
jfinis = jstol
IF(nendxl == nbound)istart = istap5
IF(nendxr == nbound)ifinis = istom5
IF(nendyl == nbound)jstart = jstap5
IF(nendyr == nbound)jfinis = jstom5

!     =========================================================================

!     INTERIOR SCHEME
!     ===============

!     TENTH ORDER EXPLICIT DIFFERENCES

DO kc = kstal,kstol
  
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
      
      fdiffa = functn(icp1,jcp1,kc) - functn(icp1,jcm1,kc)  &
             - functn(icm1,jcp1,kc) + functn(icm1,jcm1,kc)
      fdiffb = functn(icp2,jcp2,kc) - functn(icp2,jcm2,kc)  &
             - functn(icm2,jcp2,kc) + functn(icm2,jcm2,kc)
      fdiffc = functn(icp3,jcp3,kc) - functn(icp3,jcm3,kc)  &
             - functn(icm3,jcp3,kc) + functn(icm3,jcm3,kc)
      fdiffd = functn(icp4,jcp4,kc) - functn(icp4,jcm4,kc)  &
             - functn(icm4,jcp4,kc) + functn(icm4,jcm4,kc)
      fdiffe = functn(icp5,jcp5,kc) - functn(icp5,jcm5,kc)  &
             - functn(icm5,jcp5,kc) + functn(icm5,jcm5,kc)
      
      fderiv(ic,jc,kc) = acofxy*fdiffa + bcofxy*fdiffb  &
                       + ccofxy*fdiffc + dcofxy*fdiffd  &
                       + ecofxy*fdiffe
      
    END DO
    
  END DO
  
END DO

!     =========================================================================

!     LH END X-DIRECTION
!     ==================
IF(nendxl == nbound)THEN
  
!       TAKE SECOND XY-DERIVATIVE IN X-LEFT INNER HALO
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    
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
      
!           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofy1*(functn(istap1,jcp1,kc) - functn(istap1,jcm1,kc)  &
          - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
          + bcofy1*(functn(istap1,jcp2,kc) - functn(istap1,jcm2,kc)  &
          - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
      fdiffb = acofy1*(functn(istap2,jcp1,kc) - functn(istap2,jcm1,kc)  &
          - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
          + bcofy1*(functn(istap2,jcp2,kc) - functn(istap2,jcm2,kc)  &
          - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
      fdiffc = acofy1*(functn(istap3,jcp1,kc) - functn(istap3,jcm1,kc)  &
          - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
          + bcofy1*(functn(istap3,jcp2,kc) - functn(istap3,jcm2,kc)  &
          - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
      fdiffd = acofy1*(functn(istap4,jcp1,kc) - functn(istap4,jcm1,kc)  &
          - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
          + bcofy1*(functn(istap4,jcp2,kc) - functn(istap4,jcm2,kc)  &
          - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
      fderiv(istal,jc,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofy1*(functn(istal,jcp1,kc)  - functn(istal,jcm1,kc)  &
          - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
          + bcofy1*(functn(istal,jcp2,kc)  - functn(istal,jcm2,kc)  &
          - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
      fdiffb = acofy1*(functn(istap2,jcp1,kc) - functn(istap2,jcm1,kc)  &
          - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
          + bcofy1*(functn(istap2,jcp2,kc) - functn(istap2,jcm2,kc)  &
          - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
      fdiffc = acofy1*(functn(istap3,jcp1,kc) - functn(istap3,jcm1,kc)  &
          - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
          + bcofy1*(functn(istap3,jcp2,kc) - functn(istap3,jcm2,kc)  &
          - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
      fdiffd = acofy1*(functn(istap4,jcp1,kc) - functn(istap4,jcm1,kc)  &
          - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
          + bcofy1*(functn(istap4,jcp2,kc) - functn(istap4,jcm2,kc)  &
          - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
      fderiv(istap1,jc,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
          + ccf2xy*fdiffc + dcf2xy*fdiffd
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdiffa = functn(istap3,jcp1,kc) - functn(istap3,jcm1,kc)  &
          - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc)
      fdiffb = functn(istap4,jcp2,kc) - functn(istap4,jcm2,kc)  &
          - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc)
      fderiv(istap2,jc,kc) = acf3xy*fdiffa + bcf3xy*fdiffb
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdiffa = functn(istap4,jcp1,kc) - functn(istap4,jcm1,kc)  &
          - functn(istap2,jcp1,kc) + functn(istap2,jcm1,kc)
      fdiffb = functn(istap5,jcp2,kc) - functn(istap5,jcm2,kc)  &
          - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc)
      fdiffc = functn(istap6,jcp3,kc) - functn(istap6,jcm3,kc)  &
          - functn(istal,jcp3,kc)  + functn(istal,jcm3,kc)
      fderiv(istap3,jc,kc) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdiffa = functn(istap5,jcp1,kc) - functn(istap5,jcm1,kc)  &
          - functn(istap3,jcp1,kc) + functn(istap3,jcm1,kc)
      fdiffb = functn(istap6,jcp2,kc) - functn(istap6,jcm2,kc)  &
          - functn(istap2,jcp2,kc) + functn(istap2,jcm2,kc)
      fdiffc = functn(istap7,jcp3,kc) - functn(istap7,jcm3,kc)  &
          - functn(istap1,jcp3,kc) + functn(istap1,jcm3,kc)
      fdiffd = functn(istap8,jcp4,kc) - functn(istap8,jcm4,kc)  &
          - functn(istal,jcp4,kc)  + functn(istal,jcm4,kc)
      fderiv(istap4,jc,kc) = acf5xy*fdiffa + bcf5xy*fdiffb  &
          + ccf5xy*fdiffc + dcf5xy*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     RH END X-DIRECTION
!     ==================
IF(nendxr == nbound)THEN
  
!       TAKE SECOND XY-DERIVATIVE IN X-RIGHT INNER HALO
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    
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
      
!           RH POINT MINUS 4: 8TH ORDER CENTRED
      fdiffa = functn(istom3,jcp1,kc) - functn(istom3,jcm1,kc)  &
          - functn(istom5,jcp1,kc) + functn(istom5,jcm1,kc)
      fdiffb = functn(istom2,jcp2,kc) - functn(istom2,jcm2,kc)  &
          - functn(istom6,jcp2,kc) + functn(istom6,jcm2,kc)
      fdiffc = functn(istom1,jcp3,kc) - functn(istom1,jcm3,kc)  &
          - functn(istom7,jcp3,kc) + functn(istom7,jcm3,kc)
      fdiffd = functn(istol,jcp4,kc)  - functn(istol,jcm4,kc)  &
          - functn(istom8,jcp4,kc) + functn(istom8,jcm4,kc)
      fderiv(istom4,jc,kc) = acf5xy*fdiffa + bcf5xy*fdiffb  &
          + ccf5xy*fdiffc + dcf5xy*fdiffd
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdiffa = functn(istom2,jcp1,kc) - functn(istom2,jcm1,kc)  &
          - functn(istom4,jcp1,kc) + functn(istom4,jcm1,kc)
      fdiffb = functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
          - functn(istom5,jcp2,kc) + functn(istom5,jcm2,kc)
      fdiffc = functn(istol,jcp3,kc)  - functn(istol,jcm3,kc)  &
          - functn(istom6,jcp3,kc) + functn(istom6,jcm3,kc)
      fderiv(istom3,jc,kc) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdiffa = functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
          - functn(istom3,jcp1,kc) + functn(istom3,jcm1,kc)
      fdiffb = functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
          - functn(istom4,jcp2,kc) + functn(istom4,jcm2,kc)
      fderiv(istom2,jc,kc) = acf3xy*fdiffa + bcf3xy*fdiffb
      
!           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
          - functn(istol,jcp1,kc)  + functn(istol,jcm1,kc))  &
          + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
          - functn(istol,jcp2,kc)  + functn(istol,jcm2,kc))
      fdiffb = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
          - functn(istom2,jcp1,kc) + functn(istom2,jcm1,kc))  &
          + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
          - functn(istom2,jcp2,kc) + functn(istom2,jcm2,kc))
      fdiffc = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
          - functn(istom3,jcp1,kc) + functn(istom3,jcm1,kc))  &
          + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
          - functn(istom3,jcp2,kc) + functn(istom3,jcm2,kc))
      fdiffd = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
          - functn(istom4,jcp1,kc) + functn(istom4,jcm1,kc))  &
          + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
          - functn(istom4,jcp2,kc) + functn(istom4,jcm2,kc))
      fderiv(istom1,jc,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
          + ccf2xy*fdiffc + dcf2xy*fdiffd
      
!           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
          - functn(istom1,jcp1,kc) + functn(istom1,jcm1,kc))  &
          + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
          - functn(istom1,jcp2,kc) + functn(istom1,jcm2,kc))
      fdiffb = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
          - functn(istom2,jcp1,kc) + functn(istom2,jcm1,kc))  &
          + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
          - functn(istom2,jcp2,kc) + functn(istom2,jcm2,kc))
      fdiffc = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
          - functn(istom3,jcp1,kc) + functn(istom3,jcm1,kc))  &
          + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
          - functn(istom3,jcp2,kc) + functn(istom3,jcm2,kc))
      fdiffd = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
          - functn(istom4,jcp1,kc) + functn(istom4,jcm1,kc))  &
          + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
          - functn(istom4,jcp2,kc) + functn(istom4,jcm2,kc))
      fderiv(istol,jc,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     LH END Y-DIRECTION
!     ==================
IF(nendyl == nbound)THEN
  
!       TAKE SECOND XY-DERIVATIVE IN Y-LEFT INNER HALO
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    
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
      fdiffa = acofx1*(functn(icp1,jstap1,kc) - functn(icm1,jstap1,kc)  &
          - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
          + bcofx1*(functn(icp2,jstap1,kc) - functn(icm2,jstap1,kc)  &
          - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
      fdiffb = acofx1*(functn(icp1,jstap2,kc) - functn(icm1,jstap2,kc)  &
          - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
          + bcofx1*(functn(icp2,jstap2,kc) - functn(icm2,jstap2,kc)  &
          - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
      fdiffc = acofx1*(functn(icp1,jstap3,kc) - functn(icm1,jstap3,kc)  &
          - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
          + bcofx1*(functn(icp2,jstap3,kc) - functn(icm2,jstap3,kc)  &
          - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
      fdiffd = acofx1*(functn(icp1,jstap4,kc) - functn(icm1,jstap4,kc)  &
          - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
          + bcofx1*(functn(icp2,jstap4,kc) - functn(icm2,jstap4,kc)  &
          - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
      fderiv(ic,jstal,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofx1*(functn(icp1,jstal,kc)  - functn(icm1,jstal,kc)  &
          - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
          + bcofx1*(functn(icp2,jstal,kc)  - functn(icm2,jstal,kc)  &
          - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
      fdiffb = acofx1*(functn(icp1,jstap2,kc) - functn(icm1,jstap2,kc)  &
          - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
          + bcofx1*(functn(icp2,jstap2,kc) - functn(icm2,jstap2,kc)  &
          - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
      fdiffc = acofx1*(functn(icp1,jstap3,kc) - functn(icm1,jstap3,kc)  &
          - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
          + bcofx1*(functn(icp2,jstap3,kc) - functn(icm2,jstap3,kc)  &
          - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
      fdiffd = acofx1*(functn(icp1,jstap4,kc) - functn(icm1,jstap4,kc)  &
          - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
          + bcofx1*(functn(icp2,jstap4,kc) - functn(icm2,jstap4,kc)  &
          - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
      fderiv(ic,jstap1,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
          + ccf2xy*fdiffc + dcf2xy*fdiffd
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdiffa = functn(icp1,jstap3,kc) - functn(icm1,jstap3,kc)  &
          - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc)
      fdiffb = functn(icp2,jstap4,kc) - functn(icm2,jstap4,kc)  &
          - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc)
      fderiv(ic,jstap2,kc) = acf3xy*fdiffa + bcf3xy*fdiffb
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdiffa = functn(icp1,jstap4,kc) - functn(icm1,jstap4,kc)  &
          - functn(icp1,jstap2,kc) + functn(icm1,jstap2,kc)
      fdiffb = functn(icp2,jstap5,kc) - functn(icm2,jstap5,kc)  &
          - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc)
      fdiffc = functn(icp3,jstap6,kc) - functn(icm3,jstap6,kc)  &
          - functn(icp3,jstal,kc)  + functn(icm3,jstal,kc)
      fderiv(ic,jstap3,kc) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdiffa = functn(icp1,jstap5,kc) - functn(icm1,jstap5,kc)  &
          - functn(icp1,jstap3,kc) + functn(icm1,jstap3,kc)
      fdiffb = functn(icp2,jstap6,kc) - functn(icm2,jstap6,kc)  &
          - functn(icp2,jstap2,kc) + functn(icm2,jstap2,kc)
      fdiffc = functn(icp3,jstap7,kc) - functn(icm3,jstap7,kc)  &
          - functn(icp3,jstap1,kc) + functn(icm3,jstap1,kc)
      fdiffd = functn(icp4,jstap8,kc) - functn(icm4,jstap8,kc)  &
          - functn(icp4,jstal,kc)  + functn(icm4,jstal,kc)
      fderiv(ic,jstap4,kc) = acf5xy*fdiffa + bcf5xy*fdiffb  &
          + ccf5xy*fdiffc + dcf5xy*fdiffd
      
    END DO
  END DO
  
!       LH IN X LH IN Y CORNER
!       ======================
  IF(nendxl == nbound)THEN
    
    DO kc = kstal,kstol
      
!           LH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(istap1,jstap1,kc) - functn(istap1,jstal,kc)  &
          - functn(istal,jstap1,kc)  + functn(istal,jstal,kc)
      fdiffb = functn(istap2,jstap2,kc) - functn(istap2,jstal,kc)  &
          - functn(istal,jstap2,kc)  + functn(istal,jstal,kc)
      fdiffc = functn(istap3,jstap3,kc) - functn(istap3,jstal,kc)  &
          - functn(istal,jstap3,kc)  + functn(istal,jstal,kc)
      fdiffd = functn(istap4,jstap4,kc) - functn(istap4,jstal,kc)  &
          - functn(istal,jstap4,kc)  + functn(istal,jstal,kc)
      fderiv(istal,jstal,kc) = acc1xy*fdiffa + bcc1xy*fdiffb  &
          + ccc1xy*fdiffc + dcc1xy*fdiffd
      
!           LH+1 LH+1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(istal,jstal,kc)   - functn(istal,jstap1,kc)  &
          - functn(istap1,jstal,kc)  + functn(istap1,jstap1,kc)
      fdiffb = functn(istap2,jstap2,kc) - functn(istap2,jstap1,kc)  &
          - functn(istap1,jstap2,kc) + functn(istap1,jstap1,kc)
      fdiffc = functn(istap3,jstap3,kc) - functn(istap3,jstap1,kc)  &
          - functn(istap1,jstap3,kc) + functn(istap1,jstap1,kc)
      fdiffd = functn(istap4,jstap4,kc) - functn(istap4,jstap1,kc)  &
          - functn(istap1,jstap4,kc) + functn(istap1,jstap1,kc)
      fderiv(istap1,jstap1,kc) = acc2xy*fdiffa + bcc2xy*fdiffb  &
          + ccc2xy*fdiffc + dcc2xy*fdiffd
      
!           LH LH+1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xy*(functn(istap1,jstal,kc)  - functn(istap1,jstap1,kc)  &
          - functn(istal,jstal,kc)   + functn(istal,jstap1,kc))  &
          + bcf2xy*(functn(istap1,jstap2,kc) - functn(istap1,jstap1,kc)  &
          - functn(istal,jstap2,kc)  + functn(istal,jstap1,kc))  &
          + ccf2xy*(functn(istap1,jstap3,kc) - functn(istap1,jstap1,kc)  &
          - functn(istal,jstap3,kc)  + functn(istal,jstap1,kc))  &
          + dcf2xy*(functn(istap1,jstap4,kc) - functn(istap1,jstap1,kc)  &
          - functn(istal,jstap4,kc)  + functn(istal,jstap1,kc))
      fdiffb = acf2xy*(functn(istap2,jstal,kc)  - functn(istap2,jstap1,kc)  &
          - functn(istal,jstal,kc)   + functn(istal,jstap1,kc))  &
          + bcf2xy*(functn(istap2,jstap2,kc) - functn(istap2,jstap1,kc)  &
          - functn(istal,jstap2,kc)  + functn(istal,jstap1,kc))  &
          + ccf2xy*(functn(istap2,jstap3,kc) - functn(istap2,jstap1,kc)  &
          - functn(istal,jstap3,kc)  + functn(istal,jstap1,kc))  &
          + dcf2xy*(functn(istap2,jstap4,kc) - functn(istap2,jstap1,kc)  &
          - functn(istal,jstap4,kc)  + functn(istal,jstap1,kc))
      fdiffc = acf2xy*(functn(istap3,jstal,kc)  - functn(istap3,jstap1,kc)  &
          - functn(istal,jstal,kc)   + functn(istal,jstap1,kc))  &
          + bcf2xy*(functn(istap3,jstap2,kc) - functn(istap3,jstap1,kc)  &
          - functn(istal,jstap2,kc)  + functn(istal,jstap1,kc))  &
          + ccf2xy*(functn(istap3,jstap3,kc) - functn(istap3,jstap1,kc)  &
          - functn(istal,jstap3,kc)  + functn(istal,jstap1,kc))  &
          + dcf2xy*(functn(istap3,jstap4,kc) - functn(istap3,jstap1,kc)  &
          - functn(istal,jstap4,kc)  + functn(istal,jstap1,kc))
      fdiffd = acf2xy*(functn(istap4,jstal,kc)  - functn(istap4,jstap1,kc)  &
          - functn(istal,jstal,kc)   + functn(istal,jstap1,kc))  &
          + bcf2xy*(functn(istap4,jstap2,kc) - functn(istap4,jstap1,kc)  &
          - functn(istal,jstap2,kc)  + functn(istal,jstap1,kc))  &
          + ccf2xy*(functn(istap4,jstap3,kc) - functn(istap4,jstap1,kc)  &
          - functn(istal,jstap3,kc)  + functn(istal,jstap1,kc))  &
          + dcf2xy*(functn(istap4,jstap4,kc) - functn(istap4,jstap1,kc)  &
          - functn(istal,jstap4,kc)  + functn(istal,jstap1,kc))
      fderiv(istal,jstap1,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           LH+1 LH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xy*(functn(istal,jstap1,kc)  - functn(istap1,jstap1,kc)  &
          - functn(istal,jstal,kc)   + functn(istap1,jstal,kc))  &
          + bcf2xy*(functn(istap2,jstap1,kc) - functn(istap1,jstap1,kc)  &
          - functn(istap2,jstal,kc)  + functn(istap1,jstal,kc))  &
          + ccf2xy*(functn(istap3,jstap1,kc) - functn(istap1,jstap1,kc)  &
          - functn(istap3,jstal,kc)  + functn(istap1,jstal,kc))  &
          + dcf2xy*(functn(istap4,jstap1,kc) - functn(istap1,jstap1,kc)  &
          - functn(istap4,jstal,kc)  + functn(istap1,jstal,kc))
      fdiffb = acf2xy*(functn(istal,jstap2,kc)  - functn(istap1,jstap2,kc)  &
          - functn(istal,jstal,kc)   + functn(istap1,jstal,kc))  &
          + bcf2xy*(functn(istap2,jstap2,kc) - functn(istap1,jstap2,kc)  &
          - functn(istap2,jstal,kc)  + functn(istap1,jstal,kc))  &
          + ccf2xy*(functn(istap3,jstap2,kc) - functn(istap1,jstap2,kc)  &
          - functn(istap3,jstal,kc)  + functn(istap1,jstal,kc))  &
          + dcf2xy*(functn(istap4,jstap2,kc) - functn(istap1,jstap2,kc)  &
          - functn(istap4,jstal,kc)  + functn(istap1,jstal,kc))
      fdiffc = acf2xy*(functn(istal,jstap3,kc)  - functn(istap1,jstap3,kc)  &
          - functn(istal,jstal,kc)   + functn(istap1,jstal,kc))  &
          + bcf2xy*(functn(istap2,jstap3,kc) - functn(istap1,jstap3,kc)  &
          - functn(istap2,jstal,kc)  + functn(istap1,jstal,kc))  &
          + ccf2xy*(functn(istap3,jstap3,kc) - functn(istap1,jstap3,kc)  &
          - functn(istap3,jstal,kc)  + functn(istap1,jstal,kc))  &
          + dcf2xy*(functn(istap4,jstap3,kc) - functn(istap1,jstap3,kc)  &
          - functn(istap4,jstal,kc)  + functn(istap1,jstal,kc))
      fdiffd = acf2xy*(functn(istal,jstap4,kc)  - functn(istap1,jstap4,kc)  &
          - functn(istal,jstal,kc)   + functn(istap1,jstal,kc))  &
          + bcf2xy*(functn(istap2,jstap4,kc) - functn(istap1,jstap4,kc)  &
          - functn(istap2,jstal,kc)  + functn(istap1,jstal,kc))  &
          + ccf2xy*(functn(istap3,jstap4,kc) - functn(istap1,jstap4,kc)  &
          - functn(istap3,jstal,kc)  + functn(istap1,jstal,kc))  &
          + dcf2xy*(functn(istap4,jstap4,kc) - functn(istap1,jstap4,kc)  &
          - functn(istap4,jstal,kc)  + functn(istap1,jstal,kc))
      fderiv(istap1,jstal,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           LH EDGE IN Y
      DO ic = istap2,istap4
        
        icm2 = ic-2
        icm1 = ic-1
        icp1 = ic+1
        icp2 = ic+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofx1*(functn(icp1,jstap1,kc) - functn(icm1,jstap1,kc)  &
            - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
            + bcofx1*(functn(icp2,jstap1,kc) - functn(icm2,jstap1,kc)  &
            - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
        fdiffb = acofx1*(functn(icp1,jstap2,kc) - functn(icm1,jstap2,kc)  &
            - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
            + bcofx1*(functn(icp2,jstap2,kc) - functn(icm2,jstap2,kc)  &
            - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
        fdiffc = acofx1*(functn(icp1,jstap3,kc) - functn(icm1,jstap3,kc)  &
            - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
            + bcofx1*(functn(icp2,jstap3,kc) - functn(icm2,jstap3,kc)  &
            - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
        fdiffd = acofx1*(functn(icp1,jstap4,kc) - functn(icm1,jstap4,kc)  &
            - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
            + bcofx1*(functn(icp2,jstap4,kc) - functn(icm2,jstap4,kc)  &
            - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
        fderiv(ic,jstal,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofx1*(functn(icp1,jstal,kc)  - functn(icm1,jstal,kc)  &
            - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
            + bcofx1*(functn(icp2,jstal,kc)  - functn(icm2,jstal,kc)  &
            - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
        fdiffb = acofx1*(functn(icp1,jstap2,kc) - functn(icm1,jstap2,kc)  &
            - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
            + bcofx1*(functn(icp2,jstap2,kc) - functn(icm2,jstap2,kc)  &
            - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
        fdiffc = acofx1*(functn(icp1,jstap3,kc) - functn(icm1,jstap3,kc)  &
            - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
            + bcofx1*(functn(icp2,jstap3,kc) - functn(icm2,jstap3,kc)  &
            - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
        fdiffd = acofx1*(functn(icp1,jstap4,kc) - functn(icm1,jstap4,kc)  &
            - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
            + bcofx1*(functn(icp2,jstap4,kc) - functn(icm2,jstap4,kc)  &
            - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
        fderiv(ic,jstap1,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd
        
      END DO
      
!           LH EDGE IN X
      DO jc = jstap2,jstap4
        
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofy1*(functn(istap1,jcp1,kc) - functn(istap1,jcm1,kc)  &
            - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
            + bcofy1*(functn(istap1,jcp2,kc) - functn(istap1,jcm2,kc)  &
            - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
        fdiffb = acofy1*(functn(istap2,jcp1,kc) - functn(istap2,jcm1,kc)  &
            - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
            + bcofy1*(functn(istap2,jcp2,kc) - functn(istap2,jcm2,kc)  &
            - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
        fdiffc = acofy1*(functn(istap3,jcp1,kc) - functn(istap3,jcm1,kc)  &
            - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
            + bcofy1*(functn(istap3,jcp2,kc) - functn(istap3,jcm2,kc)  &
            - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
        fdiffd = acofy1*(functn(istap4,jcp1,kc) - functn(istap4,jcm1,kc)  &
            - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
            + bcofy1*(functn(istap4,jcp2,kc) - functn(istap4,jcm2,kc)  &
            - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
        fderiv(istal,jc,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofy1*(functn(istal,jcp1,kc)  - functn(istal,jcm1,kc)  &
            - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
            + bcofy1*(functn(istal,jcp2,kc)  - functn(istal,jcm2,kc)  &
            - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
        fdiffb = acofy1*(functn(istap2,jcp1,kc) - functn(istap2,jcm1,kc)  &
            - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
            + bcofy1*(functn(istap2,jcp2,kc) - functn(istap2,jcm2,kc)  &
            - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
        fdiffc = acofy1*(functn(istap3,jcp1,kc) - functn(istap3,jcm1,kc)  &
            - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
            + bcofy1*(functn(istap3,jcp2,kc) - functn(istap3,jcm2,kc)  &
            - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
        fdiffd = acofy1*(functn(istap4,jcp1,kc) - functn(istap4,jcm1,kc)  &
            - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
            + bcofy1*(functn(istap4,jcp2,kc) - functn(istap4,jcm2,kc)  &
            - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
        fderiv(istap1,jc,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      js = 0
      DO jc = jstap2,jstap4
        
        js = js+1
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
        is = 0
        DO ic = istap2,istap4
          
          is = is+1
          icm2 = ic-2
          icm1 = ic-1
          icp1 = ic+1
          icp2 = ic+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(icp1,jcp1,kc) - functn(icp1,jcm1,kc)  &
              - functn(icm1,jcp1,kc) + functn(icm1,jcm1,kc)
          fdiffb = functn(icp2,jcp2,kc) - functn(icp2,jcm2,kc)  &
              - functn(icm2,jcp2,kc) + functn(icm2,jcm2,kc)
          fderiv(ic,jc,kc) = acf3xy*fdiffa + bcf3xy*fdiffb
          fstora(is,js) = fdiffa
          fstorb(is,js) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      js = 1
      DO jc = jstap3,jstap4
        
        jsm1 = js
        js = js+1
        jcm3 = jc-3
        jcp3 = jc+3
        
        is = 1
        DO ic = istap3,istap4
          
          ism1 = is
          is = is+1
          icm3 = ic-3
          icp3 = ic+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(icp3,jcp3,kc) - functn(icp3,jcm3,kc)  &
              - functn(icm3,jcp3,kc) + functn(icm3,jcm3,kc)
          fderiv(ic,jc,kc) = acf4xy*fstora(is,js) + bcf4xy*fstorb(is,js)  &
              + ccf4xy*fdiffc
          fstorc(ism1,jsm1) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      js = 3
      is = 3
      jsm1 = 2
      ism1 = 2
      jc = jstap4
      ic = istap4
      jcm4 = jc-4
      jcp4 = jc+4
      icm4 = ic-4
      icp4 = ic+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(icp4,jcp4,kc) - functn(icp4,jcm4,kc)  &
          - functn(icm4,jcp4,kc) + functn(icm4,jcm4,kc)
      fderiv(ic,jc,kc) = acf5xy*fstora(is,js) + bcf5xy*fstorb(is,js)  &
          + ccf5xy*fstorc(ism1,jsm1) + dcf5xy*fdiffd
      
    END DO
    
  END IF
  
!       RH IN X LH IN Y CORNER
!       ======================
  IF(nendxr == nbound)THEN
    
    DO kc = kstal,kstol
      
!           RH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(istol,jstap1,kc)  - functn(istol,jstal,kc)  &
          - functn(istom1,jstap1,kc) + functn(istom1,jstal,kc)
      fdiffb = functn(istol,jstap2,kc)  - functn(istol,jstal,kc)  &
          - functn(istom2,jstap2,kc) + functn(istom2,jstal,kc)
      fdiffc = functn(istol,jstap3,kc)  - functn(istol,jstal,kc)  &
          - functn(istom3,jstap3,kc) + functn(istom3,jstal,kc)
      fdiffd = functn(istol,jstap4,kc)  - functn(istol,jstal,kc)  &
          - functn(istom4,jstap4,kc) + functn(istom4,jstal,kc)
      fderiv(istol,jstal,kc) = acc1xy*fdiffa + bcc1xy*fdiffb  &
          + ccc1xy*fdiffc + dcc1xy*fdiffd
      
!           RH-1 LH+1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(istom1,jstal,kc)  - functn(istom1,jstap1,kc)  &
          - functn(istol,jstal,kc)   + functn(istol,jstap1,kc)
      fdiffb = functn(istom1,jstap2,kc) - functn(istom1,jstap1,kc)  &
          - functn(istom2,jstap2,kc) + functn(istom2,jstap1,kc)
      fdiffc = functn(istom1,jstap3,kc) - functn(istom1,jstap1,kc)  &
          - functn(istom3,jstap3,kc) + functn(istom3,jstap1,kc)
      fdiffd = functn(istom1,jstap4,kc) - functn(istom1,jstap1,kc)  &
          - functn(istom4,jstap4,kc) + functn(istom4,jstap1,kc)
      fderiv(istom1,jstap1,kc) = acc2xy*fdiffa + bcc2xy*fdiffb  &
          + ccc2xy*fdiffc + dcc2xy*fdiffd
      
!           RH LH+1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xy*(functn(istol,jstal,kc)   - functn(istol,jstap1,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom1,jstap1,kc))  &
          + bcf2xy*(functn(istol,jstap2,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom1,jstap2,kc) + functn(istom1,jstap1,kc))  &
          + ccf2xy*(functn(istol,jstap3,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom1,jstap3,kc) + functn(istom1,jstap1,kc))  &
          + dcf2xy*(functn(istol,jstap4,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom1,jstap4,kc) + functn(istom1,jstap1,kc))
      fdiffb = acf2xy*(functn(istol,jstal,kc)   - functn(istol,jstap1,kc)  &
          - functn(istom2,jstal,kc)  + functn(istom2,jstap1,kc))  &
          + bcf2xy*(functn(istol,jstap2,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom2,jstap2,kc) + functn(istom2,jstap1,kc))  &
          + ccf2xy*(functn(istol,jstap3,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom2,jstap3,kc) + functn(istom2,jstap1,kc))  &
          + dcf2xy*(functn(istol,jstap4,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom2,jstap4,kc) + functn(istom2,jstap1,kc))
      fdiffc = acf2xy*(functn(istol,jstal,kc)   - functn(istol,jstap1,kc)  &
          - functn(istom3,jstal,kc)  + functn(istom3,jstap1,kc))  &
          + bcf2xy*(functn(istol,jstap2,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom3,jstap2,kc) + functn(istom3,jstap1,kc))  &
          + ccf2xy*(functn(istol,jstap3,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom3,jstap3,kc) + functn(istom3,jstap1,kc))  &
          + dcf2xy*(functn(istol,jstap4,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom3,jstap4,kc) + functn(istom3,jstap1,kc))
      fdiffd = acf2xy*(functn(istol,jstal,kc)   - functn(istol,jstap1,kc)  &
          - functn(istom4,jstal,kc)  + functn(istom4,jstap1,kc))  &
          + bcf2xy*(functn(istol,jstap2,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom4,jstap2,kc) + functn(istom4,jstap1,kc))  &
          + ccf2xy*(functn(istol,jstap3,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom4,jstap3,kc) + functn(istom4,jstap1,kc))  &
          + dcf2xy*(functn(istol,jstap4,kc)  - functn(istol,jstap1,kc)  &
          - functn(istom4,jstap4,kc) + functn(istom4,jstap1,kc))
      fderiv(istol,jstap1,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           RH-1 LH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xy*(functn(istom1,jstap1,kc) - functn(istol,jstap1,kc)  &
          - functn(istom1,jstal,kc)  + functn(istol,jstal,kc))  &
          + bcf2xy*(functn(istom1,jstap1,kc) - functn(istom2,jstap1,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom2,jstal,kc))  &
          + ccf2xy*(functn(istom1,jstap1,kc) - functn(istom3,jstap1,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom3,jstal,kc))  &
          + dcf2xy*(functn(istom1,jstap1,kc) - functn(istom4,jstap1,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom4,jstal,kc))
      fdiffb = acf2xy*(functn(istom1,jstap2,kc) - functn(istol,jstap2,kc)  &
          - functn(istom1,jstal,kc)  + functn(istol,jstal,kc))  &
          + bcf2xy*(functn(istom1,jstap2,kc) - functn(istom2,jstap2,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom2,jstal,kc))  &
          + ccf2xy*(functn(istom1,jstap2,kc) - functn(istom3,jstap2,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom3,jstal,kc))  &
          + dcf2xy*(functn(istom1,jstap2,kc) - functn(istom4,jstap2,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom4,jstal,kc))
      fdiffc = acf2xy*(functn(istom1,jstap3,kc) - functn(istol,jstap3,kc)  &
          - functn(istom1,jstal,kc)  + functn(istol,jstal,kc))  &
          + bcf2xy*(functn(istom1,jstap3,kc) - functn(istom2,jstap3,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom2,jstal,kc))  &
          + ccf2xy*(functn(istom1,jstap3,kc) - functn(istom3,jstap3,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom3,jstal,kc))  &
          + dcf2xy*(functn(istom1,jstap3,kc) - functn(istom4,jstap3,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom4,jstal,kc))
      fdiffd = acf2xy*(functn(istom1,jstap4,kc) - functn(istol,jstap4,kc)  &
          - functn(istom1,jstal,kc)  + functn(istol,jstal,kc))  &
          + bcf2xy*(functn(istom1,jstap4,kc) - functn(istom2,jstap4,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom2,jstal,kc))  &
          + ccf2xy*(functn(istom1,jstap4,kc) - functn(istom3,jstap4,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom3,jstal,kc))  &
          + dcf2xy*(functn(istom1,jstap4,kc) - functn(istom4,jstap4,kc)  &
          - functn(istom1,jstal,kc)  + functn(istom4,jstal,kc))
      fderiv(istom1,jstal,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           LH EDGE IN Y
      DO ic = istom4,istom2
        
        icm2 = ic-2
        icm1 = ic-1
        icp1 = ic+1
        icp2 = ic+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofx1*(functn(icp1,jstap1,kc) - functn(icm1,jstap1,kc)  &
            - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
            + bcofx1*(functn(icp2,jstap1,kc) - functn(icm2,jstap1,kc)  &
            - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
        fdiffb = acofx1*(functn(icp1,jstap2,kc) - functn(icm1,jstap2,kc)  &
            - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
            + bcofx1*(functn(icp2,jstap2,kc) - functn(icm2,jstap2,kc)  &
            - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
        fdiffc = acofx1*(functn(icp1,jstap3,kc) - functn(icm1,jstap3,kc)  &
            - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
            + bcofx1*(functn(icp2,jstap3,kc) - functn(icm2,jstap3,kc)  &
            - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
        fdiffd = acofx1*(functn(icp1,jstap4,kc) - functn(icm1,jstap4,kc)  &
            - functn(icp1,jstal,kc)  + functn(icm1,jstal,kc))  &
            + bcofx1*(functn(icp2,jstap4,kc) - functn(icm2,jstap4,kc)  &
            - functn(icp2,jstal,kc)  + functn(icm2,jstal,kc))
        fderiv(ic,jstal,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofx1*(functn(icp1,jstal,kc)  - functn(icm1,jstal,kc)  &
            - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
            + bcofx1*(functn(icp2,jstal,kc)  - functn(icm2,jstal,kc)  &
            - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
        fdiffb = acofx1*(functn(icp1,jstap2,kc) - functn(icm1,jstap2,kc)  &
            - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
            + bcofx1*(functn(icp2,jstap2,kc) - functn(icm2,jstap2,kc)  &
            - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
        fdiffc = acofx1*(functn(icp1,jstap3,kc) - functn(icm1,jstap3,kc)  &
            - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
            + bcofx1*(functn(icp2,jstap3,kc) - functn(icm2,jstap3,kc)  &
            - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
        fdiffd = acofx1*(functn(icp1,jstap4,kc) - functn(icm1,jstap4,kc)  &
            - functn(icp1,jstap1,kc) + functn(icm1,jstap1,kc))  &
            + bcofx1*(functn(icp2,jstap4,kc) - functn(icm2,jstap4,kc)  &
            - functn(icp2,jstap1,kc) + functn(icm2,jstap1,kc))
        fderiv(ic,jstap1,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd
        
      END DO
      
!           RH EDGE IN X
      DO jc = jstap2,jstap4
        
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
            - functn(istom1,jcp1,kc) + functn(istom1,jcm1,kc))  &
            + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
            - functn(istom1,jcp2,kc) + functn(istom1,jcm2,kc))
        fdiffb = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
            - functn(istom2,jcp1,kc) + functn(istom2,jcm1,kc))  &
            + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
            - functn(istom2,jcp2,kc) + functn(istom2,jcm2,kc))
        fdiffc = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
            - functn(istom3,jcp1,kc) + functn(istom3,jcm1,kc))  &
            + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
            - functn(istom3,jcp2,kc) + functn(istom3,jcm2,kc))
        fdiffd = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
            - functn(istom4,jcp1,kc) + functn(istom4,jcm1,kc))  &
            + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
            - functn(istom4,jcp2,kc) + functn(istom4,jcm2,kc))
        fderiv(istol,jc,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
            - functn(istol,jcp1,kc)  + functn(istol,jcm1,kc))  &
            + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
            - functn(istol,jcp2,kc)  + functn(istol,jcm2,kc))
        fdiffb = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
            - functn(istom2,jcp1,kc) + functn(istom2,jcm1,kc))  &
            + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
            - functn(istom2,jcp2,kc) + functn(istom2,jcm2,kc))
        fdiffc = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
            - functn(istom3,jcp1,kc) + functn(istom3,jcm1,kc))  &
            + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
            - functn(istom3,jcp2,kc) + functn(istom3,jcm2,kc))
        fdiffd = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
            - functn(istom4,jcp1,kc) + functn(istom4,jcm1,kc))  &
            + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
            - functn(istom4,jcp2,kc) + functn(istom4,jcm2,kc))
        fderiv(istom1,jc,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      js = 0
      DO jc = jstap2,jstap4
        
        js = js+1
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
        is = 0
        DO ic = istom4,istom2
          
          is = is+1
          icm2 = ic-2
          icm1 = ic-1
          icp1 = ic+1
          icp2 = ic+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(icp1,jcp1,kc) - functn(icp1,jcm1,kc)  &
              - functn(icm1,jcp1,kc) + functn(icm1,jcm1,kc)
          fdiffb = functn(icp2,jcp2,kc) - functn(icp2,jcm2,kc)  &
              - functn(icm2,jcp2,kc) + functn(icm2,jcm2,kc)
          fderiv(ic,jc,kc) = acf3xy*fdiffa + bcf3xy*fdiffb
          fstora(is,js) = fdiffa
          fstorb(is,js) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      js = 1
      DO jc = jstap3,jstap4
        
        jsm1 = js
        js = js+1
        jcm3 = jc-3
        jcp3 = jc+3
        
        is = 0
        DO ic = istom4,istom3
          
          is = is+1
          icm3 = ic-3
          icp3 = ic+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(icp3,jcp3,kc) - functn(icp3,jcm3,kc)  &
              - functn(icm3,jcp3,kc) + functn(icm3,jcm3,kc)
          fderiv(ic,jc,kc) = acf4xy*fstora(is,js) + bcf4xy*fstorb(is,js)  &
              + ccf4xy*fdiffc
          fstorc(is,jsm1) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      js = 3
      is = 1
      jsm1 = 2
      jc = jstap4
      ic = istom4
      jcm4 = jc-4
      jcp4 = jc+4
      icm4 = ic-4
      icp4 = ic+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(icp4,jcp4,kc) - functn(icp4,jcm4,kc)  &
          - functn(icm4,jcp4,kc) + functn(icm4,jcm4,kc)
      fderiv(ic,jc,kc) = acf5xy*fstora(is,js) + bcf5xy*fstorb(is,js)  &
          + ccf5xy*fstorc(is,jsm1) + dcf5xy*fdiffd
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     RH END Y-DIRECTION
!     ==================
IF(nendyr == nbound)THEN
  
!       TAKE SECOND XY-DERIVATIVE IN Y-RIGHT INNER HALO
!       EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    
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
      fdiffa = functn(icp1,jstom3,kc) - functn(icm1,jstom3,kc)  &
          - functn(icp1,jstom5,kc) + functn(icm1,jstom5,kc)
      fdiffb = functn(icp2,jstom2,kc) - functn(icm2,jstom2,kc)  &
          - functn(icp2,jstom6,kc) + functn(icm2,jstom6,kc)
      fdiffc = functn(icp3,jstom1,kc) - functn(icm3,jstom1,kc)  &
          - functn(icp3,jstom7,kc) + functn(icm3,jstom7,kc)
      fdiffd = functn(icp4,jstol,kc)  - functn(icm4,jstol,kc)  &
          - functn(icp4,jstom8,kc) + functn(icm4,jstom8,kc)
      fderiv(ic,jstom4,kc) = acf5xy*fdiffa + bcf5xy*fdiffb  &
          + ccf5xy*fdiffc + dcf5xy*fdiffd
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdiffa = functn(icp1,jstom2,kc) - functn(icm1,jstom2,kc)  &
          - functn(icp1,jstom4,kc) + functn(icm1,jstom4,kc)
      fdiffb = functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
          - functn(icp2,jstom5,kc) + functn(icm2,jstom5,kc)
      fdiffc = functn(icp3,jstol,kc)  - functn(icm3,jstol,kc)  &
          - functn(icp3,jstom6,kc) + functn(icm3,jstom6,kc)
      fderiv(ic,jstom3,kc) = acf4xy*fdiffa + bcf4xy*fdiffb  &
          + ccf4xy*fdiffc
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdiffa = functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
          - functn(icp1,jstom3,kc) + functn(icm1,jstom3,kc)
      fdiffb = functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
          - functn(icp2,jstom4,kc) + functn(icm2,jstom4,kc)
      fderiv(ic,jstom2,kc) = acf3xy*fdiffa + bcf3xy*fdiffb
      
!           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
      fdiffa = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
          - functn(icp1,jstol,kc)  + functn(icm1,jstol,kc))  &
          + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
          - functn(icp2,jstol,kc)  + functn(icm2,jstol,kc))
      fdiffb = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
          - functn(icp1,jstom2,kc) + functn(icm1,jstom2,kc))  &
          + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
          - functn(icp2,jstom2,kc) + functn(icm2,jstom2,kc))
      fdiffc = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
          - functn(icp1,jstom3,kc) + functn(icm1,jstom3,kc))  &
          + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
          - functn(icp2,jstom3,kc) + functn(icm2,jstom3,kc))
      fdiffd = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
          - functn(icp1,jstom4,kc) + functn(icm1,jstom4,kc))  &
          + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
          - functn(icp2,jstom4,kc) + functn(icm2,jstom4,kc))
      fderiv(ic,jstom1,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
          + ccf2xy*fdiffc + dcf2xy*fdiffd
      
!           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
      fdiffa = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
          - functn(icp1,jstom1,kc) + functn(icm1,jstom1,kc))  &
          + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
          - functn(icp2,jstom1,kc) + functn(icm2,jstom1,kc))
      fdiffb = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
          - functn(icp1,jstom2,kc) + functn(icm1,jstom2,kc))  &
          + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
          - functn(icp2,jstom2,kc) + functn(icm2,jstom2,kc))
      fdiffc = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
          - functn(icp1,jstom3,kc) + functn(icm1,jstom3,kc))  &
          + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
          - functn(icp2,jstom3,kc) + functn(icm2,jstom3,kc))
      fdiffd = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
          - functn(icp1,jstom4,kc) + functn(icm1,jstom4,kc))  &
          + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
          - functn(icp2,jstom4,kc) + functn(icm2,jstom4,kc))
      fderiv(ic,jstol,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
    END DO
  END DO
  
!       LH IN X RH IN Y CORNER
!       ======================
  IF(nendxl == nbound)THEN
    
    DO kc = kstal,kstol
      
!           LH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(istap1,jstol,kc) - functn(istap1,jstom1,kc)  &
          - functn(istal,jstol,kc)  + functn(istal,jstom1,kc)
      fdiffb = functn(istap2,jstol,kc) - functn(istap2,jstom2,kc)  &
          - functn(istal,jstol,kc)  + functn(istal,jstom2,kc)
      fdiffc = functn(istap3,jstol,kc) - functn(istap3,jstom3,kc)  &
          - functn(istal,jstol,kc)  + functn(istal,jstom3,kc)
      fdiffd = functn(istap4,jstol,kc) - functn(istap4,jstom4,kc)  &
          - functn(istal,jstol,kc)  + functn(istal,jstom4,kc)
      fderiv(istal,jstol,kc) = acc1xy*fdiffa + bcc1xy*fdiffb  &
          + ccc1xy*fdiffc + dcc1xy*fdiffd
      
!           LH+1 RH-1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(istal,jstom1,kc)  - functn(istal,jstol,kc)  &
          - functn(istap1,jstom1,kc) + functn(istap1,jstol,kc)
      fdiffb = functn(istap2,jstom1,kc) - functn(istap2,jstom2,kc)  &
          - functn(istap1,jstom1,kc) + functn(istap1,jstom2,kc)
      fdiffc = functn(istap3,jstom1,kc) - functn(istap3,jstom3,kc)  &
          - functn(istap1,jstom1,kc) + functn(istap1,jstom3,kc)
      fdiffd = functn(istap4,jstom1,kc) - functn(istap4,jstom4,kc)  &
          - functn(istap1,jstom1,kc) + functn(istap1,jstom4,kc)
      fderiv(istap1,jstom1,kc) = acc2xy*fdiffa + bcc2xy*fdiffb  &
          + ccc2xy*fdiffc + dcc2xy*fdiffd
      
!           LH RH-1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xy*(functn(istap1,jstom1,kc) - functn(istap1,jstol,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstol,kc))  &
          + bcf2xy*(functn(istap1,jstom1,kc) - functn(istap1,jstom2,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom2,kc))  &
          + ccf2xy*(functn(istap1,jstom1,kc) - functn(istap1,jstom3,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom3,kc))  &
          + dcf2xy*(functn(istap1,jstom1,kc) - functn(istap1,jstom4,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom4,kc))
      fdiffb = acf2xy*(functn(istap2,jstom1,kc) - functn(istap2,jstol,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstol,kc))  &
          + bcf2xy*(functn(istap2,jstom1,kc) - functn(istap2,jstom2,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom2,kc))  &
          + ccf2xy*(functn(istap2,jstom1,kc) - functn(istap2,jstom3,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom3,kc))  &
          + dcf2xy*(functn(istap2,jstom1,kc) - functn(istap2,jstom4,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom4,kc))
      fdiffc = acf2xy*(functn(istap3,jstom1,kc) - functn(istap3,jstol,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstol,kc))  &
          + bcf2xy*(functn(istap3,jstom1,kc) - functn(istap3,jstom2,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom2,kc))  &
          + ccf2xy*(functn(istap3,jstom1,kc) - functn(istap3,jstom3,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom3,kc))  &
          + dcf2xy*(functn(istap3,jstom1,kc) - functn(istap3,jstom4,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom4,kc))
      fdiffd = acf2xy*(functn(istap4,jstom1,kc) - functn(istap4,jstol,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstol,kc))  &
          + bcf2xy*(functn(istap4,jstom1,kc) - functn(istap4,jstom2,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom2,kc))  &
          + ccf2xy*(functn(istap4,jstom1,kc) - functn(istap4,jstom3,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom3,kc))  &
          + dcf2xy*(functn(istap4,jstom1,kc) - functn(istap4,jstom4,kc)  &
          - functn(istal,jstom1,kc)  + functn(istal,jstom4,kc))
      fderiv(istal,jstom1,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           LH+1 RH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xy*(functn(istal,jstol,kc)   - functn(istap1,jstol,kc)  &
          - functn(istal,jstom1,kc)  + functn(istap1,jstom1,kc))  &
          + bcf2xy*(functn(istap2,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap2,jstom1,kc) + functn(istap1,jstom1,kc))  &
          + ccf2xy*(functn(istap3,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap3,jstom1,kc) + functn(istap1,jstom1,kc))  &
          + dcf2xy*(functn(istap4,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap4,jstom1,kc) + functn(istap1,jstom1,kc))
      fdiffb = acf2xy*(functn(istal,jstol,kc)   - functn(istap1,jstol,kc)  &
          - functn(istal,jstom2,kc)  + functn(istap1,jstom2,kc))  &
          + bcf2xy*(functn(istap2,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap2,jstom2,kc) + functn(istap1,jstom2,kc))  &
          + ccf2xy*(functn(istap3,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap3,jstom2,kc) + functn(istap1,jstom2,kc))  &
          + dcf2xy*(functn(istap4,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap4,jstom2,kc) + functn(istap1,jstom2,kc))
      fdiffc = acf2xy*(functn(istal,jstol,kc)   - functn(istap1,jstol,kc)  &
          - functn(istal,jstom3,kc)  + functn(istap1,jstom3,kc))  &
          + bcf2xy*(functn(istap2,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap2,jstom3,kc) + functn(istap1,jstom3,kc))  &
          + ccf2xy*(functn(istap3,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap3,jstom3,kc) + functn(istap1,jstom3,kc))  &
          + dcf2xy*(functn(istap4,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap4,jstom3,kc) + functn(istap1,jstom3,kc))
      fdiffd = acf2xy*(functn(istal,jstol,kc)   - functn(istap1,jstol,kc)  &
          - functn(istal,jstom4,kc)  + functn(istap1,jstom4,kc))  &
          + bcf2xy*(functn(istap2,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap2,jstom4,kc) + functn(istap1,jstom4,kc))  &
          + ccf2xy*(functn(istap3,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap3,jstom4,kc) + functn(istap1,jstom4,kc))  &
          + dcf2xy*(functn(istap4,jstol,kc)  - functn(istap1,jstol,kc)  &
          - functn(istap4,jstom4,kc) + functn(istap1,jstom4,kc))
      fderiv(istap1,jstol,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           RH EDGE IN Y
      DO ic = istap2,istap4
        
        icm2 = ic-2
        icm1 = ic-1
        icp1 = ic+1
        icp2 = ic+2
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
            - functn(icp1,jstol,kc)  + functn(icm1,jstol,kc))  &
            + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
            - functn(icp2,jstol,kc)  + functn(icm2,jstol,kc))
        fdiffb = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
            - functn(icp1,jstom2,kc) + functn(icm1,jstom2,kc))  &
            + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
            - functn(icp2,jstom2,kc) + functn(icm2,jstom2,kc))
        fdiffc = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
            - functn(icp1,jstom3,kc) + functn(icm1,jstom3,kc))  &
            + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
            - functn(icp2,jstom3,kc) + functn(icm2,jstom3,kc))
        fdiffd = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
            - functn(icp1,jstom4,kc) + functn(icm1,jstom4,kc))  &
            + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
            - functn(icp2,jstom4,kc) + functn(icm2,jstom4,kc))
        fderiv(ic,jstom1,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
            - functn(icp1,jstom1,kc) + functn(icm1,jstom1,kc))  &
            + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
            - functn(icp2,jstom1,kc) + functn(icm2,jstom1,kc))
        fdiffb = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
            - functn(icp1,jstom2,kc) + functn(icm1,jstom2,kc))  &
            + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
            - functn(icp2,jstom2,kc) + functn(icm2,jstom2,kc))
        fdiffc = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
            - functn(icp1,jstom3,kc) + functn(icm1,jstom3,kc))  &
            + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
            - functn(icp2,jstom3,kc) + functn(icm2,jstom3,kc))
        fdiffd = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
            - functn(icp1,jstom4,kc) + functn(icm1,jstom4,kc))  &
            + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
            - functn(icp2,jstom4,kc) + functn(icm2,jstom4,kc))
        fderiv(ic,jstol,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd
        
      END DO
      
!           LH EDGE IN X
      DO jc = jstom4,jstom2
        
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
!             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofy1*(functn(istap1,jcp1,kc) - functn(istap1,jcm1,kc)  &
            - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
            + bcofy1*(functn(istap1,jcp2,kc) - functn(istap1,jcm2,kc)  &
            - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
        fdiffb = acofy1*(functn(istap2,jcp1,kc) - functn(istap2,jcm1,kc)  &
            - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
            + bcofy1*(functn(istap2,jcp2,kc) - functn(istap2,jcm2,kc)  &
            - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
        fdiffc = acofy1*(functn(istap3,jcp1,kc) - functn(istap3,jcm1,kc)  &
            - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
            + bcofy1*(functn(istap3,jcp2,kc) - functn(istap3,jcm2,kc)  &
            - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
        fdiffd = acofy1*(functn(istap4,jcp1,kc) - functn(istap4,jcm1,kc)  &
            - functn(istal,jcp1,kc)  + functn(istal,jcm1,kc))  &
            + bcofy1*(functn(istap4,jcp2,kc) - functn(istap4,jcm2,kc)  &
            - functn(istal,jcp2,kc)  + functn(istal,jcm2,kc))
        fderiv(istal,jc,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd
        
!             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofy1*(functn(istal,jcp1,kc)  - functn(istal,jcm1,kc)  &
            - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
            + bcofy1*(functn(istal,jcp2,kc)  - functn(istal,jcm2,kc)  &
            - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
        fdiffb = acofy1*(functn(istap2,jcp1,kc) - functn(istap2,jcm1,kc)  &
            - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
            + bcofy1*(functn(istap2,jcp2,kc) - functn(istap2,jcm2,kc)  &
            - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
        fdiffc = acofy1*(functn(istap3,jcp1,kc) - functn(istap3,jcm1,kc)  &
            - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
            + bcofy1*(functn(istap3,jcp2,kc) - functn(istap3,jcm2,kc)  &
            - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
        fdiffd = acofy1*(functn(istap4,jcp1,kc) - functn(istap4,jcm1,kc)  &
            - functn(istap1,jcp1,kc) + functn(istap1,jcm1,kc))  &
            + bcofy1*(functn(istap4,jcp2,kc) - functn(istap4,jcm2,kc)  &
            - functn(istap1,jcp2,kc) + functn(istap1,jcm2,kc))
        fderiv(istap1,jc,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      js = 0
      DO jc = jstom4,jstom2
        
        js = js+1
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
        is = 0
        DO ic = istap2,istap4
          
          is = is+1
          icm2 = ic-2
          icm1 = ic-1
          icp1 = ic+1
          icp2 = ic+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(icp1,jcp1,kc) - functn(icp1,jcm1,kc)  &
              - functn(icm1,jcp1,kc) + functn(icm1,jcm1,kc)
          fdiffb = functn(icp2,jcp2,kc) - functn(icp2,jcm2,kc)  &
              - functn(icm2,jcp2,kc) + functn(icm2,jcm2,kc)
          fderiv(ic,jc,kc) = acf3xy*fdiffa + bcf3xy*fdiffb
          fstora(is,js) = fdiffa
          fstorb(is,js) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      js = 0
      DO jc = jstom4,jstom3
        
        js = js+1
        jcm3 = jc-3
        jcp3 = jc+3
        
        is = 1
        DO ic = istap3,istap4
          
          ism1 = is
          is = is+1
          icm3 = ic-3
          icp3 = ic+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(icp3,jcp3,kc) - functn(icp3,jcm3,kc)  &
              - functn(icm3,jcp3,kc) + functn(icm3,jcm3,kc)
          fderiv(ic,jc,kc) = acf4xy*fstora(is,js) + bcf4xy*fstorb(is,js)  &
              + ccf4xy*fdiffc
          fstorc(ism1,js) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      js = 1
      is = 3
      ism1 = 2
      jc = jstom4
      ic = istap4
      jcm4 = jc-4
      jcp4 = jc+4
      icm4 = ic-4
      icp4 = ic+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(icp4,jcp4,kc) - functn(icp4,jcm4,kc)  &
          - functn(icm4,jcp4,kc) + functn(icm4,jcm4,kc)
      fderiv(ic,jc,kc) = acf5xy*fstora(is,js) + bcf5xy*fstorb(is,js)  &
          + ccf5xy*fstorc(ism1,js) + dcf5xy*fdiffd
      
    END DO
    
  END IF
  
!       RH IN X RH IN Y CORNER
!       ======================
  IF(nendxr == nbound)THEN
    
    DO kc = kstal,kstol
      
!           RH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
      fdiffa = functn(istom1,jstom1,kc) - functn(istom1,jstol,kc)  &
          - functn(istol,jstom1,kc)  + functn(istol,jstol,kc)
      fdiffb = functn(istom2,jstom2,kc) - functn(istom2,jstol,kc)  &
          - functn(istol,jstom2,kc)  + functn(istol,jstol,kc)
      fdiffc = functn(istom3,jstom3,kc) - functn(istom3,jstol,kc)  &
          - functn(istol,jstom3,kc)  + functn(istol,jstol,kc)
      fdiffd = functn(istom4,jstom4,kc) - functn(istom4,jstol,kc)  &
          - functn(istol,jstom4,kc)  + functn(istol,jstol,kc)
      fderiv(istol,jstol,kc) = acc1xy*fdiffa + bcc1xy*fdiffb  &
          + ccc1xy*fdiffc + dcc1xy*fdiffd
      
!           RH-1 RH-1 CORNER POINT: 4TH ORDER MIXED
      fdiffa = functn(istol,jstol,kc)   - functn(istol,jstom1,kc)  &
          - functn(istom1,jstol,kc)  + functn(istom1,jstom1,kc)
      fdiffb = functn(istom2,jstom2,kc) - functn(istom2,jstom1,kc)  &
          - functn(istom1,jstom2,kc) + functn(istom1,jstom1,kc)
      fdiffc = functn(istom3,jstom3,kc) - functn(istom3,jstom1,kc)  &
          - functn(istom1,jstom3,kc) + functn(istom1,jstom1,kc)
      fdiffd = functn(istom4,jstom4,kc) - functn(istom4,jstom1,kc)  &
          - functn(istom1,jstom4,kc) + functn(istom1,jstom1,kc)
      fderiv(istom1,jstom1,kc) = acc2xy*fdiffa + bcc2xy*fdiffb  &
          + ccc2xy*fdiffc + dcc2xy*fdiffd
      
!           RH RH-1 EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xy*(functn(istom1,jstol,kc)  - functn(istom1,jstom1,kc)  &
          - functn(istol,jstol,kc)   + functn(istol,jstom1,kc))  &
          + bcf2xy*(functn(istom1,jstom2,kc) - functn(istom1,jstom1,kc)  &
          - functn(istol,jstom2,kc)  + functn(istol,jstom1,kc))  &
          + ccf2xy*(functn(istom1,jstom3,kc) - functn(istom1,jstom1,kc)  &
          - functn(istol,jstom3,kc)  + functn(istol,jstom1,kc))  &
          + dcf2xy*(functn(istom1,jstom4,kc) - functn(istom1,jstom1,kc)  &
          - functn(istol,jstom4,kc)  + functn(istol,jstom1,kc))
      fdiffb = acf2xy*(functn(istom2,jstol,kc)  - functn(istom2,jstom1,kc)  &
          - functn(istol,jstol,kc)   + functn(istol,jstom1,kc))  &
          + bcf2xy*(functn(istom2,jstom2,kc) - functn(istom2,jstom1,kc)  &
          - functn(istol,jstom2,kc)  + functn(istol,jstom1,kc))  &
          + ccf2xy*(functn(istom2,jstom3,kc) - functn(istom2,jstom1,kc)  &
          - functn(istol,jstom3,kc)  + functn(istol,jstom1,kc))  &
          + dcf2xy*(functn(istom2,jstom4,kc) - functn(istom2,jstom1,kc)  &
          - functn(istol,jstom4,kc)  + functn(istol,jstom1,kc))
      fdiffc = acf2xy*(functn(istom3,jstol,kc)  - functn(istom3,jstom1,kc)  &
          - functn(istol,jstol,kc)   + functn(istol,jstom1,kc))  &
          + bcf2xy*(functn(istom3,jstom2,kc) - functn(istom3,jstom1,kc)  &
          - functn(istol,jstom2,kc)  + functn(istol,jstom1,kc))  &
          + ccf2xy*(functn(istom3,jstom3,kc) - functn(istom3,jstom1,kc)  &
          - functn(istol,jstom3,kc)  + functn(istol,jstom1,kc))  &
          + dcf2xy*(functn(istom3,jstom4,kc) - functn(istom3,jstom1,kc)  &
          - functn(istol,jstom4,kc)  + functn(istol,jstom1,kc))
      fdiffd = acf2xy*(functn(istom4,jstol,kc)  - functn(istom4,jstom1,kc)  &
          - functn(istol,jstol,kc)   + functn(istol,jstom1,kc))  &
          + bcf2xy*(functn(istom4,jstom2,kc) - functn(istom4,jstom1,kc)  &
          - functn(istol,jstom2,kc)  + functn(istol,jstom1,kc))  &
          + ccf2xy*(functn(istom4,jstom3,kc) - functn(istom4,jstom1,kc)  &
          - functn(istol,jstom3,kc)  + functn(istol,jstom1,kc))  &
          + dcf2xy*(functn(istom4,jstom4,kc) - functn(istom4,jstom1,kc)  &
          - functn(istol,jstom4,kc)  + functn(istol,jstom1,kc))
      fderiv(istol,jstom1,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           RH+1 RH EDGE POINT: 4TH ORDER MIXED
      fdiffa = acf2xy*(functn(istol,jstom1,kc)  - functn(istom1,jstom1,kc)  &
          - functn(istol,jstol,kc)   + functn(istom1,jstol,kc))  &
          + bcf2xy*(functn(istom2,jstom1,kc) - functn(istom1,jstom1,kc)  &
          - functn(istom2,jstol,kc)  + functn(istom1,jstol,kc))  &
          + ccf2xy*(functn(istom3,jstom1,kc) - functn(istom1,jstom1,kc)  &
          - functn(istom3,jstol,kc)  + functn(istom1,jstol,kc))  &
          + dcf2xy*(functn(istom4,jstom1,kc) - functn(istom1,jstom1,kc)  &
          - functn(istom4,jstol,kc)  + functn(istom1,jstol,kc))
      fdiffb = acf2xy*(functn(istol,jstom2,kc)  - functn(istom1,jstom2,kc)  &
          - functn(istol,jstol,kc)   + functn(istom1,jstol,kc))  &
          + bcf2xy*(functn(istom2,jstom2,kc) - functn(istom1,jstom2,kc)  &
          - functn(istom2,jstol,kc)  + functn(istom1,jstol,kc))  &
          + ccf2xy*(functn(istom3,jstom2,kc) - functn(istom1,jstom2,kc)  &
          - functn(istom3,jstol,kc)  + functn(istom1,jstol,kc))  &
          + dcf2xy*(functn(istom4,jstom2,kc) - functn(istom1,jstom2,kc)  &
          - functn(istom4,jstol,kc)  + functn(istom1,jstol,kc))
      fdiffc = acf2xy*(functn(istol,jstom3,kc)  - functn(istom1,jstom3,kc)  &
          - functn(istol,jstol,kc)   + functn(istom1,jstol,kc))  &
          + bcf2xy*(functn(istom2,jstom3,kc) - functn(istom1,jstom3,kc)  &
          - functn(istom2,jstol,kc)  + functn(istom1,jstol,kc))  &
          + ccf2xy*(functn(istom3,jstom3,kc) - functn(istom1,jstom3,kc)  &
          - functn(istom3,jstol,kc)  + functn(istom1,jstol,kc))  &
          + dcf2xy*(functn(istom4,jstom3,kc) - functn(istom1,jstom3,kc)  &
          - functn(istom4,jstol,kc)  + functn(istom1,jstol,kc))
      fdiffd = acf2xy*(functn(istol,jstom4,kc)  - functn(istom1,jstom4,kc)  &
          - functn(istol,jstol,kc)   + functn(istom1,jstol,kc))  &
          + bcf2xy*(functn(istom2,jstom4,kc) - functn(istom1,jstom4,kc)  &
          - functn(istom2,jstol,kc)  + functn(istom1,jstol,kc))  &
          + ccf2xy*(functn(istom3,jstom4,kc) - functn(istom1,jstom4,kc)  &
          - functn(istom3,jstol,kc)  + functn(istom1,jstol,kc))  &
          + dcf2xy*(functn(istom4,jstom4,kc) - functn(istom1,jstom4,kc)  &
          - functn(istom4,jstol,kc)  + functn(istom1,jstol,kc))
      fderiv(istom1,jstol,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
          + ccf1xy*fdiffc + dcf1xy*fdiffd
      
!           RH EDGE IN Y
      DO ic = istom4,istom2
        
        icm2 = ic-2
        icm1 = ic-1
        icp1 = ic+1
        icp2 = ic+2
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
            - functn(icp1,jstol,kc)  + functn(icm1,jstol,kc))  &
            + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
            - functn(icp2,jstol,kc)  + functn(icm2,jstol,kc))
        fdiffb = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
            - functn(icp1,jstom2,kc) + functn(icm1,jstom2,kc))  &
            + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
            - functn(icp2,jstom2,kc) + functn(icm2,jstom2,kc))
        fdiffc = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
            - functn(icp1,jstom3,kc) + functn(icm1,jstom3,kc))  &
            + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
            - functn(icp2,jstom3,kc) + functn(icm2,jstom3,kc))
        fdiffd = acofx1*(functn(icp1,jstom1,kc) - functn(icm1,jstom1,kc)  &
            - functn(icp1,jstom4,kc) + functn(icm1,jstom4,kc))  &
            + bcofx1*(functn(icp2,jstom1,kc) - functn(icm2,jstom1,kc)  &
            - functn(icp2,jstom4,kc) + functn(icm2,jstom4,kc))
        fderiv(ic,jstom1,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
            - functn(icp1,jstom1,kc) + functn(icm1,jstom1,kc))  &
            + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
            - functn(icp2,jstom1,kc) + functn(icm2,jstom1,kc))
        fdiffb = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
            - functn(icp1,jstom2,kc) + functn(icm1,jstom2,kc))  &
            + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
            - functn(icp2,jstom2,kc) + functn(icm2,jstom2,kc))
        fdiffc = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
            - functn(icp1,jstom3,kc) + functn(icm1,jstom3,kc))  &
            + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
            - functn(icp2,jstom3,kc) + functn(icm2,jstom3,kc))
        fdiffd = acofx1*(functn(icp1,jstol,kc)  - functn(icm1,jstol,kc)  &
            - functn(icp1,jstom4,kc) + functn(icm1,jstom4,kc))  &
            + bcofx1*(functn(icp2,jstol,kc)  - functn(icm2,jstol,kc)  &
            - functn(icp2,jstom4,kc) + functn(icm2,jstom4,kc))
        fderiv(ic,jstol,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd
        
      END DO
      
!           RH EDGE IN X
      DO jc = jstom4,jstom2
        
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
!             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
        fdiffa = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
            - functn(istom1,jcp1,kc) + functn(istom1,jcm1,kc))  &
            + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
            - functn(istom1,jcp2,kc) + functn(istom1,jcm2,kc))
        fdiffb = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
            - functn(istom2,jcp1,kc) + functn(istom2,jcm1,kc))  &
            + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
            - functn(istom2,jcp2,kc) + functn(istom2,jcm2,kc))
        fdiffc = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
            - functn(istom3,jcp1,kc) + functn(istom3,jcm1,kc))  &
            + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
            - functn(istom3,jcp2,kc) + functn(istom3,jcm2,kc))
        fdiffd = acofy1*(functn(istol,jcp1,kc)  - functn(istol,jcm1,kc)  &
            - functn(istom4,jcp1,kc) + functn(istom4,jcm1,kc))  &
            + bcofy1*(functn(istol,jcp2,kc)  - functn(istol,jcm2,kc)  &
            - functn(istom4,jcp2,kc) + functn(istom4,jcm2,kc))
        fderiv(istol,jc,kc) = acf1xy*fdiffa + bcf1xy*fdiffb  &
            + ccf1xy*fdiffc + dcf1xy*fdiffd
        
!             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
        fdiffa = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
            - functn(istol,jcp1,kc)  + functn(istol,jcm1,kc))  &
            + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
            - functn(istol,jcp2,kc)  + functn(istol,jcm2,kc))
        fdiffb = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
            - functn(istom2,jcp1,kc) + functn(istom2,jcm1,kc))  &
            + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
            - functn(istom2,jcp2,kc) + functn(istom2,jcm2,kc))
        fdiffc = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
            - functn(istom3,jcp1,kc) + functn(istom3,jcm1,kc))  &
            + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
            - functn(istom3,jcp2,kc) + functn(istom3,jcm2,kc))
        fdiffd = acofy1*(functn(istom1,jcp1,kc) - functn(istom1,jcm1,kc)  &
            - functn(istom4,jcp1,kc) + functn(istom4,jcm1,kc))  &
            + bcofy1*(functn(istom1,jcp2,kc) - functn(istom1,jcm2,kc)  &
            - functn(istom4,jcp2,kc) + functn(istom4,jcm2,kc))
        fderiv(istom1,jc,kc) = acf2xy*fdiffa + bcf2xy*fdiffb  &
            + ccf2xy*fdiffc + dcf2xy*fdiffd
        
      END DO
      
!           INTERIOR POINTS 4TH ORDER
      js = 0
      DO jc = jstom4,jstom2
        
        js = js+1
        jcm2 = jc-2
        jcm1 = jc-1
        jcp1 = jc+1
        jcp2 = jc+2
        
        is = 0
        DO ic = istom4,istom2
          
          is = is+1
          icm2 = ic-2
          icm1 = ic-1
          icp1 = ic+1
          icp2 = ic+2
          
!               4TH ORDER CENTRED
          fdiffa = functn(icp1,jcp1,kc) - functn(icp1,jcm1,kc)  &
              - functn(icm1,jcp1,kc) + functn(icm1,jcm1,kc)
          fdiffb = functn(icp2,jcp2,kc) - functn(icp2,jcm2,kc)  &
              - functn(icm2,jcp2,kc) + functn(icm2,jcm2,kc)
          fderiv(ic,jc,kc) = acf3xy*fdiffa + bcf3xy*fdiffb
          fstora(is,js) = fdiffa
          fstorb(is,js) = fdiffb
          
        END DO
      END DO
      
!           INTERIOR POINTS 6TH ORDER
      js = 0
      DO jc = jstom4,jstom3
        
        js = js+1
        jcm3 = jc-3
        jcp3 = jc+3
        
        is = 0
        DO ic = istom4,istom3
          
          is = is+1
          icm3 = ic-3
          icp3 = ic+3
          
!               6TH ORDER CENTRED
          fdiffc = functn(icp3,jcp3,kc) - functn(icp3,jcm3,kc)  &
              - functn(icm3,jcp3,kc) + functn(icm3,jcm3,kc)
          fderiv(ic,jc,kc) = acf4xy*fstora(is,js) + bcf4xy*fstorb(is,js)  &
              + ccf4xy*fdiffc
          fstorc(is,js) = fdiffc
          
        END DO
      END DO
      
!           INTERIOR POINT 8TH ORDER
      js = 1
      is = 1
      jc = jstom4
      ic = istom4
      jcm4 = jc-4
      jcp4 = jc+4
      icm4 = ic-4
      icp4 = ic+4
      
!           8TH ORDER CENTRED
      fdiffd = functn(icp4,jcp4,kc) - functn(icp4,jcm4,kc)  &
          - functn(icm4,jcp4,kc) + functn(icm4,jcm4,kc)
      fderiv(ic,jc,kc) = acf5xy*fstora(is,js) + bcf5xy*fstorb(is,js)  &
          + ccf5xy*fstorc(is,js) + dcf5xy*fdiffd
      
    END DO
    
  END IF
  
END IF

!     =========================================================================

!     SCALING
!     =======
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fderiv(ic,jc,kc) = fderiv(ic,jc,kc)*ovdelx*ovdely
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE d2fdxy
