SUBROUTINE dfbydx(functn,fderiv)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 12:57:47

!     *************************************************************************

!     DFBYDX
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     28-MAR-2003:  RSC MODIFIED FOR SENGA2

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES FIRST X-DERIVATIVE OF SPECIFIED FUNCTION
!     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
!     EXPLICIT 8TH,6TH,4TH,4TH ORDER END CONDITIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

real(kind=8), INTENT(IN)             :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
real(kind=8), INTENT(OUT)            :: fderiv(nxsize,nysize,nzsize)




!     LOCAL DATA
!     ==========
real(kind=8) :: fdiffa,fdiffb,fdiffc,fdiffd,fdiffe
INTEGER :: ic,jc,kc
INTEGER :: istart,ifinis
INTEGER :: icm5,icm4,icm3,icm2,icm1,iccc,icp1,icp2,icp3,icp4,icp5


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

istart = istal
ifinis = istol
IF(nendxl == nbound)istart = istap5
IF(nendxr == nbound)ifinis = istom5

!     =========================================================================

!     INTERIOR SCHEME
!     ===============

!     TENTH ORDER EXPLICIT DIFFERENCES
DO kc = kstal,kstol
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
      
      fdiffa = functn(icp1,jc,kc) - functn(icm1,jc,kc)
      fdiffb = functn(icp2,jc,kc) - functn(icm2,jc,kc)
      fdiffc = functn(icp3,jc,kc) - functn(icm3,jc,kc)
      fdiffd = functn(icp4,jc,kc) - functn(icm4,jc,kc)
      fdiffe = functn(icp5,jc,kc) - functn(icm5,jc,kc)
      
      fderiv(ic,jc,kc) = acoffx*fdiffa + bcoffx*fdiffb  &
          + ccoffx*fdiffc + dcoffx*fdiffd  &
          + ecoffx*fdiffe
      
    END DO
    
  END DO
END DO

!     =========================================================================

!     LH END
!     ======
IF(nendxl == nbound)THEN
  
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
!           LH POINT: 4TH ORDER ONE-SIDED
      fdiffa = functn(istap1,jc,kc) - functn(istal,jc,kc)
      fdiffb = functn(istap2,jc,kc) - functn(istal,jc,kc)
      fdiffc = functn(istap3,jc,kc) - functn(istal,jc,kc)
      fdiffd = functn(istap4,jc,kc) - functn(istal,jc,kc)
      fderiv(istal,jc,kc) = acof1x*fdiffa + bcof1x*fdiffb  &
          + ccof1x*fdiffc + dcof1x*fdiffd
      
!           LH POINT PLUS 1: 4TH ORDER MIXED
      fdiffa = functn(istal,jc,kc)  - functn(istap1,jc,kc)
      fdiffb = functn(istap2,jc,kc) - functn(istap1,jc,kc)
      fdiffc = functn(istap3,jc,kc) - functn(istap1,jc,kc)
      fdiffd = functn(istap4,jc,kc) - functn(istap1,jc,kc)
      fderiv(istap1,jc,kc) = acof2x*fdiffa + bcof2x*fdiffb  &
          + ccof2x*fdiffc + dcof2x*fdiffd
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdiffa = functn(istap3,jc,kc) - functn(istap1,jc,kc)
      fdiffb = functn(istap4,jc,kc) - functn(istal,jc,kc)
      fderiv(istap2,jc,kc) = acof3x*fdiffa + bcof3x*fdiffb
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdiffa = functn(istap4,jc,kc) - functn(istap2,jc,kc)
      fdiffb = functn(istap5,jc,kc) - functn(istap1,jc,kc)
      fdiffc = functn(istap6,jc,kc) - functn(istal,jc,kc)
      fderiv(istap3,jc,kc) = acof4x*fdiffa + bcof4x*fdiffb  &
          + ccof4x*fdiffc
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdiffa = functn(istap5,jc,kc) - functn(istap3,jc,kc)
      fdiffb = functn(istap6,jc,kc) - functn(istap2,jc,kc)
      fdiffc = functn(istap7,jc,kc) - functn(istap1,jc,kc)
      fdiffd = functn(istap8,jc,kc) - functn(istal,jc,kc)
      fderiv(istap4,jc,kc) = acof5x*fdiffa + bcof5x*fdiffb  &
          + ccof5x*fdiffc + dcof5x*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     RH END
!     ======
IF(nendxr == nbound)THEN
  
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
!           RH POINT MINUS 4: 8TH ORDER CENTRED
      fdiffa = functn(istom3,jc,kc) - functn(istom5,jc,kc)
      fdiffb = functn(istom2,jc,kc) - functn(istom6,jc,kc)
      fdiffc = functn(istom1,jc,kc) - functn(istom7,jc,kc)
      fdiffd = functn(istol,jc,kc)  - functn(istom8,jc,kc)
      fderiv(istom4,jc,kc) = acof5x*fdiffa + bcof5x*fdiffb  &
          + ccof5x*fdiffc + dcof5x*fdiffd
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdiffa = functn(istom2,jc,kc) - functn(istom4,jc,kc)
      fdiffb = functn(istom1,jc,kc) - functn(istom5,jc,kc)
      fdiffc = functn(istol,jc,kc)  - functn(istom6,jc,kc)
      fderiv(istom3,jc,kc) = acof4x*fdiffa + bcof4x*fdiffb  &
          + ccof4x*fdiffc
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdiffa = functn(istom1,jc,kc) - functn(istom3,jc,kc)
      fdiffb = functn(istol,jc,kc)  - functn(istom4,jc,kc)
      fderiv(istom2,jc,kc) = acof3x*fdiffa + bcof3x*fdiffb
      
!           RH POINT MINUS 1: 4TH ORDER MIXED
      fdiffa = functn(istom1,jc,kc) - functn(istol,jc,kc)
      fdiffb = functn(istom1,jc,kc) - functn(istom2,jc,kc)
      fdiffc = functn(istom1,jc,kc) - functn(istom3,jc,kc)
      fdiffd = functn(istom1,jc,kc) - functn(istom4,jc,kc)
      fderiv(istom1,jc,kc) = acof2x*fdiffa + bcof2x*fdiffb  &
          + ccof2x*fdiffc + dcof2x*fdiffd
      
!           RH POINT: 4TH ORDER ONE-SIDED
      fdiffa = functn(istol,jc,kc) - functn(istom1,jc,kc)
      fdiffb = functn(istol,jc,kc) - functn(istom2,jc,kc)
      fdiffc = functn(istol,jc,kc) - functn(istom3,jc,kc)
      fdiffd = functn(istol,jc,kc) - functn(istom4,jc,kc)
      fderiv(istol,jc,kc) = acof1x*fdiffa + bcof1x*fdiffb  &
          + ccof1x*fdiffc + dcof1x*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     SCALING
!     =======
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fderiv(ic,jc,kc) = fderiv(ic,jc,kc)*ovdelx
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE dfbydx
