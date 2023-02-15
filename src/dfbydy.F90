SUBROUTINE dfbydy(functn,fderiv)

! Code converted using TO_F90 by Alan Miller
! Date: 2022-11-09  Time: 12:26:49

!     *************************************************************************

!     DFBYDY
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     11-APR-2003:  RSC MODIFIED FOR SENGA2

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES FIRST Y-DERIVATIVE OF SPECIFIED FUNCTION
!     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
!     EXPLICIT 8TH,6TH,4TH,4TH ORDER END CONDITIONS

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
real(kind=dp) :: fdiffa,fdiffb,fdiffc,fdiffd,fdiffe
INTEGER :: ic,jc,kc
INTEGER :: jstart,jfinis
INTEGER :: jcm5,jcm4,jcm3,jcm2,jcm1,jccc,jcp1,jcp2,jcp3,jcp4,jcp5


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

jstart = jstal
jfinis = jstol
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
    
    DO ic = istal,istol
      
      fdiffa = functn(ic,jcp1,kc) - functn(ic,jcm1,kc)
      fdiffb = functn(ic,jcp2,kc) - functn(ic,jcm2,kc)
      fdiffc = functn(ic,jcp3,kc) - functn(ic,jcm3,kc)
      fdiffd = functn(ic,jcp4,kc) - functn(ic,jcm4,kc)
      fdiffe = functn(ic,jcp5,kc) - functn(ic,jcm5,kc)
      
      fderiv(ic,jc,kc) = acoffy*fdiffa + bcoffy*fdiffb  &
                       + ccoffy*fdiffc + dcoffy*fdiffd  &
                       + ecoffy*fdiffe
      
    END DO
    
  END DO
  
END DO

!     =========================================================================

!     LH END
!     ======
IF(nendyl == nbound)THEN
  
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    DO ic = istal,istol
      
!           LH POINT: 4TH ORDER ONE-SIDED
      fdiffa = functn(ic,jstap1,kc) - functn(ic,jstal,kc)
      fdiffb = functn(ic,jstap2,kc) - functn(ic,jstal,kc)
      fdiffc = functn(ic,jstap3,kc) - functn(ic,jstal,kc)
      fdiffd = functn(ic,jstap4,kc) - functn(ic,jstal,kc)
      fderiv(ic,jstal,kc) = acof1y*fdiffa + bcof1y*fdiffb  &
          + ccof1y*fdiffc + dcof1y*fdiffd
      
!           LH POINT PLUS 1: 4TH ORDER MIXED
      fdiffa = functn(ic,jstal,kc)  - functn(ic,jstap1,kc)
      fdiffb = functn(ic,jstap2,kc) - functn(ic,jstap1,kc)
      fdiffc = functn(ic,jstap3,kc) - functn(ic,jstap1,kc)
      fdiffd = functn(ic,jstap4,kc) - functn(ic,jstap1,kc)
      fderiv(ic,jstap1,kc) = acof2y*fdiffa + bcof2y*fdiffb  &
          + ccof2y*fdiffc + dcof2y*fdiffd
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdiffa = functn(ic,jstap3,kc) - functn(ic,jstap1,kc)
      fdiffb = functn(ic,jstap4,kc) - functn(ic,jstal,kc)
      fderiv(ic,jstap2,kc) = acof3y*fdiffa + bcof3y*fdiffb
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdiffa = functn(ic,jstap4,kc) - functn(ic,jstap2,kc)
      fdiffb = functn(ic,jstap5,kc) - functn(ic,jstap1,kc)
      fdiffc = functn(ic,jstap6,kc) - functn(ic,jstal,kc)
      fderiv(ic,jstap3,kc) = acof4y*fdiffa + bcof4y*fdiffb  &
          + ccof4y*fdiffc
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdiffa = functn(ic,jstap5,kc) - functn(ic,jstap3,kc)
      fdiffb = functn(ic,jstap6,kc) - functn(ic,jstap2,kc)
      fdiffc = functn(ic,jstap7,kc) - functn(ic,jstap1,kc)
      fdiffd = functn(ic,jstap8,kc) - functn(ic,jstal,kc)
      fderiv(ic,jstap4,kc) = acof5y*fdiffa + bcof5y*fdiffb  &
          + ccof5y*fdiffc + dcof5y*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     RH END
!     ======
IF(nendyr == nbound)THEN
  
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    DO ic = istal,istol
      
!           RH POINT MINUS 4: 8TH ORDER CENTRED
      fdiffa = functn(ic,jstom3,kc) - functn(ic,jstom5,kc)
      fdiffb = functn(ic,jstom2,kc) - functn(ic,jstom6,kc)
      fdiffc = functn(ic,jstom1,kc) - functn(ic,jstom7,kc)
      fdiffd = functn(ic,jstol,kc)  - functn(ic,jstom8,kc)
      fderiv(ic,jstom4,kc) = acof5y*fdiffa + bcof5y*fdiffb  &
          + ccof5y*fdiffc + dcof5y*fdiffd
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdiffa = functn(ic,jstom2,kc) - functn(ic,jstom4,kc)
      fdiffb = functn(ic,jstom1,kc) - functn(ic,jstom5,kc)
      fdiffc = functn(ic,jstol,kc)  - functn(ic,jstom6,kc)
      fderiv(ic,jstom3,kc) = acof4y*fdiffa + bcof4y*fdiffb  &
          + ccof4y*fdiffc
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdiffa = functn(ic,jstom1,kc) - functn(ic,jstom3,kc)
      fdiffb = functn(ic,jstol,kc)  - functn(ic,jstom4,kc)
      fderiv(ic,jstom2,kc) = acof3y*fdiffa + bcof3y*fdiffb
      
!           RH POINT MINUS 1: 4TH ORDER MIXED
      fdiffa = functn(ic,jstom1,kc) - functn(ic,jstol,kc)
      fdiffb = functn(ic,jstom1,kc) - functn(ic,jstom2,kc)
      fdiffc = functn(ic,jstom1,kc) - functn(ic,jstom3,kc)
      fdiffd = functn(ic,jstom1,kc) - functn(ic,jstom4,kc)
      fderiv(ic,jstom1,kc) = acof2y*fdiffa + bcof2y*fdiffb  &
          + ccof2y*fdiffc + dcof2y*fdiffd
      
!           RH POINT: 4TH ORDER ONE-SIDED
      fdiffa = functn(ic,jstol,kc) - functn(ic,jstom1,kc)
      fdiffb = functn(ic,jstol,kc) - functn(ic,jstom2,kc)
      fdiffc = functn(ic,jstol,kc) - functn(ic,jstom3,kc)
      fdiffd = functn(ic,jstol,kc) - functn(ic,jstom4,kc)
      fderiv(ic,jstol,kc) = acof1y*fdiffa + bcof1y*fdiffb  &
          + ccof1y*fdiffc + dcof1y*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     SCALING
!     =======
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fderiv(ic,jc,kc) = fderiv(ic,jc,kc)*ovdely
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE dfbydy
