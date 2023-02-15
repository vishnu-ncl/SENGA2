SUBROUTINE dfbydz(functn,fderiv)

! Code converted using TO_F90 by Alan Miller
! Date: 2022-11-09  Time: 13:11:08

!     *************************************************************************

!     DFBYDZ
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
!     EVALUATES FIRST Z-DERIVATIVE OF SPECIFIED FUNCTION
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
INTEGER :: kstart,kfinis
INTEGER :: kcm5,kcm4,kcm3,kcm2,kcm1,kccc,kcp1,kcp2,kcp3,kcp4,kcp5


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

kstart = kstal
kfinis = kstol
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
    
    DO ic = istal,istol
      
      fdiffa = functn(ic,jc,kcp1) - functn(ic,jc,kcm1)
      fdiffb = functn(ic,jc,kcp2) - functn(ic,jc,kcm2)
      fdiffc = functn(ic,jc,kcp3) - functn(ic,jc,kcm3)
      fdiffd = functn(ic,jc,kcp4) - functn(ic,jc,kcm4)
      fdiffe = functn(ic,jc,kcp5) - functn(ic,jc,kcm5)
      
      fderiv(ic,jc,kc) = acoffz*fdiffa + bcoffz*fdiffb  &
          + ccoffz*fdiffc + dcoffz*fdiffd  &
          + ecoffz*fdiffe
      
    END DO
    
  END DO
END DO

!     =========================================================================

!     LH END
!     ======
IF(nendzl == nbound)THEN
  
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
  DO jc = jstal,jstol
    DO ic = istal,istol
      
!           LH POINT: 4TH ORDER ONE-SIDED
      fdiffa = functn(ic,jc,kstap1) - functn(ic,jc,kstal)
      fdiffb = functn(ic,jc,kstap2) - functn(ic,jc,kstal)
      fdiffc = functn(ic,jc,kstap3) - functn(ic,jc,kstal)
      fdiffd = functn(ic,jc,kstap4) - functn(ic,jc,kstal)
      fderiv(ic,jc,kstal) = acof1z*fdiffa + bcof1z*fdiffb  &
          + ccof1z*fdiffc + dcof1z*fdiffd
      
!           LH POINT PLUS 1: 4TH ORDER MIXED
      fdiffa = functn(ic,jc,kstal)  - functn(ic,jc,kstap1)
      fdiffb = functn(ic,jc,kstap2) - functn(ic,jc,kstap1)
      fdiffc = functn(ic,jc,kstap3) - functn(ic,jc,kstap1)
      fdiffd = functn(ic,jc,kstap4) - functn(ic,jc,kstap1)
      fderiv(ic,jc,kstap1) = acof2z*fdiffa + bcof2z*fdiffb  &
          + ccof2z*fdiffc + dcof2z*fdiffd
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdiffa = functn(ic,jc,kstap3) - functn(ic,jc,kstap1)
      fdiffb = functn(ic,jc,kstap4) - functn(ic,jc,kstal)
      fderiv(ic,jc,kstap2) = acof3z*fdiffa + bcof3z*fdiffb
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdiffa = functn(ic,jc,kstap4) - functn(ic,jc,kstap2)
      fdiffb = functn(ic,jc,kstap5) - functn(ic,jc,kstap1)
      fdiffc = functn(ic,jc,kstap6) - functn(ic,jc,kstal)
      fderiv(ic,jc,kstap3) = acof4z*fdiffa + bcof4z*fdiffb  &
          + ccof4z*fdiffc
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdiffa = functn(ic,jc,kstap5) - functn(ic,jc,kstap3)
      fdiffb = functn(ic,jc,kstap6) - functn(ic,jc,kstap2)
      fdiffc = functn(ic,jc,kstap7) - functn(ic,jc,kstap1)
      fdiffd = functn(ic,jc,kstap8) - functn(ic,jc,kstal)
      fderiv(ic,jc,kstap4) = acof5z*fdiffa + bcof5z*fdiffb  &
          + ccof5z*fdiffc + dcof5z*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     RH END
!     ======
IF(nendzr == nbound)THEN
  
!       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
  DO jc = jstal,jstol
    DO ic = istal,istol
      
!           RH POINT MINUS 4: 8TH ORDER CENTRED
      fdiffa = functn(ic,jc,kstom3) - functn(ic,jc,kstom5)
      fdiffb = functn(ic,jc,kstom2) - functn(ic,jc,kstom6)
      fdiffc = functn(ic,jc,kstom1) - functn(ic,jc,kstom7)
      fdiffd = functn(ic,jc,kstol)  - functn(ic,jc,kstom8)
      fderiv(ic,jc,kstom4) = acof5z*fdiffa + bcof5z*fdiffb  &
          + ccof5z*fdiffc + dcof5z*fdiffd
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdiffa = functn(ic,jc,kstom2) - functn(ic,jc,kstom4)
      fdiffb = functn(ic,jc,kstom1) - functn(ic,jc,kstom5)
      fdiffc = functn(ic,jc,kstol)  - functn(ic,jc,kstom6)
      fderiv(ic,jc,kstom3) = acof4z*fdiffa + bcof4z*fdiffb  &
          + ccof4z*fdiffc
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdiffa = functn(ic,jc,kstom1) - functn(ic,jc,kstom3)
      fdiffb = functn(ic,jc,kstol)  - functn(ic,jc,kstom4)
      fderiv(ic,jc,kstom2) = acof3z*fdiffa + bcof3z*fdiffb
      
!           RH POINT MINUS 1: 4TH ORDER MIXED
      fdiffa = functn(ic,jc,kstom1) - functn(ic,jc,kstol)
      fdiffb = functn(ic,jc,kstom1) - functn(ic,jc,kstom2)
      fdiffc = functn(ic,jc,kstom1) - functn(ic,jc,kstom3)
      fdiffd = functn(ic,jc,kstom1) - functn(ic,jc,kstom4)
      fderiv(ic,jc,kstom1) = acof2z*fdiffa + bcof2z*fdiffb  &
          + ccof2z*fdiffc + dcof2z*fdiffd
      
!           RH POINT: 4TH ORDER ONE-SIDED
      fdiffa = functn(ic,jc,kstol) - functn(ic,jc,kstom1)
      fdiffb = functn(ic,jc,kstol) - functn(ic,jc,kstom2)
      fdiffc = functn(ic,jc,kstol) - functn(ic,jc,kstom3)
      fdiffd = functn(ic,jc,kstol) - functn(ic,jc,kstom4)
      fderiv(ic,jc,kstol) = acof1z*fdiffa + bcof1z*fdiffb  &
          + ccof1z*fdiffc + dcof1z*fdiffd
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     SCALING
!     =======
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fderiv(ic,jc,kc) = fderiv(ic,jc,kc)*ovdelz
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE dfbydz
