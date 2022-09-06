SUBROUTINE filtrz(functn)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 22:59:01

!     *************************************************************************

!     FILTRZ
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     30-AUG-2009:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EXPLICIT 12TH ORDER FINITE DIFFERENCE FILTER
!     WITH EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER END CONDITIONS
!     Z DIRECTION

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

REAL(KIND=8), INTENT(IN OUT)         :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)



!     LOCAL DATA
!     ==========
REAL(KIND=8) :: filter(nxsize,nysize,nzsize)
REAL(KIND=8) :: fdiffa,fdiffb,fdiffc,fdiffd,fdiffe,fdifff
INTEGER :: ic,jc,kc
INTEGER :: kstart,kfinis
INTEGER :: kcm6,kcm5,kcm4,kcm3,kcm2,kcm1,kccc
INTEGER :: kcp1,kcp2,kcp3,kcp4,kcp5,kcp6


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

kstart = kstal
kfinis = kstol
IF(nendzl == nbound)kstart = kstap6
IF(nendzr == nbound)kfinis = kstom6

!     =========================================================================

!     INTERIOR SCHEME
!     ===============

!     TWELFTH ORDER EXPLICIT FILTER
kcm5 = kstart-6
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
kcp6 = kstart+5

DO kc = kstart,kfinis
  
  kcm6 = kcm5
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
  kcp5 = kcp6
  kcp6 = kc+6
  
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      fdiffa = functn(ic,jc,kcp1) + functn(ic,jc,kcm1)
      fdiffb = functn(ic,jc,kcp2) + functn(ic,jc,kcm2)
      fdiffc = functn(ic,jc,kcp3) + functn(ic,jc,kcm3)
      fdiffd = functn(ic,jc,kcp4) + functn(ic,jc,kcm4)
      fdiffe = functn(ic,jc,kcp5) + functn(ic,jc,kcm5)
      fdifff = functn(ic,jc,kcp6) + functn(ic,jc,kcm6)
      
      filter(ic,jc,kc) = facofz*fdiffa + fbcofz*fdiffb  &
          + fccofz*fdiffc + fdcofz*fdiffd  &
          + fecofz*fdiffe + ffcofz*fdifff  &
          + fgcofz*functn(ic,jc,kc)
      
    END DO
    
  END DO
END DO

!     =========================================================================

!     LH END
!     ======
IF(nendzl == nbound)THEN
  
!       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
  DO jc = jstal,jstol
    DO ic = istal,istol
      
!           LH POINT: 6TH ORDER ONE-SIDED
      filter(ic,jc,kstal) = facf1z*functn(ic,jc,kstal)  &
          + fbcf1z*functn(ic,jc,kstap1) + fccf1z*functn(ic,jc,kstap2)  &
          + fdcf1z*functn(ic,jc,kstap3) + fecf1z*functn(ic,jc,kstap4)  &
          + ffcf1z*functn(ic,jc,kstap5) + fgcf1z*functn(ic,jc,kstap6)
      
!           LH POINT PLUS 1: 7TH ORDER MIXED
      filter(ic,jc,kstap1) = facf2z*functn(ic,jc,kstal)  &
          + fbcf2z*functn(ic,jc,kstap1) + fccf2z*functn(ic,jc,kstap2)  &
          + fdcf2z*functn(ic,jc,kstap3) + fecf2z*functn(ic,jc,kstap4)  &
          + ffcf2z*functn(ic,jc,kstap5) + fgcf2z*functn(ic,jc,kstap6)  &
          + fhcf2z*functn(ic,jc,kstap7)
      
!           LH POINT PLUS 2: 8TH ORDER MIXED
      filter(ic,jc,kstap2) = facf3z*functn(ic,jc,kstal)  &
          + fbcf3z*functn(ic,jc,kstap1) + fccf3z*functn(ic,jc,kstap2)  &
          + fdcf3z*functn(ic,jc,kstap3) + fecf3z*functn(ic,jc,kstap4)  &
          + ffcf3z*functn(ic,jc,kstap5) + fgcf3z*functn(ic,jc,kstap6)  &
          + fhcf3z*functn(ic,jc,kstap7) + ficf3z*functn(ic,jc,kstap8)
      
!           LH POINT PLUS 3: 9TH ORDER MIXED
      filter(ic,jc,kstap3) = facf4z*functn(ic,jc,kstal)  &
          + fbcf4z*functn(ic,jc,kstap1) + fccf4z*functn(ic,jc,kstap2)  &
          + fdcf4z*functn(ic,jc,kstap3) + fecf4z*functn(ic,jc,kstap4)  &
          + ffcf4z*functn(ic,jc,kstap5) + fgcf4z*functn(ic,jc,kstap6)  &
          + fhcf4z*functn(ic,jc,kstap7) + ficf4z*functn(ic,jc,kstap8)  &
          + fjcf4z*functn(ic,jc,kstap9)
      
!           LH POINT PLUS 4: 10TH ORDER MIXED
      filter(ic,jc,kstap4) = facf5z*functn(ic,jc,kstal)  &
          + fbcf5z*functn(ic,jc,kstap1) + fccf5z*functn(ic,jc,kstap2)  &
          + fdcf5z*functn(ic,jc,kstap3) + fecf5z*functn(ic,jc,kstap4)  &
          + ffcf5z*functn(ic,jc,kstap5) + fgcf5z*functn(ic,jc,kstap6)  &
          + fhcf5z*functn(ic,jc,kstap7) + ficf5z*functn(ic,jc,kstap8)  &
          + fjcf5z*functn(ic,jc,kstap9) + fkcf5z*functn(ic,jc,kstapa)
      
!           LH POINT PLUS 5: 11TH ORDER MIXED
      filter(ic,jc,kstap5) = facf6z*functn(ic,jc,kstal)  &
          + fbcf6z*functn(ic,jc,kstap1) + fccf6z*functn(ic,jc,kstap2)  &
          + fdcf6z*functn(ic,jc,kstap3) + fecf6z*functn(ic,jc,kstap4)  &
          + ffcf6z*functn(ic,jc,kstap5) + fgcf6z*functn(ic,jc,kstap6)  &
          + fhcf6z*functn(ic,jc,kstap7) + ficf6z*functn(ic,jc,kstap8)  &
          + fjcf6z*functn(ic,jc,kstap9) + fkcf6z*functn(ic,jc,kstapa)  &
          + flcf6z*functn(ic,jc,kstapb)
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     RH END
!     ======
IF(nendzr == nbound)THEN
  
!       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
  DO jc = jstal,jstol
    DO ic = istal,istol
      
!           RH POINT MINUS 5: 11TH ORDER MIXED
      filter(ic,jc,kstom5) = facf6z*functn(ic,jc,kstol)  &
          + fbcf6z*functn(ic,jc,kstom1) + fccf6z*functn(ic,jc,kstom2)  &
          + fdcf6z*functn(ic,jc,kstom3) + fecf6z*functn(ic,jc,kstom4)  &
          + ffcf6z*functn(ic,jc,kstom5) + fgcf6z*functn(ic,jc,kstom6)  &
          + fhcf6z*functn(ic,jc,kstom7) + ficf6z*functn(ic,jc,kstom8)  &
          + fjcf6z*functn(ic,jc,kstom9) + fkcf6z*functn(ic,jc,kstoma)  &
          + flcf6z*functn(ic,jc,kstomb)
      
!           RH POINT MINUS 4: 10TH ORDER MIXED
      filter(ic,jc,kstom4) = facf5z*functn(ic,jc,kstol)  &
          + fbcf5z*functn(ic,jc,kstom1) + fccf5z*functn(ic,jc,kstom2)  &
          + fdcf5z*functn(ic,jc,kstom3) + fecf5z*functn(ic,jc,kstom4)  &
          + ffcf5z*functn(ic,jc,kstom5) + fgcf5z*functn(ic,jc,kstom6)  &
          + fhcf5z*functn(ic,jc,kstom7) + ficf5z*functn(ic,jc,kstom8)  &
          + fjcf5z*functn(ic,jc,kstom9) + fkcf5z*functn(ic,jc,kstoma)
      
!           RH POINT MINUS 3: 9TH ORDER MIXED
      filter(ic,jc,kstom3) = facf4z*functn(ic,jc,kstol)  &
          + fbcf4z*functn(ic,jc,kstom1) + fccf4z*functn(ic,jc,kstom2)  &
          + fdcf4z*functn(ic,jc,kstom3) + fecf4z*functn(ic,jc,kstom4)  &
          + ffcf4z*functn(ic,jc,kstom5) + fgcf4z*functn(ic,jc,kstom6)  &
          + fhcf4z*functn(ic,jc,kstom7) + ficf4z*functn(ic,jc,kstom8)  &
          + fjcf4z*functn(ic,jc,kstom9)
      
!           RH POINT MINUS 2: 8TH ORDER MIXED
      filter(ic,jc,kstom2) = facf3z*functn(ic,jc,kstol)  &
          + fbcf3z*functn(ic,jc,kstom1) + fccf3z*functn(ic,jc,kstom2)  &
          + fdcf3z*functn(ic,jc,kstom3) + fecf3z*functn(ic,jc,kstom4)  &
          + ffcf3z*functn(ic,jc,kstom5) + fgcf3z*functn(ic,jc,kstom6)  &
          + fhcf3z*functn(ic,jc,kstom7) + ficf3z*functn(ic,jc,kstom8)
      
!           RH POINT PLUS 1: 7TH ORDER MIXED
      filter(ic,jc,kstom1) = facf2z*functn(ic,jc,kstol)  &
          + fbcf2z*functn(ic,jc,kstom1) + fccf2z*functn(ic,jc,kstom2)  &
          + fdcf2z*functn(ic,jc,kstom3) + fecf2z*functn(ic,jc,kstom4)  &
          + ffcf2z*functn(ic,jc,kstom5) + fgcf2z*functn(ic,jc,kstom6)  &
          + fhcf2z*functn(ic,jc,kstom7)
      
!           RH POINT: 6TH ORDER ONE-SIDED
      filter(ic,jc,kstol) = facf1z*functn(ic,jc,kstol)  &
          + fbcf1z*functn(ic,jc,kstom1) + fccf1z*functn(ic,jc,kstom2)  &
          + fdcf1z*functn(ic,jc,kstom3) + fecf1z*functn(ic,jc,kstom4)  &
          + ffcf1z*functn(ic,jc,kstom5) + fgcf1z*functn(ic,jc,kstom6)
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     COPY BACK
!     =========
DO kc = kstal,kstol
  DO jc = jstal,jstol
    DO ic = istal,istol
      
      functn(ic,jc,kc) = filter(ic,jc,kc)
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE filtrz
