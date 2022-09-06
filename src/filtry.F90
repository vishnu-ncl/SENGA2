SUBROUTINE filtry(functn)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 22:58:58

!     *************************************************************************

!     FILTRY
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
!     Y DIRECTION

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
INTEGER :: jstart,jfinis
INTEGER :: jcm6,jcm5,jcm4,jcm3,jcm2,jcm1,jccc
INTEGER :: jcp1,jcp2,jcp3,jcp4,jcp5,jcp6


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

jstart = jstal
jfinis = jstol
IF(nendyl == nbound)jstart = jstap6
IF(nendyr == nbound)jfinis = jstom6

!     =========================================================================

!     INTERIOR SCHEME
!     ===============

!     TWELFTH ORDER EXPLICIT FILTER
DO kc = kstal,kstol
  
  jcm5 = jstart-6
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
  jcp6 = jstart+5
  
  DO jc = jstart,jfinis
    
    jcm6 = jcm5
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
    jcp5 = jcp6
    jcp6 = jc+6
    
    DO ic = istal,istol
      
      fdiffa = functn(ic,jcp1,kc) + functn(ic,jcm1,kc)
      fdiffb = functn(ic,jcp2,kc) + functn(ic,jcm2,kc)
      fdiffc = functn(ic,jcp3,kc) + functn(ic,jcm3,kc)
      fdiffd = functn(ic,jcp4,kc) + functn(ic,jcm4,kc)
      fdiffe = functn(ic,jcp5,kc) + functn(ic,jcm5,kc)
      fdifff = functn(ic,jcp6,kc) + functn(ic,jcm6,kc)
      
      filter(ic,jc,kc) = facofy*fdiffa + fbcofy*fdiffb  &
          + fccofy*fdiffc + fdcofy*fdiffd  &
          + fecofy*fdiffe + ffcofy*fdifff  &
          + fgcofy*functn(ic,jc,kc)
      
    END DO
    
  END DO
END DO

!     =========================================================================

!     LH END
!     ======
IF(nendyl == nbound)THEN
  
!       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    DO ic = istal,istol
      
!           LH POINT: 6TH ORDER ONE-SIDED
      filter(ic,jstal,kc) = facf1y*functn(ic,jstal,kc)  &
          + fbcf1y*functn(ic,jstap1,kc) + fccf1y*functn(ic,jstap2,kc)  &
          + fdcf1y*functn(ic,jstap3,kc) + fecf1y*functn(ic,jstap4,kc)  &
          + ffcf1y*functn(ic,jstap5,kc) + fgcf1y*functn(ic,jstap6,kc)
      
!           LH POINT PLUS 1: 7TH ORDER MIXED
      filter(ic,jstap1,kc) = facf2y*functn(ic,jstal,kc)  &
          + fbcf2y*functn(ic,jstap1,kc) + fccf2y*functn(ic,jstap2,kc)  &
          + fdcf2y*functn(ic,jstap3,kc) + fecf2y*functn(ic,jstap4,kc)  &
          + ffcf2y*functn(ic,jstap5,kc) + fgcf2y*functn(ic,jstap6,kc)  &
          + fhcf2y*functn(ic,jstap7,kc)
      
!           LH POINT PLUS 2: 8TH ORDER MIXED
      filter(ic,jstap2,kc) = facf3y*functn(ic,jstal,kc)  &
          + fbcf3y*functn(ic,jstap1,kc) + fccf3y*functn(ic,jstap2,kc)  &
          + fdcf3y*functn(ic,jstap3,kc) + fecf3y*functn(ic,jstap4,kc)  &
          + ffcf3y*functn(ic,jstap5,kc) + fgcf3y*functn(ic,jstap6,kc)  &
          + fhcf3y*functn(ic,jstap7,kc) + ficf3y*functn(ic,jstap8,kc)
      
!           LH POINT PLUS 3: 9TH ORDER MIXED
      filter(ic,jstap3,kc) = facf4y*functn(ic,jstal,kc)  &
          + fbcf4y*functn(ic,jstap1,kc) + fccf4y*functn(ic,jstap2,kc)  &
          + fdcf4y*functn(ic,jstap3,kc) + fecf4y*functn(ic,jstap4,kc)  &
          + ffcf4y*functn(ic,jstap5,kc) + fgcf4y*functn(ic,jstap6,kc)  &
          + fhcf4y*functn(ic,jstap7,kc) + ficf4y*functn(ic,jstap8,kc)  &
          + fjcf4y*functn(ic,jstap9,kc)
      
!           LH POINT PLUS 4: 10TH ORDER MIXED
      filter(ic,jstap4,kc) = facf5y*functn(ic,jstal,kc)  &
          + fbcf5y*functn(ic,jstap1,kc) + fccf5y*functn(ic,jstap2,kc)  &
          + fdcf5y*functn(ic,jstap3,kc) + fecf5y*functn(ic,jstap4,kc)  &
          + ffcf5y*functn(ic,jstap5,kc) + fgcf5y*functn(ic,jstap6,kc)  &
          + fhcf5y*functn(ic,jstap7,kc) + ficf5y*functn(ic,jstap8,kc)  &
          + fjcf5y*functn(ic,jstap9,kc) + fkcf5y*functn(ic,jstapa,kc)
      
!           LH POINT PLUS 5: 11TH ORDER MIXED
      filter(ic,jstap5,kc) = facf6y*functn(ic,jstal,kc)  &
          + fbcf6y*functn(ic,jstap1,kc) + fccf6y*functn(ic,jstap2,kc)  &
          + fdcf6y*functn(ic,jstap3,kc) + fecf6y*functn(ic,jstap4,kc)  &
          + ffcf6y*functn(ic,jstap5,kc) + fgcf6y*functn(ic,jstap6,kc)  &
          + fhcf6y*functn(ic,jstap7,kc) + ficf6y*functn(ic,jstap8,kc)  &
          + fjcf6y*functn(ic,jstap9,kc) + fkcf6y*functn(ic,jstapa,kc)  &
          + flcf6y*functn(ic,jstapb,kc)
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     RH END
!     ======
IF(nendyr == nbound)THEN
  
!       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    DO ic = istal,istol
      
!           RH POINT MINUS 5: 11TH ORDER MIXED
      filter(ic,jstom5,kc) = facf6y*functn(ic,jstol,kc)  &
          + fbcf6y*functn(ic,jstom1,kc) + fccf6y*functn(ic,jstom2,kc)  &
          + fdcf6y*functn(ic,jstom3,kc) + fecf6y*functn(ic,jstom4,kc)  &
          + ffcf6y*functn(ic,jstom5,kc) + fgcf6y*functn(ic,jstom6,kc)  &
          + fhcf6y*functn(ic,jstom7,kc) + ficf6y*functn(ic,jstom8,kc)  &
          + fjcf6y*functn(ic,jstom9,kc) + fkcf6y*functn(ic,jstoma,kc)  &
          + flcf6y*functn(ic,jstomb,kc)
      
!           RH POINT MINUS 4: 10TH ORDER MIXED
      filter(ic,jstom4,kc) = facf5y*functn(ic,jstol,kc)  &
          + fbcf5y*functn(ic,jstom1,kc) + fccf5y*functn(ic,jstom2,kc)  &
          + fdcf5y*functn(ic,jstom3,kc) + fecf5y*functn(ic,jstom4,kc)  &
          + ffcf5y*functn(ic,jstom5,kc) + fgcf5y*functn(ic,jstom6,kc)  &
          + fhcf5y*functn(ic,jstom7,kc) + ficf5y*functn(ic,jstom8,kc)  &
          + fjcf5y*functn(ic,jstom9,kc) + fkcf5y*functn(ic,jstoma,kc)
      
!           RH POINT MINUS 3: 9TH ORDER MIXED
      filter(ic,jstom3,kc) = facf4y*functn(ic,jstol,kc)  &
          + fbcf4y*functn(ic,jstom1,kc) + fccf4y*functn(ic,jstom2,kc)  &
          + fdcf4y*functn(ic,jstom3,kc) + fecf4y*functn(ic,jstom4,kc)  &
          + ffcf4y*functn(ic,jstom5,kc) + fgcf4y*functn(ic,jstom6,kc)  &
          + fhcf4y*functn(ic,jstom7,kc) + ficf4y*functn(ic,jstom8,kc)  &
          + fjcf4y*functn(ic,jstom9,kc)
      
!           RH POINT MINUS 2: 8TH ORDER MIXED
      filter(ic,jstom2,kc) = facf3y*functn(ic,jstol,kc)  &
          + fbcf3y*functn(ic,jstom1,kc) + fccf3y*functn(ic,jstom2,kc)  &
          + fdcf3y*functn(ic,jstom3,kc) + fecf3y*functn(ic,jstom4,kc)  &
          + ffcf3y*functn(ic,jstom5,kc) + fgcf3y*functn(ic,jstom6,kc)  &
          + fhcf3y*functn(ic,jstom7,kc) + ficf3y*functn(ic,jstom8,kc)
      
!           RH POINT PLUS 1: 7TH ORDER MIXED
      filter(ic,jstom1,kc) = facf2y*functn(ic,jstol,kc)  &
          + fbcf2y*functn(ic,jstom1,kc) + fccf2y*functn(ic,jstom2,kc)  &
          + fdcf2y*functn(ic,jstom3,kc) + fecf2y*functn(ic,jstom4,kc)  &
          + ffcf2y*functn(ic,jstom5,kc) + fgcf2y*functn(ic,jstom6,kc)  &
          + fhcf2y*functn(ic,jstom7,kc)
      
!           RH POINT: 6TH ORDER ONE-SIDED
      filter(ic,jstol,kc) = facf1y*functn(ic,jstol,kc)  &
          + fbcf1y*functn(ic,jstom1,kc) + fccf1y*functn(ic,jstom2,kc)  &
          + fdcf1y*functn(ic,jstom3,kc) + fecf1y*functn(ic,jstom4,kc)  &
          + ffcf1y*functn(ic,jstom5,kc) + fgcf1y*functn(ic,jstom6,kc)
      
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
END SUBROUTINE filtry
