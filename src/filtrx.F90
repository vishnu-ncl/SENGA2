SUBROUTINE filtrx(functn)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 22:58:55

!     *************************************************************************

!     FILTRX
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
!     X DIRECTION

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
INTEGER :: istart,ifinis
INTEGER :: icm6,icm5,icm4,icm3,icm2,icm1,iccc
INTEGER :: icp1,icp2,icp3,icp4,icp5,icp6


!     BEGIN
!     =====

!     =========================================================================

!     END CONDITIONS
!     ==============

istart = istal
ifinis = istol
IF(nendxl == nbound)istart = istap6
IF(nendxr == nbound)ifinis = istom6

!     =========================================================================

!     INTERIOR SCHEME
!     ===============

!     TWELFTH ORDER EXPLICIT FILTER
DO kc = kstal,kstol
  DO jc = jstal,jstol
    
    icm5 = istart-6
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
    icp6 = istart+5
    
    DO ic = istart,ifinis
      
      icm6 = icm5
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
      icp5 = icp6
      icp6 = ic+6
      
      fdiffa = functn(icp1,jc,kc) + functn(icm1,jc,kc)
      fdiffb = functn(icp2,jc,kc) + functn(icm2,jc,kc)
      fdiffc = functn(icp3,jc,kc) + functn(icm3,jc,kc)
      fdiffd = functn(icp4,jc,kc) + functn(icm4,jc,kc)
      fdiffe = functn(icp5,jc,kc) + functn(icm5,jc,kc)
      fdifff = functn(icp6,jc,kc) + functn(icm6,jc,kc)
      
      filter(ic,jc,kc) = facofx*fdiffa + fbcofx*fdiffb  &
          + fccofx*fdiffc + fdcofx*fdiffd  &
          + fecofx*fdiffe + ffcofx*fdifff  &
          + fgcofx*functn(ic,jc,kc)
      
    END DO
    
  END DO
END DO

!     =========================================================================

!     LH END
!     ======
IF(nendxl == nbound)THEN
  
!       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
!           LH POINT: 6TH ORDER ONE-SIDED
      filter(istal,jc,kc) = facf1x*functn(istal,jc,kc)  &
          + fbcf1x*functn(istap1,jc,kc) + fccf1x*functn(istap2,jc,kc)  &
          + fdcf1x*functn(istap3,jc,kc) + fecf1x*functn(istap4,jc,kc)  &
          + ffcf1x*functn(istap5,jc,kc) + fgcf1x*functn(istap6,jc,kc)
      
!           LH POINT PLUS 1: 7TH ORDER MIXED
      filter(istap1,jc,kc) = facf2x*functn(istal,jc,kc)  &
          + fbcf2x*functn(istap1,jc,kc) + fccf2x*functn(istap2,jc,kc)  &
          + fdcf2x*functn(istap3,jc,kc) + fecf2x*functn(istap4,jc,kc)  &
          + ffcf2x*functn(istap5,jc,kc) + fgcf2x*functn(istap6,jc,kc)  &
          + fhcf2x*functn(istap7,jc,kc)
      
!           LH POINT PLUS 2: 8TH ORDER MIXED
      filter(istap2,jc,kc) = facf3x*functn(istal,jc,kc)  &
          + fbcf3x*functn(istap1,jc,kc) + fccf3x*functn(istap2,jc,kc)  &
          + fdcf3x*functn(istap3,jc,kc) + fecf3x*functn(istap4,jc,kc)  &
          + ffcf3x*functn(istap5,jc,kc) + fgcf3x*functn(istap6,jc,kc)  &
          + fhcf3x*functn(istap7,jc,kc) + ficf3x*functn(istap8,jc,kc)
      
!           LH POINT PLUS 3: 9TH ORDER MIXED
      filter(istap3,jc,kc) = facf4x*functn(istal,jc,kc)  &
          + fbcf4x*functn(istap1,jc,kc) + fccf4x*functn(istap2,jc,kc)  &
          + fdcf4x*functn(istap3,jc,kc) + fecf4x*functn(istap4,jc,kc)  &
          + ffcf4x*functn(istap5,jc,kc) + fgcf4x*functn(istap6,jc,kc)  &
          + fhcf4x*functn(istap7,jc,kc) + ficf4x*functn(istap8,jc,kc)  &
          + fjcf4x*functn(istap9,jc,kc)
      
!           LH POINT PLUS 4: 10TH ORDER MIXED
      filter(istap4,jc,kc) = facf5x*functn(istal,jc,kc)  &
          + fbcf5x*functn(istap1,jc,kc) + fccf5x*functn(istap2,jc,kc)  &
          + fdcf5x*functn(istap3,jc,kc) + fecf5x*functn(istap4,jc,kc)  &
          + ffcf5x*functn(istap5,jc,kc) + fgcf5x*functn(istap6,jc,kc)  &
          + fhcf5x*functn(istap7,jc,kc) + ficf5x*functn(istap8,jc,kc)  &
          + fjcf5x*functn(istap9,jc,kc) + fkcf5x*functn(istapa,jc,kc)
      
!           LH POINT PLUS 5: 11TH ORDER MIXED
      filter(istap5,jc,kc) = facf6x*functn(istal,jc,kc)  &
          + fbcf6x*functn(istap1,jc,kc) + fccf6x*functn(istap2,jc,kc)  &
          + fdcf6x*functn(istap3,jc,kc) + fecf6x*functn(istap4,jc,kc)  &
          + ffcf6x*functn(istap5,jc,kc) + fgcf6x*functn(istap6,jc,kc)  &
          + fhcf6x*functn(istap7,jc,kc) + ficf6x*functn(istap8,jc,kc)  &
          + fjcf6x*functn(istap9,jc,kc) + fkcf6x*functn(istapa,jc,kc)  &
          + flcf6x*functn(istapb,jc,kc)
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     RH END
!     ======
IF(nendxr == nbound)THEN
  
!       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
  DO kc = kstal,kstol
    DO jc = jstal,jstol
      
!           RH POINT MINUS 5: 11TH ORDER MIXED
      filter(istom5,jc,kc) = facf6x*functn(istol,jc,kc)  &
          + fbcf6x*functn(istom1,jc,kc) + fccf6x*functn(istom2,jc,kc)  &
          + fdcf6x*functn(istom3,jc,kc) + fecf6x*functn(istom4,jc,kc)  &
          + ffcf6x*functn(istom5,jc,kc) + fgcf6x*functn(istom6,jc,kc)  &
          + fhcf6x*functn(istom7,jc,kc) + ficf6x*functn(istom8,jc,kc)  &
          + fjcf6x*functn(istom9,jc,kc) + fkcf6x*functn(istoma,jc,kc)  &
          + flcf6x*functn(istomb,jc,kc)
      
!           RH POINT MINUS 4: 10TH ORDER MIXED
      filter(istom4,jc,kc) = facf5x*functn(istol,jc,kc)  &
          + fbcf5x*functn(istom1,jc,kc) + fccf5x*functn(istom2,jc,kc)  &
          + fdcf5x*functn(istom3,jc,kc) + fecf5x*functn(istom4,jc,kc)  &
          + ffcf5x*functn(istom5,jc,kc) + fgcf5x*functn(istom6,jc,kc)  &
          + fhcf5x*functn(istom7,jc,kc) + ficf5x*functn(istom8,jc,kc)  &
          + fjcf5x*functn(istom9,jc,kc) + fkcf5x*functn(istoma,jc,kc)
      
!           RH POINT MINUS 3: 9TH ORDER MIXED
      filter(istom3,jc,kc) = facf4x*functn(istol,jc,kc)  &
          + fbcf4x*functn(istom1,jc,kc) + fccf4x*functn(istom2,jc,kc)  &
          + fdcf4x*functn(istom3,jc,kc) + fecf4x*functn(istom4,jc,kc)  &
          + ffcf4x*functn(istom5,jc,kc) + fgcf4x*functn(istom6,jc,kc)  &
          + fhcf4x*functn(istom7,jc,kc) + ficf4x*functn(istom8,jc,kc)  &
          + fjcf4x*functn(istom9,jc,kc)
      
!           RH POINT MINUS 2: 8TH ORDER MIXED
      filter(istom2,jc,kc) = facf3x*functn(istol,jc,kc)  &
          + fbcf3x*functn(istom1,jc,kc) + fccf3x*functn(istom2,jc,kc)  &
          + fdcf3x*functn(istom3,jc,kc) + fecf3x*functn(istom4,jc,kc)  &
          + ffcf3x*functn(istom5,jc,kc) + fgcf3x*functn(istom6,jc,kc)  &
          + fhcf3x*functn(istom7,jc,kc) + ficf3x*functn(istom8,jc,kc)
      
!           RH POINT PLUS 1: 7TH ORDER MIXED
      filter(istom1,jc,kc) = facf2x*functn(istol,jc,kc)  &
          + fbcf2x*functn(istom1,jc,kc) + fccf2x*functn(istom2,jc,kc)  &
          + fdcf2x*functn(istom3,jc,kc) + fecf2x*functn(istom4,jc,kc)  &
          + ffcf2x*functn(istom5,jc,kc) + fgcf2x*functn(istom6,jc,kc)  &
          + fhcf2x*functn(istom7,jc,kc)
      
!           RH POINT: 6TH ORDER ONE-SIDED
      filter(istol,jc,kc) = facf1x*functn(istol,jc,kc)  &
          + fbcf1x*functn(istom1,jc,kc) + fccf1x*functn(istom2,jc,kc)  &
          + fdcf1x*functn(istom3,jc,kc) + fecf1x*functn(istom4,jc,kc)  &
          + ffcf1x*functn(istom5,jc,kc) + fgcf1x*functn(istom6,jc,kc)
      
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
END SUBROUTINE filtrx
