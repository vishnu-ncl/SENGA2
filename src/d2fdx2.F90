SUBROUTINE d2fdx2(functn,fderiv)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-26  Time: 11:25:08

!     *************************************************************************

!     D2FDX2
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     01-AUG-1996:  CREATED
!     06-APR-2003:  RSC MODIFIED FOR SENGA2

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EVALUATES SECOND X-DERIVATIVE OF SPECIFIED FUNCTION
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

real(kind=8),INTENT(IN)             :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
real(kind=8),INTENT(OUT)            :: fderiv(nxsize,nysize,nzsize)




!     LOCAL DATA
!     ==========
real(kind=8):: fdifap,fdifbp,fdifcp,fdifdp,fdifep
real(kind=8):: fdifam,fdifbm,fdifcm,fdifdm,fdifem
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
      
      fdifap = functn(icp1,jc,kc) - functn(iccc,jc,kc)
      fdifam = functn(iccc,jc,kc) - functn(icm1,jc,kc)
      fdifbp = functn(icp2,jc,kc) - functn(iccc,jc,kc)
      fdifbm = functn(iccc,jc,kc) - functn(icm2,jc,kc)
      fdifcp = functn(icp3,jc,kc) - functn(iccc,jc,kc)
      fdifcm = functn(iccc,jc,kc) - functn(icm3,jc,kc)
      fdifdp = functn(icp4,jc,kc) - functn(iccc,jc,kc)
      fdifdm = functn(iccc,jc,kc) - functn(icm4,jc,kc)
      fdifep = functn(icp5,jc,kc) - functn(iccc,jc,kc)
      fdifem = functn(iccc,jc,kc) - functn(icm5,jc,kc)
      
      fderiv(ic,jc,kc) = acofsx*(fdifap-fdifam) + bcofsx*(fdifbp-fdifbm)  &
          + ccofsx*(fdifcp-fdifcm) + dcofsx*(fdifdp-fdifdm)  &
          + ecofsx*(fdifep-fdifem)
      
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
      fdifap = functn(istap1,jc,kc) - functn(istal,jc,kc)
      fdifbp = functn(istap2,jc,kc) - functn(istal,jc,kc)
      fdifcp = functn(istap3,jc,kc) - functn(istal,jc,kc)
      fdifdp = functn(istap4,jc,kc) - functn(istal,jc,kc)
      fdifep = functn(istap5,jc,kc) - functn(istal,jc,kc)
      fderiv(istal,jc,kc) = acfs1x*fdifap + bcfs1x*fdifbp  &
          + ccfs1x*fdifcp + dcfs1x*fdifdp  &
          + ecfs1x*fdifep
      
!           LH POINT PLUS 1: 4TH ORDER MIXED
      fdifap = functn(istal,jc,kc)  - functn(istap1,jc,kc)
      fdifbp = functn(istap2,jc,kc) - functn(istap1,jc,kc)
      fdifcp = functn(istap3,jc,kc) - functn(istap1,jc,kc)
      fdifdp = functn(istap4,jc,kc) - functn(istap1,jc,kc)
      fdifep = functn(istap5,jc,kc) - functn(istap1,jc,kc)
      fderiv(istap1,jc,kc) = acfs2x*fdifap + bcfs2x*fdifbp  &
          + ccfs2x*fdifcp + dcfs2x*fdifdp  &
          + ecfs2x*fdifep
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdifap = functn(istap3,jc,kc) - functn(istap2,jc,kc)
      fdifam = functn(istap2,jc,kc) - functn(istap1,jc,kc)
      fdifbp = functn(istap4,jc,kc) - functn(istap2,jc,kc)
      fdifbm = functn(istap2,jc,kc) - functn(istal,jc,kc)
      fderiv(istap2,jc,kc) = acfs3x*(fdifap-fdifam) + bcfs3x*(fdifbp-fdifbm)
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdifap = functn(istap4,jc,kc) - functn(istap3,jc,kc)
      fdifam = functn(istap3,jc,kc) - functn(istap2,jc,kc)
      fdifbp = functn(istap5,jc,kc) - functn(istap3,jc,kc)
      fdifbm = functn(istap3,jc,kc) - functn(istap1,jc,kc)
      fdifcp = functn(istap6,jc,kc) - functn(istap3,jc,kc)
      fdifcm = functn(istap3,jc,kc) - functn(istal,jc,kc)
      fderiv(istap3,jc,kc) = acfs4x*(fdifap-fdifam)  &
          + bcfs4x*(fdifbp-fdifbm) + ccfs4x*(fdifcp-fdifcm)
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdifap = functn(istap5,jc,kc) - functn(istap4,jc,kc)
      fdifam = functn(istap4,jc,kc) - functn(istap3,jc,kc)
      fdifbp = functn(istap6,jc,kc) - functn(istap4,jc,kc)
      fdifbm = functn(istap4,jc,kc) - functn(istap2,jc,kc)
      fdifcp = functn(istap7,jc,kc) - functn(istap4,jc,kc)
      fdifcm = functn(istap4,jc,kc) - functn(istap1,jc,kc)
      fdifdp = functn(istap8,jc,kc) - functn(istap4,jc,kc)
      fdifdm = functn(istap4,jc,kc) - functn(istal,jc,kc)
      fderiv(istap4,jc,kc) = acfs5x*(fdifap-fdifam)  &
          + bcfs5x*(fdifbp-fdifbm) + ccfs5x*(fdifcp-fdifcm)  &
          + dcfs5x*(fdifdp-fdifdm)
      
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
      fdifap = functn(istom3,jc,kc) - functn(istom4,jc,kc)
      fdifam = functn(istom4,jc,kc) - functn(istom5,jc,kc)
      fdifbp = functn(istom2,jc,kc) - functn(istom4,jc,kc)
      fdifbm = functn(istom4,jc,kc) - functn(istom6,jc,kc)
      fdifcp = functn(istom1,jc,kc) - functn(istom4,jc,kc)
      fdifcm = functn(istom4,jc,kc) - functn(istom7,jc,kc)
      fdifdp = functn(istol,jc,kc)  - functn(istom4,jc,kc)
      fdifdm = functn(istom4,jc,kc) - functn(istom8,jc,kc)
      fderiv(istom4,jc,kc) = acfs5x*(fdifap-fdifam)  &
          + bcfs5x*(fdifbp-fdifbm) + ccfs5x*(fdifcp-fdifcm)  &
          + dcfs5x*(fdifdp-fdifdm)
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdifap = functn(istom2,jc,kc) - functn(istom3,jc,kc)
      fdifam = functn(istom3,jc,kc) - functn(istom4,jc,kc)
      fdifbp = functn(istom1,jc,kc) - functn(istom3,jc,kc)
      fdifbm = functn(istom3,jc,kc) - functn(istom5,jc,kc)
      fdifcp = functn(istol,jc,kc)  - functn(istom3,jc,kc)
      fdifcm = functn(istom3,jc,kc) - functn(istom6,jc,kc)
      fderiv(istom3,jc,kc) = acfs4x*(fdifap-fdifam)  &
          + bcfs4x*(fdifbp-fdifbm) + ccfs4x*(fdifcp-fdifcm)
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdifap = functn(istom1,jc,kc) - functn(istom2,jc,kc)
      fdifam = functn(istom2,jc,kc) - functn(istom3,jc,kc)
      fdifbp = functn(istol,jc,kc)  - functn(istom2,jc,kc)
      fdifbm = functn(istom2,jc,kc) - functn(istom4,jc,kc)
      fderiv(istom2,jc,kc) = acfs3x*(fdifap-fdifam) + bcfs3x*(fdifbp-fdifbm)
      
!           RH POINT MINUS 1: 4TH ORDER MIXED
      fdifap = functn(istol,jc,kc)  - functn(istom1,jc,kc)
      fdifbp = functn(istom2,jc,kc) - functn(istom1,jc,kc)
      fdifcp = functn(istom3,jc,kc) - functn(istom1,jc,kc)
      fdifdp = functn(istom4,jc,kc) - functn(istom1,jc,kc)
      fdifep = functn(istom5,jc,kc) - functn(istom1,jc,kc)
      fderiv(istom1,jc,kc) = acfs2x*fdifap + bcfs2x*fdifbp  &
          + ccfs2x*fdifcp + dcfs2x*fdifdp  &
          + ecfs2x*fdifep
      
!           RH POINT: 4TH ORDER ONE-SIDED
      fdifap = functn(istom1,jc,kc)  - functn(istol,jc,kc)
      fdifbp = functn(istom2,jc,kc)  - functn(istol,jc,kc)
      fdifcp = functn(istom3,jc,kc)  - functn(istol,jc,kc)
      fdifdp = functn(istom4,jc,kc)  - functn(istol,jc,kc)
      fdifep = functn(istom5,jc,kc)  - functn(istol,jc,kc)
      fderiv(istol,jc,kc) = acfs1x*fdifap + bcfs1x*fdifbp  &
          + ccfs1x*fdifcp + dcfs1x*fdifdp  &
          + ecfs1x*fdifep
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     SCALING
!     =======
DO kc = kstal, kstol
  DO jc = jstal, jstol
    DO ic = istal, istol
      
      fderiv(ic,jc,kc) = fderiv(ic,jc,kc)*ovdlx2
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE d2fdx2
