SUBROUTINE d2fdy2(functn,fderiv)

! Code converted using TO_F90 by Alan Miller
! Date: 2022-11-09  Time: 13:16:49

!     *************************************************************************

!     D2FDY2
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
!     EVALUATES SECOND Y-DERIVATIVE OF SPECIFIED FUNCTION
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

REAL(kind=8), INTENT(IN)             :: functn(nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
REAL(kind=8), INTENT(OUT)            :: fderiv(nxsize,nysize,nzsize)




!     LOCAL DATA
!     ==========
REAL(kind=8) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep
REAL(kind=8) :: fdifam,fdifbm,fdifcm,fdifdm,fdifem
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
      
      fdifap = functn(ic,jcp1,kc) - functn(ic,jccc,kc)
      fdifam = functn(ic,jccc,kc) - functn(ic,jcm1,kc)
      fdifbp = functn(ic,jcp2,kc) - functn(ic,jccc,kc)
      fdifbm = functn(ic,jccc,kc) - functn(ic,jcm2,kc)
      fdifcp = functn(ic,jcp3,kc) - functn(ic,jccc,kc)
      fdifcm = functn(ic,jccc,kc) - functn(ic,jcm3,kc)
      fdifdp = functn(ic,jcp4,kc) - functn(ic,jccc,kc)
      fdifdm = functn(ic,jccc,kc) - functn(ic,jcm4,kc)
      fdifep = functn(ic,jcp5,kc) - functn(ic,jccc,kc)
      fdifem = functn(ic,jccc,kc) - functn(ic,jcm5,kc)
      
      fderiv(ic,jc,kc) = acofsy*(fdifap-fdifam) + bcofsy*(fdifbp-fdifbm)  &
                       + ccofsy*(fdifcp-fdifcm) + dcofsy*(fdifdp-fdifdm)  &
                       + ecofsy*(fdifep-fdifem)
      
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
      fdifap = functn(ic,jstap1,kc) - functn(ic,jstal,kc)
      fdifbp = functn(ic,jstap2,kc) - functn(ic,jstal,kc)
      fdifcp = functn(ic,jstap3,kc) - functn(ic,jstal,kc)
      fdifdp = functn(ic,jstap4,kc) - functn(ic,jstal,kc)
      fdifep = functn(ic,jstap5,kc) - functn(ic,jstal,kc)
      fderiv(ic,jstal,kc) = acfs1y*fdifap + bcfs1y*fdifbp  &
                          + ccfs1y*fdifcp + dcfs1y*fdifdp  &
                          + ecfs1y*fdifep
      
!           LH POINT PLUS 1: 4TH ORDER MIXED
      fdifap = functn(ic,jstal,kc)  - functn(ic,jstap1,kc)
      fdifbp = functn(ic,jstap2,kc) - functn(ic,jstap1,kc)
      fdifcp = functn(ic,jstap3,kc) - functn(ic,jstap1,kc)
      fdifdp = functn(ic,jstap4,kc) - functn(ic,jstap1,kc)
      fdifep = functn(ic,jstap5,kc) - functn(ic,jstap1,kc)
      fderiv(ic,jstap1,kc) = acfs2y*fdifap + bcfs2y*fdifbp  &
                           + ccfs2y*fdifcp + dcfs2y*fdifdp  &
                           + ecfs2y*fdifep
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdifap = functn(ic,jstap3,kc) - functn(ic,jstap2,kc)
      fdifam = functn(ic,jstap2,kc) - functn(ic,jstap1,kc)
      fdifbp = functn(ic,jstap4,kc) - functn(ic,jstap2,kc)
      fdifbm = functn(ic,jstap2,kc) - functn(ic,jstal,kc)
      fderiv(ic,jstap2,kc) = acfs3y*(fdifap-fdifam) + bcfs3y*(fdifbp-fdifbm)
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdifap = functn(ic,jstap4,kc) - functn(ic,jstap3,kc)
      fdifam = functn(ic,jstap3,kc) - functn(ic,jstap2,kc)
      fdifbp = functn(ic,jstap5,kc) - functn(ic,jstap3,kc)
      fdifbm = functn(ic,jstap3,kc) - functn(ic,jstap1,kc)
      fdifcp = functn(ic,jstap6,kc) - functn(ic,jstap3,kc)
      fdifcm = functn(ic,jstap3,kc) - functn(ic,jstal,kc)
      fderiv(ic,jstap3,kc) = acfs4y*(fdifap-fdifam)  &
          + bcfs4y*(fdifbp-fdifbm) + ccfs4y*(fdifcp-fdifcm)
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdifap = functn(ic,jstap5,kc) - functn(ic,jstap4,kc)
      fdifam = functn(ic,jstap4,kc) - functn(ic,jstap3,kc)
      fdifbp = functn(ic,jstap6,kc) - functn(ic,jstap4,kc)
      fdifbm = functn(ic,jstap4,kc) - functn(ic,jstap2,kc)
      fdifcp = functn(ic,jstap7,kc) - functn(ic,jstap4,kc)
      fdifcm = functn(ic,jstap4,kc) - functn(ic,jstap1,kc)
      fdifdp = functn(ic,jstap8,kc) - functn(ic,jstap4,kc)
      fdifdm = functn(ic,jstap4,kc) - functn(ic,jstal,kc)
      fderiv(ic,jstap4,kc) = acfs5y*(fdifap-fdifam)  &
          + bcfs5y*(fdifbp-fdifbm) + ccfs5y*(fdifcp-fdifcm)  &
          + dcfs5y*(fdifdp-fdifdm)
      
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
      fdifap = functn(ic,jstom3,kc) - functn(ic,jstom4,kc)
      fdifam = functn(ic,jstom4,kc) - functn(ic,jstom5,kc)
      fdifbp = functn(ic,jstom2,kc) - functn(ic,jstom4,kc)
      fdifbm = functn(ic,jstom4,kc) - functn(ic,jstom6,kc)
      fdifcp = functn(ic,jstom1,kc) - functn(ic,jstom4,kc)
      fdifcm = functn(ic,jstom4,kc) - functn(ic,jstom7,kc)
      fdifdp = functn(ic,jstol,kc)  - functn(ic,jstom4,kc)
      fdifdm = functn(ic,jstom4,kc) - functn(ic,jstom8,kc)
      fderiv(ic,jstom4,kc) = acfs5y*(fdifap-fdifam)  &
          + bcfs5y*(fdifbp-fdifbm) + ccfs5y*(fdifcp-fdifcm)  &
          + dcfs5y*(fdifdp-fdifdm)
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdifap = functn(ic,jstom2,kc) - functn(ic,jstom3,kc)
      fdifam = functn(ic,jstom3,kc) - functn(ic,jstom4,kc)
      fdifbp = functn(ic,jstom1,kc) - functn(ic,jstom3,kc)
      fdifbm = functn(ic,jstom3,kc) - functn(ic,jstom5,kc)
      fdifcp = functn(ic,jstol,kc)  - functn(ic,jstom3,kc)
      fdifcm = functn(ic,jstom3,kc) - functn(ic,jstom6,kc)
      fderiv(ic,jstom3,kc) = acfs4y*(fdifap-fdifam)  &
          + bcfs4y*(fdifbp-fdifbm) + ccfs4y*(fdifcp-fdifcm)
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdifap = functn(ic,jstom1,kc) - functn(ic,jstom2,kc)
      fdifam = functn(ic,jstom2,kc) - functn(ic,jstom3,kc)
      fdifbp = functn(ic,jstol,kc)  - functn(ic,jstom2,kc)
      fdifbm = functn(ic,jstom2,kc) - functn(ic,jstom4,kc)
      fderiv(ic,jstom2,kc) = acfs3y*(fdifap-fdifam) + bcfs3y*(fdifbp-fdifbm)
      
!           RH POINT MINUS 1: 4TH ORDER MIXED
      fdifap = functn(ic,jstol,kc)  - functn(ic,jstom1,kc)
      fdifbp = functn(ic,jstom2,kc) - functn(ic,jstom1,kc)
      fdifcp = functn(ic,jstom3,kc) - functn(ic,jstom1,kc)
      fdifdp = functn(ic,jstom4,kc) - functn(ic,jstom1,kc)
      fdifep = functn(ic,jstom5,kc) - functn(ic,jstom1,kc)
      fderiv(ic,jstom1,kc) = acfs2y*fdifap + bcfs2y*fdifbp  &
                           + ccfs2y*fdifcp + dcfs2y*fdifdp  &
                           + ecfs2y*fdifep
      
!           RH POINT: 4TH ORDER ONE-SIDED
      fdifap = functn(ic,jstom1,kc) - functn(ic,jstol,kc)
      fdifbp = functn(ic,jstom2,kc) - functn(ic,jstol,kc)
      fdifcp = functn(ic,jstom3,kc) - functn(ic,jstol,kc)
      fdifdp = functn(ic,jstom4,kc) - functn(ic,jstol,kc)
      fdifep = functn(ic,jstom5,kc) - functn(ic,jstol,kc)
      fderiv(ic,jstol,kc) = acfs1y*fdifap + bcfs1y*fdifbp  &
                          + ccfs1y*fdifcp + dcfs1y*fdifdp  &
                          + ecfs1y*fdifep
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     SCALING
!     =======
DO kc = kstal, kstol
  DO jc = jstal, jstol
    DO ic = istal, istol
      
      fderiv(ic,jc,kc) = fderiv(ic,jc,kc)*ovdly2
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE d2fdy2
