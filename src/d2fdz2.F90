SUBROUTINE d2fdz2(functn,fderiv)

! Code converted using TO_F90 by Alan Miller
! Date: 2022-11-09  Time: 13:29:05

!     *************************************************************************

!     D2FDZ2
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
!     EVALUATES SECOND Z-DERIVATIVE OF SPECIFIED FUNCTION
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
real(kind=dp) :: fdifap,fdifbp,fdifcp,fdifdp,fdifep
real(kind=dp) :: fdifam,fdifbm,fdifcm,fdifdm,fdifem
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
      
      fdifap = functn(ic,jc,kcp1) - functn(ic,jc,kccc)
      fdifam = functn(ic,jc,kccc) - functn(ic,jc,kcm1)
      fdifbp = functn(ic,jc,kcp2) - functn(ic,jc,kccc)
      fdifbm = functn(ic,jc,kccc) - functn(ic,jc,kcm2)
      fdifcp = functn(ic,jc,kcp3) - functn(ic,jc,kccc)
      fdifcm = functn(ic,jc,kccc) - functn(ic,jc,kcm3)
      fdifdp = functn(ic,jc,kcp4) - functn(ic,jc,kccc)
      fdifdm = functn(ic,jc,kccc) - functn(ic,jc,kcm4)
      fdifep = functn(ic,jc,kcp5) - functn(ic,jc,kccc)
      fdifem = functn(ic,jc,kccc) - functn(ic,jc,kcm5)
      
      fderiv(ic,jc,kc) = acofsz*(fdifap-fdifam) + bcofsz*(fdifbp-fdifbm)  &
          + ccofsz*(fdifcp-fdifcm) + dcofsz*(fdifdp-fdifdm)  &
          + ecofsz*(fdifep-fdifem)
      
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
      fdifap = functn(ic,jc,kstap1) - functn(ic,jc,kstal)
      fdifbp = functn(ic,jc,kstap2) - functn(ic,jc,kstal)
      fdifcp = functn(ic,jc,kstap3) - functn(ic,jc,kstal)
      fdifdp = functn(ic,jc,kstap4) - functn(ic,jc,kstal)
      fdifep = functn(ic,jc,kstap5) - functn(ic,jc,kstal)
      fderiv(ic,jc,kstal) = acfs1z*fdifap + bcfs1z*fdifbp  &
          + ccfs1z*fdifcp + dcfs1z*fdifdp  &
          + ecfs1z*fdifep
      
!           LH POINT PLUS 1: 4TH ORDER MIXED
      fdifap = functn(ic,jc,kstal)  - functn(ic,jc,kstap1)
      fdifbp = functn(ic,jc,kstap2) - functn(ic,jc,kstap1)
      fdifcp = functn(ic,jc,kstap3) - functn(ic,jc,kstap1)
      fdifdp = functn(ic,jc,kstap4) - functn(ic,jc,kstap1)
      fdifep = functn(ic,jc,kstap5) - functn(ic,jc,kstap1)
      fderiv(ic,jc,kstap1) = acfs2z*fdifap + bcfs2z*fdifbp  &
          + ccfs2z*fdifcp + dcfs2z*fdifdp  &
          + ecfs2z*fdifep
      
!           LH POINT PLUS 2: 4TH ORDER CENTRED
      fdifap = functn(ic,jc,kstap3) - functn(ic,jc,kstap2)
      fdifam = functn(ic,jc,kstap2) - functn(ic,jc,kstap1)
      fdifbp = functn(ic,jc,kstap4) - functn(ic,jc,kstap2)
      fdifbm = functn(ic,jc,kstap2) - functn(ic,jc,kstal)
      fderiv(ic,jc,kstap2) = acfs3z*(fdifap-fdifam) + bcfs3z*(fdifbp-fdifbm)
      
!           LH POINT PLUS 3: 6TH ORDER CENTRED
      fdifap = functn(ic,jc,kstap4) - functn(ic,jc,kstap3)
      fdifam = functn(ic,jc,kstap3) - functn(ic,jc,kstap2)
      fdifbp = functn(ic,jc,kstap5) - functn(ic,jc,kstap3)
      fdifbm = functn(ic,jc,kstap3) - functn(ic,jc,kstap1)
      fdifcp = functn(ic,jc,kstap6) - functn(ic,jc,kstap3)
      fdifcm = functn(ic,jc,kstap3) - functn(ic,jc,kstal)
      fderiv(ic,jc,kstap3) = acfs4z*(fdifap-fdifam)  &
          + bcfs4z*(fdifbp-fdifbm) + ccfs4z*(fdifcp-fdifcm)
      
!           LH POINT PLUS 4: 8TH ORDER CENTRED
      fdifap = functn(ic,jc,kstap5) - functn(ic,jc,kstap4)
      fdifam = functn(ic,jc,kstap4) - functn(ic,jc,kstap3)
      fdifbp = functn(ic,jc,kstap6) - functn(ic,jc,kstap4)
      fdifbm = functn(ic,jc,kstap4) - functn(ic,jc,kstap2)
      fdifcp = functn(ic,jc,kstap7) - functn(ic,jc,kstap4)
      fdifcm = functn(ic,jc,kstap4) - functn(ic,jc,kstap1)
      fdifdp = functn(ic,jc,kstap8) - functn(ic,jc,kstap4)
      fdifdm = functn(ic,jc,kstap4) - functn(ic,jc,kstal)
      fderiv(ic,jc,kstap4) = acfs5z*(fdifap-fdifam)  &
          + bcfs5z*(fdifbp-fdifbm) + ccfs5z*(fdifcp-fdifcm)  &
          + dcfs5z*(fdifdp-fdifdm)
      
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
      fdifap = functn(ic,jc,kstom3) - functn(ic,jc,kstom4)
      fdifam = functn(ic,jc,kstom4) - functn(ic,jc,kstom5)
      fdifbp = functn(ic,jc,kstom2) - functn(ic,jc,kstom4)
      fdifbm = functn(ic,jc,kstom4) - functn(ic,jc,kstom6)
      fdifcp = functn(ic,jc,kstom1) - functn(ic,jc,kstom4)
      fdifcm = functn(ic,jc,kstom4) - functn(ic,jc,kstom7)
      fdifdp = functn(ic,jc,kstol)  - functn(ic,jc,kstom4)
      fdifdm = functn(ic,jc,kstom4) - functn(ic,jc,kstom8)
      fderiv(ic,jc,kstom4) = acfs5z*(fdifap-fdifam)  &
          + bcfs5z*(fdifbp-fdifbm) + ccfs5z*(fdifcp-fdifcm)  &
          + dcfs5z*(fdifdp-fdifdm)
      
!           RH POINT MINUS 3: 6TH ORDER CENTRED
      fdifap = functn(ic,jc,kstom2) - functn(ic,jc,kstom3)
      fdifam = functn(ic,jc,kstom3) - functn(ic,jc,kstom4)
      fdifbp = functn(ic,jc,kstom1) - functn(ic,jc,kstom3)
      fdifbm = functn(ic,jc,kstom3) - functn(ic,jc,kstom5)
      fdifcp = functn(ic,jc,kstol)  - functn(ic,jc,kstom3)
      fdifcm = functn(ic,jc,kstom3) - functn(ic,jc,kstom6)
      fderiv(ic,jc,kstom3) = acfs4z*(fdifap-fdifam)  &
          + bcfs4z*(fdifbp-fdifbm) + ccfs4z*(fdifcp-fdifcm)
      
!           RH POINT MINUS 2: 4TH ORDER CENTRED
      fdifap = functn(ic,jc,kstom1) - functn(ic,jc,kstom2)
      fdifam = functn(ic,jc,kstom2) - functn(ic,jc,kstom3)
      fdifbp = functn(ic,jc,kstol)  - functn(ic,jc,kstom2)
      fdifbm = functn(ic,jc,kstom2) - functn(ic,jc,kstom4)
      fderiv(ic,jc,kstom2) = acfs3z*(fdifap-fdifam) + bcfs3z*(fdifbp-fdifbm)
      
!           RH POINT MINUS 1: 4TH ORDER MIXED
      fdifap = functn(ic,jc,kstol)  - functn(ic,jc,kstom1)
      fdifbp = functn(ic,jc,kstom2) - functn(ic,jc,kstom1)
      fdifcp = functn(ic,jc,kstom3) - functn(ic,jc,kstom1)
      fdifdp = functn(ic,jc,kstom4) - functn(ic,jc,kstom1)
      fdifep = functn(ic,jc,kstom5) - functn(ic,jc,kstom1)
      fderiv(ic,jc,kstom1) = acfs2z*fdifap + bcfs2z*fdifbp  &
          + ccfs2z*fdifcp + dcfs2z*fdifdp  &
          + ecfs2z*fdifep
      
!           RH POINT: 4TH ORDER ONE-SIDED
      fdifap = functn(ic,jc,kstom1)  - functn(ic,jc,kstol)
      fdifbp = functn(ic,jc,kstom2)  - functn(ic,jc,kstol)
      fdifcp = functn(ic,jc,kstom3)  - functn(ic,jc,kstol)
      fdifdp = functn(ic,jc,kstom4)  - functn(ic,jc,kstol)
      fdifep = functn(ic,jc,kstom5)  - functn(ic,jc,kstol)
      fderiv(ic,jc,kstol) = acfs1z*fdifap + bcfs1z*fdifbp  &
          + ccfs1z*fdifcp + dcfs1z*fdifdp  &
          + ecfs1z*fdifep
      
    END DO
  END DO
  
END IF

!     =========================================================================

!     SCALING
!     =======
DO kc = kstal, kstol
  DO jc = jstal, jstol
    DO ic = istal, istol
      
      fderiv(ic,jc,kc) = fderiv(ic,jc,kc)*ovdlz2
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE d2fdz2
