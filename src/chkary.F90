SUBROUTINE chkary(array,  &
        nchkxl,nchkxr,nchkyl,nchkyr,nchkzl,nchkzr,ncpcmx,  &
        istac,istoc,jstac,jstoc,kstac,kstoc,icpec)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:15:42

!     *************************************************************************

!     CHKARY
!     ======

!     AUTHOR
!     ------
!     R.S.CANT - CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     04-JUL-2004:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     DIAGNOSTIC ROUTINE
!     CHECKS UNIFORMITY OF THE SPECIFIED ARRAY
!     VARIANT FOR SPECIES ARRAYS

!     *************************************************************************
use data_types

!     ARGUMENTS
!     =========

REAL(KIND=dp), INTENT(IN)             :: array(nchkxl:nchkxr,nchkyl:nchkyr,nchkzl:nchkzr,ncpcmx)
INTEGER, INTENT(IN OUT)                  :: nchkxl
INTEGER, INTENT(IN OUT)                  :: nchkxr
INTEGER, INTENT(IN OUT)                  :: nchkyl
INTEGER, INTENT(IN OUT)                  :: nchkyr
INTEGER, INTENT(IN OUT)                  :: nchkzl
INTEGER, INTENT(IN OUT)                  :: nchkzr
INTEGER, INTENT(IN OUT)                  :: ncpcmx
INTEGER, INTENT(IN)                      :: istac
INTEGER, INTENT(IN)                      :: istoc
INTEGER, INTENT(IN)                      :: jstac
INTEGER, INTENT(IN)                      :: jstoc
INTEGER, INTENT(IN)                      :: kstac
INTEGER, INTENT(IN)                      :: kstoc
INTEGER, INTENT(IN OUT)                  :: icpec


!     LOCAL DATA
!     ==========
REAL(KIND=dp) :: arrmin,arrmax
INTEGER :: ic,jc,kc
INTEGER :: icmin,jcmin,kcmin
INTEGER :: icmax,jcmax,kcmax


!     BEGIN
!     =====

WRITE(6,*)nchkxl,nchkxr,nchkyl,nchkyr,nchkzl,nchkzr,ncpcmx
WRITE(6,*)istac,istoc,jstac,jstoc,kstac,kstoc,icpec

!     =========================================================================

arrmin = array(istac,jstac,kstac,icpec)
icmin = istac
jcmin = jstac
kcmin = kstac
arrmax = array(istac,jstac,kstac,icpec)
icmax = istac
jcmax = jstac
kcmax = kstac

DO kc = kstac,kstoc
  DO jc = jstac,jstoc
    DO ic = istac,istoc
      
      IF(array(ic,jc,kc,icpec) < arrmin)THEN
        arrmin = array(ic,jc,kc,icpec)
        icmin = ic
        jcmin = jc
        kcmin = kc
      END IF
      IF(array(ic,jc,kc,icpec) > arrmax)THEN
        arrmax = array(ic,jc,kc,icpec)
        icmax = ic
        jcmax = jc
        kcmax = kc
      END IF
      
    END DO
  END DO
END DO

WRITE(6,'(2(1PE12.4,3I5))')arrmin,icmin,jcmin,kcmin, arrmax,icmax,jcmax,kcmax

!     =========================================================================


RETURN
END SUBROUTINE chkary
