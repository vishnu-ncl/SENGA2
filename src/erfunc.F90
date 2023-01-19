FUNCTION erfunc(argmnt)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:16:24

!     *************************************************************************

!     ERFUNC
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     09-JAN-2005:  CREATED

!     DESCRIPTION
!     -----------
!     COMPUTES THE ERROR FUNCTION
!     REF: ABRAMOWITZ AND STEGUN p299 sec 7.1.26

!     *************************************************************************
use data_types

!     PARAMETERS
!     ==========

REAL(KIND=8), INTENT(IN OUT)         :: argmnt

REAL(KIND=8), PARAMETER :: zero=0.0_8
REAL(KIND=8), PARAMETER :: one=1.0_8
REAL(KIND=8), PARAMETER :: half=0.50_8


INTEGER, PARAMETER :: ncoeff = 5
INTEGER, PARAMETER :: ncofm1 = ncoeff-1


!     FUNCTION
!     ========
REAL(KIND=8) :: erfunc


!     ARGUMENT
!     ========



!     LOCAL DATA
!     ==========
REAL(KIND=8) :: ecoeff(ncoeff)
REAL(KIND=8) :: pcoeff
REAL(KIND=8) :: etotal,zvalue,tvalue
INTEGER :: icoeff


!     BEGIN
!     =====

!     =========================================================================

!     SET THE COEFFICIENTS
pcoeff = 0.3275911_8
ecoeff(1) = 0.254829592_8
ecoeff(2) =-0.284496736_8
ecoeff(3) = 1.421413741_8
ecoeff(4) =-1.453152027_8
ecoeff(5) = 1.061405429_8

!     EVALUATE ERROR FUNCTION
zvalue = ABS(argmnt)
tvalue = one/(one+pcoeff*zvalue)

etotal = ecoeff(ncoeff)
DO icoeff = ncofm1,1,-1
  etotal = ecoeff(icoeff) + etotal*tvalue
END DO
etotal = etotal*tvalue

erfunc = one - etotal*EXP(-zvalue*zvalue)
IF(argmnt < zero)erfunc = -erfunc

!     =========================================================================


RETURN
END FUNCTION erfunc
