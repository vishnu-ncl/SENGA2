FUNCTION espovk(waveno)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:16:44

!     *************************************************************************

!     ESPOVK
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     20-APR-1997:  CREATED
!     24-MAY-2003:  RSC UPDATED FOR SENGA2

!     DESCRIPTION
!     -----------
!     EVALUATES THE TURBULENT ENERGY SPECTRUM DIVIDED BY WAVENUMBER
!     AT THE GIVEN WAVENUMBER MAGNITUDE

!     *************************************************************************


!     EXTERNAL FUNCTION
!     =================

REAL(kind=8), INTENT(IN)             :: waveno
REAL(kind=8) :: espect
EXTERNAL espect


!     PARAMETER
!     =========

REAL(kind=8), PARAMETER :: zero = 0.0_8


!     FUNCTION
!     ========
REAL(kind=8) :: espovk


!     ARGUMENT
!     ========



!     BEGIN
!     =====

!     =========================================================================

espovk = zero
IF(waveno > zero)espovk = espect(waveno)/waveno

!     =========================================================================


RETURN
END FUNCTION espovk
