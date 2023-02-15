FUNCTION espksq(waveno)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:16:38

!     *************************************************************************

!     ESPKSQ
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
!     EVALUATES THE TURBULENT ENERGY SPECTRUM TIMES WAVENUMBER SQUARED
!     AT THE GIVEN WAVENUMBER MAGNITUDE

!     *************************************************************************
use data_types

!     EXTERNAL FUNCTION
!     =================

REAL(KIND=dp), INTENT(IN)             :: waveno
REAL(KIND=dp) :: espect
EXTERNAL espect


!     FUNCTION
!     ========
REAL(KIND=dp) :: espksq


!     ARGUMENT
!     ========



!     BEGIN
!     =====

!     =========================================================================

espksq = waveno*waveno*espect(waveno)

!     =========================================================================


RETURN
END FUNCTION espksq
