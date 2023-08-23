FUNCTION espect(waveno)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-14  Time: 11:16:31

!     *************************************************************************

!     ESPECT
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  DEPARTMENT OF MECHANICAL ENGINEERING, UMIST

!     CHANGE RECORD
!     -------------
!     28-JUN-1994:  CREATED
!     13-JUL-2001:  RSC NEW VERSION FOR BATCHELOR-TOWNSEND SPECTRUM
!     24-MAY-2003:  RSC UPDATED FOR SENGA2

!     DESCRIPTION
!     -----------
!     EVALUATES THE TURBULENT ENERGY SPECTRUM
!     AT THE GIVEN WAVENUMBER MAGNITUDE

!     REFERENCES
!     ----------
!     1) BATCHELOR, G.K., TOWNSEND, A.A.: JFM 88(4), 685-709, 1948.

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------

use com_espect
!     -------------------------------------------------------------------------


!     PARAMETERS
!     ==========

real(kind=8),INTENT(IN)             :: waveno

real(kind=8),PARAMETER :: two = 2.0_8


!     FUNCTION
!     ========
real(kind=8):: espect


!     ARGUMENT
!     ========



!     LOCAL DATA
!     ==========
real(kind=8):: wavrat


!     BEGIN
!     =====

!     =========================================================================

wavrat = waveno*ovk0
wavrat = wavrat*wavrat
espect = covk0*wavrat*wavrat*EXP(-two*wavrat)

!     =========================================================================


RETURN
END FUNCTION espect
