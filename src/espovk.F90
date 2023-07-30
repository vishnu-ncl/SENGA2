FUNCTION espovk(waveno)

!   *************************************************************************

!   ESPOVK
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   20-APR-1997:  CREATED
!   24-MAY-2003:  RSC UPDATED FOR SENGA2

!   DESCRIPTION
!   -----------
!   EVALUATES THE TURBULENT ENERGY SPECTRUM DIVIDED BY WAVENUMBER
!   AT THE GIVEN WAVENUMBER MAGNITUDE

!   *************************************************************************

!   EXTERNAL FUNCTION
!   =================
    real(kind=8), intent(IN)             :: waveno
    real(kind=8) :: espect
    EXTERNAL espect

!   PARAMETER
!   =========
    real(kind=8), parameter :: zero = 0.0_8

!   FUNCTION
!   ========
    real(kind=8) :: espovk

!   ARGUMENT
!   ========

!   BEGIN
!   =====

!   =========================================================================

    espovk = zero
    IF(waveno > zero) espovk = espect(waveno)/waveno

!   =========================================================================

END FUNCTION espovk
