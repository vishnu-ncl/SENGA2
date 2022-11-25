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
    use data_types

!   EXTERNAL FUNCTION
!   =================
    real(kind=dp), intent(IN)             :: waveno
    real(kind=dp) :: espect
    EXTERNAL espect

!   PARAMETER
!   =========
    real(kind=dp), parameter :: zero = 0.0_dp

!   FUNCTION
!   ========
    real(kind=dp) :: espovk

!   ARGUMENT
!   ========

!   BEGIN
!   =====

!   =========================================================================

    espovk = zero
    IF(waveno > zero) espovk = espect(waveno)/waveno

!   =========================================================================

END FUNCTION espovk
