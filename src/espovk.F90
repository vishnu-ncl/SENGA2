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
    real(8), intent(IN)             :: waveno
    real(8) :: espect
    EXTERNAL espect

!   PARAMETER
!   =========
    real(8), parameter :: zero = 0.0_8

!   FUNCTION
!   ========
    real(8) :: espovk

!   ARGUMENT
!   ========

!   BEGIN
!   =====

!   =========================================================================

    espovk = zero
    IF(waveno > zero) espovk = espect(waveno)/waveno

!   =========================================================================

END FUNCTION espovk
