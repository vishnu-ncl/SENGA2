FUNCTION espksq(waveno)

!   *************************************************************************

!   ESPKSQ
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
!   EVALUATES THE TURBULENT ENERGY SPECTRUM TIMES WAVENUMBER SQUARED
!   AT THE GIVEN WAVENUMBER MAGNITUDE

!   *************************************************************************

!   EXTERNAL FUNCTION
!   =================
    real(kind=8), intent(IN)             :: waveno
    real(kind=8) :: espect
    EXTERNAL espect

!   FUNCTION
!   ========
    real(kind=8) :: espksq

!   ARGUMENT
!   ========

!   BEGIN
!   =====

!   =========================================================================

    espksq = waveno*waveno*espect(waveno)

!   =========================================================================

END FUNCTION espksq
