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
    use data_types

!   EXTERNAL FUNCTION
!   =================
    real(kind=dp), intent(IN)             :: waveno
    real(kind=dp) :: espect
    EXTERNAL espect

!   FUNCTION
!   ========
    real(kind=dp) :: espksq

!   ARGUMENT
!   ========

!   BEGIN
!   =====

!   =========================================================================

    espksq = waveno*waveno*espect(waveno)

!   =========================================================================

END FUNCTION espksq
