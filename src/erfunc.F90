FUNCTION erfunc(argmnt)

!   *************************************************************************

!   ERFUNC
!   ======

!   AUTHOR
!   ------
!   R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!   CHANGE RECORD
!   -------------
!   09-JAN-2005:  CREATED

!   DESCRIPTION
!   -----------
!   COMPUTES THE ERROR FUNCTION
!   REF: ABRAMOWITZ AND STEGUN p299 sec 7.1.26

!   *************************************************************************

!   PARAMETERS
!   ==========
    real(kind=8), intent(in out)         :: argmnt

    real(kind=8), parameter:: zero=0.0_8
    real(kind=8), parameter:: one=1.0_8
    real(kind=8), parameter:: half=0.50_8


    integer(kind=4), parameter:: ncoeff = 5
    integer(kind=4), parameter:: ncofm1 = ncoeff-1

!   FUNCTION
!   ========
    real(kind=8) :: erfunc

!   ARGUMENT
!   ========

!   LOCAL DATA
!   ==========
    real(kind=8) :: ecoeff(ncoeff)
    real(kind=8) :: pcoeff
    real(kind=8) :: etotal,zvalue,tvalue
    integer(kind=4) :: icoeff

!   BEGIN
!   =====

!   =========================================================================

!   SET THE COEFFICIENTS
    pcoeff = 0.3275911_8
    ecoeff(1) = 0.254829592_8
    ecoeff(2) =-0.284496736_8
    ecoeff(3) = 1.421413741_8
    ecoeff(4) =-1.453152027_8
    ecoeff(5) = 1.061405429_8

!   EVALUATE ERROR FUNCTION
    zvalue = ABS(argmnt)
    tvalue = one/(one+pcoeff*zvalue)

    etotal = ecoeff(ncoeff)
    DO icoeff = ncofm1,1,-1
        etotal = ecoeff(icoeff) + etotal*tvalue
    END DO
    etotal = etotal*tvalue

    erfunc = one - etotal*EXP(-zvalue*zvalue)
    IF(argmnt < zero) erfunc = -erfunc

!   =========================================================================

END FUNCTION erfunc
