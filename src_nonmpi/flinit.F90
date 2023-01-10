SUBROUTINE flinit
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-04  Time: 21:03:01

!     *************************************************************************

!     FLINIT
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT

!     CHANGE RECORD
!     -------------
!     30-AUG-2009:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     INITIALISES SPATIAL FILTERING
!     12TH ORDER EXPLICIT FILTERING
!     6TH,7TH,8TH,9TH,10TH,11TH ORDER EXPLICIT BOUNDARY SCHEMES

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     BEGIN
!     =====

!     =========================================================================

!     INTERIOR SCHEME
!     ---------------

!     TWELFTH ORDER EXPLICIT CENTRED DIFFERENCES
twfact = -EXP(-12.0_dp*LOG(two))
facofx =  792.0_dp*twfact
fbcofx = -495.0_dp*twfact
fccofx =  220.0_dp*twfact
fdcofx = -66.0_dp*twfact
fecofx =  12.0_dp*twfact
ffcofx = -1.00_dp*twfact
fgcofx = -924.0_dp*twfact
fgcofx =  1.00_dp + fgcofx

facofy = facofx
fbcofy = fbcofx
fccofy = fccofx
fdcofy = fdcofx
fecofy = fecofx
ffcofy = ffcofx
fgcofy = fgcofx

facofz = facofx
fbcofz = fbcofx
fccofz = fccofx
fdcofz = fdcofx
fecofz = fecofx
ffcofz = ffcofx
fgcofz = fgcofx


!     BOUNDARY TREATMENT
!     ------------------

!     FIRST POINT SCHEME (6TH ORDER ONE SIDED)
facf1x = -1.00_dp*twfact
fbcf1x =  6.00_dp*twfact
fccf1x = -15.0_dp*twfact
fdcf1x =  20.0_dp*twfact
fecf1x = -15.0_dp*twfact
ffcf1x =  6.00_dp*twfact
fgcf1x = -1.00_dp*twfact
facf1x =  1.00_dp + facf1x

facf1y = facf1x
fbcf1y = fbcf1x
fccf1y = fccf1x
fdcf1y = fdcf1x
fecf1y = fecf1x
ffcf1y = ffcf1x
fgcf1y = fgcf1x

facf1z = facf1x
fbcf1z = fbcf1x
fccf1z = fccf1x
fdcf1z = fdcf1x
fecf1z = fecf1x
ffcf1z = ffcf1x
fgcf1z = fgcf1x

!     SECOND POINT SCHEME (7TH ORDER MIXED)
facf2x =  6.00_dp*twfact
fbcf2x = -37.0_dp*twfact
fccf2x =  96.0_dp*twfact
fdcf2x = -135.0_dp*twfact
fecf2x =  110.0_dp*twfact
ffcf2x = -51.0_dp*twfact
fgcf2x =  12.0_dp*twfact
fhcf2x = -1.00_dp*twfact
fbcf2x =  1.00_dp + fbcf2x

facf2y = facf2x
fbcf2y = fbcf2x
fccf2y = fccf2x
fdcf2y = fdcf2x
fecf2y = fecf2x
ffcf2y = ffcf2x
fgcf2y = fgcf2x
fhcf2y = fhcf2x

facf2z = facf2x
fbcf2z = fbcf2x
fccf2z = fccf2x
fdcf2z = fdcf2x
fecf2z = fecf2x
ffcf2z = ffcf2x
fgcf2z = fgcf2x
fhcf2z = fhcf2x


!     THIRD POINT SCHEME  (8TH ORDER MIXED)
facf3x = -15.0_dp*twfact
fbcf3x =  96.0_dp*twfact
fccf3x = -262.0_dp*twfact
fdcf3x =  396.0_dp*twfact
fecf3x = -360.0_dp*twfact
ffcf3x =  200.0_dp*twfact
fgcf3x = -66.0_dp*twfact
fhcf3x =  12.0_dp*twfact
ficf3x = -1.00_dp*twfact
fccf3x =  1.00_dp + fccf3x

facf3y = facf3x
fbcf3y = fbcf3x
fccf3y = fccf3x
fdcf3y = fdcf3x
fecf3y = fecf3x
ffcf3y = ffcf3x
fgcf3y = fgcf3x
fhcf3y = fhcf3x
ficf3y = ficf3x

facf3z = facf3x
fbcf3z = fbcf3x
fccf3z = fccf3x
fdcf3z = fdcf3x
fecf3z = fecf3x
ffcf3z = ffcf3x
fgcf3z = fgcf3x
fhcf3z = fhcf3x
ficf3z = ficf3x

!     FOURTH POINT SCHEME (9TH ORDER MIXED)
facf4x =  20.0_dp*twfact
fbcf4x = -135.0_dp*twfact
fccf4x =  396.0_dp*twfact
fdcf4x = -662.0_dp*twfact
fecf4x =  696.0_dp*twfact
ffcf4x = -480.0_dp*twfact
fgcf4x =  220.0_dp*twfact
fhcf4x = -66.0_dp*twfact
ficf4x =  12.0_dp*twfact
fjcf4x = -1.00_dp*twfact
fdcf4x =  1.00_dp + fdcf4x

facf4y = facf4x
fbcf4y = fbcf4x
fccf4y = fccf4x
fdcf4y = fdcf4x
fecf4y = fecf4x
ffcf4y = ffcf4x
fgcf4y = fgcf4x
fhcf4y = fhcf4x
ficf4y = ficf4x
fjcf4y = fjcf4x

facf4z = facf4x
fbcf4z = fbcf4x
fccf4z = fccf4x
fdcf4z = fdcf4x
fecf4z = fecf4x
ffcf4z = ffcf4x
fgcf4z = fgcf4x
fhcf4z = fhcf4x
ficf4z = ficf4x
fjcf4z = fjcf4x

!     FIFTH POINT SCHEME  (10TH ORDER MIXED)
facf5x = -15.0_dp*twfact
fbcf5x =  110.0_dp*twfact
fccf5x = -360.0_dp*twfact
fdcf5x =  696.0_dp*twfact
fecf5x = -887.0_dp*twfact
ffcf5x =  786.0_dp*twfact
fgcf5x = -495.0_dp*twfact
fhcf5x =  220.0_dp*twfact
ficf5x = -66.0_dp*twfact
fjcf5x =  12.0_dp*twfact
fkcf5x = -1.00_dp*twfact
fecf5x =  1.00_dp + fecf5x

facf5y = facf5x
fbcf5y = fbcf5x
fccf5y = fccf5x
fdcf5y = fdcf5x
fecf5y = fecf5x
ffcf5y = ffcf5x
fgcf5y = fgcf5x
fhcf5y = fhcf5x
ficf5y = ficf5x
fjcf5y = fjcf5x
fkcf5y = fkcf5x

facf5z = facf5x
fbcf5z = fbcf5x
fccf5z = fccf5x
fdcf5z = fdcf5x
fecf5z = fecf5x
ffcf5z = ffcf5x
fgcf5z = fgcf5x
fhcf5z = fhcf5x
ficf5z = ficf5x
fjcf5z = fjcf5x
fkcf5z = fkcf5x

!     SIXTH POINT SCHEME  (11TH ORDER MIXED)
facf6x =  6.00_dp*twfact
fbcf6x = -51.0_dp*twfact
fccf6x =  200.0_dp*twfact
fdcf6x = -480.0_dp*twfact
fecf6x =  786.0_dp*twfact
ffcf6x = -923.0_dp*twfact
fgcf6x =  792.0_dp*twfact
fhcf6x = -495.0_dp*twfact
ficf6x =  220.0_dp*twfact
fjcf6x = -66.0_dp*twfact
fkcf6x =  12.0_dp*twfact
flcf6x = -1.00_dp*twfact
ffcf6x =  1.00_dp + ffcf6x

facf6y = facf6x
fbcf6y = fbcf6x
fccf6y = fccf6x
fdcf6y = fdcf6x
fecf6y = fecf6x
ffcf6y = ffcf6x
fgcf6y = fgcf6x
fhcf6y = fhcf6x
ficf6y = ficf6x
fjcf6y = fjcf6x
fkcf6y = fkcf6x
flcf6y = flcf6x

facf6z = facf6x
fbcf6z = fbcf6x
fccf6z = fccf6x
fdcf6z = fdcf6x
fecf6z = fecf6x
ffcf6z = ffcf6x
fgcf6z = fgcf6x
fhcf6z = fhcf6x
ficf6z = ficf6x
fjcf6z = fjcf6x
fkcf6z = fkcf6x
flcf6z = flcf6x

!     =========================================================================


RETURN
END SUBROUTINE flinit
