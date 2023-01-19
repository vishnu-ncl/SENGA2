SUBROUTINE exhals(bigarr,ispec,buffer,indexl,indexr,  &
        jndexl,jndexr,  &
        kndexl,kndexr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-05  Time: 14:07:55

!     *************************************************************************

!     EXHALS
!     ======

!     AUTHOR
!     ------
!     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPERTMENT

!     CHANGE RECORD
!     -------------
!     16-MAY-2003:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     EXTRACTS HALO DATA FOR PARALLEL TRANSFER BUFFER
!     FOR SPECIES MASS FRACTIONS

!     *************************************************************************


!     GLOBAL DATA
!     ===========
!     -------------------------------------------------------------------------
use data_types
use com_senga
!     -------------------------------------------------------------------------


!     ARGUMENTS
!     =========

REAL(KIND=8), INTENT(IN)             :: bigarr(nspcmx,nxbigl:nxbigr,nybigl:nybigr,nzbigl:nzbigr)
INTEGER, INTENT(IN OUT)                  :: ispec
REAL(KIND=8), INTENT(OUT)            :: buffer(nparay)
INTEGER, INTENT(IN)                      :: indexl
INTEGER, INTENT(IN)                      :: indexr
INTEGER, INTENT(IN)                      :: jndexl
INTEGER, INTENT(IN)                      :: jndexr
INTEGER, INTENT(IN)                      :: kndexl
INTEGER, INTENT(IN)                      :: kndexr






!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc
INTEGER :: icount


!     BEGIN
!     =====

!     =========================================================================

icount = 0
DO kc = kndexl,kndexr
  DO jc = jndexl,jndexr
    DO ic = indexl,indexr
      
      icount = icount + 1
      buffer(icount) = bigarr(ispec,ic,jc,kc)
      
    END DO
  END DO
END DO

!     =========================================================================


RETURN
END SUBROUTINE exhals
