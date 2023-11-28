SUBROUTINE ardump(array,nxphys,nyphys,nzphys)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 21:15:34

!     *************************************************************************

!     ARDUMP
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     23-JUL-1992:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     DIAGNOSTIC ROUTINE
!     PRINTS THE CONTENTS OF THE SPECIFIED ARRAY

!     *************************************************************************



!     ARGUMENTS
!     =========

REAL(kind=8), INTENT(IN OUT)         :: array(nxphys,nyphys,nzphys)
INTEGER, INTENT(IN OUT)                  :: nxphys
INTEGER, INTENT(IN)                      :: nyphys
INTEGER, INTENT(IN)                      :: nzphys




!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc


!     BEGIN
!     =====

WRITE(6,*)

DO kc = nzphys,1,-1
  WRITE(6,'(" Z-PLANE ",I10)')kc
  DO jc = nyphys,1,-1
    WRITE(6,'(I5,21(1PE10.3))')jc,(array(ic,jc,kc),ic = 1,nxphys)
  END DO
END DO


RETURN
END SUBROUTINE ardump
