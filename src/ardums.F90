SUBROUTINE ardums(array,nxphys,nyphys,nzphys,nspec,ispec)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2022-09-13  Time: 20:58:16

!     *************************************************************************

!     ARDUMS
!     ======

!     AUTHOR
!     ------
!     R.S.CANT

!     CHANGE RECORD
!     -------------
!     11-SEP-2005:  CREATED

!     DESCRIPTION
!     -----------
!     DNS CODE SENGA2
!     DIAGNOSTIC ROUTINE
!     PRINTS THE CONTENTS OF THE SPECIFIED SPECIES ARRAY

!     *************************************************************************


!     ARGUMENTS
!     =========

INTEGER, INTENT(IN OUT)                  :: nxphys
INTEGER, INTENT(IN)                      :: nyphys
INTEGER, INTENT(IN)                      :: nzphys
INTEGER, INTENT(IN OUT)                  :: nspec
INTEGER, INTENT(IN OUT)                  :: ispec
real(kind=8),INTENT(IN OUT)         :: array(nxphys,nyphys,nzphys,nspec)




!     LOCAL DATA
!     ==========
INTEGER :: ic,jc,kc


!     BEGIN
!     =====

WRITE(6,*)

DO kc = nzphys,1,-1
  WRITE(6,'(" Z-PLANE ",I10)')kc
  DO jc = nyphys,1,-1
    WRITE(6,'(I5,21(1PE10.3))') jc,(array(ic,jc,kc,ispec),ic = 1,nxphys)
  END DO
END DO


RETURN
END SUBROUTINE ardums
