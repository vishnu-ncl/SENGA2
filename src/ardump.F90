SUBROUTINE ardump(array,nxphys,nyphys,nzphys)

!   *************************************************************************

!   ARDUMP
!   ======

!   AUTHOR
!   ------
!   R.S.CANT

!   CHANGE RECORD
!   -------------
!   23-JUL-1992:  CREATED

!   DESCRIPTION
!   -----------
!   DNS CODE SENGA2
!   DIAGNOSTIC ROUTINE
!   PRINTS THE CONTENTS OF THE SPECIFIED ARRAY

!   *************************************************************************

!   ARGUMENTS
!   =========

    integer(kind=4), intent(IN OUT)             :: nxphys
    integer(kind=4), intent(IN)                 :: nyphys
    integer(kind=4), intent(IN)                 :: nzphys
    real(kind=8), intent(IN OUT)       :: array(nxphys,nyphys,nzphys)

!   LOCAL DATA
!   ==========
    integer(kind=4) :: ic,jc,kc

!   BEGIN
!   =====

    WRITE(6,*)

    DO kc = nzphys,1,-1
        WRITE(6,'(" Z-PLANE ",I10)')kc
        DO jc = nyphys,1,-1
            WRITE(6,'(I5,21(1PE10.3))')jc,(array(ic,jc,kc),ic = 1,nxphys)
        END DO
    END DO

!   ==========================================================================

END SUBROUTINE ardump
